#!/usr/bin/env python3
"""
Interface web (Flask) para o executável Multigrid.

Estrutura esperada:

    pasta/
      app.py
      templates/index.html
      Multigrid/            <- código-fonte C (Makefile, *.c, *.h ...)

Para apontar para outra pasta de fontes:  MULTIGRID_DIR=/caminho python3 app.py

O app NÃO mexe no seu config.mk nem no seu ./multigrid: para cada combinação
(problema, dimensão) ele compila um executável próprio (com -DPLOT) numa
pasta de cache (_cache/), usando links simbólicos para os seus fontes, e
executa cada requisição num diretório temporário isolado.
"""
import base64
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import threading
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from flask import Flask, jsonify, render_template, request

app = Flask(__name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.environ.get("MULTIGRID_DIR") or next(
    (d for d in (os.path.join(BASE_DIR, "Multigrid"), BASE_DIR)
     if os.path.isfile(os.path.join(d, "Makefile"))),
    os.path.join(BASE_DIR, "Multigrid"),
)
CACHE_DIR = os.path.join(BASE_DIR, "_cache")

TIMEOUT_EXEC = 180          # segundos por execução
TIMEOUT_BUILD = 240
_sem_exec = threading.Semaphore(2)          # no máximo 2 execuções simultâneas
_trava_build = {}                           # uma trava por configuração
_trava_global = threading.Lock()

# k máximo da malha n = 16 * 2**k, por (problema, dimensão)
KMAX = {
    ("POISSON", "1D"): 10,   # 16384
    ("POISSON", "2D"): 6,    # 1024
    ("POISSON", "3D"): 3,    # 128
    ("CALOR",   "1D"): 9,    # 8192
    ("CALOR",   "2D"): 5,    # 512  (o modo PLOT grava vários quadros .dat)
    ("CALOR",   "3D"): 2,    # 64
}
PASSOS_MAX = 200
TOL_MAX = 14


# ----------------------------------------------------------------------
# Utilidades
# ----------------------------------------------------------------------
def _define(arquivo, nome, padrao):
    """Lê '#define NOME valor' de um header."""
    try:
        with open(os.path.join(SRC_DIR, arquivo), errors="replace") as f:
            for linha in f:
                p = linha.split()
                if len(p) >= 3 and p[0] == "#define" and p[1] == nome:
                    return float(p[2])
    except OSError:
        pass
    return padrao


# valores padrão = os que estão nos headers originais
PADRAO = {
    "K": int(_define("problema.h", "K", 4)),          # Poisson
    "EXP": int(_define("problema.h", "EXP", 2)),      # Poisson
    "KC": int(_define("calor.h", "K", 10)),           # calor
    "ALFA": _define("calor.h", "ALFA", 4.0),          # calor
    "TF": _define("calor.h", "TEMPO_FINAL", 0.1),     # calor
}
# faixas dos sliders
FAIXAS = {"K": (1, 8), "EXP": (1, 3), "KC": (1, 20), "ALFA": (0.0, 20.0)}
ALFA_PASSO = 0.5
# tempos finais permitidos (slider logarítmico: índice -> valor)
TF_VALORES = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
if PADRAO["TF"] not in TF_VALORES:
    TF_VALORES = sorted(TF_VALORES + [PADRAO["TF"]])


def soma(x, K, EXP):
    s, freq = 0.0, K
    for _ in range(EXP):
        s = s + np.sin(freq * np.pi * x)
        freq *= K
    return s


def _num(x):
    """Formata número como literal C e como pedaço de nome de pasta."""
    return repr(float(x))


def _slug(x):
    return _num(x).replace("-", "m").replace(".", "p")


def _int(nome, valor, minimo, maximo):
    try:
        n = int(valor)
    except (TypeError, ValueError):
        raise ValueError(f"{nome} deve ser inteiro.")
    if not (minimo <= n <= maximo):
        raise ValueError(f"{nome} deve estar entre {minimo} e {maximo}.")
    return n


def _escolha(nome, valor, opcoes):
    if valor not in opcoes:
        raise ValueError(f"{nome} inválido.")
    return valor


# ----------------------------------------------------------------------
# Compilação (uma vez por problema/dimensão)
# ----------------------------------------------------------------------
def _headers_modificados(problema, p):
    """Devolve {arquivo: texto} com os #define de K, EXP, ALFA, TEMPO_FINAL trocados."""
    if problema == "POISSON":
        arq, troca = "problema.h", {"K": str(int(p["K"])), "EXP": str(int(p["EXP"]))}
    else:
        arq, troca = "calor.h", {"K": _num(p["KC"]), "ALFA": _num(p["ALFA"]),
                                 "TEMPO_FINAL": _num(p["TF"])}
    with open(os.path.join(SRC_DIR, arq), errors="replace") as f:
        texto = f.read()
    for nome, valor in troca.items():
        texto, n = re.subn(rf"(?m)^(\s*#define\s+{nome}\s+).*?(\r?)$",
                           lambda m: f"{m.group(1)}{valor}{m.group(2)}", texto)
        if n != 1:
            raise RuntimeError(f"#define {nome} não encontrado em {arq}.")
    return {arq: texto}


def _garante_executavel(problema, dim, p):
    if problema == "POISSON":
        chave = f"{problema}_{dim}_K{int(p['K'])}_E{int(p['EXP'])}"
    else:
        chave = (f"{problema}_{dim}_K{_slug(p['KC'])}_A{_slug(p['ALFA'])}"
                 f"_T{_slug(p['TF'])}")
    modificados = _headers_modificados(problema, p)
    with _trava_global:
        trava = _trava_build.setdefault(chave, threading.Lock())
    pasta = os.path.join(CACHE_DIR, chave)
    with trava:
        os.makedirs(pasta, exist_ok=True)
        # links simbólicos para os fontes (make detecta alterações neles)
        for arq in glob.glob(os.path.join(SRC_DIR, "*.c")) + \
                   glob.glob(os.path.join(SRC_DIR, "*.h")) + \
                   [os.path.join(SRC_DIR, "Makefile")]:
            nome = os.path.basename(arq)
            if nome in modificados:
                continue
            destino = os.path.join(pasta, nome)
            if not os.path.lexists(destino):
                os.symlink(arq, destino)
        # headers com parâmetros: cópia própria desta configuração
        for nome, texto in modificados.items():
            destino = os.path.join(pasta, nome)
            if os.path.islink(destino):
                os.remove(destino)
            antigo = open(destino, newline="").read() if os.path.exists(destino) else None
            if antigo != texto:
                with open(destino, "w", newline="") as f:
                    f.write(texto)
        with open(os.path.join(pasta, "config.mk"), "w") as f:
            f.write(f"PROBLEMA = {problema}\nDIMENSAO = {dim}\n"
                    "PRINT = ON\nDEBUG = OFF\n")
        r = subprocess.run(["make", "-j4", "multigridplot"], cwd=pasta,
                           capture_output=True, text=True,
                           timeout=TIMEOUT_BUILD)
        if r.returncode != 0:
            raise RuntimeError("Falha ao compilar:\n" + (r.stderr or r.stdout)[-2000:])
    return os.path.join(pasta, "multigridplot")


def _roda(exe, args, cwd):
    with _sem_exec:
        r = subprocess.run([exe] + [str(a) for a in args], cwd=cwd,
                           capture_output=True, timeout=TIMEOUT_EXEC)
    out = r.stdout.decode("utf-8", errors="replace")
    err = r.stderr.decode("utf-8", errors="replace")
    return r.returncode, out, err


# ----------------------------------------------------------------------
# Leitura dos dados e figuras
# ----------------------------------------------------------------------
def _le(caminho):
    if not os.path.isfile(caminho) or os.path.getsize(caminho) == 0:
        raise RuntimeError(f"O executável não gerou {os.path.basename(caminho)}.")
    try:
        import pandas as pd
        return pd.read_csv(caminho, sep=r"\s+", header=None,
                           skip_blank_lines=True).to_numpy(dtype=float)
    except ImportError:
        return np.loadtxt(caminho)


def _png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _passo(n, alvo=400):
    return max(1, int(np.ceil(n / alvo)))


def _mapa(ax, x, y, Z, titulo, cmap="viridis", log=False, rotulo=""):
    """Z[iy, ix]."""
    if log:
        pos = Z[Z > 0]
        vmax = float(Z.max()) if Z.size else 1.0
        vmin = max(float(pos.min()) if pos.size else 1e-16, vmax * 1e-8, 1e-300)
        norma = LogNorm(vmin=vmin, vmax=max(vmax, vmin * 10))
        m = ax.pcolormesh(x, y, np.ma.masked_less_equal(Z, 0), shading="auto",
                          cmap="inferno", norm=norma)
    else:
        m = ax.pcolormesh(x, y, Z, shading="auto", cmap=cmap)
    ax.set_title(titulo, fontsize=10)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    cb = plt.colorbar(m, ax=ax, fraction=0.046, pad=0.04)
    if rotulo:
        cb.set_label(rotulo)


def figs_poisson(dim, n, pasta, K, EXP):
    arq = os.path.join(pasta, "grafico.dat")
    d = _le(arq)
    figuras = []

    if dim == "1D":
        x, u, e = d[:, 0], d[:, 1], d[:, 2]
        xf = np.linspace(0, 1, 4000)
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.6))
        a1.plot(xf, soma(xf, K, EXP), color="#2459a6", lw=2, label="Analítica")
        a1.plot(x, u, "--", color="#d62728", lw=1.4, label="Numérica")
        a1.set_xlabel("x"); a1.set_ylabel("u(x)"); a1.grid(alpha=.3)
        a1.legend(); a1.set_title("Solução", fontsize=10)
        ep = np.where(e > 0, e, np.nan)
        a2.semilogy(x, ep, color="#d62728", lw=1.2)
        a2.set_xlabel("x"); a2.set_ylabel("|u_h − u|"); a2.grid(alpha=.3, which="both")
        a2.set_title("Erro absoluto", fontsize=10)
        fig.suptitle(f"Poisson 1D — n = {n} — K = {K}, EXP = {EXP}")
        figuras.append(("Solução e erro", _png(fig)))

    elif dim == "2D":
        # arquivo: laço externo em x (índice i), interno em y (j)
        N = n
        X = d[:, 0].reshape(N, N)
        Y = d[:, 1].reshape(N, N)
        U = d[:, 2].reshape(N, N)
        E = d[:, 3].reshape(N, N)
        s = _passo(N)
        xs, ys = X[::s, 0], Y[0, ::s]
        U, E = U[::s, ::s].T, E[::s, ::s].T            # -> [iy, ix]
        A = np.outer(soma(ys, K, EXP), soma(xs, K, EXP))
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
        _mapa(ax[0], xs, ys, A, "Analítica", "RdBu_r")
        _mapa(ax[1], xs, ys, U, "Numérica", "RdBu_r")
        _mapa(ax[2], xs, ys, E, "Erro absoluto", log=True)
        fig.suptitle(f"Poisson 2D — n = {n} — K = {K}, EXP = {EXP}")
        figuras.append(("Solução analítica, numérica e erro", _png(fig)))

    else:  # 3D: corte z ≈ 0.5
        N = n
        D = d.reshape(N, N, N, d.shape[1])            # [i, j, k, col]
        # z = 0,5 é nó de sen(Kπz): escolhe o corte de maior amplitude
        k = int(np.argmax(np.abs(D[:, :, :, 5]).sum(axis=(0, 1))))
        xs = D[:, 0, k, 0]; ys = D[0, :, k, 1]
        U = D[:, :, k, 3].T; E = D[:, :, k, 4].T; A = D[:, :, k, 5].T
        zc = D[0, 0, k, 2]
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
        _mapa(ax[0], xs, ys, A, "Analítica", "RdBu_r")
        _mapa(ax[1], xs, ys, U, "Numérica", "RdBu_r")
        _mapa(ax[2], xs, ys, E, "Erro absoluto", log=True)
        fig.suptitle(f"Poisson 3D — n = {n} — corte z = {zc:.3f} — K = {K}, EXP = {EXP}")
        figuras.append(("Corte do volume", _png(fig)))
    return figuras


def figs_calor(dim, n, pasta, tf, KC, ALFA):
    instantes = [("inicio", "t = 0"), ("meio", f"t = T/2 = {tf / 2:g}"),
                 ("fim", f"t = T = {tf:g}")]
    dados = [_le(os.path.join(pasta, f"grafico_{i}.dat")) for i, _ in instantes]
    N = n + 1
    # (rótulo da figura, coluna do arquivo por dimensão, log?)
    col = {"1D": {"num": 1, "ex": 2, "err": 3},
           "2D": {"num": 2, "ex": 3, "err": 4},
           "3D": {"num": 3, "ex": 4, "err": 5}}[dim]
    itens = [("num", "Solução numérica", "u_h"),
             ("ex", "Solução exata", "u"),
             ("err", "Erro absoluto", "|u_h − u|")]
    figuras = []

    kcorte = 0
    if dim == "3D":
        D0 = dados[0].reshape(N, N, N, dados[0].shape[1])      # [k, j, i, col]
        # a exata em t=0 define o corte de maior amplitude
        kcorte = int(np.argmax(np.abs(D0[:, :, :, col["ex"]]).sum(axis=(1, 2))))

    def fatia(d, c):
        if dim == "2D":                       # laço externo y, interno x
            Z = d[:, c].reshape(N, N)         # [iy, ix]
            xs = d[:N, 0]; ys = d[::N, 1]
            return xs, ys, Z
        # 3D: laços k (z) > j (y) > i (x)
        D = d.reshape(N, N, N, d.shape[1])    # [k, j, i, col]
        k = kcorte
        return D[k, 0, :, 0], D[k, :, 0, 1], D[k, :, :, c]

    # escala de cor comum aos 3 instantes (para comparação)
    for chave, nome, rot in itens:
        c = col[chave]
        if dim == "1D":
            fig, ax = plt.subplots(1, 3, figsize=(15, 3.6), sharey=(chave != "err"))
            for a, d, (_, ti) in zip(ax, dados, instantes):
                y = d[:, c]
                if chave == "err":
                    a.semilogy(d[:, 0], np.where(y > 0, y, np.nan), color="#d62728")
                else:
                    a.plot(d[:, 0], y, color="#2459a6", lw=1.6)
                a.set_title(ti, fontsize=10); a.set_xlabel("x")
                a.set_ylabel(rot); a.grid(alpha=.3, which="both")
        else:
            fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
            todos = []
            for d in dados:
                xs, ys, Z = fatia(d, c)
                s = _passo(N)
                todos.append((xs[::s], ys[::s], Z[::s, ::s]))
            if chave != "err":
                vmin = min(float(t[2].min()) for t in todos)
                vmax = max(float(t[2].max()) for t in todos)
            for a, (xs, ys, Z), (_, ti) in zip(ax, todos, instantes):
                if chave == "err":
                    _mapa(a, xs, ys, Z, ti, log=True, rotulo=rot)
                else:
                    _mapa(a, xs, ys, Z, ti, "RdBu_r", rotulo=rot)
                    a.collections[0].set_clim(vmin, vmax)
        extra = ""
        if dim == "3D":
            zc = dados[0].reshape(N, N, N, dados[0].shape[1])[kcorte, 0, 0, 2]
            extra = f" — corte z = {zc:.3f}"
        fig.suptitle(f"Equação do calor {dim} — {nome} (n = {n}, K = {KC:g}, α = {ALFA:g}){extra}")
        figuras.append((nome, _png(fig)))
    return figuras


def fig_residuos(texto):
    """Histórico de resíduo a partir das linhas impressas (PRINT=ON)."""
    vals = []
    for lin in texto.splitlines():
        if lin.startswith("Ciclo") or lin.startswith("Tempo:"):
            m = re.search(r"[Rr]es[ií]duo:\s*([0-9.eE+\-]+)", lin)
            if m:
                try:
                    vals.append(float(m.group(1)))
                except ValueError:
                    pass
    if len(vals) < 2:
        return None
    fig, ax = plt.subplots(figsize=(7, 3.4))
    ax.semilogy(range(1, len(vals) + 1), vals, "o-", ms=3, color="#2459a6")
    ax.set_xlabel("ciclo V (acumulado)")
    ax.set_ylabel("resíduo")
    ax.grid(alpha=.3, which="both")
    ax.set_title("Histórico do resíduo", fontsize=10)
    return ("Convergência do resíduo", _png(fig))


def _resultados(texto):
    """Extrai o bloco RESULTADOS em pares (rótulo, valor)."""
    pares, dentro = [], False
    for lin in texto.splitlines():
        if lin.strip() == "RESULTADOS":
            dentro = True
            continue
        if dentro and ":" in lin:
            k, v = lin.split(":", 1)
            pares.append([k.strip(), v.strip()])
        elif dentro and lin.strip() and not pares:
            continue
    return pares


def _corta_log(texto, cabeca=40, cauda=60):
    linhas = texto.splitlines()
    if len(linhas) <= cabeca + cauda + 5:
        return texto
    omit = len(linhas) - cabeca - cauda
    return "\n".join(linhas[:cabeca] + [f"... ({omit} linhas omitidas) ..."]
                     + linhas[-cauda:])


# ----------------------------------------------------------------------
# Montagem da linha de comando
# ----------------------------------------------------------------------
def _prefixo(algoritmo):
    return {"fmg": [], "v": ["v"], "v1": ["v", "1"]}[algoritmo]


def _parametros(fonte):
    problema = _escolha("problema", fonte.get("problema", ""), ("POISSON", "CALOR"))
    dim = _escolha("dimensão", fonte.get("dimensao", ""), ("1D", "2D", "3D"))
    algoritmo = _escolha("algoritmo", fonte.get("algoritmo", "fmg"), ("fmg", "v", "v1"))
    kmax = KMAX[(problema, dim)]
    k = _int("malha", fonte.get("k", ""), 1 if algoritmo == "v1" else 0, kmax)
    n = 16 * 2 ** k
    passos = _int("passos", fonte.get("passos", "10"), 1, PASSOS_MAX)
    exp_tol = fonte.get("tol", "")
    exp_tol = None if exp_tol in ("", "auto") else _int("tolerância", exp_tol, 0, TOL_MAX)

    p = {}
    if problema == "POISSON":
        p["K"] = _int("K", fonte.get("K", PADRAO["K"]), *FAIXAS["K"])
        p["EXP"] = _int("EXP", fonte.get("EXP", PADRAO["EXP"]), *FAIXAS["EXP"])
        p["TF"] = PADRAO["TF"]
    else:
        p["KC"] = _int("K", fonte.get("KC", PADRAO["KC"]), *FAIXAS["KC"])
        try:
            alfa = float(fonte.get("ALFA", PADRAO["ALFA"]))
        except ValueError:
            raise ValueError("α deve ser numérico.")
        if not (FAIXAS["ALFA"][0] <= alfa <= FAIXAS["ALFA"][1]):
            raise ValueError("α fora da faixa.")
        p["ALFA"] = round(alfa / ALFA_PASSO) * ALFA_PASSO
        i = _int("tempo final", fonte.get("tf", TF_VALORES.index(PADRAO["TF"])),
                 0, len(TF_VALORES) - 1)
        p["TF"] = TF_VALORES[i]
    return problema, dim, algoritmo, k, n, passos, exp_tol, p


# ----------------------------------------------------------------------
# Rotas
# ----------------------------------------------------------------------
@app.route("/")
def index():
    return render_template(
        "index.html", kmax_js=json.dumps({f"{p}_{d}": v for (p, d), v in KMAX.items()}),
        passos_max=PASSOS_MAX, tol_max=TOL_MAX,
        padrao_js=json.dumps(PADRAO), faixas_js=json.dumps(FAIXAS),
        tf_js=json.dumps(TF_VALORES), tf_idx=TF_VALORES.index(PADRAO["TF"]), tf_max=len(TF_VALORES) - 1,
        alfa_passo=ALFA_PASSO)


@app.route("/api/executar")
def api_executar():
    try:
        problema, dim, alg, k, n, passos, exp_tol, p = _parametros(request.args)
        exe = _garante_executavel(problema, dim, p)

        if problema == "POISSON":
            args = _prefixo(alg) + ([0, n, exp_tol] if exp_tol is not None else [n])
        else:
            args = _prefixo(alg) + ["n", n, passos]

        tmp = tempfile.mkdtemp(prefix="mg_")
        try:
            rc, out, err = _roda(exe, args, tmp)
            if rc != 0 or "RESULTADOS" not in out:
                raise RuntimeError((out + err).strip() or "execução falhou")
            figuras = (figs_poisson(dim, n, tmp, p["K"], p["EXP"])
                       if problema == "POISSON"
                       else figs_calor(dim, n, tmp, p["TF"], p["KC"], p["ALFA"]))
            fr = fig_residuos(out)
            if fr:
                figuras.append(fr)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        return jsonify({
            "comando": "./multigrid " + " ".join(str(a) for a in args),
            "resultados": _resultados(out),
            "figuras": [{"titulo": t, "img": i} for t, i in figuras],
            "log": _corta_log(out),
        })
    except ValueError as e:
        return jsonify({"erro": str(e)}), 400
    except subprocess.TimeoutExpired:
        return jsonify({"erro": "Tempo limite excedido; reduza a malha ou os passos."}), 504
    except RuntimeError as e:
        return jsonify({"erro": str(e)}), 500


@app.route("/api/ordem")
def api_ordem():
    """Estudo de convergência (equivalente ao 'make valid'): n e 2n."""
    try:
        problema, dim, alg, k, n, passos, exp_tol, p = _parametros(request.args)
        tipo = _escolha("tipo", request.args.get("tipo", "espacial"),
                        ("espacial", "temporal"))
        if problema == "POISSON" and tipo == "temporal":
            raise ValueError("Poisson é estacionário: só há ordem espacial.")
        if k + 1 > KMAX[(problema, dim)]:
            raise ValueError("Para comparar n e 2n, escolha uma malha menor "
                             "(2n precisa caber no limite).")
        exe = _garante_executavel(problema, dim, p)
        pre = _prefixo(alg)

        if problema == "POISSON":
            args = pre + ([0, n, 2 * n, exp_tol] if exp_tol is not None else [n, 2 * n])
        elif tipo == "espacial":
            dt = p["TF"] / passos
            args = pre + ["s", n, 2 * n, repr(dt)]
        else:
            if 2 * passos > PASSOS_MAX * 2:
                raise ValueError("Passos demais.")
            args = pre + ["n", n, passos, 2 * passos]

        tmp = tempfile.mkdtemp(prefix="mg_")
        try:
            rc, out, err = _roda(exe, args, tmp)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        if rc != 0:
            raise RuntimeError((out + err).strip() or "execução falhou")
        # só o final: a(s) linha(s) da ordem e o bloco de resultados
        return jsonify({
            "comando": "./multigrid " + " ".join(str(a) for a in args),
            "saida": _corta_log(out, 10, 40),
        })
    except ValueError as e:
        return jsonify({"erro": str(e)}), 400
    except subprocess.TimeoutExpired:
        return jsonify({"erro": "Tempo limite excedido."}), 504
    except RuntimeError as e:
        return jsonify({"erro": str(e)}), 500


if __name__ == "__main__":
    if not os.path.isfile(os.path.join(SRC_DIR, "Makefile")):
        raise SystemExit(f"Makefile não encontrado em {SRC_DIR}. "
                         "Coloque a pasta Multigrid/ ao lado de app.py "
                         "ou defina MULTIGRID_DIR.")
    print(f"Fontes: {SRC_DIR}")
    app.run(host="0.0.0.0", port=5000, debug=False, threaded=True)
