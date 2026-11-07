# Multigrid na web (Flask)

    pip install -r requirements.txt
    python3 app.py            # abre em http://localhost:5000

Requer `gcc` e `make`. O app compila sozinho um executável por (problema, dimensão)
em `_cache/`, sem tocar no seu `config.mk` nem no seu `./multigrid`.
Se os fontes estiverem em outra pasta: `MULTIGRID_DIR=/caminho/Multigrid python3 app.py`.

Controles: problema, dimensão, algoritmo (completo / V-cycle / V-cycle 1 nível),
malha n = 16·2^k, passos de tempo (calor), tolerância (Poisson) e o botão
"Ordem de convergência" (equivalente ao `make valid`: roda n e 2n).
Limites de malha por configuração: dicionário KMAX em app.py.

Parâmetros do problema (sliders): K e EXP (Poisson); K, α e tempo final T (calor).
Eles são `#define` nos headers, então cada combinação é compilada uma vez
(alguns segundos) em `_cache/<problema>_<dim>_<parâmetros>/`. Os headers originais
não são alterados. Apague `_cache/` quando quiser liberar espaço.
Faixas e valores de T: `FAIXAS` e `TF_VALORES` em app.py.
