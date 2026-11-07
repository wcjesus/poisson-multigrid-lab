COMECE A LEITURA POR AQUI, MAS DEPOIS VEJA PROBLEMAS.MD, QUE ENUMERA LIMITAÇÕES CONHECIDAS DO CÓDIGO E EXTENSOES.MD, QUE LISTA AS POSSÍVEIS MELHORIAS ENXERGADAS
USE O CÓDIGO DO PROJETO MENOR (lab.c) COMO LABORATÓRIO PARA A IMPLEMENTAÇÃO DESSE MAIOR

Esse projeto implementa um Multigrid geométrico para resolver equação do calor ou problema de Poisson em diferentes dimensões.
Existem diversos arquivos c e headers, além de Scripts bash, cuja utilidade varia. make clean sempre limpa todos os artefatos de texto, ou saídas que o programa gera, o que deve ser levado em conta antes de acionar.
O programa sempre compila um único executável. Isso é importante: a configuração do problema é definida na compilação.
Mas as configurações dele podem ser mudadas antes

Por exemplo:
POISSON 2D
Pode mudar para
EQUACAO DO CALOR 3D

Deve-se saber como sempre qual a configuração atual antes de compilar. 
Também como mudar a configuração.

Use make help para esse fim.

Todos os alvos do Makefile também serão explicados com esse comando.

Deve-se saber como validar o algoritmo. Existem três tipos de validacao

- Validação visual: make plot faz o gráfico da solução numérica do problema, da analítica e do erro para comparação 

- Validação da convergência: mais importante que o gráfico, make valid faz a tabela de convergência, que mostra a ordem de convergência assintotica, buscando demonstrar que a discretizacao do domínio está corretamente implementada. Para espaço, é diferença finita central (segunda ordem: p -> 2), para tempo, Euler implícito (primeira ordem: p -> 1). 

- O próprio algoritmo calcula e imprime o erro final obtido apos a convergência usando duas normas.

O make help imprime instruções de execução do programa. Além disso, executar o programa sem argumentos da seguinte forma:

./multigrid 

Mostra a sintaxe correta para cada modo.

Adicionalmente, para se acompanhar a convergência, pode-se ligar e desligar as impressões 

"Ciclo %d: Resíduo: %.12e"

Na tela, ou ainda fazer com que a execução grave todo o caminho percorrido pelo código em um arquivo.

Para medição de desempenho, o correto é executar sem essas saídas. Para acompanhamento e manutenção do funcionamento, aí as impressões são desejáveis.

Na hora da execução, pode-se escolher com que o programa faça Two Grid e Vcycle puro além do ciclo completo do Multigrid. 

Para experimentos, recomenda-se variar os parâmetros do problema, a malha espacial, o intervalo de tempo ou número de passos de tempo (quando aplicável, na Equação do calor) e o algoritmo entre as execuções. 

IMPORTANTE: Ao realizar testes demorados, deve-se verificar se a demora é devida ao dimensionamento do programa e não a paginação/swap de memoria
