Aqui buscamos documentar aspectos possivelmente problemáticos do código.

- No modo DEBUG, existe abertura de arquivo que não é fechado corretamente. Além disso, ainda que o texto gerado diga "erro suavizado", o erro calculado ali refere-se à solução do problema Ae=r na malha grossa, enquanto o resultado impresso é obtido a partir do cálculo do módulo da diferença entre a solução obtida e a solução fabricada. Essa comparação não tem sentido.
O modo foi concebido como uma forma de se acompanhar toda a execução do programa para depuração. No passado, já teve serventia em apontar a causa exata de comportamentos indevidos do algoritmo observados nas execuções. Foi mantido como legado, mas recomenda-se evitar seu uso, ou então corrigir seus problemas antes de usar. Ao menos tomar ciência deles.

- Muitas instruções de pré processador espalhadas comprometem a legibilidade do código. Foi almejado, mas não concluído, mudar a lógica para que diferentes executáveis fossem compilados e a configuração fosse trocado a partir das execução, em vez da compilação. 

- Execuções de make plot, make valid e outros apagam a compilação do usuário. Portanto deve-se tomar cuidado, ou corrigir para que eles gerem executáveis diferentes para fazerem suas tarefas e ao final apaguem apenas estes, sem eliminar uma compilação precedente. 

- Possivelmente novas otimizações no gerenciamento de memória são ainda possíveis.
