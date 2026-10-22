# Proveniência e licenças

Este núcleo de pesquisa reaproveita conceitos, decisões estruturais e o cálculo de primeiras derivadas MLS do repositório próprio [EFGMI-H(curl)](https://github.com/Lucas-Sperotto/efgmi-hcurl-diffuse-whitney), especialmente `src/mls/mls_evaluator.cpp`, sob **licença MIT**, copyright (c) 2026 Sperotto.

A presente adaptação foi reescrita na API `al_maleh::MLSEvaluator`, acrescentando derivadas de segunda ordem e peso Wendland C2. Na auditoria de 2026-10-22 a base passou de escala global fixa para base deslocada ao ponto de avaliação (escala = maior raio de suporte) e foi acrescentado o peso spline quártica transcrito da Eq. (4-2) de Al Maleh (2009). Respeitar o `LICENSE` de ambos os repositórios.

Os outros repositórios NASA TP-3485 e `simple-waveguide-efg-fem-solver` são referências metodológicas. **Nenhuma rotina de montagem FEM/EFG foi copiada deles nesta etapa.** O `mls_shape.cpp` legado não deve ser usado sem correção das derivadas.

Este código **não** é identificado como código original de Hussam Al Maleh. As equações computadas são uma construção independente destinada à sua reprodução crítica.
