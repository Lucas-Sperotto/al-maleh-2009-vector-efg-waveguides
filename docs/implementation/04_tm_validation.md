# M4 — Correspondência modal e modos espúrios (experimental)

Esta etapa compara o espectro calculado para o campo magnético transversal com as soluções analíticas TM de guia retangular PEC. **Ainda não é uma reprodução integral verificada de Al Maleh (2009).** A revisão da teoria original ficou para uma etapa posterior.

## Critério de forma

Para o domínio [0,2] × [0,1], os valores analíticos são:

$$
\lambda_{mn} = \pi^2(m^2/4+n^2), \quad m,n\geq1.
$$

O rotacional do campo H é proporcional ao campo elétrico longitudinal Ez de um modo TM; portanto comparamos curl(H) com a função escalar senoidal:

$$
f_{mn}(x,y)=\sin(m\pi x/2)\sin(n\pi y).
$$

A sobreposição é o cosseno em L² entre estas funções:

$$
\mathrm{overlap}_{mn}(h)=
\frac{|\int_\Omega f_{mn}\,\mathrm{curl}(H_h)\,d\Omega|}
{\|f_{mn}\|_{L^2}\,\|\mathrm{curl}(H_h)\|_{L^2}}.
$$

A forma numérica é obtida diretamente das derivadas MLS. O produto interno e a matriz de rigidez usam a **mesma quadratura**, para consistência. Em autovalores degenerados, consideramos a norma da **projeção no subespaço** de autovetores quase degenerados, não só um vetor isolado. Isso evita falsos negativos decorrentes da rotação arbitrária de autovetores degenerados.

## Comparações e critérios

Para cada TM analítico, registramos separadamente: índice espectral mais próximo, autovalor e erro relativo; índice com **maior similaridade de forma** e seu erro relativo; sobreposição individual e de cluster; resíduo algébrico; divergência relativa; indicador de fluxo normal nas paredes; condicionamento da massa e dimensão numérica do núcleo curl-curl.

Uma triagem exploratória chamada shape_supported_candidate exige simultaneamente erro relativo de lambda <= 0,15, sobreposição de cluster >= 0,75, divergência relativa <= 0,25, indicador normal <= 0,25 e resíduo <= 1e-7. Os limiares são **escolhas de pesquisa**: não são regras da dissertação, nem provas de presença ou ausência de modos espúrios. Os casos restantes são rotulados como unconfirmed_candidate, **nunca spurious_proven**.

Uma aproximação pode apresentar pequeno resíduo algébrico e ainda não representar um modo físico. Também é possível que um modo válido esteja mal resolvido numa nuvem grosseira; a classificação precisa ser analisada em conjunto com refinamento e comparação por subespaços.

## Campanha automatizada

Compilar e executar:

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
python3 scripts/m4_campaign.py --binary build/tm_validation --out out/m4 --profile ci
~~~

Casos executados: três refinamentos (5×5, 7×7, 9×7, MLS grau 1, suporte 3,3 e Gauss 3); teste de suportes 2,8 e 3,8 com 7×7; teste MLS quadrático (7×7, suporte 3,8); quadratura Gauss 2 e 4 versus Gauss 3 de referência.

Cada execução guarda tabelas de correspondência, espectro, metadados, stdout e stderr. O relatório da campanha agrega runs.csv, tm_matches_all.csv e REPORT.md. Falhas numéricas **são reportadas**, inclusive posto ou condicionamento insuficiente da massa. A CI exige as três execuções de refinamento concluídas, mas falhas dos ensaios de sensibilidade são registradas como resultados válidos de diagnóstico.

Limitações: o solucionador ainda é denso, limitado a 600 DOFs; a validação inicial utiliza geometria retangular e campos TM com índice m,n pequeno; os indicadores de fronteira e divergência são discretos e dependem da quadratura. Uma conclusão sobre convergência física ou modos espúrios exige campanhas mais amplas e comparação independente com FEM.

Código: src/tm_validation.cpp, include/al_maleh/tm_validation.hpp; aplicativo apps/tm_validation.cpp; teste tests/test_tm_validation.cpp; orquestrador scripts/m4_campaign.py.
