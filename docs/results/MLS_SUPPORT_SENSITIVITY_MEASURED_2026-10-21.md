# Varredura controlada de suporte MLS M2/M5 — resultados medidos

**Fonte:** [GitHub Actions 39642651230](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39642651230) no commit 02eb257f (CI: **28/28 CTests passaram**, campanha 12/12 execuções de solver, 24 linhas de modos). [Artefato com todos os espectros, matrizes e CSVs](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39642651230/artifacts/12365842067).

## Contrato

Somente o **fator radial de suporte MLS** é alterado: 3,0 / 3,3 / 3,6. Constantes: nuvens 3×5 e 4×7, geometria 2×1 m, MLS linear, Wendland C², quadratura 2×2 por subcélula e **quatro subcélulas por retângulo nodal**. Dois métodos: M2 cartesiano e M5 gradiente/J-gradiente. Nenhum valor é atribuído ao Capítulo 4 ainda não auditado.

Erro percentual relativo **assinado** em lambda=k². rho é a correlação absoluta de curl(H) com a referência TM. Para M5 as duas colunas de massa e núcleo contam direções **numéricas na redução de posto**, não modos físicos propagantes.

| Método | Modo | Nós | Fator de suporte | Erro autor (%) | Erro nosso melhor curl (%) | Diferença (pp) | rho curl | rank massa | zeros curl |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| M2 | TM11 | 3×5 | 3.0 | +17.199 | +17.874 | +0.675 | 0.896232 | — | — |
| M2 | TM21 | 3×5 | 3.0 | +11.846 | +12.375 | +0.529 | 0.885607 | — | — |
| M5 | TM11 | 3×5 | 3.0 | +36.183 | +19.314 | -16.869 | 0.827131 | 26 | 14 |
| M5 | TM21 | 3×5 | 3.0 | +38.224 | +30.566 | -7.658 | 0.582856 | 26 | 14 |
| M2 | TM11 | 3×5 | 3.3 | +17.199 | +17.887 | +0.688 | 0.896722 | — | — |
| M2 | TM21 | 3×5 | 3.3 | +11.846 | +12.456 | +0.610 | 0.884854 | — | — |
| M5 | TM11 | 3×5 | 3.3 | +36.183 | +19.287 | -16.896 | 0.867775 | 26 | 14 |
| M5 | TM21 | 3×5 | 3.3 | +38.224 | +31.465 | -6.759 | 0.635449 | 26 | 14 |
| M2 | TM11 | 3×5 | 3.6 | +17.199 | +17.863 | +0.664 | 0.897137 | — | — |
| M2 | TM21 | 3×5 | 3.6 | +11.846 | +12.440 | +0.594 | 0.884878 | — | — |
| M5 | TM11 | 3×5 | 3.6 | +36.183 | +18.647 | -17.536 | 0.857820 | 26 | 14 |
| M5 | TM21 | 3×5 | 3.6 | +38.224 | -8.503 | -46.727 | 0.718638 | 26 | 14 |
| M2 | TM11 | 4×7 | 3.0 | +3.608 | +0.049 | -3.559 | 0.999326 | — | — |
| M2 | TM21 | 4×7 | 3.0 | +2.038 | +0.037 | -2.001 | 0.994919 | — | — |
| M5 | TM11 | 4×7 | 3.0 | +34.120 | +2.157 | -31.963 | 0.989868 | 52 | 27 |
| M5 | TM21 | 4×7 | 3.0 | +33.723 | +1.672 | -32.051 | 0.992172 | 52 | 27 |
| M2 | TM11 | 4×7 | 3.3 | +3.608 | -0.004 | -3.612 | 0.999779 | — | — |
| M2 | TM21 | 4×7 | 3.3 | +2.038 | -0.049 | -2.087 | 0.999291 | — | — |
| M5 | TM11 | 4×7 | 3.3 | +34.120 | +2.612 | -31.508 | 0.987131 | 52 | 27 |
| M5 | TM21 | 4×7 | 3.3 | +33.723 | +1.604 | -32.119 | 0.992118 | 52 | 27 |
| M2 | TM11 | 4×7 | 3.6 | +3.608 | -0.018 | -3.626 | 0.999932 | — | — |
| M2 | TM21 | 4×7 | 3.6 | +2.038 | -0.051 | -2.089 | 0.999734 | — | — |
| M5 | TM11 | 4×7 | 3.6 | +34.120 | +3.839 | -30.281 | 0.979483 | 52 | 27 |
| M5 | TM21 | 4×7 | 3.6 | +33.723 | +1.532 | -32.191 | 0.991325 | 52 | 27 |

## Resultados interpretáveis

1. **Estabilidade do M2 grosseiro:** TM11 em 3×5 varia apenas de +17,863% a +17,887% nos três fatores; os valores publicados são +17,199%. Não há evidência nessa faixa de que o suporte sozinho explique diferenças importantes nas malhas maiores.
2. **Não monotonicidade M5:** TM11 em 4×7 passa de +2,157% a +3,839%, enquanto o erro publicado pelo autor é +34,120%. Essa pequena varredura não elimina a diferença.
3. **Troca abrupta no TM21 M5 3×5:** candidato de melhor curl altera erro de **+30,566% / +31,465% / −8,503%**, enquanto a correlação rho **aumenta** (0,582856 / 0,635449 / 0,718638). A interpretação preferencial é **alerta de possível troca de ramo/modal matching**; não concluir uma bifurcação física do autovalor sem seguir autovetores por sobreposição entre parametrizações.
4. **Nulidades do M5 preservadas:** nas três parametrizações, a malha 3×5 permanece com posto de massa 26, 4 dependências e 14 zeros de curl; a malha 4×7 permanece com posto 52, 4 dependências e 27 zeros de curl. Estabilidade do posto não implica estabilidade da correspondência modal.
5. **Não foi testado** o peso MLS original da dissertação, nem base quadrática, variação de quadratura, tolerância de posto, divergência ou contorno. O estudo examina só o suporte de Wendland C² do nosso programa.

## Próximo procedimento discriminante

Recalcular os dois candidatos (próximo em lambda e melhor curl) por fator; projetar *todos* os campos em um grid físico comum; montar matriz cruzada de sobreposição com massa comum e identificar a trajetória de cada modo por correlação absoluta e subespaços degenerados. Registrar o índice original do autovetor e seus resíduos, curl e divergência. Não selecionar automaticamente um novo índice pelo máximo de rho sem deixar evidência de troca.

Arquivos nativos e log brutos na URL do artefato indicada. O script fonte é scripts/mls_support_sensitivity.py; nenhum espectro histórico publicado anteriormente foi modificado.
