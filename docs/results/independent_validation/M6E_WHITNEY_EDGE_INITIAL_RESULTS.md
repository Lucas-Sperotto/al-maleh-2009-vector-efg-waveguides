# M6e — FEM de arestas Whitney versus EFG cartesiano: primeira campanha

Resultados medidos em 21/10/2026 pela
[execução GitHub Actions 39623279639](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39623279639),
com 20/20 CTests aprovados. O FEM independente foi compilado diretamente
do repositório Lucas-Sperotto/fem-eigen-em-2d3d-whitney no commit
dc1d9545c1238554c2292aa0f64a3c1b883d4118.

**Ramo utilizado:** TM vetor transversal H com funções Whitney de aresta;
nunca confundir com o ramo TE vetor transversal E, que tem condições de
contorno diferentes.

**Escala:** o domínio FEM nativo é 1×0,5 m e o domínio EFG é
2×1 m; para comparar usou-se lambda(edge)/4.
O programa não modificou o solver FEM nem o problema EFG original.

## Raízes abaixo da referência analítica TM11

| Nuvem EFG | Mesh edge (células) | EFG positivos abaixo de TM11 | FEM edge positivos abaixo de TM11 | FEM edge abaixo de 0,5 TM11 |
|---|---|---:|---:|---:|
| 9×7 | 8×4 | 49 | 1 | 0 |
| 13×11 | 12×6 | 94 | 1 | 0 |
| 17×15 | 16×8 | 144 | 1 | 0 |

O método Whitney apresentou uma raiz positiva ligeiramente **abaixo**
do autovalor exato de TM11 em cada refinamento; ela foi atribuída
experimentalmente ao próprio TM11, por sua proximidade espectral.
NÃO é correto dizer que Whitney não apresentou qualquer autovalor
positivo abaixo da referência analítica. O que foi observado é
**ausência de dezenas de raízes adicionais muito pequenas** para
os cortes numéricos e meshes empregados.

O cutoff numérico explícito do edge está em
out/m6e/summary.csv e nos artefatos nativos da execução.

## Comparação de autovalores simples

Erros relativos percentuais em relação aos autovalores TM analíticos
do guia PEC 2×1. Pareamento por proximidade espectral **não** certifica
correspondência vetorial de campos.

| Modo | EFG 9×7 | Whitney 8×4 | EFG 13×11 | Whitney 12×6 | EFG 17×15 | Whitney 16×8 |
|---|---:|---:|---:|---:|---:|---:|
| TM11 | 0,00392% | 0,10220% | 0,00031% | 0,04105% | 0,00030% | 0,02232% |
| TM21 | 0,00641% | 1,59507% | 0,00296% | 0,73881% | 0,00181% | 0,42126% |
| TM31 | 0,04825% | 1,20467% | 0,01320% | 0,62449% | 0,00651% | 0,36743% |
| TM12 | 0,01359% | 4,93572% | 0,00129% | 2,15203% | 0,00091% | 1,21016% |

No refinamento fino, lambda(TM11)_analítica=12,3370055,
lambda_EFG=12,3370425 e lambda_Whitney_reescalado=12,3342513.

**Interpretação:** o EFG demonstra erros espectrais baixos nos modos
simples escolhidos, mesmo com poluição de autovalores abaixo de TM11.
A discrepância FEM Whitney diminui sistematicamente entre as três
malhas para os quatro modos simples. A comparação usa graus de
liberdade, escalas e malhas diferentes, portanto não pode ser lida
como análise de eficiência por custo.

## Dupla degenerada TM22/TM41

A referência possui lambda=5π²≈49,3480220 com multiplicidade 2.

| Mesh | EFG A | EFG B | FEM Whitney A (reescalado) | FEM Whitney B (reescalado) |
|---|---:|---:|---:|---:|
| 9×7 / 8×4 | 49,3592571 | 49,4393821 | 48,7770345 | 48,9917712 |
| 13×11 / 12×6 | 49,3476601 | 49,3605466 | 49,1997120 | 49,2449322 |
| 17×15 / 16×8 | 49,3486325 | 49,3575232 | 49,2831947 | 49,2980266 |

Essas duas raízes foram identificadas por proximidade com 5π²,
não por correlação vetorial ou projeção conjunta dos campos.
**A existência de duas raízes próximas NÃO prova que o espaço
invariante corresponda aos dois modos TM corretos**.
O defeito de divergência da dupla EFG fina descoberto em M6b,
com máximo 0,415325, continua em aberto.

## Qual interpretação é fundamentada?

O fato de o Whitney edge (espaço H(curl)) não exibir dezenas de
raízes numéricas baixas nos meshes testados é consistente com
a hipótese M6d: o EFG vetorial cartesiano não representa de forma
exata e compatível o núcleo gradiente da discretização de Maxwell.
Contudo, associação não equivale a prova matemática de que somente
esse mecanismo causa a poluição.

O programa guardou os CSVs nativos completos do FEM de aresta,
incluindo S/T, todos os autovetores, logs, modos TE/TM e campos
exportados, mais o espectro EFG arquivado e hashes de proveniência.

## Auditoria independente das matrizes nativas FEM

A campanha posterior incluiu uma segunda verificação matemática,
executada no [GitHub Actions run 39623940672](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39623940672).
O script scripts/verify_m6e_edge.py reabriu TODOS os coeficientes
de autovetores TM da rotina LAPACK, reconstruindo S/T dos arquivos
CRS nativos, e verificou simultaneamente a forma generalizada
S v = lambda T v e o produto V^T T V.

| Malha FEM Whitney | Autovetores TM | Erro máximo V^T T V-I | Resíduo EVP relativo máximo |
|---|---:|---:|---:|
| 8×4 | 108 | 1,9317881e-14 | 5,9739309e-16 |
| 12×6 | 234 | 2,5091040e-14 | 1,1387772e-15 |
| 16×8 | 408 | 6,5725203e-14 | 1,1071835e-15 |

Todos os três casos passaram. Isso elimina como explicação óbvia
da diferença de espectros um problema de diagonalização ou
reconstrução dos coeficientes do FEM edge no arquivo publicado.

## Etapa seguinte

Comparar campos H em um grid de quadratura **comum** usando os
coeficientes EFG e as bases Whitney, com normalização L2, orientação
de sinal e ângulos principais para o dubleto degenerado. Também
comparar métricas de divergência e fluxo normal, e não apenas
autovalores.

Ver: docs/implementation/09_whitney_edge_cross_validation.md
