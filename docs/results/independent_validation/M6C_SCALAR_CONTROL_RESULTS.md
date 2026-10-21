# M6c — Primeira campanha do controle escalar MLS/EFG

Resultados efetivamente produzidos em 21/10/2026 no
[GitHub Actions run 39613312752](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39613312752),
a partir da branch `feat/m6c-scalar-efg-control`.
A reabertura de espectros e autovetores completos pela rotina independente
`scripts/verify_scalar_control.py` passou na
[run 39613541043](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39613541043).

**Natureza do experimento:** esta é uma **terceira formulação** para
controle, escalar MLS/EFG com Dirichlet aproximado por colocalização de
`Ez=0` nos nós de contorno. Não substitui o problema vetorial EFG
original nem o FEM P1 independente NASA.

## Autovalor fundamental

Referência analítica para o guia PEC 2×1:
`lambda(TM11)=5*pi²/4=12.3370055013617`.

| Nuvem | DOFs escalares reduzidos | lambda mínimo | Erro relativo TM11 (%) | Sobreposição do campo Ez |
|---|---:|---:|---:|---:|
| 9×7 | 35 | 12,3342563890303 | 0,022283 | 0,99999639 |
| 13×11 | 99 | 12,3359404348179 | 0,008633 | 0,99999895 |
| 17×15 | 195 | 12,3361161171425 | 0,007209 | 0,99999957 |

**Importante sobre os modos abaixo de TM11:** em termos estritos,
o primeiro autovalor **em cada refinamento** está ligeiramente abaixo
da referência analítica. Esses três valores são as aproximações do
**próprio TM11**, identificadas pela forma modal e pela convergência, e
não constituem dezenas de modos adicionais de baixa frequência.
Não se deve escrever que o espectro escalar não contém nenhum valor
positivo numericamente inferior a lambda(TM11).

Ao contrário do FEM P1 H¹₀ conforme, a colocalização MLS não impõe
`Ez=0` em **todos os pontos da parede**. Portanto, a aproximação escalar
não tem garantia de cota superior de Rayleigh e os valores podem estar
ligeiramente abaixo da referência analítica.

## Precisão dos quatro modos TM simples

| Nuvem | TM11 erro % | TM21 erro % | TM31 erro % | TM12 erro % |
|---|---:|---:|---:|---:|
| 9×7 | 0,022283 | 0,028285 | 0,054200 | 0,074978 |
| 13×11 | 0,008633 | 0,008055 | 0,003244 | 0,009430 |
| 17×15 | 0,007209 | 0,005396 | 0,000991 | 0,006770 |

Todos os quatro erros espectrais tabelados diminuem de 9×7
até 17×15. A atribuição foi feita por proximidade espectral, de forma
exploratória: a tabela não valida degenerescências nem o restante do
espectro.

## Violar as condições de parede entre nós ainda é possível

| Nuvem | TM11 traço PEC normalizado | TM11 fluxo normal H normalizado |
|---|---:|---:|
| 9×7 | 0,001198 | 0,005573 |
| 13×11 | 0,000440 | 0,003185 |
| 17×15 | 0,000227 | 0,002204 |

Há indícios de melhoria da condição de contorno física para o TM11,
**não** de sua satisfação exata. A divergência nula de
`H=z×grad(Ez)` é uma identidade da própria representação, não uma
verificação independente de como o espaço vetorial EFG comporta-se.

## Comparação com M6a/M6b

O espectro vetorial EFG cartesiano não restrito exibiu **49, 94 e 144**
valores positivos abaixo da referência TM11. No controle escalar, os
primeiros valores são já aproximações do próprio TM11, com boa correlação
de campo. Este resultado **sustenta uma hipótese** de que a proliferação
das raízes pequenas vem das propriedades do espaço vetorial/curl-curl,
mas não determina isoladamente qual mecanismo discreto a causa.

Para resolver a questão de forma rigorosa, a próxima etapa deverá
examinar o gradiente/quase-núcleo ou uma discretização vetorial
compatível, e não substituir o problema de Maxwell pelo problema
escalar quando o objetivo for reproduzir o método de Al Maleh.

## Reabrir a campanha

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
/usr/bin/ctest --test-dir build --output-on-failure
python3 scripts/m6c_scalar_campaign.py --binary build/scalar_efg_control --out out/m6c
OPENBLAS_NUM_THREADS=1 python3 scripts/verify_scalar_control.py --campaign out/m6c
~~~

Saídas integrais em `out/m6c/{9x7,13x11,17x15}/`, incluindo
autovalores/autovetores completos e matrizes S/M/C/Q, mais
`out/m6c/REPORT.md` e `out/m6c/summary.csv`.

Os dados e manifestos da campanha M6a original permanecem intocados.
