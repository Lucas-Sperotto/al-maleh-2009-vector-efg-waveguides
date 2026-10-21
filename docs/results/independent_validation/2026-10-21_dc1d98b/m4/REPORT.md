# M4 — Campanha espectral EFG cartesiano (diagnóstico exploratório)

**Não constitui certificação de reprodução da dissertação nem identificação definitiva de modos TM físicos.**

A geometria é um retângulo 2×1. O autovalor analítico TM11 é 1.25π² ≈ 12.3370055. Cada caso monta K/M novamente e resolve o autoproblema.

| Caso | Índice próximo TM11 | Erro relativo λ11 | Sobreposição curl TM11 | Diagnóstico |
|---|---:|---:|---:|---|
| local_9x7 | 54 | 0.000039245 | 0.999962611 | shape_supported_candidate |
| local_13x11 | 100 | 0.000003146 | 0.999993772 | shape_supported_candidate |
| local_17x15 | 153 | 0.000002997 | 0.999992794 | shape_supported_candidate |

## Correspondência dos modos TM21 e TM12 (refinamentos)

| Caso | Referência | Erro relativo lambda | Projeção de curl | Divergência | Fluxo normal |
|---|---|---:|---:|---:|---:|
| local_9x7 | TM21 | 0.00006409 | 0.99993272 | 0.031259 | 0.014696 |
| local_9x7 | TM12 | 0.00013588 | 0.99950497 | 0.069497 | 0.039969 |
| local_13x11 | TM21 | 0.00002959 | 0.99998616 | 0.012313 | 0.003987 |
| local_13x11 | TM12 | 0.00001294 | 0.99997458 | 0.022206 | 0.010385 |
| local_17x15 | TM21 | 0.00001813 | 0.99999301 | 0.013238 | 0.004166 |
| local_17x15 | TM12 | 0.00000912 | 0.99999494 | 0.012726 | 0.004871 |

## Diagnóstico do núcleo e massa

| Caso | Estado | Núcleo numérico | Condição da massa |
|---|---|---:|---:|
| local_9x7 | success | 5 | 72373.616629218581 |
| local_13x11 | success | 6 | 285083.85639685922 |
| local_17x15 | success | 9 | 443914.0966216285 |

Interpretação: a sobreposição mede a correlação normalizada entre curl(H) e sin(mπx/a) sin(nπy/b), com tratamento de clusters espectrais degenerados. Erro pequeno no autovalor **sozinho** não é evidência de um modo físico. Falhas de condicionamento ou de posto não são ocultadas; consulte runs.csv e os stderr.txt.

A convergência deve ser estudada em sequência com diagnóstico de divergência, fluxo normal, estabilidade de suporte e convergência de autofunções.

## Ambiente

- Python: 3.12.3
- Plataforma: Linux-6.18.33.2-microsoft-standard-WSL2-x86_64-with-glibc2.39
- Commit CI: não registrado
- Execução: local

Arquivos-fonte: tm_matches_all.csv, runs.csv e pasta por caso.
