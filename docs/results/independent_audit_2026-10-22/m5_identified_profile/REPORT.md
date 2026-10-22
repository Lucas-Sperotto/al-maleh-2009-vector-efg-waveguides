# M5 — Al Maleh (2009) Capítulo 6: campanha vetorial rotacional/irrotacional

**Valores publicados e valores reproduzidos são dados distintos.**
Eqs. (6-1)..(6-6) implementadas: grad Ni e z×grad Ni;
o autoproblema usa as formas fracas curl-curl/massa do Capítulo 5.
A massa singular é reduzida por autovetores positivos acima de
um limiar declarado; o espectro retido inclui todos os autovalores
negativos e quase zero, sem reclassificação seletiva.

| Modo | Grade | Gauss | Autor (%) | Candidato λ próximo (%) | Melhor forma curl (%) | Diferença forma/autor (pp) |
|---|---|---|---:|---:|---:|---:|
| TM11 | 3×5 | 1×1 | -3.572 | -3.572 | -25.447 | -21.875 |
| TM11 | 3×5 | 2×2 | +36.183 | +4.833 | +36.183 | +0.000 |
| TM11 | 3×5 | 3×3 | +28.258 | +2.674 | +28.256 | -0.002 |
| TM11 | 3×5 | 4×4 | +26.008 | -0.320 | +26.011 | +0.003 |
| TM11 | 4×7 | 1×1 | -13.716 | -13.943 | -13.943 | -0.227 |
| TM11 | 4×7 | 2×2 | +34.120 | -0.473 | +34.128 | +0.008 |
| TM11 | 4×7 | 3×3 | +25.177 | -11.596 | +25.173 | -0.004 |
| TM11 | 4×7 | 4×4 | +23.504 | -0.287 | +23.509 | +0.005 |
| TM11 | 5×9 | 1×1 | -13.382 | -13.383 | -13.383 | -0.001 |
| TM11 | 5×9 | 2×2 | +34.360 | +34.362 | +34.362 | +0.002 |
| TM11 | 6×11 | 1×1 | -13.273 | -13.273 | -13.273 | +0.000 |
| TM11 | 7×13 | 1×1 | -13.257 | -13.255 | -13.255 | +0.002 |
| TM21 | 3×5 | 1×1 | -1.509 | +0.948 | +1.509 | +3.018 |
| TM21 | 3×5 | 2×2 | +38.224 | -14.885 | +38.225 | +0.001 |
| TM21 | 3×5 | 3×3 | +29.853 | -19.840 | +29.855 | +0.002 |
| TM21 | 3×5 | 4×4 | +27.490 | -21.243 | +27.491 | +0.001 |
| TM21 | 4×7 | 1×1 | -14.353 | +1.942 | -14.354 | -0.001 |
| TM21 | 4×7 | 2×2 | +33.723 | -16.170 | +33.725 | +0.002 |
| TM21 | 4×7 | 3×3 | +24.864 | -21.767 | +24.866 | +0.002 |
| TM21 | 4×7 | 4×4 | +23.227 | -22.807 | +23.228 | +0.001 |
| TM21 | 5×9 | 1×1 | -12.985 | -12.986 | -12.986 | -0.001 |
| TM21 | 5×9 | 2×2 | +33.415 | +7.327 | +33.420 | +0.005 |
| TM21 | 6×11 | 1×1 | -12.840 | -12.840 | -12.840 | -0.000 |
| TM21 | 7×13 | 1×1 | -12.831 | -12.830 | -12.830 | +0.001 |

## Rank e número de modos

| Malha | Gauss | Rank M | Nullity M | Zeros curl | Autovalores negativos |
|---|---:|---:|---:|---:|---:|
| 3×5 | 1 | 26 | 4 | 14 | 0 |
| 3×5 | 2 | 26 | 4 | 14 | 0 |
| 3×5 | 3 | 26 | 4 | 14 | 0 |
| 3×5 | 4 | 26 | 4 | 14 | 0 |
| 4×7 | 1 | 52 | 4 | 27 | 0 |
| 4×7 | 2 | 52 | 4 | 27 | 0 |
| 4×7 | 3 | 52 | 4 | 27 | 0 |
| 4×7 | 4 | 52 | 4 | 27 | 0 |
| 5×9 | 1 | 86 | 4 | 44 | 0 |
| 5×9 | 2 | 86 | 4 | 44 | 0 |
| 6×11 | 1 | 128 | 4 | 65 | 0 |
| 7×13 | 1 | 178 | 4 | 90 | 0 |

Os zeros de massa (dependências entre funções) NÃO são os mesmos
que os autovalores zero de rigidez no problema reduzido.
Erro pequeno em λ sozinho não identifica um modo físico.

## Referência da dissertação

Tabela 4 descreve muitos autovalores zero no método rotacional/irrotacional;
e resultados negativos na quadratura 2×2, reduzidos por ordens mais altas.
Para TM11, as três raízes destacadas são aproximadamente 16,8, 15,8 e
15,5 (Gauss 2×2, 3×3, 4×4) no caso 3×5, conforme Tabela 4.
As Tabelas 5 e 6 apresentam erros sem convergência sistemática;
as Figuras 18–26 exibem campos irregulares e modos espúrios.

Nenhum parâmetro foi ajustado para reproduzir o resultado negativo
da dissertação. Diferenças quantitativas exigem cotejo com o capítulo 4
(MLS original e integrais), seleção modal e sensibilidade da redução
de posto. O driver preserva no CSV bruto o espectro retido completo
e todo o conjunto de coeficientes dos dois potenciais.

As saídas tm11_{irrot,rot}_field.csv/.svg mostram separadamente
os dois termos da Eq. (6-6), nunca copiam as Figuras 22–23.
