# Auditoria independente teoria × código — campanha 2026-10-22

Campanha NOVA. Nenhum dado histórico foi regravado. Números de Al Maleh
vêm das Tabelas 2–6 da dissertação (PDF McGill, SHA-256
`ba356899e516894301a907f5ead06392a2132dd15e5c7b3811733dbae06cc506`); números do programa vêm dos CSV desta pasta.

## 1. Identificação inversa do perfil MLS do Capítulo 4 (inferência)

A Eq. (4-2) fixa a spline quártica; o valor de r NÃO consta na tese.
Varredura do fator r/h contra os 15 autovalores não nulos cartesianos e
os 5 rotacionais/irrotacionais da Tabela 4 (malha 3×5, Gauss 2–4):

| Peso | Grau | Formulação | Melhor r/h | média |log10| G2/G3/G4 | pior erro relativo |
|---|---:|---|---:|---|---:|
| quartic_spline | 1 | cartesian | 1.5000 | 0.0006/0.0009/0.0007 | 0.46% |
| quartic_spline | 1 | rot_irrot | 1.5000 | 0.0007/0.0007/0.0007 | 0.38% |
| quartic_spline | 2 | cartesian | 2.2100 | 3.7684/3.7291/3.7215 | >1000% (sem correspondência) |
| quartic_spline | 2 | rot_irrot | 3.5000 | 0.1164/0.1372/0.1278 | 67.39% |
| wendland_c2 | 1 | cartesian | 1.7400 | 0.2140/0.1635/0.1594 | 650.39% |
| wendland_c2 | 1 | rot_irrot | 3.5000 | 0.0228/0.0284/0.0266 | 19.22% |

Os valores do autor têm 3 algarismos (arredondamento ≤0,5%, |log10| ≤ 0,0022).
Melhor r/h no limite da faixa varrida indica ausência de mínimo interno, não identificação.

## 2. Tabelas 2, 3, 5 e 6 com o perfil identificado

| Formulação | Classificação | Linhas |
|---|---|---:|
| M2 | match_best_curl | 28 |
| M5 | magnitude_match_opposite_sign | 1 |
| M5 | match_best_curl | 21 |
| M5 | match_nearest_lambda_only | 1 |
| M5 | unexplained_difference | 1 |

| Form. | Modo | Malha | Gauss | Autor % | Melhor curl % | λ próximo % | Δ curl (pp) | Classe |
|---|---|---|---:|---:|---:|---:|---:|---|
| M2 | TM11 | 3×5 | 1 | +23.767 | +23.7666 | +0.2930 | -0.0004 | match_best_curl |
| M2 | TM11 | 3×5 | 2 | +17.199 | +17.1994 | -3.7810 | +0.0004 | match_best_curl |
| M2 | TM11 | 3×5 | 3 | +17.210 | +17.2100 | -3.7282 | +0.0000 | match_best_curl |
| M2 | TM11 | 3×5 | 4 | +17.206 | +17.2057 | -3.8639 | -0.0003 | match_best_curl |
| M2 | TM11 | 3×5 | 5 | +17.206 | +17.2055 | -3.8561 | -0.0005 | match_best_curl |
| M2 | TM11 | 4×7 | 1 | +5.978 | +5.9787 | -4.3823 | +0.0007 | match_best_curl |
| M2 | TM11 | 4×7 | 2 | +3.608 | +3.6078 | +3.6078 | -0.0002 | match_best_curl |
| M2 | TM11 | 4×7 | 3 | +3.763 | +3.7627 | +3.7627 | -0.0003 | match_best_curl |
| M2 | TM11 | 4×7 | 4 | +3.759 | +3.7590 | +3.7590 | +0.0000 | match_best_curl |
| M2 | TM11 | 4×7 | 5 | +3.756 | +3.7559 | +3.7559 | -0.0001 | match_best_curl |
| M2 | TM11 | 5×9 | 1 | +2.805 | +2.8057 | +2.8057 | +0.0007 | match_best_curl |
| M2 | TM11 | 5×9 | 2 | +1.617 | +1.6175 | +1.6175 | +0.0005 | match_best_curl |
| M2 | TM11 | 6×11 | 1 | +1.567 | +1.5674 | +1.5674 | +0.0004 | match_best_curl |
| M2 | TM11 | 6×11 | 2 | +0.839 | +0.8390 | +0.8390 | -0.0000 | match_best_curl |
| M2 | TM21 | 3×5 | 1 | +18.557 | +18.5576 | +18.5576 | +0.0006 | match_best_curl |
| M2 | TM21 | 3×5 | 2 | +11.846 | +11.8466 | +11.8466 | +0.0006 | match_best_curl |
| M2 | TM21 | 3×5 | 3 | +12.103 | +12.1035 | +12.1035 | +0.0005 | match_best_curl |
| M2 | TM21 | 3×5 | 4 | +12.089 | +12.0894 | +12.0894 | +0.0004 | match_best_curl |
| M2 | TM21 | 3×5 | 5 | +12.085 | +12.0853 | +12.0853 | +0.0003 | match_best_curl |
| M2 | TM21 | 4×7 | 1 | +4.508 | +4.5086 | +3.3265 | +0.0006 | match_best_curl |
| M2 | TM21 | 4×7 | 2 | +2.038 | +2.0381 | +2.0381 | +0.0001 | match_best_curl |
| M2 | TM21 | 4×7 | 3 | +2.174 | +2.1747 | +2.1747 | +0.0007 | match_best_curl |
| M2 | TM21 | 4×7 | 4 | +2.164 | +2.1645 | -2.0804 | +0.0005 | match_best_curl |
| M2 | TM21 | 4×7 | 5 | +2.161 | +2.1614 | -2.1179 | +0.0004 | match_best_curl |
| M2 | TM21 | 5×9 | 1 | +2.313 | +2.3131 | +2.3131 | +0.0001 | match_best_curl |
| M2 | TM21 | 5×9 | 2 | +1.071 | +1.0713 | +1.0713 | +0.0003 | match_best_curl |
| M2 | TM21 | 6×11 | 1 | +1.307 | +1.3072 | +1.3072 | +0.0002 | match_best_curl |
| M2 | TM21 | 6×11 | 2 | +0.539 | +0.5390 | +0.5390 | -0.0000 | match_best_curl |
| M5 | TM11 | 3×5 | 1 | -3.572 | -25.4467 | -3.5715 | -21.8747 | match_nearest_lambda_only |
| M5 | TM11 | 3×5 | 2 | +36.183 | +36.1833 | +4.8327 | +0.0003 | match_best_curl |
| M5 | TM11 | 3×5 | 3 | +28.258 | +28.2561 | +2.6744 | -0.0019 | match_best_curl |
| M5 | TM11 | 3×5 | 4 | +26.008 | +26.0106 | -0.3199 | +0.0026 | match_best_curl |
| M5 | TM11 | 4×7 | 1 | -13.716 | -13.9430 | -13.9430 | -0.2270 | unexplained_difference |
| M5 | TM11 | 4×7 | 2 | +34.120 | +34.1277 | -0.4727 | +0.0077 | match_best_curl |
| M5 | TM11 | 4×7 | 3 | +25.177 | +25.1730 | -11.5958 | -0.0040 | match_best_curl |
| M5 | TM11 | 4×7 | 4 | +23.504 | +23.5091 | -0.2872 | +0.0051 | match_best_curl |
| M5 | TM11 | 5×9 | 1 | -13.382 | -13.3825 | -13.3825 | -0.0005 | match_best_curl |
| M5 | TM11 | 5×9 | 2 | +34.360 | +34.3616 | +34.3616 | +0.0016 | match_best_curl |
| M5 | TM11 | 6×11 | 1 | -13.273 | -13.2730 | -13.2730 | +0.0000 | match_best_curl |
| M5 | TM11 | 7×13 | 1 | -13.257 | -13.2555 | -13.2555 | +0.0015 | match_best_curl |
| M5 | TM21 | 3×5 | 1 | -1.509 | +1.5086 | +0.9480 | +3.0176 | magnitude_match_opposite_sign |
| M5 | TM21 | 3×5 | 2 | +38.224 | +38.2250 | -14.8854 | +0.0010 | match_best_curl |
| M5 | TM21 | 3×5 | 3 | +29.853 | +29.8549 | -19.8399 | +0.0019 | match_best_curl |
| M5 | TM21 | 3×5 | 4 | +27.490 | +27.4911 | -21.2434 | +0.0011 | match_best_curl |
| M5 | TM21 | 4×7 | 1 | -14.353 | -14.3536 | +1.9419 | -0.0006 | match_best_curl |
| M5 | TM21 | 4×7 | 2 | +33.723 | +33.7245 | -16.1702 | +0.0015 | match_best_curl |
| M5 | TM21 | 4×7 | 3 | +24.864 | +24.8660 | -21.7669 | +0.0020 | match_best_curl |
| M5 | TM21 | 4×7 | 4 | +23.227 | +23.2282 | -22.8068 | +0.0012 | match_best_curl |
| M5 | TM21 | 5×9 | 1 | -12.985 | -12.9859 | -12.9859 | -0.0009 | match_best_curl |
| M5 | TM21 | 5×9 | 2 | +33.415 | +33.4197 | +7.3268 | +0.0047 | match_best_curl |
| M5 | TM21 | 6×11 | 1 | -12.840 | -12.8401 | -12.8401 | -0.0001 | match_best_curl |
| M5 | TM21 | 7×13 | 1 | -12.831 | -12.8298 | -12.8298 | +0.0012 | match_best_curl |

## 3. Consistência de 2ª ordem (defeito de Δ para x²+y², valor exato 4)

| Peso | Grau | r/h | Malha | cond. máx. (base local) | defeito máx. | RMS interior |
|---|---:|---:|---|---:|---:|---:|
| quartic_spline | 1 | 1.50 | 3×5 | 26.3 | 4.720e+00 | — |
| quartic_spline | 1 | 1.50 | 5×9 | 26.3 | 4.720e+00 | 1.812e+00 |
| quartic_spline | 1 | 1.50 | 9×17 | 26.3 | 4.720e+00 | 1.812e+00 |
| quartic_spline | 1 | 1.50 | 17×33 | 26.3 | 4.720e+00 | 1.812e+00 |
| quartic_spline | 1 | 2.50 | 3×5 | 26.4 | 2.090e+00 | — |
| quartic_spline | 1 | 2.50 | 5×9 | 26.4 | 2.039e+00 | — |
| quartic_spline | 1 | 2.50 | 9×17 | 26.4 | 2.039e+00 | 4.897e-01 |
| quartic_spline | 1 | 2.50 | 17×33 | 26.4 | 2.039e+00 | 4.897e-01 |
| quartic_spline | 1 | 3.30 | 3×5 | 28.7 | 2.667e+00 | — |
| quartic_spline | 1 | 3.30 | 5×9 | 28.5 | 2.447e+00 | — |
| quartic_spline | 1 | 3.30 | 9×17 | 28.5 | 2.447e+00 | 8.864e-02 |
| quartic_spline | 1 | 3.30 | 17×33 | 28.5 | 2.447e+00 | 1.007e-01 |
| wendland_c2 | 1 | 3.30 | 3×5 | 33.4 | 2.601e+00 | — |
| wendland_c2 | 1 | 3.30 | 5×9 | 33.3 | 2.624e+00 | — |
| wendland_c2 | 1 | 3.30 | 9×17 | 33.3 | 2.624e+00 | 6.218e-02 |
| wendland_c2 | 1 | 3.30 | 17×33 | 33.3 | 2.624e+00 | 7.089e-02 |
| quartic_spline | 2 | 2.50 | 3×5 | 1040.8 | 9.166e-13 | — |
| quartic_spline | 2 | 2.50 | 5×9 | 1040.8 | 3.798e-12 | — |
| quartic_spline | 2 | 2.50 | 9×17 | 1040.8 | 1.395e-11 | 2.001e-13 |
| quartic_spline | 2 | 2.50 | 17×33 | 1040.8 | 5.309e-11 | 8.961e-13 |
| quartic_spline | 2 | 3.30 | 3×5 | 1023.5 | 3.446e-13 | — |
| quartic_spline | 2 | 3.30 | 5×9 | 988.2 | 1.606e-12 | — |
| quartic_spline | 2 | 3.30 | 9×17 | 988.2 | 6.125e-12 | 1.087e-13 |
| quartic_spline | 2 | 3.30 | 17×33 | 988.2 | 2.512e-11 | 5.071e-13 |
| wendland_c2 | 2 | 3.30 | 3×5 | 1337.5 | 3.499e-13 | — |
| wendland_c2 | 2 | 3.30 | 5×9 | 1327.7 | 1.148e-12 | — |
| wendland_c2 | 2 | 3.30 | 9×17 | 1327.7 | 4.547e-12 | 1.391e-13 |
| wendland_c2 | 2 | 3.30 | 17×33 | 1327.7 | 2.220e-11 | 6.282e-13 |

## 4. Convergência M5 (TM11, melhor curl; quadratura de produção G2 e convergida 16 subcélulas×G5)

| Peso | Grau | r/h | Malha | erro G2 (%) | erro convergido (%) | ρ curl | ‖KZ‖/‖K‖ | resíduo completo |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| quartic_spline | 1 | 1.5 | 3×5 | +36.183 | +26.388 | 0.8208 | 2.0e-15 | 3.1e-16 |
| quartic_spline | 1 | 1.5 | 5×9 | +34.362 | +24.045 | 0.8950 | 1.8e-15 | 2.9e-16 |
| quartic_spline | 1 | 1.5 | 7×13 | +34.162 | +23.836 | 0.8975 | 1.9e-15 | 2.7e-16 |
| quartic_spline | 1 | 1.5 | 9×17 | +33.948 | +23.619 | 0.8989 | 1.9e-15 | 3.6e-16 |
| quartic_spline | 1 | 1.5 | 13×25 | +33.669 | +23.338 | 0.9002 | 5.5e-15 | 6.5e-16 |
| quartic_spline | 1 | 2.5 | 3×5 | +22.620 | +21.864 | 0.8786 | 1.2e-14 | 1.5e-16 |
| quartic_spline | 1 | 2.5 | 5×9 | +2.785 | +2.607 | 0.9868 | 6.2e-15 | 4.5e-16 |
| quartic_spline | 1 | 2.5 | 7×13 | +2.103 | +2.038 | 0.9899 | 5.9e-15 | 6.3e-16 |
| quartic_spline | 1 | 2.5 | 9×17 | +1.890 | +1.864 | 0.9908 | 4.9e-15 | 5.9e-16 |
| quartic_spline | 1 | 2.5 | 13×25 | +1.754 | +1.751 | 0.9914 | 6.0e-15 | 4.9e-16 |
| quartic_spline | 1 | 3.3 | 3×5 | +9.213 | +9.117 | 0.9262 | 8.6e-14 | 2.1e-16 |
| quartic_spline | 1 | 3.3 | 5×9 | +0.737 | +0.331 | 0.9982 | 2.9e-14 | 4.0e-16 |
| quartic_spline | 1 | 3.3 | 7×13 | +0.871 | +0.161 | 0.9992 | 4.3e-14 | 8.8e-16 |
| quartic_spline | 1 | 3.3 | 9×17 | +0.907 | +0.114 | 0.9994 | 3.4e-14 | 5.0e-16 |
| quartic_spline | 1 | 3.3 | 13×25 | +0.945 | +0.096 | 0.9995 | 2.5e-14 | 1.1e-15 |
| wendland_c2 | 1 | 3.3 | 3×5 | +19.287 | +18.892 | 0.8699 | 1.6e-14 | 4.9e-16 |
| wendland_c2 | 1 | 3.3 | 5×9 | +0.420 | +0.226 | 0.9988 | 1.7e-14 | 4.4e-16 |
| wendland_c2 | 1 | 3.3 | 7×13 | +0.314 | +0.127 | 0.9994 | 6.2e-14 | 5.0e-16 |
| wendland_c2 | 1 | 3.3 | 9×17 | +0.274 | +0.082 | 0.9996 | 1.1e-13 | 5.3e-16 |
| wendland_c2 | 1 | 3.3 | 13×25 | +0.250 | +0.053 | 0.9997 | 2.1e-13 | 7.3e-16 |
| quartic_spline | 2 | 2.5 | 3×5 | +23.103 | +11.141 | 0.8520 | 1.6e-15 | 2.3e-16 |
| quartic_spline | 2 | 2.5 | 5×9 | +5.580 | +5.531 | 0.9728 | 1.5e-15 | 3.8e-16 |
| quartic_spline | 2 | 2.5 | 7×13 | +2.474 | +2.622 | 0.9870 | 1.7e-15 | 2.5e-16 |
| quartic_spline | 2 | 2.5 | 9×17 | +1.386 | +1.543 | 0.9923 | 1.5e-15 | 2.2e-16 |
| quartic_spline | 2 | 2.5 | 13×25 | +0.615 | +0.724 | 0.9964 | 2.8e-15 | 2.1e-16 |
| wendland_c2 | 2 | 3.3 | 3×5 | +16.186 | +16.275 | 0.8943 | 1.9e-15 | 2.8e-16 |
| wendland_c2 | 2 | 3.3 | 5×9 | +0.541 | +0.634 | 0.9963 | 8.6e-15 | 1.6e-16 |
| wendland_c2 | 2 | 3.3 | 7×13 | +0.261 | +0.302 | 0.9985 | 1.4e-14 | 2.2e-16 |
| wendland_c2 | 2 | 3.3 | 9×17 | +0.159 | +0.186 | 0.9991 | 3.6e-14 | 1.8e-16 |
| wendland_c2 | 2 | 3.3 | 13×25 | +0.078 | +0.093 | 0.9995 | 7.9e-14 | 2.6e-16 |

## 5. Limiar de posto da massa M5

| Caso | limiar | posto | nulidade | resíduo projetado | resíduo completo | ‖KZ‖/‖K‖ |
|---|---:|---:|---:|---:|---:|---:|
| wendland_c2 1 3.3 3×5 G2 | 1e-02 | 16 | 14 | 4.0e-16 | 2.6e-02 | 9.1e-02 |
| wendland_c2 1 3.3 3×5 G2 | 1e-04 | 24 | 6 | 9.2e-16 | 5.2e-04 | 7.0e-03 |
| wendland_c2 1 3.3 3×5 G2 | 1e-06 | 26 | 4 | 5.3e-15 | 1.7e-16 | 9.2e-15 |
| wendland_c2 1 3.3 3×5 G2 | 1e-08 | 26 | 4 | 5.3e-15 | 1.7e-16 | 9.2e-15 |
| wendland_c2 1 3.3 3×5 G2 | 1e-10 | 26 | 4 | 5.3e-15 | 1.7e-16 | 9.2e-15 |
| wendland_c2 1 3.3 3×5 G2 | 1e-12 | 26 | 4 | 5.3e-15 | 1.7e-16 | 9.2e-15 |
| wendland_c2 1 3.3 3×5 G2 | 1e-14 | 26 | 4 | 5.3e-15 | 1.7e-16 | 9.2e-15 |
| wendland_c2 1 3.3 3×5 G2 | 1e-15 | 26 | 4 | 5.3e-15 | 1.7e-16 | 9.2e-15 |
| quartic_spline 1 1.5 3×5 G2 | 1e-02 | 24 | 6 | 1.9e-16 | 1.1e-02 | 3.2e-02 |
| quartic_spline 1 1.5 3×5 G2 | 1e-04 | 26 | 4 | 2.7e-16 | 3.3e-16 | 6.0e-16 |
| quartic_spline 1 1.5 3×5 G2 | 1e-06 | 26 | 4 | 2.7e-16 | 3.3e-16 | 6.0e-16 |
| quartic_spline 1 1.5 3×5 G2 | 1e-08 | 26 | 4 | 2.7e-16 | 3.3e-16 | 6.0e-16 |
| quartic_spline 1 1.5 3×5 G2 | 1e-10 | 26 | 4 | 2.7e-16 | 3.3e-16 | 6.0e-16 |
| quartic_spline 1 1.5 3×5 G2 | 1e-12 | 26 | 4 | 2.7e-16 | 3.3e-16 | 6.0e-16 |
| quartic_spline 1 1.5 3×5 G2 | 1e-14 | 26 | 4 | 2.7e-16 | 3.3e-16 | 6.0e-16 |
| quartic_spline 1 1.5 3×5 G2 | 1e-15 | 26 | 4 | 2.7e-16 | 3.3e-16 | 6.0e-16 |
| wendland_c2 1 3.3 5×9 G2 | 1e-02 | 58 | 32 | 4.4e-16 | 5.8e-03 | 7.2e-02 |
| wendland_c2 1 3.3 5×9 G2 | 1e-04 | 84 | 6 | 2.2e-15 | 7.8e-05 | 1.9e-03 |
| wendland_c2 1 3.3 5×9 G2 | 1e-06 | 86 | 4 | 5.5e-15 | 2.7e-16 | 9.3e-15 |
| wendland_c2 1 3.3 5×9 G2 | 1e-08 | 86 | 4 | 5.5e-15 | 2.7e-16 | 9.3e-15 |
| wendland_c2 1 3.3 5×9 G2 | 1e-10 | 86 | 4 | 5.5e-15 | 2.7e-16 | 9.3e-15 |
| wendland_c2 1 3.3 5×9 G2 | 1e-12 | 86 | 4 | 5.5e-15 | 2.7e-16 | 9.3e-15 |
| wendland_c2 1 3.3 5×9 G2 | 1e-14 | 86 | 4 | 5.5e-15 | 2.7e-16 | 9.3e-15 |
| wendland_c2 1 3.3 5×9 G2 | 1e-15 | 86 | 4 | 5.5e-15 | 2.7e-16 | 9.3e-15 |
| wendland_c2 2 3.3 5×9 G2 | 1e-02 | 68 | 22 | 4.5e-16 | 4.2e-03 | 4.0e-02 |
| wendland_c2 2 3.3 5×9 G2 | 1e-04 | 82 | 8 | 1.0e-15 | 6.5e-05 | 1.3e-03 |
| wendland_c2 2 3.3 5×9 G2 | 1e-06 | 84 | 6 | 2.1e-15 | 1.7e-16 | 4.4e-15 |
| wendland_c2 2 3.3 5×9 G2 | 1e-08 | 84 | 6 | 2.1e-15 | 1.7e-16 | 4.4e-15 |
| wendland_c2 2 3.3 5×9 G2 | 1e-10 | 84 | 6 | 2.1e-15 | 1.7e-16 | 4.4e-15 |
| wendland_c2 2 3.3 5×9 G2 | 1e-12 | 84 | 6 | 2.1e-15 | 1.7e-16 | 4.4e-15 |
| wendland_c2 2 3.3 5×9 G2 | 1e-14 | 84 | 6 | 2.1e-15 | 1.7e-16 | 4.4e-15 |
| wendland_c2 2 3.3 5×9 G2 | 1e-15 | 84 | 6 | 2.1e-15 | 1.7e-16 | 4.4e-15 |
| wendland_c2 1 3.3 7×13 G1 | 1e-02 | 116 | 66 | 6.1e-16 | 4.7e-03 | 7.3e-02 |
| wendland_c2 1 3.3 7×13 G1 | 1e-04 | 172 | 10 | 1.5e-15 | 2.7e-04 | 4.9e-03 |
| wendland_c2 1 3.3 7×13 G1 | 1e-06 | 178 | 4 | 1.3e-14 | 5.5e-16 | 3.0e-14 |
| wendland_c2 1 3.3 7×13 G1 | 1e-08 | 178 | 4 | 1.3e-14 | 5.5e-16 | 3.0e-14 |
| wendland_c2 1 3.3 7×13 G1 | 1e-10 | 178 | 4 | 1.3e-14 | 5.5e-16 | 3.0e-14 |
| wendland_c2 1 3.3 7×13 G1 | 1e-12 | 178 | 4 | 1.3e-14 | 5.5e-16 | 3.0e-14 |
| wendland_c2 1 3.3 7×13 G1 | 1e-14 | 178 | 4 | 1.3e-14 | 5.5e-16 | 3.0e-14 |
| wendland_c2 1 3.3 7×13 G1 | 1e-15 | 178 | 4 | 1.3e-14 | 5.5e-16 | 3.0e-14 |
| quartic_spline 1 1.5 7×13 G2 | 1e-02 | 168 | 14 | 5.1e-16 | 3.2e-03 | 9.7e-03 |
| quartic_spline 1 1.5 7×13 G2 | 1e-04 | 178 | 4 | 3.0e-15 | 3.0e-16 | 1.1e-15 |
| quartic_spline 1 1.5 7×13 G2 | 1e-06 | 178 | 4 | 3.0e-15 | 3.0e-16 | 1.1e-15 |
| quartic_spline 1 1.5 7×13 G2 | 1e-08 | 178 | 4 | 3.0e-15 | 3.0e-16 | 1.1e-15 |
| quartic_spline 1 1.5 7×13 G2 | 1e-10 | 178 | 4 | 3.0e-15 | 3.0e-16 | 1.1e-15 |
| quartic_spline 1 1.5 7×13 G2 | 1e-12 | 178 | 4 | 3.0e-15 | 3.0e-16 | 1.1e-15 |
| quartic_spline 1 1.5 7×13 G2 | 1e-14 | 178 | 4 | 3.0e-15 | 3.0e-16 | 1.1e-15 |
| quartic_spline 1 1.5 7×13 G2 | 1e-15 | 178 | 4 | 3.0e-15 | 3.0e-16 | 1.1e-15 |

## 6. Dubleto TM22/TM41 cartesiano (grandezas invariantes de subespaço)

| Nuvem | Objeto | Índices | cos. principais c/ analítico | div mín–máx | fração gradiente máx |
|---|---|---|---|---|---:|
| 9x7 | eigenvector | 69 | 0.01386, 0.00000 | 2.8967–2.8967 | 0.0106 |
| 9x7 | eigenvector | 70 | 0.99939, 0.00000 | 0.0862–0.0862 | 0.0012 |
| 9x7 | eigenvector | 71 | 0.99886, 0.00000 | 0.0901–0.0901 | 0.0014 |
| 9x7 | computed_doublet_span | 70;71 | 0.99939, 0.99886 | 0.0862–0.0901 | 0.0014 |
| 9x7 | l2_projection_of_analytic_doublet | P(TM22);P(TM41) | 0.99998, 0.99977 | 0.0159–0.0373 | 0.0000 |
| 9x7 | l2_projection_TM22 | single | 0.99998, 0.00000 | 0.0159–0.0159 | 0.0000 |
| 9x7 | l2_projection_TM41 | single | 0.99977, 0.00000 | 0.0373–0.0373 | 0.0000 |
| 9x7 | best_2d_subspace_of_window | 69;70;71 | 0.99939, 0.99895 | 0.0773–0.0862 | 0.0015 |
| 13x11 | eigenvector | 126 | 0.00000, 0.00000 | 3.8170–3.8170 | 0.9806 |
| 13x11 | eigenvector | 127 | 1.00000, 0.00000 | 0.0065–0.0065 | 0.0000 |
| 13x11 | eigenvector | 128 | 0.99974, 0.00000 | 0.0652–0.0652 | 0.0004 |
| 13x11 | computed_doublet_span | 127;128 | 1.00000, 0.99974 | 0.0065–0.0652 | 0.0004 |
| 13x11 | l2_projection_of_analytic_doublet | P(TM22);P(TM41) | 1.00000, 0.99999 | 0.0042–0.0103 | 0.0000 |
| 13x11 | l2_projection_TM22 | single | 1.00000, 0.00000 | 0.0042–0.0042 | 0.0000 |
| 13x11 | l2_projection_TM41 | single | 0.99999, 0.00000 | 0.0103–0.0103 | 0.0000 |
| 13x11 | best_2d_subspace_of_window | 126;127;128 | 1.00000, 0.99974 | 0.0065–0.0652 | 0.0004 |
| 17x15 | eigenvector | 190 | 0.00000, 0.00000 | 4.3036–4.3036 | 0.9918 |
| 17x15 | eigenvector | 191 | 0.12782, 0.00000 | 3.2273–3.2273 | 0.9645 |
| 17x15 | eigenvector | 192 | 1.00000, 0.00000 | 0.0075–0.0075 | 0.0000 |
| 17x15 | eigenvector | 193 | 0.99179, 0.00000 | 0.4153–0.4153 | 0.0160 |
| 17x15 | computed_doublet_span | 192;193 | 1.00000, 0.99179 | 0.0075–0.4153 | 0.0160 |
| 17x15 | l2_projection_of_analytic_doublet | P(TM22);P(TM41) | 1.00000, 1.00000 | 0.0022–0.0045 | 0.0000 |
| 17x15 | l2_projection_TM22 | single | 1.00000, 0.00000 | 0.0022–0.0022 | 0.0000 |
| 17x15 | l2_projection_TM41 | single | 1.00000, 0.00000 | 0.0045–0.0045 | 0.0000 |
| 17x15 | best_2d_subspace_of_window | 190;191;192;193 | 1.00000, 0.99999 | 0.0075–0.0145 | 0.0000 |

## 7. Rastreamento M7 (M5, 3×5, Gauss 2×2, Wendland, r/h 3,00→3,60)

- Fatores: 31; transições com λ<30: 210;
  ambíguas: 0; menor sobreposição L² rastreada: 0.9561.
- Trocas do candidato 'melhor curl' TM21 (não sucessor rastreado): 3.30, 3.46, 3.50.
- Erro TM21 do candidato melhor curl varia de -20.013% a +34.000%; ρ entre 0.5834 e 0.7425.

## 8. Regressão dos espectros arquivados M6a após a mudança de base MLS

| Nuvem | modos | máx |Δλ| | máx Δλ/λ (1<λ<60) | cond. MLS arquivado (base global) | atual (base local) | núcleo arq./atual |
|---|---:|---:|---:|---:|---:|---|
| 9x7 | 126 | 2.67e-10 | 7.12e-13 | 265.6 | 33.9 | 5/5 |
| 13x11 | 286 | 3.50e-09 | 1.20e-12 | 629.3 | 34.7 | 6/6 |
| 17x15 | 510 | 2.35e-08 | 9.24e-12 | 1148.9 | 35.0 | 9/9 |

Interpretação científica: ver docs/results/INDEPENDENT_THEORY_CODE_AUDIT.md.
