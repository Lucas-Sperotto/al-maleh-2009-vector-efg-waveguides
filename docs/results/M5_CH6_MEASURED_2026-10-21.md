# M5 — Primeira campanha medida com as bases originais do capítulo 6

**Fonte da execução:** [CI 39636316113](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39636316113) (commit `7c7f7d4f`): CTest aprovado, **12/12 configurações** e **24 comparações**. Valores de referência publicados foram transcritos das Tabelas 5/6 das imagens fornecidas; dados de execução provêm do código desta campanha.

Todos os números abaixo são **erros percentuais assinados em lambda=k²**. Candidato por autovalor mais próximo e por máxima correlação de curl(H) são registrados separadamente; nenhum é automaticamente certificado como o campo físico correto.

## Tabelas 5 e 6 — referência e simulação

| Modo | Nós | Gauss | Autor (%) | Nosso λ próximo (%) | Nosso melhor curl (%) | Δ melhor curl − autor (pp) |
|---|---|---|---:|---:|---:|---:|
| TM11 | 3×5 | 1×1 | -3.572 | +13.274 | +22.689 | +26.261 |
| TM11 | 3×5 | 2×2 | +36.183 | +19.287 | +19.287 | -16.896 |
| TM11 | 3×5 | 3×3 | +28.258 | +18.932 | +18.932 | -9.326 |
| TM11 | 3×5 | 4×4 | +26.008 | +18.901 | +18.901 | -7.107 |
| TM11 | 4×7 | 1×1 | -13.716 | +0.352 | +0.352 | +14.068 |
| TM11 | 4×7 | 2×2 | +34.120 | +2.612 | +2.612 | -31.508 |
| TM11 | 4×7 | 3×3 | +25.177 | +2.278 | +2.278 | -22.899 |
| TM11 | 4×7 | 4×4 | +23.504 | +2.251 | +2.251 | -21.253 |
| TM11 | 5×9 | 1×1 | -13.382 | -0.725 | -0.725 | +12.657 |
| TM11 | 5×9 | 2×2 | +34.360 | +0.420 | +0.420 | -33.940 |
| TM11 | 6×11 | 1×1 | -13.273 | -0.696 | -0.696 | +12.577 |
| TM11 | 7×13 | 1×1 | -13.257 | -0.880 | -0.880 | +12.377 |
| TM21 | 3×5 | 1×1 | -1.509 | -7.571 | +41.354 | +42.863 |
| TM21 | 3×5 | 2×2 | +38.224 | -6.363 | +31.465 | -6.759 |
| TM21 | 3×5 | 3×3 | +29.853 | -9.215 | +31.042 | +1.189 |
| TM21 | 3×5 | 4×4 | +27.490 | -9.478 | +30.981 | +3.491 |
| TM21 | 4×7 | 1×1 | -14.353 | -0.902 | -0.902 | +13.451 |
| TM21 | 4×7 | 2×2 | +33.723 | +1.604 | +1.604 | -32.119 |
| TM21 | 4×7 | 3×3 | +24.864 | +1.243 | +1.243 | -23.621 |
| TM21 | 4×7 | 4×4 | +23.227 | +1.217 | +1.217 | -22.010 |
| TM21 | 5×9 | 1×1 | -12.985 | -0.839 | -0.839 | +12.146 |
| TM21 | 5×9 | 2×2 | +33.415 | +0.481 | +0.481 | -32.934 |
| TM21 | 6×11 | 1×1 | -12.840 | -0.859 | -0.859 | +11.981 |
| TM21 | 7×13 | 1×1 | -12.831 | -0.933 | -0.933 | +11.898 |

## Posto e quase-núcleo (confronto com Tabela 4)

| Nós | Gauss | Rank massa | Massa nula | Curl próximo de zero | Autovalores materialmente negativos |
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

**Resultado estrutural:** com 3×5 nós, as 30 funções vetoriais originais geram uma matriz de massa de posto 26, com **4 direções nulas**, e a rigidez reduzida apresenta **14 autovalores próximos de zero** em todas as quadraturas testadas. O autor também apresenta numerosos autovalores zero na Tabela 4, cuja interpretação exige distinguir a nulidade de massa da nulidade física de curl.

## Comparações especialmente importantes

- TM11 3×5, Gauss 2×2: autor **+36,183%**; nosso melhor curl **+19,287%**, divergência de −16,896 pontos percentuais.
- TM21 3×5, Gauss 2×2: autor **+38,224%**; nosso melhor curl **+31,465%**, divergência de −6,759 pontos percentuais.
- TM21 3×5, Gauss 3×3: autor **+29,853%**; nosso melhor curl **+31,042%**, divergência de apenas +1,189 ponto percentual.
- TM11 4×7, Gauss 2×2: autor **+34,120%**; nosso melhor curl **+2,612%**, divergência de −31,508 pontos percentuais.
- TM11 5×9, Gauss 1×1: autor **−13,382%**; nosso melhor curl **−0,725%**, sinal coincidente e magnitude diferente.
- TM21 3×5, Gauss 1×1: o candidato espectral mais próximo dá **−7,571%** e o candidato por campo dá **+41,354%**, enquanto o autor informa **−1,509%**. São modos distintos; nenhum foi silenciosamente selecionado para forçar acordo.

A execução registra **zero autovalores materialmente negativos** com o cutoff de classificação padrão. Isso não contradiz os resultados negativos relatados pelo autor em 2×2: autovalores de magnitude inferior ao cutoff aparecem como near_zero e seu **sinal exato continua preservado** no eigenvalues.csv do pacote bruto.

**Não é reprodução quantitativa integral.** As equações vetoriais (6-1)–(6-6) conferem, mas a escolha de peso/suporte MLS e detalhes de tratamento de singularidade no trabalho original ainda precisam ser auditados. Diferenças grandes nos resultados de malhas refinadas não demonstram, por si, superioridade física do programa presente.

Ver `docs/implementation/11_rot_irrot_original_chapter6.md`, `docs/results/m5_author_ch6_tables_5_6.csv`, os testes M5 e a issue #11. Dados M2/M6 históricos permanecem inalterados.
