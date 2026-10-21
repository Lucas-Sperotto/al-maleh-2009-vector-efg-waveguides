# M5 — Formulação vetorial rotacional/irrotacional de Al Maleh (2009), capítulo 6

**Fonte primária:** páginas fornecidas do Capítulo 6 da dissertação, “Irrotational and Rotational Shape Functions”, Equações (6-1)–(6-6), Figuras 16–26, Tabelas 4–6. As funções vetoriais foram agora efetivamente cotejadas com a fonte. Os detalhes de pesos e raio MLS da seção 4.1 e da redução de matriz singular original ainda precisam de confirmação. Logo, não reivindicamos coincidência quantitativa integral.

## Equações originais conferidas

O autor escreve H=grad(phi)+R (6-1), grad(phi)=sum_g grad(N_g) phi_g (6-2), R=sum_g rho_g c_g (6-3), rho_g=z×grad(N_g) (6-4), e intercala funções de cada nó (6-5). A Eq. (6-6) é a soma das duas famílias:

\[
\mathbf H_h=\sum_g [\varphi_g\nabla N_g+c_g\,\hat z\times\nabla N_g].
\]

Em coordenadas do guia, as funções vetoriais são

\[
G_i=(N_{i,x},N_{i,y}),\quad R_i=(-N_{i,y},N_{i,x}).
\]

Conferem as identidades diferenciais

\[
\operatorname{curl} G_i=0,\quad \nabla\cdot R_i=0,\quad
\operatorname{curl} R_i=\Delta N_i=N_{i,xx}+N_{i,yy}.
\]

Nossa MLS disponibiliza Hessianas consistentes; a formulação utiliza realmente as segundas derivadas, sem confundir a representação rotacional com a base cartesiana de funções N_i xhat e N_i yhat do capítulo 5.

## Formas bilineares resultantes

Usamos a mesma forma fraca do capítulo 5, S_ij=∫curl(N_i)curl(N_j), T_ij=∫N_i·N_j e autoproblema S h=lambda T h. A redação do capítulo 6 remete às Eqs. (5-10) e (5-11) como matrizes, embora a (5-11) mostrada no capítulo 5 seja a referência analítica para k; registramos a discrepância de numeração sem alterá-la silenciosamente.

Para a ordenação por blocos, primeiro coeficientes do potencial escalar e depois coeficientes rotacionais:

\[
S=\begin{bmatrix}0&0\\0&L\end{bmatrix},
\quad L_{ij}=\int_\Omega\Delta N_i\Delta N_j\,d\Omega.
\]

\[
T=\begin{bmatrix}A&B\\B^\top&A\end{bmatrix},\quad
A_{ij}=\int_\Omega\nabla N_i\cdot\nabla N_j\,d\Omega,
\]

\[
B_{ij}=\int_\Omega(-N_{i,x}N_{j,y}+N_{i,y}N_{j,x})\,d\Omega.
\]

O bloco cruzado B é antissimétrico, mas a massa inteira T é simétrica. O código computa esses blocos diretamente na quadratura das funções MLS. Os testes verificam identidades diferenciais, formas quadráticas, simetria e reprodução afim/quadrática.

## Massa singular e zeros do operador

Por partição de unidade, ∑N_i=1, logo ∑grad(N_i)=0 e ∑(z×grad N_i)=0. Há dependências lineares necessárias entre as funções vetoriais. Como a MLS reproduz x e y, também existe sobreposição entre funções de campos constantes produzidos pelas duas famílias. Portanto, a matriz de massa T é em geral **singular**.

A implementação diagonaliza T, retém autovetores de massa com autovalores positivos acima de max(eps_abs,eps_rel*mu_max), e forma Q com QᵀTQ=I. Resolve QᵀSQ u=lambda u e reconstitui h=Qu.

**Não confundir** a nulidade de massa (funções dependentes descartadas apenas na construção de uma base independente) com os autovalores **zero de S** no espaço fisicamente representado. A campanha reporta separadamente rank, mass_nullity, número de zero de curl, todos os valores negativos e todos os valores positivos. Autovalores negativos não são apagados por uma regra pós-solução. A sensibilidade ao limiar de redução deve ser investigada antes de reproduzir a contagem de zeros do autor na Tabela 4.

## Dados primários do capítulo 6

A Tabela 4 confronta o espectro cartesiano e o rotacional/irrotacional numa grade 3×5, com Gauss 2×2, 3×3 e 4×4. Muitos autovalores rotacionais tornam-se nulos; o autor informa autovalores negativos em 2×2, que desaparecem ou se aproximam de zero com quadraturas maiores. O autovalor físico TM11 destacado fica **aproximadamente 16,8, 15,8 e 15,5** nas três quadraturas, respectivamente.

**Nota sobre a contagem:** a Tabela 4 abrevia os primeiros autovalores com a indicação “14 zeros”, mas a prosa do capítulo menciona quinze zeros exatos. Não devemos resolver essa diferença por suposição; além disso, a contagem de autovalores curl-curl nulos no espaço de massa reduzido difere da nulidade da própria matriz de massa. Nosso experimento 3×5 reproduz 14 raízes curl-curl próximas de zero no espaço reduzido de rank 26, mas a correspondência linha a linha com os modos do autor ainda deve ser investigada.

Os erros percentuais **assinados** das Tabelas 5 e 6 são:

| Modo | Grade | Gauss 1×1 | 2×2 | 3×3 | 4×4 |
|---|---|---:|---:|---:|---:|
| TM11 | 3×5 | −3,572 | +36,183 | +28,258 | +26,008 |
| TM11 | 4×7 | −13,716 | +34,120 | +25,177 | +23,504 |
| TM11 | 5×9 | −13,382 | +34,360 | — | — |
| TM11 | 6×11 | −13,273 | — | — | — |
| TM11 | 7×13 | −13,257 | — | — | — |
| TM21 | 3×5 | −1,509 | +38,224 | +29,853 | +27,490 |
| TM21 | 4×7 | −14,353 | +33,723 | +24,864 | +23,227 |
| TM21 | 5×9 | −12,985 | +33,415 | — | — |
| TM21 | 6×11 | −12,840 | — | — | — |
| TM21 | 7×13 | −12,831 | — | — | — |

Os traços são valores não apresentados pelo autor; não foram preenchidos. O CSV original com 24 observações fica em docs/results/m5_author_ch6_tables_5_6.csv.

As Figuras 16–17 mostram campo radial e campo circular associados a uma função escalar MLS central. Figuras 18–21 mostram os campos TM11, 22–23 a contribuição separada irrotacional/rotacional do TM11, 24–25 os campos TM21, 26 um modo espúrio. Segundo o autor, os modos TM físicos são pouco precisos, não têm padrão de convergência e há modos espúrios remanescentes. **Essa é uma observação sobre o método testado na dissertação, não uma prova geral de impossibilidade de métodos rotacional/irrotacional.**

## Reproduzir

O executável rot_irrot_rect recebe nx,ny,Gauss,fator_de_suporte,subdivisao,saida e limiar relativo opcional de redução da massa. A grade 3×5 do autor corresponde a nx=5, ny=3; quatro células por retângulo nodal significam subdivisao=2. O padrão da massa relativo é 1e-10 e todos os parâmetros são exportados.

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
/usr/bin/ctest --test-dir build --output-on-failure
./build/rot_irrot_rect 5 3 2 3.3 2 out/m5/manual_3x5_gauss2

python3 scripts/m5_ch6_campaign.py \
  --binary build/rot_irrot_rect --out out/m5
~~~

Saídas por caso: spectrum completo eigenvalues.csv, eigenvectors.csv, matrix_elements.csv, tm_matches.csv (candidato por lambda e por correlação de curl separados), run_metadata.txt, campos originais TM11/TM21, e parcelas irrotacional/rotacional de TM11 em CSV e SVG. No diretório geral: REPORT.md, comparison.csv, spectrum_summary.csv, metadata.json e manifest_sha256.csv.

**Validação independente:** `python3 scripts/verify_m5_ch6.py --campaign out/m5` verifica 12 espectros completos, a quantidade de autovetores/DOFs, resíduos algébricos, nulidade de massa, os 14 autovalores de curl quase nulos no caso 3×5, coordenadas e campos exportados e todos os hashes SHA256. Os diagramas radiais/circulares do nó central da grade 5×9 correspondentes às Figuras 16–17 são gerados como `irrotational_shape_central.svg` e `rotational_shape_central.svg`, a partir da nossa própria MLS.

**Resultados efetivamente calculados:** [relatório das 24 comparações e posto do operador](../results/M5_CH6_MEASURED_2026-10-21.md), [artefato GitHub Actions com todos os autovetores e matrizes](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39636814699/artifacts/12360864291). O método foi executado com os parâmetros declarados, não com os parâmetros da dissertação ainda não cotejados.

O fator de suporte 3,3, peso Wendland C² e redução de posto são escolhas experimentais declaradas desta implementação. Não houve ajuste artificial dos autovalores à dissertação; campos próximos da forma analítica não ficam automaticamente certificados como físicos. Não modificar dados M6a–M6f.

**Dependência de branches:** esta etapa M5 depende dos novos pontos de Gauss e subcélulas incorporados à PR #15 do M2. O merge da PR M5 deverá ocorrer depois da PR M2.
