# M2 — Formulação vetorial cartesiana EFG, conferida com Al Maleh (2009), capítulo 5

**Estado:** montagem K/M implementada e testada; as equações (5-1)–(5-11) e resultados do capítulo 5 foram finalmente confrontados com as páginas fornecidas pelo usuário em 21/10/2026. A reprodução exata dos erros percentuais das Tabelas 2 e 3 **não está demonstrada**. A implementação da representação rotacional/irrotacional é outra etapa (M5) e será conferida contra as páginas do capítulo 6.

## 1. Equações do autor: o que agora está realmente conferido

**Capítulo 5, “Element-Free Galerkin for Vector Fields: Cartesian Shape Functions”**:

| Fonte primária | Formulação conforme as páginas do autor | Implementação |
|---|---|---|
| (5-1a), (5-2a)–(5-2d) | H(s)=sum_i (Hxi xhat + Hyi yhat + Hzi zhat) Ni; em 2D, Ni*xhat, Ni*yhat por nó | MLS escalar para cada componente transversal Hx, Hy |
| (5-3) | condição PEC natural: (curl H) × n = 0 | condição natural da forma fraca; sem anulação artificial dos Hx/Hy dos nós na parede |
| (5-4), (5-5) | Maxwell harmônico com meio homogêneo e livre de fontes | hipótese do problema de autovalores |
| (5-6)–(5-8) | problema de Helmholtz vetorial e forma fraca; [S−k²T]h=0 | diagonalização K h = λ M h, com λ=k² |
| (5-9) | T_ij=∫Ω Ni·Nj dΩ | massa M |
| (5-10) | S_ij=∫Ω curl(Ni)·curl(Nj) dΩ | rigidez curl-curl K |
| (5-11) | k=√[(mπ/a)²+(nπ/b)²] | referências TM exatas para a=2 m, b=1 m |
| Fig. 9 | grade **3×5 nós** com **quatro células de integração por retângulo entre quatro nós vizinhos** | novo parâmetro explícito subdivision=2, mantém subdivision=1 como padrão histórico |
| Tabelas 2–3 | grades 3×5, 4×7, 5×9, 6×11 e quadraturas 1×1 até 5×5 | campanha M2-Ch5 de 14 combinações e 28 valores de referência originais |
| Figs. 10–13 | vetores de H transversal nos modos TM11 e TM21 | amostragem original Hx,Hy e desenho SVG de quiver, sem reutilizar figuras da dissertação |
| Fig. 15 | modo espúrio no caso 6×11, quadratura 2×2, com k²≈12,1456 | diagnóstico explícito; proximidade numérica de λ(TM11) não basta para classificação física |

A base vetorial do autor é enumerada de modo intercalado (N1 xhat, N1 yhat, N2 xhat, N2 yhat, ...). **Nossa ordenação de coeficientes difere apenas por permutação:** primeiro todos Hx_i, depois todos Hy_i. Essa diferença não altera os autovalores se as matrizes forem permutadas consistentemente.

### Observação crítica sobre sinais das equações originais

Na reprodução visual disponibilizada, a Eq. (5-6) aparece com termo **+k²H**, mas a forma fraca (5-7) contém **−k²(Ni·H)** e o sistema (5-8) é claramente **[S−k²T]h=0**. Isso constitui uma aparente diferença de sinal **dentro do texto fornecido**; não a corrigimos silenciosamente nem concluímos que o autor errou sem verificar a íntegra e convenções de rotacional. **Nossa implementação segue precisamente o autoproblema generalizado positivo das Eqs. (5-7)–(5-10): S h=λT h.**

## 2. Blocos cartesianos da discretização implementada

Se phi_i é uma função MLS escalar e N_i^x=(phi_i,0), N_i^y=(0,phi_i), então:

\[
\mathrm{curl}_t\,\mathbf N_i^x=-\partial_y\phi_i,\quad
\mathrm{curl}_t\,\mathbf N_i^y=+\partial_x\phi_i.
\]

\[
K_{xx}^{ij}=\int_\Omega \phi_{i,y}\phi_{j,y},\quad
K_{xy}^{ij}=-\int_\Omega\phi_{i,y}\phi_{j,x},
\quad K_{yx}=K_{xy}^\top,\quad
K_{yy}^{ij}=\int_\Omega\phi_{i,x}\phi_{j,x}.
\]

\[
M=\mathrm{diag}(M_s,M_s),\qquad
(M_s)_{ij}=\int_\Omega\phi_i\phi_j.
\]

As integrais seguem diretamente (5-9) e (5-10). O núcleo gradiente do curl-curl é inerente ao operador e **não** é apagado por penalização ou por exclusão automática de autovalores.

## 3. Geometria e células: importante para confrontar a tese

O guia do autor tem a=2 m, b=1 m. Na nomenclatura das imagens, **3×5 nós significa ny=3, nx=5**: cinco posições ao longo de x e três ao longo de y (Fig. 9). O mesmo vale para 4×7 => (nx=7, ny=4), 5×9 => (nx=9, ny=5), 6×11 => (nx=11, ny=6).

A Figura 9 descreve quatro células de integração entre quatro nós adjacentes, ou subdivisão 2 em cada eixo. Nosso protótipo inicial integrava em apenas uma célula por intervalo. Agora ambos os esquemas são possíveis, sem alterar o padrão das campanhas antigas:

- rotina original: BackgroundQuadrature(rect, nx−1, ny−1, gauss);
- configuração do capítulo 5: BackgroundQuadrature(rect, 2(nx−1), 2(ny−1), gauss).

Em cada célula, a quadratura Gauss-Legendre tensorial recebe de 1 até 5 pontos por eixo, ou seja, **1, 4, 9, 16 ou 25 pontos por célula**. Testamos áreas, momentos, quadratura de grau alto e as contagens de pontos explicitamente. A malha de integração não é a base EFG e não cria funções FEM.

O fator de suporte MLS **3,3 e o peso Wendland C² continuam escolhas deste nosso protótipo**; as imagens não documentam os detalhes completos dos pesos, estratégias de suporte, tolerâncias e solver da tese. Logo, não podemos exigir igualdade percentual rígida com os resultados do autor. O script deixa as divergências visíveis.

## 4. Referências numéricas transcritas das páginas

Pela Eq. (5-11), o autor escreve:

- TM11: k²≈12,3370 (rad/m)² — valor com maior precisão em nosso código 12,337005501;
- TM21: k²≈19,7392 (rad/m)² — 19,739208802.

**Tabela 2: erro percentual com sinal publicado para TM11**

| Nós autor | Gauss 1×1 | 2×2 | 3×3 | 4×4 | 5×5 |
|---|---:|---:|---:|---:|---:|
| 3×5 | +23,767 | +17,199 | +17,210 | +17,206 | +17,206 |
| 4×7 | +5,978 | +3,608 | +3,763 | +3,759 | +3,756 |
| 5×9 | +2,805 | +1,617 | — | — | — |
| 6×11 | +1,567 | +0,839 | — | — | — |

**Tabela 3: erro percentual com sinal publicado para TM21**

| Nós autor | Gauss 1×1 | 2×2 | 3×3 | 4×4 | 5×5 |
|---|---:|---:|---:|---:|---:|
| 3×5 | +18,557 | +11,846 | +12,103 | +12,089 | +12,085 |
| 4×7 | +4,508 | +2,038 | +2,174 | +2,164 | +2,161 |
| 5×9 | +2,313 | +1,071 | — | — | — |
| 6×11 | +1,307 | +0,539 | — | — | — |

Traços representam células **não publicadas**, não zero nem resultados estimados. Os dados de referência legíveis e sua proveniência estão em [m2_author_ch5_tables_2_3.csv](../results/m2_author_ch5_tables_2_3.csv); testes impedem preencher silenciosamente essas lacunas.

As imagens também registram Fig. 10 (H para TM11, 3×5, 5×5), Fig. 11 (TM11, 6×11, 2×2), Fig. 12 (TM21, 3×5, 5×5), Fig. 13 (TM21, 4×7, 2×2) e Fig. 15 (campo espúrio 6×11, 2×2, k²≈12,1456). **Não confundir a Figura 15 com TM11** apenas porque 12,1456 é próximo a 12,3370. Nossas saídas gráficas serão calculadas com nossos próprios autovetores, sem republicar as imagens da dissertação.

## 5. Testes e execução da nova campanha

Os testes originais cobrem reprodução de funções constantes/lineares, energia de rotação, gradiente no núcleo curl-curl, blocos K/M, simetria e massa. Novos testes cobrem quadraturas 1×1 e 5×5 e topologia de quatro células da Figura 9. O espectro é resolvido pelo módulo M3/M4 já implementado: o **M2 isolado monta K/M**, mas a comparação de autovalores precisa do solver generalizado.

No repositório:

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
/usr/bin/ctest --test-dir build --output-on-failure

# Fonte M2-Ch5: 14 combinações de grade/quadratura e 28 valores originais
python3 scripts/m2_ch5_benchmark.py \
  --binary build/tm_validation --out out/m2_ch5
~~~

O script gera: out/m2_ch5/REPORT.md, comparison.csv, metadata.json, 14 pastas com espectros completos e matrizes exportadas, campos Hx/Hy para TM11 e TM21 e **SVGs originais** correspondentes. Cada linha compara tanto o candidato com lambda mais próximo quanto o candidato que maximiza a sobreposição de curl(H) com o modo analítico. São critérios diferentes para evitar ocultar candidatos espúrios.

**O resultado de uma campanha concluída deve ser citado como valor experimental, não confundido com os percentuais publicados do autor.** As fontes M6a–M6e e seus manifestos imutáveis permanecem intactos.

### Teste manual da Figura 9

~~~bash
./build/tm_validation 5 3 1 3.3 5 out/m2_ch5/manual_3x5 2
~~~

O argumento final 2 controla a subdivisão de células. Sem ele o executável preserva o comportamento histórico (subdivisão 1). O quinto argumento 5 indica 5×5 pontos de Gauss **por célula**, não 5×5 células.

## 6. Limitações

- A conclusão do autor de que a representação cartesiana introduz modos espúrios é consistente com as campanhas M6a–M6f, mas ainda não foi reproduzido quantitativamente o exemplo de k²=12,1456 da Figura 15.
- Mesmo quando o módulo de diferença de autovalores é pequeno, o campo pode estar fisicamente incorreto; usar sobreposição, divergência, contorno e subespaços invariantes, especialmente nos modos degenerados.
- O capítulo 6 da dissertação, necessário para o **M5 rotacional/irrotacional**, será documentado separadamente quando fornecido. O conteúdo fornecido agora é o **capítulo 5**, que fundamenta o M2/M3/M4 cartesiano.


## Auditoria M2 pós-merges: qualidade do pareamento

A comparação das 28 linhas transcritas das Tabelas 2/3 mostrou diferença média absoluta dos **erros assinados de lambda**, em pontos percentuais, de **4,998 pp** se o modo for escolhido pelo autovalor mais próximo, contra **1,859 pp** se escolhido pela maior correlação espacial de curl(H). Na malha 3×5 a discrepância média cai de **9,432** para **0,645 pp**. Isso comprova que a seleção espectral sozinha pode escolher um ramo modal incorreto, mas **não certifica** a identidade física por correlação de curl isoladamente.

Uma varredura controlada de suporte (3,0/3,3/3,6), com Gauss 2×2 e quatro subcélulas, encontrou para TM11 em 3×5: +17,874%, +17,887%, +17,863%; o autor publica +17,199%. O efeito do suporte **nessa faixa** é pequeno no M2 grosseiro. Os resultados originais e limitações estão no [relatório de suporte MLS](../results/MLS_SUPPORT_SENSITIVITY_MEASURED_2026-10-21.md); relatório geral em [auditoria pós-merges](../results/POST_MERGE_SCIENTIFIC_AUDIT_2026-10-21.md). **Não extrapolar a faixa testada ou atribuir Wendland C² ao autor.**
