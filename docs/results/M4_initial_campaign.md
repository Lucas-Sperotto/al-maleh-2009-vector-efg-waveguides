# M4 — Primeira campanha numérica de autovalores vetoriais (21/10/2026)

**Status:** resultados computados no GitHub Actions usando C++/Eigen do projeto, não valores sintéticos. [Execução CI com resultados](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39589873251). Artefato compactado no GitHub Actions: m4-tm-validation, incluindo tabelas completas e configurações individuais.

**Atenção:** uma correspondência numérica com modos TM não é certificação física de todos os autovalores; os indicadores de campo e as séries de refinamento são obrigatórios. Estas avaliações usam as bases MLS e pesos escolhidos para esta reprodução crítica, e ainda não provam coincidência exata com os parâmetros de Al Maleh (2009).

## Problema de referência

Retângulo PEC com lados a=2, b=1. Autovalor analítico TM11:

$$
\lambda_{11}=\pi^2(1/4+1)=1,25\pi^2\approx12,3370055.
$$

A tabela mede erro relativo em lambda, e sobreposição do rotacional calculado com a função seno TM11 (projeção de subespaço em caso de degenerescência). Os números aqui estão arredondados dos logs de CI; valores com maior precisão constam nos CSV do artefato.

| Caso | Erro relativo TM11 | Sobreposição TM11 | Núcleo numérico | Condição da massa |
|---|---:|---:|---:|---:|
| Refinamento 5×5, grau 1, suporte 3,3, Gauss 3 | 0,00013 | 0,99959 | 5 | 122952 |
| Refinamento 7×7, grau 1, suporte 3,3, Gauss 3 | 0,00003 | 0,99993 | 6 | 398000 |
| Refinamento 9×7, grau 1, suporte 3,3, Gauss 3 | 0,00004 | 0,99996 | 5 | 72374 |
| 7×7, suporte 2,8 | 0,00015 | 0,99961 | 5 | 87737 |
| 7×7, suporte 3,8 | 0,00010 | 0,99985 | 6 | 1326688 |
| 7×7, grau polinomial 2, suporte 3,8 | 0,00006 | 0,99923 | 9 | 302028 |
| 7×7, Gauss 2 | 0,00016 | 0,99996 | 6 | 328619 |
| 7×7, Gauss 4 | 0,00004 | 0,99989 | 5 | 407661 |

As **oito execuções finalizaram**. O erro espectral TM11 é pequeno mas *não monotônico*. A dimensão do núcleo numérico varia com suporte, ordem MLS e quadratura: isso requer investigação e não pode ser tratado como mero detalhe de arredondamento.

## Modos TM21 e TM12 — diferença entre autovalor e estrutura física

| Nuvem | Referência | Erro relativo lambda | Sobreposição | Divergência relativa | Fluxo normal relativo |
|---|---|---:|---:|---:|---:|
| 5×5 | TM21 | 0,00048894 | 0,99758116 | 0,172378 | 0,104247 |
| 5×5 | TM12 | 0,00963703 | **0,00000000** | **2,133391** | **1,593198** |
| 7×7 | TM21 | 0,00131734 | 0,97515439 | **0,915179** | **0,408939** |
| 7×7 | TM12 | 0,00226941 | 0,97467166 | **0,566112** | **0,414811** |
| 9×7 | TM21 | 0,00006409 | 0,99993272 | 0,031259 | 0,014696 |
| 9×7 | TM12 | 0,00013588 | 0,99950497 | 0,069497 | 0,039969 |

**Interpretação:** a nuvem 5×5 encontra um candidato com erro espectral TM12 inferior a 1%, mas cuja forma e divergência não são compatíveis com um bom modo físico. A nuvem 7×7 encontra TM21 com erro de lambda de 0,13%, mas apresenta divergência e normal de parede elevados. A nuvem 9×7 melhora ambas as famílias.

A classificação adequada nesta fase é **candidato confirmado por forma e indicadores exploratórios** ou **candidato não confirmado**; não afirmar que uma classificação negativa estabelece matematicamente um modo espúrio. Um modo físico insuficientemente resolvido pode falhar nos mesmos testes.

## Limitações e próximos experimentos

1. Variar o tamanho da nuvem mantendo critérios uniformes de escala do suporte; ampliar refinamento além de 9×7 com solucionador espectral esparso.
2. Usar estimativa de divergência e normal de fronteira e examinar sua convergência; a integração de parede agora é segmentada por células de fundo.
3. Investigar os subespaços do núcleo e a estabilidade de sua dimensão. O núcleo de discretização cartesiana não tem necessariamente a mesma estrutura de uma discretização compatível H(curl).
4. Validar correspondência entre modos numericamente quase degenerados usando projeções de subespaços e, em etapa futura, autovetores/funções FEM independentes.
5. Reconciliar os parâmetros da implementação com as equações originais da dissertação somente após obter e verificar o texto primário.

A documentação de metodologia está em [04_tm_validation.md](../implementation/04_tm_validation.md). Nenhum capítulo de teoria foi alterado nesta campanha.
