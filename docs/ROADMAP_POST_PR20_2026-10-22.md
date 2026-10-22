# Próximos passos científicos após a PR #20 — 22/10/2026

**Repositório:** \`Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides\`  
**Estado-base:** PR [#20](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/pull/20) **integrada à main** no commit \`0a728d4d293c7c3b7c8687c45e1d1f6e1421bf62\`. As duas CIs finais da PR no head \`4b3e6b6\` foram aprovadas: [pull_request](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39820787944) e [push](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39820781132). **A CI específica do commit de merge na main deve ser conferida separadamente** — sucesso do head não é substituto automático de sucesso da main.

## I. O que está efetivamente demonstrado

- **MLS / Capítulo 4:** o peso **spline quártica**, Eq. (4-2), foi verificado contra a cópia integral da dissertação. Grau polinomial **linear** e fator radial \`r/h = 1,5\` foram **inferidos** da Tabela 4 e conferidos sem reajustar contra outras tabelas; não são parâmetros explicitamente prescritos pelo autor. O baseline histórico Wendland C²/3,3h continua uma variante distinta.
- **M2 / Tabelas 2–3:** 28/28 desvios assinados \`|Δ| < 0,001 pp\`; **20/28** ficam dentro de meia unidade da última casa impressa (±0,0005 pp). Não declarar igualdade exata em todas as linhas.
- **M5 / Tabela 4:** **45/45** autovalores cartesianos não nulos dentro do arredondamento individual; **12/15** valores da representação gradiente/rotacional dentro dele. Reproduzida apenas aproximadamente para M5.
- **M5 / Tabelas 5–6:** **21/24** aproximações por melhor correlação \`curl(H)\` com \`|Δ| < 0,008 pp\`; apenas **3/24** dentro de ±0,0005 pp desse critério. Falhas F6/F7/F8 continuam abertas.
- **M5 / singularidade:** o solver tem certificados de \`ker M ⊆ ker K\` (via \`KZ\`), resíduos completos e independência do subespaço com rank correto. **Esses certificados não demonstram completude quando um cutoff elimina direções de massa positivas**. Rank e tolerância devem ser reportados.
- **M5 / refinamento:** no perfil quártico linear inferido, houve **estagnação do erro no conjunto estudado**; com MLS de grau 2, o erro decresceu nos refinamentos testados. Isso não prova uma taxa assintótica ou uma formulação universalmente livre de espúrios.
- **M7:** biblioteca para sobreposição em campo \`H\`, associação Hungarian e ângulos principais já existe e foi testada. **O driver rastreia suporte com nuvem, quadratura e rank fixos**, mas ainda *não* rastreia automaticamente o subespaço inteiro em degenerescência; marca clusters ambíguos.
- **M6:** o dubleto TM22/TM41 permanece sem certificação física. A explicação por hibridização com espúrios é uma **hipótese**, não causa demonstrada. Misturar vetores com \`λ\` distintos pode melhorar um subespaço geométrico, mas não cria autopares exatos. A comparação M6f por **um centroide de cada triângulo** ainda não é quadratura espacial convergida.

Referências primárias versionadas: [auditoria teoria × código](results/INDEPENDENT_THEORY_CODE_AUDIT.md), [revisão final da PR #20](results/PR20_FINAL_REVIEW.md), [campanha independente arquivada](results/independent_audit_2026-10-22/REPORT.md), [manual M7](implementation/12_mode_tracking_and_audit_tools.md). A auditoria prévia [21/10](results/POST_MERGE_SCIENTIFIC_AUDIT_2026-10-21.md) é **histórica**, não deve ser confundida com os resultados sob o novo peso.

## II. Ordem proposta de investigação

| Prioridade | Etapa | Questão falsificável | Entrega verificável | Issues |
|---|---|---|---|---|
| **P0** | **M8a — F6/F7/F8** | Qual hipótese explica as três exceções da reprodução das Tabelas 5/6, sem escolher arbitrariamente raízes? | Por caso, varredura modal completa, mapa de erros assinados, índices, campo, norma, contorno e residual; comparação explícita com a fonte sem substituição de dados | [#11](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/11), [#17](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/17) |
| **P0** | **M8b — rastreamento modal extensivo** | As trocas de candidatos persistem quando o *mesmo* ramo é seguido ao variar quadratura, cutoff de massa e malha? | Driver em geometria/quadratura comum; scores ponderados H/curl, seleção um-a-um, ambiguidades e clusters, log de correspondência | [#18](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/18) |
| **P1** | **M8c — M6f integral** | Correlações e ângulos principais EFG × Whitney convergem com a regra de integração espacial? E a violação de divergência/PEC? | Campos reconstruídos no mesmo Gauss físico 2×2/3×3/5×5, múltiplos refinamentos, métricas invariantes e comparação com os centroides | [#12](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/12), [#6](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/6) |
| **P1** | **M8d — consistência de ordem 2 e convergência M5** | A melhora com base quadrática persiste com quadratura efetivamente convergida, cortes de posto e malhas não idênticas? | Curvas \`h\` × \`λ\`, defeito de \`Δ\`, condicionamento, rank, \`KZ\`, resíduos completos, campos e índices físicos | [#11](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/11), [#17](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/17) |
| **P2** | **M9 — formulação vetorial compatível** | Um subespaço discreto compatível com a sequência de de Rham reduz modos espúrios e satisfaz PEC/divergência melhor que o cartesiano atual? | Proposta matemática e experimento comparativo *independente*, sem retirar zeros artificiosamente nem confundir com o escalar TM | [#6](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/6) |

### M8a — três casos concretos

- **F6**: Tabela 6, TM21, 3×5, Gauss 1×1: valor publicado **−1,509%**; candidatos atuais de erro **+1,508603%** (maior curl) e **+0,948004%** (mais próximo). Coincidência aproximada de módulo é compatível com erro de sinal, **não prova erro tipográfico** na dissertação. Investigar modo, convenção, campo e possibilidade de outra seleção.
- **F7**: Tabela 5, TM11, 3×5, Gauss 1×1: publicado **−3,572%**, candidato mais próximo **−3,571510%** (índice 18), mas candidato de melhor curl **−25,446678%** (índice 15). Esta linha demonstra que o autor pode ter usado outro critério de seleção; não apagar nenhum dos candidatos.
- **F8**: Tabela 5, TM11, 4×7, Gauss 1×1: **diferença de −0,227 pp** para o melhor candidato documentado. Causa não explicada; testar discretização, suporte e quadratura *isoladamente*, mantendo a fonte imutável.

**Critérios de aceite M8a:** reproduzir deterministamente cada caso com perfil spline quártica/grau linear/\`r/h=1,5\`, comparar contra os dados efetivamente publicados, reportar todo o espectro (incluindo valores pequenos ou negativos), índices e formas de campo; confirmar por testes sintéticos que a escolha modal não é manipulada em função do número de referência; elaborar tabela “suportado / refutado / indeterminado” para cada explicação; testes, SHA-256 e CI.

### M8b — critérios de identidade

Usar produtos internos físicos de \`H\` e \`curl(H)\` em quadratura comum, com massas e normalizações claramente documentadas; guardar **índices de coluna**, \`λ\`, sobreposições, margem de confiança e casos ambíguos. Para multiplicidades, rastrear o **subespaço inteiro** e comparar ângulos principais invariantes a rotações internas. Variar suporte, Gauss, tolerância de rank e refinamento separadamente. Nunca usar overlap alto para alegar que o vetor é físico sem também verificar divergência, fluxo normal e residual.

**Aceite:** testes sintéticos com troca de sinal/escala, cruzamento de eigenvalues, rotação interna de dubleto e casos com número diferente de colunas; relatório de trajetória com ambiguidades *sem seleção forçada*; matriz completa de candidatos e manifestos.

### M8c — referência independente e contorno

Fixar o repositório FEM Whitney no SHA \`dc1d9545c1238554c2292aa0f64a3c1b883d4118\` até decisão documentada de atualização. Reconstruir os campos de Whitney *dentro dos elementos* em pontos de Gauss físicos, respeitando orientações de arestas. Transformar corretamente geometria 1×0,5 → 2×1, \`λ_edge/4\`, normalizações e pesos de área; evitar comparar vetores normalizados apenas por pico como se fossem soluções em norma de massa. Medir \`div H\`, \`H·n\` e o pior vetor do dubleto, além de correlação e ângulo principal. **Não concluir que os espúrios desapareceram só porque campos selecionados ficaram próximos.**

### M8d — quadratura e consistência

Distinguir “quadratura mais fina entre opções avaliadas” de “quadratura convergida”: testar no mínimo mais de uma regra e registrar diferença entre integrais/métricas sucessivas, não declarar convergência a partir de uma só regra. Para MLS linear \`r/h=1,5\` o defeito de \`Δ(x²+y²)\` e a estagnação do autovalor TM11 são **observados nas famílias de malhas do relatório**, não resultados universais. O grau quadrático reproduz esse Laplaciano polinomial, mas pode demandar suportes maiores e sofrer condicionamento; publicar sempre métricas de \`cond(A)\`.

## III. Protocolo de execução e governança

1. **Uma issue de pesquisa e uma branch por hipótese separável.** Não misturar correção de solver com troca experimental de peso/raio. Preservar sempre o baseline Wendland C² e o perfil quártico *inferido* como experimentos distintos.
2. **CI e ciência separadas:** 33 CTests da PR #20 verificam código e regressões, mas não são certificação de Maxwell nem das Tabelas publicadas. Toda campanha científica nova deve apresentar \`command_line\`, \`git_sha\`, versões, perfil, dados brutos, tolerâncias e manifesto SHA-256.
3. **Arquivos históricos imutáveis:** dados anteriores permanecem na história do Git; não reescrever CSVs ou relatórios de 21/10 por resultados de 22/10. Relatórios compactados podem omitir matrizes grandes, mas devem especificar como regenerá-las e verificar seus hashes.
4. **Higiene Git:** a PR #20 está integrada; branches remotas antigas podem ser avaliadas para remoção **depois de revisão manual**. Não apagá-las automaticamente. Artefatos de Actions com retenção limitada não devem ser a única fonte de dados irrecuperáveis.
5. **Encerramento científico:** uma issue permanece aberta até seus critérios físicos/quantitativos serem atingidos ou até hipótese ser refutada por experimento rastreável. Um merge sem falha na CI não é conclusão científica.

## IV. Comandos de ponto de partida (MyHoneyBee)

\`\`\`bash
cd ~/al-maleh-2009-vector-efg-waveguides
git status --short              # preserve mls_probe.csv or other local outputs
git fetch origin
git switch main
git pull --ff-only origin main
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
/usr/bin/ctest --test-dir build --output-on-failure

# Read the consolidated evidence without rewriting experimental archives:
less docs/results/PR20_FINAL_REVIEW.md
less docs/results/INDEPENDENT_THEORY_CODE_AUDIT.md
\`\`\`

Para o próximo desenvolvimento, criar \`feat/m8a-f6-f7-f8-modal-investigation\` desde a \`main\`, com testes dos três casos e relatório em um **novo** diretório de \`docs/results/\`, sem manipular os arquivos históricos.
