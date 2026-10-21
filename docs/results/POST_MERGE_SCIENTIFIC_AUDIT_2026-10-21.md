# Auditoria científica pós-integração — Al Maleh EFG (2026-10-21)

**Escopo:** main após integração de M1–M6f; comparação de equações e tabelas extraídas das **imagens fornecidas dos Capítulos 5 e 6** da dissertação (2009); código C++/Python atual; campanhas efetivamente executadas. O Capítulo 4 **não está disponível em texto integral validado** nesta auditoria: embora links públicos para a tese existam, o acesso automático ao PDF retornou erro. Nenhum peso, raio ou solver será atribuído a Al Maleh sem fonte primária.

**Conclusão executiva:** a construção cartesiana do Capítulo 5 e as bases gradiente/rotacional do Capítulo 6 estão implementadas e passaram por verificações algébricas. **Não se reproduziu integralmente os números do autor**, e boa aproximação de autovalores não demonstra ausência de espúrios nem satisfação das condições de Maxwell.

## 1. Software auditado e proveniência

PRs já integradas, na ordem aplicável: #5 (M6a), #7 (M6b), #8 (M6c), #9 (M6d), #10 (M6e), #15 (M2 Cap. 5), #16 (M5 Cap. 6), #13 (M6f). Commit da main conferido: **b2454e992c50bc7b9b2414383ddbcb1c4b898858**. O [GitHub Actions da main](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39641176860) terminou com sucesso.

M2 original: [28 dados transcritos](../results/m2_author_ch5_tables_2_3.csv), [28 medições assinadas](../results/m2_ch5_observed_2026-10-21.csv). M5 original: [24 dados transcritos](../results/m5_author_ch6_tables_5_6.csv), [24 medições e tabela de rank](M5_CH6_MEASURED_2026-10-21.md). Os scripts desta auditoria recalculam desvios a partir dos **dados já versionados**, não de números estimados.

## 2. Identidade matemática e diferenças metodológicas

**Capítulo 5, Eqs. (5-1)–(5-10):** o autor usa funções escalares MLS para componentes independentes Hx e Hy; massa pela integral de produtos internos, rigidez curl-curl, guia PEC de 2×1 m. A formulação S h = k² T h implementada corresponde a (5-8). Registra-se aparente diferença de sinal entre a escrita (5-6) e (5-8) no material fornecido. Figura 9: quatro células de integração por quadrilátero de nós; tabelas com Gauss 1×1 a 5×5 por célula. A malha denominada 3×5 no autor corresponde a nx=5, ny=3.

**Capítulo 6, Eqs. (6-1)–(6-6):** as funções vetoriais são grad(Ni) e J grad(Ni), com J(dx,dy)=(-dy,dx); campo H é a soma de ambas as famílias, com coeficientes independentes. A identidade curl(grad Ni)=0 e curl(J grad Ni)=Delta(Ni) é usada corretamente na montagem. A escolha de redução de posto da massa pelo espectro positivo é nossa e **não foi verificada como idêntica à usada por Al Maleh**.

**Lacuna de Capítulo 4:** peso Wendland C², base polinomial linear, suporte radial igual a 3,3 vezes o maior espaçamento nodal, fatoração numérica e cortes de rank são escolhas da nossa implementação. Faltam páginas originais com funções peso, suporte, forma da MLS, número de vizinhos, quadratura, fronteira, solver e procedimento de classificação modal; nenhuma equivalência quantitativa pode ser pressuposta.

## 3. Auditoria quantitativa M2 cartesiano

As 28 linhas comparam os erros **percentuais assinados de lambda = k²** das Tabelas 2 e 3, não erros percentuais de frequência. A média de diferença absoluta entre os percentuais do autor e da implementação é medida em **pontos percentuais (pp)**:

| Grupo | N | Desvio médio com lambda mais próximo | Desvio médio com melhor correlação de curl(H) |
|---|---:|---:|---:|
| Todos M2 | 28 | 4,998 pp | **1,859 pp** |
| TM11 | 14 | 8,490 pp | 2,213 pp |
| TM21 | 14 | 1,506 pp | 1,506 pp |
| Malha 3×5 | 10 | 9,432 pp | **0,645 pp** |
| Outras malhas | 18 | 2,534 pp | 2,534 pp |

Em **5/28 linhas**, o pareamento por autovalor difere materialmente do pareamento por campo (mais de 1 pp); são cinco linhas TM11 na malha 3×5. No exemplo TM11, 3×5, Gauss 2×2: autor **+17,199%**, candidato mais próximo **+3,377%**, candidato de melhor correlação de curl(H) **+17,887%**. São candidatos diferentes: a coincidência de 17,887% com a ordem de grandeza publicada **não prova** que o campo é o mesmo autovetor do autor.

Na malha 4×7, Gauss 2×2, TM11: autor **+3,608%**, nosso melhor candidato **−0,004%**. Na 6×11, Gauss 2×2, TM21: autor **+0,539%**, nosso resultado **+0,004%**. Isso demonstra divergência de reprodução, **não superioridade física comprovada**.

Com pareamento por curl(H), **11 das 28 comparações** diferem menos de 1 pp da fonte; as restantes 17 ainda divergem 1 pp ou mais. O comportamento quantitativo precisa ser auditado por peso, suporte e seleção de modos, evitando escolher artificialmente valores próximos aos publicados.

## 4. Auditoria quantitativa M5 rotacional/irrotacional

As 24 linhas das Tabelas 5 e 6 foram reproduzidas com duas regras de seleção sem alterar os autoproblemas:

| Grupo | N | Desvio médio com lambda mais próximo | Desvio médio com melhor curl(H) |
|---|---:|---:|---:|
| Todos M5 | 24 | 20,762 pp | **18,139 pp** |
| TM11 | 12 | 17,621 pp | 18,406 pp |
| TM21 | 12 | 23,904 pp | 17,872 pp |
| Malha 3×5 | 8 | 22,107 pp | 14,237 pp |
| Outras malhas | 16 | 20,090 pp | 20,090 pp |

**5/24** linhas mudam materialmente entre os critérios de seleção. **Nenhuma das 24 comparações** fica a menos de 1 pp dos erros relatados pelo autor usando o candidato de melhor curl(H).

Exemplos:
- TM11 3×5 Gauss 2×2: autor +36,183%; nossa medição +19,287% (diferença −16,896 pp).
- TM21 3×5 Gauss 2×2: autor +38,224%; nossa medição +31,465% (diferença −6,759 pp).
- TM21 3×5 Gauss 3×3: autor +29,853%; nossa medição +31,042% (diferença +1,189 pp).
- TM11 4×7 Gauss 2×2: autor +34,120%; nossa medição +2,612% (diferença −31,508 pp).
- TM21 3×5 Gauss 1×1: o candidato mais próximo tem −7,571%, mas o de melhor curl(H) tem +41,354%; o autor informa −1,509%. Não escolher arbitrariamente uma raiz para simular fidelidade.

A implementação atual mostra autovalores numericamente mais próximos das referências analíticas em vários casos refinados; isso não reproduz o comportamento irregular das Tabelas 5/6, nem autoriza concluir que a decomposição resolveu a poluição espectral.

## 5. Descoberta estrutural: nulidade de massa e núcleo de curl

Medições completas em cinco nuvens e nas ordens de Gauss publicadas:

| Número de nós n | DOFs totais 2n | Posto M | Nulidade M | Autovalores próximos de zero de curl no espaço reduzido | Resto do espectro reduzido |
|---:|---:|---:|---:|---:|---:|
| 15 | 30 | 26 | 4 | 14 | 12 |
| 28 | 56 | 52 | 4 | 27 | 25 |
| 45 | 90 | 86 | 4 | 44 | 42 |
| 66 | 132 | 128 | 4 | 65 | 63 |
| 91 | 182 | 178 | 4 | 90 | 88 |

**Padrão observado nas nuvens verificadas:** rank(M)=2n−4; dim(ker M)=4; dim(ker K_reduzido)=n−1. O espaço rotacional/irrotacional contém duas dependências evidentes pela partição da unidade (somas das funções gradiente e das rotacionadas), e mais duas porque os campos constantes são representáveis por ambas as famílias quando MLS reproduz x e y. Isso demonstra **pelo menos quatro dependências**; a igualdade em quatro foi observada nessas configurações, **não provada para toda MLS**.

O conjunto de gradientes escalares tem dimensão n−1 após remover a constante e gera naturalmente campos de rotacional nulo; a contagem observada n−1 é coerente. O complemento do núcleo inclui n−3 autovalores, **sem garantia de que todos sejam modos físicos**. Na Tabela 4, o texto do autor menciona quinze zeros exatos enquanto o quadro abrevia 14 zeros: não resolver a diferença textual sem esclarecer o significado de posto, redundância e tolerância.

### Autovalores negativos

Nossa rigidez é construída de contribuições positivas ao quadrado do Laplaciano e a massa de produtos internos: ambas são formas de Gram semidefinidas positivas com quadratura de pesos positivos. Depois de reduzir a massa a um subespaço positivo, **autovalores materialmente negativos não são esperados matematicamente**. Autovalores negativos minúsculos podem decorrer de arredondamento; a referência relata negativos na quadratura 2×2. Não introduzir penalização, sinais ou filtros artificiais apenas para reproduzir tal observação.

**Risco residual de verificação:** o solucionador atual grava tanto o resíduo completo K h−lambda M h quanto sua projeção Qᵀ(K h−lambda M h), mas a CI bloqueia explicitamente apenas resíduos projetados. O resíduo completo e a sensibilidade à tolerância relativa da massa precisam ser examinados antes de tratar as raízes reduzidas como equivalentes ao problema integral original.

## 6. O que as validações independentes M6 realmente comprovam

M6a: boa convergência de quatro modos simples TM e comparação independente FEM; não certifica todo o espectro. M6b: divergência invariante máxima do dubleto TM22/TM41 fino aproximadamente **0,415325**, acima do limite exploratório 0,25; giro de autovetores não corrige essa falha. M6d: 49, 94 e 144 raízes positivas abaixo de TM11 em três nuvens cartesianas e fração gradiente média crescente — forte indício de poluição numérica, não classificação absoluta.

M6e: o Whitney FEM de arestas, ramo TM-H, não apresentou dezenas de raízes positivas muito pequenas, embora um autovalor físico estivesse ligeiramente abaixo de TM11 analítico. M6f: correlações vetoriais TM11 EFG×Whitney **0,988801 / 0,995026 / 0,997201**; TM12 **0,933912 / 0,972125 / 0,984560**; dubleto fino correlação principal mínima **0,983058**, ângulo máximo **10,562°**. A comparação ainda usa **um centroide por triângulo**, não quadratura espacial convergida. Alta correlação modal não apaga a reprovação de divergência do M6b.

O problema escalar TM de um guia homogêneo é referência matemática legítima; a discussão científica é sobre um espaço discreto **vetorial** que respeite Maxwell e generalize o tratamento a problemas mais difíceis.

## 7. Hipóteses concorrentes, prioridade e testes discriminantes

| Priorização | Questão | Ensaio verificável |
|---|---|---|
| P0 | Peso/suporte MLS da tese versus perfil atual | Obter Cap. 4 e transcrever literalmente pesos, raio, parâmetros e páginas; implementar perfil distinto sem trocar silenciosamente o baseline |
| P0 | Pareamento de autovetores físicos | Comparar campos H e curl(H), ângulos invariantes de subespaços degenerados e condições PEC, nunca só autovalor próximo |
| P0 | Massas singulares e projeção M5 | Varrer tolerâncias de posto 1e−8, 1e−10 e 1e−12; comparar rank, lambda, resíduos completos e projetados |
| P1 | Influência da integração | Variar Gauss e subdivisões separadamente, verificar positividade de pesos e definitude das formas |
| P1 | Espúrios e inconsistência de fronteira | Medir div(H), fluxo normal H·n, resíduos de Maxwell em quadratura comum |
| P2 | M6f definitivo | Reconstruir campo Whitney e campo EFG nos mesmos pontos de Gauss, refinando grade física e comparando ângulos principais |

**Novo experimento em desenvolvimento nesta PR:** scripts/mls_support_sensitivity.py muda **somente** o suporte (fatores 3,0 / 3,3 / 3,6) em 3×5 e 4×7, mantém MLS linear, Wendland C², Gauss 2×2 e quatro subcélulas. Executa as duas formulações e registra seleção espectral, correlação de curl, massa e quase-núcleo, sem alterar arquivos históricos. A primeira execução CI define seus resultados; **não são parâmetros confirmados do autor**.

## 8. Comandos de reprodução

Depois de atualizar a branch de auditoria e compilar:

    /usr/bin/ctest --test-dir build --output-on-failure
    python3 scripts/audit_m2_m5_reproduction.py --out out/m2_m5_audit.json
    PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 python3 scripts/mls_support_sensitivity.py --cartesian-binary build/tm_validation --rot-irrot-binary build/rot_irrot_rect --out out/support_sensitivity

O primeiro script só requer Python padrão e referências versionadas. O segundo realiza 12 execuções de solver e publica REPORT.md, support_sensitivity.csv e dados brutos por caso; não sobrepõe M2–M6 arquivados.

## 9. Decisão científica

Aprovar ferramentas de auditoria após CI verde. **Não encerrar** as issues de fidelidade M2/M5 e de validação integral M6f com base apenas em compilação e testes algébricos. A próxima intervenção que pode explicar as discrepâncias quantitativas é conferir o **Capítulo 4** na fonte primária, formalizar um perfil fiel ao autor e comparar fator por fator contra o baseline. Preservar explicitamente resultados desfavoráveis, incluindo modos espúrios e falhas de divergência.
