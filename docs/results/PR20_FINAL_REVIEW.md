# Revisão final da PR #20 — 22/10/2026

> **Atualização pós-revisão — 22/10/2026:** a [PR #20](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/pull/20) foi **integrada à `main`** no commit `0a728d4d293c7c3b7c8687c45e1d1f6e1421bf62`, após a aprovação das CIs [39820787944](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39820787944) e [39820781132](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39820781132) no head `4b3e6b6`. As frases abaixo sobre “PR aberta”, “sem merge” ou “CI futura” documentam o **estado histórico anterior** à integração e foram mantidas como trilha da revisão. O [roteiro científico pós-merge](../ROADMAP_POST_PR20_2026-10-22.md) define as pendências. Nenhuma issue de fidelidade/física foi encerrada automaticamente por esse merge.


Revisão dos commits `b74932d`, `4b113bb`, `83d53b4` e `ae8de93`, na branch
`audit/teoria-codigo-efg`. Complementa e corrige a interpretação da campanha
histórica; não altera seus CSV, manifesto ou `REPORT.md`.
Correções de código: commit `583f13e`.
A PR permanece aberta para decisão do mantenedor, sem merge.

## Falha de CI e correções

As execuções `39807691022` e `39807681821` passaram nos 33 CTests e nas
campanhas anteriores. O executável CMake `build/independent_audit` foi
passado como `--out` e `Path.iterdir()` lançou `NotADirectoryError`.
O workflow executa e publica agora `build/independent_audit_results/`,
verifica seu manifesto e publica também um caso M5 completo verificado
independentemente por SciPy e diferenças finitas.

| Arquivos | Correção e evidência |
|---|---|
| `.github/workflows/ci.yml` | Separação entre executável e diretório; verificação dos hashes; caso M5 independente na CI |
| `scripts/independent_audit_campaign.py` | Erro explícito para arquivo e diretório não vazio; aceita vazio/inexistente; precisão por última casa impressa; hashes dos executáveis e referências |
| `tests/test_independent_audit_campaign.py` | Preservação de conteúdo/permissões, três estados da saída, arredondamento e zeros finais da fonte |
| `src/mls.cpp`, `include/al_maleh/mls.hpp`, `tests/test_mls.cpp` | Spline estável `(1−q)^3(1+3q)`; escala apenas dos suportes ativos; referência QR independente e invariâncias |
| `src/mode_tracking.cpp`, `include/al_maleh/mode_tracking.hpp`, `apps/mode_tracking.cpp`, `tests/test_mode_tracking.cpp` | Normalização robusta; confiança considera disputa por destino; clusters degenerados ambíguos; índices originais e status explícitos |
| `include/al_maleh/rot_irrot.hpp`, `tests/test_rot_irrot.cpp` | Hipóteses precisas da equivalência e três lápis analíticos independentes; solver mantido |
| `scripts/verify_m5_quartic_independent.py` | Recálculo SciPy de K/M exportadas e reconstrução MLS/Laplaciano por diferenças finitas, sem importar rotinas C++ |
| `INDEPENDENT_THEORY_CODE_AUDIT.md`, este documento e `pr20_review_2026-10-22/` | Estado dos commits, conclusões qualificadas e novas evidências, preservando os arquivos experimentais anteriores |

Problemas adicionais demonstrados: a spline expandida podia produzir peso
negativo (−4,44e−16 em vez de +4,00e−18 em q=1−1e−6); um nó inativo
remoto podia tornar a matriz local artificialmente mal condicionada; overlaps
com amplitudes 1e200/1e−200 falhavam e o QR perdia direções apenas por
escala. A margem só na linha podia ignorar competição pelo mesmo destino.
Todos esses casos têm regressões. Hungarian foi comparado com enumeração
exaustiva de 25 matrizes retangulares de dimensões 1…5.

## Fonte, perfil inferido e precisão

Conferida a cópia local integral do PDF McGill, SHA-256
`ba356899e516894301a907f5ead06392a2132dd15e5c7b3811733dbae06cc506`.
A Eq. (4-2) e a Tabela 4 foram conferidas na fonte; as 52 porcentagens das
Tabelas 2/3/5/6 foram extraídas do texto e comparadas aos CSV: zero diferenças.
Valores publicados continuam separados dos simulados. A consulta web do
PDF retornou HTTP 403 nesta revisão; usou-se a cópia local de hash idêntico.

A spline quártica é da fonte. Grau linear e r/h=1,5 são identificação inversa
pela Tabela 4, pesquisando fatores em passo 0,01 e famílias específicas;
não são parâmetros documentalmente confirmados. Tabelas 2/3/5/6 foram
avaliadas com esse perfil fixo, sem ajuste adicional. Isso oferece validação
cruzada parcial, sem identificar unicamente o solver ou todos os parâmetros.

| Tabelas | Resultado da nova execução |
|---|---|
| 2–3, M2 | 28/28 com diferença <0,001 pp; máximo 0,000728142 pp; 20/28 dentro de ±0,0005 pp da precisão impressa |
| 4, cartesiano | 45/45 valores não nulos dentro da meia unidade da última casa impressa |
| 4, M5 | 12/15 dentro desse intervalo; reprodução aproximada; os três primeiros valores não nulos diferem além do arredondamento |
| 5–6, M5 | 21/24 pelo melhor curl com diferença <0,008 pp; máximo 0,007658697 pp; somente 3/24 dentro de ±0,0005 pp por esse critério |

Exceções da Tabela 4 (autor → programa): G2 5,91→5,91877675;
G3 5,62→5,59897062; G4 5,99→6,01255361. O intervalo é ±0,005
nessas três entradas; o corte uniforme de 0,5% não certificava arredondamento.
A classificação operacional `match_best_curl` usa 0,01 pp e é explicitamente
mais larga que a precisão das porcentagens impressas. Seus zeros finais
vêm dos CSV fonte: converter +34,360 em 34,36 relaxaria indevidamente o teste.

As três exceções M5 foram preservadas:

- **F7:** TM11 3×5 G1: autor −3,572%; candidato por λ −3,571510%, índice 18;
  melhor curl −25,446678%, índice 15, ρ≈0,742419. Compatível com outra seleção,
  sem comprovar o critério do autor.
- **F6:** TM21 3×5 G1: autor −1,509%; melhor curl +1,508603%, índice 21;
  nearest +0,948004%, índice 20. Coincidência de módulo apoia hipótese de
  sinal; não comprova erro tipográfico.
- **F8:** TM11 4×7 G1: diferença de 0,227 pp permanece sem explicação.

A seleção de maior correlação é feita com todos os candidatos positivos;
não equivale a certificação física nem a rastrear o mesmo ramo entre regras.
Na convergência a correlação usa a regra comum fina; na reprodução das
tabelas usa a regra da montagem. Por isso, TM21 3×5 Wendland 3,3 G2 pode
selecionar −20,013% no primeiro estudo e +31,465% no segundo, sem mudança
no espectro. Não se deve misturar essas séries como trajetória de um modo.

## Massa singular, convergência e dubleto

A equivalência exige remover o núcleo verdadeiro de M e verificar
ker M⊆ker K. KZ≈0 e resíduos completos pequenos não provam completude
quando o corte elimina massa positiva. Um teste sintético perde um modo
zero físico mesmo com ambos os certificados nulos. Outro lápis PSD com
acoplamento ao núcleo tem λ finitos {1,5;3}, enquanto a projeção dá {2;3};
o resíduo completo detecta o problema. A fórmula de Schur requer inversão
ou pseudoinversão com compatibilidade, não vale para todo pencil simétrico.

No caso novo M5 quártico 3×5 G1, SciPy independente confirmou posto 26,
nulidade 4, ‖KZ‖/‖K‖≈4,41e−16, resíduo completo máximo≈1,63e−16 e
diferença espectral≈4,41e−13 frente ao C++. A reconstrução WLS afim com
Laplaciano por diferenças de quarta ordem reconstituiu K com erro relativo
≈7,47e−10 e confirmou índices/correlações de F6/F7 em três passos.
Na campanha completa, os máximos foram 2,13e−13 para KZ relativo e
1,31e−15 para resíduo completo. Cortes excessivos mantêm resíduos projetados
pequenos enquanto os completos sobem a 1e−4…3e−2; não há aceitação cega.

O defeito do Laplaciano linear persiste em fases homólogas de malhas
uniformes e r/h fixo; foi confirmado também fora da rotina de auditoria.
Isso não é prova geral de não convergência espectral. Para o perfil inferido,
TM11 na regra fina dá erros 26,388→24,045→23,836→23,619→23,338% em
3×5,5×9,7×13,9×17,13×25, evidência de estagnação nessa família. A regra
4 subdivisões por eixo/G5 é a referência mais fina, não quadratura
certificadamente convergida: 2/G4 ainda difere cerca de 0,29 pp no caso fino.
Com Wendland linear 3,3h os erros continuam caindo; não se demonstrou patamar.

Com grau quadrático há estudo de refino efetivo: Wendland 3,3h,
TM11 0,634→0,302→0,186→0,093%; quártico 2,5h,
5,531→2,622→1,543→0,724%, de 5×9 a 13×25. Taxas efetivas de cerca de
1,7–1,9 após 5×9 apoiam melhora; não provam taxa assintótica, compactness
discreta ou ausência geral de espúrios.

No dubleto TM22/TM41, o conteúdo analítico distribuído na janela e a
proximidade do espúrio apoiam a hipótese de hibridização. A combinação de
autovetores de λ distintos pode aproximar melhor o subespaço analítico,
mas não é novo autopar exato. A análise não corrige o resíduo/divergência
dos autopares originais nem demonstra causalidade sem rastreamento adicional.

M7 preservou os índices, λ, associações e seleção do arquivo histórico:
31 fatores; 210 transições com λ<30; zero ambíguas; menor overlap 0,956145;
trocas do candidato TM21 em 3,30/3,46/3,50. O executável varia suporte com
nuvem/quadratura/posto fixos. A biblioteca mede ângulos principais; o driver
sinaliza clusters como ambíguos, sem rastrear automaticamente o subespaço
inteiro. A flag antiga de seleção registra concordância com Hungarian;
o novo campo de status registra a confiança nessa associação.

As contagens cartesianas M6a/M6d são conciliadas sem alterar dados:
49/94/144 positivos abaixo de TM11 com corte relativo 1e−9 λmax;
49/95/148 com corte absoluto 1e−9. Em 13×11 um autovalor adicional de
1,42e−7 e em 17×15 quatro de 6,40e−8…1,98e−6 entram no conjunto positivo.
SciPy sobre K/M históricos reproduziu ambas as contagens.
Dimensão exatamente cinco do núcleo cartesiano não é teorema para toda nuvem.

## Validação, arquivos e pendências

- Release principal: configuração, build `--parallel 2` e `/usr/bin/ctest`:
  **33/33**, zero falhas.
- Release limpo `build-review-clean`, no Ubuntu/WSL disponível, mesmo build
  e CTest: **33/33**, zero falhas. Não é execução do runner remoto.
- Campanha completa em `build/independent_audit_results`: **121 hashes** válidos.
- Manifestos históricos: **120 e 168 hashes** válidos; nenhum dado regravado.
- Rejeição CLI de `--out build/independent_audit`: código 2, mensagem explícita,
  sem traceback; hash e permissões do executável preservados.
- Os 17 testes Python da campanha e a verificação SciPy/FD também passaram.
- As duas execuções remotas antigas falharam pelo conflito de caminho.
  O resultado da nova CI após push será registrado na descrição/checks da
  [PR #20](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/pull/20),
  sem antecipar aprovação neste registro local.

[pr20_review_2026-10-22/](pr20_review_2026-10-22/) contém resumo, comparações
com precisão explícita, proveniência e inventário da execução local, JSON da
verificação independente e um caso M5 **completo** com K/M, coeficientes e
campos. O inventário `local_campaign_manifest.json` descreve os 121 arquivos
na pasta de build, não afirma que todos estejam nesta pasta compacta.
O manifesto raiz desta pasta inventaria somente os arquivos nela presentes.
Como a validação precedeu o commit, a proveniência local registra `ae8de93`
com árvore alterada e hashes dos executáveis; não simula uma execução em
árvore limpa. A compilação limpa refere-se ao diretório de build novo.
A campanha histórica compacta mantém espectros/diagnósticos/logs; matrizes e
campos completos dos demais casos são regeneráveis pelos drivers. O workflow
publica a campanha e o caso independente como artefatos com seus manifestos.
`mls_probe.csv` não rastreado foi preservado com SHA-256
`75984091a49f31755bb5b7d5cac4ebd855dc1ed1e98efdbf623b53427a48c684`.

Pendências cobertas pelas issues #11/#17: F6/F7/F8, parâmetros/solver da
fonte e convergência/quadratura; #12: hibridização e quadratura comum
EFG×Whitney; #18: rastreamento entre quadraturas/postos e de subespaços
degenerados. Não há novo problema crítico/alto demonstrado sem correção ou
registro. Os limites científicos acima permanecem explícitos; essas issues
não são concluídas pela PR. Recomenda-se integrar as correções e a auditoria
apenas com CI final verde e aceitação desse escopo parcial. O merge é decisão
do mantenedor e não foi executado.

Para repetir o controle independente, use uma pasta nova no comando do driver:

```bash
build/rot_irrot_rect 5 3 1 1.5 2 build/m5_quartic_review_local 1e-10 1 quartic_spline
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 python3 scripts/verify_m5_quartic_independent.py build/m5_quartic_review_local
python3 scripts/validation_manifest.py --verify --campaign docs/results/pr20_review_2026-10-22
```
