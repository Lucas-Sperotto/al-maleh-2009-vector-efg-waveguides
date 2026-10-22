# Auditoria científica e técnica independente — teoria × código (2026-10-22)

Branch `audit/teoria-codigo-efg`, criada a partir de `main` em
`c217aa5b57e92fd9c5be54db06d4458f4d2d0dbc`. A auditoria original foi
publicada nos commits `b74932d`, `4b113bb`, `83d53b4` e `ae8de93`, na
[PR #20](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/pull/20).
A revisão final e suas ressalvas estão em
[PR20_FINAL_REVIEW.md](PR20_FINAL_REVIEW.md). Campanha original preservada e verificável:
[`independent_audit_2026-10-22/`](independent_audit_2026-10-22/REPORT.md)
(manifesto SHA-256 de 120 arquivos; `python3 scripts/validation_manifest.py --verify --campaign docs/results/independent_audit_2026-10-22`).

## Convenção de proveniência usada em todo o texto

| Rótulo | Significado |
|---|---|
| **[FONTE]** | publicado por Al Maleh (2009), conferido no PDF |
| **[PROGRAMA]** | produzido pelo código deste repositório (arquivo citado) |
| **[TEOREMA]** | propriedade matemática demonstrada neste documento |
| **[INFERÊNCIA]** | conclusão deduzida dos dados, não afirmada pela fonte |
| **[HIPÓTESE]** | explicação plausível ainda não testada |
| **[DECISÃO]** | escolha metodológica nossa |

---

## 1. Escopo e commit inicial

- Commit inicial: `c217aa5` (merge da PR #19). Árvore limpa, exceto o
  arquivo não rastreado `mls_probe.csv`, preservado intacto
  (SHA-256 `75984091…c684`).
- Baseline: Release, `cmake -S . -B build`, `/usr/bin/ctest`: **28/28
  testes aprovados**.
- Ambiente: Ubuntu 24.04.5 (WSL2, MyHoneyBee), g++ 13.3.0, CMake 3.28.3,
  Eigen 3.4.0, Python 3.12.3. Com `PYTHONNOUSERSITE=1`: NumPy 1.26.4,
  SciPy 1.11.4, Matplotlib 3.6.3 (o NumPy 2.4.6 do site do usuário é
  incompatível com o SciPy do sistema; os testes do CMake já isolam isso).
- Auditados: `include/`, `src/`, `apps/`, `scripts/`, `tests/`, `docs/`,
  `CMakeLists.txt` e `.github/workflows/ci.yml`. O repositório FEM
  independente não foi modificado nem clonado; a análise do M6e/M6f usou
  apenas os dados e scripts deste repositório.

## 2. Fontes científicas efetivamente consultadas

1. **Dissertação completa** — H. Al Maleh, *A Vector Element-Free Galerkin
   Method for Waveguide Eigenvalue Problems*, M.Eng., McGill, 2009. PDF
   obtido em `mcgill.scholaris.ca` (59 páginas, SHA-256
   `ba356899e516894301a907f5ead06392a2132dd15e5c7b3811733dbae06cc506`).
   Texto extraído com `mutool`; páginas 18 (Eqs. 4-1/4-2), 28 (Eqs.
   5-6–5-10) e 39 (Tabela 4) foram renderizadas e conferidas visualmente.
   A auditoria anterior registrava o Capítulo 4 como indisponível: **esta
   lacuna foi fechada**.
2. As 52 entradas já transcritas das Tabelas 2, 3, 5 e 6 (CSV do
   repositório) foram comparadas com o texto do PDF: **0 divergências**.
3. Nova transcrição: [Tabela 4](author_ch6_table4_spectra.csv) (78
   linhas: 60 cartesianas, 18 rotacionais/irrotacionais), conferida com a
   imagem da página 39.

Nenhuma outra referência bibliográfica foi usada como fonte de números.

## 3. Equações confirmadas

| Item | Fonte | Situação no código |
|---|---|---|
| Maxwell harmônico, convenção e^{jωt}: ∇×E=−jωμ₀H (5-4), ∇×H=jωε₀E (5-5) | [FONTE] | implica ∇×∇×H = k²H, k²=ω²μ₀ε₀ |
| Condição PEC natural (∇×H)×aₙ=0 (5-3) | [FONTE] | forma fraca sem termo de contorno: `assemble_cartesian`, `assemble_rot_irrot` |
| Forma fraca a(Nᵢ,H)=(∇×Nᵢ)·(∇×H)−k²Nᵢ·H (5-7); [S−k²T]h=0 (5-8) | [FONTE] | `K h = λ M h`, λ=k² ≥ 0 |
| Tᵢⱼ=∫Nᵢ·Nⱼ (5-9), Sᵢⱼ=∫(∇×Nᵢ)·(∇×Nⱼ) (5-10) | [FONTE] | blocos cartesianos conferidos termo a termo; teste de energia independente pontual |
| k²=(mπ/a)²+(nπ/b)², a=2, b=1; TM11=12,3370, TM21=19,7392 (5-11) | [FONTE] | `rectangular_tm_reference` |
| H=∇φ+R, ρ_g=ẑ×∇N_g, Eqs. (6-1)–(6-6) | [FONTE] | Gᵢ=(Nᵢ,ₓ,Nᵢ,ᵧ), Rᵢ=(−Nᵢ,ᵧ,Nᵢ,ₓ); M_GG=M_RR=∫∇Nᵢ·∇Nⱼ, M_GR=∫(−Nᵢ,ₓNⱼ,ᵧ+Nᵢ,ᵧNⱼ,ₓ), K_RR=∫ΔNᵢΔNⱼ |
| MLS (4-3)–(4-12): A=Σw p pᵀ, N_g=pᵀA⁻¹B_g | [FONTE] | `MLSEvaluator` (derivadas 1ª/2ª por diferenciação de Aγ=p) |
| **Peso spline quártica** w(l/r)=1−6(l/r)²+8(l/r)³−3(l/r)⁴ (4-2), l=‖s−x_g‖ (4-1) | [FONTE] | **antes ausente**; agora `MLSWeight::QuarticSpline` |
| Quadratura Gauss-Legendre por célula; 4 células por retângulo nodal (Fig. 9) | [FONTE] | `BackgroundQuadrature`, `subdivision=2` |
| Base linear {1,x,y} ou quadrática (Tabela 1) | [FONTE] lista ambas, sem dizer qual usou | identificada como **linear** (§9) [INFERÊNCIA] |

Identidades diferenciais verificadas: curl Gᵢ = Nᵢ,ᵧₓ−Nᵢ,ₓᵧ = 0;
div Rᵢ = −Nᵢ,ᵧₓ+Nᵢ,ₓᵧ = 0; curl Rᵢ = ΔNᵢ (a MLS devolve uma única Hessiana
mista, logo as duas identidades são exatas por construção).

### 3.1 Problema contínuo [TEOREMA]

Para H ∈ H(curl,Ω): (curl H, curl v) = k²(H,v) ∀v. Tomando v=∇q, q∈H¹(Ω)
sem condição de contorno: k²(H,∇q)=0. Logo, para k²≠0, div H=0 em Ω e
n·H=0 em ∂Ω no sentido fraco. Com ψ = curl H, a equação dá −Δψ=k²ψ e a
condição natural é ψ=0 em ∂Ω (Ez=0). Espectro: 0 com núcleo de dimensão
infinita (∇H¹) e os valores TM (mπ/a)²+(nπ/b)², m,n ≥ 1. Simetria e
semidefinitude de K e M não excluem modos espúrios do problema discreto.

## 4. Inconsistências teóricas identificadas na fonte

| ID | Achado | Evidência | Classe |
|---|---|---|---|
| F1 | (1-1) e (5-6) escrevem ∇×∇×H **+** k²H = 0; de (5-4)/(5-5) segue ∇×∇×H **−** k²H = 0. (5-7) e (5-8) estão corretas e foram as implementadas | PDF p. 3 e p. 28 (imagem) | [FONTE] erro tipográfico |
| F2 | Cabeçalho da Tabela 2 diz k²=12,3770; o texto diz 12,3370 | Reproduzimos as 14 linhas TM11 com 5π²/4=12,33700550 a <0,001 pp; com 12,3770 o desvio seria ≈0,38 pp | [INFERÊNCIA] tipografia no cabeçalho |
| F3 | Cap. 6 cita as matrizes como "(5-10) e (5-11)"; devem ser (5-9) e (5-10) | PDF p. 36 | [FONTE] numeração |
| F4 | Texto: "Fifteen of these are exact zeros"; Tabela 4: "14 zeros" | O problema reduzido tem exatamente n−1=14 zeros (§5.2); não é possível obter 15 sem contar a parte singular do pencil | [INFERÊNCIA] inconsistência interna não resolvível |
| F5 | "In the 2×2 integration scheme, some negative eigenvalues appear" | Teorema §5.3: autovalores materialmente negativos são impossíveis no problema deflacionado com pesos de Gauss positivos | [HIPÓTESE] artefato do solver do autor (ex.: pencil singular não deflacionado) |
| F6 | Tabela 6, TM21 3×5, Gauss 1×1: −1,509% | Com o perfil identificado o espectro tem +1,5086% e +0,948% perto de TM21, nenhum valor em −1,509% | [HIPÓTESE] possível erro de sinal, não comprovado |
| F7 | Tabela 5, TM11 3×5, Gauss 1×1: −3,572% | É exatamente nosso candidato pelo λ mais próximo (−3,5715%); o candidato de maior correlação de curl (ρ=0,742) é outro autovetor (−25,447%) | [INFERÊNCIA] compatível com critério de seleção diferente |
| F8 | Tabela 5, TM11 4×7, Gauss 1×1: −13,716% | Nosso −13,943% (Δ=0,227 pp), mesmo autovetor por curl e por λ | **não explicado** |
| F9 | Cap. 4 não informa o raio r | — | lacuna da fonte; inferido em §9 |

## 5. Análise matemática da formulação rotacional/irrotacional (M5)

### 5.1 Pencil singular [TEOREMA]

Seja f = φ + iψ um polinômio holomorfo reproduzido pela MLS (grau 1:
f ∈ {1, i, z, iz}; grau 2: também z², iz²). Por Cauchy–Riemann,
∇φ + ẑ×∇ψ ≡ 0, e ψ é harmônica, logo ΔNᵢ aplicado aos coeficientes
nodais de ψ é zero. O vetor de coeficientes (φ(xᵢ), ψ(xᵢ)) está em
**ker M ∩ ker K** para qualquer nuvem e qualquer quadratura. Portanto
rank(M) ≤ 2n−4 (grau 1) ou 2n−6 (grau 2) e **det(K−λM) ≡ 0**: o problema
generalizado original é um pencil singular e precisa ser deflacionado.
Verificado em `test_rot_irrot.cpp::singular_pencil_theorem` (‖Mz‖,‖Kz‖
relativos < 1e-13, graus 1 e 2).

### 5.2 Correspondência entre problema reduzido e original [TEOREMA]

Seja Z uma base ortonormal do **núcleo verdadeiro** de M e Q uma base
M-ortonormal de seu complemento. Se KZ = 0 (ker M ⊆ ker K), então
range K ⊆ range M e todo autopar de QᵀKQ y=λy é autopar de Kh=λMh com
**resíduo completo nulo**, h=Qy, com liberdade de adicionar o núcleo comum.
Trata-se do espectro no quociente pelo núcleo comum: det(K−λM)≡0 não
define, por si só, um espectro regular. Se KZ ≠ 0 e ZᵀKZ é inversível,
o complemento de Schur é QᵀKQ − QᵀKZ(ZᵀKZ)⁻¹ZᵀKQ. Para K semidefinida
positiva, a pseudoinversa também pode ser usada com a compatibilidade do
bloco cruzado. Isso não é uma fórmula geral para todo pencil simétrico.
Quando o corte remove direções de massa positiva, essas hipóteses deixam
de valer: KZ pequeno e resíduos pequenos, sozinhos, não certificam a
completude do espectro. O resíduo projetado
Qᵀ(Kh−λMh) é ~0 **por construção** e não distingue os dois casos.

Implementado: `mass_kernel_stiffness_coupling = ‖KZ‖_F/‖K‖_F`,
`max_full_residual`, `max_reduced_residual`, `mass_eigen_max_removed`.
[PROGRAMA] Em todos os 120 casos M5 da campanha: acoplamento ≤ 2,1e-13 e
resíduo completo ≤ 1,1e-15. Com limiar 1e-2 (sobre-redução), o resíduo
projetado continua ≤ 1e-15 enquanto o completo sobe a 1e-4…3e-2
(`m5_rank_sweep.csv`) — demonstração de que o critério projetado não basta.

Dimensões: dim ker(K_red) = dim ker K − dim ker M quando ker M ⊆ ker K.
Com grau 1, dim ker K = n + d_Δ, d_Δ = dim{b : Σbᵢ ΔNᵢ(x_q)=0 ∀q} ≥ 3
(1, x, y). **Se dim ker M=4**, segue dim ker(K_red) ≥ n−1.
A igualdade n−1 requer também d_Δ=3; subintegração pode criar dependências
adicionais. Os gradientes de V_h pertencem ao núcleo de curl por construção. A igualdade rank(M)=2n−4 e
dim ker(K_red)=n−1 é **observação** [PROGRAMA] em todas as nuvens
testadas (3×5 a 13×25, perfis Wendland e spline quártica), não teorema.

### 5.3 Sinal dos autovalores [TEOREMA]

K = Σ_q w_q c_q c_qᵀ com w_q > 0 é semidefinida; QᵀKQ também. Todo
autovalor do problema reduzido é ≥ 0; valores negativos só podem ter
magnitude de arredondamento, ~ε‖K‖/μ_min (observado: |λ| ≤ 7e-15).
Nada foi filtrado: os sinais brutos continuam nos CSV.

### 5.4 Limiar de posto [PROGRAMA]

Há um patamar entre ~1e-15 e μ_min/μ_max (8,6e-5 em 3×5; 1,1e-6 em
7×13) no qual posto, espectro e certificados não mudam. O padrão 1e-10
está no patamar em todos os casos medidos, mas a lacuna diminui com o
refino; a decisão deve sempre ser acompanhada dos certificados acima.

## 6. Problemas de codificação encontrados

| ID | Grav. | Arquivo:linha (base `c217aa5`) | Problema | Evidência |
|---|---|---|---|---|
| C1 | MÉDIO | `src/mls.cpp:110-112,130,157` | Base polinomial centrada no centro da nuvem e escalada pelo seu tamanho: cond(A) cresce como (L/h)^{2p}, a métrica `condition_number` exportada mede a escala da nuvem e não a geometria local, e grau 2 perde precisão e é rejeitado (>1e12) em nuvens finas | grau 2, f=3,3: cond 3,9e3 → 8,0e4 → 1,5e6 (5×3→17×9); spline quártica f=2,0, 15×17: defeito de Δ(x²+y²) 2,4e-3 |
| C2 | MÉDIO | `src/rot_irrot.cpp:76-144` | Redução de posto sem certificado de equivalência com o pencil original; testes só cobriam resíduo projetado | §5.2, `m5_rank_sweep.csv` |
| C3 | MÉDIO | `src/spectral_export.cpp:272` | JSON de metadados grava `"weight": "Wendland C2"` e `"centered scaled polynomial"` fixos, independentemente do perfil | leitura do código |
| C4 | MÉDIO | `src/mls.cpp` (ausência) | Peso da Eq. (4-2) da fonte não implementado; impossível testar o perfil do autor | §9 |
| C5 | BAIXO | `src/rot_irrot.cpp:12`, `apps/rot_irrot_rect.cpp:53,59` | M5 limitado a 120 nós e grau 1/Wendland fixos no executável: impedia estudos de convergência e de perfil | — |
| C6 | MÉDIO | `src/spectral.cpp:62-63,75` | Corte de zero 1e-9·λ_max cresce com o refino; no cartesiano rotula como `near_zero` modos espúrios quase-gradiente (λ 6e-8…2e-6) além do núcleo exato (5): 6 em 13×11 e 9 em 17×15 | `efg/*/eigenvalues.csv` da campanha M6a |
| C7 | BAIXO | `scripts/independent_validation.py:127-138` | Dubleto procurado apenas entre **pares** de autovetores; em 17×15 o conteúdo físico está em 3 vetores (§8.2) | `doublet_17x15.csv` |
| C8 | BAIXO | `src/scalar_control.cpp:190` | `div=-s.dxy+s.dxy` é zero por construção (já documentado como identidade, não validação) | leitura do código |
| C9 | BAIXO | `apps/constraint_subspace.cpp:94` | razão normal com constante 3,0 fixa (=perímetro/área do guia 2×1) | leitura do código |
| C10 | BAIXO | `src/mls.cpp:126` | `evaluate` percorre todos os nós (sem busca espacial): custo O(N) por ponto | desempenho, não correção |

Não foram encontrados erros de sinal, índice ou peso nos blocos K/M
cartesianos e rotacionais, nas derivadas MLS (conferidas por diferenças
finitas com estudo de passo) nem na quadratura.

## 7. Correções implementadas

| ID | Correção | Arquivos |
|---|---|---|
| C1 | Base deslocada para o ponto de avaliação e escalada pelo maior raio de suporte. Pela invariância afim da MLS os valores e derivadas são idênticos em aritmética exata; cond(A) torna-se geométrica (≈1,3e3 em grau 2, ≈26–35 em grau 1) e independente do refino | `src/mls.cpp`, `include/al_maleh/mls.hpp` |
| C2 | Certificados `mass_kernel_stiffness_coupling`, `max_full_residual`, `max_reduced_residual`, `mass_eigen_max_removed` (cálculo e exportação), sem filtrar autovalores | `src/rot_irrot.cpp`, `include/al_maleh/rot_irrot.hpp`, `apps/rot_irrot_rect.cpp` |
| C3 | JSON grava o peso real e a nova descrição da base | `src/spectral_export.cpp` |
| C4 | `MLSWeight::QuarticSpline` (Eq. 4-2) com derivadas até 2ª ordem; padrão continua Wendland C² | `src/mls.cpp`, `include/al_maleh/mls.hpp` |
| C5 | Limite M5 da biblioteca 600 nós; argumentos opcionais `degree` e `weight` em `rot_irrot_rect` e `weight` em `tm_validation`; `--weight` em `m2_ch5_benchmark.py` e `m5_ch6_campaign.py`. Padrões inalterados | `apps/*.cpp`, `scripts/*.py` |
| M7 | Rastreamento modal por sobreposição L² em regra espacial comum, diagnóstico de ângulos principais e atribuição Hungarian; clusters degenerados são sinalizados como ambíguos, sem rastreamento automático de subespaços | `include/al_maleh/mode_tracking.hpp`, `src/mode_tracking.cpp`, `apps/mode_tracking.cpp` |
| — | Estudos reprodutíveis e campanha | `apps/independent_audit.cpp`, `scripts/independent_audit_campaign.py` |

C6–C10 **não** foram alterados (mudariam critérios de classificação
históricos ou são de baixo impacto); ver §10 e issues sugeridas.

## 8. Comparação antes/depois e testes

### 8.1 Impacto numérico da correção C1 [PROGRAMA]

`archive_regression_m6a.csv` recalcula as nuvens da campanha M6a
arquivada e compara com os espectros publicados:

| Nuvem | máx Δλ/λ (1<λ<60) | máx |Δλ| (todo o espectro) | cond. MLS antes → depois | núcleo |
|---|---:|---:|---:|---:|
| 9×7 | 7,1e-13 | 2,7e-10 | 265,6 → 33,9 | 5 = 5 |
| 13×11 | 1,2e-12 | 3,5e-09 | 629,3 → 34,7 | 6 = 6 |
| 17×15 | 9,2e-12 | 2,3e-08 | 1148,9 → 35,0 | 9 = 9 |

As diferenças estão no nível de arredondamento (os maiores |Δλ| ocorrem
nos autovalores ~10³, relativos ≤ 1e-11). Conclusões científicas das
campanhas M2–M6f não mudam. O teste de condicionamento novo falha com a
base antiga (cond. 4,2e3 → 2,7e7) e passa com a nova.

### 8.2 Testes

- **Antes:** 28/28 CTest aprovados.
- **Depois:** **33/33** CTest aprovados (Release). Novos: `m7_mode_tracking_invariance`,
  `m7_mode_tracking_smoke`, `audit_doublet_smoke`, `audit_rank_sweep_smoke`,
  `independent_audit_campaign_archive`. Testes ampliados:
  `mls_consistency` (spline quártica, estudo de passo de diferenças
  finitas com regime de truncamento O(h²) e de arredondamento, cantos e
  bordas da nuvem retangular, condicionamento independente do refino,
  consistência de 2ª ordem), `m5_ch6_basis_and_rank` (teorema do pencil
  singular graus 1 e 2, certificados, resíduo projetado insuficiente,
  Tabela 4 rot/irrot), `generalized_spectrum` (Tabela 4 cartesiana e
  Tabela 2).
- **Sanitizers:** build `RelWithDebInfo` com
  `-fsanitize=address,undefined`: 20/20 testes C++ e os executáveis de
  auditoria sem erros.
- **Determinismo:** duas execuções completas da campanha produziram os
  117 arquivos de dados bit a bit idênticos (diferem apenas 2 logs, pelo
  caminho de saída, e `provenance.json`, pela data).
- **Integridade:** manifesto M6a histórico verificado (168 arquivos);
  manifesto da nova campanha verificado (120 arquivos).
- Transcrições do autor: 52/52 entradas das Tabelas 2, 3, 5, 6 conferem
  com o PDF.

## 9. Impacto sobre os resultados científicos

### 9.1 Perfil MLS do Capítulo 4 [INFERÊNCIA forte]

Varredura de r/h com o peso da Eq. (4-2) contra os 15 autovalores não
nulos cartesianos e os 5 rotacionais/irrotacionais da Tabela 4 (três
quadraturas, 60 valores): mínimo único em **r = 1,5 × espaçamento
nodal, base linear**, com erro relativo máximo 0,46% (cartesiano) e
0,38% (rot/irrot). A conferência individual da última casa impressa
aprova 45/45 entradas cartesianas e 12/15 M5; os três primeiros
autovalores não nulos M5 ficam fora do intervalo de arredondamento.
A Tabela 4 é, portanto, reprodução **aproximada**, não identidade numérica;
f=1,49 ou 1,51 já pioram ~20×. Wendland C² e base quadrática não
reproduzem a tabela na faixa de raios e passo pesquisados. O raio **não está escrito na
dissertação**: é uma identificação inversa, registrada como tal.

### 9.2 Reprodução das Tabelas 2, 3, 5 e 6 [PROGRAMA]

Com esse perfil e **o nosso código** (M2 e M5 com redução de posto):

- M2: **28/28** entradas reproduzidas com |Δ| < 0,001 pp.
- M5: **21/24** com |Δ| ≤ 0,008 pp pelo critério de maior correlação de
  curl; 1 coincide pelo λ mais próximo (F7); 1 coincide em módulo com
  sinal oposto (F6, hipótese de erro de sinal, não comprovada); 1 diferença
  de 0,227 pp sem explicação (F8).

A classificação operacional usa 0,01 pp, que é mais larga que meia unidade
da última casa publicada (0,0005 pp). Os limites <0,001 pp e ≤0,008 pp
descrevem proximidade, não igualdade na precisão impressa. O perfil foi
fixado com a Tabela 4 antes de conferir as demais tabelas, constituindo
validação cruzada parcial; não identifica unicamente todos os parâmetros
ou o solver do autor.

Consequências: (i) as equações (5-1)–(5-11) e (6-1)–(6-6) estão
implementadas como na fonte; (ii) a redução de posto da massa aproxima o
espectro não nulo e recupera os 14 zeros M5 da Tabela 4; (iii) **as grandes
discrepâncias das campanhas M2/M5 de 21/10 (até 34 pp) eram devidas ao
perfil MLS (Wendland C², 3,3h), não a erros de formulação**. Os números
históricos permanecem válidos para o perfil em que foram calculados.

### 9.3 Evidência de não convergência M5 no perfil inferido [PROGRAMA + TEOREMA parcial]

A rigidez do M5 usa ∫ΔNᵢΔNⱼ: requer consistência de **segunda ordem**.
A MLS linear reproduz só polinômios de grau 1, e o defeito
Σ q(xᵢ)ΔNᵢ − Δq para q=x²+y² é **independente de h** no interior
[TEOREMA: para malha uniforme, Σ|xᵢ−x|²Nᵢ(x) = h²ĝ((x−x_c)/h) e seu
Laplaciano não depende de h]. [PROGRAMA] RMS interior: 1,81 (45% de 4)
para spline quártica 1,5h; 0,49 para 2,5h; ~0,09 para 3,3h; ~0,07 para
Wendland 3,3h; ~1e-12 com base quadrática.

Com a quadratura de referência mais fina (16 subcélulas × Gauss 5×5), o erro de TM11 do
perfil do autor **satura em ≈ +23%** (26,4 / 24,0 / 23,8 / 23,6 / 23,3%
de 3×5 a 13×25), com ρ de curl estagnado ≈ 0,90; com Gauss 2×2 satura em
≈ +34% e com 1×1 em ≈ −13%. A persistência com quadratura mais fina apoia a hipótese de defeito
de consistência de 2ª ordem, agravado por suporte pequeno. A regra
16 subcélulas × Gauss 5 é a referência mais fina testada, não uma prova
de convergência da quadratura. O defeito pontual, isoladamente, não prova
não convergência espectral; os dados mostram estagnação nesta família. Com base linear e suporte maior os erros
medidos caem (≈1,7% em 2,5h; ≈0,05% em 3,3h na malha mais fina);
a queda ainda observada em 3,3h não estabelece um patamar assintótico.
Com **base quadrática o erro
decresce sistematicamente no refino medido** (Wendland 3,3h: 0,634 → 0,302 → 0,186 → 0,093%).
O mecanismo é coerente com os dados, mas a afirmação geral "base linear
nunca converge" não foi demonstrada para qualquer perfil (só medida).

### 9.4 Modos espúrios

- **Cartesiano:** a reprodução linear garante **ao menos** 5
  campos afins de rotacional nulo. O corte absoluto 1e-9 identifica 5
  autovalores quase nulos nas famílias testadas; isso não prova dimensão
  exatamente 5 para toda nuvem. Os demais gradientes discretos não estão
  no espaço, gerando dezenas de autovalores positivos pequenos (14 abaixo
  de TM11 em 3×5; 49/95/148 em 9×7/13×11/17×15, somando os rotulados
  `near_zero` pelo corte relativo, achado C6).
- **Rot/irrot [TEOREMA + PROGRAMA]:** os gradientes discretos são exatos
  (n−1 zeros). A afirmação do autor de que "certain other spurious modes
  still remain" foi **reproduzida em 3×5** (Tabela 4: 5,91; 12,9; 15,0;
  15,5 abaixo do TM11=16,8) e também ocorre em 4×7 (Gauss 2×2: λ=12,279
  abaixo do TM11 identificado em 16,547). Nas grades de 5×9 a 13×25 da
  campanha, em todos os perfis e regras, nenhum autovalor não nulo fica
  abaixo do TM11 identificado pelo campo (que coincide com o mais
  próximo de 5π²/4); nelas o problema principal é a imprecisão (§9.3).

### 9.5 Dubleto TM22/TM41 (M6b) — conclusão invariante [PROGRAMA]

`doublet_17x15.csv`: o autovetor 192 é físico (cosseno 1,00000 com o
subespaço analítico, div 0,0075); o 193 tem div 0,415 e fração gradiente
0,016. O vizinho **espúrio 191** (λ=49,204, fração gradiente 0,96, gap
relativo 3,1e-3) contém cosseno 0,128 (≈1,6% de energia) do dubleto
físico: os dados apoiam a **hipótese de hibridização por quase-degenerescência** entre um modo
quase-gradiente e o dubleto. A melhor aproximação L² do dubleto analítico
no mesmo espaço MLS tem div ≤ 0,0045, e o melhor subespaço 2D de
{190,…,193} recupera o dubleto com cossenos (1,00000; 0,99999) e
div ≤ 0,0145. O espaço MLS tem capacidade de aproximação melhor que
a dos autopares selecionados; isso é compatível com a diferença entre
13×11 e 17×15, mas não demonstra sua causa. As combinações da janela
com λ distintos e as projeções L² são aproximações geométricas; não
são novas soluções exatas do autoproblema nem corrigem os autopares.

Invariante à base do subespaço {192,193}: autovalores, ângulos principais
com o analítico, extremos de div (0,0075–0,4153) e de fração gradiente.
Dependente de representação: rótulos "TM22"/"TM41" por vetor e as
correlações individuais com cada modo analítico degenerado. Como
λ192≠λ193 (separação 1,8e-4), os dois autovetores são individualmente
determinados (a menos de sinal) pelo operador discreto. A alta correlação
EFG×Whitney do M6f (mín. 0,983) é compatível com 1,6% de contaminação e
**não** elimina a falha de divergência.

### 9.6 Identidade modal M7 (TM21, 3×5, f 3,00→3,60) [PROGRAMA]

Varredura com passo 0,02 e rastreamento por L² de H em regra comum
(32×16 células, Gauss 3): 210 transições com λ<30, **nenhuma ambígua**,
sobreposição mínima 0,956 — todos os ramos variam continuamente. O
"candidato de maior correlação de curl" muda de ramo em f=3,30, 3,46 e
3,50; seu erro varia de −20,0% a +34,0% com ρ entre 0,58 e 0,74. Há um
cruzamento de ordem perto de λ≈14,7 (f 3,06→3,12) e troca de caráter
(quase-cruzamento) entre os ramos ≈16,4 e ≈17,2 (f≈3,42–3,48). O conteúdo
TM21 está repartido entre 2–3 autovetores (cosseno principal do par
0,82–0,99). Conclusão: o salto +31,465% → −8,503% da campanha de 21/10 é
**troca do autovetor selecionado**, não descontinuidade espectral nem
efeito do posto (constante 26); na malha 3×5 nenhum autovetor representa
bem o TM21, e o erro "do TM21" depende da regra de seleção.

### 9.7 Comparação FEM Whitney (M6e/M6f)

Conferidos: escala λ_EFG = λ_edge/4 (guia 1×0,5 → 2×1), mapeamento de
coordenadas ×2, campo analítico H = (−∂ᵧψ, ∂ₓψ), pesos de área e ângulos
principais por SVD de bases L²-ortonormalizadas. O valor Whitney no
centroide é exato como média do elemento (campo afim por triângulo); o
EFG é amostra pontual. A documentação já declara M6f como exploratório
(1 ponto por triângulo, não quadratura convergida).

## 10. Problemas ainda não resolvidos

1. F8 (TM11 4×7 Gauss 1×1, 0,227 pp) sem explicação.
2. Origem dos "negative eigenvalues" e dos "fifteen exact zeros" do autor
   (F4/F5): exige o código original; hipótese do pencil singular não
   deflacionado não testada.
3. C6: corte de zero relativo a λ_max no cartesiano mistura núcleo exato
   e espúrios. Alterar o critério muda classificações históricas: não foi
   feito; proposta em issue.
4. C7: buscar o dubleto em subespaços de janela (não só pares) mudaria o
   veredito do M6a; mantido, documentado.
5. Convergência M5 com base quadrática foi medida até 13×25; não há prova
   de taxa nem de ausência de espúrios em geral (compactness discreta não
   demonstrada).
6. M6f continua por centroides; quadratura comum convergida pendente.

## 11. Próximos experimentos propostos

1. Rodar as campanhas M2/M5 (Tabelas 2–6) com perfil quártico 1,5h e
   base quadrática para isolar o efeito da consistência de 2ª ordem sobre
   os modos espúrios de 3×5.
2. Critério de núcleo cartesiano baseado em lacuna espectral/arredondamento
   (ε‖K‖/μ_min) reportado **ao lado** do atual.
3. Busca de dubleto por melhor subespaço 2D de janela espectral no M6a/M6b
   (nova campanha, sem alterar a histórica).
4. Rastreamento M7 do perfil do autor em grades 4×7 e 5×9 e de TM22/TM41
   cartesiano entre 13×11 e 17×15 para medir a aproximação espúrio-físico.
5. M6f com Gauss comum em cada triângulo e refino independente da malha FEM.
6. Avaliar a hipótese F5 resolvendo o pencil singular completo com QZ
   (LAPACK `dggev`) e comparando os 4 autovalores indeterminados.

## Arquivo compacto e validação da revisão

A pasta histórica arquiva espectros, diagnósticos, correlações, varreduras
comprimidas, logs e proveniência, conforme seu manifesto. Matrizes K/M,
coeficientes e campos completos M2/M5 foram removidos pela compactação
e são regeneráveis pelos drivers. Os relatórios internos M2/M5 descrevem
a saída completa antes dessa remoção. A revisão final e a nova execução
CI/local são registradas separadamente em [PR20_FINAL_REVIEW.md](PR20_FINAL_REVIEW.md);
o REPORT.md histórico mantém o texto original, inclusive a interpretação
de arredondamento corrigida acima.

## Reprodução

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
/usr/bin/ctest --test-dir build --output-on-failure
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 python3 scripts/independent_audit_campaign.py --out build/independent_audit_results
python3 scripts/validation_manifest.py --verify --campaign build/independent_audit_results
```

A campanha leva ~1,5 min em 2 núcleos. Ela recusa diretórios de saída
não vazios e caminhos que são arquivos para nunca sobrescrever resultados.
