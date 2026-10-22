# M7 — Rastreamento de identidade modal e ferramentas da auditoria de 2026-10-22

Contexto científico e resultados: [relatório da auditoria independente](../results/INDEPENDENT_THEORY_CODE_AUDIT.md).

## Por que rastrear modos

A posição de um autovetor no espectro ordenado não identifica o modo:
quando um parâmetro (raio MLS, quadratura, nuvem) varia, autovalores
cruzam ou quase cruzam e a ordem muda. O critério "maior correlação de
curl com TM(m,n)" também pode mudar de ramo sem que nenhum autovalor
salte. O M7 segue cada autovetor pelo **campo físico** H.

## Método (`include/al_maleh/mode_tracking.hpp`)

- `normalized_overlaps(A,B,w)`: |⟨aᵢ,bⱼ⟩_w|/(‖aᵢ‖_w‖bⱼ‖_w), invariante a
  sinal e escala de cada coluna.
- `principal_cosines(A,B,w)`: cossenos dos ângulos principais entre
  span(A) e span(B) no produto L² ponderado (QR com pivotamento + SVD);
  invariante a A→AR e B→BS não singulares; colunas dependentes são
  descartadas. É a comparação correta para autoespaços (quase)
  degenerados.
- `max_overlap_assignment(O)`: associação um-a-um de máxima sobreposição
  total (algoritmo Hungarian), retangular quando necessário.
- `eigenvalue_clusters(λ,tol)`: agrupamento por lacunas relativas.

O executável `mode_tracking` resolve cada parâmetro, amostra H e curl(H)
de todos os modos com 0<λ≤λ_max numa **regra espacial fixa** (32×16
células, Gauss 3×3, independente da nuvem e do raio), associa modos
entre fatores consecutivos por L² de H e registra:

- `spectrum.csv`: λ e sobreposições com os campos analíticos TM11, TM21,
  TM31, TM12, TM22, TM41 (H e curl);
- `tracking.csv`: associação, sobreposição, segunda melhor (margem) e
  status `tracked` (≥0,9 e margem ≥0,2) ou `ambiguous`;
- `selection.csv`: candidato de maior ρ de curl e de λ mais próximo, se o
  candidato de ρ é o sucessor rastreado do anterior (`0` marca troca de
  ramo) e o cosseno principal do campo analítico com o par dos dois
  melhores candidatos (mistura modal).

```bash
./build/mode_tracking m5 5 3 2 2 1 wendland_c2 40 out/m7 3.0 3.1 3.2 3.3 3.4 3.5 3.6
```

Argumentos: formulação (`m2`|`m5`), nx, ny, Gauss, subdivisão, grau,
peso, λ_max, saída, fatores de suporte.

Testes: `tests/test_mode_tracking.cpp` (invariância a sinal/escala,
produto ponderado, invariância a rotação/escala de base dos ângulos
principais, ângulo conhecido, colunas dependentes, Hungarian vs. guloso,
caso retangular, agrupamento) e `m7_mode_tracking_smoke`.

## `independent_audit` (estudos numéricos reprodutíveis)

| Subcomando | Saída |
|---|---|
| `ch4-scan weight degree fmin fmax fstep out.csv` | 20 primeiros autovalores cartesianos e rot/irrot da malha 3×5 (Gauss 2–4) por fator r/h |
| `laplacian out.csv` | defeito Σq(xᵢ)ΔNᵢ−4 para q=x²+y², cond(A), por perfil e malha |
| `convergence out.csv` | M5 TM11/TM21 com quadratura de produção e convergida, certificados |
| `rank-sweep out.csv` | posto, nulidade, resíduos projetado/completo e ‖KZ‖/‖K‖ por limiar |
| `doublet nx ny out.csv` | dubleto TM22/TM41 cartesiano: vetores, subespaço, projeção L² do analítico e melhor subespaço 2D da janela espectral |

`scripts/independent_audit_campaign.py --out DIR` executa todos os
estudos, as campanhas M2/M5 com o perfil identificado (`--weight
quartic_spline --support-factor 1.5`), o rastreamento M7 e a regressão
dos espectros arquivados do M6a; grava `REPORT.md`, `provenance.json` e
`manifest.json` (verificável por `scripts/validation_manifest.py
--verify`). Recusa diretório não vazio.

## Perfis MLS

`MLSEvaluator(nodes, degree, max_condition, weight)` aceita
`MLSWeight::WendlandC2` (padrão histórico do projeto) e
`MLSWeight::QuarticSpline` (Eq. 4-2 da dissertação). O raio r=1,5h do
perfil do autor é **inferido** pela varredura contra a Tabela 4; a tese
não o informa.

Executáveis com perfil configurável (argumentos opcionais finais, padrões
inalterados):

```bash
./build/tm_validation 5 3 1 1.5 2 out/m2_q 2 quartic_spline
./build/rot_irrot_rect 5 3 2 1.5 2 out/m5_q 1e-10 1 quartic_spline
```
