---
title: "05 — Arquitetura C++, reutilização de código e rastreabilidade científica"
author: "Projeto de reprodução de Al Maleh (2009)"
date: "20 de outubro de 2026"
lang: pt-BR
---

# 1. Política do projeto

Este repositório tem objetivo específico: uma **reprodução crítica** da dissertação de Hussam Al Maleh (2009). Não alterar a implementação principal do projeto NASA TP-3485 nem fundir o experimento à base Whitney difusa antes de obter resultados independentes. Registrar separadamente (A) comportamento demonstrado no resumo da tese, (B) matemática reconstruída, (C) escolhas próprias de engenharia e (D) implementações verificadas contra páginas do original.

A edição de 20/10/2026 é uma **especificação inicial sem evidência de execução numérica**. Não afirmar testes aprovados sem CI e artefatos correspondentes.

## 2. Repositórios inspecionados

| Projeto | Arquivos ou módulos úteis | Tipo de aproveitamento | Risco ou adaptação |
|:--|:--|:--|:--|
| `Lucas-Sperotto/fem-eigen-em-2d3d-whitney` | `src/core`, `src/edge`, `src/helmvec`, `src/explicit`, `scripts`, `docs/results` | geometria, quadratura, solucionador generalizado, exportação e benchmarks FEM | não misturar base de arestas FEM com EFG nodal; equivalência física TE/TM antes de comparar |
| `Lucas-Sperotto/simple-waveguide-efg-fem-solver` | `mls_shape.cpp`, `efgm_solver.cpp`, `fem_solver.cpp`, `waveguide_modes.py`, `plot_eigenvalues.py` | protótipos, geração de nuvem, comparação escalar e gráficos | `mls_shape.cpp` omite derivadas de coeficientes MLS; `efgm_solver.cpp` penaliza ambas matrizes $K$ e $M$; peso gaussiano truncado apresenta salto |
| `Lucas-Sperotto/efgmi-hcurl-diffuse-whitney` | `docs/theory`, `docs/validation`, implementação de MLS/Whitney e `TODO.md` | testes algébricos, diagnósticos de núcleo, experimentos de funções vetoriais compatíveis | política `all_active_pairs` redundante e diagnóstico de núcleo ainda não aprovado; base conceitualmente distinta |
| `Lucas-Sperotto/efg-3d-electromagnetics-2006` | `include/efg`, `src`, estrutura CMake/CTest | padrão didático de organização C | documento do repositório descreve projeto ainda em esqueleto; não copiar solver inexistente |

Os quatro repositórios possuem arquivo `LICENSE` MIT (inspecionado). Em qualquer cópia real de código, preservar avisos de copyright, texto de licença e documentação da origem do módulo. Dependências terceirizadas permanecem submetidas às suas próprias licenças.

## 3. Organização sugerida

```text
al-maleh-2009-vector-efg-waveguides/
  README.md
  LICENSE
  CMakeLists.txt
  include/al_maleh/{geometry,mls,quadrature,vector,eigen,boundary,io}/
  src/{geometry,mls,quadrature,vector,eigen,boundary,io}/
  apps/{tm_reference,tm_cartesian,tm_rot_irrot,compare}/
  tests/{mls,derivatives,quadrature,vector_basis,spectral,benchmarks}/
  docs/theory/{01...,02...,03...,04...,05...}.md
  docs/references/
  benchmarks/rectangular/
  scripts/
  results/                  # fora do controle de versão, exceto resumos curados
```

Recomendações: C++17, CMake/CTest, Eigen para vetores/matrizes locais, LAPACK ou Spectra/ARPACK para problemas maiores, Gmsh apenas para **integração e visualização** quando adequado, VTK/CSV para campos. Matriz esparsa em CSR/CSC; montagens densas reservadas a testes pequenos.

## 4. Interfaces a definir antes de copiar código

```cpp
struct Point2 { double x{}, y{}; };
struct Node { int id{}; Point2 x; };
struct IntegrationPoint { Point2 x; double weight{}; };

struct MLSConfig {
  int polynomial_degree{};
  double support_radius{};
  WeightType weight{};
};

struct MLSData {
  std::vector<int> active;
  std::vector<double> phi, dx, dy, dxx, dxy, dyy;
  double condition_number{};
};

struct Eigenpair {
  double lambda{};
  Vector coefficients;
  double relative_residual{};
};
```

A função `evaluate_mls(point, cloud, config)` deve calcular a avaliação completa (ou retornar status explícito) e ser reutilizada por todas as formulações. A classe de montagem não deve recalcular derivadas de outra maneira. Fornecer políticas de DOFs vetoriais distintas: duas componentes por nó em cartesiano; outras funções após cotejo do texto original.

## 5. Mapa de rastreabilidade: matemática para código e teste

| Relação teórica | Módulo proposto | Teste obrigatório |
|:--|:--|:--|
| $Aq_I=B_I$ | `src/mls/mls.cpp` | posto, estabilidade, PU e reprodução polinomial |
| $Aq_{I,\alpha}=B_{I,\alpha}-A_{,\alpha}q_I$ | `src/mls/mls_derivatives.cpp` | diferenças finitas e patch linear |
| $Aq_{I,\alpha\beta}=B_{I,\alpha\beta}-A_{,\alpha\beta}q_I-A_{,\alpha}q_{I,\beta}-A_{,\beta}q_{I,\alpha}$ | `src/mls/mls_hessian.cpp` | derivadas mistas e patch quadrático |
| $\mathbf N_i^{x,y}$ e rotacionais | `src/vector/cartesian_basis.cpp` | curl por diferenciação independente |
| $K_{IJ}=\int\operatorname{curl}N_I\operatorname{curl}N_J$ | `src/assembly/curlcurl.cpp` | simetria, PSD e quadratura |
| $M_{IJ}=\int N_I\cdot N_J$ | `src/assembly/mass.cpp` | simetria, posto/definitude, escala geométrica |
| $Kq=\lambda Mq$ | `src/eigen/generalized.cpp` | resíduos e autoproblema manufaturado |
| $\lambda_{mn}=(m\pi/a)^2+(n\pi/b)^2$ | `src/validation/rect_tm.cpp` | TM11, TM21 e degenerados |
| $\nabla\cdot H=0$ para modos positivos | `src/validation/fields.cpp` | norma da divergência e fluxo normal |
| $N^I=\nabla\phi$, $N^R=J\nabla\phi$ | `src/vector/rot_irrot_basis.cpp` | Hessiana, núcleo, rank de massa — **experimental até cotejo do original** |

## 6. Estratégia de cópia dos projetos existentes

**Pode ser reaproveitado com testes de integração:** gerador e leitor de geometrias, quad Gauss, infraestrutura CMake/CTest, exportação VTK/CSV, solver LAPACK e scripts de comparação de espectros. Utilizar cópia seletiva com referência ao caminho e commit de origem, não importar uma árvore inteira por conveniência.

**Deve ser reimplementado ou corrigido:** avaliação MLS e derivadas; escolha de pesos compactos e suaves; tratamento das condições de contorno do campo magnético TM; montagem vetorial cartesiana e rotacional/irrotacional. O `mls_shape.cpp` legado não implementa derivada total de $\phi_I$; o peso `weight_gaussian()` truncado não é contínuo no raio de corte. A `efgm_solver.cpp` antiga modifica também a matriz de massa para tratar fronteira e deve permanecer apenas como referência histórica.

**Não reutilizar como se fosse equivalente:** a base Whitney FEM com DOFs de aresta e a base Whitney difusa experimental não são as funções nodais cartesianas da primeira formulação de Al Maleh. Elas pertencem à etapa de comparação cruzada, não à definição do método reproduzido.

## 7. Plano de entregas e critérios de bloqueio

1. **M0 — Fonte bibliográfica.** Registrar em `docs/references/` o endereço, metadados da tese, páginas, numeração de equações e tabela de diferenças entre o original e nossa reconstrução. Critério: não usar a expressão “idêntico ao autor” sem citação verificável.
2. **M1 — MLS.** Base, pesos, suportes, derivadas de primeira/segunda ordem e testes de reprodução. Critério: todos os patch tests aprovados e diagnóstico de condicionamento.
3. **M2 — EFG cartesiano.** Montagem esparsa de K/M, tratamento natural/essencial consistente, solver e exportação de espectro. Critério: simetria e massa de posto suficiente; execução com três refinamentos.
4. **M3 — Validação TM.** Solução retangular, campos, divergência, fronteiras, resíduos e modos múltiplos. Critério: erro decrescente e classificação documentada de modos não físicos.
5. **M4 — Segunda formulação.** Somente após conferir a base do capítulo da tese. Critério: fórmulas da dissertação mapeadas e testes de posto antes de resolver autovalores.
6. **M5 — Comparação independente.** Confrontar com projeto NASA e EFGMI-H(curl) em geometria, normalização e quadratura equivalentes. Critério: relatório reprodutível e gráficos auditáveis.

## 8. Relatório científico automatizado

Para cada execução salvar `config.json`, `run.log`, `eigenvalues.csv`, `modes.csv`, `fields_*.csv`, `matrices_summary.json`, `diagnostics.json`, `timings.csv` e, quando viável, campos em VTK. Incluir commit, versão do compilador, BLAS/LAPACK, ordem da quadratura, base MLS, peso, raio de suporte, tolerância de posto, restrições de fronteira e dimensionalidade. Comparar `lambda`, `kc` e `fc` em colunas distintas.

## 9. Referências de repositório

- https://github.com/Lucas-Sperotto/fem-eigen-em-2d3d-whitney
- https://github.com/Lucas-Sperotto/efgmi-hcurl-diffuse-whitney
- https://github.com/Lucas-Sperotto/simple-waveguide-efg-fem-solver
- https://github.com/Lucas-Sperotto/efg-3d-electromagnetics-2006
- Al Maleh, H. (2009), tese McGill. URL: https://www.collectionscanada.gc.ca/obj/thesescanada/vol2/002/MR61527.PDF