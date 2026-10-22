# Vector Element-Free Galerkin — Autovalores de guias de onda

Reprodução crítica da dissertação de **Hussam Al Maleh (2009)**, *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems* (McGill University). O projeto estuda aproximações vetoriais sem malha para modos **TM** de guias PEC, modos próprios espúrios e comparação com soluções analíticas e FEM.

> **Rigor científico:** a documentação inclui reconstruções matemáticas voltadas à implementação. A [auditoria independente de 2026-10-22](docs/results/INDEPENDENT_THEORY_CODE_AUDIT.md) cotejou o **PDF completo** da dissertação (Capítulos 4–6): o peso MLS do autor é a spline quártica da Eq. (4-2); o raio não consta na tese e foi **inferido** (r = 1,5 × espaçamento, base linear) por identificação inversa contra a Tabela 4. Testes algébricos e boa concordância de autovalores não certificam ausência de modos espúrios.

**Auditoria independente teoria × código (2026-10-22, integrada pela PR #20):** [relatório](docs/results/INDEPENDENT_THEORY_CODE_AUDIT.md) e [campanha reprodutível com manifesto SHA-256](docs/results/independent_audit_2026-10-22/REPORT.md). Com o perfil identificado, o código reproduz **28/28** entradas das Tabelas 2–3 (<0,001 pp), o espectro da Tabela 4 nas duas formulações e **21/24** entradas das Tabelas 5–6 (<0,008 pp pelo critério de maior correlação curl). **A coincidência dentro da precisão impressa é menor: 20/28 nas Tabelas 2–3 e 3/24 nas Tabelas 5–6.** Auditoria anterior (histórica): [pós-merges de 21/10](docs/results/POST_MERGE_SCIENTIFIC_AUDIT_2026-10-21.md).

**Estado pós-merge da PR #20:** a auditoria e os testes estão integrados na `main` (commit `0a728d4d`). O [roteiro científico priorizado de 22/10/2026](docs/ROADMAP_POST_PR20_2026-10-22.md) define as próximas etapas M8a–M8d e M9, com issues, riscos, critérios de aceite e reprodutibilidade. **CI verde valida execução, não a ausência de modos espúrios.**

## Estado do desenvolvimento

| Etapa | Estado |
|---|---|
| M1 — MLS linear/quadrático, gradientes e Hessianas | **Implementado e testado** (CTest e GitHub Actions); pesos Wendland C² e spline quártica (Eq. 4-2); base local com condicionamento independente do refino |
| M2 — Formulação cartesiana, quadratura, matrizes de massa e curl-curl | **Implementado e confrontado com o capítulo 5**; Tabelas 2–3 com 28/28 desvios inferiores a 0,001 pp; Tabela 4 com 45/45 raízes cartesianas não nulas dentro do arredondamento, usando perfil MLS parcialmente inferido |
| M3 — Solucionador generalizado, referências TM e exportação de campos | **Protótipo implementado e testado**; convergência física ainda pendente |
| M4 — Comparação modal, núcleo curl-curl e campanha de refinamento | **Implementado e testado**; classificação física ainda exploratória |
| M5 — Formulação rotacional/irrotacional | **Equações (6-1)–(6-6) verificadas; Tabela 4 M5 aproximada (12/15 valores), Tabelas 5–6 com 21/24 desvios <0,008 pp**; pencil singular certificado sob hipóteses de núcleo válido; estagnação MLS linear observada no perfil inferido, sem prova universal |
| M6a–M6f — Validação FEM independente, núcleo gradiente e campos Whitney | **Integrados à main**; falha do dubleto TM22/TM41 17×15 ainda sem correção; hibridização com espúrio é hipótese; M6f por centroides continua exploratório |
| M7 — Identidade modal entre parâmetros | **Implementado e testado para variação de suporte com nuvem/quadratura/posto fixos**: sobreposição L², ângulos principais e Hungarian; clusters degenerados são sinalizados, mas rastreamento automático de subespaços permanece pendente ([documentação](docs/implementation/12_mode_tracking_and_audit_tools.md)) |

O M1 usa C++17, Eigen3, bases de grau 1 e 2, base polinomial deslocada ao ponto de avaliação e dois pesos: Wendland C² (padrão histórico, decisão própria) e spline quártica (Eq. 4-2 da dissertação). Veja [implementação e resultados dos testes](docs/implementation/01_mls_core.md), [proveniência](docs/implementation/CODE_PROVENANCE.md), [API](include/al_maleh/mls.hpp), [código](src/mls.cpp) e [testes](tests/test_mls.cpp).

## Auditoria independente de 2026-10-22 — reprodução

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
/usr/bin/ctest --test-dir build --output-on-failure          # 33 testes
PYTHONNOUSERSITE=1 python3 scripts/independent_audit_campaign.py --out out/independent_audit
python3 scripts/validation_manifest.py --verify --campaign docs/results/independent_audit_2026-10-22
```

Perfil do autor nos executáveis (argumentos opcionais; padrões inalterados):

```bash
./build/tm_validation 5 3 1 1.5 2 out/m2_autor 2 quartic_spline
./build/rot_irrot_rect 5 3 2 1.5 2 out/m5_autor 1e-10 1 quartic_spline
./build/mode_tracking m5 5 3 2 2 1 wendland_c2 40 out/m7 3.0 3.3 3.6
```

## Compilar e testar

Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libeigen3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
/usr/bin/ctest --test-dir build --output-on-failure
./build/mls_probe > mls_probe.csv
```

`mls_probe` exporta funções MLS e derivadas. O executável `cartesian_rect` monta as matrizes K/M; `cartesian_spectrum` resolve o autoproblema generalizado em dimensões pequenas e exporta espectro, referências TM e campos preliminares:

```bash
./build/cartesian_spectrum 5 5 1 3.3 out/m3/rect
```

**Alguns modos TM simples foram identificados e validados por comparação independente, mas o espectro completo contém poluição numérica e a condição física do dubleto TM22/TM41 continua reprovada em divergência.** Veja [implementação M3](docs/implementation/03_generalized_eigenproblem.md), [solucionador](src/spectral.cpp), [testes](tests/test_spectral.cpp) e [exemplos CSV no GitHub Actions](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions).

## M2 — Capítulo 5 original: equações e Tabelas 2–3

As imagens da dissertação de Al Maleh (2009) permitiram confrontar
as **equações (5-1) a (5-11)** com nossa montagem curl-curl/massa,
e transcrever os erros publicados para TM11 e TM21 nas grades
3×5, 4×7, 5×9 e 6×11. O autor usa **quatro células de integração
por retângulo entre nós** e quadratura Gauss de 1×1 até 5×5.
O programa agora suporta essas opções **sem modificar o padrão anterior**.
A auditoria posterior da PR #20 **verificou o peso spline quártica da Eq. (4-2)**, mas inferiu o grau linear e o raio de suporte `r/h=1,5` por identificação inversa. Com esse perfil, 28/28 linhas das Tabelas 2–3 diferem menos de 0,001 pp; 20/28 ficam na precisão impressa. **A campanha histórica abaixo com Wendland C²/3,3h continua preservada e não representa o perfil inferido da fonte.**

```bash
python3 scripts/m2_ch5_benchmark.py --binary build/tm_validation --out out/m2_ch5
```

O relatório compara os valores originais e os nossos erros
de autovalores, preserva o espectro completo e gera campos
vetoriais próprios TM11 e TM21 em CSV/SVG (sem copiar
as figuras originais).

- [Equações originais, Tabelas 2/3 e procedimentos M2](docs/implementation/02_cartesian_assembly.md)
- [Tabela numérica transcrita do capítulo 5](docs/results/m2_author_ch5_tables_2_3.csv)
- [Resultados numéricos efetivamente executados em 14 casos](docs/results/M2_CH5_MEASURED_2026-10-21.md)
- [Artefatos completos M2: matrizes, espectros e SVGs](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39630223717/artifacts/12358927216)

**Ambiente MyHoneyBee:** caso um `ctest` de `~/.local/bin`
produza `ModuleNotFoundError: No module named 'cmake'`,
use `/usr/bin/ctest` (o executável oficial CMake/CTest).
## M5 — Formulação rotacional/irrotacional original (Capítulo 6)

Com as imagens do capítulo 6 foi possível verificar que Al Maleh usa
**grad(N_i)** e **z×grad(N_i)** para cada nó, exatamente como mostram
as equações (6-1)–(6-6). A rigidez envolve o Laplaciano das MLS;
a massa pode ser singular por dependências entre bases, e por isso
o solver M5 faz redução explícita de posto, preservando e documentando
o espectro no subespaço reduzido (inclusive autovalores negativos).

Os **24 erros numéricos originais** das Tabelas 5 e 6 foram
transcritos. As nossas simulações são publicadas separadamente,
com distinção explícita entre o peso quártico **documentado** e o grau/raio **inferidos**. A PR #20 encontrou 21/24 desvios inferiores a 0,008 pp nas Tabelas 5–6 pelo melhor curl, não igualdade na última casa impressa; F6/F7/F8 permanecem em investigação.

```bash
python3 scripts/m5_ch6_campaign.py --binary build/rot_irrot_rect --out out/m5
```

Os CSV e SVG produzidos incluem TM11, TM21 e os campos irrotacional
e rotacional de TM11, além de espectros completos e relatórios de posto.
Veja [implementação e limitações M5](docs/implementation/11_rot_irrot_original_chapter6.md)
e [dados originais do autor](docs/results/m5_author_ch6_tables_5_6.csv). Veja ainda [os 24 resultados realmente calculados, postos e divergências](docs/results/M5_CH6_MEASURED_2026-10-21.md) e [o pacote integral CI de matrizes, espectros e campos](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39636316113/artifacts/12363260079).
## Validacao espectral M4

O executavel `tm_validation` compara formas de campo (`curl(H)`) e autovalores com os modos analiticos TM do guia retangular. A campanha reproduz tres refinamentos nodais e variacoes de raio de suporte, grau MLS e quadratura:

```bash
python3 scripts/m4_campaign.py --binary build/tm_validation --out out/m4 --profile ci
```

Na primeira campanha de oito casos, **todos foram executados**. O modo TM11 apresentou erro relativo de aproximadamente 0,013%, 0,003% e 0,004% nas nuvens 5x5, 7x7 e 9x7. Entretanto, modos superiores em nuvens grosseiras podem ter autovalores proximos aos analiticos e campos fisicamente inconsistentes. **Os resultados sao diagnósticos, nao certificacao de ausencia de modos espurios.**

- [Metodologia, criterios e limites de M4](docs/implementation/04_tm_validation.md)
- [Resultados medidos do M4](docs/results/M4_initial_campaign.md)
- [Codigo e correspondencia modal](src/tm_validation.cpp)
- [Artefatos reproduziveis no GitHub Actions](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions)

## Teste local com refinamento elevado

Para comparar o espectro e as formas TM do guia retangular PEC 2×1 em três
refinamentos (9×7, 13×11, 17×15 nós), com CMake Release e CTest:

```bash
JOBS=2 PROFILE=local bash scripts/run_local_reference.sh
```

Para o master antigo do Stargate, inicie com as nuvens menores (9×7 e 13×11):

```bash
JOBS=1 PROFILE=stargate bash scripts/run_local_reference.sh
```

Resultados em `out/m4/local/REPORT.md` ou `out/m4/stargate/REPORT.md`,
além de tabelas CSV detalhadas por modo e indicadores de divergência/contorno.
**A comparação é feita com soluções TM analíticas e ainda não equivale à
reprodução numérica de uma tabela específica da dissertação.**
Consulte [o guia completo de execução](docs/implementation/05_local_refinement.md).
Veja também os [resultados medidos em três refinamentos (até 510 DOFs)](docs/results/M4_high_refinement.md).

## Validação independente M6a e espectros completos

A [campanha FEM × EFG de 21/10/2026](docs/results/independent_validation/2026-10-21_dc1d98b/REPORT.md)
executou os dois programas independentes e preservou 2174 autovetores completos,
suas matrizes, mapas de DOFs e campos. Quatro modos simples passam nos critérios
usados; o subespaço TM22/TM41 fino reprova por divergência, apesar da proximidade
espectral. A verificação adicional Eigen × LAPACK nas mesmas matrizes EFG passou.
O [guia de reprodução e formatos](docs/results/independent_validation/README.md)
inclui comandos, testes, manifesto SHA-256 e limites científicos.

## M6b — Diagnóstico de subespaços degenerados (experimental)

As novas formas quadráticas de divergência `D` e fluxo normal `B`
permitem calcular extremos invariantes à escolha de base em TM22/TM41,
sem modificar a matriz K/M nem os arquivos históricos da campanha M6a.

```bash
python3 scripts/m6b_campaign.py --binary build/constraint_subspace --out out/m6b
```

A execução analisa os dubletos das nuvens 9×7, 13×11 e 17×15,
e registra também os candidatos espectrais positivos abaixo de TM11.
**A validação física continua parcial e não declaramos ausência de espúrios.**
Consulte [M6b — método, testes e limites](docs/implementation/06_constraint_subspaces.md).

[Resultados medidos dos três refinamentos do M6b](docs/results/independent_validation/M6B_INITIAL_RESULTS.md).
## M6c — Controle escalar MLS/EFG para modos TM

Uma terceira formulação, **escalar e independente da discretização vetorial**,
aproxima `Ez` com MLS, impõe `Ez=0` por colocalização em nós de contorno
e diagonaliza o problema generalizado reduzido à base de núcleo dessa
restrição. Esse experimento não corrige nem substitui o EFG vetorial.

```bash
python3 scripts/m6c_scalar_campaign.py --binary build/scalar_efg_control --out out/m6c
OPENBLAS_NUM_THREADS=1 python3 scripts/verify_scalar_control.py --campaign out/m6c
```

O autovalor fundamental calculado no refinamento 17×15 foi
`12,336116117`, contra referência TM11 `12,337005501`
(erro `0,007209%`). O autovalor está **ligeiramente abaixo** da referência,
porque a restrição nodal não equivale a Dirichlet contínuo conforme.
Não há proliferação de dezenas de autovalores positivos de baixa frequência
como na formulação cartesiana não restrita. O campo e os contornos
precisam de validação independente antes de alegações mais gerais.

Veja [método e limitações M6c](docs/implementation/07_scalar_efg_control.md) e
[primeiros resultados efetivamente medidos](docs/results/independent_validation/M6C_SCALAR_CONTROL_RESULTS.md).
## M6d — Auditoria do quase-núcleo gradiente

O M6d calcula a projeção L2 dos gradientes escalares MLS no espaço
vetorial cartesiano e mede a fração de energia gradiente de **todos**
os autovetores, sem alterar K/M nem filtrar o espectro.

Nos três refinamentos, os 49/94/144 candidatos positivos abaixo de TM11
apresentaram frações gradientes médias de **91,64%, 97,05% e 98,29%**.
A coluna B do dubleto TM22/TM41 fino, embora reprovada por divergência
no M6b, tem fração gradiente de apenas 1,60%; portanto sua falha
permanece sem explicação completa.

```bash
python3 scripts/m6d_gradient_campaign.py --binary build/gradient_kernel_audit --out out/m6d
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 python3 scripts/verify_m6d_gradient.py --campaign out/m6d
```

Consulte [formulação e limites M6d](docs/implementation/08_gradient_kernel_audit.md)
e [resultados medidos](docs/results/independent_validation/M6D_GRADIENT_KERNEL_RESULTS.md).
## M6e — FEM de aresta Whitney independente

A comparação vetorial agora inclui o executável `edge_rect` (ramo TM
em H_t) do repositório independente [fem-eigen-em-2d3d-whitney](https://github.com/Lucas-Sperotto/fem-eigen-em-2d3d-whitney).
A referência é fixada ao commit `dc1d9545` e compilada separadamente
na CI, sem alterações no código FEM.

O driver FEM edge usa retângulo 1×0,5, enquanto o EFG usa 2×1.
Os autovalores são portanto comparados mediante λ_edge/4.
A etapa compara espectros e referências TM; campos vetoriais completos
e o dubleto degenerado ainda exigem matching por forma/ângulos principais.

Veja [protocolo M6e, fontes, reprodução e limitações](docs/implementation/09_whitney_edge_cross_validation.md) e [resultados medidos](docs/results/independent_validation/M6E_WHITNEY_EDGE_INITIAL_RESULTS.md).
## M6f — campos vetoriais H EFG × Whitney (diagnóstico por centroides)

O experimento M6f avalia os campos TM-H dos dois métodos nos **mesmos**
centroides de triângulos, após transformar as coordenadas do guia FEM
1×0,5 para o guia EFG 2×1. Compara correlações L2 dos modos simples
e ângulos principais invariantes do dubleto TM22/TM41.

**Esta é uma campanha de amostragem exploratória**, não uma quadratura
com convergência nem certificação física de divergência e contorno.
Para reproduzir, execute M6e antes e siga
[o protocolo de comparação dos campos](docs/implementation/10_whitney_efg_vector_overlap.md).

## Documentação

O README serve como índice. As deduções detalhadas, equações, critérios de validação e discussões de arquitetura estão organizados nos arquivos Markdown:

| Tema | Documento técnico | Notas completas extraídas do README antigo |
|---|---|---|
| Maxwell, modos TM e forma fraca | [01 — Maxwell e TM](docs/theory/01_maxwell_tm_forma_fraca.md) | [Notas 01](docs/notes/01_maxwell_e_forma_fraca.md) |
| MLS/EFG e derivadas | [02 — MLS](docs/theory/02_mls_efg_derivadas.md) | [Notas 02](docs/notes/02_mls_e_derivadas.md) |
| Formulação cartesiana e rotacional/irrotacional | [03 — Vetorial](docs/theory/03_formulacoes_vetoriais.md) | [Notas 03](docs/notes/03_formulacoes_vetoriais.md) |
| Modos espúrios e validação | [04 — Espectro](docs/theory/04_autovalores_espurios_validacao.md) | [Notas 04](docs/notes/04_espurios_e_validacao.md) |
| Arquitetura, reuso e plano | [05 — Implementação](docs/theory/05_arquitetura_reuso_rastreabilidade.md) | [Notas 05](docs/notes/05_arquitetura_e_programa_de_pesquisa.md) |

Veja também o [índice teórico](docs/theory/README.md), o [contexto e bibliografia](docs/notes/00_contexto_e_fontes.md), o [mapa de migração](docs/notes/README.md) e a [cópia integral e inalterada do README anterior](docs/archive/README_antes_reorganizacao_2026-10-20.md). Os documentos técnicos existentes não foram modificados nesta etapa.

## Referências e projetos relacionados

- [Al Maleh (2009) — McGill University, PDF](https://mcgill.scholaris.ca/bitstreams/4d79b125-540f-4b38-8c06-722d14eb60bb/download) e [Library and Archives Canada, PDF](https://www.collectionscanada.gc.ca/obj/thesescanada/vol2/002/MR61527.PDF).
- [NASA TP-3485: FEM independente](https://github.com/Lucas-Sperotto/fem-eigen-em-2d3d-whitney).
- [EFGMI-H(curl): Whitney difusas](https://github.com/Lucas-Sperotto/efgmi-hcurl-diffuse-whitney).
- [EFGM/FEM escalar anterior](https://github.com/Lucas-Sperotto/simple-waveguide-efg-fem-solver).
- [EFG eletromagnético 3D (2006)](https://github.com/Lucas-Sperotto/efg-3d-electromagnetics-2006).

Licença [MIT](LICENSE). Preservar direitos autorais e rastreabilidade ao adaptar código.
