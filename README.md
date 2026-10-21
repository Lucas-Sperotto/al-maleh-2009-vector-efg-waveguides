# Vector Element-Free Galerkin — Autovalores de guias de onda

Reprodução crítica da dissertação de **Hussam Al Maleh (2009)**, *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems* (McGill University). O projeto estuda aproximações vetoriais sem malha para modos **TM** de guias PEC, modos próprios espúrios e comparação com soluções analíticas e FEM.

> **Rigor científico:** a documentação inclui reconstruções matemáticas voltadas à implementação, não uma transcrição integral verificada da dissertação. A identificação exata da segunda formulação rotacional/irrotacional depende de cotejo com o PDF. Testes MLS aprovados não validam ainda autovalores eletromagnéticos.

## Estado do desenvolvimento

| Etapa | Estado |
|---|---|
| M1 — MLS linear/quadrático, gradientes e Hessianas | **Implementado e testado** (CTest e GitHub Actions) |
| M2 — Formulação cartesiana, quadratura, matrizes de massa e curl-curl | **Implementado e testado** (protótipo de montagem) |
| M3 — Solucionador generalizado, referências TM e exportação de campos | **Protótipo implementado e testado**; convergência física ainda pendente |
| M4 — Comparação modal, núcleo curl-curl e campanha de refinamento | **Implementado e testado**; classificação física ainda exploratória |
| M5 — Formulação rotacional/irrotacional | Pendente de verificação bibliográfica |
| M6 — Comparação independente com NASA TP-3485 e EFGMI-H(curl) | Pendente |

O M1 usa C++17, Eigen3, bases de grau 1 e 2 e peso Wendland C² como decisões próprias de engenharia. Veja [implementação e resultados dos testes](docs/implementation/01_mls_core.md), [proveniência](docs/implementation/CODE_PROVENANCE.md), [API](include/al_maleh/mls.hpp), [código](src/mls.cpp) e [testes](tests/test_mls.cpp).

## Compilar e testar

Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libeigen3-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/mls_probe > mls_probe.csv
```

`mls_probe` exporta funções MLS e derivadas. O executável `cartesian_rect` monta as matrizes K/M; `cartesian_spectrum` resolve o autoproblema generalizado em dimensões pequenas e exporta espectro, referências TM e campos preliminares:

```bash
./build/cartesian_spectrum 5 5 1 3.3 out/m3/rect
```

**A identificação física e a convergência dos modos TM ainda não foram estabelecidas.** Veja [implementação M3](docs/implementation/03_generalized_eigenproblem.md), [solucionador](src/spectral.cpp), [testes](tests/test_spectral.cpp) e [exemplos CSV no GitHub Actions](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions).

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
