# Vector Element-Free Galerkin — Autovalores de guias de onda

Reprodução crítica da dissertação de **Hussam Al Maleh (2009)**, *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems* (McGill University). O projeto estuda aproximações vetoriais sem malha para modos **TM** de guias PEC, modos próprios espúrios e comparação com soluções analíticas e FEM.

> **Rigor científico:** a documentação inclui reconstruções matemáticas voltadas à implementação, não uma transcrição integral verificada da dissertação. A identificação exata da segunda formulação rotacional/irrotacional depende de cotejo com o PDF. Testes MLS aprovados não validam ainda autovalores eletromagnéticos.

## Estado do desenvolvimento

| Etapa | Estado |
|---|---|
| M1 — MLS linear/quadrático, gradientes e Hessianas | **Implementado e testado** (CTest e GitHub Actions) |
| M2 — Formulação cartesiana, quadratura, matrizes de massa e curl-curl | Pendente |
| M3 — Autovalores TM e campos de guia retangular PEC | Pendente |
| M4 — Diagnóstico dos modos espúrios | Pendente |
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

`mls_probe` exporta as funções MLS e suas derivadas. **Ainda não há solver de Maxwell/guia de ondas implementado.**

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
