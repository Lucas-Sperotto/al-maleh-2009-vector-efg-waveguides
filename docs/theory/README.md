# Documentação teórica — Reprodução de Al Maleh (2009)

Edição: 20 de outubro de 2026. A documentação é uma **reconstrução matemática própria** orientada à implementação, e **não** uma transcrição das equações internas da dissertação: o resumo foi conferido, mas ainda não se obteve acesso integral verificável ao PDF para cotejamento página a página.

1. `01_maxwell_tm_forma_fraca.md`: campos TM, fronteiras PEC, formulação fraca e autovalores retangulares.
2. `02_mls_efg_derivadas.md`: MLS com derivadas totais, Hessiana, pesos, estabilidade e testes.
3. `03_formulacoes_vetoriais.md`: matriz curl-curl cartesiana e modelo rotacional/irrotacional **a verificar no original**.
4. `04_autovalores_espurios_validacao.md`: núcleo, massa singular, análise de modos espúrios e benchmarks.
5. `05_arquitetura_reuso_rastreabilidade.md`: organização C++, auditoria de código já existente e plano de implementação.

Cada fonte Markdown é adequada para controle de versões; os correspondentes documentos `.docx` foram gerados com equações editáveis por conversão Pandoc.

Fonte principal: Hussam Al Maleh. *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems*. McGill University, 2009. https://www.collectionscanada.gc.ca/obj/thesescanada/vol2/002/MR61527.PDF