# M6f — primeira comparação vetorial TM-H EFG × Whitney

**Estado experimental:** correlações sobre centroides de triângulos Whitney. Esta amostragem barata não é quadratura física convergida; os ângulos obtidos são diagnósticos preliminares. A issue #12 define a validação completa.

## Método

O método compara os campos Hx,Hy nos mesmos pontos físicos, sem alterar K, M ou quaisquer autovetores. A fonte FEM é o executável edge_rect independente, com modelo TM de H transversal, compilado no commit dc1d9545c1238554c2292aa0f64a3c1b883d4118 da campanha M6e. Os vetores EFG são produzidos com a mesma solução cartesiana M3–M6a.

O guia FEM nativo tem dimensões 1×0,5 e o EFG 2×1. Transformam-se ambas as coordenadas FEM segundo x_EFG=2*x_FEM, y_EFG=2*y_FEM. As áreas dos triângulos são multiplicadas por quatro. Isso preserva o posicionamento relativo dos campos; as amplitudes de ambas as formulações são renormalizadas em uma mesma norma discreta L2. Os campos nativos FEM são normalizados por pico, mas isso não afeta as correlações normalizadas.

Para TM11, TM21, TM31 e TM12 calculam-se as correlações absolutas e erros L2 orientados, além da concordância de cada método com a referência H=(-∂y psi,∂x psi), psi=sin(mπx/2)sin(nπy). Para TM22/TM41, os dois vetores de cada método são QR-ortonormalizados na quadratura discreta, e são obtidos os valores singulares do produto cruzado e dois ângulos principais, invariantes à rotação arbitrária de base. Nunca rotular individualmente as duas colunas do dubleto como TM22 ou TM41.

A correspondência inicial dos autovetores vem da campanha M6e, por proximidade de autovalores; não foi validada ainda por busca global de sobreposição de campos. Os resultados completos, arquivos selecionados, coordenadas e hashes constam dos artefatos de M6e/M6f.

## Arquivos e execução

Aplicativo C++: apps/efg_sample_vector_fields.cpp — avalia todos os seis campos EFG nos centroides efetivos do FEM edge. Script Python: scripts/m6f_centroid_overlap.py. Testes: tests/test_m6f_centroid_overlap.py.

A campanha grava REPORT.md, simple_modes.csv, doublet_subspaces.csv, efg_samples.csv por caso, logs, índices de modos e hashes SHA256 dos arquivos FEM originais.

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel 2
    /usr/bin/ctest --test-dir build --output-on-failure
    PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 python3 scripts/m6f_centroid_overlap.py --efg-binary build/efg_sample_vector_fields --m6e-campaign out/m6e --out out/m6f

## Limitações e próximo passo

Não se deve confundir amostragem nos centroides com validação do campo contínuo, nem uma boa correlação espacial com a satisfação das equações de Maxwell. Faltam: reconstrução de Whitney/EFG em grade de quadratura Gauss comum com convergência, divergência H, fluxo normal H nas paredes, e avaliação da direção física pior do dubleto. O M6b permanece reprovado no caso fino por divergência 0,415325.

A investigação M5 rotacional/irrotacional é uma frente separada, prioritária após cotejo com a dissertação original, conforme issue #11.
