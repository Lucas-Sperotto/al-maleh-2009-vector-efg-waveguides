# M6d — Projeção L2 do quase-núcleo gradiente

**Estado: diagnóstico experimental, sem alteração do autoproblema original.**

## Motivação

A campanha M6a identificou 49/94/144 autovalores EFG positivos abaixo
da referência analítica TM11 nas três nuvens. M6b demonstrou a falha
de divergência invariante do dubleto TM22/TM41 fino. M6c testou uma
formulação escalar separada, sem a proliferação dessas raízes adicionais.

No contínuo, curl(grad f)=0. Contudo, o gradiente das funções escalares
MLS não pertence necessariamente ao espaço vetorial formado pelos dois
blocos cartesianos das mesmas funções. Não há prova de que estes
espaços constituam um complexo discreto compatível.

## Formulação diagnóstica

Se Ni são as funções vetoriais e phi_j as funções MLS escalares:

- H_ij = integral_Omega Ni dot grad(phi_j) dOmega
- S_ij = integral_Omega grad(phi_i) dot grad(phi_j) dOmega
- M_ij = integral_Omega Ni dot Nj dOmega (matriz vetorial original)

A melhor aproximação vetorial L2 dos gradientes escalares tem
coeficientes B = M^-1 H, calculados por Cholesky. A Gramiana
dos gradientes projetados é G = H^T M^-1 H.

Diagonalizamos G. Para cada autovalor sigma maior que
tau * sigma_max, construímos q = (M^-1 H)u / sqrt(sigma).
O resultado Q satisfaz Q^T M Q aproximadamente igual a I.
A imagem de Q define o projetor M-ortogonal P = Q Q^T M.

O indicador de participação gradiente de um autovetor v é:

    f_grad(v) = ||P v||_M^2 / ||v||_M^2

Também arquivamos 1 - f_grad, o quociente curl-curl da
componente projetada e o erro de representação da função gradiente
original, dado por:

    fit(c) = sqrt(1 - (c^T G c) / (c^T S c))

Nenhuma dessas medidas modifica os autovalores do EFG.

O corte tau=1e-8 é uma decisão numérica de diagnóstico,
**não** o corte de autovalores perto de zero do problema curl-curl.
A campanha inclui tau=1e-6 e tau=1e-10 na nuvem 17x15
para testar a sensibilidade de posto e indicadores.

## Arquivos e integridade

Por nuvem:

- eigenmode_gradient.csv: todos os autovetores e frações gradiente;
- selected_modes.csv: TM11, TM21, TM31, TM12 e o dubleto;
- gradient_fit_probes.csv: gradientes MLS de x, y, x²+y² e seno TM11;
- gradient_gram_eigenvalues.csv: espectro completo da Gramiana G;
- gradient_coupling.mtx, scalar_gradient_energy.mtx;
- projected_scalar_gram.mtx, gradient_projector_basis.mtx;
- vector_mass.mtx, eigenvectors.mtx (todas as colunas);
- metadata.txt e logs com convenções e limiares.

Matrizes em Matrix Market com 17 algarismos significativos.
A rotina Python de reabertura verifica G = H^T M^-1 H,
Q^T M Q, energia de projeção por autovetor, limites
de positividade e autovalores arquivados da campanha M6a.
Os dados históricos originais não são sobrescritos.

## Reprodução

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel 2
    /usr/bin/ctest --test-dir build --output-on-failure
    python3 scripts/m6d_gradient_campaign.py --binary build/gradient_kernel_audit --out out/m6d
    PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 python3 scripts/verify_m6d_gradient.py --campaign out/m6d

Os resultados são gravados em out/m6d/REPORT.md, summary.csv e
rank_sensitivity.csv.

No MyHoneyBee use PYTHONNOUSERSITE=1 para evitar o conflito
entre NumPy instalado no diretório do usuário e SciPy do Ubuntu.

## Limites da conclusão

Uma grande fração gradiente sugere associação ao quase-núcleo,
mas não é uma prova de que o autovetor seja espúrio.
A subtração dessas componentes **não foi aplicada** à formulação
e não deve ser realizada sem justificar uma estrutura compatível
e verificar contorno e divergência dos campos resultantes.
A formulação M5 rotacional/irrotacional continua pendente
de cotejo bibliográfico com a dissertação original.
