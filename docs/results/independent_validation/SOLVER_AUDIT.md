# Auditoria dos sistemas espectrais EFG e FEM

Auditoria em 21 de outubro de 2026. Esta página registra os programas, matrizes,
solucionadores e convenções. As tabelas de comparação física da campanha ficam
no `REPORT.md` e nos CSVs da campanha identificada por data e commit.

## Proveniência e preservação dos projetos

| Projeto | Caminho local auditado | Origem científica inicial |
|---|---|---|
| EFG cartesiano | `/home/sperotto/al-maleh-2009-vector-efg-waveguides` | `33fbf75` |
| FEM TP-3485 | `/home/sperotto/fem-eigen-em-2d3d-whitney` | `4c56a526b1a308edd2e17e3beec306f76dd5dce5` |

Os remotos são, respectivamente,
`https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides.git` e
`https://github.com/Lucas-Sperotto/fem-eigen-em-2d3d-whitney.git`.
O commit exato do executável EFG ampliado e os parâmetros de cada execução
constam em `run_metadata.json` e no manifesto da campanha; o hash inicial
acima não é apresentado como hash do novo exportador.

O FEM estava na branch `main`, com um único worktree no caminho indicado.
Havia saídas locais não rastreadas em `out/2d/`, `out/RESULTS_REPORT.md`,
`out/img_all/`, `out/run_all.log` e `out/validation/`. Elas foram preservadas.
A compilação e a execução desta campanha usaram diretórios separados e
`TP3485_OUT_DIR`. O código FEM permaneceu intacto; não há commit FEM desta
validação. A branch EFG de trabalho é
`feat/m6a-independent-fem-validation`.

## Sistema EFG

A montagem em `src/cartesian.cpp` integra

\[
 K_{ij}=\int_\Omega \operatorname{curl}N_i\operatorname{curl}N_j\,dA,
 \qquad M_{ij}=\int_\Omega N_i\cdot N_j\,dA.
\]

As funções vetoriais são `[phi_i,0]` e `[0,phi_i]`, com `phi_i` MLS.
A ordem dos coeficientes é **todos Hx, depois todos Hy**. O rotacional
reconstruído é `dHy/dx - dHx/dy`. A campanha usa domínio 2 × 1, MLS linear,
fator de suporte 3,3 e quadratura 3 × 3 por célula. Os espaços têm
126, 286 e 510 DOFs para nuvens 9 × 7, 13 × 11 e 17 × 15.

A montagem M3/M4 mantém todos os coeficientes de contorno. Não há eliminação
forte do fluxo normal nem penalização da divergência. A condição natural é
`curl(H)=0` no contorno. O operador curl-curl tem um núcleo de campos com
rotacional nulo, que deve ser distinguido de modos TM positivos; autovalores
positivos, isoladamente, não demonstram ausência de campos espúrios.

`src/spectral.cpp` valida dimensões, finitude e simetria, usando tolerância
relativa padrão `1e-10` para as normas de Frobenius. O sistema efetivamente
entregue ao Eigen usa `A=(K+K^T)/2` e `B=(M+M^T)/2`. O exportador guarda
as matrizes montadas originais e seus indicadores de assimetria. O espectro
de B deve ser estritamente positivo; o limite padrão de condicionamento é
`1e13`. O limite do solucionador denso de referência é 600 DOFs.

A chamada é `Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd>(A,B)`.
A instalação auditada é Eigen 3.4.0, em aritmética `double` IEEE 754.
Os cabeçalhos locais
`Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h` e
`SelfAdjointEigenSolver.h` mostram Cholesky `B=LL^T`, redução
`L^-1 A L^-T`, tridiagonalização e QR simétrico implícito com deslocamento
de Wilkinson, seguidos da transformação inversa dos autovetores.
As colunas são preservadas com `V^T M V ≈ I`.

O limiar de núcleo é
`max(zero_absolute, zero_relative * max(lambda_max,0))`, com ambos os
parâmetros padrão iguais a `1e-9`. Um autovalor inferior ao negativo desse
limiar causa erro; pequenos valores negativos compatíveis com arredondamento
são mantidos. Todos os valores `lambda <= zero_cutoff` são classificados
como `near_zero`, sem descartar colunas na persistência. A exportação mantém
o sinal de lambda; `kc=sqrt(max(0,lambda))` é apenas uma coluna derivada.

## Sistema FEM independente

`helm10_rect` gera dois triângulos P1 por célula, com diagonal
`n00 -> n11`. A malha é numerada por linhas, com x variando mais rapidamente.
A família escolhida para esta validação é **TM escalar com Dirichlet**:
`Ez=0` em todos os nós da fronteira, eliminados antes da solução. Os DOFs
livres seguem a ordem dos nós interiores. O comando também resolve TE com
Neumann, mas seus autovetores não entram na comparação TM.

As matrizes independentes, montadas por
`src/core/helm10_scalar_system.cpp` com o backend `closed-form`, são

\[
 S_{ij}=\int_\Omega \nabla N_i\cdot\nabla N_j\,dA,
 \qquad T_{ij}=\int_\Omega N_i N_j\,dA,
 \qquad S u=\lambda T u.
\]

Para um triângulo de área A, `S_e(i,j)=(b_i b_j+c_i c_j)/(4A)` e
`T_e=A/12 * [[2,1,1],[1,2,1],[1,1,2]]`, com `eps_r=mu_r=1`.
S e T são simétricas; T é definida positiva no espaço interior. A matriz S
TM também é definida positiva, ao contrário do problema TE, cujo modo
constante é nulo. O pós-processamento nativo TE remove valores menores que
`1e-9` ao listar modos físicos. A exportação espectral nativa contém todos
os autovalores e todas as colunas, independentemente desse filtro de exibição.
O empacotamento FEM classifica `abs(lambda) <= 1e-10*max(1,max(abs(lambda)))`
como `near_zero` e valores abaixo do negativo desse limiar como `negative`;
nenhuma dessas classes é removida dos arquivos. Nas três soluções TM todos
os valores são estritamente positivos. Esse limiar de diagnóstico FEM não
substitui o limiar configurável de núcleo EFG.

`src/core/lapack_eig.hpp` copia S/T para vetores coluna-major e chama

```cpp
LAPACKE_dsygv(LAPACK_COL_MAJOR, 1, 'V', 'U', n,
              A.data(), n, B.data(), n, w.data());
```

`itype=1` representa `S u=lambda T u`; `jobz=V` pede todos os autovetores;
`uplo=U` usa o triângulo superior. `info != 0` é erro fatal. A implementação
instalada, LAPACK/LAPACKE 3.12.0 com OpenBLAS 0.3.26, opera em binary64.
A inspeção de `dsygv_` e `dsyev_` na biblioteca local mostrou as chamadas
`DPOTRF`, `DSYGST`, `DSYEV`, `DSYTRD`, `DORGTR`, `DSTEQR` e `DTRSM`:
Cholesky, redução simétrica, tridiagonalização, iteração QL/QR e transformação
inversa. A opção usada fornece `Z^T T Z=I`. O programa não substitui esse
solve por Eigen ou por valores analíticos.

## Execuções reais e prova de captura transparente

Compilação isolada executada:

```bash
cmake -S /home/sperotto/fem-eigen-em-2d3d-whitney \
  -B /tmp/m6a-fem/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/m6a-fem/build --target helm10_rect -j2
cc -O2 -shared -fPIC scripts/fem_dsygv_capture.c \
  -o /tmp/m6a-fem/capture_lapack.so -ldl -lm
```

Cada rodada usou `OPENBLAS_NUM_THREADS=1`, `TP3485_OUT_DIR`,
`M6A_LAPACK_CAPTURE_DIR` e `LD_PRELOAD`, com a CLI:

```bash
helm10_rect --ar-m 2.0 --nx NX --ny NY --nmodos 12 --backend closed-form
```

O parâmetro `--ar-m 2.0` produz `a=2 m`, `b=1 m`. O EFG usa as mesmas
coordenadas numéricas, interpretadas aqui com **uma unidade de comprimento
igual a um metro**. Assim, lambda tem unidade `m^-2` e `kc` tem `m^-1`.
Nenhuma escala geométrica automática é aplicada para melhorar concordância.

O interceptor de `scripts/fem_dsygv_capture.c` resolve a função real por
`dlsym(RTLD_NEXT,...)` e encaminha os argumentos originais. Ele registra
S/T antes da chamada, os autovetores/autovalores depois dela e o retorno
`info`; não altera as matrizes ou os resultados. Os logs guardados mostram
as chamadas TE e TM em cada execução, com `layout=102`, `itype=1`,
`jobz=V`, `uplo=U` e `info=0`. A segunda chamada é o sistema TM usado.

| Rodada | nx × ny | Nós | Triângulos | DOFs TM / autovalores | DOFs TE registrados |
|---|---:|---:|---:|---:|---:|
| coarse | 16 × 8 | 153 | 256 | 105 | 153 |
| medium | 30 × 15 | 496 | 900 | 406 | 496 |
| fine | 40 × 20 | 861 | 1600 | 741 | 861 |

Um controle coarse adicional executou o mesmo binário, com os mesmos
argumentos e número de threads, sem `LD_PRELOAD`. As quatro saídas nativas
TM — autovalores, autovetores, S e T em CSV — são idênticas byte a byte.
Os hashes estão no metadado coarse. Para os autovalores, SHA-256:
`4b13453f025addf2c73f13c5516a4fac75bce07d3d7d98b9ea409b51c9b8b005`;
para os autovetores:
`5933947b45db874ccc7c4d925d3a1bc7505e050abd37d880894af057009a081a`.
Essa verificação controla a transparência da instrumentação; não é um
segundo solucionador numérico.

## Persistência e verificações numéricas

Os CSVs nativos FEM usam 16 algarismos significativos. A captura publica
`eigenvectors.mtx.gz`, `mass_matrix.mtx.gz` e `stiffness_matrix.mtx.gz`
com 17 algarismos significativos, no formato Matrix Market `array real
general`; todas as colunas TM são preservadas. O índice zero de
`eigenvalues.csv` corresponde à primeira coluna matemática do arquivo
Matrix Market. O gzip é determinístico, com `mtime=0`.

`dof_map.csv`, `mesh_nodes.csv` e `mesh_elements.csv` permitem reconstruir
Ez por interpolação P1. O VTK nativo do primeiro modo é preservado em
`native_exports/tm11_rect_sv.vtk.gz`. A conectividade VTK coincide
exatamente com a conectividade exportada. As coordenadas do VTK nativo
medium têm arredondamento de exibição de até `3.334e-6 m`; as coordenadas
precisas do mapa seguem a fórmula da malha e são verificadas pela remontagem
de S/T. O erro máximo dessa remontagem foi `3.553e-15` em S e `4.337e-18`
em T. A diferença entre S/T capturadas e os CSVs nativos de 16 dígitos
foi no máximo `4.441e-16` e `4.337e-19`, respectivamente.

Além dos coeficientes completos, `fields_nodes.csv.gz` salva Ez dos
primeiros 12 modos TM, sem mudar sua normalização em T.
`fields_elements.csv.gz` salva seus gradientes P1 exatos por triângulo.
Os campos e CSVs nativos também são preservados comprimidos. O campo
transversal didático nativo FEM usa `Ex=-dEz/dx`, `Ey=-dEz/dy`; metadados
nativos de frequência, beta e Ztm não devem ser confundidos com a
normalização usada na comparação longitudinal.

Para tornar explícitas as escalas dos resíduos, a rotina EFG original usa

\[
 r_F=\frac{\|Kv-\lambda Mv\|_2}
 {\max(10^{-300},(\|K\|_F+|\lambda|\|M\|_F)\|v\|_2)}.
\]

O empacotamento FEM registra um resíduo de ação das matrizes,

\[
 r_A=\frac{\|Sv-\lambda Tv\|_2}
 {\max(\mathrm{tiny},\|Sv\|_2+|\lambda|\|Tv\|_2)}.
\]

As duas definições têm denominadores diferentes e seus números não devem
ser comparados como se usassem a mesma escala. A comparação entre
solucionadores sobre as matrizes EFG deve ainda tratar os subespaços, pois
um núcleo ou um autovalor degenerado pode mudar sua base sem mudar a solução.
Precisamente, o numerador calculado pelo EFG usa as matrizes simetrizadas
A/B; as normas do denominador usam K/M montadas. A fórmula com K/M acima
pressupõe simetria exata. A leitura dos arquivos persistidos verifica também
o resíduo nas matrizes originais exportadas.

| Rodada FEM | max r_A no espectro completo | max abs(ZᵀTZ−I) | Menor autovalor de T |
|---|---:|---:|---:|
| coarse | 4.333e-14 | 4.441e-15 | 4.13924e-3 |
| medium | 3.348e-13 | 7.994e-15 | 1.13050e-3 |
| fine | 6.453e-13 | 1.066e-14 | 6.31164e-4 |

As 14 verificações FEM por rodada passaram: finitude, dimensões completas,
simetrias S/T, positividade de T, resíduos, ortogonalidade em T, remontagem
S/T pela malha, positividade do espectro TM, conectividade e coordenadas do
VTK nativo e concordância de S/T com os exports nativos. São 42 verificações
aprovadas e nenhuma falha neste empacotamento. Nenhum autovalor TM foi
classificado como quase nulo ou negativo; esse resultado não afirma que o
núcleo EFG seja vazio.

O fato de um par resolver bem seu sistema matricial demonstra precisão
algébrica. A validação física exige, separadamente, correlação de Ez FEM com
curl(H) EFG na mesma malha física, erros de campo após normalização e ajuste
de sinal e indicadores EFG de divergência e fluxo normal. Para TM22/TM41,
cujo autovalor exato é `5*pi^2`, a unidade de comparação física é o subespaço
de dimensão dois; a orientação de uma coluna isolada é arbitrária.

A preparação dos arquivos, incluindo essas verificações, é reproduzida por
`scripts/prepare_fem_exports.py`. O comando completo de campanha é fornecido
em `scripts/run_independent_validation.sh`. Os logs das execuções e os
metadados versionados são a evidência dos resultados; as pastas `/tmp` acima
descrevem somente o ambiente da execução original.
