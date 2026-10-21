# Validação independente FEM × EFG

Campanhas novas são armazenadas em diretórios `AAAA-MM-DD_<commit>`. O arquivo
[`LATEST`](LATEST) identifica a campanha publicada; seu `REPORT.md` apresenta
os valores medidos, limites e veredictos. Resultados M4 anteriores foram preservados.
Consulte também a [auditoria dos solucionadores](SOLVER_AUDIT.md).

## Reproduzir os programas e todos os resultados

Dependências Ubuntu/Debian: `cmake`, `g++`, `libeigen3-dev`, `zlib1g-dev`,
`liblapacke-dev`, `python3-numpy`, `python3-scipy` e `python3-matplotlib`.
Alternativamente, instale os pacotes Python de `requirements-validation.txt`
num ambiente virtual; Eigen, zlib e LAPACKE continuam sendo dependências nativas.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
FEM_REPO=../fem-eigen-em-2d3d-whitney JOBS=2 \
  bash scripts/run_independent_validation.sh
```

O wrapper recusa sobrescrever campanhas existentes. Defina `CAMPAIGN_DIR` para
escolher um novo destino. Ele executa a campanha pública M4 com `JOBS=2
PROFILE=local`, compila o FEM original e executa `helm10_rect --ar-m 2.0
--nx NX --ny NY --nmodos 12 --backend closed-form` para 16×8, 30×15 e 40×20
células. Os DOFs livres TM são 105, 406 e 741. O parâmetro `--nmodos` limita
os campos nativos ilustrativos; o solucionador e a captura preservam o espectro
integral. O executável também resolve TE; somente a segunda chamada, TM com
Dirichlet, entra na análise. Os logs permitem conferir essa separação.

O registrador `fem_dsygv_capture.c` encaminha a chamada original a
`LAPACKE_dsygv`, sem substituir matrizes, resultados ou solver. Uma execução
adicional sem registrador verifica identidade dos arquivos FEM do caso grosso.
Os resultados nativos ficam recuperáveis no diretório temporário informado
pelo wrapper. Os arquivos necessários à análise são comprimidos e versionados.

`PYTHONNOUSERSITE=1` evita interferência de pacotes Python do usuário;
`OPENBLAS_NUM_THREADS=1` e `OMP_NUM_THREADS=1` controlam os threads.
`CTEST_BIN` pode selecionar uma instalação funcional do CTest; o wrapper usa
o executável que acompanha o `cmake` selecionado.

## Reabrir os resultados sem executar as discretizações

```bash
CAMPAIGN="docs/results/independent_validation/$(cat docs/results/independent_validation/LATEST)"
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 \
  python3 scripts/independent_validation.py --campaign "$CAMPAIGN" --check-only
python3 scripts/validation_manifest.py --campaign "$CAMPAIGN" --verify
```

Para gerar novamente tabelas e cinco figuras a partir das matrizes e campos
persistidos, omita `--check-only`. Essa etapa executa a verificação adicional
Eigen × LAPACK nas matrizes EFG salvas, sem montar novamente EFG ou FEM.
Depois de qualquer regeneração deliberada, recrie o manifesto com
`scripts/validation_manifest.py --campaign "$CAMPAIGN" --fem-repo CAMINHO_FEM`.

## Formatos e convenções

Cada diretório `efg/9x7`, `efg/13x11`, `efg/17x15` e `fem/coarse`,
`fem/medium`, `fem/fine` contém o espectro e os autovetores completos:

- `eigenvalues.csv`: índice iniciado em zero, λ assinado, kc, classe e resíduo.
- `eigenvectors.mtx.gz`: Matrix Market `array real general`, valores por coluna;
  a coluna matemática `i+1` pertence ao índice CSV `i`.
- `stiffness_matrix.mtx.gz`, `mass_matrix.mtx.gz`: matrizes originais; formato
  `coordinate` no EFG e `array` na captura FEM, ambos Matrix Market completo.
- `dof_map.csv`: EFG primeiro Hx, depois Hy, coordenadas e raio de cada nó;
  FEM um coeficiente Ez por nó interior, com valores de fronteira eliminados.
- `mode_diagnostics.csv` e `run_metadata.json`: medidas, parâmetros, commit,
  método numérico e convenções. Nenhum autovetor recebe arredondamento ou
  renormalização após sair do solver. A escrita autoritativa usa 17 algarismos.

No EFG, `grid.csv` e `field_curl.mtx.gz` guardam o rotacional MLS de todas as
colunas na grade comum 80×40, nos centros de células, em ordem y externa/x
interna. Os DOFs permitem reconstituir Hx, Hy e divergência por avaliação MLS.
No FEM, `mesh_nodes.csv` e `mesh_elements.csv` preservam malha/conectividade;
`fields_nodes.csv.gz` e `fields_elements.csv.gz` guardam Ez e seus gradientes
P1 para os primeiros 12 modos. Todos os outros campos podem ser reconstruídos
dos autovetores completos. `native_exports` preserva saídas do executável
original; sua precisão histórica é 16 algarismos, e a do VTK é menor.
`common_grid_ez.mtx.gz` e seu mapa de colunas guardam os candidatos usados
na comparação na mesma grade física do EFG.

Exemplo Python para reabrir as matrizes:

```python
import gzip
from scipy.io import mmread
with gzip.open("eigenvectors.mtx.gz", "rb") as f:
    V = mmread(f)  # V[:, i] corresponde à linha de índice i de eigenvalues.csv
```

A escala comum é uma unidade de comprimento, interpretada como 1 m para
comparar ao CLI FEM: a=2 m e b=1 m, λ em m⁻² e kc em m⁻¹. O EFG usa
coordenadas adimensionais; nenhum fator geométrico está implícito na comparação.
Em outra escala física L, λ deve ser dividido por L².

A comparação usa L2 físico ponderado, ajuste de sinal e subespaços para o
par TM22/TM41. Esses dois nomes identificam referências analíticas no mesmo
espaço de dimensão dois; as colunas discretas não recebem nomes individuais
arbitrários. Próximidade espectral, correlação espacial e condições físicas
de fronteira são verificados separadamente. Os modos quase nulos permanecem
nos arquivos, com λ negativo de arredondamento quando produzido pelo solver.

O manifesto enumera caminhos relativos, tamanhos e SHA-256 de todos os
arquivos da campanha, exceto ele próprio. Os commits dos dois programas e os
parâmetros de cada caso também constam nele. O arquivo não garante validade
física: garante que os dados auditados correspondem aos arquivos publicados.
