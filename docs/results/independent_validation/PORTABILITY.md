# Reabertura com implementações BLAS/LAPACK diferentes

A campanha `2026-10-21_dc1d98b` foi calculada com LAPACK ligado ao OpenBLAS.
O primeiro CI da PR usou a implementação BLAS/LAPACK de referência instalada
pelo Ubuntu. Compilação, testes de exportação, testes dos algoritmos e SHA-256
passaram, mas a reabertura falhou ao exigir igualdade quase exata do **erro
observado entre os dois solvers**, em
`solver_comparison.9x7.maximum_eigenvalue_absolute_difference`.

Esse erro observado é uma medida de arredondamento da execução arquivada,
e muda com a implementação BLAS, mesmo quando os novos autovalores, resíduos,
ortogonalidade e subespaços passam nos limites originais. Exigir igualdade
desse erro em outra implementação não testa a mesma propriedade que exigir
concordância dos autovalores dentro da tolerância científica declarada.

A correção verifica os resultados do solver adicional pelos limites já
publicados: tolerância absoluta `200*n*eps*max(1,max|lambda|)` mais termo
relativo `1e-9`; resíduos `<=1e-7`, erro de Gram em massa `<=1e-5` e correlação
mínima dos subespaços `>=1-1e-6`. Os erros derivados arquivados são conferidos
contra os próprios autovalores/diagnósticos arquivados, e os novos solves
continuam sendo verificados de forma independente. O nome da biblioteca
arquivada permanece como proveniência, sem exigir que todo ambiente futuro
use o mesmo nome de arquivo.

Os autovetores autoritativos Eigen e FEM estão arquivados; os autovetores do
solve LAPACK adicional sobre as matrizes EFG são recalculados pelo verificador.
Por isso, resíduos e correlações históricos desse solve adicional são
validados pelos limites e pela consistência interna, enquanto a proteção
contra qualquer mudança de bytes no histórico é feita separadamente pelo
manifesto SHA-256. A reabertura não reescreve medidas históricas para imitar
o backend novo.

As tabelas de comparação de campos, índices, dimensões, classificações e
geometria mantêm a verificação de reconstrução. Nenhum autovalor, autovetor,
matriz ou figura da campanha foi recalculado para corrigir esse teste.
Os limites físicos M4 permaneceram iguais, inclusive a reprovação de
TM22/TM41 por divergência `0,415285 > 0,25`.

Em Debian/Ubuntu com as bibliotecas instaladas, é possível reabrir a mesma
campanha usando os dois backends, sem alterar as alternativas do sistema:

```bash
CAMPAIGN="docs/results/independent_validation/$(cat docs/results/independent_validation/LATEST)"
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 \
  python3 scripts/independent_validation.py --campaign "$CAMPAIGN" --check-only
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 \
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3:/usr/lib/x86_64-linux-gnu/blas/libblas.so.3 \
  python3 scripts/independent_validation.py --campaign "$CAMPAIGN" --check-only
```

Essa segunda chamada serve apenas para verificar portabilidade do
pós-processamento. A execução FEM original e seus autovetores completos
continuam sendo os arquivos capturados da chamada real a `LAPACKE_dsygv`.

A correção passou na suíte local completa de 12 CTests e em 20 testes Python
dos algoritmos (incluindo seis novos casos de regressão de portabilidade),
além do teste Python do manifesto. A reabertura integral das 1896 linhas
passou tanto com OpenBLAS quanto com o preload de BLAS/LAPACK de referência.
Os 168 hashes da campanha permaneceram iguais.
