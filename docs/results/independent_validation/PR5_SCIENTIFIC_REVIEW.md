# Auditoria científica complementar da PR #5 (21/10/2026)

**Estado:** revisão adicional do espectro arquivado; **não altera** autovalores, autovetores, matrizes, figuras, critérios originais ou manifesto SHA-256 da campanha `2026-10-21_dc1d98b`.

## Conclusões verificadas

1. **Concordância algébrica forte:** Eigen e LAPACK concordam sobre as mesmas matrizes K/M EFG, com diferença máxima em lambda de aproximadamente 1,03e-8. Isso valida a diagonalização, **não a física de todos os pares**.
2. **Quatro modos TM simples sustentados** na comparação independente FEM × EFG, com boa concordância de campos e erros analíticos pequenos. Os arquivos originais continuam sendo a referência autoritativa.
3. **Falha física do subespaço TM22/TM41:** os autovetores arquivados de índices 192/193 apresentam razões de divergência 0,007521 / 0,415285, respectivamente. O limite exploratório é 0,25. A direção 193 pertence ao subespaço de dois vetores; logo o **supremo da razão de divergência nesse subespaço** deve ser pelo menos tão grande quanto o valor dessa direção (após a normalização comum por sqrt(lambda_T M) e norma L2), independentemente de como uma futura análise escolha a base. Um valor calculado apenas como o máximo de duas colunas *não é ele próprio invariante sob rotações*: é somente **limite inferior** para esse supremo.
4. **Contaminação do espectro de baixa frequência:** para os refinamentos EFG de 126/286/510 graus de liberdade, há 49/94/144 entradas positivas `lambda>zero_cutoff` abaixo da primeira referência TM `lambda_11=5*pi^2/4`. Nenhuma dessas entradas é um modo TM positivo do problema retangular contínuo de Dirichlet; suas propriedades físicas individuais, e eventual natureza de quase-núcleo, necessitam estudo. A contagem depende do corte numérico, logo não mede a dimensão exata de um núcleo.
5. **Convergência física ainda incompleta:** campo curl(H) correlacionado com Ez e pequeno resíduo matricial **não** garantem campo H solenoidal e condições de contorno; a coluna 193 ilustra isso. A grade comum 80×40 de comparação dos campos ainda não teve estudo de refinamento separado.

## Interpretação matemática

No operador contínuo `curl-curl` em 2D com condição natural `curl H = 0` no contorno, qualquer autovetor de autovalor positivo obedece, formalmente, `div H=0`, enquanto `s=curl H` satisfaz `-Delta s=lambda s` com `s=0` nas paredes. Logo o espectro positivo ideal corresponde ao escalar TM Dirichlet. O operador contínuo também tem um núcleo gradiente grande. O espaço MLS cartesiano de dois blocos **não foi demonstrado** compatível com uma sequência discreta gradiente-rotacional, portanto valores próprios pequenos não nulos e modos sem convergência física são possíveis.

Isso motiva validar uma formulação discretamente solenoidal/compatível, não modificar artificialmente a tolerância de classificação para ocultar candidatos.

## Como reabrir esta evidência

~~~bash
python3 scripts/audit_spectral_pollution.py \
  --campaign docs/results/independent_validation/$(cat docs/results/independent_validation/LATEST)
~~~

O programa utiliza apenas bibliotecas padrão Python e os CSV/JSON publicados. Confere ordenação, índices, limiar, contagens, distribuição das raízes de baixa frequência e um **limite inferior** da pior divergência no subespaço degenerado. Para salvar relatório JSON separadamente, utilize `--output` apontando para outro arquivo. O script não altera resultados arquivados por padrão.

Há três testes de regressão em `tests/test_spectral_pollution.py`, integrados ao CTest, que verificam as contagens de 49/94/144, a reprovação do par e entradas inválidas.

## Próxima campanha recomendada (M6b)

- Construir as formas discretas **independentes** de divergência e fluxo normal: `D_ij=int(div Ni)(div Nj)`, `B_ij=int_boundary(n.Ni)(n.Nj)`, com testes de quadratura e invariância.
- Diagonalizar **dentro do subespaço de duas dimensões** as formas de divergência e de parede com gramiana L2, produzindo indicadores máximo e mínimo invariantes a mudanças de base. Comparar com os limites inferiores desta auditoria.
- Examinar projeção de campos quase irrotacionais/gradientes e estratégias compatíveis de restrição ou formulação mista; evitar penalização ad hoc que apenas melhore a tabela.
- Repetir os três refinamentos EFG e adicionar estudo isolado da grade física comum (por exemplo, 40×20, 80×40 e 160×80) sobre os **mesmos** autovetores persistidos.
- Preservar a primeira campanha como baseline; qualquer mudança física na formulação deve gerar **nova campanha e novo manifesto**, com comparação antes/depois e sem substituir a reprovação histórica.

**Importante:** Não considerar a validação física como completa até que a contaminação de baixa frequência e o par degenerado sejam investigados.
