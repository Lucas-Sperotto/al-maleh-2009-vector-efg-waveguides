# M6b — Formas de divergência e fronteira, métricas invariantes

**Estado:** nova instrumentação matemática de diagnóstico. **Não altera a
formulação K h = lambda M h**, não introduz penalização e não recalcula nem
modifica os arquivos M6a de 21/10/2026.

## Por que esta etapa existe?

Na campanha M6a, a dupla TM22/TM41 do retângulo PEC 2×1 foi identificada
no espaço de coeficientes EFG 17×15 por índices 192 e 193. Seus indicadores
individuais de divergência foram 0,007521 e 0,415285 (limite exploratório
0,25). Porém, qualquer base ortogonal diferente de um subespaço degenerado
pode produzir outras duas colunas. **O maior valor entre duas colunas não é
invariante sob rotações de base.** O valor que importa para uma afirmação
sobre todo o espaço é seu **supremo sobre todas as combinações não nulas**.

## Novas formas quadráticas

Para as bases vetoriais cartesianas do EFG, `N_i`, calculamos de modo
independente das matrizes M e K:

```text
D_ij = integral_Omega div(N_i) div(N_j) dOmega
B_ij = integral_boundary (n.N_i) (n.N_j) dGamma
M_ij = integral_Omega N_i.N_j dOmega  [já existente no M2]
```

A quadratura volumétrica D usa a mesma regra de integração de M/K.
A regra de fronteira B usa Gauss 2, 3 ou 4, aplicada **a cada lado de
cada célula de integração na fronteira**, para não subintegrar campos
com muitas oscilações. A normal externa vale (±1,0) e (0,±1).

Para uma matriz `V` com `k` colunas independentes, descrevendo o
subespaço físico `S = range(V)`, definimos:

```text
G = V^T M V  (Gram L2 de S)
D_S = V^T D V
B_S = V^T B V
```

Resolvemos dois autoproblemas simétricos-definidos **pequenos**:

```text
D_S c = mu G c
B_S c = nu G c
```

Os indicadores físico-numéricos extremos dentro de S são:

```text
div_min = sqrt(mu_min / lambda_ref)
div_max = sqrt(mu_max / lambda_ref)
normal_min = sqrt((area/perimeter)*nu_min)
normal_max = sqrt((area/perimeter)*nu_max)
```

Usa-se `lambda_ref=5*pi^2` para o dubleto TM22/TM41. Esta
normalização fixa, em vez do autovalor individual de cada componente,
é necessária para discutir todas as combinações do subespaço.

Sob uma troca de base não singular `V -> V R`, G, D_S e B_S sofrem
transformação por congruência; **os autovalores generalizados não mudam**.
A rotina retorna também a combinação `c` que maximiza cada indicador.
A rotina rejeita bases mal condicionadas/dependentes em produto M.

## Execução

Após compilar a branch M6b:

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure

python3 scripts/m6b_campaign.py --binary build/constraint_subspace --out out/m6b
~~~

O executável também pode ser chamado diretamente para reproduzir um
caso, por exemplo:

~~~bash
./build/constraint_subspace 17 15 192 193 out/m6b/17x15
~~~

Casos usados no experimento inicial:

| Nuvem | DOFs | Índices EFG do dubleto |
|---|---:|---|
| 9×7 | 126 | 70, 71 |
| 13×11 | 286 | 127, 128 |
| 17×15 | 510 | 192, 193 |

**Os índices vêm da campanha M6a e só podem ser reutilizados porque
a geometria, parâmetros e formulários K/M são os mesmos.** Se os
autovetores forem reordenados, esses índices deixam de representar os
mesmos modos. A comparação deve ser feita por projeção modal/subespaços,
não por índice cego.

O script gera `REPORT.md`, `summary.csv`, além de
`subspace_constraints.csv`, `positive_below_tm11.csv`,
`run_metadata.txt` e logs individuais. A CI guarda o artefato
`m6b-constraint-subspaces`.

## Testes

`test_constraints.cpp` verifica:

- D e B simétricas; energias de divergência e de fronteira;
- Campos constante, rotacional linear e gradiente linear, para MLS
  linear e quadrático e três regras de Gauss no contorno;
- Extremos calculados antes/depois de rotacionar a base;
- Correspondência da direção extremante com seu quociente de Rayleigh;
- Rejeição de subespaços dependentes e argumentos inválidos.

**O M6b ainda não resolve a origem dos autovalores espúrios.** A etapa
seguinte será calcular restrições compatíveis (ou projeções solenoidais),
investigar o núcleo gradiente e repetir a campanha integral de
validação física, preservando o baseline M6a.
