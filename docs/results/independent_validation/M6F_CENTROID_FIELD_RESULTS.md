# M6f — Comparação experimental de campos TM-H nos centroides de Whitney

Campanha executada em 21/10/2026 no [GitHub Actions run 39626802076](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39626802076) da branch experimental M6f. Nesta execução, **21/21 CTests passaram**, inclusive testes de sinal, escala e invariância de subespaços degenerados.

## Método: o que foi e não foi comparado

Foram confrontados os campos H=(Hx,Hy) do EFG cartesiano original e do FEM Whitney do executável edge_rect, **ramo TM-H**, nos mesmos centroides. O guia FEM foi reescalado de 1×0,5 para 2×1 por x'=2x e y'=2y; foram usados pesos de área física iguais por triângulo. Os campos foram renormalizados em L2 discreto, removendo a diferença entre as normalizações nativas EFG (massa) e FEM (pico unitário para exibição).

O EFG foi resolvido com as mesmas matrizes e parâmetros M6a e seus autovetores foram reconstruídos por MLS nos centroides; o FEM já disponibiliza os campos nativos nos centroides, sem alteração do código externo. A escolha inicial das colunas foi espectral via M6e, não pareamento global por forma.

**Isto é amostragem com um ponto por triângulo, NÃO uma quadratura convergida.** Não representa medição independente do erro de divergência ou das condições PEC.

## Correlação dos modos simples

A coluna rho representa a correlação absoluta na norma L2 amostrada sobre vetores TM-H completos; o erro associado, alinhando sinal, é sqrt(2−2*rho).

| Nuvem | rho(TM11) | rho(TM21) | rho(TM31) | rho(TM12) |
|---|---:|---:|---:|---:|
| 9×7 | 0,98880076 | 0,99121466 | 0,97767109 | 0,93391244 |
| 13×11 | 0,99502610 | 0,99613305 | 0,99063332 | 0,97212456 |
| 17×15 | **0,99720134** | **0,99784700** | **0,99484453** | **0,98455959** |

A concordância melhora nas três nuvens para os quatro modos testados. O modo TM12 segue com maior desajuste que TM11 e TM21. A coincidência espectral, por si, teria ocultado essa diferença de forma espacial.

Na malha fina, os erros L2 normalizados calculados foram: TM11 0,074815, TM21 0,065620, TM31 0,101543 e TM12 0,175729. São erros de **campos vetoriais em centroides**, não iguais aos erros anteriores de curl(H)×Ez no M6a.

## Dubleto degenerado TM22/TM41

Para cada discretização, formou-se uma base ortonormal L2 de dimensão 2 com os campos TM próximos da referência 5π². A SVD do produto cruzado mede ângulos principais invariantes à rotação interna dos vetores próprios.

| Nuvem | Menor cos(ângulo principal) EFG×Whitney | Maior ângulo (graus) | Menor cos EFG×analítico | Menor cos Whitney×analítico |
|---|---:|---:|---:|---:|
| 9×7 | 0,94339168 | 19,370823 | 0,99931025 | 0,94270874 |
| 13×11 | 0,97759909 | 12,150230 | 0,99985479 | 0,97768563 |
| 17×15 | **0,98305818** | **10,561673** | 0,99488196 | 0,98804999 |

O par EFG fino **permanece reprovado** na análise invariante de divergência de M6b (máximo 0,41532472 para o limite exploratório 0,25). Boa correlação por centroides não anula essa falha física. Também não permite concluir que Whitney ou EFG obedeçam exatamente às condições de contorno.

## Próximas condições de aceitação

- Reconstruir ambos os campos em quadratura Gauss comum e aumentar a resolução para testar convergência dos ângulos/correlações.
- Auditar divergência e H·n de campos vetoriais reconstruídos, inclusive a direção pior do dubleto.
- Confirmar correspondência modal por campo, não por proximidade de eigenvalores.
- Preservar espectros M6a–M6e, scripts, hashes e dados originais.

Protocolo: [10_whitney_efg_vector_overlap.md](../../implementation/10_whitney_efg_vector_overlap.md). A investigação principal M5 rotacional/irrotacional do autor segue na issue #11.
