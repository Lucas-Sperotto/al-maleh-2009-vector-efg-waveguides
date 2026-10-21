# M6e — Comparação independente FEM Whitney edge × EFG cartesiano

**Propósito:** confrontar o problema vetorial EFG curl-curl no retângulo
2×1, que contém raízes positivas adicionais abaixo do modo físico TM11,
com uma implementação FEM de **elementos de aresta Whitney 1-form**
já disponível no repositório independente NASA TP-3485.

**Esta etapa não substitui o EFG por FEM e não altera os autovalores M6a–M6d.**

## Fonte FEM independente, fixada por SHA

Repositório: https://github.com/Lucas-Sperotto/fem-eigen-em-2d3d-whitney

Commit: dc1d9545c1238554c2292aa0f64a3c1b883d4118

Arquivos existentes efetivamente usados:
- src/helmvec/main_edge_rect.cpp: executável edge_rect;
- src/edge/edge_dofs.cpp: arestas globais e orientação;
- src/edge/edge_basis.cpp: bases Whitney locais;
- src/edge/edge_assembly.cpp: formas curl-curl e massa;
- src/edge/mode_match_rect_edge_tm.hpp: campo transversal TM em H;
- src/core/spectral_csv.hpp: exportação de TODAS as colunas/autovalores.
- src/core/output_paths.hpp: TP3485_OUT_DIR escolhe saída isolada.

O FEM forma S_ij = ∫ curl Wi curl Wj e T_ij = ∫ Wi·Wj,
com diagonalização LAPACKE em modo TM PEC (contorno natural,
arestas de fronteira mantidas). O sistema TE utiliza outro
tratamento de fronteira; **não deve ser confundido com o TM**
físico em H_t implementado pelo EFG.

## Geometria e escala

O executável edge_rect original define a=1, b=0,5.
O EFG usa a=2, b=1. São retângulos similares com uma ampliação
isotrópica de fator 2, portanto, com propriedades materiais iguais:

    lambda_EFG = lambda_edge / 4
    kc_EFG = kc_edge / 2

Não se deve comparar diretamente os autovalores numéricos brutos.

A campanha M6e usa:

| EFG (nós) | FEM edge (células) | FEM edge (triângulos) |
|---|---:|---:|
| 9×7 | 8×4 | 64 |
| 13×11 | 12×6 | 144 |
| 17×15 | 16×8 | 256 |

Os graus de liberdade EFG e edge não são os mesmos; comparar
precisão diretamente apenas pelo número de elementos seria
inadequado. A campanha separa ambos e cita explicitamente os
parâmetros de malha.

## Produtos do experimento

Por nuvem ficam arquivados os artefatos NATIVOS do FEM edge,
inclusive espectro e autovetores de TE/TM (nenhuma coluna apagada),
S/T em CRS, campos nos centroides e logs; no projeto EFG são
gravados summary.csv, comparison.csv e REPORT.md.

O arquivo provenance.json registra o SHA de cada CSV original
e do espectro EFG M6a. O script verifica a revisão Git exata e a
ausência de edições locais nas fontes FEM relevantes.

Atribuímos preliminarmente modos pelos autovalores mais próximos
de TM11, TM21, TM31, TM12 e duas raízes no entorno de 5π².
A correspondência é **apenas espectral**: o dubleto TM22/TM41
pode sofrer rotações de base. O modo efetivamente físico depende
de correlação vetorial/ângulos principais sobre uma mesma
quadratura espacial; essa comparação é um trabalho posterior.

O relatório distingue:
- número de positivos abaixo da referência analítica TM11;
- número abaixo de 0,5×lambda(TM11), mais distante da referência;
- corte numérico explícito de classificação positiva;
- erros de autovalores EFG versus FEM edge versus analítico.

Não se classificam todos os valores abaixo da referência como
espúrios: uma aproximação numérica do próprio TM11 pode estar
ligeiramente abaixo de seu autovalor exato.

## Reproduzir com duas cópias dos repositórios

Na árvore EFG, com o segundo repositório FEM clonado em um
diretório irmão, construa o executável FEM:

    cmake -S ../fem-eigen-em-2d3d-whitney -B ../fem-eigen-em-2d3d-whitney/build -DCMAKE_BUILD_TYPE=Release
    cmake --build ../fem-eigen-em-2d3d-whitney/build --target edge_rect --parallel 2

**O script exige a revisão FEM pinned**. Para não mexer na
branch de trabalho do FEM, obtenha uma segunda cópia/worktree
fixada no commit dc1d9545..., compile-a e passe seus caminhos
nos argumentos --edge-source-dir/--edge-binary.

    python3 scripts/m6e_edge_compare.py \
        --edge-binary ../fem-eigen-em-2d3d-whitney/build/edge_rect \
        --edge-source-dir ../fem-eigen-em-2d3d-whitney \
        --out out/m6e

No GitHub Actions, o clone independente usa ref SHA explícito,
sparse checkout e compila exclusivamente o target edge_rect;
não faz commit no repositório FEM.

## Próxima etapa

Após checar estabilidade e precisão espectral, implementar
avaliação dos campos Whitney e dos campos EFG no **mesmo domínio
e pontos de quadratura**, com alinhamento de sinal e de escala
física de H. Para o dubleto 5π², comparar **subespaços de dimensão 2**
e ângulos principais em vez de colunas arbitrárias. Registrar
divergência, H·n, resíduos e efeitos de orientação dos DOFs.
Isso é necessário para afirmar mais do que concordância espectral.

Fontes M6a/M6b/M6c/M6d permanecem congeladas.
