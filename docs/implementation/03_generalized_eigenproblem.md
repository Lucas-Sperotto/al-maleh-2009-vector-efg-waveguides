# M3 — Autoproblema generalizado, diagnóstico e referências TM

**Estado:** protótipo espectral numérico para o sistema vetorial cartesiano montado no M2. **Não é ainda uma reprodução verificada das tabelas da dissertação de Al Maleh.**

## Sistema de autovalores

Para o campo magnético transversal aproximado por funções nodais MLS, o sistema discretizado é

$$
K\,h = \lambda\,M\,h,\qquad \lambda=k_c^2.
$$

O solucionador atual usa Eigen, com matrizes densas **somente para até 600 graus de liberdade**. Recebe as matrizes esparsas simétricas de M2 e verifica dimensões, finitude, simetria, positividade definida da matriz de massa e seu condicionamento. Caso a massa esteja singular ou muito mal condicionada, a execução **falha explicitamente** em vez de calcular autovalores enganadores.

O Eigen GeneralizedSelfAdjointEigenSolver retorna todos os autovalores (inclusive valores próximos de zero) e autovetores normalizados na massa. O programa registra, por modo, o resíduo relativo da equação matricial e a tolerância utilizada para classificar o núcleo numérico. A classificação numérica não é prova de que um modo seja físico.

## Referência TM: guia retangular PEC

Para um retângulo homogêneo com lados a e b, a referência escalar analítica dos modos TM é

$$
\lambda_{mn}=\pi^2\left[\left(\frac{m}{a}\right)^2+
\left(\frac{n}{b}\right)^2\right],\quad m,n\ge 1.
$$

Para a=2 e b=1, TM11 vale 1.25 pi²; TM21 vale 2 pi². O programa compara cada referência ao autovalor positivo numericamente mais próximo, sem repetir índice. **Isso é apenas uma associação preliminar de valores**: modos espúrios e degenerações podem invalidar a interpretação.

## Critérios do campo

O pós-processamento avalia campos por funções MLS, seus rotacionais e suas divergências. Registra o indicador adimensional:

$$
\eta_\mathrm{div} =
\frac{\|\nabla_t\cdot H_t\|_{L^2(\Omega)}}
{\sqrt{\lambda}\,\|H_t\|_{L^2(\Omega)}}.
$$

Também avalia um indicador integrado de fluxo normal sobre as quatro paredes PEC. Para um modo TM convergente, esperamos redução da divergência e do fluxo normal sob refinamento. Os autovalores positivos não são automaticamente modos TM.

## Arquivos e execução

Código-fonte: [src/spectral.cpp](../../src/spectral.cpp), API em [include/al_maleh/spectral.hpp](../../include/al_maleh/spectral.hpp), testes em [tests/test_spectral.cpp](../../tests/test_spectral.cpp), aplicativo em [apps/cartesian_spectrum.cpp](../../apps/cartesian_spectrum.cpp).

Comandos:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/cartesian_spectrum 5 5 1 3.3 out/m3/rect
```

Saídas: eigenvalues.csv (espectro completo e resíduos), tm_reference.csv (associação analítica), mode_diagnostics.csv (métricas dos primeiros candidatos positivos), field_candidate_1.csv até field_candidate_4.csv (amostras de campos) e run_metadata.txt (parâmetros/limiares). A pasta out/ é ignorada pelo Git.

## Limitações e M4

- M3 não incorpora projeção explícita no subespaço solenoidal nem uma técnica de eliminação de modos espúrios.
- O cálculo do espectro é denso, adequado a problemas didáticos pequenos.
- Valores de convergência física, sensibilidade a raio de influência/ordem/quadratura, correspondência de autovetores e fronteira PEC exigem campanhas adicionais e comparação independente com FEM.
- A revisão dos demais documentos teóricos continua pendente, conforme solicitado pelo projeto.
