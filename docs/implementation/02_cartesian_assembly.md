# M2 — Montagem vetorial cartesiana EFG (protótipo)

**Status:** implementado na branch de desenvolvimento, sujeito a testes CTest/CI. O trabalho **não** é ainda uma reprodução completa dos experimentos de Al Maleh nem resolve o autoproblema de Maxwell.

## Incógnita e espaço vetorial

O campo aproximado é a dupla transversal de componentes de campo magnético

$$
\mathbf H_t^h=(\sum_i\phi_i h_{x,i},\;\sum_i\phi_i h_{y,i})^T.
$$

As funções de base são $\mathbf N_i^x=(\phi_i,0)^T$ e $\mathbf N_i^y=(0,\phi_i)^T$. Definimos o rotacional escalar 2D por $\mathrm{curl}_t\mathbf H_t=\partial_xH_y-\partial_yH_x$. Portanto:

$$
\mathbf K=\begin{bmatrix} K_{xx}&K_{xy}\\K_{yx}&K_{yy}\end{bmatrix},
\qquad
K_{xx}^{ij}=\int_\Omega\phi_{i,y}\phi_{j,y},\quad
K_{xy}^{ij}=-\int_\Omega\phi_{i,y}\phi_{j,x},
$$

$$
K_{yx}^{ij}=-\int_\Omega\phi_{i,x}\phi_{j,y},\quad
K_{yy}^{ij}=\int_\Omega\phi_{i,x}\phi_{j,x}.
$$

A matriz de massa tem blocos diagonais idênticos:

$$
\mathbf M=\begin{bmatrix}M_s&0\\0&M_s\end{bmatrix},
\qquad M_s^{ij}=\int_\Omega\phi_i\phi_j.
$$

A ordenação global de DOFs é **todas as componentes $H_x$ seguidas de todas as $H_y$**. O modelo inicial toma coeficientes materiais unitários e constantes.

## Geometria, quadratura e suporte

- Nuvem regular de nós `make_regular_cloud` e domínio retangular com limites explícitos.
- Quadratura gaussiana tensorial em células de fundo, ordens 2, 3, 4.
- **A malha de fundo serve apenas à integração; não produz funções de forma FEM.**
- Integração de funções MLS existentes de grau 1/2, com derivadas consistentes, por suporte circular de Wendland C².
- Montagem por `Eigen::SparseMatrix` e tripletos (didática; otimização futura).

## Fronteira e espectro

No modelo TM homogêneo em campo magnético transversal, a condição natural associada à forma fraca curl-curl é `curl_t H_t = 0` na fronteira PEC. Nesta **etapa de montagem**, nenhuma restrição nodal artificial é aplicada ao campo, e o operador discreto inclui um amplo núcleo de campos tipo gradiente. **Não classificar autovalores deste sistema como modos TM validados antes de implementar e analisar adequadamente a seleção dos modos positivos.**

Em particular, não penalizamos a matriz de massa e não impomos `H_t=0` na parede: ambas as decisões produziriam outro problema físico.

## Testes analíticos M2

Para $\Omega=[0,2]\times[0,1]$, a quadratura deve integrar área e momentos polinomiais $1,x,y,x^2,xy,y^2$ exatamente. Para campos reproduzidos pela base MLS linear/quadrática:

- Campo uniforme $\mathbf H=(1,0)$: massa $2$, rigidez zero.
- Campo uniforme $\mathbf H=(0,1)$: massa $2$, rigidez zero.
- Campo rotacional $\mathbf H=(-y,x)$: $\mathrm{curl}\mathbf H=2$, energia curl-curl $8$, massa $10/3$.
- Campo gradiente $\mathbf H=(x,y)$: rotacional zero e energia curl-curl zero.
- Simetria de $K$/$M$, bloco de massa cruzado zero e formas quadráticas não negativas em vetores de teste.

## Arquivos

- [API](../../include/al_maleh/cartesian.hpp)
- [Implementação](../../src/cartesian.cpp)
- [Testes](../../tests/test_cartesian.cpp)
- [Aplicativo de diagnóstico](../../apps/cartesian_rect.cpp)

Compilação:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/cartesian_rect 7 5 1 4.0
```

## Pendente para M3

Solver generalizado $Kq=\lambda Mq$, diagnóstico de posto e modos nulos, convergência dos autovalores TM positivos para casos PEC, projeção/divergência de campos e comparação independente com NASA TP-3485. A interpretação física de modos próximos de zero deve ser feita por critérios espectrais, não pela filtragem indiscriminada dos valores pequenos.
