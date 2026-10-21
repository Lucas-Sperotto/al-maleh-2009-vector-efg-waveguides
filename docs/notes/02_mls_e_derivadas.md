# MLS e derivadas

> Análise inicial migrada do antigo README (linhas 258–444), mantida para consulta e rastreabilidade. As equações são reconstruções de estudo, não transcrições verificadas de Al Maleh. Veja a [referência técnica](../theory/02_mls_efg_derivadas.md) e o [estado atual](../../README.md). Os delimitadores matemáticos foram adequados ao GitHub.

---

## 3. Construção das funções de forma pelo EFG

Antes de construir funções vetoriais, precisamos implementar corretamente as funções escalares EFG.

O método utiliza aproximações por mínimos quadrados móveis, ou MLS (*Moving Least Squares*).

Considere um conjunto de $N$ nós:

$$
\mathcal N=
\{\mathbf x_1,\mathbf x_2,\ldots,\mathbf x_N\}.
$$

A aproximação escalar é:

$$
u^h(\mathbf x)
=
\mathbf p^T(\mathbf x)\mathbf a(\mathbf x).
$$

Para uma base polinomial linear:

$$
\mathbf p(\mathbf x)
=
\begin{bmatrix}
1\\x\\y
\end{bmatrix}.
$$

Os coeficientes $\mathbf a(\mathbf x)$ são calculados minimizando:

$$
J(\mathbf a)
=
\sum_{i=1}^{N}
w_i(\mathbf x)
\left[
\mathbf p^T(\mathbf x_i)\mathbf a
-
u_i
\right]^2.
$$

A condição de mínimo fornece:

$$
\mathbf A(\mathbf x)\mathbf a(\mathbf x)
=
\mathbf B(\mathbf x)\mathbf u,
$$

com:

$$
\mathbf A(\mathbf x)
=
\sum_{i=1}^{N}
w_i(\mathbf x)
\mathbf p(\mathbf x_i)
\mathbf p^T(\mathbf x_i)
$$

e

$$
\mathbf B(\mathbf x)
=
\begin{bmatrix}
w_1\mathbf p(\mathbf x_1)&
\cdots&
w_N\mathbf p(\mathbf x_N)
\end{bmatrix}.
$$

Portanto:

$$
\mathbf a(\mathbf x)
=
\mathbf A^{-1}(\mathbf x)
\mathbf B(\mathbf x)\mathbf u.
$$

Substituindo na aproximação:

$$
u^h(\mathbf x)
=
\sum_{i=1}^{N}
\phi_i(\mathbf x)u_i,
$$

onde:

$$
\boxed{
\boldsymbol\phi^T(\mathbf x)
=
\mathbf p^T(\mathbf x)
\mathbf A^{-1}(\mathbf x)
\mathbf B(\mathbf x).
}
$$

Essas são as funções de forma MLS.

Na implementação, entretanto, não devemos calcular explicitamente a inversa de $\mathbf A$. É preferível resolver os sistemas lineares por fatoração.

### 3.1. Derivadas das funções MLS

Este é um dos pontos que considero mais importantes para a implementação.

As funções de forma dependem da posição tanto por meio da base polinomial quanto pelas funções peso.

Logo:

$$
\phi_i(\mathbf x)
=
\mathbf p^T(\mathbf x)\mathbf q_i(\mathbf x),
$$

com:

$$
\mathbf q_i=\mathbf A^{-1}\mathbf B_i.
$$

A derivada em relação a uma coordenada $\alpha$, com $\alpha=x$ ou $y$, é:

$$
\boxed{
\phi_{i,\alpha}
=
\mathbf p_{,\alpha}^{T}\mathbf q_i
+
\mathbf p^T\mathbf q_{i,\alpha}.
}
$$

Por sua vez:

$$
\boxed{
\mathbf q_{i,\alpha}
=
\mathbf A^{-1}
\left(
\mathbf B_{i,\alpha}
-
\mathbf A_{,\alpha}\mathbf q_i
\right).
}
$$

Observe que **não podemos considerar os coeficientes MLS constantes durante a diferenciação**.

Isso gera derivadas incorretas das funções de forma.

#### Uma descoberta no seu repositório existente

Examinei o arquivo `mls_shape.cpp` do seu repositório `simple-waveguide-efg-fem-solver`.

Atualmente, ele calcula as derivadas essencialmente por:

```cpp
result.dphidx[id] =
    (dpxdx.transpose() * coeff.col(j))(0);

result.dphidy[id] =
    (dpxdy.transpose() * coeff.col(j))(0);
```

Essas expressões não incluem as derivadas de $\mathbf A^{-1}\mathbf B$.

Portanto, **não representam as derivadas completas das funções MLS**.

Encontrei também outro problema na montagem: a penalização de contorno é adicionada tanto à matriz de rigidez quanto à matriz de massa. Isso altera o autoproblema e pode introduzir autovalores artificiais.

Arquivos examinados: [mls_shape.cpp](https://github.com/Lucas-Sperotto/simple-waveguide-efg-fem-solver/blob/main/mls_shape.cpp?utm_source=chatgpt.com) e [efgm_solver.cpp](https://github.com/Lucas-Sperotto/simple-waveguide-efg-fem-solver/blob/main/efgm_solver.cpp?utm_source=chatgpt.com).

Isso significa que podemos aproveitar a organização desse projeto, mas não devemos reutilizar seu núcleo MLS sem uma correção matemática e novos testes.

---
