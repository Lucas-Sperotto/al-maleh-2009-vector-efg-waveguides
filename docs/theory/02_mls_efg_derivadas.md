---
title: "02 — Aproximação MLS, método EFG e derivadas consistentes"
author: "Projeto de reprodução de Al Maleh (2009)"
date: "20 de outubro de 2026"
lang: pt-BR
---

# 1. Objetivo numérico

A estrutura MLS é a única dependência matemática comum às duas aproximações vetoriais planejadas. Erros nas primeiras derivadas contaminam a matriz curl-curl; erros nas segundas derivadas comprometem bases formadas a partir de gradientes. **Este texto apresenta uma formulação implementável, não afirma reproduzir a escolha exata de pesos ou suportes da dissertação antes de conferência no original.**

## 2. Aproximação por mínimos quadrados móveis

Considere pontos $\mathbf x_I=(x_I,y_I)$, coeficientes nodais generalizados $u_I$ e funções peso compactas $w_I(\mathbf x)$. Uma aproximação polinomial local é $u^h(\mathbf x)=\mathbf p^T(\mathbf x)\mathbf a(\mathbf x)$. Minimiza-se

$$
J(\mathbf a;\mathbf x)=\sum_{I\in\mathcal S(\mathbf x)} w_I(\mathbf x)
\bigl(\mathbf p^T(\mathbf x_I)\mathbf a(\mathbf x)-u_I\bigr)^2.
$$

Na reprodução inicial, usar a base linear $\mathbf p=(1,x,y)^T$; para avaliar segundas derivadas com precisão, experimentar também $\mathbf p=(1,x,y,x^2,xy,y^2)^T$. É preferível utilizar coordenadas centralizadas/escaladas localmente para reduzir condicionamento.

Definindo

$$
\mathbf A=\sum_I w_I\mathbf p_I\mathbf p_I^T,
\qquad \mathbf B_I=w_I\mathbf p_I,
\qquad \mathbf A\mathbf q_I=\mathbf B_I,
$$

obtemos as funções MLS

$$
\phi_I(\mathbf x)=\mathbf p^T(\mathbf x)\mathbf q_I(\mathbf x),
\qquad u^h(\mathbf x)=\sum_I \phi_I(\mathbf x)u_I.
$$

A soma só envolve os nós cujo suporte contém o ponto de quadratura. A matriz $\mathbf A$ precisa ter posto igual à dimensão da base; a contagem de vizinhos é necessária, mas não suficiente para garantir posto completo.

## 3. Derivadas de primeira ordem: regra do produto integral

Como $\mathbf A$ e $\mathbf B_I$ dependem de $\mathbf x$, diferenciar somente a base $\mathbf p$ está errado. Da identidade $\mathbf A\mathbf q_I=\mathbf B_I$ segue

$$
\mathbf A\mathbf q_{I,\alpha}=\mathbf B_{I,\alpha}-\mathbf A_{,\alpha}\mathbf q_I,
\quad\alpha\in\{x,y\}.
$$

Então

$$
\boxed{\phi_{I,\alpha}=\mathbf p_{,\alpha}^T\mathbf q_I+\mathbf p^T\mathbf q_{I,\alpha}.}
$$

As derivadas de $\mathbf A$ podem ser montadas de modo exato a partir das derivadas da função peso:

$$
\mathbf A_{,\alpha}=\sum_I w_{I,\alpha}\mathbf p_I\mathbf p_I^T,
\qquad\mathbf B_{I,\alpha}=w_{I,\alpha}\mathbf p_I,
$$

caso $\mathbf p_I$ esteja na base física fixa. Se a base estiver centrada no ponto de avaliação, incluir também as derivadas explícitas da mudança de coordenadas; não misturar as duas convenções.

## 4. Derivadas de segunda ordem

Diferenciando novamente a identidade do sistema MLS, para $\alpha,\beta\in\{x,y\}$:

$$
\mathbf A\mathbf q_{I,\alpha\beta}=
\mathbf B_{I,\alpha\beta}-\mathbf A_{,\alpha\beta}\mathbf q_I
-\mathbf A_{,\alpha}\mathbf q_{I,\beta}-\mathbf A_{,\beta}\mathbf q_{I,\alpha}.
$$

A função de forma tem derivada

$$
\boxed{\phi_{I,\alpha\beta}=\mathbf p_{,\alpha\beta}^T\mathbf q_I
+\mathbf p_{,\alpha}^T\mathbf q_{I,\beta}
+\mathbf p_{,\beta}^T\mathbf q_{I,\alpha}
+\mathbf p^T\mathbf q_{I,\alpha\beta}.}
$$

Para bases rotacional/irrotacional formadas por $\nabla\phi_I$, o laplaciano $\Delta\phi_I=\phi_{I,xx}+\phi_{I,yy}$ entra diretamente no rotacional. As funções peso precisam permitir diferenciação de segunda ordem e continuidade compatível com a integração escolhida.

## 5. Pesos, suporte e suavidade

Uma opção experimental é o peso de Wendland $C^2$ em coordenada radial $q=r/h$, para $0\le q<1$:

$$
w(q)=(1-q)^4(4q+1),\qquad w(q)=0\text{ para }q\ge1.
$$

Essa escolha **não é atribuída a Al Maleh**. O peso, escala $h$, forma do suporte, regularidade e base devem ser expostos na configuração e registrados nos metadados de cada execução.

No projeto legado `simple-waveguide-efg-fem-solver`, o peso é implementado por $e^{-6{,}25q^2}$ para $q<1$ e zero abruptamente fora do suporte. Há um salto em $q=1$ ($e^{-6{,}25}\neq0$); esse peso truncado não é sequer contínuo na fronteira do suporte. Não reutilizar essa função na versão que requer derivadas suaves sem corrigir a construção.

## 6. Estabilidade da solução linear

Nunca montar $\mathbf A^{-1}$ explicitamente. Resolver $\mathbf A\mathbf q=\mathbf B$ por LDLT/Cholesky quando a matriz for simétrica definida positiva, ou QR com pivoteamento/SVD nos diagnósticos. Registrar $\kappa(\mathbf A)$ e a menor singularidade por quadratura. A regularização deve ser uma decisão controlada e comparada com estratégia adaptativa de suportes.

## 7. Testes que bloqueiam o assembly

Para cada ponto de quadratura, verificar:

$$
\sum_I\phi_I=1,\qquad \sum_I\nabla\phi_I=\mathbf0,
\qquad\sum_I\phi_Ix_I=x,\qquad\sum_I\phi_Iy_I=y.
$$

Para base linear, também:

$$
\sum_I x_I\phi_{I,x}=1,\quad \sum_I y_I\phi_{I,y}=1,
\quad \sum_I x_I\phi_{I,y}=0,\quad\sum_I y_I\phi_{I,x}=0.
$$

Se a implementação usa segundas derivadas, verificar $\sum_I\phi_{I,xx}=0$, $\sum_I\phi_{I,xy}=0$ e $\sum_I\phi_{I,yy}=0$ e a reprodução de polinômios quadráticos quando essa base é adotada. Comparar gradientes/Hessianas analíticas com diferenciação numérica em pontos internos onde o conjunto de nós ativos permaneça constante.

## 8. Contrato recomendado de API

```cpp
struct MLSResult {
    std::vector<int> active_nodes;
    Vector phi, dphidx, dphidy;
    Vector d2phidx2, d2phidxdy, d2phidy2;
    double condition_A;
};

MLSResult evaluate_mls(
    const Point2& x,
    const NodeCloud& cloud,
    const MLSConfig& config);
```

`MLSResult` deve informar erro explícito se o posto ou a continuidade necessária falhar; não devolver um vetor de zeros que será silenciosamente interpretado como resultado físico.

## 9. Auditoria de código legado

Foi verificado que `simple-waveguide-efg-fem-solver/mls_shape.cpp` calcula gradientes usando apenas a derivada de $\mathbf p$, omitindo o termo $\mathbf p^T\mathbf q_{I,\alpha}$. Portanto esse arquivo serve como referência histórica, **não como código científico a copiar sem correção**. O arquivo `efgm_solver.cpp` também adiciona uma penalização aos termos diagonais de **ambas** as matrizes $K$ e $M$ para nós marcados como fronteira; isso altera o espectro e deve ser reavaliado. O tratamento correto de fronteira não decorre da marcação nodal, pois MLS não apresenta em geral a propriedade de Kronecker.

## Referências

Belytschko, T.; Lu, Y. Y.; Gu, L. Element-free Galerkin methods. *International Journal for Numerical Methods in Engineering*, 37, 229–256, 1994. DOI 10.1002/nme.1620370205.

Al Maleh, H. *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems*. McGill University, 2009. A escolha exata de peso da dissertação deve ser conferida no texto primário.