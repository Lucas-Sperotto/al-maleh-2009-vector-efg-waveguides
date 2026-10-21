# Formulações vetoriais

> Análise inicial migrada do antigo README (linhas 445–912), mantida para consulta e rastreabilidade. As equações são reconstruções de estudo, não transcrições verificadas de Al Maleh. Veja a [referência técnica](../theory/03_formulacoes_vetoriais.md) e o [estado atual](../../README.md). Os delimitadores matemáticos foram adequados ao GitHub.

---

## 4. Primeira formulação: aproximação pelas componentes cartesianas

Essa é a primeira abordagem descrita por Al Maleh.

O campo magnético é aproximado separadamente em suas componentes:

$$
H_x^h(x,y)
=
\sum_{i=1}^{N}
\phi_i(x,y)h_{xi}
$$

e

$$
H_y^h(x,y)
=
\sum_{i=1}^{N}
\phi_i(x,y)h_{yi}.
$$

Assim:

$$
\boxed{
\mathbf H_t^h
=
\sum_{i=1}^{N}
\begin{bmatrix}
\phi_i&0\\
0&\phi_i
\end{bmatrix}
\begin{bmatrix}
h_{xi}\\
h_{yi}
\end{bmatrix}.
}
$$

Temos $2N$ incógnitas.

Podemos definir duas funções vetoriais para cada nó:

$$
\mathbf N_i^x=
\begin{bmatrix}
\phi_i\\
0
\end{bmatrix}
$$

e

$$
\mathbf N_i^y=
\begin{bmatrix}
0\\
\phi_i
\end{bmatrix}.
$$

O rotacional dessas funções é:

$$
\operatorname{curl}_t\mathbf N_i^x
=
-\frac{\partial\phi_i}{\partial y}
$$

e

$$
\operatorname{curl}_t\mathbf N_i^y
=
\frac{\partial\phi_i}{\partial x}.
$$

Agora podemos construir as matrizes globais.

### 4.1. Matriz de rigidez

A matriz de rigidez é:

$$
K_{IJ}
=
\int_\Omega
\operatorname{curl}_t\mathbf N_I
\operatorname{curl}_t\mathbf N_J
\,d\Omega.
$$

Organizando por componentes cartesianas:

$$
\boxed{
\mathbf K=
\begin{bmatrix}
\mathbf K_{xx}&\mathbf K_{xy}\\
\mathbf K_{yx}&\mathbf K_{yy}
\end{bmatrix}.
}
$$

Os quatro blocos são:

$$
\boxed{
K_{xx}^{ij}
=
\int_\Omega
\phi_{i,y}\phi_{j,y}\,d\Omega
}
$$

$$
\boxed{
K_{xy}^{ij}
=
-\int_\Omega
\phi_{i,y}\phi_{j,x}\,d\Omega
}
$$

$$
\boxed{
K_{yx}^{ij}
=
-\int_\Omega
\phi_{i,x}\phi_{j,y}\,d\Omega
}
$$

$$
\boxed{
K_{yy}^{ij}
=
\int_\Omega
\phi_{i,x}\phi_{j,x}\,d\Omega.
}
$$

Essas expressões já são suficientes para escrever uma rotina de montagem por quadratura numérica.

### 4.2. Matriz de massa

A matriz de massa é:

$$
M_{IJ}
=
\int_\Omega
\mathbf N_I\cdot\mathbf N_J\,d\Omega.
$$

Como as bases cartesianas são ortogonais:

$$
\boxed{
\mathbf M=
\begin{bmatrix}
\mathbf M_s&0\\
0&\mathbf M_s
\end{bmatrix},
}
$$

em que:

$$
\boxed{
M_s^{ij}
=
\int_\Omega
\phi_i\phi_j\,d\Omega.
}
$$

O problema resultante é:

$$
\boxed{
\mathbf K\mathbf h
=
k_c^2\mathbf M\mathbf h.
}
$$

Essa é a primeira formulação que eu implementaria.

Ela tem uma grande vantagem: precisamos apenas das funções MLS e de suas primeiras derivadas.

### 4.3. Por que aparecem modos espúrios?

Observe o rotacional da aproximação cartesiana:

$$
\operatorname{curl}_t\mathbf H_t^h
=
\sum_i
\left(
\phi_{i,x}h_{yi}
-
\phi_{i,y}h_{xi}
\right).
$$

Entretanto, a divergência é:

$$
\nabla_t\cdot\mathbf H_t^h
=
\sum_i
\left(
\phi_{i,x}h_{xi}
+
\phi_{i,y}h_{yi}
\right).
$$

Não existe nenhuma restrição explícita sobre os coeficientes $h_{xi}$ e $h_{yi}$ que garanta:

$$
\nabla_t\cdot\mathbf H_t^h=0.
$$

No problema contínuo, a divergência nula decorre da equação vetorial para autovalores positivos.

Na discretização, contudo, essa relação pode não ser preservada adequadamente.

A dificuldade envolve a compatibilidade entre os espaços discretos de gradiente e rotacional.

Por isso, uma aproximação vetorial construída simplesmente duplicando funções nodais escalares pode produzir espectros contaminados.

O resumo da dissertação confirma precisamente esse comportamento geral: a formulação cartesiana identifica os modos físicos, mas também apresenta modos espúrios. [ResearchGate](https://www.researchgate.net/publication/38437582_A_vector_Element-Free_Galerkin_method_for_waveguide_eigenvalue_problems?utm_source=chatgpt.com)

É importante distinguir duas situações: os autovalores nulos associados ao núcleo matemático do operador curl-curl não são necessariamente erros; já os autovalores positivos que não convergem para modos físicos representam contaminação espectral.

---

## 5. Segunda formulação: componentes rotacional e irrotacional

A segunda abordagem descrita pelo autor tenta melhorar a representação vetorial utilizando componentes rotacionais e irrotacionais.

Matematicamente, isso nos conduz à decomposição de Helmholtz.

Em duas dimensões, um campo vetorial suficientemente regular pode ser representado, sob hipóteses adequadas, por componentes da forma:

$$
\mathbf H_t
=
\nabla\psi+
\mathbf J\nabla\chi,
$$

onde:

$$
\mathbf J=
\begin{bmatrix}
0&-1\\
1&0
\end{bmatrix}.
$$

A primeira parcela é irrotacional:

$$
\operatorname{curl}_t(\nabla\psi)=0.
$$

A segunda é solenoidal:

$$
\nabla_t\cdot(\mathbf J\nabla\chi)=0.
$$

Assim, temos uma maneira natural de separar componentes de gradiente e componentes relacionadas ao rotacional.

**Esta é uma construção matemática adequada para investigar a segunda abordagem; ainda precisamos conferir a definição específica das funções vetoriais empregadas no capítulo 6 para afirmar equivalência exata com a implementação de Al Maleh.**

### 5.1. Funções de forma vetoriais

Podemos construir funções irrotacionais:

$$
\boxed{
\mathbf N_i^I=
\begin{bmatrix}
\phi_{i,x}\\
\phi_{i,y}
\end{bmatrix}
}
$$

e funções solenoidais:

$$
\boxed{
\mathbf N_i^R=
\begin{bmatrix}
-\phi_{i,y}\\
\phi_{i,x}
\end{bmatrix}.
}
$$

A aproximação passa a ser:

$$
\mathbf H_t^h
=
\sum_i a_i\mathbf N_i^I
+
\sum_i b_i\mathbf N_i^R.
$$

O rotacional da função irrotacional é nulo:

$$
\operatorname{curl}_t\mathbf N_i^I=0.
$$

Para a função rotacional:

$$
\operatorname{curl}_t\mathbf N_i^R
=
\phi_{i,xx}+\phi_{i,yy}.
$$

Portanto:

$$
\boxed{
\operatorname{curl}_t\mathbf N_i^R
=
\nabla_t^2\phi_i.
}
$$

Essa expressão revela uma dificuldade importante.

**Agora precisamos calcular segundas derivadas das funções MLS.**

Isso aumenta a sensibilidade do método à regularidade da função peso, ao condicionamento da aproximação local e à quadratura.

Para essa implementação eu consideraria inicialmente uma base polinomial quadrática:

$$
\mathbf p=
\begin{bmatrix}
1&x&y&x^2&xy&y^2
\end{bmatrix}^T,
$$

acompanhada de uma função peso suficientemente suave, pelo menos $C^2$.

### 5.2. Matriz de rigidez

Para a construção rotacional/irrotacional apresentada, a matriz possui a estrutura:

$$
\mathbf K=
\begin{bmatrix}
\mathbf K_{II}&\mathbf K_{IR}\\
\mathbf K_{RI}&\mathbf K_{RR}
\end{bmatrix}.
$$

Como:

$$
\operatorname{curl}_t\mathbf N_i^I=0,
$$

temos:

$$
\mathbf K_{II}=0,
$$

$$
\mathbf K_{IR}=0,
$$

$$
\mathbf K_{RI}=0.
$$

O bloco rotacional é:

$$
\boxed{
K_{RR}^{ij}
=
\int_\Omega
(\nabla_t^2\phi_i)
(\nabla_t^2\phi_j)
\,d\Omega.
}
$$

A matriz de massa, por sua vez, não é necessariamente diagonal por blocos.

Por exemplo:

$$
M_{II}^{ij}
=
\int_\Omega
\nabla\phi_i\cdot\nabla\phi_j\,d\Omega
$$

e:

$$
M_{RR}^{ij}
=
\int_\Omega
\nabla\phi_i\cdot\nabla\phi_j\,d\Omega.
$$

Existem também termos de acoplamento:

$$
M_{IR}^{ij}
=
\int_\Omega
\nabla\phi_i\cdot
\mathbf J\nabla\phi_j\,d\Omega.
$$

O autoproblema continua sendo:

$$
\boxed{
\mathbf K\mathbf q
=
k_c^2\mathbf M\mathbf q.
}
$$

Mas surge uma questão adicional: a independência linear das funções vetoriais.

Como as funções MLS reproduzem a partição da unidade:

$$
\sum_i\phi_i=1,
$$

temos:

$$
\sum_i\nabla\phi_i=0.
$$

Consequentemente, as funções vetoriais construídas exclusivamente a partir desses gradientes apresentam dependências lineares.

A matriz de massa poderá ser singular.

Isso exige uma escolha cuidadosa dos graus de liberdade ou uma redução de base antes da solução espectral.

Não devemos simplesmente enviar o sistema completo para `LAPACKE_dsygv` sem verificar a definitude da matriz de massa.

---
