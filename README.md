# Análise matemática e proposta de implementação do método de Al Maleh (2009)

**Minha recomendação é reproduzir primeiro as duas formulações de Al Maleh em um repositório próprio, utilizando o projeto NASA TP-3485 como referência de validação.** Depois, poderemos integrar os resultados ao projeto EFGMI-H(curl), no qual você já está investigando funções de Whitney difusas.

Examinei também seus repositórios existentes. Encontrei alguns pontos matemáticos e computacionais que merecem atenção antes de iniciarmos a implementação.

O trabalho de Hussam Al Maleh é especialmente interessante porque aborda um problema fundamental do eletromagnetismo computacional: **como aproximar campos vetoriais por um método sem malha sem introduzir modos próprios espúrios**.

Há, entretanto, uma limitação nesta análise: localizei o PDF original e confirmei o resumo, a organização dos capítulos e os resultados gerais, mas os servidores não permitiram recuperar integralmente as páginas nesta sessão. Assim, as equações abaixo são uma reconstrução matemática do problema e das duas classes de aproximação descritas pelo autor, não uma transcrição verificada de todas as equações da dissertação. Para uma reprodução fiel, será necessário conferir especialmente as expressões dos capítulos 5 e 6. [ResearchGate](https://www.researchgate.net/publication/38437582_A_vector_Element-Free_Galerkin_method_for_waveguide_eigenvalue_problems?utm_source=chatgpt.com)

Documento original: [A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems — McGill University](https://mcgill.scholaris.ca/bitstreams/4d79b125-540f-4b38-8c06-722d14eb60bb/download?utm_source=chatgpt.com).

---

## 1. Qual problema Al Maleh pretende resolver?

A dissertação considera um guia de ondas com paredes perfeitamente condutoras, cuja seção transversal ocupa um domínio bidimensional:

\[
\Omega\subset\mathbb{R}^{2}.
\]

O objetivo é determinar os **modos TM (Transverse Magnetic)** e seus respectivos números de onda de corte.

O campo magnético transversal é representado por:

\[
\mathbf H_t(x,y)=
\begin{bmatrix}
H_x(x,y)\\
H_y(x,y)
\end{bmatrix}.
\]

A escolha desse campo é importante.

Em vez de resolver diretamente o autoproblema escalar para a componente longitudinal do campo elétrico \(E_z\), o autor utiliza uma formulação vetorial na qual as incógnitas estão associadas às componentes transversais do campo magnético.

Isso permite investigar as dificuldades da aproximação vetorial pelo método EFG.

### 1.1. Partindo das equações de Maxwell

Para campos harmônicos, temos:

\[
\nabla\times\mathbf E=-j\omega\mu\mathbf H
\]

e

\[
\nabla\times\mathbf H=j\omega\varepsilon\mathbf E.
\]

Em um meio homogêneo, sem fontes, obtemos a equação vetorial:

\[
\nabla\times\nabla\times\mathbf H
=
k^2\mathbf H,
\]

em que

\[
k^2=\omega^2\mu\varepsilon.
\]

Para a formulação transversal de um guia uniforme, o autoproblema de corte pode ser escrito como:

\[
\boxed{
\nabla_t\times\nabla_t\times\mathbf H_t
=
k_c^2\mathbf H_t
}
\]

com

\[
k_c^2=k^2-\beta^2.
\]

Aqui, \(\beta\) é a constante de propagação longitudinal.

O objetivo numérico é determinar os autovalores:

\[
\lambda_i=k_{c,i}^2
\]

e os respectivos autovetores, que representam os modos eletromagnéticos.

### 1.2. A condição de divergência

Pelas equações de Maxwell:

\[
\nabla\cdot\mathbf B=0.
\]

Para permeabilidade constante:

\[
\nabla\cdot\mathbf H=0.
\]

Como estamos trabalhando com modos TM, para os quais \(H_z=0\), temos:

\[
\boxed{
\frac{\partial H_x}{\partial x}
+
\frac{\partial H_y}{\partial y}=0.
}
\]

Essa relação é fundamental.

Um campo magnético transversal físico deve ser solenoidal.

É justamente a capacidade da aproximação numérica de preservar essa estrutura que torna o problema interessante.

### 1.3. Condições de contorno

Para um guia PEC, o campo elétrico tangencial à superfície condutora deve desaparecer.

No problema TM, isso implica:

\[
E_z=0
\qquad \text{em }\partial\Omega.
\]

Em termos do campo magnético transversal, as soluções físicas satisfazem:

\[
\mathbf H_t\cdot\mathbf n=0
\qquad \text{em }\partial\Omega.
\]

Também podemos formular o problema vetorial com a condição natural:

\[
\operatorname{curl}_t\mathbf H_t=0
\qquad \text{em }\partial\Omega.
\]

Para os autovalores positivos da formulação curl-curl homogênea, essas relações são compatíveis.

**Atenção:** não devemos aplicar automaticamente a condição \(\mathbf n\times\mathbf H_t=0\). Essa não é a condição PEC apropriada à incógnita magnética transversal escolhida por Al Maleh.

Essa distinção será importante ao comparar a implementação com o módulo `helmvec` do repositório NASA, que também possui uma formulação em campo elétrico.

---

# 2. Formulação variacional do problema

A etapa seguinte é transformar o problema diferencial em uma formulação integral.

Considere:

\[
\nabla_t\times\nabla_t\times\mathbf H_t
=
k_c^2\mathbf H_t.
\]

Multiplicando por uma função teste vetorial \(\mathbf v\) e integrando no domínio:

\[
\int_\Omega
\mathbf v\cdot
(\nabla_t\times\nabla_t\times\mathbf H_t)
\,d\Omega
=
k_c^2
\int_\Omega
\mathbf v\cdot\mathbf H_t\,d\Omega.
\]

Aplicando a identidade de integração por partes e considerando a condição natural apropriada, obtemos:

\[
\boxed{
\int_\Omega
(\nabla_t\times\mathbf v)
(\nabla_t\times\mathbf H_t)
\,d\Omega
=
k_c^2
\int_\Omega
\mathbf v\cdot\mathbf H_t
\,d\Omega.
}
\]

Em duas dimensões, o rotacional de um campo transversal é um escalar:

\[
\operatorname{curl}_t\mathbf H_t
=
\frac{\partial H_y}{\partial x}
-
\frac{\partial H_x}{\partial y}.
\]

Essa é a formulação básica que queremos discretizar.

Definindo:

\[
a(\mathbf H,\mathbf v)
=
\int_\Omega
\operatorname{curl}_t\mathbf H
\operatorname{curl}_t\mathbf v\,d\Omega
\]

e

\[
m(\mathbf H,\mathbf v)
=
\int_\Omega
\mathbf H\cdot\mathbf v\,d\Omega,
\]

o problema variacional consiste em encontrar \(\mathbf H\neq0\) e \(\lambda\) tais que:

\[
\boxed{
a(\mathbf H,\mathbf v)
=
\lambda m(\mathbf H,\mathbf v).
}
\]

A partir daqui podemos utilizar diferentes aproximações para o campo vetorial.

Essa é precisamente a questão investigada na dissertação.

---

# 3. Construção das funções de forma pelo EFG

Antes de construir funções vetoriais, precisamos implementar corretamente as funções escalares EFG.

O método utiliza aproximações por mínimos quadrados móveis, ou MLS (*Moving Least Squares*).

Considere um conjunto de \(N\) nós:

\[
\mathcal N=
\{\mathbf x_1,\mathbf x_2,\ldots,\mathbf x_N\}.
\]

A aproximação escalar é:

\[
u^h(\mathbf x)
=
\mathbf p^T(\mathbf x)\mathbf a(\mathbf x).
\]

Para uma base polinomial linear:

\[
\mathbf p(\mathbf x)
=
\begin{bmatrix}
1\\x\\y
\end{bmatrix}.
\]

Os coeficientes \(\mathbf a(\mathbf x)\) são calculados minimizando:

\[
J(\mathbf a)
=
\sum_{i=1}^{N}
w_i(\mathbf x)
\left[
\mathbf p^T(\mathbf x_i)\mathbf a
-
u_i
\right]^2.
\]

A condição de mínimo fornece:

\[
\mathbf A(\mathbf x)\mathbf a(\mathbf x)
=
\mathbf B(\mathbf x)\mathbf u,
\]

com:

\[
\mathbf A(\mathbf x)
=
\sum_{i=1}^{N}
w_i(\mathbf x)
\mathbf p(\mathbf x_i)
\mathbf p^T(\mathbf x_i)
\]

e

\[
\mathbf B(\mathbf x)
=
\begin{bmatrix}
w_1\mathbf p(\mathbf x_1)&
\cdots&
w_N\mathbf p(\mathbf x_N)
\end{bmatrix}.
\]

Portanto:

\[
\mathbf a(\mathbf x)
=
\mathbf A^{-1}(\mathbf x)
\mathbf B(\mathbf x)\mathbf u.
\]

Substituindo na aproximação:

\[
u^h(\mathbf x)
=
\sum_{i=1}^{N}
\phi_i(\mathbf x)u_i,
\]

onde:

\[
\boxed{
\boldsymbol\phi^T(\mathbf x)
=
\mathbf p^T(\mathbf x)
\mathbf A^{-1}(\mathbf x)
\mathbf B(\mathbf x).
}
\]

Essas são as funções de forma MLS.

Na implementação, entretanto, não devemos calcular explicitamente a inversa de \(\mathbf A\). É preferível resolver os sistemas lineares por fatoração.

## 3.1. Derivadas das funções MLS

Este é um dos pontos que considero mais importantes para a implementação.

As funções de forma dependem da posição tanto por meio da base polinomial quanto pelas funções peso.

Logo:

\[
\phi_i(\mathbf x)
=
\mathbf p^T(\mathbf x)\mathbf q_i(\mathbf x),
\]

com:

\[
\mathbf q_i=\mathbf A^{-1}\mathbf B_i.
\]

A derivada em relação a uma coordenada \(\alpha\), com \(\alpha=x\) ou \(y\), é:

\[
\boxed{
\phi_{i,\alpha}
=
\mathbf p_{,\alpha}^{T}\mathbf q_i
+
\mathbf p^T\mathbf q_{i,\alpha}.
}
\]

Por sua vez:

\[
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
\]

Observe que **não podemos considerar os coeficientes MLS constantes durante a diferenciação**.

Isso gera derivadas incorretas das funções de forma.

### Uma descoberta no seu repositório existente

Examinei o arquivo `mls_shape.cpp` do seu repositório `simple-waveguide-efg-fem-solver`.

Atualmente, ele calcula as derivadas essencialmente por:

```cpp
result.dphidx[id] =
    (dpxdx.transpose() * coeff.col(j))(0);

result.dphidy[id] =
    (dpxdy.transpose() * coeff.col(j))(0);
```

Essas expressões não incluem as derivadas de \(\mathbf A^{-1}\mathbf B\).

Portanto, **não representam as derivadas completas das funções MLS**.

Encontrei também outro problema na montagem: a penalização de contorno é adicionada tanto à matriz de rigidez quanto à matriz de massa. Isso altera o autoproblema e pode introduzir autovalores artificiais.

Arquivos examinados: [mls_shape.cpp](https://github.com/Lucas-Sperotto/simple-waveguide-efg-fem-solver/blob/main/mls_shape.cpp?utm_source=chatgpt.com) e [efgm_solver.cpp](https://github.com/Lucas-Sperotto/simple-waveguide-efg-fem-solver/blob/main/efgm_solver.cpp?utm_source=chatgpt.com).

Isso significa que podemos aproveitar a organização desse projeto, mas não devemos reutilizar seu núcleo MLS sem uma correção matemática e novos testes.

---

# 4. Primeira formulação: aproximação pelas componentes cartesianas

Essa é a primeira abordagem descrita por Al Maleh.

O campo magnético é aproximado separadamente em suas componentes:

\[
H_x^h(x,y)
=
\sum_{i=1}^{N}
\phi_i(x,y)h_{xi}
\]

e

\[
H_y^h(x,y)
=
\sum_{i=1}^{N}
\phi_i(x,y)h_{yi}.
\]

Assim:

\[
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
\]

Temos \(2N\) incógnitas.

Podemos definir duas funções vetoriais para cada nó:

\[
\mathbf N_i^x=
\begin{bmatrix}
\phi_i\\
0
\end{bmatrix}
\]

e

\[
\mathbf N_i^y=
\begin{bmatrix}
0\\
\phi_i
\end{bmatrix}.
\]

O rotacional dessas funções é:

\[
\operatorname{curl}_t\mathbf N_i^x
=
-\frac{\partial\phi_i}{\partial y}
\]

e

\[
\operatorname{curl}_t\mathbf N_i^y
=
\frac{\partial\phi_i}{\partial x}.
\]

Agora podemos construir as matrizes globais.

## 4.1. Matriz de rigidez

A matriz de rigidez é:

\[
K_{IJ}
=
\int_\Omega
\operatorname{curl}_t\mathbf N_I
\operatorname{curl}_t\mathbf N_J
\,d\Omega.
\]

Organizando por componentes cartesianas:

\[
\boxed{
\mathbf K=
\begin{bmatrix}
\mathbf K_{xx}&\mathbf K_{xy}\\
\mathbf K_{yx}&\mathbf K_{yy}
\end{bmatrix}.
}
\]

Os quatro blocos são:

\[
\boxed{
K_{xx}^{ij}
=
\int_\Omega
\phi_{i,y}\phi_{j,y}\,d\Omega
}
\]

\[
\boxed{
K_{xy}^{ij}
=
-\int_\Omega
\phi_{i,y}\phi_{j,x}\,d\Omega
}
\]

\[
\boxed{
K_{yx}^{ij}
=
-\int_\Omega
\phi_{i,x}\phi_{j,y}\,d\Omega
}
\]

\[
\boxed{
K_{yy}^{ij}
=
\int_\Omega
\phi_{i,x}\phi_{j,x}\,d\Omega.
}
\]

Essas expressões já são suficientes para escrever uma rotina de montagem por quadratura numérica.

## 4.2. Matriz de massa

A matriz de massa é:

\[
M_{IJ}
=
\int_\Omega
\mathbf N_I\cdot\mathbf N_J\,d\Omega.
\]

Como as bases cartesianas são ortogonais:

\[
\boxed{
\mathbf M=
\begin{bmatrix}
\mathbf M_s&0\\
0&\mathbf M_s
\end{bmatrix},
}
\]

em que:

\[
\boxed{
M_s^{ij}
=
\int_\Omega
\phi_i\phi_j\,d\Omega.
}
\]

O problema resultante é:

\[
\boxed{
\mathbf K\mathbf h
=
k_c^2\mathbf M\mathbf h.
}
\]

Essa é a primeira formulação que eu implementaria.

Ela tem uma grande vantagem: precisamos apenas das funções MLS e de suas primeiras derivadas.

## 4.3. Por que aparecem modos espúrios?

Observe o rotacional da aproximação cartesiana:

\[
\operatorname{curl}_t\mathbf H_t^h
=
\sum_i
\left(
\phi_{i,x}h_{yi}
-
\phi_{i,y}h_{xi}
\right).
\]

Entretanto, a divergência é:

\[
\nabla_t\cdot\mathbf H_t^h
=
\sum_i
\left(
\phi_{i,x}h_{xi}
+
\phi_{i,y}h_{yi}
\right).
\]

Não existe nenhuma restrição explícita sobre os coeficientes \(h_{xi}\) e \(h_{yi}\) que garanta:

\[
\nabla_t\cdot\mathbf H_t^h=0.
\]

No problema contínuo, a divergência nula decorre da equação vetorial para autovalores positivos.

Na discretização, contudo, essa relação pode não ser preservada adequadamente.

A dificuldade envolve a compatibilidade entre os espaços discretos de gradiente e rotacional.

Por isso, uma aproximação vetorial construída simplesmente duplicando funções nodais escalares pode produzir espectros contaminados.

O resumo da dissertação confirma precisamente esse comportamento geral: a formulação cartesiana identifica os modos físicos, mas também apresenta modos espúrios. [ResearchGate](https://www.researchgate.net/publication/38437582_A_vector_Element-Free_Galerkin_method_for_waveguide_eigenvalue_problems?utm_source=chatgpt.com)

É importante distinguir duas situações: os autovalores nulos associados ao núcleo matemático do operador curl-curl não são necessariamente erros; já os autovalores positivos que não convergem para modos físicos representam contaminação espectral.

---

# 5. Segunda formulação: componentes rotacional e irrotacional

A segunda abordagem descrita pelo autor tenta melhorar a representação vetorial utilizando componentes rotacionais e irrotacionais.

Matematicamente, isso nos conduz à decomposição de Helmholtz.

Em duas dimensões, um campo vetorial suficientemente regular pode ser representado, sob hipóteses adequadas, por componentes da forma:

\[
\mathbf H_t
=
\nabla\psi+
\mathbf J\nabla\chi,
\]

onde:

\[
\mathbf J=
\begin{bmatrix}
0&-1\\
1&0
\end{bmatrix}.
\]

A primeira parcela é irrotacional:

\[
\operatorname{curl}_t(\nabla\psi)=0.
\]

A segunda é solenoidal:

\[
\nabla_t\cdot(\mathbf J\nabla\chi)=0.
\]

Assim, temos uma maneira natural de separar componentes de gradiente e componentes relacionadas ao rotacional.

**Esta é uma construção matemática adequada para investigar a segunda abordagem; ainda precisamos conferir a definição específica das funções vetoriais empregadas no capítulo 6 para afirmar equivalência exata com a implementação de Al Maleh.**

## 5.1. Funções de forma vetoriais

Podemos construir funções irrotacionais:

\[
\boxed{
\mathbf N_i^I=
\begin{bmatrix}
\phi_{i,x}\\
\phi_{i,y}
\end{bmatrix}
}
\]

e funções solenoidais:

\[
\boxed{
\mathbf N_i^R=
\begin{bmatrix}
-\phi_{i,y}\\
\phi_{i,x}
\end{bmatrix}.
}
\]

A aproximação passa a ser:

\[
\mathbf H_t^h
=
\sum_i a_i\mathbf N_i^I
+
\sum_i b_i\mathbf N_i^R.
\]

O rotacional da função irrotacional é nulo:

\[
\operatorname{curl}_t\mathbf N_i^I=0.
\]

Para a função rotacional:

\[
\operatorname{curl}_t\mathbf N_i^R
=
\phi_{i,xx}+\phi_{i,yy}.
\]

Portanto:

\[
\boxed{
\operatorname{curl}_t\mathbf N_i^R
=
\nabla_t^2\phi_i.
}
\]

Essa expressão revela uma dificuldade importante.

**Agora precisamos calcular segundas derivadas das funções MLS.**

Isso aumenta a sensibilidade do método à regularidade da função peso, ao condicionamento da aproximação local e à quadratura.

Para essa implementação eu consideraria inicialmente uma base polinomial quadrática:

\[
\mathbf p=
\begin{bmatrix}
1&x&y&x^2&xy&y^2
\end{bmatrix}^T,
\]

acompanhada de uma função peso suficientemente suave, pelo menos \(C^2\).

## 5.2. Matriz de rigidez

Para a construção rotacional/irrotacional apresentada, a matriz possui a estrutura:

\[
\mathbf K=
\begin{bmatrix}
\mathbf K_{II}&\mathbf K_{IR}\\
\mathbf K_{RI}&\mathbf K_{RR}
\end{bmatrix}.
\]

Como:

\[
\operatorname{curl}_t\mathbf N_i^I=0,
\]

temos:

\[
\mathbf K_{II}=0,
\]

\[
\mathbf K_{IR}=0,
\]

\[
\mathbf K_{RI}=0.
\]

O bloco rotacional é:

\[
\boxed{
K_{RR}^{ij}
=
\int_\Omega
(\nabla_t^2\phi_i)
(\nabla_t^2\phi_j)
\,d\Omega.
}
\]

A matriz de massa, por sua vez, não é necessariamente diagonal por blocos.

Por exemplo:

\[
M_{II}^{ij}
=
\int_\Omega
\nabla\phi_i\cdot\nabla\phi_j\,d\Omega
\]

e:

\[
M_{RR}^{ij}
=
\int_\Omega
\nabla\phi_i\cdot\nabla\phi_j\,d\Omega.
\]

Existem também termos de acoplamento:

\[
M_{IR}^{ij}
=
\int_\Omega
\nabla\phi_i\cdot
\mathbf J\nabla\phi_j\,d\Omega.
\]

O autoproblema continua sendo:

\[
\boxed{
\mathbf K\mathbf q
=
k_c^2\mathbf M\mathbf q.
}
\]

Mas surge uma questão adicional: a independência linear das funções vetoriais.

Como as funções MLS reproduzem a partição da unidade:

\[
\sum_i\phi_i=1,
\]

temos:

\[
\sum_i\nabla\phi_i=0.
\]

Consequentemente, as funções vetoriais construídas exclusivamente a partir desses gradientes apresentam dependências lineares.

A matriz de massa poderá ser singular.

Isso exige uma escolha cuidadosa dos graus de liberdade ou uma redução de base antes da solução espectral.

Não devemos simplesmente enviar o sistema completo para `LAPACKE_dsygv` sem verificar a definitude da matriz de massa.

---

# 6. O problema central: preservar as propriedades matemáticas do campo

Aqui está o ponto que mais me interessa nesse trabalho.

O problema não consiste apenas em obter bons autovalores.

Precisamos verificar se as funções vetoriais discretas preservam as propriedades fundamentais do operador eletromagnético.

No espaço contínuo, temos a sequência:

\[
H^1(\Omega)
\xrightarrow{\nabla}
H(\operatorname{curl};\Omega)
\xrightarrow{\operatorname{curl}}
L^2(\Omega).
\]

Uma propriedade importante é:

\[
\operatorname{curl}(\nabla u)=0.
\]

O espaço discreto deve reproduzir essa relação adequadamente.

Caso contrário, o núcleo do operador discreto pode conter campos que não correspondem à estrutura física e matemática desejada.

Essa é a motivação para utilizar elementos vetoriais compatíveis, como os elementos de Whitney/Nédélec.

### Uma ligação direta com seu projeto EFGMI-H(curl)

No repositório:

[efgmi-hcurl-diffuse-whitney](https://github.com/Lucas-Sperotto/efgmi-hcurl-diffuse-whitney?utm_source=chatgpt.com)

você já está investigando funções da forma:

\[
\boxed{
\mathbf W_{ij}
=
\phi_i\nabla\phi_j
-
\phi_j\nabla\phi_i.
}
\]

Essa construção é muito interessante porque busca associar a aproximação MLS a uma representação vetorial inspirada nos elementos de Whitney.

Ela também tem uma identidade importante:

\[
\nabla\phi_i
=
-\sum_j\mathbf W_{ij},
\]

quando os pares necessários estão presentes e a partição da unidade é satisfeita.

Portanto, existe uma conexão matemática direta entre essa formulação e o problema de Al Maleh.

Entretanto, encontrei uma questão importante na documentação atual do seu projeto.

No arquivo `TODO.md`, o diagnóstico registrado para o núcleo do operador apresenta:

\[
\dim(\ker K)=236,
\]

enquanto o valor esperado naquele experimento era:

\[
N-1=24.
\]

Isso significa que o teste de núcleo ainda não foi aprovado.

O resultado não prova, isoladamente, que todos esses graus de liberdade adicionais sejam modos espúrios físicos. Parte da discrepância pode decorrer da redundância das funções vetoriais, da política de pares virtuais ou de deficiências de posto do espaço discreto.

Mas mostra que ainda precisamos estabilizar a estrutura algébrica da formulação.

É por isso que considero importante reproduzir primeiro a formulação de Al Maleh separadamente.

Ela oferecerá uma terceira referência além dos seus métodos FEM e EFGMI.

---

# 7. Como validar matematicamente a implementação

O melhor caso inicial é um guia retangular homogêneo.

Considere:

\[
\Omega=(0,a)\times(0,b).
\]

Para os modos TM, a componente longitudinal elétrica satisfaz:

\[
\nabla_t^2 E_z+k_c^2E_z=0,
\]

com:

\[
E_z=0
\qquad\text{em }\partial\Omega.
\]

As autofunções são:

\[
E_z^{mn}(x,y)
=
\sin\left(\frac{m\pi x}{a}\right)
\sin\left(\frac{n\pi y}{b}\right),
\]

onde:

\[
m,n=1,2,3,\ldots
\]

Os autovalores analíticos são:

\[
\boxed{
k_{c,mn}^2
=
\left(\frac{m\pi}{a}\right)^2
+
\left(\frac{n\pi}{b}\right)^2.
}
\]

Podemos também construir o campo magnético transversal correspondente:

\[
\mathbf H_t^{mn}
\propto
\mathbf J\nabla_t E_z^{mn}.
\]

Logo, desprezando um fator multiplicativo constante:

\[
H_x^{mn}
=
-\frac{n\pi}{b}
\sin\left(\frac{m\pi x}{a}\right)
\cos\left(\frac{n\pi y}{b}\right)
\]

e:

\[
H_y^{mn}
=
\frac{m\pi}{a}
\cos\left(\frac{m\pi x}{a}\right)
\sin\left(\frac{n\pi y}{b}\right).
\]

Essa solução permite validar não apenas o autovalor, mas também a distribuição espacial do campo.

## 7.1. Primeiro experimento numérico

Sugiro utilizar inicialmente:

\[
a=2,\qquad b=1.
\]

Os primeiros modos TM apresentam:

| Modo | Autovalor analítico \(k_c^2\) |
|---|---|
| TM₁₁ | \(1{,}25\pi^2\) |
| TM₂₁ | \(2\pi^2\) |
| TM₃₁ | \(3{,}25\pi^2\) |
| TM₁₂ | \(4{,}25\pi^2\) |
| TM₂₂ | \(5\pi^2\) |
| TM₄₁ | \(5\pi^2\) |

Observe que TM₂₂ e TM₄₁ apresentam autovalores degenerados.

Isso é útil para testar o tratamento de autovalores múltiplos.

Também precisamos lembrar que **TM₁₁ é o modo TM de menor frequência de corte** nesse guia, enquanto TE₁₀ é o modo dominante considerando a família TE.

Como o problema original de Al Maleh utiliza \(\mathbf H_t\) para encontrar os modos TM, a validação inicial deve se concentrar nessa família.

## 7.2. Critérios de validação

Não considero suficiente comparar apenas os primeiros autovalores.

Precisamos acompanhar também:

\[
\varepsilon_{\lambda,i}
=
\frac{
|\lambda_i^h-\lambda_i^{exata}|
}{
|\lambda_i^{exata}|
}
\]

e um indicador para a divergência:

\[
\boxed{
\eta_{\operatorname{div},i}
=
\frac{
\|\nabla_t\cdot\mathbf H_{t,i}^h\|_{L^2(\Omega)}
}{
k_{c,i}\|\mathbf H_{t,i}^h\|_{L^2(\Omega)}
}.
}
\]

Para uma solução física bem aproximada, esperamos:

\[
\eta_{\operatorname{div},i}\rightarrow0
\]

com o refinamento.

Também devemos medir o resíduo do autoproblema:

\[
\boxed{
r_i=
\frac{
\|\mathbf K\mathbf h_i-\lambda_i\mathbf M\mathbf h_i\|_2
}{
(\|\mathbf K\|_2+
|\lambda_i|\|\mathbf M\|_2)\|\mathbf h_i\|_2
}.
}
\]

Esse indicador distingue problemas de solução do sistema algébrico de erros físicos da discretização.

Um autovalor espúrio pode apresentar resíduo algébrico pequeno. Por isso, precisamos analisar simultaneamente convergência, divergência, condições de contorno e semelhança com os campos analíticos.

---

# 8. Onde colocar essa implementação?

Agora podemos discutir os repositórios.

Examinei quatro projetos seus relacionados a esse problema.

| Repositório | Situação | Minha avaliação |
|---|---|---|
| `fem-eigen-em-2d3d-whitney` | Reprodução do NASA TP-3485, com os 14 casos documentados | Excelente referência para validação |
| `efgmi-hcurl-diffuse-whitney` | Investigação de funções vetoriais Whitney difusas | Excelente para comparação e evolução |
| `simple-waveguide-efg-fem-solver` | Implementação escalar FEM/EFGM | Pode fornecer componentes reutilizáveis após correções |
| `efg-3d-electromagnetics-2006` | Reprodução do trabalho de Parreira et al. | Menos relacionado diretamente ao autoproblema vetorial |

### Minha proposta

Criaria um repositório:

**`al-maleh-2009-vector-efg-waveguides`**

Ele teria um objetivo bem definido:

> Reproduzir e avaliar numericamente as formulações vetoriais do método Element-Free Galerkin propostas por Hussam Al Maleh (2009) para problemas de autovalores em guias de onda, com ênfase na identificação dos modos físicos TM, na análise dos modos espúrios e na comparação com soluções analíticas e pelo método dos elementos finitos.

A organização seria:

```text
al-maleh-2009-vector-efg-waveguides/
│
├── CMakeLists.txt
├── README.md
├── AGENTS.md
├── CITATION.cff
│
├── include/
│   ├── geometry/
│   ├── mls/
│   ├── vector_basis/
│   ├── assembly/
│   ├── boundary/
│   ├── eigen/
│   └── validation/
│
├── src/
│   ├── geometry/
│   ├── mls/
│   ├── vector_basis/
│   │   ├── cartesian.cpp
│   │   └── rotational_irrotational.cpp
│   ├── assembly/
│   ├── boundary/
│   └── eigen/
│
├── apps/
│   ├── scalar_reference.cpp
│   ├── cartesian_waveguide.cpp
│   ├── rotational_waveguide.cpp
│   └── compare_methods.cpp
│
├── tests/
│   ├── test_mls.cpp
│   ├── test_mls_derivatives.cpp
│   ├── test_vector_basis.cpp
│   ├── test_matrix_assembly.cpp
│   ├── test_rectangular_modes.cpp
│   └── test_spurious_modes.cpp
│
├── docs/
│   ├── theory/
│   ├── equations/
│   ├── implementation/
│   ├── validation/
│   ├── comparison/
│   └── references/
│
├── benchmarks/
│   ├── rectangular/
│   └── nasa_tp3485/
│
├── scripts/
│
└── results/
```

Manteria C++17/CMake, seguindo seus projetos de computação científica, com Eigen ou LAPACK para as operações lineares.

Também adotaria a mesma filosofia de rastreabilidade que utilizamos no projeto NASA:

\[
\text{Equação da dissertação}
\longrightarrow
\text{Formulação discreta}
\longrightarrow
\text{Função C++}
\longrightarrow
\text{Teste}
\longrightarrow
\text{Resultado}.
\]

Não misturaria inicialmente o código de Al Maleh ao núcleo de reprodução da NASA.

O repositório `fem-eigen-em-2d3d-whitney` já possui uma organização muito clara em torno do trabalho de Jin e dos 14 casos originais. Introduzir ali outra dissertação como método principal dificultaria a separação da autoria e dos resultados científicos.

**O NASA TP-3485 deve funcionar como referência independente, não como recipiente de todas as implementações.**

---

# 9. Planejamento da implementação

Minha proposta seria trabalhar em seis etapas, com testes e commits pequenos.

**Etapa 1 — Fundamentos matemáticos e MLS.** Reconstruir a formulação da dissertação, corrigir a implementação de MLS, implementar primeiras e segundas derivadas e executar testes de reprodução polinomial. Nenhuma matriz de autovalores seria montada antes desses testes.

**Etapa 2 — Formulação cartesiana.** Implementar as funções vetoriais, montar as matrizes \(\mathbf K\) e \(\mathbf M\) e resolver o primeiro autoproblema. A montagem deve passar por testes de simetria, integração e comparação com integrais conhecidas.

**Etapa 3 — Reprodução dos resultados físicos.** Utilizar o guia retangular, comparar os modos TM com a solução analítica e gerar gráficos dos campos magnéticos. Refinar a nuvem nodal e estudar a convergência.

**Etapa 4 — Investigação dos modos espúrios.** Analisar o núcleo de \(\mathbf K\), a condição de divergência, os resíduos espectrais e a sensibilidade aos parâmetros MLS. Reproduzir, na medida do possível, as condições numéricas da dissertação.

**Etapa 5 — Formulação rotacional/irrotacional.** Depois de conferir as expressões originais do capítulo 6, implementar a segunda representação vetorial, suas matrizes e o tratamento das dependências lineares.

**Etapa 6 — Comparação com FEM e EFGMI-H(curl).** Executar benchmarks equivalentes e investigar quais propriedades matemáticas e numéricas são preservadas por cada aproximação.

Ao final teremos três caminhos independentes de resolução vetorial para comparação: FEM com Whitney/Nédélec, EFG cartesiano de Al Maleh e EFGMI com Whitney difusas. A quarta implementação, rotacional/irrotacional de Al Maleh, ampliará esse conjunto.

---

# 10. O que eu investigaria além da dissertação

O resultado mais interessante não seria apenas reproduzir o trabalho de 2009.

Seria utilizar essa reprodução para estudar por que algumas discretizações sem malha são capazes de aproximar corretamente os campos eletromagnéticos, enquanto outras introduzem espectros contaminados.

Em particular, considero três questões matemáticas importantes.

**Primeira: a estrutura do núcleo do operador.** Precisamos determinar se a dimensão elevada do núcleo observada no projeto EFGMI-H(curl) está associada à redundância das funções de Whitney difusas, à seleção de pares virtuais ou à falta de compatibilidade discreta.

**Segunda: o efeito das derivadas MLS.** A correção das primeiras derivadas e o cálculo consistente das segundas derivadas podem afetar substancialmente os resultados numéricos. Vale investigar quanto do erro vem da representação vetorial e quanto vem da própria aproximação MLS.

**Terceira: a compatibilidade entre os espaços discretos.** A aproximação cartesiana utiliza diretamente duas famílias de funções escalares; a formulação rotacional utiliza derivadas dessas funções; e sua abordagem Whitney difusa constrói uma família vetorial por combinações antissimétricas. Podemos investigar a capacidade de cada espaço representar gradientes, campos solenoidais e autofunções de Maxwell.

Isso pode produzir uma comparação científica mais interessante do que simplesmente apresentar tabelas de erros de autovalores.

---

## 11. Minha conclusão

Eu **não começaria implementando diretamente a formulação rotacional/irrotacional**.

Começaria corrigindo e validando o núcleo MLS e depois desenvolveria a aproximação cartesiana de Al Maleh. Assim, conseguiremos reproduzir o comportamento principal descrito na dissertação antes de introduzir uma base vetorial mais complexa.

Também não considero adequado modificar agora o núcleo do repositório NASA para comportar esse método. Ele já é uma boa referência FEM, e preservar essa independência permitirá comparações mais confiáveis.

O caminho que recomendo é:

\[
\boxed{
\text{MLS validado}
\longrightarrow
\text{EFG cartesiano}
\longrightarrow
\text{Modos TM}
\longrightarrow
\text{Espúrios}
\longrightarrow
\text{EFG rotacional}
}
\]

Depois:

\[
\boxed{
\text{Al Maleh}
\;\;\longleftrightarrow\;\;
\text{NASA TP-3485}
\;\;\longleftrightarrow\;\;
\text{EFGMI-H(curl)}
}
\]

**Há uma oportunidade de pesquisa bastante concreta aqui:** investigar a estrutura e a compatibilidade dos espaços vetoriais sem malha, utilizando uma implementação independente de Al Maleh como referência para avaliar e aperfeiçoar suas funções de Whitney difusas.

Antes de escrevermos o código definitivo, a prioridade será obter a conferência integral das equações dos capítulos 5 e 6. Isso evitará transformar uma implementação matematicamente plausível em uma reprodução que atribua ao autor escolhas que ele não fez.

Não alterei nenhum dos seus repositórios nesta análise. A decisão arquitetural que considero mais sólida é criar `al-maleh-2009-vector-efg-waveguides` e reservar o NASA TP-3485 para validação cruzada.
