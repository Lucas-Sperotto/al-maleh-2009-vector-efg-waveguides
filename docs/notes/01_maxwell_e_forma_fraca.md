# Maxwell e formulação fraca

> Análise inicial migrada do antigo README (linhas 27–257), mantida para consulta e rastreabilidade. As equações são reconstruções de estudo, não transcrições verificadas de Al Maleh. Veja a [referência técnica](../theory/01_maxwell_tm_forma_fraca.md) e o [estado atual](../../README.md). Os delimitadores matemáticos foram adequados ao GitHub.

---

### 1. Qual problema Al Maleh pretende resolver?

A dissertação considera um guia de ondas com paredes perfeitamente condutoras, cuja seção transversal ocupa um domínio bidimensional:

$$
\Omega\subset\mathbb{R}^{2}.
$$

O objetivo é determinar os **modos TM (Transverse Magnetic)** e seus respectivos números de onda de corte.

O campo magnético transversal é representado por:

$$
\mathbf H_t(x,y)=
\begin{bmatrix}
H_x(x,y)\\
H_y(x,y)
\end{bmatrix}.
$$

A escolha desse campo é importante.

Em vez de resolver diretamente o autoproblema escalar para a componente longitudinal do campo elétrico $E_z$, o autor utiliza uma formulação vetorial na qual as incógnitas estão associadas às componentes transversais do campo magnético.

Isso permite investigar as dificuldades da aproximação vetorial pelo método EFG.

#### 1.1. Partindo das equações de Maxwell

Para campos harmônicos, temos:

$$
\nabla\times\mathbf E=-j\omega\mu\mathbf H
$$

e

$$
\nabla\times\mathbf H=j\omega\varepsilon\mathbf E.
$$

Em um meio homogêneo, sem fontes, obtemos a equação vetorial:

$$
\nabla\times\nabla\times\mathbf H
=
k^2\mathbf H,
$$

em que

$$
k^2=\omega^2\mu\varepsilon.
$$

Para a formulação transversal de um guia uniforme, o autoproblema de corte pode ser escrito como:

$$
\boxed{
\nabla_t\times\nabla_t\times\mathbf H_t
=
k_c^2\mathbf H_t
}
$$

com

$$
k_c^2=k^2-\beta^2.
$$

Aqui, $\beta$ é a constante de propagação longitudinal.

O objetivo numérico é determinar os autovalores:

$$
\lambda_i=k_{c,i}^2
$$

e os respectivos autovetores, que representam os modos eletromagnéticos.

#### 1.2. A condição de divergência

Pelas equações de Maxwell:

$$
\nabla\cdot\mathbf B=0.
$$

Para permeabilidade constante:

$$
\nabla\cdot\mathbf H=0.
$$

Como estamos trabalhando com modos TM, para os quais $H_z=0$, temos:

$$
\boxed{
\frac{\partial H_x}{\partial x}
+
\frac{\partial H_y}{\partial y}=0.
}
$$

Essa relação é fundamental.

Um campo magnético transversal físico deve ser solenoidal.

É justamente a capacidade da aproximação numérica de preservar essa estrutura que torna o problema interessante.

#### 1.3. Condições de contorno

Para um guia PEC, o campo elétrico tangencial à superfície condutora deve desaparecer.

No problema TM, isso implica:

$$
E_z=0
\qquad \text{em }\partial\Omega.
$$

Em termos do campo magnético transversal, as soluções físicas satisfazem:

$$
\mathbf H_t\cdot\mathbf n=0
\qquad \text{em }\partial\Omega.
$$

Também podemos formular o problema vetorial com a condição natural:

$$
\operatorname{curl}_t\mathbf H_t=0
\qquad \text{em }\partial\Omega.
$$

Para os autovalores positivos da formulação curl-curl homogênea, essas relações são compatíveis.

**Atenção:** não devemos aplicar automaticamente a condição $\mathbf n\times\mathbf H_t=0$. Essa não é a condição PEC apropriada à incógnita magnética transversal escolhida por Al Maleh.

Essa distinção será importante ao comparar a implementação com o módulo `helmvec` do repositório NASA, que também possui uma formulação em campo elétrico.

---

## 2. Formulação variacional do problema

A etapa seguinte é transformar o problema diferencial em uma formulação integral.

Considere:

$$
\nabla_t\times\nabla_t\times\mathbf H_t
=
k_c^2\mathbf H_t.
$$

Multiplicando por uma função teste vetorial $\mathbf v$ e integrando no domínio:

$$
\int_\Omega
\mathbf v\cdot
(\nabla_t\times\nabla_t\times\mathbf H_t)
\,d\Omega
=
k_c^2
\int_\Omega
\mathbf v\cdot\mathbf H_t\,d\Omega.
$$

Aplicando a identidade de integração por partes e considerando a condição natural apropriada, obtemos:

$$
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
$$

Em duas dimensões, o rotacional de um campo transversal é um escalar:

$$
\operatorname{curl}_t\mathbf H_t
=
\frac{\partial H_y}{\partial x}
-
\frac{\partial H_x}{\partial y}.
$$

Essa é a formulação básica que queremos discretizar.

Definindo:

$$
a(\mathbf H,\mathbf v)
=
\int_\Omega
\operatorname{curl}_t\mathbf H
\operatorname{curl}_t\mathbf v\,d\Omega
$$

e

$$
m(\mathbf H,\mathbf v)
=
\int_\Omega
\mathbf H\cdot\mathbf v\,d\Omega,
$$

o problema variacional consiste em encontrar $\mathbf H\neq0$ e $\lambda$ tais que:

$$
\boxed{
a(\mathbf H,\mathbf v)
=
\lambda m(\mathbf H,\mathbf v).
}
$$

A partir daqui podemos utilizar diferentes aproximações para o campo vetorial.

Essa é precisamente a questão investigada na dissertação.

---
