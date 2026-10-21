# Modos espúrios e validação

> Análise inicial migrada do antigo README (linhas 913–1165), mantida para consulta e rastreabilidade. As equações são reconstruções de estudo, não transcrições verificadas de Al Maleh. Veja a [referência técnica](../theory/04_autovalores_espurios_validacao.md) e o [estado atual](../../README.md). Os delimitadores matemáticos foram adequados ao GitHub.

---

## 6. O problema central: preservar as propriedades matemáticas do campo

Aqui está o ponto que mais me interessa nesse trabalho.

O problema não consiste apenas em obter bons autovalores.

Precisamos verificar se as funções vetoriais discretas preservam as propriedades fundamentais do operador eletromagnético.

No espaço contínuo, temos a sequência:

$$
H^1(\Omega)
\xrightarrow{\nabla}
H(\operatorname{curl};\Omega)
\xrightarrow{\operatorname{curl}}
L^2(\Omega).
$$

Uma propriedade importante é:

$$
\operatorname{curl}(\nabla u)=0.
$$

O espaço discreto deve reproduzir essa relação adequadamente.

Caso contrário, o núcleo do operador discreto pode conter campos que não correspondem à estrutura física e matemática desejada.

Essa é a motivação para utilizar elementos vetoriais compatíveis, como os elementos de Whitney/Nédélec.

#### Uma ligação direta com seu projeto EFGMI-H(curl)

No repositório:

[efgmi-hcurl-diffuse-whitney](https://github.com/Lucas-Sperotto/efgmi-hcurl-diffuse-whitney?utm_source=chatgpt.com)

você já está investigando funções da forma:

$$
\boxed{
\mathbf W_{ij}
=
\phi_i\nabla\phi_j
-
\phi_j\nabla\phi_i.
}
$$

Essa construção é muito interessante porque busca associar a aproximação MLS a uma representação vetorial inspirada nos elementos de Whitney.

Ela também tem uma identidade importante:

$$
\nabla\phi_i
=
-\sum_j\mathbf W_{ij},
$$

quando os pares necessários estão presentes e a partição da unidade é satisfeita.

Portanto, existe uma conexão matemática direta entre essa formulação e o problema de Al Maleh.

Entretanto, encontrei uma questão importante na documentação atual do seu projeto.

No arquivo `TODO.md`, o diagnóstico registrado para o núcleo do operador apresenta:

$$
\dim(\ker K)=236,
$$

enquanto o valor esperado naquele experimento era:

$$
N-1=24.
$$

Isso significa que o teste de núcleo ainda não foi aprovado.

O resultado não prova, isoladamente, que todos esses graus de liberdade adicionais sejam modos espúrios físicos. Parte da discrepância pode decorrer da redundância das funções vetoriais, da política de pares virtuais ou de deficiências de posto do espaço discreto.

Mas mostra que ainda precisamos estabilizar a estrutura algébrica da formulação.

É por isso que considero importante reproduzir primeiro a formulação de Al Maleh separadamente.

Ela oferecerá uma terceira referência além dos seus métodos FEM e EFGMI.

---

## 7. Como validar matematicamente a implementação

O melhor caso inicial é um guia retangular homogêneo.

Considere:

$$
\Omega=(0,a)\times(0,b).
$$

Para os modos TM, a componente longitudinal elétrica satisfaz:

$$
\nabla_t^2 E_z+k_c^2E_z=0,
$$

com:

$$
E_z=0
\qquad\text{em }\partial\Omega.
$$

As autofunções são:

$$
E_z^{mn}(x,y)
=
\sin\left(\frac{m\pi x}{a}\right)
\sin\left(\frac{n\pi y}{b}\right),
$$

onde:

$$
m,n=1,2,3,\ldots
$$

Os autovalores analíticos são:

$$
\boxed{
k_{c,mn}^2
=
\left(\frac{m\pi}{a}\right)^2
+
\left(\frac{n\pi}{b}\right)^2.
}
$$

Podemos também construir o campo magnético transversal correspondente:

$$
\mathbf H_t^{mn}
\propto
\mathbf J\nabla_t E_z^{mn}.
$$

Logo, desprezando um fator multiplicativo constante:

$$
H_x^{mn}
=
-\frac{n\pi}{b}
\sin\left(\frac{m\pi x}{a}\right)
\cos\left(\frac{n\pi y}{b}\right)
$$

e:

$$
H_y^{mn}
=
\frac{m\pi}{a}
\cos\left(\frac{m\pi x}{a}\right)
\sin\left(\frac{n\pi y}{b}\right).
$$

Essa solução permite validar não apenas o autovalor, mas também a distribuição espacial do campo.

### 7.1. Primeiro experimento numérico

Sugiro utilizar inicialmente:

$$
a=2,\qquad b=1.
$$

Os primeiros modos TM apresentam:

| Modo | Autovalor analítico $k_c^2$ |
|---|---|
| TM₁₁ | $1{,}25\pi^2$ |
| TM₂₁ | $2\pi^2$ |
| TM₃₁ | $3{,}25\pi^2$ |
| TM₁₂ | $4{,}25\pi^2$ |
| TM₂₂ | $5\pi^2$ |
| TM₄₁ | $5\pi^2$ |

Observe que TM₂₂ e TM₄₁ apresentam autovalores degenerados.

Isso é útil para testar o tratamento de autovalores múltiplos.

Também precisamos lembrar que **TM₁₁ é o modo TM de menor frequência de corte** nesse guia, enquanto TE₁₀ é o modo dominante considerando a família TE.

Como o problema original de Al Maleh utiliza $\mathbf H_t$ para encontrar os modos TM, a validação inicial deve se concentrar nessa família.

### 7.2. Critérios de validação

Não considero suficiente comparar apenas os primeiros autovalores.

Precisamos acompanhar também:

$$
\varepsilon_{\lambda,i}
=
\frac{
|\lambda_i^h-\lambda_i^{exata}|
}{
|\lambda_i^{exata}|
}
$$

e um indicador para a divergência:

$$
\boxed{
\eta_{\operatorname{div},i}
=
\frac{
\|\nabla_t\cdot\mathbf H_{t,i}^h\|_{L^2(\Omega)}
}{
k_{c,i}\|\mathbf H_{t,i}^h\|_{L^2(\Omega)}
}.
}
$$

Para uma solução física bem aproximada, esperamos:

$$
\eta_{\operatorname{div},i}\rightarrow0
$$

com o refinamento.

Também devemos medir o resíduo do autoproblema:

$$
\boxed{
r_i=
\frac{
\|\mathbf K\mathbf h_i-\lambda_i\mathbf M\mathbf h_i\|_2
}{
(\|\mathbf K\|_2+
|\lambda_i|\|\mathbf M\|_2)\|\mathbf h_i\|_2
}.
}
$$

Esse indicador distingue problemas de solução do sistema algébrico de erros físicos da discretização.

Um autovalor espúrio pode apresentar resíduo algébrico pequeno. Por isso, precisamos analisar simultaneamente convergência, divergência, condições de contorno e semelhança com os campos analíticos.

---
