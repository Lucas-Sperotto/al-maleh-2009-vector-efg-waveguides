---
title: "01 — Equações de Maxwell, modos TM e formulação variacional"
author: "Projeto de reprodução de Al Maleh (2009)"
date: "20 de outubro de 2026"
lang: pt-BR
---

# Escopo e estatuto das afirmações

Este documento reconstrói a matemática de um autoproblema de guia PEC homogêneo com o campo magnético transversal como incógnita. O **resumo publicado** de Al Maleh (2009) confirma: modos TM, incógnita magnética transversal, EFG cartesiano e EFG rotacional/irrotacional, validação em guia retangular e ocorrência de modos espúrios. **As equações numeradas aqui são deduções próprias para implementação, não transcrições ou numeração da dissertação**. Antes de alegar reprodução fiel, confrontar capítulo, sinais, normalizações e imposição de fronteira com o PDF integral.

## 1. Convenções físicas

Seja $\Omega=(0,a)\times(0,b)$ a seção de um guia de paredes perfeitamente condutoras (PEC) e eixo longitudinal $z$. Usamos campos fasoriais com dependência temporal $e^{+j\omega t}$ e longitudinal $e^{-j\beta z}$, no meio linear isotrópico homogêneo. Denotamos $k_0=\omega\sqrt{\mu\varepsilon}$ (número de onda do material) e $k_c^2=k_0^2-\beta^2$ (autovalor de corte). Em material geral, preferir $k=\omega\sqrt{\mu\varepsilon}$ e reservar $k_0$ ao vácuo.

As equações de Maxwell, sem fontes, são

$$
\nabla\times\mathbf E=-j\omega\mu\mathbf H,\qquad
\nabla\times\mathbf H=j\omega\varepsilon\mathbf E,
\qquad \nabla\cdot(\mu\mathbf H)=0.
$$

Em um modo TM em relação a $z$, $H_z=0$ e $E_z\neq0$. Para parâmetros materiais constantes, obtemos

$$
\nabla_t\cdot\mathbf H_t=0,\qquad
\mathbf H_t=\gamma\,\mathbf J\nabla_t E_z,
\quad \mathbf J=\begin{bmatrix}0&-1\\1&0\end{bmatrix},
$$

em que $\gamma$ é um fator não nulo dependente da convenção fasorial e de $k_c$. O sinal e a fase de $\gamma$ não afetam os autovalores, mas precisam ser fixados ao reconstruir campos físicos.

## 2. Formulações escalar e vetorial

O campo longitudinal elétrico satisfaz

$$
-\Delta_t E_z=k_c^2E_z\quad\text{em }\Omega,
\qquad E_z=0\quad\text{em }\partial\Omega.
$$

O problema vetorial transversal que servirá de referência de implementação é

$$
\operatorname{curl}_t\operatorname{curl}_t\mathbf H_t=k_c^2\mathbf H_t,
\qquad \operatorname{curl}_t\mathbf H_t=\partial_xH_y-\partial_yH_x.
$$

A primeira operação de rotacional produz um escalar em 2D; a segunda deve ser interpretada como o rotacional vetorial do escalar: $\operatorname{curl}_t q=(\partial_yq,-\partial_xq)$, de modo que o operador composto é positivo. Essa definição evita ambiguidades de sinal; de maneira equivalente, escreva $-\mathbf J\nabla_t q$.

Para autovalores $k_c^2>0$, aplicar divergência à equação fornece $\nabla_t\cdot\mathbf H_t=0$. **Para $k_c^2=0$, essa inferência é inválida**. A existência do núcleo de gradientes deve ser tratada explicitamente.

## 3. Fronteira e espaços funcionais

A parede PEC impõe $\mathbf n\times\mathbf E=0$. Para TM, $E_z=0$ na fronteira; em consequência, o campo magnético físico tem $\mathbf H_t\cdot\mathbf n=0$ (salvo singularidades de canto). Pela relação entre $E_z$ e $\operatorname{curl}_t\mathbf H_t$ obtém-se também $\operatorname{curl}_t\mathbf H_t=0$ na fronteira para soluções regulares.

A formulação fraca curl-curl com condição natural $\operatorname{curl}_t\mathbf H_t=0$ é conveniente computacionalmente, mas **não se deve confundir o espaço fraco vetorial irrestrito com o subespaço físico TM**. Devemos testar divergência e fluxo normal dos autovetores calculados. Para usar condição normal como essencial, é necessário espaço conforme $H(\mathrm{div})$ ou uma imposição fraca adequada; o espaço nodal vetorial MLS é formalmente mais regular onde as funções são suaves, mas não tem propriedade interpolatória nodal.

O espaço de teste natural é $V=H(\mathrm{curl};\Omega)$. Defina

$$
a(\mathbf u,\mathbf v)=\int_\Omega\mu_r^{-1}
(\operatorname{curl}_t\mathbf u)(\operatorname{curl}_t\mathbf v)\,d\Omega,
\quad m(\mathbf u,\mathbf v)=\int_\Omega\varepsilon_r\mathbf u\cdot\mathbf v\,d\Omega.
$$

Para a reprodução inicial, assumir $\mu_r=\varepsilon_r=1$ evita confundir o autoproblema TM homogêneo com um problema de materiais não homogêneos; o tratamento dos pesos materiais deve ser derivado novamente quando necessário.

O problema fraco é: encontrar $(\lambda,\mathbf H)\in\mathbb R\times V$, $\mathbf H\ne0$, tal que

$$
a(\mathbf H,\mathbf v)=\lambda m(\mathbf H,\mathbf v)
\quad\forall\mathbf v\in V,\qquad\lambda=k_c^2.
$$

A integração por partes produz um termo de fronteira envolvendo $\operatorname{curl}_t\mathbf H$ e a componente tangencial da função teste; a condição natural o anula. Na implementação devem ser checados sinais e orientação do contorno por teste de manufatura, não por mera transcrição de fórmulas.

## 4. Solução analítica do guia retangular

Para o TM$_{mn}$, $m,n\ge1$:

$$
E_z^{mn}=E_0\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b},
\qquad\lambda_{mn}=k_{c,mn}^2=\left(\frac{m\pi}{a}\right)^2+\left(\frac{n\pi}{b}\right)^2.
$$

Tomando $\mathbf H_t=\mathbf J\nabla_tE_z$ sem fator físico de escala:

$$
H_x=-\frac{n\pi}{b}\sin\frac{m\pi x}{a}\cos\frac{n\pi y}{b},
\quad H_y=\frac{m\pi}{a}\cos\frac{m\pi x}{a}\sin\frac{n\pi y}{b}.
$$

O campo satisfaz $\nabla_t\cdot\mathbf H_t=0$, $\mathbf n\cdot\mathbf H_t=0$ e, na fronteira, $\operatorname{curl}_t\mathbf H_t=-\lambda E_z=0$ (com a convenção definida). Para $a=2$ e $b=1$, os primeiros autovalores TM incluem $1{,}25\pi^2$ (TM$_{11}$), $2\pi^2$ (TM$_{21}$), $3{,}25\pi^2$ (TM$_{31}$), $4{,}25\pi^2$ (TM$_{12}$), e $5\pi^2$ (TM$_{22}$ e TM$_{41}$, degenerados).

## 5. Conversões e falsos positivos

- $\lambda$ tem unidade $\mathrm{m}^{-2}$; $k_c=\sqrt\lambda$ tem $\mathrm{m}^{-1}$; $f_c=k_c/(2\pi\sqrt{\mu\varepsilon})$ tem $\mathrm{Hz}$.
- Não comparar um autovalor quadrático diretamente com a frequência ou com $k_c$.
- O menor modo do guia completo é TE$_{10}$, mas o menor **TM** é TM$_{11}$; não misturar famílias ao validar Al Maleh.
- Degenerescências exigem comparação de subespaços e não somente um-para-um de autovetores.

## Referências

Al Maleh, H. *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems*. Master of Engineering thesis, McGill University, 2009. URL: https://www.collectionscanada.gc.ca/obj/thesescanada/vol2/002/MR61527.PDF . Resumo verificado; equações internas pendentes de cotejo.

Jin, J. M. *Finite Element Method for Eigenvalue Problems in Electromagnetics*. NASA TP-3485, 1994, usado como referência FEM externa e independente.