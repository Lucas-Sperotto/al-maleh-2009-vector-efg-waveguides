---
title: "03 — Formulações vetoriais EFG: cartesiana e rotacional/irrotacional"
author: "Projeto de reprodução de Al Maleh (2009)"
date: "20 de outubro de 2026"
lang: pt-BR
---

# Nota de rastreabilidade

O resumo de Al Maleh (2009) afirma que foram exploradas duas representações: (A) componentes cartesianas; (B) decomposição rotacional e irrotacional. **A formulação (A), cartesiana, foi agora conferida com as imagens do capítulo 5 da dissertação, em particular (5-1), (5-2), (5-7)–(5-10): a estrutura de funções vetoriais e das formas bilineares coincide com o que foi montado.** Uma possível diferença de sinal na equação (5-6) e os parâmetros originais de MLS são registrados em [M2, conferência do capítulo 5](../implementation/02_cartesian_assembly.md). Os valores percentuais efetivamente medidos pelo nosso código **não coincidem integralmente** com as Tabelas 2/3 (veja [campanha calculada](../results/M2_CH5_MEASURED_2026-10-21.md)). **A formulação (B) foi agora conferida diretamente com as imagens do Capítulo 6, Eqs. (6-1)–(6-6): o autor usa exatamente grad(N_i) e z×grad(N_i) como funções de forma independentes por nó.** A implementação M5 segue essas funções e as formas de massa/rigidez do capítulo 5; o tratamento original de degenerescências de massa, o peso MLS e parâmetros ainda exigem cotejo, de modo que não alegamos igualdade experimental sem medição. Ver [M5 Capítulo 6](../implementation/11_rot_irrot_original_chapter6.md).

## 1. Formulação A: base cartesiana

Sejam $\phi_i\in H^1(\Omega)$ as funções MLS; usar em cada nó duas funções vetoriais

$$
\mathbf N_i^{x}=\begin{bmatrix}\phi_i\\0\end{bmatrix},\qquad
\mathbf N_i^{y}=\begin{bmatrix}0\\\phi_i\end{bmatrix}.
$$

O campo aproximado é

$$
\mathbf H_h=\sum_{i=1}^{N}\left(h_{xi}\mathbf N_i^x+h_{yi}\mathbf N_i^y\right),\qquad
\mathbf h=(h_{x1},\ldots,h_{xN},h_{y1},\ldots,h_{yN})^T.
$$

Com $\operatorname{curl}_t(u_x,u_y)=\partial_xu_y-\partial_yu_x$:

$$
\operatorname{curl}_t\mathbf N_i^x=-\phi_{i,y},\qquad
\operatorname{curl}_t\mathbf N_i^y=\phi_{i,x}.
$$

## 2. Matrizes cartesianas

Para $\mu_r=\varepsilon_r=1$, a matriz $K_{IJ}=\int_\Omega\operatorname{curl}N_I\operatorname{curl}N_J$ e a matriz $M_{IJ}=\int_\Omega N_I\cdot N_J$ têm os blocos:

$$
\mathbf K=\begin{bmatrix}K_{xx}&K_{xy}\\K_{yx}&K_{yy}\end{bmatrix},
\qquad \mathbf M=\begin{bmatrix}M_s&0\\0&M_s\end{bmatrix}.
$$

$$
(K_{xx})_{ij}=\int_\Omega\phi_{i,y}\phi_{j,y}\,d\Omega,
\quad (K_{xy})_{ij}=-\int_\Omega\phi_{i,y}\phi_{j,x}\,d\Omega,
$$

$$
(K_{yx})_{ij}=-\int_\Omega\phi_{i,x}\phi_{j,y}\,d\Omega,
\quad (K_{yy})_{ij}=\int_\Omega\phi_{i,x}\phi_{j,x}\,d\Omega,
$$

$$
(M_s)_{ij}=\int_\Omega\phi_i\phi_j\,d\Omega,
\qquad K_{yx}=K_{xy}^T.
$$

O autoproblema é $K\mathbf h=\lambda M\mathbf h$ com $\lambda=k_c^2$. É necessário distinguir o núcleo de $K$ de autovalores positivos espúrios; o núcleo de gradientes é parte da matemática do operador.

## 3. Montagem por quadratura

Escolha células de integração de fundo $C$ e pontos gaussianos $\mathbf x_q$ com pesos físicos $W_q=w_q|J_C(\xi_q)|$. A malha de integração **não é malha de interpolação**; os graus de liberdade pertencem à nuvem de nós EFG. Em cada ponto, calcular $\phi_i$, $\phi_{i,x}$ e $\phi_{i,y}$ para o conjunto local de nós ativos. Para cada par $i,j$:

$$
K_{ij}\gets K_{ij}+W_q\,c_i c_j,\qquad
M_{ij}\gets M_{ij}+W_q\,\mathbf N_i\cdot\mathbf N_j,
$$

onde $c_i=\operatorname{curl}\mathbf N_i$ para cada DOF cartesiano. Montar os blocos locais por produtos externos reduz a duplicação de laços e facilita testes de simetria. Agregar índices esparsos consistentemente; registrar quadratura e raio de suporte em cada execução.

### Pseudocódigo da base cartesiana

```text
for cell in integration_cells:
  for q in gauss_points(cell):
    sf = mls.evaluate(q.position, first_derivatives=true)
    dofs = cartesian_dofs(sf.active_nodes)
    c = [-sf.dphidy, +sf.dphidx]     # mesmo ordenamento de dofs
    B = [ [phi_i, 0] / [0, phi_i] ]  # duas componentes por nó
    K[dofs,dofs] += q.weight * outer(c,c)
    M[dofs,dofs] += q.weight * (B^T B)
assert symmetry(K,M)
solve_generalized_eigenproblem(K,M)
```

## 4. Formulação B: hipótese de decomposição de Helmholtz

Em domínio simplesmente conexo, com escolhas de fronteira e regularidade apropriadas, considere a decomposição formal

$$
\mathbf H=\nabla\psi+\mathbf J\nabla\chi,
\qquad \mathbf J=\begin{bmatrix}0&-1\\1&0\end{bmatrix}.
$$

A parcela $\nabla\psi$ é irrotacional e a parcela $\mathbf J\nabla\chi$ é solenoidal. Construir funções experimentais

$$
\mathbf N_i^{I}=\nabla\phi_i,\qquad
\mathbf N_i^{R}=\mathbf J\nabla\phi_i
=\begin{bmatrix}-\phi_{i,y}\\\phi_{i,x}\end{bmatrix}.
$$

Suas identidades são

$$
\operatorname{curl}\mathbf N_i^I=0,\qquad
\nabla\cdot\mathbf N_i^R=0,\qquad
\operatorname{curl}\mathbf N_i^R=\Delta\phi_i.
$$

Portanto, as segundas derivadas MLS tornam-se necessárias. Essa decomposição não equivale automaticamente a uma decomposição ortogonal em qualquer condição de contorno; termos de fronteira e topologia importam.

## 5. Blocos da formulação B proposta

Escreva $\mathbf H_h=\sum_i a_i\mathbf N_i^I+\sum_i b_i\mathbf N_i^R$ e $\mathbf q=(\mathbf a,\mathbf b)^T$. A matriz de rigidez satisfaz

$$
K=\begin{bmatrix}0&0\\0&K_{RR}\end{bmatrix},\qquad
(K_{RR})_{ij}=\int_\Omega(\Delta\phi_i)(\Delta\phi_j)\,d\Omega.
$$

A matriz de massa é

$$
M=\begin{bmatrix}M_{II}&M_{IR}\\M_{RI}&M_{RR}\end{bmatrix},
\qquad (M_{II})_{ij}=(M_{RR})_{ij}=\int_\Omega\nabla\phi_i\cdot\nabla\phi_j.
$$

Os blocos cruzados são

$$
(M_{IR})_{ij}=\int_\Omega\nabla\phi_i\cdot\mathbf J\nabla\phi_j\,d\Omega,
\quad M_{RI}=M_{IR}^T.
$$

A identidade $\sum_i\phi_i=1$ implica $\sum_i\nabla\phi_i=0$. Isso introduz dependências lineares entre funções vetoriais; pode tornar a matriz de massa singular. Solução espectral requer identificar o espaço nulo, construir projeção de posto completo ou formular redução compatível, sem simplesmente descartar autovalores desconfortáveis.

## 6. Condições de contorno: não aplicar eliminação nodal ingênua

Funções MLS não possuem, em geral, $\phi_i(\mathbf x_j)=\delta_{ij}$, portanto zerar coeficientes associados a nós geométricos de fronteira não impõe automaticamente condição Dirichlet. Para o escalar $E_z=0$, estratégias robustas incluem multiplicadores de Lagrange, Nitsche consistente ou penalização **somente com justificativa da forma fraca**. A formulação vetorial $H_t$ deve respeitar a condição física apropriada; termo de penalização normal $\eta\int_{\partial\Omega}(\mathbf n\cdot\mathbf H)^2$ requer estudo do espectro alterado e de dependência de $\eta$.

## 7. Testes matemáticos mínimos

1. Simetria relativa de $K$ e $M$, incluindo blocos cruzados.
2. Não negatividade da forma de energia $\mathbf v^TK\mathbf v$; definitude de massa no espaço reduzido.
3. Consistência de $\operatorname{curl}\mathbf N_i^x$, $\operatorname{curl}\mathbf N_i^y$, $\operatorname{curl}\mathbf N_i^I$ e $\operatorname{curl}\mathbf N_i^R$ por diferenciação independente.
4. Dimensão do núcleo de $K$ e classificação de direções linearmente dependentes; nunca assumir que todo autovalor pequeno é espúrio.
5. Convergência em refinamentos e comparação dos campos por correlação modal, além dos autovalores.

## Referência primária

Al Maleh, H. 2009, *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems*. As formulações de funções de forma foram verificadas contra os capítulos 5 e 6, incluindo (6-1)–(6-6). Ainda faltam o cotejo dos pesos/raios MLS e a abordagem original para a matriz de massa singular. As Tabelas 5/6 originais estão transcritas, mas os resultados numéricos precisam ser confrontados separadamente.
## Auditoria de 2026-10-22: propriedades demonstradas da formulação rotacional/irrotacional

1. **Pencil singular (teorema).** Se f=φ+iψ é polinômio holomorfo
   reproduzido pela MLS, ∇φ+ẑ×∇ψ≡0 (Cauchy–Riemann) e Δψ≡0; o vetor de
   coeficientes pertence a ker M ∩ ker K. Grau 1: f∈{1,i,z,iz} (4 vetores);
   grau 2: mais z², iz² (6). Logo det(K−λM)≡0 e o problema exige deflação.
2. **Equivalência da redução (teorema).** Com Z base do núcleo descartado de
   M: se KZ=0, os autopares de QᵀKQ são autopares de Kh=λMh com resíduo
   completo nulo; caso contrário os autovalores finitos são os do
   complemento de Schur e a projeção dá cotas superiores.
3. **Núcleo de rotacional (teorema + observação).** dim ker(K_red)=dim ker K −
   dim ker M ≥ n−1 (gradientes discretos ∇V_h); igualdade observada.
4. **Positividade (teorema).** Com pesos de Gauss positivos, QᵀKQ é
   semidefinida: autovalores negativos só em nível de arredondamento.
5. **Ortogonalidade aos gradientes (teorema).** Autovetores com λ≠0 são
   M-ortogonais a ∇V_h (divergência fraca discreta nula e n·H=0 natural).

Validação numérica dessas propriedades e reprodução das Tabelas 4–6:
[relatório](../results/INDEPENDENT_THEORY_CODE_AUDIT.md).
