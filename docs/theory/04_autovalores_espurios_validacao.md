---
title: "04 — Autovalores, modos espúrios e protocolo de validação"
author: "Projeto de reprodução de Al Maleh (2009)"
date: "20 de outubro de 2026"
lang: pt-BR
---

# 1. Por que há modos espúrios?

O autoproblema curl-curl não é coercivo no espaço vetorial completo, pois qualquer gradiente suave tem rotacional zero. Esquematicamente,

$$
\nabla\phi\in\ker(\operatorname{curl}),\qquad
\nabla\phi\in\ker(\operatorname{curl}\operatorname{curl}).
$$

Os autovalores nulos do núcleo contínuo não são por si mesmos modos espúrios: são propriedade estrutural do operador. Chamaremos de **modo espúrio positivo** uma autofunção discretizada associada a autovalor positivo que não converge ao espectro físico alvo sob refinamento, ou viola as restrições físicas sem desaparecer numericamente. Diferenciar ainda vetores de coeficientes redundantes, em que $M$ é singular, de modos eletromagnéticos de energia não nula.

## 2. Compatibilidade de espaços

No contínuo existe o complexo

$$
H^1(\Omega)\xrightarrow{\nabla}H(\operatorname{curl};\Omega)
\xrightarrow{\operatorname{curl}}L^2(\Omega).
$$

O operador curl anula exatamente a imagem do gradiente. Aproximações com funções vetoriais arbitrárias podem não reproduzir uma sequência discreta exata com propriedades de aproximação espectral adequadas. A definição da base, o quadraturamento, as restrições de contorno, a redundância e os projetores influenciam o espectro.

No projeto EFGMI-H(curl) existente, a função Whitney difusa proposta é $W_{ij}=\phi_i\nabla\phi_j-\phi_j\nabla\phi_i$. O repositório registra um teste com $\dim\ker K=236$ versus valor esperado de 24, ainda não aprovado. **Esse número não deve ser importado como prova de 212 modos espúrios físicos**: é indispensável separar dependência linear da base, núcleo real de curl e espaços de gradiente.

## 3. Solucionador generalizado

Resolvemos

$$
K\mathbf h=\lambda M\mathbf h.
$$

Para usar uma rotina simétrica definida como LAPACK `dsygv`, $K$ e $M$ devem ser simétricas e $M$ deve ser definida positiva no subespaço usado. Caso $M$ seja semidefinida ou singular, identificar seu posto, projetar a um espaço de massa positiva ou escolher um solver generalizado apropriado; não usar deslocamentos diagonais arbitrários sem documentar seu efeito no espectro.

Os resíduos relativos são

$$
\rho_i=\frac{\|K\mathbf h_i-\lambda_iM\mathbf h_i\|_2}
{(\|K\|_2+|\lambda_i|\|M\|_2)\|\mathbf h_i\|_2}.
$$

Para avaliação física de modos positivos,

$$
\eta_{{\mathrm{div}},i}=\frac{\|\nabla_t\cdot\mathbf H_i\|_{L^2(\Omega)}}
{k_{c,i}\|\mathbf H_i\|_{L^2(\Omega)}},
\qquad
\eta_{n,i}=\frac{\|\mathbf n\cdot\mathbf H_i\|_{L^2(\partial\Omega)}}
{\|\mathbf H_i\|_{L^2(\Omega)}}.
$$

Os indicadores dependem de norma, dimensão e discretização; no artigo experimental fixar escalas, quadratura e tolerâncias explicitamente. O valor pequeno de $\rho_i$ **não garante** que o modo seja físico.

## 4. Benchmarks analíticos

No guia retangular PEC homogêneo, $\lambda_{mn}=(m\pi/a)^2+(n\pi/b)^2$ com $m,n\ge1$ para TM. Para $a=2$ e $b=1$:

| Modo TM | $\lambda_{mn}$ | Aproximação numérica |
|:--|--:|--:|
| TM11 | $1.25\pi^2$ | 12.3370055 |
| TM21 | $2\pi^2$ | 19.7392088 |
| TM31 | $3.25\pi^2$ | 32.0762143 |
| TM12 | $4.25\pi^2$ | 41.9458182 |
| TM22 e TM41 | $5\pi^2$ | 49.3480220 |

A duplicidade TM22/TM41 exige casamento modal por **subespaço**. Em cada campo reconstruído compare correlação normalizada $|\langle H_h,H_{\mathrm{ref}}\rangle|/(\|H_h\|\|H_{\mathrm{ref}}\|)$; os sinais/fases de autovetores são arbitrários.

## 5. Protocolo de refinamento

Executar, no mínimo, três nuvens regulares ou quase regulares com $N_x\times N_y$ crescendo, mantendo um critério transparente para o raio de influência em unidades de espaçamento nodal. Repetir com quadraturas de ordens diferentes para separar erro de integração de erro de aproximação. Em cada rodada exportar:

- configurações: $a$, $b$, $N_x$, $N_y$, suporte, base, peso, quadratura, estratégia de fronteira, versão do código;
- medidas numéricas: posto de MLS, $\kappa(A)$ máximo, simetria de $K/M$, posto de $M$, dimensão numérica do núcleo de $K$, tolerâncias;
- espectro: $\lambda_i$, $k_{c,i}$, resíduos, etiqueta de modo, erro analítico relativo, estabilidade sob refinamento;
- campos: $H_x$, $H_y$, $\nabla\cdot H$, $\operatorname{curl}H$, fluxo normal nas paredes; visualização VTK/CSV;
- custo: montagem, solução, memória, tamanho esparso e densidade de não nulos.

## 6. Classes de falhas a distinguir

| Observação | Hipóteses a testar | Evidência requerida |
|:--|:--|:--|
| Massa singular | base redundante; suporte degenerado; quadratura insuficiente | posto numérico e vetores nulos de $M$ |
| Muitos zeros em $K$ | núcleo físico de gradientes; redundância; erro de rank | comparação $\ker K$ versus $\operatorname{im}\nabla_h$ |
| Modo positivo sem correspondência TM | espaço incompatível; BC incorreta; integração fraca | convergência, campo, divergência e fronteira |
| Autovalor convergente mas campo errado | matching degenerado; normalização; modo misturado | projeção sobre subespaço analítico |
| Sensibilidade ao suporte | MLS mal condicionado; peso não suave | estudo paramétrico de suporte e $\kappa(A)$ |

## 7. Critérios de aceitação

O primeiro milestone só deve ser anunciado como **reprodução numérica preliminar** após: (a) MLS reproduzir polinômios com derivadas corretas; (b) simetria e massa válida; (c) pelo menos três níveis de refinamento; (d) primeiros modos TM com erro decrescente; (e) análise explícita dos modos que não convergem; (f) documentação de diferenças entre este código e a tese original. Definir os limiares de erro antes do experimento, não após visualizar os dados.

## 8. Conexão com NASA TP-3485

O projeto `fem-eigen-em-2d3d-whitney` já documenta matrizes, casos TE/TM, exportação VTK/CSV e comparações. Aqui ele será referência externa em configuração física equivalente. Não tomar coincidência numérica de autovalor como confirmação de que a base vetor EFG é compatible com Maxwell. Preservar resultados dos dois códigos separados e publicar script de comparação reprodutível.

## Referências

Al Maleh, H. 2009. *A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems*. McGill University.

Arnold, D. N.; Falk, R. S.; Winther, R. Finite element exterior calculus: from Hodge theory to numerical stability. *Bulletin of the AMS*, 47(2), 281–354, 2010. DOI: 10.1090/S0273-0979-10-01278-4.