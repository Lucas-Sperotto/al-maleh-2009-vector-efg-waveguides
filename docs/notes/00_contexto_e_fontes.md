# Contexto e fontes

> Análise inicial migrada do antigo README (linhas 13–26), mantida para consulta e rastreabilidade. As equações são reconstruções de estudo, não transcrições verificadas de Al Maleh. Veja a [referência técnica](../theory/README.md) e o [estado atual](../../README.md). Os delimitadores matemáticos foram adequados ao GitHub.

---

## Análise matemática e proposta de implementação do método de Al Maleh (2009)

**Minha recomendação é reproduzir primeiro as duas formulações de Al Maleh em um repositório próprio, utilizando o projeto NASA TP-3485 como referência de validação.** Depois, poderemos integrar os resultados ao projeto EFGMI-H(curl), no qual você já está investigando funções de Whitney difusas.

Examinei também seus repositórios existentes. Encontrei alguns pontos matemáticos e computacionais que merecem atenção antes de iniciarmos a implementação.

O trabalho de Hussam Al Maleh é especialmente interessante porque aborda um problema fundamental do eletromagnetismo computacional: **como aproximar campos vetoriais por um método sem malha sem introduzir modos próprios espúrios**.

Há, entretanto, uma limitação nesta análise: localizei o PDF original e confirmei o resumo, a organização dos capítulos e os resultados gerais, mas os servidores não permitiram recuperar integralmente as páginas nesta sessão. Assim, as equações abaixo são uma reconstrução matemática do problema e das duas classes de aproximação descritas pelo autor, não uma transcrição verificada de todas as equações da dissertação. Para uma reprodução fiel, será necessário conferir especialmente as expressões dos capítulos 5 e 6. [ResearchGate](https://www.researchgate.net/publication/38437582_A_vector_Element-Free_Galerkin_method_for_waveguide_eigenvalue_problems?utm_source=chatgpt.com)

Documento original: [A Vector Element-Free Galerkin Method for Waveguide Eigenvalue Problems — McGill University](https://mcgill.scholaris.ca/bitstreams/4d79b125-540f-4b38-8c06-722d14eb60bb/download?utm_source=chatgpt.com).

---
