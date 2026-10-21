# M4 — Refinamento local reproduzível (9×7 → 13×11 → 17×15)

Origem dos números: [GitHub Actions, execução 39591730349](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39591730349), usando a montagem com acumuladores compactos e a comparação com soluções analíticas TM. O artefato da execução contém CSV de espectros, métricas e metadados. **Esses números foram medidos em CI; a execução local do usuário será independente.**

## Geometria e autovalores analíticos

Retângulo PEC com a=2, b=1; funções MLS de grau 1, raio de suporte 3,3 vezes o maior espaçamento da nuvem, quadratura 3×3 Gauss por célula de fundo.

O modo TM11 tem lambda = 1,25 pi² ≈ 12,337005501. A tabela a seguir contém **erros relativos em lambda** e **sobreposições normalizadas do campo curl(H) com o perfil TM11**.

| Nuvem | Nós | DOFs | Erro relativo TM11 | Sobreposição TM11 | Condição(M) | Núcleo numérico |
|---|---:|---:|---:|---:|---:|---:|
| 9×7 | 63 | 126 | 0,000039245 | 0,999962611 | 72373,6 | 5 |
| 13×11 | 143 | 286 | 0,000003146 | 0,999993772 | 285083,9 | 6 |
| 17×15 | 255 | 510 | 0,000002997 | 0,999992794 | 443914,1 | 9 |

Em porcentagem, os erros TM11 são aproximadamente **0,0039245%**, **0,0003146%** e **0,0002997%**. O refinamento melhora o autovalor TM11, mas a sobreposição não é estritamente monotônica.

## Modos superiores — forma e condição física

| Nuvem | Modo | Erro lambda | Sobreposição curl | Divergência relativa | Indicador normal PEC |
|---|---|---:|---:|---:|---:|
| 9×7 | TM21 | 0,00006409 | 0,99993272 | 0,031259 | 0,014696 |
| 13×11 | TM21 | 0,00002959 | 0,99998616 | 0,012313 | 0,003987 |
| 17×15 | TM21 | 0,00001813 | 0,99999301 | 0,013238 | 0,004166 |
| 9×7 | TM12 | 0,00013588 | 0,99950497 | 0,069497 | 0,039969 |
| 13×11 | TM12 | 0,00001294 | 0,99997458 | 0,022206 | 0,010385 |
| 17×15 | TM12 | 0,00000912 | 0,99999494 | 0,012726 | 0,004871 |

**Interpretação:** as três primeiras famílias TM são bem aproximadas por esse protótipo cartesiano nas nuvens maiores. A dimensão do núcleo numérico aumenta de 5 para 9; não se deve inferir automaticamente a dimensão exata do núcleo funcional, nem a ausência de modos espúrios. A variação da condição da massa com refinamento também exige acompanhamento.

O código e a campanha ampliada agora incluem perfis analíticos TM31, TM41, TM22 e outros, permitindo analisar a **degenerescência TM22–TM41** do retângulo 2×1. A tabela acima foi registrada antes dessa ampliação, portanto não atribui números que não foram observados nesse experimento.

## Como repetir

~~~bash
JOBS=2 PROFILE=local bash scripts/run_local_reference.sh
~~~

O relatório e todos os CSV serão salvos em out/m4/local. Consulte [o guia de execução](../implementation/05_local_refinement.md).

## Limitações da reprodução

Esses valores confirmam a aproximação às **soluções analíticas de um guia retangular PEC** com as decisões de implementação adotadas. **Não se pode afirmar que reproduzimos uma tabela específica da dissertação de 2009** sem conferir diretamente a fonte original (domínio, unidades, pesos, nós, fronteira, quadratura e observáveis). O software permanece experimental.
