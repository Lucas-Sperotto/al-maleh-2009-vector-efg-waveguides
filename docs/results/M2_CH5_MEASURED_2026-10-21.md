# M2 — Resultados executados frente às Tabelas 2/3 de Al Maleh (2009)

> **Nota posterior (auditoria independente de 2026-10-22):** Com o perfil do Capítulo 4 identificado (spline quártica, base linear, r=1,5h inferido), as 28 entradas das Tabelas 2–3 são reproduzidas a <0,001 pp; as diferenças abaixo refletem o perfil Wendland C²/3,3h usado nesta campanha. Ver [relatório da auditoria](INDEPENDENT_THEORY_CODE_AUDIT.md). O conteúdo abaixo é preservado sem alterações como registro histórico.

**Data: 21/10/2026**. Campanha integral [GitHub Actions #39630223717](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39630223717), executada no commit `ad51043d` (a CI verificou o programa e os dados de fonte sem alterar o solver original de M3/M4).

**Resultado:** 23/23 testes CTest aprovados; todas as **14 combinações** de nuvem e quadratura concluídas, produzindo **28 comparações** com a transcrição do capítulo 5. [Artefato integral M2](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39630223717/artifacts/12358927216) preserva todos os autovalores/autovetores, matrizes, diagnósticos, CSVs originais Hx/Hy, SVGs originais e logs.

## Configuração deste experimento

- Guia PEC 2×1 m; nós `3×5` (ny=3,nx=5), `4×7`, `5×9`, `6×11`; quatro células de integração por retângulo de nós, conforme Figura 9.
- Quadraturas Gauss-Legendre 1×1 a 5×5 pontos por célula; apenas combinações com valores publicados foram executadas.
- MLS linear, Wendland C², raio de suporte 3,3 vezes o maior espaçamento nodal, que são decisões da implementação atual; **não há comprovação ainda de identidade com os parâmetros originais do autor**.
- Dois pareamentos por modo: candidato de λ mais próximo da referência analítica e candidato com maior sobreposição normalizada de `curl(H)` com a forma TM analítica.
- Valores e erros aqui são percentuais **assinados** relativos a λ=k². Diferenças da última coluna são **pontos percentuais**, não erros percentuais relativos entre percentuais.

## Comparação completa

| Modo | Nós autor | Gauss | Autor (%) | Nosso mais próximo (%) | Nosso melhor `curl(H)` (%) | Diferença absoluta do mais próximo (pp) |
|---|---|---|---:|---:|---:|---:|
| TM11 | 3×5 | 1×1 | +23.767 | -12.417 | +24.647 | 36.184 |
| TM11 | 3×5 | 2×2 | +17.199 | +3.377 | +17.887 | 13.822 |
| TM11 | 3×5 | 3×3 | +17.210 | +3.370 | +17.887 | 13.840 |
| TM11 | 3×5 | 4×4 | +17.206 | +3.386 | +17.887 | 13.820 |
| TM11 | 3×5 | 5×5 | +17.206 | +3.387 | +17.887 | 13.819 |
| TM11 | 4×7 | 1×1 | +5.978 | -0.045 | -0.045 | 6.023 |
| TM11 | 4×7 | 2×2 | +3.608 | -0.004 | -0.004 | 3.612 |
| TM11 | 4×7 | 3×3 | +3.763 | +0.006 | +0.006 | 3.757 |
| TM11 | 4×7 | 4×4 | +3.759 | +0.006 | +0.006 | 3.753 |
| TM11 | 4×7 | 5×5 | +3.756 | +0.006 | +0.006 | 3.750 |
| TM11 | 5×9 | 1×1 | +2.805 | +0.221 | +0.221 | 2.584 |
| TM11 | 5×9 | 2×2 | +1.617 | +0.082 | +0.082 | 1.535 |
| TM11 | 6×11 | 1×1 | +1.567 | +0.031 | +0.031 | 1.536 |
| TM11 | 6×11 | 2×2 | +0.839 | +0.011 | +0.011 | 0.828 |
| TM21 | 3×5 | 1×1 | +18.557 | +17.461 | +17.461 | 1.096 |
| TM21 | 3×5 | 2×2 | +11.846 | +12.456 | +12.456 | 0.610 |
| TM21 | 3×5 | 3×3 | +12.103 | +12.470 | +12.470 | 0.367 |
| TM21 | 3×5 | 4×4 | +12.089 | +12.470 | +12.470 | 0.381 |
| TM21 | 3×5 | 5×5 | +12.085 | +12.470 | +12.470 | 0.385 |
| TM21 | 4×7 | 1×1 | +4.508 | -0.042 | -0.042 | 4.550 |
| TM21 | 4×7 | 2×2 | +2.038 | -0.049 | -0.049 | 2.087 |
| TM21 | 4×7 | 3×3 | +2.174 | -0.030 | -0.030 | 2.204 |
| TM21 | 4×7 | 4×4 | +2.164 | -0.030 | -0.030 | 2.194 |
| TM21 | 4×7 | 5×5 | +2.161 | -0.030 | -0.030 | 2.191 |
| TM21 | 5×9 | 1×1 | +2.313 | +0.143 | +0.143 | 2.170 |
| TM21 | 5×9 | 2×2 | +1.071 | +0.044 | +0.044 | 1.027 |
| TM21 | 6×11 | 1×1 | +1.307 | +0.024 | +0.024 | 1.283 |
| TM21 | 6×11 | 2×2 | +0.539 | +0.004 | +0.004 | 0.535 |

## Achados

1. **TM11 na nuvem 3×5 com Gauss 2×2:** o autor publica 17.199% e o candidato de maior correlação de forma da implementação atual apresenta +17.887%, diferença de 0.688 ponto percentual. O candidato espectral mais próximo, porém, está em +3.377%. Isso demonstra que a seleção apenas pelo autovalor pode estar escolhendo outro campo.
2. **TM11 na nuvem 4×7 com Gauss 2×2:** o autor informa +3.608%, contra -0.004% no nosso solver. Não é reprodução numérica da Tabela 2, mesmo com formas matemáticas compatíveis.
3. **TM21 na nuvem 6×11 com Gauss 2×2:** o autor publica +0.539%, e o nosso experimento apresenta +0.004%. A tendência de erro menor em malhas mais refinadas permanece visível, mas os valores não são numericamente iguais.
4. **Modo espúrio da Figura 15:** o autor informa k²≈12,1456 para uma malha 6×11 com Gauss 2×2. A figura é uma referência histórica para uma verificação de forma física; não se reivindica reprodução desse autovalor individual ou de seu autovetor nesta primeira campanha.

## Conclusão e pendências

Há consistência estrutural com (5-1)–(5-11), em especial com a forma generalizada [S−k²T]h=0, mas **não** há correspondência quantitativa geral com as Tabelas 2 e 3. Diferenças de MLS, peso/suporte, imposição efetiva de fronteira, integração e identificação modal precisam ser confrontadas com os detalhes completos do capítulo 4 e do capítulo 5 antes de ajustar a implementação.

Os SVGs são **desenhos originais dos autovetores desta campanha**, não cópias das Figuras 10–13. O modo selecionado para desenho é o de melhor correlação `curl(H)`, conforme o relatório bruto; sua qualidade física não é certificada por essa correlação isoladamente.

Referência para números do autor: imagens fornecidas do Capítulo 5 de Al Maleh (2009), Tabelas 2 e 3; transcrição verificável no repositório em `docs/results/m2_author_ch5_tables_2_3.csv`. Nenhum arquivo histórico M6a–M6f foi reescrito.
