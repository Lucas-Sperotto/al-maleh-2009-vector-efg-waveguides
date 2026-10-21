# M6d — Resultados do diagnóstico do quase-núcleo gradiente

Campanha executada em 21/10/2026 via [GitHub Actions
run 39620895466](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39620895466),
na branch feat/m6d-gradient-kernel-audit, baseada em M6c.
19/19 CTests passaram. A verificação independente Python reabriu
as matrizes e o espectro completo das três nuvens e das duas
configurações adicionais de posto.

## Fração da energia projetada no espaço gradiente

A fração f_grad=||P_grad v||_M²/||v||_M² mede a energia do
autovetor que está no espaço formado pelas melhores aproximações
L2 vetoriais dos gradientes MLS escalares.

| Nuvem | DOFs vetoriais | Posto gradiente | Positivos abaixo TM11 | Fração gradiente média | Mediana | f_grad>=90% |
|---|---:|---:|---:|---:|---:|---:|
| 9×7 | 126 | 62 | 49 | 0,916422 | 0,998846 | 43 |
| 13×11 | 286 | 142 | 94 | 0,970516 | 0,999563 | 87 |
| 17×15 | 510 | 254 | 144 | 0,982869 | 0,999712 | 139 |

Os postos 62/142/254 são iguais a n-1 para n=63/143/255,
o que é consistente com a direção escalar constante de gradiente nulo.
Esse acordo **não prova** por si só a existência de um complexo discreto
compatível com Maxwell.

A elevada fração gradiente dos valores positivos abaixo de TM11
constitui evidência numérica direta de sua associação com o
quase-núcleo gradiente no EFG cartesiano. Os demais candidatos
precisam de análise modal individual.

## Modos selecionados

| Nuvem | Fração gradiente TM11 | Fração gradiente coluna A do dubleto | Coluna B do dubleto |
|---|---:|---:|---:|
| 9×7 | aprox. 0,000000 | 0,001156 | 0,001445 |
| 13×11 | aprox. 0,000001 | 0,000004 | 0,000424 |
| 17×15 | aprox. 0,000010 | 0,000003 | 0,016017 |

**Distinção importante:** a coluna B do dubleto fino foi
reprovada no M6b por divergência de aproximadamente 0,415325,
mas sua fração gradiente projetada foi somente 0,016017
(1,60%). Portanto a falha do dubleto **não foi explicada**
por simples pertença majoritária ao espaço gradiente projetado.
O defeito do dubleto e as raízes pequenas podem necessitar
mecanismos ou restrições discretas adicionais.

## Verificação independente e estabilidade

Todas as reaberturas passaram e reconstruíram
G = H^T M^-1 H, Q^T M Q, as participações
gradientes por autovetor e os autovalores vetoriais originais
contra os arquivos da campanha M6a.
Os piores erros reportados na nuvem 17×15 foram:

- Erro relativo da reconstrução da Gramiana: 1,757e-15;
- Erro máximo de Q^T M Q - I: 3,755e-12;
- Diferença máxima de participações gravadas versus recalculadas: 2,717e-12.

A varredura adicional de cortes relativos
tau=1e-6, 1e-8 e 1e-10 manteve o posto 254
na nuvem 17×15.

## Limitações

A projeção é estritamente diagnóstica: o método M6d **não filtrou**
modos, **não modificou K/M** e **não alterou** a campanha
M6a nem seus manifestos SHA-256. A quase-coincidência
numérica entre gradientes e valores próprios pequenos
não equivale ainda a uma prova de convergência de um
espaço de Whitney/H(curl) ou de uma sequência exata discreta.

Para investigação seguinte: experimentar uma formulação
mista ou um espaço gradiente-rotacional de compatibilidade
demonstrável, com provas de posto, consistência de contorno
e comparações FEM×EFG nos modos simples e degenerados.

## Reprodução

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel 2
    /usr/bin/ctest --test-dir build --output-on-failure
    python3 scripts/m6d_gradient_campaign.py --binary build/gradient_kernel_audit --out out/m6d
    PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 python3 scripts/verify_m6d_gradient.py --campaign out/m6d

Dados completos da campanha M6d estão disponíveis como
artefato no GitHub Actions. Ver também
docs/implementation/08_gradient_kernel_audit.md.
