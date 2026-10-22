# M6b — Extremos invariantes do subespaço TM22/TM41

> **Nota posterior (auditoria independente de 2026-10-22):** A falha de divergência 0,415 do dubleto 17×15 foi explicada: hibridização com o modo espúrio quase-gradiente 191 (λ=49,204, a 0,31% do dubleto), que contém ≈1,6% da energia física. A projeção L² do dubleto analítico no mesmo espaço tem div ≤ 0,0045; não é piso de aproximação. Ver [relatório da auditoria](../INDEPENDENT_THEORY_CODE_AUDIT.md). O conteúdo abaixo é preservado sem alterações como registro histórico.

Resultados computados em 21/10/2026 com o executável
`constraint_subspace` e publicados em
[GitHub Actions run 39602371572](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/actions/runs/39602371572),
artefato **m6b-constraint-subspaces**. A forma discreta K/M e os critérios
históricos **não foram alterados**. Os resultados provêm de nova integração
das formas D/B, não de uma cópia textual da tabela antiga.

## Resultados dos três refinamentos

| Nuvem EFG | DOFs | Positivos abaixo de TM11 | Div min | Div max | Normal min | Normal max |
|---|---:|---:|---:|---:|---:|---:|
| 9×7 | 126 | 49 | 0,086237391 | 0,090135263 | 0,061217306 | 0,070964743 |
| 13×11 | 286 | 94 | 0,006457085 | 0,065223689 | 0,003712919 | 0,040078167 |
| 17×15 | 510 | 144 | 0,007521385 | **0,415324720** | 0,002980709 | 0,243585507 |

As métricas são calculadas pelas razões extremas de Rayleigh projetadas
no **espaço de dimensão dois** selecionado no M6a, utilizando norma L2
da massa e `lambda_ref=5*pi^2` para todo o dubleto.

**Conclusão objetiva:** a divergência máxima do subespaço EFG 17×15
continua acima do limiar exploratório 0,25, mesmo depois de remover a
dependência da escolha da base. Trata-se de uma falha física numérica
não explicada pelo posicionamento dos autovetores na diagonalização.

### O que este experimento NÃO demonstrou

- Não foi alterada a matriz de rigidez nem imposta restrição solenoidal.
- Não se explicou por que o refinamento fino degrada o subespaço degenerado.
- Não foram classificados individualmente todos os 49/94/144 candidatos.
- A proximidade da condição normal máxima 0,243586 ao limite 0,25 **não**
  valida o subespaço: a divergência máxima 0,415325 reprova.
- A comparação FEM×EFG da campanha M6a é um baseline parcial, não uma
  certificação de inexistência de espúrios.

### Experimentos seguintes

Montar a projeção do subespaço gradiente/quase-núcleo e uma formulação
solenoidal ou mista fundamentada na identidade diferencial contínua;
repetir a campanha em branch separada e comparar os espectros inteiros
com os dados históricos preservados, nunca apagando a falha de M6a.

## Execução local

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
python3 scripts/m6b_campaign.py --binary build/constraint_subspace --out out/m6b
~~~

Saídas: `out/m6b/summary.csv`, `out/m6b/REPORT.md`, e CSVs
de subespaços e espectros inferiores a TM11 nas três subpastas.

Os resultados são evidência de execução na CI, e não de ensaio de
desempenho no Stargate.
