# Arquitetura e pesquisa

> Análise inicial migrada do antigo README (linhas 1166–1355), mantida para consulta e rastreabilidade. As equações são reconstruções de estudo, não transcrições verificadas de Al Maleh. Veja a [referência técnica](../theory/05_arquitetura_reuso_rastreabilidade.md) e o [estado atual](../../README.md). Os delimitadores matemáticos foram adequados ao GitHub.

---

## 8. Onde colocar essa implementação?

Agora podemos discutir os repositórios.

Examinei quatro projetos seus relacionados a esse problema.

| Repositório | Situação | Minha avaliação |
|---|---|---|
| `fem-eigen-em-2d3d-whitney` | Reprodução do NASA TP-3485, com os 14 casos documentados | Excelente referência para validação |
| `efgmi-hcurl-diffuse-whitney` | Investigação de funções vetoriais Whitney difusas | Excelente para comparação e evolução |
| `simple-waveguide-efg-fem-solver` | Implementação escalar FEM/EFGM | Pode fornecer componentes reutilizáveis após correções |
| `efg-3d-electromagnetics-2006` | Reprodução do trabalho de Parreira et al. | Menos relacionado diretamente ao autoproblema vetorial |

#### Minha proposta

Criaria um repositório:

**`al-maleh-2009-vector-efg-waveguides`**

Ele teria um objetivo bem definido:

> Reproduzir e avaliar numericamente as formulações vetoriais do método Element-Free Galerkin propostas por Hussam Al Maleh (2009) para problemas de autovalores em guias de onda, com ênfase na identificação dos modos físicos TM, na análise dos modos espúrios e na comparação com soluções analíticas e pelo método dos elementos finitos.

A organização seria:

```text
al-maleh-2009-vector-efg-waveguides/
│
├── CMakeLists.txt
├── README.md
├── AGENTS.md
├── CITATION.cff
│
├── include/
│   ├── geometry/
│   ├── mls/
│   ├── vector_basis/
│   ├── assembly/
│   ├── boundary/
│   ├── eigen/
│   └── validation/
│
├── src/
│   ├── geometry/
│   ├── mls/
│   ├── vector_basis/
│   │   ├── cartesian.cpp
│   │   └── rotational_irrotational.cpp
│   ├── assembly/
│   ├── boundary/
│   └── eigen/
│
├── apps/
│   ├── scalar_reference.cpp
│   ├── cartesian_waveguide.cpp
│   ├── rotational_waveguide.cpp
│   └── compare_methods.cpp
│
├── tests/
│   ├── test_mls.cpp
│   ├── test_mls_derivatives.cpp
│   ├── test_vector_basis.cpp
│   ├── test_matrix_assembly.cpp
│   ├── test_rectangular_modes.cpp
│   └── test_spurious_modes.cpp
│
├── docs/
│   ├── theory/
│   ├── equations/
│   ├── implementation/
│   ├── validation/
│   ├── comparison/
│   └── references/
│
├── benchmarks/
│   ├── rectangular/
│   └── nasa_tp3485/
│
├── scripts/
│
└── results/
```

Manteria C++17/CMake, seguindo seus projetos de computação científica, com Eigen ou LAPACK para as operações lineares.

Também adotaria a mesma filosofia de rastreabilidade que utilizamos no projeto NASA:

$$
\text{Equação da dissertação}
\longrightarrow
\text{Formulação discreta}
\longrightarrow
\text{Função C++}
\longrightarrow
\text{Teste}
\longrightarrow
\text{Resultado}.
$$

Não misturaria inicialmente o código de Al Maleh ao núcleo de reprodução da NASA.

O repositório `fem-eigen-em-2d3d-whitney` já possui uma organização muito clara em torno do trabalho de Jin e dos 14 casos originais. Introduzir ali outra dissertação como método principal dificultaria a separação da autoria e dos resultados científicos.

**O NASA TP-3485 deve funcionar como referência independente, não como recipiente de todas as implementações.**

---

## 9. Planejamento da implementação

Minha proposta seria trabalhar em seis etapas, com testes e commits pequenos.

**Etapa 1 — Fundamentos matemáticos e MLS.** Reconstruir a formulação da dissertação, corrigir a implementação de MLS, implementar primeiras e segundas derivadas e executar testes de reprodução polinomial. Nenhuma matriz de autovalores seria montada antes desses testes.

**Etapa 2 — Formulação cartesiana.** Implementar as funções vetoriais, montar as matrizes $\mathbf K$ e $\mathbf M$ e resolver o primeiro autoproblema. A montagem deve passar por testes de simetria, integração e comparação com integrais conhecidas.

**Etapa 3 — Reprodução dos resultados físicos.** Utilizar o guia retangular, comparar os modos TM com a solução analítica e gerar gráficos dos campos magnéticos. Refinar a nuvem nodal e estudar a convergência.

**Etapa 4 — Investigação dos modos espúrios.** Analisar o núcleo de $\mathbf K$, a condição de divergência, os resíduos espectrais e a sensibilidade aos parâmetros MLS. Reproduzir, na medida do possível, as condições numéricas da dissertação.

**Etapa 5 — Formulação rotacional/irrotacional.** Depois de conferir as expressões originais do capítulo 6, implementar a segunda representação vetorial, suas matrizes e o tratamento das dependências lineares.

**Etapa 6 — Comparação com FEM e EFGMI-H(curl).** Executar benchmarks equivalentes e investigar quais propriedades matemáticas e numéricas são preservadas por cada aproximação.

Ao final teremos três caminhos independentes de resolução vetorial para comparação: FEM com Whitney/Nédélec, EFG cartesiano de Al Maleh e EFGMI com Whitney difusas. A quarta implementação, rotacional/irrotacional de Al Maleh, ampliará esse conjunto.

---

## 10. O que eu investigaria além da dissertação

O resultado mais interessante não seria apenas reproduzir o trabalho de 2009.

Seria utilizar essa reprodução para estudar por que algumas discretizações sem malha são capazes de aproximar corretamente os campos eletromagnéticos, enquanto outras introduzem espectros contaminados.

Em particular, considero três questões matemáticas importantes.

**Primeira: a estrutura do núcleo do operador.** Precisamos determinar se a dimensão elevada do núcleo observada no projeto EFGMI-H(curl) está associada à redundância das funções de Whitney difusas, à seleção de pares virtuais ou à falta de compatibilidade discreta.

**Segunda: o efeito das derivadas MLS.** A correção das primeiras derivadas e o cálculo consistente das segundas derivadas podem afetar substancialmente os resultados numéricos. Vale investigar quanto do erro vem da representação vetorial e quanto vem da própria aproximação MLS.

**Terceira: a compatibilidade entre os espaços discretos.** A aproximação cartesiana utiliza diretamente duas famílias de funções escalares; a formulação rotacional utiliza derivadas dessas funções; e sua abordagem Whitney difusa constrói uma família vetorial por combinações antissimétricas. Podemos investigar a capacidade de cada espaço representar gradientes, campos solenoidais e autofunções de Maxwell.

Isso pode produzir uma comparação científica mais interessante do que simplesmente apresentar tabelas de erros de autovalores.

---

### 11. Minha conclusão

Eu **não começaria implementando diretamente a formulação rotacional/irrotacional**.

Começaria corrigindo e validando o núcleo MLS e depois desenvolveria a aproximação cartesiana de Al Maleh. Assim, conseguiremos reproduzir o comportamento principal descrito na dissertação antes de introduzir uma base vetorial mais complexa.

Também não considero adequado modificar agora o núcleo do repositório NASA para comportar esse método. Ele já é uma boa referência FEM, e preservar essa independência permitirá comparações mais confiáveis.

O caminho que recomendo é:

$$
\boxed{
\text{MLS validado}
\longrightarrow
\text{EFG cartesiano}
\longrightarrow
\text{Modos TM}
\longrightarrow
\text{Espúrios}
\longrightarrow
\text{EFG rotacional}
}
$$

Depois:

$$
\boxed{
\text{Al Maleh}
\;\;\longleftrightarrow\;\;
\text{NASA TP-3485}
\;\;\longleftrightarrow\;\;
\text{EFGMI-H(curl)}
}
$$

**Há uma oportunidade de pesquisa bastante concreta aqui:** investigar a estrutura e a compatibilidade dos espaços vetoriais sem malha, utilizando uma implementação independente de Al Maleh como referência para avaliar e aperfeiçoar suas funções de Whitney difusas.

Antes de escrevermos o código definitivo, a prioridade será obter a conferência integral das equações dos capítulos 5 e 6. Isso evitará transformar uma implementação matematicamente plausível em uma reprodução que atribua ao autor escolhas que ele não fez.

Não alterei nenhum dos seus repositórios nesta análise. A decisão arquitetural que considero mais sólida é criar `al-maleh-2009-vector-efg-waveguides` e reservar o NASA TP-3485 para validação cruzada.


---
