# M1 — Núcleo MLS consistente (protótipo de pesquisa)

Status: **implementado e testado localmente**, mas **não é ainda um solver de guias de ondas**. O peso Wendland C2, as bases e os suportes são **decisões próprias**, ainda não confrontadas página a página com Al Maleh (2009).

## Origem e algoritmo

Adaptação independente do núcleo de primeiras derivadas em `Lucas-Sperotto/efgmi-hcurl-diffuse-whitney/src/mls/mls_evaluator.cpp` (MIT, copyright 2026 Sperotto). O novo núcleo amplia a derivação para segundas derivadas sem inversão explícita de matrizes.

A matriz de momentos é `A = sum_I w_I p_I p_I^T`; para uma base fixa de polinômios centralizados/escalados no domínio, defina `g=A^{-1} p(x)`. Como `phi_I = w_I p_I^T g`, resolve-se pela mesma fatoração Cholesky:

```text
A gx  = px  - Ax g
A gy  = py  - Ay g
A gxx = pxx - Axx g - 2 Ax gx
A gxy = pxy - Axy g - Ax gy - Ay gx
A gyy = pyy - Ayy g - 2 Ay gy
```

As derivadas das funções são obtidas pela regra do produto (incluindo derivadas dos pesos). Não formar `inverse(A)`.

## Peso radial

`w(q) = (1-q)^4 (4q+1)` para `0<=q<1`; `w(q)=0` fora; `q=|x-x_I|/h_I`. A implementação avalia gradiente e Hessiana analiticamente, incluindo o limite no centro do suporte.

## Interface e estabilidade

- `MLSEvaluator(nodes, degree=1|2, max_condition=1e12).evaluate(point)` retorna `phi, dx, dy, dxx, dxy, dyy`, lista de nós ativos e condição da matriz de momentos.
- A origem e a escala polinomial são definidas uma vez pela caixa delimitadora da nuvem, não dependem do ponto de avaliação; isso mantém as expressões de derivadas válidas.
- Suportes são circulares e podem variar por nó; vizinhos são encontrados por busca linear (otimização posterior).
- Exceções explícitas para suportes insuficientes, matriz não definida positiva, condição elevada e valores não finitos.
- Não há regularização silenciosa nem imposição de condições PEC nesta etapa.

## Compilação e diagnóstico

Ubuntu/Debian: `sudo apt install -y build-essential cmake libeigen3-dev`.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/mls_probe > mls_probe.csv
```

`mls_probe` usa a base quadrática de 25 nós em um exemplo pequeno e imprime valores e derivadas por nó ativo. Testes verificam partição da unidade, derivadas dessa identidade, reprodução exata de polinômios lineares/quadráticos e suas derivadas, diferenças finitas de `phi` e `grad phi`, singularidade e avaliação exatamente em um nó.

O teste local de 20/10/2026 foi compilado com g++ e Eigen headers e retornou PASS. A validação em CI/repositório remoto deverá ser verificada independentemente.

## Limitações e próximo bloqueio

A reprodução polinomial MLS *não* valida o autoproblema vetorial de Maxwell. Faltam quadratura em domínio, discretização cartesiana de `H_t`, fronteira PEC, montagem K/M e solução espectral. Aumentar o suporte para garantir unisolvência não substitui análise de condicionamento. Antes de atribuir essas escolhas à dissertação, conferir texto original e registrar diferenças.


## Auditoria pós-integração — identificação bibliográfica pendente

A implementação e os testes do M1 são numericamente consistentes, mas o peso Wendland C², a base polinomial e o fator de suporte usados **não foram ainda atribuídos à dissertação de Al Maleh (2009)**. As imagens verificadas pelo projeto cobrem os capítulos 5 e 6; o capítulo 4 original, com os parâmetros da construção MLS, ainda não foi cotejado na íntegra. A [auditoria científica pós-merges](../results/POST_MERGE_SCIENTIFIC_AUDIT_2026-10-21.md) distingue explicitamente implementação atual e valores publicados. A [issue #17](https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides/issues/17) especifica como comparar peso, suporte, derivadas e solver sem mudar silenciosamente os experimentos M2–M6.

A [primeira varredura controlada de suporte](../results/MLS_SUPPORT_SENSITIVITY_MEASURED_2026-10-21.md) utilizou os fatores 3,0/3,3/3,6, mantendo MLS linear e Wendland C², e demonstrou estabilidade parcial do M2, mas possível troca de modo no M5. **É um estudo da nossa implementação, não um ajuste dos parâmetros da tese.**

Reprodução no MyHoneyBee: use /usr/bin/ctest caso um lançador Python em ~/.local/bin/ctest apresente erro do módulo cmake.

## Auditoria independente de 2026-10-22 — Capítulo 4 cotejado e base local

O PDF completo da dissertação foi consultado (SHA-256 `ba356899…c506`).
A Eq. (4-2) define o peso **spline quártica** w(l/r)=1−6(l/r)²+8(l/r)³−3(l/r)⁴,
com l=‖s−x_g‖ (4-1); a Tabela 1 lista bases linear e quadrática e o
texto **não informa o raio r**. Mudanças desta auditoria:

- `MLSWeight::QuarticSpline` implementa a Eq. (4-2) com gradiente
  −12(1−q)²/r²·(x−x_g) e Hessiana radial 24(1−q)/(r³l); o padrão continua
  `WendlandC2` (escolha própria do projeto, nunca atribuída ao autor).
- **Correção de condicionamento:** a base polinomial passou a ser
  deslocada para o ponto de avaliação e escalada pelo maior raio. Pela
  invariância afim da MLS, valores e derivadas não mudam em aritmética
  exata (espectros arquivados M6a reproduzidos com Δλ/λ ≤ 9,2e-12 na faixa
  física), mas cond(A) deixa de crescer como (L/h)^{2·grau}: em grau 2,
  f=3,3, passou de 3,9e3/8,0e4/1,5e6 (5×3/9×5/17×9) para ≈1,3e3 em todas.
  `condition_number` agora mede apenas a geometria local dos vizinhos.
- Novos testes: estudo de passo das diferenças finitas (regime O(h²) e de
  arredondamento), cantos e bordas da nuvem retangular 5×3 com r=1,5h e
  3,3h, condicionamento independente do refino (falha com a base antiga) e
  consistência de 2ª ordem: grau 2 reproduz Δ(x²+y²)=4 a ~1e-11; grau 1
  tem defeito O(1) **independente de h** (mecanismo da não convergência do
  M5 do autor).

Pela identificação inversa contra a Tabela 4, o perfil do autor é
spline quártica, base linear, r=1,5h ([relatório](../results/INDEPENDENT_THEORY_CODE_AUDIT.md),
[campanha](../results/independent_audit_2026-10-22/REPORT.md)). Esse
raio é **inferência**, não valor impresso na tese.
