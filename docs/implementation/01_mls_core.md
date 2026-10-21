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
