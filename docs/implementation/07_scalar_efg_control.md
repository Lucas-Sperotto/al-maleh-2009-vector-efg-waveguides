# M6c — Controle independente de EFG escalar para modos TM

> **Experimento paralelo de diagnóstico.** Este programa **não implementa**
> a formulação vetorial de Al Maleh e **não corrige** a contaminação
> espectral detectada nas campanhas M6a e M6b.
> A finalidade é testar o impacto da escolha do espaço aproximante
> utilizando as mesmas funções MLS, geometria e quadratura.

## Formulação

Para o guia retangular PEC `Ω=(0,2)×(0,1)`, o potencial longitudinal TM
satisfaz:

```text
-Δ Ez = λ Ez            em Ω
Ez = 0                  em ∂Ω
λ = kc²
```

A aproximação `Ez_h(x)=sum_i phi_i(x) c_i` usa MLS linear e peso Wendland
C². As matrizes são montadas independentemente das matrizes vetoriais:

```text
S_ij = ∫Ω ∇phi_i · ∇phi_j dΩ
M_ij = ∫Ω phi_i phi_j dΩ
```

**MLS não é interpolatória.** Portanto, não basta zerar coeficientes de
nós na fronteira. Avaliamos as funções de forma em todos os nós de
fronteira da nuvem para formar a matriz:

```text
C_(q,i) = phi_i(x_q)        x_q em ∂Ω
C c = 0
```

Uma decomposição SVD `C = U Σ Vᵀ` fornece uma base `Q` para o núcleo
numérico de C. O autoproblema admissível é:

```text
(Qᵀ S Q) y = λ (Qᵀ M Q) y
c = Q y
```

O código usa `Eigen::GeneralizedSelfAdjointEigenSolver` com checagens
de positividade definida, condição da massa, resíduos **projetados**
no subespaço admissível, ortogonalidade M e restrições Cc.

É incorreto exigir `Sc−λMc=0` no espaço irrestrito, já que a solução
com restrições de contorno admite forças do multiplicador de Lagrange.
O teste relevante é `Qᵀ(Sc−λMc)=0`.

## Campo magnético de controle

```text
H_t = (-∂Ez/∂y, ∂Ez/∂x)
div H_t = 0
curl H_t = ΔEz
```

Por construção, a divergência é nula em aritmética matemática para
funções MLS com derivadas mistas consistentes. Portanto, **o zero
de divergência nesta formulação não é uma validação independente da
discretização vetorial cartesiana.**

Também medimos por quadratura segmentada:

- Sobreposição L2 de `Ez_h` com `sin(mπx/2)sin(nπy)`;
- Traço de `Ez` nas paredes, **inclusive entre os pontos de colocalização**;
- Fluxo normal de `H_t` nas paredes.

**Limitação:** Cc=0 impõe `Ez=0` exatamente apenas nos pontos de
colocalização; a condição PEC pode ser violada entre eles. O estudo da
convergência dessas violações é indispensável.

## Dados preservados

O executável exporta:

- `eigenvalues.csv`: espectro reduzido integral, números de corte
  e resíduos projetados / de colocalização;
- `eigenvectors.mtx`: todos os autovetores escalares completos
  nos coeficientes originais (não apenas a solução reduzida);
- `stiffness_matrix.mtx` e `mass_matrix.mtx`: matrizes escalares
  originais do EFG;
- `boundary_collocation.mtx` e `boundary_nullspace.mtx`: C e Q;
- `dof_map.csv`: coordenadas originais e marcação da fronteira;
- `tm_reference.csv`: comparações exploratórias por proximidade
  espectral, correlações e resíduos físicos de contorno;
- `run_metadata.txt`: geometria, raio de suporte, posto,
  grau de liberdade, condicionamento, integrações etc.

Todos os `.mtx` utilizam **Matrix Market array real general** em
precisão textual de 17 algarismos significativos. O índice do espectro
começa em zero; a coluna matemática no Matrix Market começa em um.

É proibido rotular o par TM22/TM41 como resolvido com base apenas em
autovalores próximos. A identificação de um subespaço degenerado exige
análise conjunta por ângulos principais.

## Execução

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
/usr/bin/ctest --test-dir build --output-on-failure
python3 scripts/m6c_scalar_campaign.py --binary build/scalar_efg_control --out out/m6c
python3 scripts/verify_scalar_control.py --campaign out/m6c
~~~

O teste `verify_scalar_control.py` reabre todas as matrizes e autovetores,
compara dimensões, ordenação, o produto de Gram `VᵀMV`, `CV`,
`CQ`, `QᵀQ` e os resíduos do autoproblema projetado.

A comparação M6a FEM × EFG vetorial **não é sobrescrita**: M6c cria
uma terceira solução, potencialmente útil para identificar efeitos
introduzidos exclusivamente pelo espaço EFG vetorial cartesiano.
