# Guia de execução local — reprodução numérica do guia retangular

**Estado:** protocolo experimental para a discretização vetorial EFG cartesiana do repositório. O objetivo imediato é **confrontar modos TM e frequências de corte com os valores analíticos do guia retangular PEC de lados 2×1**. Isso **não é ainda** o cotejamento das mesmas tabelas, parâmetros e figuras da dissertação de Al Maleh (2009).

## Preparação (Ubuntu, Ubuntu no WSL ou Stargate)

~~~bash
sudo apt update
sudo apt install -y git build-essential cmake libeigen3-dev python3
git clone https://github.com/Lucas-Sperotto/al-maleh-2009-vector-efg-waveguides.git
cd al-maleh-2009-vector-efg-waveguides
~~~

Se já existe a pasta, apenas acesse-a e execute git pull após atualizar para o commit principal mais recente.

## Campanha completa no computador local (recomendada)

~~~bash
JOBS=2 PROFILE=local bash scripts/run_local_reference.sh
~~~

O script configura CMake em Release, compila usando duas tarefas simultâneas, executa os testes CTest e roda as seguintes nuvens nodais:

| Configuração | Número de nós | DOFs vetoriais | Grau MLS | Gauss |
|---|---:|---:|---:|---:|
| local_9x7 | 63 | 126 | 1 | 3 |
| local_13x11 | 143 | 286 | 1 | 3 |
| local_17x15 | 255 | 510 | 1 | 3 |

Todas usam domínio 2×1 e raio de suporte igual a 3,3 vezes o maior espaçamento da nuvem. São resolvidas as matrizes generalizadas K h = lambda M h e medidas a divergência de H, o fluxo normal nas paredes e a correlação entre curl(H) e os perfis analíticos de Ez.

Artefatos em out/m4/local:

- REPORT.md: tabela sintética dos erros, sobreposições, núcleo e condicionamento;
- tm_matches_all.csv: comparação por TM11, TM21, TM31, TM12, TM22, TM41 etc., incluindo pares degenerados;
- runs.csv: status de cada execução e falhas numéricas, caso existam;
- local_9x7/, local_13x11/, local_17x15/: tabelas originais e stdout/stderr de cada caso.

Também é possível executar somente uma configuração:

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure

./build/tm_validation 17 15 1 3.3 3 out/m4/manual_17x15
~~~

## Perfil para o Stargate

O master antigo do Stargate tem recursos mais limitados; para iniciar sem executar a matriz de 510 DOFs, use:

~~~bash
JOBS=1 PROFILE=stargate bash scripts/run_local_reference.sh
~~~

Esse perfil executa inicialmente apenas 9×7 e 13×11 (126 e 286 DOFs). Para uma campanha mais refinada, usar PROFILE=local depois de examinar RAM, troca (swap), tempo de execução e condicionamento. O solucionador atual converte matrizes para **densas** e sua implementação limita o tamanho a 600 DOFs; **não** tentar aumentar indefinidamente a resolução sem mudar o autossolucionador.

O desempenho no master Xeon legado não deve ser estimado a partir da execução no GitHub Actions. O CI é um experimento independente, não uma medida de velocidade do Stargate.

## Valores analíticos esperados

Com a=2, b=1 e lambda=kc²:

| Modo | lambda analítico | kc analítico |
|---|---:|---:|
| TM11 | 12,337005501 | 3,5124074 |
| TM21 | 19,739208802 | 4,4428829 |
| TM31 | 32,076214304 | 5,6635867 |
| TM12 | 41,945818705 | 6,4765592 |
| TM22 | 49,348022005 | 7,0248147 |
| TM41 | 49,348022005 | 7,0248147 |

TM22 e TM41 são **degenerados** nesse domínio; identificar dois autovetores isolados pelo índice pode ser enganoso. A validação M4 compara projeções dos subespaços espectrais e mede simultaneamente erro no autovalor, similaridade de curl(H), divergência, fluxo normal e resíduo.

Para converter kc em frequência de corte física (Hz), é preciso especificar as dimensões **em metros** e as propriedades do meio; o modelo atual usa a=2 e b=1 em unidades consistentes e trata kc em unidades inversas de comprimento. Não atribuir frequências reais sem essa conversão.

## Resultado anterior, validado em GitHub Actions

Em uma campanha com 9×7, 13×11 e 17×15 nós, todos os casos foram executados. As sobreposições de TM11 reportadas foram maiores que 0,9999. Para TM21, os erros relativos em lambda foram aproximadamente 6,409e-5, 2,959e-5 e 1,813e-5; para TM12, 1,3588e-4, 1,294e-5 e 9,12e-6. A condição de M depende fortemente dos suportes e deve ser monitorada.

Esses resultados foram obtidos no runner de CI, não na máquina pessoal nem no Stargate. A campanha revisada também foi ampliada para incluir TM31, TM41 e outros modos, e os valores desta ampliação devem ser consultados em seus próprios artefatos de execução.

## Como aproximar a dissertação original

A dissertação apresenta duas representações do campo (cartesiana e rotacional/irrotacional), testadas em guia retangular PEC. Para afirmar que **reproduzimos uma figura ou tabela específica**, é obrigatório conhecer seus valores de dimensões, densidade de nós, funções peso, bases MLS, imposição de fronteiras, quadratura e espectro originalmente publicado. Esses elementos ainda **não foram cotejados integralmente** com a fonte primária; os números acima são uma validação contra **soluções analíticas**, não uma reprodução fiel comprovada de todas as tabelas da dissertação. 
