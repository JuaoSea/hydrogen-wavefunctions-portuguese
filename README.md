# 🚧 Atenção 🚧

Esse código foi traduzido do inglês:

Funções de onda do átomo de hidrogênio - PT-BR
https://github.com/ssebastianmag/hydrogen-wavefunctions


# Funções de onda do hidrogênio e gráficos de densidade eletrônica

Modelagem e visualização de funções de onda de átomos de hidrogênio e densidade de probabilidade de elétrons.

* Python 3.11.4
* Matplotlib 3.7.2
* Seaborn 0.12.2
* NumPy 1.25.2
* SciPy 1.11.1

### 1. Mecânica Quântica e Sistemas Atômicos: Uma breve visão geral

A mecânica quântica (QM) é a teoria fundamental da física que fornece uma
descrição das propriedades físicas da natureza na escala de átomos e partículas subatômicas.
Ao contrário da mecânica clássica, que descreve fenômenos macroscópicos, a QM aborda o comportamento
da matéria e da energia no nível quântico (discreto).

O átomo de hidrogênio é especialmente significativo, pois é o átomo mais simples, contendo apenas um elétron.
Sua função de onda pode ser tratada analiticamente, fornecendo insights profundos sobre a natureza dos sistemas quânticos.

<br>

<p align='center'>
  <img src='https://github.com/ssebastianmag/hydrogen-wavefunctions/blob/main/img/hydrogen_probability_densities.png' width=85% />
</p>
<p align='center'>
    <i>Densidade de probabilidade eletrônica para orbitais de átomos de hidrogênio mostrados como seções transversais</i>
</p>

---

#### 1.1 Funções de onda

Uma função de onda, frequentemente denotada como ($\psi$), representa o estado quântico de uma partícula em um sistema.
Ela fornece informações sobre a amplitude de probabilidade dos estados de posição e momento da partícula. 



#### 1.2 Densidade eletrônica | Densidade de probabilidade

A magnitude quadrada da função de onda $|\psi|^2$, fornece a densidade de probabilidade para
a posição da partícula no espaço. Para um elétron em um átomo, ela descreve a distribuição espacial
da probabilidade de localizar o elétron.




#### 1.3 Orbitais atômicos

Essas são funções matemáticas que descrevem o comportamento ondulatório de um elétron ou de um par de elétrons em um átomo. Essas funções podem ser usadas para determinar a probabilidade de encontrar um elétron em qualquer região específica ao redor do núcleo do átomo.

<br>

<p align='center'>
    <img src='https://github.com/ssebastianmag/hydrogen-wavefunctions/blob/edda6d746cbe2163f2e92e1191126d0fe7d6488a/img/(3%2C2%2C1)%5Blt%5D.png' width=50% />
</p>
<p align='center'>
    <i>Gráfico de densidade eletrônica mostrando regiões de probabilidade eletrônica variável</i>
</p>

---

#### 1.4 Números quânticos

Podemos descrever números quânticos como um conjunto de valores numéricos que fornecem uma descrição completa do estado de uma partícula quântica. Para elétrons em um átomo, há tipicamente quatro números quânticos:
<br>

- Número quântico principal ($n$): `( 1 <= n )`<br>
Representa o nível de energia do elétron e o tamanho relativo do orbital.


- Número quântico azimutal ($l$): `( 0 <= l <= n-1 )`<br>
Relaciona-se com a forma do orbital atômico.


- Número quântico magnético ($m_l$): `( -l <= m <= l )`<br>
Especifica a orientação do orbital no espaço.


- Número quântico de spin ($m_s$): `( +1/2 or -1/2 )`<br>
Descreve o spin intrínseco do elétron.

<br>

> [!NOTE]
> No átomo de hidrogênio, ou qualquer átomo com um único elétron (como hélio ionizado, lítio, etc.),
o spin do elétron não interage com mais nada para afetar sua distribuição espacial.
>
> Para nossa aplicação específica com o átomo de hidrogênio, focaremos nos três primeiros números quânticos.
Como o spin do elétron não influencia a forma ou distribuição da nuvem de elétrons.

---

### 2. Equação de Chrödinger para funções de onda do átomo de hidrogênio

A equação de Schrödinger serve como base da mecânica quântica,
é uma equação diferencial que determina as funções de onda de um sistema quântico.
Para o átomo de hidrogênio, usamos a seguinte representação da equação de Schrödinger independente do tempo:


$\large \hat{H} \psi = E \psi$

$H$ é o operador hamiltoniano, que representa a energia total (cinética + potencial) do sistema,
e $E$ é a energia total do sistema.

Dada a simetria esférica do átomo de hidrogênio, podemos expressá-la em termos de
coordenadas esféricas $(r, \theta, \varphi)$ em vez de coordenadas retangulares $(x, y, z)$.
Onde $r$ é a coordenada radial, $\theta$ é o ângulo polar (relativo ao eixo z vertical),
e $\varphi$ é o ângulo azimutal (relativo ao eixo x).

<p align='center'>
  <img src='https://github.com/ssebastianmag/hydrogen-wavefunctions/blob/main/img/coordinate_system.png' width=38% />
</p>
<p align='center'>
    <i>Relação entre os sistemas de coordenadas esféricas e retangulares</i>
</p>

A função de onda $\psi(r, \theta, \varphi)$ pode ser representada como um produto de funções radiais e angulares:

$\large \psi(r, \theta, \varphi) = R(r) Y(\theta, \varphi)$

Quando o hamiltoniano é expresso em coordenadas esféricas, ele contém partes radiais e angulares.
Ao substituir isso na equação de Schrödinger, separamos a equação em duas partes:
uma que depende apenas de $r$ (a parte radial) e outra que depende de $\theta$ e $\varphi$ (a parte angular).


---

#### 2.1 Componente radial

$\large R_{n \ell}(r) = \sqrt{\left( \frac{2}{n a_0} \right)^3 \frac{(n-\ell-1)!}{2n(n+\ell)!}} e^{-\frac{r}{n a_0}} \left( \frac{2r}{n a_0} \right)^{\ell} L_{n-\ell-1}^{2\ell+1}\left(\frac{2r}{n a_0}\right)$

A função de onda radial nos dá informações sobre a distribuição de probabilidade
do elétron como uma função da distância $r$ do
núcleo. Sua forma abrange três termos principais:

**2.1.1 Decaimento Exponencial**: Significa o decaimento de probabilidade de encontrar um
elétron conforme nos afastamos do núcleo. Aqui, $a_0$ é o raio de Bohr
que define uma escala característica para dimensões atômicas:

$\large e^{-\frac{r}{n a_0}}$

<br>

**2.1.2 Termo de potência**: Determina como a probabilidade muda com $r$.
O número quântico azimutal $\ell$ desempenha um papel significativo na determinação
do número de nós na distribuição radial:

$\large \left( \frac{2r}{n a_0} \right)^{\ell}$

<br>

**2.1.3 Polinômios de Laguerre associados**: Esses polinômios contribuem para a estrutura mais fina da parte radial,
especialmente definindo nós (regiões onde a probabilidade é zero):

$\large L_{n-\ell-1}^{2\ell+1}\left(\frac{2r}{n a_0}\right)$

---

#### 2.2 Componente angular

$\large Y_{\ell}^{m}(\theta, \varphi) = (-1)^m \sqrt{\frac{(2\ell+1)}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}} P_{\ell}^{m}(\cos\theta) e^{im\varphi}$

A função de onda angular produz os harmônicos esféricos, que fornecem a dependência angular da função de onda em
termos dos ângulos polar ($\theta$) e azimutal ($\varphi$).

Esses harmônicos esféricos fornecem um relato detalhado das formas e orientações dos orbitais atômicos,
caracterizando como as distribuições de probabilidade de elétrons são espalhadas no espaço.
Ele tem dois componentes:

**2.2.1 Polinômios de Legendre Associados**: Eles ditam a forma do orbital na direção polar ($\theta$),
ajudando a definir as formas características (s, p, d, etc.) que frequentemente associamos aos orbitais atômicos:

$\large P_{\ell}^{m}(\cos\theta)$

<br>

**2.2.2 Termo Azimutal Exponencial**: Este termo fornece a orientação do orbital no plano azimutal, conforme
determinado pelo número quântico magnético $m$:

$\large e^{im\varphi}$

---

#### 2.3 Função de onda normalizada

TA função de onda normalizada resultante para o átomo de hidrogênio é o produto das soluções dos componentes radial e angular:

$\large \psi_{n \ell m}(r, \theta, \varphi) = R_{n \ell}(r) Y_{\ell}^{m}(\theta, \varphi)$

<br>

Para determinar a densidade de probabilidade do elétron estar em um determinado local,
integramos a magnitude quadrada da função de onda sobre todo o espaço: $|\psi_{n \ell m}|^2$

$\large P(r, \theta, \varphi) = |\psi_{n \ell m}(r, \theta, \varphi)|^2$

---
> Por meio da análise do modelo de função de onda do átomo de hidrogênio, o comportamento e a distribuição da densidade de elétrons
dentro dos sistemas atômicos se tornam aparentes, lançando luz sobre a incerteza inerente da mecânica quântica.
---

## Implementação

#### Command line arguments:

```
$ python hydrogen_wavefunction_cli.py --help
```

```   
usage: main.py [-h] [--dark_theme] [--colormap COLORMAP] [n] [l] [m] [a0_scale_factor]

Átomo de hidrogênio - Visualização da função de onda e densidade eletrônica
para estados quânticos específicos (n, l, m).  

positional arguments:
  n                     (n) Número quântico principal (int)
  l                     (l) Número quântico azimutal (int)
  m                     (m) Número quântico magnético (int)
  a0_scale_factor       Fator de escala do raio de Bohr (float)

options:
  -h, --help            Mostra essa mensagem de ajuda e sai
  --dark_theme          Caso setado, será plotado um gráfico com tema escuro
  --colormap COLORMAP   Seaborn plot colormap

```

---

#### Input args:
    $ python main.py 3 2 1 0.3

|    Parameter    |            Description            | Value |  Constraint   |
|:---------------:|:---------------------------------:|:-----:|:-------------:|
|        n        |  Número quântico principal ($n$)   |   3   |    1 <= n     |
|        l        | Número quântico azimutal ($\ell$) |   2   | 0 <= l <= n-1 |
|        m        |   Número quântico magnético ($m$)   |   1   | -l <= m <= l  |
| a0_scale_factor | Fator de escala do raio de Bohr ($a_0$)  |  0.3  |               |
|   dark_theme    |      Habilitar tema escuro do gráfico       |       |               |
|    colormap     |       Seaborn plot colormap       |       |               |

#### Output:

<p align='left'>
  <img src='https://github.com/ssebastianmag/hydrogen-wavefunctions/blob/edda6d746cbe2163f2e92e1191126d0fe7d6488a/img/(3%2C2%2C1)%5Blt%5D.png' width=60% />
</p>

---

#### Input args:
    $ python main.py 3 2 1 0.3 --dark_theme

|    Parameter    |            Description            | Value |  Constraint   |
|:---------------:|:---------------------------------:|:-----:|:-------------:|
|        n        |  Número quântico principal ($n$)   |   3   |    1 <= n     |
|        l        | Número quântico azimutal ($\ell$) |   2   | 0 <= l <= n-1 |
|        m        |   Número quântico magnético ($m$)   |   1   | -l <= m <= l  |
| a0_scale_factor | Fator de escala do raio de Bohr ($a_0$)  |  0.3  |               |
|   dark_theme    |      Habilitar tema escuro do gráfico       |   --dark_theme   |               |
|    colormap     |       Seaborn plot colormap       |       |               |
#### Output:

<p align='left'>
  <img src='https://github.com/ssebastianmag/hydrogen-wavefunctions/blob/edda6d746cbe2163f2e92e1191126d0fe7d6488a/img/(3%2C2%2C1)%5Bdt%5D.png' width=60% />
</p>

---
#### Input args:
    $ python main.py 20 10 5 0.01 --dark_theme --colormap "mako"

|    Parameter    |            Description            |    Value     |  Constraint   |
|:---------------:|:---------------------------------:|:------------:|:-------------:|
|        n        |  Número quântico principal ($n$)   |      20      |    1 <= n     |
|        l        | Número quântico azimutal ($\ell$) |      10      | 0 <= l <= n-1 |
|        m        |   Número quântico magnético ($m$)   |      5       | -l <= m <= l  |
| a0_scale_factor | Fator de escala do raio de Bohr ($a_0$)  |     0.01     |               |
|   dark_theme    |      Habilitar tema escuro do gráfico       | --dark_theme |               |
|    colormap     |       Seaborn plot colormap       |    "mako"    |               |

#### Output:

<p align='left'>
  <img src='https://github.com/ssebastianmag/hydrogen-wavefunctions/blob/edda6d746cbe2163f2e92e1191126d0fe7d6488a/img/(20%2C10%2C5)%5Bdt%5D.png' width=60% />
</p>

Para números quânticos extremamente altos, os seguintes efeitos podem ser observados:

- A complexidade aumenta ainda mais, resultando em vários nós e padrões intrincados.
- Avaliar a função de onda em um vasto domínio espacial se torna computacionalmente intensivo.
- A visualização pode se tornar confusa, dificultando o discernimento de detalhes ou características específicas.

---
