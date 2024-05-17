from scipy.constants import physical_constants
import matplotlib.pyplot as plt
import scipy.special as sp
import seaborn as sns
import numpy as np
import argparse


def radial_function(n, l, r, a0):
    """ Cálculo da parte radial da função de onda utilizando polinômios de 
    Laguerre e um fator de decaimento exponencial.

    Args:
        n (int): Número quântico principal
        l (int): Númeria quântico azimutal
        r (numpy.ndarray): Coordenada radial
        a0 (float): Raio de Bohr em escala
    Returns:
        numpy.ndarray: Componente radial da função de onda
    """

    laguerre = sp.genlaguerre(n - l - 1, 2 * l + 1)
    p = 2 * r / (n * a0)

    constant_factor = np.sqrt(
        ((2 / n * a0) ** 3 * (sp.factorial(n - l - 1))) /
        (2 * n * (sp.factorial(n + l)))
    )
    return constant_factor * np.exp(-p / 2) * (p ** l) * laguerre(p)


def angular_function(m, l, theta, phi):
    """ Cálculo da parte radial da função de onda utilizando polinômios de 
    Laguerre e um fator de decaimento exponencial.

    Args:
        n (int): Número quântico principal
        l (int): Númeria quântico azimutal
        r (numpy.ndarray): Coordenada radial
        a0 (float): Raio de Bohr em escala
    Returns:
        numpy.ndarray: Componente radial da função de onda
    """

    legendre = sp.lpmv(m, l, np.cos(theta))

    constant_factor = ((-1) ** m) * np.sqrt(
        ((2 * l + 1) * sp.factorial(l - np.abs(m))) /
        (4 * np.pi * sp.factorial(l + np.abs(m)))
    )
    return constant_factor * legendre * np.real(np.exp(1.j * m * phi))


def compute_wavefunction(n, l, m, a0_scale_factor):
    """ Cálculo da parte radial da função de onda utilizando polinômios de 
    Laguerre e um fator de decaimento exponencial.

    Args:
        n (int): Número quântico principal
        l (int): Númeria quântico azimutal
        r (numpy.ndarray): Coordenada radial
        a0 (float): Raio de Bohr em escala
    Returns:
        numpy.ndarray: Componente radial da função de onda
    """

    # Raio de Bohr em escala para melhor visualização 
    a0 = a0_scale_factor * physical_constants['Bohr radius'][0] * 1e+12

    # z-x plano para representar a distribuição espacial de elétrons
    grid_extent = 480
    grid_resolution = 680
    z = x = np.linspace(-grid_extent, grid_extent, grid_resolution)
    z, x = np.meshgrid(z, x)

    # Use epsilon para evitar divisão por zero durante cáculo de ângulos 
    eps = np.finfo(float).eps

    # Ψnlm(r,θ,φ) = Rnl(r).Ylm(θ,φ)
    psi = radial_function(
        n, l, np.sqrt((x ** 2 + z ** 2)), a0
    ) * angular_function(
        m, l, np.arctan(x / (z + eps)), 0
    )
    return psi


def compute_probability_density(psi):
    """ Calculo de densidade de probabilidade de uma determinada função de onda.
    Args:
        psi (numpy.ndarray): Função de onda
    Returns:
        numpy.ndarray: Densidade de probabilidade de uma função de onda
    """
    return np.abs(psi) ** 2


def plot_wf_probability_density(n, l, m, a0_scale_factor, dark_theme=False, colormap='rocket'):
    """Plotar a densidade de probabilidade da função de onda do átomo de 
    hidrogênio para determinado estado quântico (n, l, m)


    Args:
        n (int): Número quântico principal, determina o nível de energia e o tamanho do orbital
        l (int): Número quântico azimutal, determina o formato do orbital
        m (int): Número quântico magnético, define a orientação do orbital
        a0_scale_factor (float): Raio de Bohr em escala
        dark_theme (bool): Se VERDADEIRO, usa um fundo preto para a plotagem, o padrão é FALSO
        colormap (str): Mapa de cor Seaborn, padrão é 'rocket'
    """

    # Validação dos números quânticos
    if not isinstance(n, int) or n < 1:
        raise ValueError('n deve ser um número inteiro que satisfaça a condição: n >= 1')
    if not isinstance(l, int) or not (0 <= l < n):
        raise ValueError('l deve ser um número inteiro que satisfaça a condição: 0 <= l < n')
    if not isinstance(m, int) or not (-l <= m <= l):
        raise ValueError('m deve ser um número inteiro que satisfaça a condição: -l <= m <= l')

    # Validação do mapa de cor
    try:
        sns.color_palette(colormap)
    except ValueError:
        raise ValueError(f'{colormap} não é um mapa de cor reconhecido.')

    # Configura a estética da plotagem usando matplotlib rcParams
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['xtick.major.width'] = 4
    plt.rcParams['ytick.major.width'] = 4
    plt.rcParams['xtick.major.size'] = 15
    plt.rcParams['ytick.major.size'] = 15
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    plt.rcParams['axes.linewidth'] = 4

    fig, ax = plt.subplots(figsize=(16, 16.5))
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(right=0.905)
    plt.subplots_adjust(left=-0.1)

    # Computa e visualiza a densidade de probabilidade da função de onda
    psi = compute_wavefunction(n, l, m, a0_scale_factor)
    prob_density = compute_probability_density(psi)

    # Aqui transpomos a matriz para alinhar o plano z-x calculado com a exibição y-x imshow do Matplotlib
    im = ax.imshow(np.sqrt(prob_density).T, cmap=sns.color_palette(colormap, as_cmap=True))

    cbar = plt.colorbar(im, fraction=0.046, pad=0.03)
    cbar.set_ticks([])

    # Aplicar parâmetros de fundo escuro
    if dark_theme:
        theme = 'dt'
        background_color = sorted(
            sns.color_palette(colormap, n_colors=100),
            key=lambda color: 0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2]
        )[0]
        plt.rcParams['text.color'] = '#dfdfdf'
        title_color = '#dfdfdf'
        fig.patch.set_facecolor(background_color)
        cbar.outline.set_visible(False)
        ax.tick_params(axis='x', colors='#c4c4c4')
        ax.tick_params(axis='y', colors='#c4c4c4')
        for spine in ax.spines.values():
            spine.set_color('#c4c4c4')

    else:  # Aplicar parâmetros de fundo claro
        theme = 'lt'
        plt.rcParams['text.color'] = '#000000'
        title_color = '#000000'
        ax.tick_params(axis='x', colors='#000000')
        ax.tick_params(axis='y', colors='#000000')

    ax.set_title('Hydrogen Atom - Wavefunction Electron Density', pad=130, fontsize=44, loc='left', color=title_color)
    ax.text(0, 722, (
        r'$|\psi_{n \ell m}(r, \theta, \varphi)|^{2} ='
        r' |R_{n\ell}(r) Y_{\ell}^{m}(\theta, \varphi)|^2$'
    ), fontsize=36)
    ax.text(30, 615, r'$({0}, {1}, {2})$'.format(n, l, m), color='#dfdfdf', fontsize=42)
    ax.text(770, 140, 'Probabilidade de distribuição do elétron', rotation='vertical', fontsize=30)
    ax.text(705, 700, 'Maior\nprobabilidade', fontsize=24)
    ax.text(705, -60, 'Menor\nprobabilidade', fontsize=24)
    ax.text(775, 590, '+', fontsize=34)
    ax.text(769, 82, '−', fontsize=34, rotation='vertical')
    ax.invert_yaxis()

    # Salvar e mostrar a plotagem
    plt.savefig(f'({n},{l},{m})[{theme}].png')
    # plt.show()


# - - - Execução:
if __name__ == '__main__':

    # Configurando as linhas de comando dos argumentos
    parser = argparse.ArgumentParser(
        description='Hydrogen Atom - Wavefunction and Electron Density Visualization '
                    'for specific quantum states (n, l, m).'
    )

    # Argumentos requeridos
    parser.add_argument('n', type=int, help='(n) Número quântico principal (int)', nargs='?')
    parser.add_argument('l', type=int, help='(l) Número quântico azimutal (int)', nargs='?')
    parser.add_argument('m', type=int, help='(m) Número quântico magnético (int)', nargs='?')
    parser.add_argument('a0_scale_factor', type=float, help='Raio de Bohr (float)', nargs='?')

    # Argumentos opcionais 
    parser.add_argument('--dark_theme', action='store_true', help='If set, the plot uses a dark theme')
    parser.add_argument('--colormap', type=str, default='rocket', help='Seaborn plot colormap')
    args = parser.parse_args()

    #Essa seção vai executar se você clicar no script ou rodar sem qualquer argumento
    if args.n is None:
        print('--- --- --- --- --- --- --- --- ---')
        print('Átomo de Hidrogênio - Visualizador de Função de onda e Densidade do elétron')
        print('\nParâmetros requeridos /')

        while True:
            try:
                args.n = int(input('Número quântico principal (n): '))
                if args.n < 1:
                    raise ValueError
                break
            except ValueError:
                print('(!) n deve ser um número inteiro que satisfaça a condição: n >= 1')

        while True:
            try:
                args.l = int(input('Número quântico azimmutal (l): '))
                if not (0 <= args.l < args.n):
                    raise ValueError
                break
            except ValueError:
                print('(!) l deve ser um número inteiro que satisfaça a condição: 0 <= l < n')

        while True:
            try:
                args.m = int(input('Número quântico magnético (m): '))
                if not (-args.l <= args.m <= args.l):
                    raise ValueError
                break
            except ValueError:
                print('(!) m deve ser um número inteiro que satisfaça a condição: -l <= m <= l')

        while True:
            try:
                args.a0_scale_factor = float(input('Raio de Bohr: '))
                if args.a0_scale_factor <= 0:
                    raise ValueError
                break
            except ValueError:
                print('(!) Insira um número válido maior que  0')

        print('\nParâmetros opcionais  /')
        dark_theme_choice = input(
            '[Pressione ENTER para pular] Você quer usar o tema escuro? [yes/y] (padrão é no): '
        ).lower()
        args.dark_theme = True if dark_theme_choice in ['yes', 'y'] else False

        while True:
            args.colormap = input('[Pressione ENTER para pular] Mapa de cor Seaborn (padrão é "rocket"): ')
            if not args.colormap:
                args.colormap = 'rocket'
                break
            try:
                sns.color_palette(args.colormap)
                break
            except ValueError:
                print(f'{args.colormap} não é um mapa de cor reconhecido.')
        print('\n--- --- --- --- --- --- --- --- ---')

    # Plotar densidade de probabilidade da função de onda
    plot_wf_probability_density(args.n, args.l, args.m, args.a0_scale_factor, args.dark_theme, args.colormap)
    