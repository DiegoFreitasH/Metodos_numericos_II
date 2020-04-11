import math
from algorithms import *

def function(x):
    '''Function f(x) to be integrated'''
    return ((math.sin(2*x) + 4*(x**2) + 3*x)**2)

integral_precisa = 17.8764703
'''Integrais da função de teste seguindo as especificações das integrais de
Gauss-Hermite -> ∫ e^(-(x^2))f(x)dx, intervalo (-∞, +∞)
Gauss-Laguerre -> ∫ e^(-x)f(x)dx, intervalo (0, +∞)
Gauss-Chebyshev -> ∫ (1/sqrt(1-x^2))f(x)dx, intervalo (-1, +1) '''
integral_precisa_gauss_hermite = 34.0277796460935
integral_precisa_gauss_laguerre = 547.1745882352782
integral_precisa_gauss_chebyshev = 46.05236716716156

def output_newton_cotes(f, output_to_file=True, path="newton-cotes"):
    trapezio = integrate(f, 0, 1, math.pow(10, -6), regra_trapezio, output_to_file, path)
    trapezio_aberto = integrate(f, 0, 1, math.pow(10, -6), regra_trapezio_aberto, output_to_file, path)
    simpson = integrate(f, 0, 1, math.pow(10, -6), regra_simpson, output_to_file, path)
    milne = integrate(f, 0, 1, math.pow(10, -6), regra_milne, output_to_file, path)
    simpson_3_8 = integrate(f, 0, 1, math.pow(10, -6), regra_simpson_3_8, output_to_file, path)
    poli_3_aberta = integrate(f, 0, 1, math.pow(10, -6), regra_poli_3_aberta, output_to_file, path)
    poli_4_fechada = integrate(f, 0, 1, math.pow(10, -6), regra_poli_4_fechada, output_to_file, path)
    poli_4_aberta = integrate(f, 0, 1, math.pow(10, -6), regra_poli_4_aberta, output_to_file, path)

    print(f"Trapezio Fechado \t= {trapezio:.7f}")
    print(f"Trapezio Aberto \t= {trapezio_aberto:.7f}")
    print(f"Simpson 1/8 \t\t= {simpson:.7f}")
    print(f"Milne \t\t\t= {milne:.7f}")
    print(f"Simpson 3/8 \t\t= {simpson_3_8:.7f}")
    print(f"Polinomio 3 Aberta \t= {poli_3_aberta:.7f}")
    print(f"Polinomio 4 Fechada \t= {poli_4_fechada:.7f}")
    print(f"Polinomio 4 Aberta \t= {poli_4_aberta:.7f}")
    print(f"Integral \t\t= {integral_precisa}\n")


def output_gauss_legendre(f, output_to_file=True, path="gauss-legendre"):
    gl_2p = integrate(f, 0, 1, math.pow(10, -6), gauss_legendre_2p, output_to_file, path)
    gl_3p = integrate(f, 0, 1, math.pow(10, -6), gauss_legendre_3p, output_to_file, path)
    gl_4p = integrate(f, 0, 1, math.pow(10, -6), gauss_legendre_4p, output_to_file, path)

    print(f"Gauss-Legendre 2 pontos = {gl_2p:.7f}")
    print(f"Gauss-Legendre 3 pontos = {gl_3p:.7f}")
    print(f"Gauss-Legendre 4 pontos = {gl_4p:.7f}")
    print(f"Integral \t\t= {integral_precisa}\n")

def output_gauss_hermite(f):
    gh_2 = gauss_hermite_n2(f)
    gh_3 = gauss_hermite_n3(f)
    gh_4 = gauss_hermite_n4(f)
    print(f"Gauss Hermite 2 \t= {gh_2:.7f}")
    print(f"Gauss Hermite 3 \t= {gh_3:.7f}")
    print(f"Gauss Hermite 4 \t= {gh_4:.7f}")
    print(f"Gauss Hermite Precisa \t= {integral_precisa_gauss_hermite:.7f}\n")

def output_gauss_laguerre(f):
    gl_2 = gauss_laguerre_n2(f)
    gl_3 = gauss_laguerre_n3(f)
    print(f"Gauss Laguerre 2 \t= {gl_2:.7f}")
    print(f"Gauss Laguerre 3 \t= {gl_3:.7f}")
    print(f"Gauss Laguerre Precisa \t= {integral_precisa_gauss_laguerre:.7f}\n")

def output_gauss_chebyshev(f):
    gc_2 = gauss_chebyshev_n2(f)
    gc_3 = gauss_chebyshev_n3(f)
    gc_4 = gauss_chebyshev_n4(f)
    print(f"Gauss Chebyshev 2 \t= {gc_2:.7f}")
    print(f"Gauss Chebyshev 3 \t= {gc_3:.7f}")
    print(f"Gauss Chebyshev 4 \t= {gc_4:.7f}")
    print(f"Gauss Chebyshev Precisa = {integral_precisa_gauss_chebyshev:.7f}\n")

def main():
    output_newton_cotes(function)
    output_gauss_legendre(function)
    output_gauss_hermite(function)
    output_gauss_laguerre(function)
    output_gauss_chebyshev(function)

if __name__ == "__main__":
    main()