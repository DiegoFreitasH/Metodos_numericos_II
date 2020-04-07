import math
from algorithms import *

def function(x):
    '''Function f(x) to be integrated'''
    return ((math.sin(2*x) + 4*(x**2) + 3*x)**2)

integral_precisa = 17.8764703

def output_newton_cotes(output_to_file=True):
    trapezio = integrate(function, 0, 1, math.pow(10, -6), regra_trapezio, output_to_file)
    trapezio_aberto = integrate(function, 0, 1, math.pow(10, -6), regra_trapezio_aberto, output_to_file)
    simpson = integrate(function, 0, 1, math.pow(10, -6), regra_simpson, output_to_file)
    milne = integrate(function, 0, 1, math.pow(10, -6), regra_milne, output_to_file)
    simpson_3_8 = integrate(function, 0, 1, math.pow(10, -6), regra_simpson_3_8, output_to_file)
    poli_3_aberta = integrate(function, 0, 1, math.pow(10, -6), regra_poli_3_aberta, output_to_file)
    poli_4_fechada = integrate(function, 0, 1, math.pow(10, -6), regra_poli_4_fechada, output_to_file)
    poli_4_aberta = integrate(function, 0, 1, math.pow(10, -6), regra_poli_4_aberta, output_to_file)

    print(f"Trapezio Fechado \t= {trapezio:.7f}")
    print(f"Trapezio Aberto \t= {trapezio_aberto:.7f}")
    print(f"Simpson 1/8 \t\t= {simpson:.7f}")
    print(f"Milne \t\t\t= {milne:.7f}")
    print(f"Simpson 3/8 \t\t= {simpson_3_8:.7f}")
    print(f"Polinomio 3 Aberta \t= {poli_3_aberta:.7f}")
    print(f"Polinomio 4 Fechada \t= {poli_4_fechada:.7f}")
    print(f"Polinomio 4 Aberta \t= {poli_4_aberta:.7f}")
    print(f"Integral \t\t= {integral_precisa}")

def output_gauss_legendre(output_to_file=True):
    gl_2p = integrate(function, 0, 1, math.pow(10, -6), gauss_legendre_2p, output_to_file)
    gl_3p = integrate(function, 0, 1, math.pow(10, -6), gauss_legendre_3p, output_to_file)
    gl_4p = integrate(function, 0, 1, math.pow(10, -6), gauss_legendre_4p, output_to_file)

    print(f"Gauss-Legendre 2 pontos = {gl_2p:.7f}")
    print(f"Gauss-Legendre 3 pontos = {gl_3p:.7f}")
    print(f"Gauss-Legendre 4 pontos = {gl_4p:.7f}")

    print(f"Integral \t\t= {integral_precisa}")

def main():
    # output_newton_cotes()
    output_gauss_legendre()

if __name__ == "__main__":
    main()