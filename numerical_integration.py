import math
from algorithms import *

def function(x):
    '''Function f(x) to be integrated'''
    return ((math.sin(2*x) + 4*(x**2) + 3*x)**2)

# Funções para o problema de singularidade
def function_1(x):
    '''Função 1 da tarefa da aula 14'''
    return 1/((x*x)**(1/3.0))

function_1_precisa = 6.0

def function_2(x):
    '''Função 2 da tarefa da aula 14'''
    return 1/(math.sqrt(4-x**2))

function_2_precisa = math.pi / 2.0

integral_precisa = 17.8764703
'''Integrais da função de teste seguindo as especificações das integrais de
Gauss-Hermite -> ∫ e^(-(x^2))f(x)dx, intervalo (-∞, +∞)
Gauss-Laguerre -> ∫ e^(-x)f(x)dx, intervalo (0, +∞)
Gauss-Chebyshev -> ∫ (1/sqrt(1-x^2))f(x)dx, intervalo (-1, +1) '''
integral_precisa_gauss_hermite = 34.0277796460935
integral_precisa_gauss_laguerre = 547.1745882352782
integral_precisa_gauss_chebyshev = 46.05236716716156

def output_newton_cotes(function, output_to_file=True, path="newton-cotes"):
    trapezio = integrate(function, 0, 1, math.pow(10, -6), regra_trapezio, output_to_file, path)
    trapezio_aberto = integrate(function, 0, 1, math.pow(10, -6), regra_trapezio_aberto, output_to_file, path)
    simpson = integrate(function, 0, 1, math.pow(10, -6), regra_simpson, output_to_file, path)
    milne = integrate(function, 0, 1, math.pow(10, -6), regra_milne, output_to_file, path)
    simpson_3_8 = integrate(function, 0, 1, math.pow(10, -6), regra_simpson_3_8, output_to_file, path)
    poli_3_aberta = integrate(function, 0, 1, math.pow(10, -6), regra_poli_3_aberta, output_to_file, path)
    poli_4_fechada = integrate(function, 0, 1, math.pow(10, -6), regra_poli_4_fechada, output_to_file, path)
    poli_4_aberta = integrate(function, 0, 1, math.pow(10, -6), regra_poli_4_aberta, output_to_file, path)

    print(f"{'-'*6}Formulas de Newton-Cotes{'-'*6}")
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

    print(f"{'-'*6}Quadratura de Gauss-Legendre{'-'*6}")
    print(f"Gauss-Legendre 2 pontos = {gl_2p:.7f}")
    print(f"Gauss-Legendre 3 pontos = {gl_3p:.7f}")
    print(f"Gauss-Legendre 4 pontos = {gl_4p:.7f}")
    print(f"Integral \t\t= {integral_precisa}\n")

def output_gauss_hermite(f):
    gh_2 = gauss_hermite_n2(f)
    gh_3 = gauss_hermite_n3(f)
    gh_4 = gauss_hermite_n4(f)

    print(f"{'-'*6}Quadratura de Gauss-Hermite{'-'*6}")
    print(f"Gauss Hermite 2 \t= {gh_2:.7f}")
    print(f"Gauss Hermite 3 \t= {gh_3:.7f}")
    print(f"Gauss Hermite 4 \t= {gh_4:.7f}")
    print(f"Gauss Hermite Precisa \t= {integral_precisa_gauss_hermite:.7f}\n")

def output_gauss_laguerre(f):
    gl_2 = gauss_laguerre_n2(f)
    gl_3 = gauss_laguerre_n3(f)
    gl_4 = gauss_laguerre_n4(f)

    print(f"{'-'*6}Quadratura de Gauss-Leguerre{'-'*6}")
    print(f"Gauss Laguerre 2 \t= {gl_2:.7f}")
    print(f"Gauss Laguerre 3 \t= {gl_3:.7f}")
    print(f"Gauss Laguerre 4 \t= {gl_4:.7f}")
    print(f"Gauss Laguerre Precisa \t= {integral_precisa_gauss_laguerre:.7f}\n")

def output_gauss_chebyshev(f):
    gc_2 = gauss_chebyshev_n2(f)
    gc_3 = gauss_chebyshev_n3(f)
    gc_4 = gauss_chebyshev_n4(f)

    print(f"{'-'*6}Quadratura de Gauss-Chebyshev{'-'*6}")
    print(f"Gauss Chebyshev 2 \t= {gc_2:.7f}")
    print(f"Gauss Chebyshev 3 \t= {gc_3:.7f}")
    print(f"Gauss Chebyshev 4 \t= {gc_4:.7f}")
    print(f"Gauss Chebyshev Precisa = {integral_precisa_gauss_chebyshev:.7f}\n")

def output_singularidades(f_1, f_2):
    expo_simples_1 = 2 * exponencial_simples(f_1, 0, 1, pow(10, -7))
    expo_dupla_1 = 2 * exponencial_dupla(f_1, 0, 1, pow(10, -7))
    expo_simples_2 = exponencial_simples(f_2, -2, 0, pow(10, -7))
    expo_dupla_2 = exponencial_dupla(f_2, -2, 0, pow(10, -7))

    print(f"{'-'*6}Integrais com Singularidade{'-'*6}")
    print(f"Exponencial Simples da Função 1 = {expo_simples_1:.7f}")
    print(f"Exponencial Dupla da Função 1 \t= {expo_dupla_1:.7f}")
    print(f"Integral Precisa 1 \t\t= {function_1_precisa:.7f}")
    print(f"Exponencial Simples da Função 2 = {expo_simples_2:.7f}")
    print(f"Exponencial Dupla da Função 2 \t= {expo_dupla_2:.7f}")
    print(f"Integral Precisa 2 \t\t= {function_2_precisa:.7f}")

def output_surface_area():
    w_list = [5/9, 8/9, 5/9]
    x_list = [-math.sqrt(3/5), 0, math.sqrt(3/5)]

    def g(r, s):
        alpha = (1 + r) / 2.0
        beta = (1 + s) * math.pi
        x = 40 * alpha * math.cos(beta)
        y = 40 * alpha * math.sin(beta)
        return math.sqrt((0.4*x)**2 + (0.4*y)**2 + 1) * 1600 *  alpha * math.pi/2
    
    I = 0
    for i, b in enumerate(x_list):
        for j, a in enumerate(x_list):
            I += w_list[j]*w_list[i]*g(a,b)

    print(f"\n{'-'*6}Calculo de área da superfice{'-'*6}")
    print(f"Area da Surperfice problema 2 \t= {I:.7f}")

def output_volume():
    w_list = [5.0/9, 8.0/9, 5.0/9]
    x_list = [-math.sqrt(3/5), 0, math.sqrt(3/5)]

    def g(r, s):
        alpha = (r + 1) / 2.0
        beta = (s + 1) * math.pi/4
        x = 40 * alpha * math.cos(beta) 
        y = 20 * alpha * math.sin(beta)
        return 0.2 * (x**2 - y**2) * 800 * alpha * math.pi/8
   
    I = 0
    for i, b in enumerate(x_list):
        for j, a in enumerate(x_list):
            value = w_list[j]*w_list[i]*g(a,b)
            I += value
    
    I *= 4
    print(f"\n{'-'*6}Calculo do Volume{'-'*6}")
    print(f"Volume do Solido problema 2 \t= {I:.7f}")

def main():
    # output_newton_cotes(function)
    # output_gauss_legendre(function)
    # output_gauss_hermite(function)
    # output_gauss_laguerre(function)
    # output_gauss_chebyshev(function)
    # output_singularidades(function_1, function_2)
    # output_surface_area()
    # output_volume()
    A = [
        [5, 2, 1],
        [2, 3, 1],
        [1, 1, 2]
    ]
    v0 = [1, 1, 1]
    # decomLU(A)
    potencia_inversa(A, v0, 10**-7)
    # potencia_regular(A, v0, 10**-7)

if __name__ == "__main__":
    main()
