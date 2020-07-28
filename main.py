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

def output_potencia_regular():
    A1 = [
        [5, 2, 1],
        [2, 3, 1],
        [1, 1, 2]
    ]
    A2 = [
        [40, 8, 4, 2, 1],
        [8, 30, 12, 6, 2],
        [4, 12, 20, 1, 2],
        [2, 6, 1, 25, 4],
        [1, 2, 2, 4, 5]
    ]
    v0 = [1, 1, 1]
    autovalor_A1, autovetor_A1 = potencia_regular(A1, v0, 10**-7)
    autovalor_A2, autovetor_A2 = potencia_regular(A2, [1,1,1,1,1], 10**-7)
    
    print(f"\n{'-'*6}Potencia Regular{'-'*6}")
    print(f'Autovalor Dominante de A1 = {autovalor_A1}')
    print(f'Autovetor Dominante de A1 = {autovetor_A1}\n')
    print(f'Autovalor Dominante de A2 = {autovalor_A2}')
    print(f'Autovetor Dominante de A2 = {autovetor_A2}')
    

def output_potencia_inversa():
    A1 = [
        [5, 2, 1],
        [2, 3, 1],
        [1, 1, 2]
    ]
    A2 = [
        [-14, 1, -2],
        [1, -1, 1],
        [-2, 1, -11]
    ]
    A3 = [
        [40, 8, 4, 2, 1],
        [8, 30, 12, 6, 2],
        [4, 12, 20, 1, 2],
        [2, 6, 1, 25, 4],
        [1, 2, 2, 4, 5]
    ]
    v0 = [1, 1, 1]
    autovalor_A1, autovetor_A1 = potencia_inversa(A1, v0, 10**-7)
    autovalor_A2, autovetor_A2 = potencia_inversa(A2, v0, 10**-7)
    autovalor_A3, autovetor_A3 = potencia_inversa(A3, [1,1,1,1,1], 10**-7)
    
    print(f"\n{'-'*6}Potencia Inversa{'-'*6}")
    print(f'Menor Autovalor de A1 = {autovalor_A1}')
    print(f'Menor Autovetor de A1 = {autovetor_A1}\n')
    print(f'Menor Autovalor de A2 = {autovalor_A2}')
    print(f'Menor Autovetor de A2 = {autovetor_A2}\n')
    print(f'Menor Autovalor de A3 = {autovalor_A3}')
    print(f'Menor Autovetor de A3 = {autovetor_A3}\n')

def output_potencia_com_deslocamento(u1, u2, u3):
    A1 = [
        [5, 2, 1],
        [2, 3, 1],
        [1, 1, 2]
    ]
    A2 = [
        [-14, 1, -2],
        [1, -1, 1],
        [-2, 1, -11]
    ]
    A3 = [
        [40, 8, 4, 2, 1],
        [8, 30, 12, 6, 2],
        [4, 12, 20, 1, 2],
        [2, 6, 1, 25, 4],
        [1, 2, 2, 4, 5]
    ]
    v0 = [1, 1, 1]
    autovalor_A1, autovetor_A1 = potencia_com_deslocamento(A1, [1, 1, 1], 10**-7, u1)
    autovalor_A2, autovetor_A2 = potencia_com_deslocamento(A2, [1, 1, 1], 10**-7, u2)
    autovalor_A3, autovetor_A3 = potencia_com_deslocamento(A3, [1,1,1,1,1], 10**-7, u3)
    
    print(f"\n{'-'*6}Potencia Com Deslocamento{'-'*6}")
    print(f'Autovalor de A1 = {autovalor_A1}')
    print(f'Autovetor de A1 = {autovetor_A1}')
    print(f'Deslocamento de A1 = {u1}\n')
    print(f'Autovalor de A2 = {autovalor_A2}')
    print(f'Autovetor de A2 = {autovetor_A2}')
    print(f'Deslocamento de A2 = {u2}\n')
    print(f'Autovalor de A3 = {autovalor_A3}')
    print(f'Autovetor de A3 = {autovetor_A3}')
    print(f'Deslocamento de A3 = {u3}\n')

def output_metodo_de_householder():
    A = [
        [40, 8, 4, 2, 1],
        [8, 30, 12, 6, 2],
        [4, 12, 20, 1, 2],
        [2, 6, 1, 25, 4],
        [1, 2, 2, 4, 5]
    ]

    A_barra, H = metodo_de_householder(A, len(A))
    
    max_autovalor, max_autovetor = potencia_regular(A_barra, [1, 1, 1, 1, 1], 10**-7)
    min_autovalor, min_autovetor = potencia_inversa(A_barra, [1, 1, 1, 1, 1], 10**-7)

    autovalores_A_barra = [min_autovalor]
    autovetores_A_barra = [min_autovetor]

    for i in range(math.floor(min_autovalor), math.floor(max_autovalor), 5):
        autovalor, autovetor = potencia_com_deslocamento(A_barra, [1, 1, 1, 1, 1], 10**-7, i)
        if abs(autovalor - autovalores_A_barra[-1]) > 1:
            autovalores_A_barra.append(autovalor)
            autovetores_A_barra.append(autovetor)
    
    autovetores_A = []

    for autovetor in autovetores_A_barra:
        autovetores_A.append(tranpose_matrix(mult_matrix_matrix(H, tranpose_matrix([autovetor])))[0])

    resultados_A = [(n,m) for n,m in zip(autovalores_A_barra, autovetores_A)]
    
    print(f"\n{'-'*6}Metodo de Householder{'-'*6}")
    print('Matriz A Complementar = {')
    print_matrix(A_barra)
    print('}')
    print('Matriz H Acumulada = {')
    print_matrix(H)
    print('}')
    print('\n-> Autovalores e Autovalores\n')
    for resultado in resultados_A:
        print(f"Autovalor: {resultado[0]}")
        print(f"Autovetor: {resultado[1]}")
        print()
        
def output_metodo_de_jacobi():
    A = [
        [40, 8, 4, 2, 1],
        [8, 30, 12, 6, 2],
        [4, 12, 20, 1, 2],
        [2, 6, 1, 25, 4],
        [1, 2, 2, 4, 5]
    ]
    n = len(A)

    print(f"\n{'-'*6}Metodo de Jacobi{'-'*6}")
    P, autovalores = metodoDeJacobi(A, n, 10**-7, True)
    A_barra = [[autovalores[i] if i == j else 0 for i in range(n)] for j in range(n)]

    print("Matriz A barra = {")
    print_matrix(A_barra)
    print("}")

    print("Matriz Acumulada P = {")
    print_matrix(P)
    print("}")

    paresVetoresValores = [i for i in zip(autovalores, tranpose_matrix(P))]
    for i in range(n):
        print(f"Autovalor: {paresVetoresValores[i][0]}")
        print(f"Autovetor: {paresVetoresValores[i][1]}\n")

    print("Varredura de Jacobi com Matriz Tridiagonal")

    # TODO Item 2.1
    A_barra, H = metodo_de_householder(A, n)
    P, autovalores = metodoDeJacobi(A_barra, n, 10**-7)
    
    print("Matriz Acumulada P = {")
    print_matrix(P)
    print("}")
    
    P = mult_matrix_matrix(H, P)

    print("Matriz Acumulada P*H = {")
    print_matrix(P)
    print("}")
    
def output_metodo_QR():
    A = [
        [40, 8, 4, 2, 1],
        [8, 30, 12, 6, 2],
        [4, 12, 20, 1, 2],
        [2, 6, 1, 25, 4],
        [1, 2, 2, 4, 5]
    ]
    n = len(A)

    print(f"\n{'-'*6} Metodo QR {'-'*6}")
    P, autovalores = metodoQR(A, n, 10**-7)

    print("Matriz Acumulada P = {")
    print_matrix(P)
    print("}")
    
    # TODO Item 1.3
    paresVetoresValores = [i for i in zip(autovalores, tranpose_matrix(P))]
    for i in range(n):
        print(f"Autovalor: {paresVetoresValores[i][0]}")
        print(f"Autovetor: {paresVetoresValores[i][1]}\n")

    # TODO Item 2.1 
    A_barra, H = metodo_de_householder(A, n)
    P, autovalores = metodoQR(A_barra, n, 10**-7)
    
    print("Metodo QR com matriz Tridiagonal")

    print("Matriz Acumulada P = {")
    print_matrix(P)
    print("}")
    
    P = mult_matrix_matrix(H, P)

    print("Matriz Acumulada P*H = {")
    print_matrix(P)
    print("}")


def main():
    # output_newton_cotes(function)
    # output_gauss_legendre(function)
    # output_gauss_hermite(function)
    # output_gauss_laguerre(function)
    # output_gauss_chebyshev(function)
    # output_singularidades(function_1, function_2)
    # output_surface_area()
    # output_volume()
    # output_potencia_regular()
    # output_potencia_inversa()
    # output_potencia_com_deslocamento(2, -6, 10)
    # output_metodo_de_householder()
    output_metodo_de_jacobi()
    output_metodo_QR()

if __name__ == "__main__":
    main()
