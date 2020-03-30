import math


def function(x):
    return ((math.sin(2*x) + 4*(x**2) + 3*x)**2)

def regra_trapezio(f, xi, xf, k, delta_x):
    return (delta_x/2)*(f(xi) + f(xf)) 

def regra_trapezio_aberto(f, xi, xf, k, delta_x):
    h = delta_x / 3
    return (delta_x/2) * (f(xi + h) + f(xi + 2*h))

def regra_simpson(f, xi, xf, k, delta_x):
    h = delta_x / 2
    return (h/3)*(f(xi) + 4*f(xi+h) + f(xf))

def regra_milne(f, xi, xf, k, delta_x): 
    h = delta_x / 4
    return (delta_x/3) * (2*f(xi + h) - f(xi + 2*h) + 2*f(xi + 3*h))

def regra_simpson_3_8(f, xi, xf, k, delta_x):
    h = delta_x / 3
    return (3*h/8) * (f(xi) + 3*f(xi + h) + 3*f(xi + 2*h) + f(xf))

def regra_poli_3_aberta(f, xi, xf, k, delta_x):
    h = delta_x / 5
    return (5*h/24) * (11*f(xi + h) + f(xi + 2*h) + f(xi + 3*h) + 11*f(xi + 4*h))

def regra_poli_4_fechada(f, xi, xf, k, delta_x):
    h = delta_x / 4
    return (2*h/45) * (7*f(xi) + 32*f(xi + h) + 12*f(xi + 2*h) + 32*f(xi + 3*h) + 7*f(xf))

def regra_poli_4_aberta(f, xi, xf, k, delta_x):
    h = delta_x / 6
    return (3*h/10) * (11*f(xi + h) - 14*f(xi + 2*h) + 26*f(xi + 3*h) - 14*f(xi + 4*h) + 11*f(xi + 5*h))

def integrate_n(f, a, b, n, subs_func = regra_trapezio):
    I = 0
    h = (b-a)/n

    for k in range(n):
        x0 = a + k*h
        x1 = x0 + h
        I += subs_func(f, x0, x1, k, h)
    
    return I

def integrate(f, a, b, E, subs_func=regra_trapezio, debug=False):
    I = 0
    N = 1
    erro = 1
    i = 0
    
    if(debug): file = open(f"{subs_func.__name__}.txt", "w")

    while erro > E:
        N *= 2
        In = integrate_n(f, a, b, N, subs_func)
        erro = abs((In - I)/In)
        if(debug):
            i += 1
            file.write(f"Iteração {i}\n")
            file.write(f"I{i-1} \t= {I:.7f}\n")
            file.write(f"I{i} \t= {In:.7f}\n")
            file.write(f"Erro \t= {erro}\n")
            file.write("----------------\n")
        I = In
    if(debug):
        file.write(f"Regra: {subs_func.__name__}\n")
        file.write(f"Num. Iterações: {i}\n")
        file.write(f"Resultado: {I:.7f}\n")
        file.close()
    return I

trapezio = integrate(function, 0, 1, math.pow(10, -6), regra_trapezio)
trapezio_aberto = integrate(function, 0, 1, math.pow(10, -6), regra_trapezio_aberto)
simpson = integrate(function, 0, 1, math.pow(10, -6), regra_simpson)
milne = integrate(function, 0, 1, math.pow(10, -6), regra_milne)
simpson_3_8 = integrate(function, 0, 1, math.pow(10, -6), regra_simpson_3_8)
poli_3_aberta = integrate(function, 0, 1, math.pow(10, -6), regra_poli_3_aberta)
poli_4_fechada = integrate(function, 0, 1, math.pow(10, -6), regra_poli_4_fechada)
poli_4_aberta = integrate(function, 0, 1, math.pow(10, -6), regra_poli_4_aberta)

# print(f"Trapezio Fechado = {trapezio:.7f}")
# print(f"Trapezio Aberto = {trapezio_aberto:.7f}")
# print(f"Simpson 1/8 = {simpson:.7f}")
# print(f"Milne = {milne:.7f}")
# print(f"Simpson 3/8 = {simpson_3_8:.7f}")
# print(f"Polinomio 3 Aberta = {poli_3_aberta:.7f}")
# print(f"Polinomio 4 Fechada = {poli_4_fechada:.7f}")
# print(f"Polinomio 4 Aberta = {poli_4_aberta:.7f}")
# print(f"Integral = 17.8764703")
