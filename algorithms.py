import math
import os
# ======================================= #
#       Fórmulas de newton-colos          #
# ======================================= #

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

# ======================================= #
#       Quadratura de Gauss-Legendre      #
# ======================================= #

def x(xi, xf, a):
    return ((xi + xf)/2) + ((xf - xi)/2)*a

def gauss_legendre_2p(f, xi, xf, k, delta_x):
    x_a1 = x(xi, xf, -1*math.sqrt(1/3))
    x_a2 = x(xi, xf, math.sqrt(1/3)) 
    return ((xf-xi)/2)*(f(x_a1) + f(x_a2))

def gauss_legendre_3p(f, xi, xf, k, delta_x):
    x_a1 = x(xi, xf, -1*math.sqrt(3/5))
    x_a2 = x(xi, xf, 0) 
    x_a3 = x(xi, xf, math.sqrt(3/5))
    w1 = w3 = 5/9
    w2 = 8/9
    return ((xf-xi)/2)*(f(x_a1)*w1 + f(x_a2)*w2 + f(x_a3)*w3)

def gauss_legendre_4p(f, xi, xf, k, delta_x):
    x_a1 = x(xi, xf, -1*math.sqrt( (15-2*math.sqrt(30))/35 ))
    x_a2 = x(xi, xf, -1*math.sqrt( (15+2*math.sqrt(30))/35 ))
    x_a3 = x(xi, xf, math.sqrt( (15-2*math.sqrt(30))/35 ))
    x_a4 = x(xi, xf, math.sqrt( (15+2*math.sqrt(30))/35 ))
    w1 = w3 = (18 + math.sqrt(30))/36
    w2 = w4 = (18 - math.sqrt(30))/36
    return ((xf-xi)/2)*(f(x_a1)*w1 + f(x_a2)*w2 + f(x_a3)*w3 + f(x_a4)*w4)

# ======================================= #
#     Funções de integração numerica      #
# ======================================= #

def integrate_n(f, a, b, n, subs_func = regra_trapezio):
    I = 0
    h = (b-a)/n

    for k in range(n):
        x0 = a + k*h
        x1 = x0 + h
        I += subs_func(f, x0, x1, k, h)
    
    return I

def integrate(f, a, b, E, subs_func=regra_trapezio, debug=False, dir=""):
    I = 0
    N = 1
    erro = 1
    i = 0
    
    if(debug):
        path_to_output = os.path.join("./out", dir, f"{subs_func.__name__}.txt") 
        file = open(path_to_output, "w")

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

# ======================================= #
#    Gauss-{Hermite,Laguerre,Chebyshev}   #
# ======================================= #

def gauss_hermite_n2(f):
    x1 = -1/math.sqrt(2)
    x2 = 1/math.sqrt(2)
    w1 = w2 = math.sqrt(math.pi)/2
    return w1*f(x1) + w2*f(x2)

def gauss_hermite_n3(f):
    x1 = -1*math.sqrt(3/2)
    x2 = 0
    x3 = math.sqrt(3/2)
    w1 = w3 = math.sqrt(math.pi)/6
    w2 = 2*math.sqrt(math.pi)/3
    return w1*f(x1) + w2*f(x2) + w3*f(x3)

def gauss_hermite_n4(f):
    x1 = -1*math.sqrt( (3+math.sqrt(6))/2 )
    x2 = -1*math.sqrt( (3-math.sqrt(6))/2 )
    x3 = math.sqrt( (3-math.sqrt(6))/2 )
    x4 = math.sqrt( (3+math.sqrt(6))/2 )
    w1 = w4 = (math.sqrt(9*math.pi) - math.sqrt(6*math.pi))/12
    w2 = w3 = (math.sqrt(9*math.pi) + math.sqrt(6*math.pi))/12
    return w1*f(x1) + w2*f(x2) + w3*f(x3) + w1*f(x4)

def gauss_laguerre_n2(f):
    x1 = 2 - math.sqrt(2)
    x2 = 2 + math.sqrt(2)
    w1 = (1/4) * (2 + math.sqrt(2))
    w2 = (1/4) * (2 - math.sqrt(2))
    return w1*f(x1) + w2*f(x2)

def gauss_laguerre_n3(f):
    x1 = 0.4157745568
    x2 = 2.2942803603
    x3 = 6.2899450829
    w1 = 0.7110930099
    w2 = 0.2785177336
    w3 = 0.0103892565
    return w1*f(x1) + w2*f(x2) + w3*f(x3)

def gauss_laguerre_n4(f): 
    x1 = 0.322548
    x2 = 1.74576
    x3 = 4.53662
    x4 = 9.39507
    w1 = 0.603154
    w2 = 0.357419
    w3 = 0.0388879
    w4 = 0.000539295
    return w1*f(x1) + w2*f(x2) + w3*f(x3) + w4*f(x4)

def gauss_chebyshev_n2(f):
    x1 = -1/math.sqrt(2)
    x2 = 1/math.sqrt(2)
    w = math.pi/2
    return w*f(x1) + w*f(x2)

def gauss_chebyshev_n3(f):
    x1 = -math.sqrt(3)/2
    x2 = 0
    x3 = math.sqrt(3)/2
    w = math.pi/3
    return w*f(x1) + w*f(x2) + w*f(x3)

def gauss_chebyshev_n4(f):
    x1 = -1*math.sqrt( (2+math.sqrt(2))/4 )
    x2 = -1*math.sqrt( (2-math.sqrt(2))/4 )
    x3 = math.sqrt( (2-math.sqrt(2))/4 )
    x4 = math.sqrt( (2+math.sqrt(2))/4 )
    w = math.pi/4
    return w*f(x1) + w*f(x2) + w*f(x3) + w*f(x4)