import math
import os
import numpy as np
# ======================================= #
#       Fórmulas de newton-cotes          #
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
    h = delta_x / 4.0
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

def integrate_n(f, a, b, n, subs_func=regra_trapezio):
    I = 0
    delta = float((b-a)/float(n))

    for k in range(n):
        x0 = a + k*delta
        x1 = x0 + delta
        I += subs_func(f, x0, x1, k, delta)
    
    return I

def integrate(f, a, b, E, subs_func=regra_trapezio, debug=False, dir=""):
    I = 0
    N = 1
    erro = 1
    i = 0
    
    if(debug):
        path_to_output = os.path.join("./out", dir, f"{subs_func.__name__}.txt") 
        file = open(path_to_output, "w")

    while erro >= E:
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

# ======================================= #
#        Problema de singularidade        #
# ======================================= #

def x_expo_simples(a, b, s):
    x = ((a+b)/2.0) + (((b-a)/2.0)*math.tanh(s))
    return x

def ds_expo_simples(a, b, s):
    return (b-a)/2 * (1/(math.cosh(s)**2))

def x_expo_dupla(a, b, s):
    return ((a+b)/2.0) + ((b-a)/2.0 * math.tanh(math.pi/2 * math.sinh(s)))

def ds_expo_dupla(a, b, s):
    return (b-a)/2.0 * (math.pi/2 * math.cosh(s)/math.pow(math.cosh(math.pi/2*math.sinh(s)),2.0))


def exponencial_simples(f, a, b, E):
    erro = 1
    I = 0
    c = 1

    def g(x):
        return f(x_expo_simples(a, b, x)) * ds_expo_simples(a, b, x)
    
    # c é no máximo 19 para evitar o erro de precisão por conta do limite de ponto flutuante
    while erro >= E and c <= 19:
        In = integrate(g, -c, c, 10**(-7), subs_func=regra_poli_4_fechada)
        print(In)
        erro = abs((In - I)/ In)
        I = In
        c += 1
    return I

def exponencial_dupla(f, a, b, E):
    erro = 1
    I = 0
    c = 1
    
    def g(x):
        return f(x_expo_dupla(a, b, x)) * ds_expo_dupla(a, b, x)
    
    # c é no máximo 3 para evitar o erro de precisão por conta do limite de ponto flutuante
    while erro >= E and c <= 3:
        In = integrate(g, -c, c, 10**(-7), subs_func=gauss_legendre_4p)
        erro = abs((In - I)/ In)
        I = In
        c += .5
    return I


# ======================================= #
#        Autovalores e Autovetores        #
# ======================================= #

def normalizar(v):
    out = []
    div = math.sqrt(mult_vector(v,v))
    for n in v:
        out.append(n/div)
    return out

def mult_vector(v0, v1):
    out = 0
    for n,m in zip(v0,v1):
        out += n*m
    return out

def mult_matrix_vector(A, v):
    n = len(A)
    out = []
    for i in range(n):
        s = 0.0
        for j in range(n):
            s += v[j]*A[i][j]
        out.append(s)
    return out

def mult_matrix_matrix(A, B):
    out = [[0 for i in range(len(B[0]))] for j in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                out[i][j] += A[i][k] * B[k][j]
    return out

def getIdentityMatrix(n):
    return [[1 if i == j else 0 for i in range(n)] for j in range(n)]

def print_matrix(matrix):
    n = len(matrix)
    for i in range(n):
        print('\t', end='')
        for j in range(len(matrix[0])):
            print(f'{matrix[i][j]:.5f}'.zfill(8), end=' ')
        print()

def decomLU(A):
    n = len(A)
    L = [[0 for x in range(n)] for y in range(n)]
    U = [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        
        for k in range(i,n):
            soma = 0
            for j in range(i):
                soma += (L[i][j] * U[j][k])
            U[i][k] = A[i][k] - soma
            if(U[i][i] == 0):
                raise ZeroDivisionError
        
        for k in range(i,n):
            if(i == k): 
                L[i][i] = 1
            else:
                soma = 0
                for j in range(i):
                    soma += L[k][j] * U[j][i]
                L[k][i] = (A[k][i] - soma) / U[i][i]
    return L, U

def subs_sucessivas(matrix, vector):
    n = len(matrix)
    out = [0 for i in range(n)]
    soma = 0
    out[0] = vector[0] / matrix[0][0]

    for i in range(1, n):
        soma = 0
        for j in range(i):
            soma += matrix[i][j] * out[j]
        out[i] = (vector[i] - soma) / matrix[i][i]
    
    return out

def subs_retroativas(matrix, vector):
    n = len(matrix)
    out = [0 for i in range(n)]
    soma = 0
    out[n-1] = vector[n-1]/matrix[n-1][n-1]
    
    for i in range(n-2, -1, -1):
        soma = 0
        for j in range(i+1, n):
            soma += matrix[i][j] * out[j]
        out[i] = (vector[i] - soma) / matrix[i][i]
    
    return out

def solveLU(L, U, x):
    y = subs_sucessivas(L, x)
    return subs_retroativas(U, y)

def potencia_regular(A, v0, e):
    gamma_novo = 0
    vk_novo = v0
    error = 10
    while(error > e):
        
        gamma_velho = gamma_novo
        vk_velho = vk_novo

        x = normalizar(vk_velho)
        vk_novo = mult_matrix_vector(A, x)
        gamma_novo = mult_vector(x, vk_novo)
        error = abs((gamma_novo - gamma_velho)/gamma_novo)
    
    return gamma_novo, x

def potencia_inversa(A, v0, e):
    
    L, U = decomLU(A)
    gamma_novo = 0
    vk_novo = v0
    error = 10
    
    while(error > e):
        gamma_velho = gamma_novo
        vk_velho = vk_novo
        x = normalizar(vk_velho)

        vk_novo = solveLU(L, U, x)
        gamma_novo = mult_vector(x, vk_novo)
        error = abs((gamma_novo - gamma_velho)/gamma_novo)
    
    gamma = float(1.0/gamma_novo)

    return gamma, x

def potencia_com_deslocamento(A, v0, e, u):
    n = len(A)
    A_temp = [[A[i][j] if i != j else A[i][j] - u for j in range(n)] for i in range(n)]
                
    gamma, x = potencia_inversa(A_temp, v0, e)
    gamma += u
    return gamma, x

def select_pivo(matrix, k):
    n = len(matrix)
    pv = matrix[k][k]
    index = k
    for i in range(k+1, n):
        pv_aux = matrix[i][k]
        if(abs(pv_aux) > pv):
            pv = pv_aux
            index = i
    return pv, index

def permute(matrix, p, k, r):
    n = len(matrix)
    p[k], p[r] = p[r], p[k]
    for i in range(n):
        matrix[k][i], matrix[r][i] = matrix[r][i], matrix[k][i]

def fact_LU(A):
    n = len(A)
    U = A[:]
    L = [[0 for i in range(n)] for j in range(n)]
    permutations = [i for i in range(n)]
    L[0][0] = 1

    for k in range(n-1):
        
        pv, index = select_pivo(U, k)
        if(pv == 0): raise ZeroDivisionError
        if(index != k): permute(U, permutations, k, index)
   
        for i in range(k+1, n):
            m = -(U[i][k] / U[k][k])
            U[i][k] = 0
            L[i][k] = -m
            
            for j in range(k+1, n):
                U[i][j] = U[i][j] + (m * U[k][j])
                L[i][j] = 1 if (i == j) else 0
    return L, U 

def tranpose_matrix(A):
    n = len(A) 
    return [[A[i][j] for i in range(n)] for j in range(len(A[0]))]

def get_lenght(v):
    return math.sqrt(mult_vector(v, v))

def subract_vector(v1, v2):
    n = len(v1)
    out = []
    for n, m in zip(v1, v2):
        out.append(n - m)
    return out

# ======================================= #
#          Metodo de Householder          #
# ======================================= #

def get_matriz_householder(A, i):
    n = len(A)
    w = [0 if k <= i else A[k][i] for k in range(n)]
    w2 = [0 for i in range(n)]

    w2[i + 1] = get_lenght(w)
    
    N = subract_vector(w, w2)
    
    n_norm = normalizar(N)
    
    temp_n = [2*n_norm[i] for i in range(n)]
    temp = mult_matrix_matrix(tranpose_matrix([temp_n]), [n_norm])
    
    H = [[1 - temp[i][j] if i == j else -temp[i][j] for j in range(n)] for i in range(n)]
    
    return H

def metodo_de_householder(A, n):
    H = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    A_old = A[:]
    for i in range(n-2):
        
        H_new = get_matriz_householder(A_old, i)

        temp = mult_matrix_matrix(H_new, A_old)
        A_new = mult_matrix_matrix(temp, H_new)
        
        A_old = A_new
        
        H = mult_matrix_matrix(H, H_new)

    return A_new, H

# ======================================= #
#            Metodo de Jacobi             #
# ======================================= #

def getMatrizJacobiBaseadaEmIJ(A, i, j, n):
    e = 10**-6

    J = getIdentityMatrix(n)

    if(abs(A[i][j]) <= e):
        return J
    elif(abs(A[i][i] - A[j][j]) <= e):
        theta = math.pi / 4
    else:
        theta = 0.5*math.atan(-2*A[i][j] / (A[i][i] - A[j][j]))
    
    J[i][i] = math.cos(theta)
    J[j][j] = math.cos(theta)
    J[i][j] = math.sin(theta)
    J[j][i] = -math.sin(theta)

    return J

def varreduraDeJacobi(A, n, DEBUG=False):
    J = getIdentityMatrix(n)

    A_old = A

    if(DEBUG):
        print("\nVarredura de Jacobi")

    for j in range(n-1):
        for i in range(j+1, n):
            J_new = getMatrizJacobiBaseadaEmIJ(A_old, i, j, n)

            J_newt = tranpose_matrix(J_new)
            temp = mult_matrix_matrix(J_newt, A_old)
            A_new = mult_matrix_matrix(temp, J_new)

            A_old = A_new
            J = mult_matrix_matrix(J, J_new)
            
            if(DEBUG):
                print("Matriz A_nova = {")
                print_matrix(A_new)
                print("}")
    
    A_barra = A_new
    return A_barra, J

def somaAbaixoDaDiagonal(A, n):
    soma = 0
    for j in range(n-1):
        for i in range(j+1, n):
            soma += A[i][j]**2
    return soma

def metodoDeJacobi(A, n, e, DEBUG=False, PRINT_TRI=False):
    P = getIdentityMatrix(n)
    A_old = A[:]
    val = 100.0

    while(val > e):
        A_new, J = varreduraDeJacobi(A_old, n, PRINT_TRI)
        if(DEBUG):
            print("Varredura = {")
            print_matrix(A_new)
            print("}")
        A_old = A_new
        P = mult_matrix_matrix(P, J)
        val = somaAbaixoDaDiagonal(A_new, n)

    Lamb = [A_new[i][i] for i in range(n)]
    return P, Lamb

# ======================================= #
#                Metodo QR                #
# ======================================= #

def getMatrizJacobiBaseadaEmIJ_De_R(R, i, j, n):
    e = 10**-6
    J = getIdentityMatrix(n)
    if(abs(R[i][j]) <= e):
        return J
    if(abs(R[j][j] <= e)):
        if(R[i][j] < 0):
            theta = math.pi/2
        else:
            theta = -math.pi/2
    else:
        theta = math.atan(-R[i][j]/R[j][j])
    
    J[i][i] = math.cos(theta)
    J[j][j] = math.cos(theta)
    J[i][j] = math.sin(theta)
    J[j][i] = -math.sin(theta)

    return J

def decompQR(A, n):
    QT = getIdentityMatrix(n)
    R_old = A
    for j in range(n-1):
        for i in range(j+1, n):
            J_new = getMatrizJacobiBaseadaEmIJ_De_R(R_old, i, j , n)
            R_new = mult_matrix_matrix(J_new, R_old)
            R_old = R_new
            QT = mult_matrix_matrix(J_new, QT)
    Q = tranpose_matrix(QT)
    R = R_new
    return Q, R

def metodoQR(A, n, e, DEBUG=False):
    val = 100

    P = getIdentityMatrix(n)
    A_old = A[:]
    
    while(val > e):
        Q, R = decompQR(A_old, n)
        A_new = mult_matrix_matrix(R, Q)
        A_old = A_new
        if(DEBUG):
            print("Matriz QR = {")
            print_matrix(A_new)
            print("}")
        P = mult_matrix_matrix(P, Q)
        val = somaAbaixoDaDiagonal(A_new, n)
    autovetores = [A_new[i][i] for i in range(n)]
    return P, autovetores

# ======================================= #
#          Problema Valor Inicial         #
# ======================================= #

def euler_explicito(S, F, delta_t):
    return S + delta_t * F(S)

def runge_kutta(S, F, delta_t):
    F1 = F(S[0])
    S_1 = S + delta_t/2 * F1

    F2 = F(S_1[0])
    S_2 = S + delta_t * (-1*F1 + 2*F2)

    F3 = F(S_2[0])

    S_out = S + delta_t * (1/6*F1 + 4/6*F2 + 1/6*F3)
    return S_out

def runge_kutta_4_ordem(S, F, delta_t):
    F1 = F(S[0])
    S_1 = S + delta_t/2 * F1

    F2 = F(S_1[0])
    S_2 = S + delta_t/2 * F2

    F3 = F(S_2[0])
    S_3 = S + delta_t * F3

    F4 = F(S_3[0])

    S_out = S + delta_t/6 * (F1 + 2*F2 + 2*F3 + F4)
    return S_out

def preditor_corretor_4_ordem(S, F, delta_t, values):

    error = 10**-7
    F_1 = F(values[0][0])
    F_2 = F(values[1][0])
    F_3 = F(values[2][0])
    F_4 = F(S[0])

    # Predição
    S_pred = S + delta_t/24 * (-9*F_1 + 37*F_2 - 59*F_3 + 55*F_4)

    # Correção
    F_new = F(S_pred[0])
    S_correction = S + delta_t/24 * (F_2 - 5*F_3 + 19*F_4 + 9*F_new)

    S_error = abs((S_correction - S_pred)/S_correction)
    while(S_error[0] > error and S_error[1] > error):
        S_pred = S_correction[:]
        F_new = F(S_pred[0])
        S_correction = S + delta_t/24 * (F_2 - 5*F_3 + 19*F_4 + 9*F_new)
        S_error = abs((S_correction - S_pred)/S_correction)
    
    return S_correction