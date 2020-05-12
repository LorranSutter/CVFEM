import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#--------- Escrito em Python 3 ---------#

# Le o arquivo que contem os seguintes dados:
#    Numero de nos
#    Coordendadas dos nos
#    Vetor contendo o numero de suportes por no
#    Numero maximo de suportes por no
#    Matriz de suportes
#    Numero maximo de nos por contorno
#    Vetor contendo o numero de nos por contorno
#    Matriz de contorno


def read_file(arq):
    print('Iniciada a leitura do arquivo de entrada')

    # Abre o arquivo
    f = open(arq, 'r')

    # Armazena o numero de nos
    n = int(f.readline())
    print('Numero total de nos: ', n)

    # Inicializacao dos vetores de coordenadas
    x = np.zeros(n)
    y = np.zeros(n)

    for k in range(n):
        x[k], y[k] = list(map(float, f.readline().split()))

    # Inicializacao do vetor que contem o numero de suportes por no
    n_s = np.zeros(n, dtype=np.int64)

    for k in range(n):
        n_s[k] = int(f.readline())

    # Numero maximo de suportes por no
    num_suport = int(f.readline()) + 1

    # Inicializacao da matriz de suportes
    S = np.zeros((n, num_suport))

    print('Iniciada a criacao da matriz de suportes (' +
          str(n) + ',' + str(num_suport) + ')')
    for k in range(n):
        S[k, :] = list(map(int, (f.readline().split())))

    num_seg = int(f.readline())  # boundary segments

    n_b = np.zeros(num_seg)

    for k in range(num_seg):
        n_b[k] = int(f.readline())

    print('Iniciada a criacao da matriz de contorno (' +
          str(num_seg) + ',' + str(int(max(n_b))) + ')')
    # Global node number of j_th node on i_th boundary
    B = np.zeros((num_seg, int(max(n_b))))

    for k in range(num_seg):
        aux = list(map(int, f.readline().split()))
        for w in range(len(aux)):
            B[k, w] = aux[w]

    print('Finalizada a criacao da matriz de contorno')
    f.close()
    print('Finilizada a leitura do arquivo')

    return [n, x, y, n_s, num_suport, S, num_seg, n_b, B]

# Define os coeficientes a_i e a_i,j


def set_coefficients(n, x, y, vx, vy, S, n_s):
    ap = np.zeros(n)  # Coeficiente do no da iteracao
    a = np.zeros((n, num_suport))  # Coeficiente dos suportes do no da iteracao
    V = np.zeros(n)  # Volume do VC
    xL = np.zeros(3)  # Coordenada Local x dos 3 nos no elemento da iteracao
    yL = np.zeros(3)  # Coordenada Local y dos 3 nos no elemento da iteracao
    # Difusividade nodal Local dos 3 nos no elemento da iteracao
    dif_nodal_L = np.zeros(3)
    # Velocidade nodal Local x dos 3 nos no elemento da iteracao
    vxL = np.zeros(3)
    # Velocidade nodal Local y dos 3 nos no elemento da iteracao
    vyL = np.zeros(3)
    Vele = 0  # Volume do elemento da iteracao
    V = np.zeros(n)  # Volume de controle do no da iteracao
    # Derivada em x da funcao de forma dos 3 nos no elemento da iteracao
    Nx = np.zeros(3)
    # Derivada em y da funcao de forma dos 3 nos no elemento da iteracao
    Ny = np.zeros(3)
    dif_face = 0  # Difusividade no ponto medio da face
    vx_face = 0  # Velocidade x no ponto medio da face
    vy_face = 0  # Velocidade y no ponto medio da face
    delta_x = 0  # Componente x normal da area na face
    delta_y = 0  # Componente y normal da area na face
    q_face = 0  # Fluxo de volume na face

    for k in range(n):
        if np.mod(k, 100) == 0:
            print('Iteracao: ', k, '\t Faltam: ', n-k)

        loopcount = int(n_s[k])  # Numero de suportes por no
        print(k)
        if S[k, n_s[k]] == 0:  # Detecta se eh contorno
            loopcount = int(n_s[k]-1)

        xL[0] = x[k]  # Coordenda Local x do no da iteracao
        yL[0] = y[k]  # Coordenda Local y do no da iteracao
        dif_nodal_L[0] = dif_nodal[k]  # Difusividade nodal do no da iteracao
        vxL[0] = vx[k]
        vyL[0] = vy[k]

        # No da iteracao e suportes
        #           k,w+2 0------0 k,w+1
        #               / |     /|
        #              /  |    / |
        #             /   |   /  |
        #            /    |  /   |
        #           /     i /    |
        #    k,w+3 0------0------0 k,w; k,w+6
        #          |     /|     /
        #          |    / |    /
        #          |   /  |   /
        #          |  /   |  /
        #          | /    | /
        #    k,w+4 0------0 k,w+5

        # Loop para cada suporte do no da iteracao
        for w in range(loopcount):
            # Itera em pares de nos suportes no sentido anti-horario
            k1 = int(S[k, w])
            k2 = int(S[k, w+1])

            # Coordendas locais dos nos suportes do no da iteracao
            xL[1] = x[k1-1]
            yL[1] = y[k1-1]
            xL[2] = x[k2-1]
            yL[2] = y[k2-1]

            # Difusividade nodal Local dos nos suportes do no da iteracao
            dif_nodal_L[1] = dif_nodal[k1-1]
            dif_nodal_L[2] = dif_nodal[k2-1]

            # Velocidade nodal Local dos nos suportes do no da iteracao
            vxL[1] = vx[k1-1]
            vyL[1] = vy[k1-1]
            vxL[2] = vx[k2-1]
            vyL[2] = vy[k2-1]

            # Volume do elemento da iteracao
            Vele = (xL[1]*yL[2]-xL[2]*yL[1]-xL[0]*yL[2] +
                    xL[0]*yL[1]+yL[0]*xL[2]-yL[0]*xL[1])/2

            # Contribuicao para o VC, usa-se em situacoes que ha fonte
            V[k] = V[k] + Vele/3

            # Derivadas das funcoes de forma do elemento
            Nx[0] = (yL[1]-yL[2])/(2*Vele)
            Nx[1] = (yL[2]-yL[0])/(2*Vele)
            Nx[2] = (yL[0]-yL[1])/(2*Vele)
            Ny[0] = -(xL[1]-xL[2])/(2*Vele)
            Ny[1] = -(xL[2]-xL[0])/(2*Vele)
            Ny[2] = -(xL[0]-xL[1])/(2*Vele)

            # -------------------- FACE 1 --------------------

            # Difusividade no ponto medio da face
            dif_face = (5*dif_nodal_L[0] + 5 *
                        dif_nodal_L[1] + 2*dif_nodal_L[2])/12

            # Velocidade no ponto medio da face
            vx_face = (5*vxL[0] + 5*vxL[1] + 2*vxL[2])/12
            vy_face = (5*vyL[0] + 5*vyL[1] + 2*vyL[2])/12

            # Componente normal da area na face
            delta_x = ((xL[0]+xL[1]+xL[2])/3) - ((xL[0]+xL[1])/2)
            delta_y = ((yL[0]+yL[1]+yL[2])/3) - ((yL[0]+yL[1])/2)

            # Fluxo de volume na face
            q_face = vx_face*delta_y - vy_face*delta_x

            # Coeficientes de difusao
            ap[k] = ap[k] + dif_face*(-Nx[0] * delta_y + Ny[0] * delta_x)
            a[k, w] = a[k, w] + dif_face*(Nx[1] * delta_y - Ny[1] * delta_x)
            a[k, w+1] = a[k, w+1] + dif_face * \
                (Nx[2] * delta_y - Ny[2] * delta_x)

            # Coeficientes de adveccao
            ap[k] = ap[k] + max(q_face, 0)
            a[k, w] = a[k, w] + max(-q_face, 0)

            # -------------------- FACE 2 --------------------

            # Difusividade no ponto medio da face
            dif_face = (5*dif_nodal_L[0] + 2 *
                        dif_nodal_L[1] + 5*dif_nodal_L[2])/12

            # Velocidade no ponto medio da face
            vx_face = (5*vxL[0] + 2*vxL[1] + 5*vxL[2])/12
            vy_face = (5*vyL[0] + 2*vyL[1] + 5*vyL[2])/12

            # Componente normal da area na face
            delta_x = ((xL[0]+xL[2])/2) - ((xL[0]+xL[1]+xL[2])/3)
            delta_y = ((yL[0]+yL[2])/2) - ((yL[0]+yL[1]+yL[2])/3)

            # Fluxo de volume na face
            q_face = vx_face*delta_y - vy_face*delta_x

            # Coeficientes de difusao
            ap[k] = ap[k] + dif_face*(-Nx[0] * delta_y + Ny[0] * delta_x)
            a[k, w] = a[k, w] + dif_face*(Nx[1] * delta_y - Ny[1] * delta_x)
            a[k, w+1] = a[k, w+1] + dif_face * \
                (Nx[2] * delta_y - Ny[2] * delta_x)

            # Coeficientes de adveccao
            ap[k] = ap[k] + max(q_face, 0)
            a[k, w+1] = a[k, w+1] + max(-q_face, 0)

        # endfor loopcount

        a[k, 0] = a[k, 0] + a[k, n_s[k]]
        a[k, n_s[k]] = 0

    # endfor n

    return [a, ap, V]

# Define os contornos


def set_boundary(B, n_b):
    BC = np.zeros(n)
    BB = np.zeros(n)

    # Para este caso especifico
    for k in range(int(n_b[0])):
        BC[int(B[2, k])-1] = 1e18

    for k in range(int(n_b[2])):
        BC[int(B[0, k])-1] = 1e18
        BB[int(B[0, k])-1] = 1e18

    return [BC, BB]

# Define a fonte, para este caso


def set_source():
    QC = np.zeros(n)
    QB = np.zeros(n)
    return [QC, QB]


def solver(n, S, n_s, a, ap, BC, BB, QC, QB, maxIt, tolerancia):
    '''
    Equacao 5.28, Voller
    '''
    phi_old = np.zeros(n)
    phi_new = np.zeros(n)
    RHS = 0.0

    for it in range(maxIt):
        if(np.mod(it, 100) == 0):
            print('Solver, iteracao: ', it, '\t Faltam: ', maxIt-it)
        for k in range(n):
            RHS = BB[k] + QB[k]
            for w in range(int(n_s[k])):
                RHS = RHS + a[k, w]*phi_new[int(S[k, w])-1]
            phi_new[k] = float(RHS/(ap[k]+BC[k]+QC[k]))

        if(max(abs(phi_new-phi_old)) <= tolerancia):
            break
        phi_old = phi_new.copy()
    # endfor it

    return phi_new


if __name__ == "__main__":
    arq = sys.argv[1]
    maxIt = 1500
    tolerancia = 0.001

    [n, x, y, n_s, num_suport, S, num_seg, n_b, B] = read_file(arq)

    # ------------ Parte, por hora, necessaria de conversao em coordenada polar ------------
    radius = np.zeros(n)
    theta = np.zeros(n)
    for k in range(n):
        radius[k] = math.sqrt(x[k]**2 + y[k]**2)
        theta[k] = np.arctan2(y[k], x[k])

    vx = np.zeros(n)
    vy = np.zeros(n)
    for k in range(n):
        vx[k] = math.cos(theta[k])/radius[k]
        vy[k] = math.sin(theta[k])/radius[k]

    dif_nodal = np.zeros(n)
    for k in range(n):
        dif_nodal[k] = 1/math.sqrt(x[k]**2 + y[k]**2)
    # -------------- Fim da parte necessaria de conversao em coordenada polar --------------

    print('Iniciando o set dos coeficientes de adveccao e difusao')
    [a, ap, V] = set_coefficients(n, x, y, vx, vy, S, n_s)
    print('Finalizado o set dos coeficientes de adveccao e difusao')

    print('Iniciado o set das condicoes de contorno')
    [BC, BB] = set_boundary(B, n_b)
    print('Finalizado o set das condicoes de contorno')

    print('Iniciado o set da fonte')
    [QC, QB] = set_source()
    print('Finalizado o set da fonte')

    print('Iniciado o SOLVER')
    phi = solver(n, S, n_s, a, ap, BC, BB, QC, QB, maxIt, tolerancia)
    print('Finalizado o SOLVER')

    pout = 0
    xout = []
    phiout = []
    for k in range(n):
        if y[k] < 1e-5:
            pout = pout + 1
            xout.append(x[k])
            phiout.append(phi[k])

    ana = np.zeros(pout)
    for k in range(pout):
        ana[k] = (math.exp(xout[k]) - math.exp(2))/(math.exp(1) - math.exp(2))

    plt.plot(xout, ana, 'ro', label='Analitico')
    plt.plot(xout, phiout, 'bo', label='Numerico')
    plt.legend(loc='best')
    plt.show()
