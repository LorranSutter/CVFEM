import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def read_file(arq):
    '''
    Read the file containing the following data
        - Number of nodes
        - Nodes coordinates
        - Array of number of supports by node
        - Max number of supports by node
        - Support matrix
        - Max number of nodes by boundary
        - Array of number of nodes by boundary
        - Boundary matrix

    '''
    print('\nStarting input file reading')

    # Open the file
    f = open(arq, 'r')

    # Store number of nodes
    n = int(f.readline())
    print('Total number of nodes: ', n)

    # Init coordinates array
    x = np.zeros(n)
    y = np.zeros(n)

    for k in range(n):
        x[k], y[k] = list(map(float, f.readline().split()))

    # Init the array that contains the number of supports by node
    n_s = np.zeros(n, dtype=np.int64)

    for k in range(n):
        n_s[k] = int(f.readline())

    # Max number of supports by node
    num_suport = int(f.readline()) + 1

    # Init support matrix
    S = np.zeros((n, num_suport))

    print('Started creating the support matrix (' +
          str(n) + ',' + str(num_suport) + ')')
    for k in range(n):
        S[k, :] = list(map(int, (f.readline().split())))
    print('Finished creating the support matrix')

    # Boundary segments
    num_seg = int(f.readline())

    n_b = np.zeros(num_seg)

    for k in range(num_seg):
        n_b[k] = int(f.readline())

    print('Started creating boundary matrix (' +
          str(num_seg) + ',' + str(int(max(n_b))) + ')')

    # Global node number of j_th node on i_th boundary
    B = np.zeros((num_seg, int(max(n_b))))

    for k in range(num_seg):
        aux = list(map(int, f.readline().split()))
        for w in range(len(aux)):
            B[k, w] = aux[w]
    print('Finished creating the boundary matrix')

    f.close()
    print('Finished reading the file\n')

    return [n, x, y, n_s, num_suport, S, num_seg, n_b, B]


def set_coefficients(n, x, y, vx, vy, S, n_s):
    '''
    Defines the coefficients a_i and a_i,j
    '''
    ap = np.zeros(n)              # Iteration node coefficient
    a = np.zeros((n, num_suport)) # Coefficient of iteration node supports
    V = np.zeros(n)               # Volume in the CV
    xL = np.zeros(3)              # Local coordinate x of the 3 nodes in the iteration element
    yL = np.zeros(3)              # Local coordinate y of the 3 nodes in the iteration element
    dif_nodal_L = np.zeros(3)     # Local nodal diffusivity of the 3 nodes in the iteration element
    vxL = np.zeros(3)             # Local nodal velocity x of the 3 nodes in the iteration element
    vyL = np.zeros(3)             # Local nodal velocity y of the 3 nodes in the iteration element
    Vele = 0                      # Volume of the iteration element
    V = np.zeros(n)               # Control volume of the iteration node
    Nx = np.zeros(3)              # Derivate in x of the shape function of the 3 nodes in the iteration element
    Ny = np.zeros(3)              # Derivate in y of the shape function of the 3 nodes in the iteration element
    diff_face = 0                 # Diffusivity in the midpoint of the face
    vx_face = 0                   # Velocity x in the midpoint of the face
    vy_face = 0                   # Velocity y in the midpoint of the face
    delta_x = 0                   # Face area normal component x
    delta_y = 0                   # Face area normal component y
    q_face = 0                    # Face volume flow

    for k in range(n):
        if np.mod(k, 100) == 0:
            print('Iteration: ', k, '\t To finish: ', n-k)

        loopcount = int(n_s[k])  # Number of supports by node
        if S[k, n_s[k]] == 0:  # Detects if is a boundary
            loopcount = int(n_s[k]-1)

        xL[0] = x[k]                   # Local coordinate x of the iteration node
        yL[0] = y[k]                   # Local coordinate y of the iteration node
        dif_nodal_L[0] = dif_nodal[k]  # Nodal diffusivity in the iteration node
        vxL[0] = vx[k]
        vyL[0] = vy[k]

        # Iteration node and supports
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

        # Loop for each iteration node support
        for w in range(loopcount):
            # Iterate in pairs of support nodes counterclockwise
            k1 = int(S[k, w])
            k2 = int(S[k, w+1])

            # Local coordinates of the supporting nodes of the iteration node
            xL[1] = x[k1-1]
            yL[1] = y[k1-1]
            xL[2] = x[k2-1]
            yL[2] = y[k2-1]

            # Local nodal diffusivity of the nodes supporting the iteration node
            dif_nodal_L[1] = dif_nodal[k1-1]
            dif_nodal_L[2] = dif_nodal[k2-1]

            # Local nodal speed of iteration node support nodes
            vxL[1] = vx[k1-1]
            vyL[1] = vy[k1-1]
            vxL[2] = vx[k2-1]
            vyL[2] = vy[k2-1]

            # Volume of the iteration element
            Vele = (xL[1]*yL[2]-xL[2]*yL[1]-xL[0]*yL[2] +
                    xL[0]*yL[1]+yL[0]*xL[2]-yL[0]*xL[1])/2

            # Contribution to CV, used in situations that have a source
            V[k] = V[k] + Vele/3

            # Derivates of the element shape functions
            Nx[0] = (yL[1]-yL[2])/(2*Vele)
            Nx[1] = (yL[2]-yL[0])/(2*Vele)
            Nx[2] = (yL[0]-yL[1])/(2*Vele)
            Ny[0] = -(xL[1]-xL[2])/(2*Vele)
            Ny[1] = -(xL[2]-xL[0])/(2*Vele)
            Ny[2] = -(xL[0]-xL[1])/(2*Vele)

            # -------------------- FACE 1 --------------------

            # Diffusivity at the midpoint of the face
            diff_face = (5*dif_nodal_L[0] + 5 *
                        dif_nodal_L[1] + 2*dif_nodal_L[2])/12

            # Velocity at the midpoint of the face
            vx_face = (5*vxL[0] + 5*vxL[1] + 2*vxL[2])/12
            vy_face = (5*vyL[0] + 5*vyL[1] + 2*vyL[2])/12

            # Normal area component on the face
            delta_x = ((xL[0]+xL[1]+xL[2])/3) - ((xL[0]+xL[1])/2)
            delta_y = ((yL[0]+yL[1]+yL[2])/3) - ((yL[0]+yL[1])/2)

            # Face volume flow
            q_face = vx_face*delta_y - vy_face*delta_x

            # Diffusion coefficients
            ap[k] = ap[k] + diff_face*(-Nx[0] * delta_y + Ny[0] * delta_x)
            a[k, w] = a[k, w] + diff_face*(Nx[1] * delta_y - Ny[1] * delta_x)
            a[k, w+1] = a[k, w+1] + diff_face * \
                (Nx[2] * delta_y - Ny[2] * delta_x)

            # Advection coefficients
            ap[k] = ap[k] + max(q_face, 0)
            a[k, w] = a[k, w] + max(-q_face, 0)

            # -------------------- FACE 2 --------------------

            # Diffusivity at the midpoint of the face
            diff_face = (5*dif_nodal_L[0] + 2 *
                        dif_nodal_L[1] + 5*dif_nodal_L[2])/12

            # Velocity at the midpoint of the face
            vx_face = (5*vxL[0] + 2*vxL[1] + 5*vxL[2])/12
            vy_face = (5*vyL[0] + 2*vyL[1] + 5*vyL[2])/12

            # Normal area component on the face
            delta_x = ((xL[0]+xL[2])/2) - ((xL[0]+xL[1]+xL[2])/3)
            delta_y = ((yL[0]+yL[2])/2) - ((yL[0]+yL[1]+yL[2])/3)

            # Face volume flow
            q_face = vx_face*delta_y - vy_face*delta_x

            # Diffusion coefficients
            ap[k] = ap[k] + diff_face*(-Nx[0] * delta_y + Ny[0] * delta_x)
            a[k, w] = a[k, w] + diff_face*(Nx[1] * delta_y - Ny[1] * delta_x)
            a[k, w+1] = a[k, w+1] + diff_face * \
                (Nx[2] * delta_y - Ny[2] * delta_x)

            # Advection coefficients
            ap[k] = ap[k] + max(q_face, 0)
            a[k, w+1] = a[k, w+1] + max(-q_face, 0)

        a[k, 0] = a[k, 0] + a[k, n_s[k]]
        a[k, n_s[k]] = 0

    return [a, ap, V]

def set_boundary(B, n_b):
    '''
    Define boundaries
    '''
    BC = np.zeros(n)
    BB = np.zeros(n)

    # For this specific case
    for k in range(int(n_b[0])):
        BC[int(B[2, k])-1] = 1e18

    for k in range(int(n_b[2])):
        BC[int(B[0, k])-1] = 1e18
        BB[int(B[0, k])-1] = 1e18

    return [BC, BB]

def set_source():
    '''
    Defines the source, for this case
    '''
    QC = np.zeros(n)
    QB = np.zeros(n)
    return [QC, QB]


def solver(n, S, n_s, a, ap, BC, BB, QC, QB, maxIt, tolerance):
    '''
    Equation 5.28, Voller
    '''
    phi_old = np.zeros(n)
    phi_new = np.zeros(n)
    RHS = 0.0

    for it in range(maxIt):
        if(np.mod(it, 100) == 0):
            print('Solver, iteration: ', it, '\t To finish: ', maxIt-it)
        for k in range(n):
            RHS = BB[k] + QB[k]
            for w in range(int(n_s[k])):
                RHS = RHS + a[k, w]*phi_new[int(S[k, w])-1]
            phi_new[k] = float(RHS/(ap[k]+BC[k]+QC[k]))

        if(max(abs(phi_new-phi_old)) <= tolerance):
            break
        phi_old = phi_new.copy()

    return phi_new


if __name__ == "__main__":
    arq = sys.argv[1]
    maxIt = 1500
    tolerance = 0.001

    [n, x, y, n_s, num_suport, S, num_seg, n_b, B] = read_file(arq)

    # ------------ Convert to polar coordinates ------------
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
    # --------- End of polar coordinates convertion --------

    print('Starting to set advection and diffusion coefficients')
    [a, ap, V] = set_coefficients(n, x, y, vx, vy, S, n_s)
    print('Finished to set advection and diffusion coefficients')

    print('Starting to set boundary conditions')
    [BC, BB] = set_boundary(B, n_b)
    print('Finished to set boundary conditions')

    print('Starting to set the source')
    [QC, QB] = set_source()
    print('Finished to set the source')

    print('\nStarting SOLVER')
    phi = solver(n, S, n_s, a, ap, BC, BB, QC, QB, maxIt, tolerance)
    print('Finished SOLVER\n')

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

    plt.plot(xout, ana, 'rd', label='Analitic')
    plt.plot(xout, phiout, 'bo', label='Numerical')
    plt.title('$\phi$ value as a function of radial position')
    plt.legend(loc='best')
    plt.show()
