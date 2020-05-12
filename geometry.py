import sys
import math


def distancia(x1, x2, y1, y2, z1=0, z2=0):
    """
    distancia(x1,x2,y1,y2,z1=0,z2=0) -> float

    Determina a distancia entre dois pontos

    Parametros
    ----------
    x1 -> coordenada x do primeiro ponto
    x2 -> coordenada x do segundo ponto
    y1 -> coordenada y do primeiro ponto
    y2 -> coordenada y do segundo ponto
    z1 -> coordenada z do primeiro ponto (opcional)
    z2 -> coordenada z do segundo ponto (opcional)

    Retorno
    -------
    float -> numero real
    """
    return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)


def convex_hull_index(points, points_index):
    """
    convex_hull_index(points,points_index) -> list

    Encontra os pontos que determinam a caixa minima de 
    um conjunto de pontos.
    Retorna uma lista de pontos no sentido anti-horario,
    juntamente com uma lista de seus respectivos indices.

    Parametros
    ----------
    points       -> lista de coordenadas dos pontos 
                    [[x0,y0],[x1,y1],...,[xn,yn]]
    points_index -> lista dos indices dos pontos 
                    [[idx0,idy1],[idx0,idy1],...,[idxn,idyn]]

    Retorno
    -------
    list_points       -> lista de coordenadas dos pontos 
                         ordenados no sentido anti-horario
    list_points_index -> lista de indices dos pontos 
                         ordenados no sentido anti-horario
    """
    points_original = points.copy()

    points = sorted(set(points))

    if len(points) <= 1:
        return points

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    lower_index = []
    for p in points:
        if(len(points) > 20):
            print(lower)
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for i in lower:
        lower_index.append(points_index[points_original.index(i)])

    upper = []
    upper_index = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    for i in upper:
        upper_index.append(points_index[points_original.index(i)])

    res_index = lower_index
    for i in upper_index:
        if(i not in res_index):
            res_index.append(i)

    return [lower[:-1] + upper[:-1], res_index]


def convex_hull(points):
    """
    convex_hull(points) -> list

    Encontra os pontos que determinam a caixa minima de 
    um conjunto de pontos.
    Retorna uma lista de pontos no sentido anti-horario

    Parametros
    -----------------------------------------
    points -> lista de coordenadas dos pontos 
              [[x0,y0],[x1,y1],...,[xn,yn]]

    Retorno
    ----------------------------------------------
    list_points -> lista de coordenadas dos pontos 
                   ordenados no sentido anti-horario
    """
    print('here')

    points = sorted(points)
    print(points)
    #points = sorted(set(points))

    if len(points) <= 1:
        return points

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in points:
        print(lower)
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(points):
        print(upper)
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return [lower[:-1] + upper[:-1]]


def poly_clockwise(x, y, poly_points):
    """
    poly_clockwise(x, y, poly_points) -> bool

    Determina se os pontos do dado poligono estao
    orenados no sentido horario

    Parametros
    ----------
    x             -> lista de coordendas x
    y             -> lista de coordendas y
    poly_points   -> lista de pontos
                     [p0,p1,p2,...,pn]

    Retorno
    -------
    bool -> Esta ou nao ordenado no sentido horario
    """

    # print(len(x))
    # print(len(y))
    # print(len(poly_points))
    soma = 0
    for k in range(len(poly_points)-1):
        # print(k,x[poly_points[k-1]],x[poly_points[k]],x[poly_points[k+1]],y[poly_points[k-1]],y[poly_points[k]],y[poly_points[k+1]])
        soma += (x[poly_points[k]] - x[poly_points[k-1]]) * \
            (y[poly_points[k+1]] - y[poly_points[k-1]])
        soma -= (x[poly_points[k+1]] - x[poly_points[k-1]]) * \
            (y[poly_points[k]] - y[poly_points[k-1]])
        # print(soma)

    if(soma == 0):
        print('\n Nao eh um poligono \n')
        sys.exit()

    return True if soma < 0 else False
