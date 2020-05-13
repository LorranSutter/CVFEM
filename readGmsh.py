import sys
import numpy as np
import matplotlib.pyplot as plt

import geometry


def read_msh(arq):
    '''
    file.msh - nodes
    read two-dimensional .msg file data
    '''
    f = open(arq, 'r')

    line = ''
    while line.find('$MeshFormat') != 0:
        line = f.readline()

    mesh_format = f.readline()
    mesh_format = mesh_format[:-1]

    while line.find('$Nodes') != 0:
        line = f.readline()

    # Number of nodes
    line = f.readline()
    num_nodes = int(line)

    # Coordinates
    x = np.zeros(num_nodes)
    y = np.zeros(num_nodes)
    for i in range(num_nodes):
        line = f.readline().split()

        x[i] = float(line[1])
        y[i] = float(line[2])

    line = ''
    while line.find('$Elements') != 0:
        line = f.readline()

    # Number of elements
    line = f.readline()
    num_elem_total = int(line)

    # Coordinates
    is_contour = [False for i in range(num_nodes)]
    num_contours = 0
    point_elements = []
    line_elements = []
    tri_elements = []
    num_elem_point = 0
    num_elem_line = 0
    num_elem_tri = 0
    for i in range(num_elem_total):
        line = f.readline().split()

        if(line[1] == '15'):
            num_elem_point += 1
            point_elements.append(list(map(int, line[4:])))
            for k in point_elements[-1]:
                if(not is_contour[k-1]):
                    is_contour[k-1] = True
                    num_contours += 1

        elif(line[1] == '1'):
            num_elem_line += 1
            line_elements.append(list(map(int, line[4:])))
            for k in line_elements[-1]:
                if(not is_contour[k-1]):
                    is_contour[k-1] = True
                    num_contours += 1

        elif(line[1] == '2'):
            num_elem_tri += 1
            tri_elements.append(list(map(int, line[5:])))

    f.close()

    return [mesh_format,
            num_nodes,
            x,
            y,
            num_contours,
            is_contour,
            num_elem_total,
            num_elem_point,
            point_elements,
            num_elem_line,
            line_elements,
            num_elem_tri,
            tri_elements]


def remove_disconnected(num_nodes, x, y, num_contours, is_contour, num_elem_total, num_elem_point, point_elements, num_elem_line, line_elements, num_elem_tri, tri_elements):
    '''
    Remove all points that does not have connectivity
    with an element from the domain
    They were created, possibly, only as an aux
    Return True in case is necessary to create an updated file
    '''

    # Array of possibly disconnected points
    disconnected_points = [k[0] for k in point_elements]

    # Determines if each possibly disconnected point is really disconnected
    for k in tri_elements:
        for w in disconnected_points:
            if(w in k):
                disconnected_points.remove(w)  # Remove connected point
                if(disconnected_points == []):  # If no point is connected
                    # There is no need to create an updated file
                    return [False, num_nodes, x, y, num_contours, is_contour, num_elem_total, num_elem_point, point_elements, line_elements, tri_elements]

    # Update the quantity of nodes removing them from disconnected points
    num_nodes -= len(disconnected_points)
    num_contours -= len(disconnected_points)
    num_elem_total -= len(disconnected_points)
    num_elem_point -= len(disconnected_points)

    x = list(x)
    y = list(y)

    # Remove disconnected points from coordinates array
    for k in disconnected_points:
        x.pop(k-1)
        y.pop(k-1)
        is_contour.pop(k-1)

    x = np.array(x)
    y = np.array(y)

    # Update line elements after remove disconnected points
    for k in range(num_elem_line):
        for w in disconnected_points:
            for p in range(len(line_elements[k][1:])):
                if(line_elements[k][p+1] >= w):
                    line_elements[k][p+1] -= 1

    # Update triangle elements after remove disconnected points
    for k in range(num_elem_tri):
        for w in disconnected_points:
            for p in range(len(tri_elements[k])):
                if(tri_elements[k][p] >= w):
                    tri_elements[k][p] -= 1

    # There is a need to create an updated file
    return [num_nodes, x, y, num_contours, is_contour, num_elem_total, num_elem_point, point_elements, line_elements, tri_elements]


def element_to_support(num_elem_tri, line_elements, tri_elements, x, y, num_nodes, is_contour):
    '''
    Transform the list of elements into a support matrix
    '''
    S = [set() for i in range(1, num_nodes+1)]

    for i in tri_elements:
        for j in i:
            S[int(j-1)] = S[int(j-1)].union(set(list(i)))

    for i in range(len(S)):
        S[i] = S[i] - {i+1}

    S_temp = S.copy()

    S = organizes_supports(S, x, y, num_nodes, is_contour)
    B = contour_matrix(line_elements, x, y)

    return [S, S_temp, B]


def organizes_supports(S, x, y, num_nodes, is_contour):
    '''
    Sorts the supports of the nodes counterclockwise
    '''
    S_list = []
    for k, i in enumerate(S):
        points = [(x[int(j)-1], y[int(j)-1]) for j in i]
        points_index = list(i)

        [points, points_index] = geometry.convex_hull_index(
            points, points_index)

        if(is_contour[k]):
            points_index = sort_contour_supports(points_index)
            points_index.append(0)
        else:
            points_index.append(points_index[0])

        S_list.append(points_index)

    return S_list


def sort_contour_supports(points_index):
    '''
    Sorts the supports of the nodes counterclockwise according to Voller criteria.
    You need to start counting from a point on the contour and the path must be within the domain
    '''
    if(len(points_index) == 2):
        return points_index
    elif(is_contour[int(points_index[0]-1)] and is_contour[int(points_index[-1]-1)]):
        return points_index
    else:
        for _ in range(len(points_index)+1):
            points_index.append(points_index.pop(0))
            if(is_contour[int(points_index[0]-1)] and is_contour[int(points_index[-1]-1)]):
                return points_index
        print('There is something wrong with the file')
        print('Stopped in: ', points_index)
        sys.exit(1)


def contour_matrix(line_elements, x, y):
    '''
    Creates contour matrix
    '''
    B = []                # Contour matrix
    B_aux = []            # Aux contour matrix
    num_entity_arr = []   # Array that stores the entity number of each contour

    # Init arrays
    for k in line_elements:
        if(k[0] not in num_entity_arr):
            num_entity_arr.append(k[0])
            B.append([])
            B_aux.append([])

    # Fill B_aux with the pairs of nodes of each contour line elements
    for k in line_elements:
        B_aux[num_entity_arr.index(k[0])].append(k[1:])

    # Find out what the extreme points of the first segment are
    extreme_points = []
    is_extreme1 = is_extreme2 = True
    for k in B_aux[0]:
        for w in B_aux[0]:
            if(k is not w):
                if(k[0] == w[0] or k[0] == w[1]):
                    is_extreme1 = False
                if(k[1] == w[0] or k[1] == w[1]):
                    is_extreme2 = False
                if(not is_extreme1 and not is_extreme2):
                    break
        if(is_extreme1):
            extreme_points.append(k[0])
        elif(is_extreme2):
            extreme_points.append(k[1])
        is_extreme1 = is_extreme2 = True

    # ---------------------------------
    # Choose arbitrarily an initial extreme point,
    # sort all and then check
    # if is counterclockwise
    # ---------------------------------

    point_aux = extreme_points[0]
    index_B = 0
    index_B_aux = 0
    while B_aux != []:
        found_point = False
        segment = B_aux[index_B_aux]
        while segment != []:
            for k in segment:
                if(point_aux == k[0]):
                    found_point = True
                    B[index_B].append(k[0])
                    point_aux = k[1]
                    segment.remove(k)
                elif(point_aux == k[1]):
                    found_point = True
                    B[index_B].append(k[1])
                    point_aux = k[0]
                    segment.remove(k)
            if(not found_point):
                if(index_B_aux + 1 >= len(B_aux)):
                    index_B_aux = 0
                    segment = B_aux[0]
                else:
                    index_B_aux += 1
                    segment = B_aux[index_B_aux]
        B[index_B].append(point_aux)
        B_aux.pop(index_B_aux)
        index_B_aux = 0
        index_B += 1

    # Turns the contour matrix into an one-dimensional array
    contour_points = [w for k in B for w in k[1:]]

    # Inverts the contour if clockwise
    if(geometry.poly_clockwise(x, y, contour_points)):
        for k in B:
            k.reverse()
        B.reverse()

    return B


def write_cvfem_file(arq, n, x, y, n_s, S, n_b, B):
    '''
    Write file to be read by cvfem
    '''
    f = open(arq, 'w')

    f.write(str(n) + '\n')
    for i, j in zip(x, y):
        f.write(str(i) + ' ' + str(j) + '\n')
    for i in n_s:
        f.write(str(i) + '\n')
    f.write(str(max(n_s)) + '\n')
    for i in S:
        for j in i:
            f.write(str(j) + ' ')
        f.write('\n')
    f.write(str(len(n_b)) + '\n')
    for i in n_b:
        f.write(str(i) + '\n')
    for i in B:
        for j in i:
            f.write(str(j) + ' ')
        f.write('\n')
    f.close()


#--------- Main ---------#
arq = sys.argv[1]
[mesh_format,
 num_nodes,
 x,
 y,
 num_contours,
 is_contour,
 num_elem_total,
 num_elem_point,
 point_elements,
 num_elem_line,
 line_elements,
 num_elem_tri,
 tri_elements
 ] = read_msh(arq)

[num_nodes,
 x,
 y,
 num_contours,
 is_contour,
 num_elem_total,
 num_elem_point,
 point_elements,
 line_elements,
 tri_elements
 ] = remove_disconnected(num_nodes, x, y, num_contours, is_contour, num_elem_total, num_elem_point, point_elements, num_elem_line, line_elements, num_elem_tri, tri_elements)

[S,
 S_temp,
 B
 ] = element_to_support(num_elem_tri, line_elements, tri_elements, x, y, num_nodes, is_contour)

n_s = [len(i)-1 for i in S]
n_b = [len(i) for i in B]

aux = []
for i in range(len(S)):
    if(len(S[i]) < max(n_s)+1):
        aux = [0 for j in range(max(n_s)+1-len(S[i]))]
        S[i].extend(aux)
        aux = []

write_cvfem_file('output' + arq.split('.')
                 [0].capitalize() + '.dat', num_nodes, x, y, n_s, S, n_b, B)
