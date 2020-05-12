import sys
import numpy as np
import matplotlib.pyplot as plt

import geometry

#--------- Escrito em Python 3 ---------#

#--------- file.msh - nos ---------#
# --------- le os dados de um arquivo .msh bidimensional


def read_msh(arq):
    f = open(arq, 'r')

    line = ''
    while line.find('$MeshFormat') != 0:
        line = f.readline()

    mesh_format = f.readline()
    mesh_format = mesh_format[:-1]

    while line.find('$Nodes') != 0:
        line = f.readline()

    # Numero de nos
    line = f.readline()
    num_nos = int(line)

    # Coordenadas
    x = np.zeros(num_nos)
    y = np.zeros(num_nos)
    for i in range(num_nos):
        line = f.readline().split()

        x[i] = float(line[1])
        y[i] = float(line[2])

    line = ''
    while line.find('$Elements') != 0:
        line = f.readline()

    # Numero de elementos
    line = f.readline()
    num_elem_total = int(line)

    # Coordenadas
    eh_contorno = [False for i in range(num_nos)]
    num_contornos = 0
    elementos_ponto = []
    elementos_line = []
    elementos_tri = []
    num_elem_ponto = 0
    num_elem_line = 0
    num_elem_tri = 0
    for i in range(num_elem_total):
        line = f.readline().split()

        if(line[1] == '15'):
            num_elem_ponto += 1
            elementos_ponto.append(list(map(int, line[4:])))
            for k in elementos_ponto[-1]:
                if(not eh_contorno[k-1]):
                    eh_contorno[k-1] = True
                    num_contornos += 1

        elif(line[1] == '1'):
            num_elem_line += 1
            elementos_line.append(list(map(int, line[4:])))
            for k in elementos_line[-1]:
                if(not eh_contorno[k-1]):
                    eh_contorno[k-1] = True
                    num_contornos += 1

        elif(line[1] == '2'):
            num_elem_tri += 1
            elementos_tri.append(list(map(int, line[5:])))

    f.close()

    return [mesh_format,
            num_nos,
            x,
            y,
            num_contornos,
            eh_contorno,
            num_elem_total,
            num_elem_ponto,
            elementos_ponto,
            num_elem_line,
            elementos_line,
            num_elem_tri,
            elementos_tri]

# --------- Remove pontos que nao possuem conectividade
# --------- com nenhum elemento do dominio
# --------- Foram criados, possivelmente, somente como auxilio
# --------- Retorna True caso seja necessaria a criacao de um arquivo atualizado


def remove_desconexos(num_nos, x, y, num_contornos, eh_contorno, num_elem_total, num_elem_ponto, elementos_ponto, num_elem_line, elementos_line, num_elem_tri, elementos_tri):
    # Vetor de pontos possivelmente desconexos
    pontos_desconexos = [k[0] for k in elementos_ponto]
    num_desconexos = len(pontos_desconexos)

    # Determina se cada ponto possivelmente desconexos
    # eh realmente desconexo
    for k in elementos_tri:
        for w in pontos_desconexos:
            if(w in k):
                pontos_desconexos.remove(w)  # Remove ponto conexo
                if(pontos_desconexos == []):  # Caso nenhum ponto seja desconexo
                    # Nao ha a necessidade da criacao de um arquivo atualizado
                    return [False, num_nos, x, y, num_contornos, eh_contorno, num_elem_total, num_elem_ponto, elementos_ponto, elementos_line, elementos_tri]

    # Atualiza a quantidade de nos removendo-se os pontos desconexos
    num_nos -= len(pontos_desconexos)
    num_contornos -= len(pontos_desconexos)
    num_elem_total -= len(pontos_desconexos)
    num_elem_ponto -= len(pontos_desconexos)

    x = list(x)
    y = list(y)
    # Remove os pontos desconexos dos vetores de coordenadas
    for k in pontos_desconexos:
        x.pop(k-1)
        y.pop(k-1)
        eh_contorno.pop(k-1)

    x = np.array(x)
    y = np.array(y)

    # Atualiza os elementos linha apos remocao dos pontos desconexos
    for k in range(num_elem_line):
        for w in pontos_desconexos:
            for p in range(len(elementos_line[k][1:])):
                # print(k,w,p)
                if(elementos_line[k][p+1] >= w):
                    #print('here', elementos_line[k][p],elementos_line[k])
                    # input()
                    elementos_line[k][p+1] -= 1

    # Atualiza os elementos triangulo apos remocao dos pontos desconexos
    for k in range(num_elem_tri):
        for w in pontos_desconexos:
            for p in range(len(elementos_tri[k])):
                if(elementos_tri[k][p] >= w):
                    elementos_tri[k][p] -= 1

    # Ha a necessidade da criacao de um arquivo atualizado
    return [True, num_nos, x, y, num_contornos, eh_contorno, num_elem_total, num_elem_ponto, elementos_ponto, elementos_line, elementos_tri]

# --------- Transforma a relacao dos elementos
# --------- em matriz de suporte


def elemento_para_suporte(num_elem_tri, elementos_line, elementos_tri, x, y, num_nos, eh_contorno):
    S = [set() for i in range(1, num_nos+1)]

    for i in elementos_tri:
        for j in i:
            S[int(j-1)] = S[int(j-1)].union(set(list(i)))

    for i in range(len(S)):
        S[i] = S[i] - {i+1}

    S_temp = S.copy()

    [points, S] = organiza_suportes(S, x, y, num_nos, eh_contorno)
    B = matriz_contorno(elementos_line, x, y)

    return [S, S_temp, B]

#--------- Ordena os suportes dos nos no sentido anti-horario ---------#


def organiza_suportes(S, x, y, num_nos, eh_contorno):
    S_lista = []
    S_points = []
    X = list(x)
    Y = list(y)
    for k, i in enumerate(S):
        points = [(x[int(j)-1], y[int(j)-1]) for j in i]
        points_index = list(i)

        [points, points_index] = geometry.convex_hull_index(points, points_index)

        if(eh_contorno[k]):
            points_index = ordena_suportes_contorno(points_index)
            points_index.append(0)
        else:
            points_index.append(points_index[0])

        S_lista.append(points_index)
        S_points.append(points)

    return [S_points, S_lista]

# --------- Ordena os suportes de um contorno dos nos no sentido
# --------- anti-horario de acordo com os criterios do Voller.
# --------- Precisa-se comecar a contar de um ponto no contorno e
# --------- o caminho precisa estar dentro do dominio


def ordena_suportes_contorno(points_index):
    if(len(points_index) == 2):
        return points_index
    elif(eh_contorno[int(points_index[0]-1)] and eh_contorno[int(points_index[-1]-1)]):
        return points_index
    else:
        for i in range(len(points_index)+1):
            points_index.append(points_index.pop(0))
            if(eh_contorno[int(points_index[0]-1)] and eh_contorno[int(points_index[-1]-1)]):
                return points_index
        print('Ha algo errado com o arquivo')
        print('Parou em: ', points_index)
        sys.exit(1)

#--------- Cria matriz de contorno ---------#


def matriz_contorno(elementos_line, x, y):
    B = []                # Matriz de contornos
    B_aux = []            # Matriz de contornos auxiliar
    vet_num_entidade = []  # Vetor que armazena o numero da entidade de cada contorno

    # Inicializa os vetores
    for k in elementos_line:
        if(k[0] not in vet_num_entidade):
            vet_num_entidade.append(k[0])
            B.append([])
            B_aux.append([])

    # Preenche B_aux com os pares de nos de cada elementos linha do contorno
    for k in elementos_line:
        B_aux[vet_num_entidade.index(k[0])].append(k[1:])

    # Descobre quais sao os pontos extremos do primeiro segmento
    pontos_extremos = []
    eh_extremo1 = eh_extremo2 = True
    for k in B_aux[0]:
        for w in B_aux[0]:
            if(k is not w):
                if(k[0] == w[0] or k[0] == w[1]):
                    eh_extremo1 = False
                if(k[1] == w[0] or k[1] == w[1]):
                    eh_extremo2 = False
                if(not eh_extremo1 and not eh_extremo2):
                    break
            # endif
        # endfor
        if(eh_extremo1):
            pontos_extremos.append(k[0])
        elif(eh_extremo2):
            pontos_extremos.append(k[1])
        eh_extremo1 = eh_extremo2 = True
        # endfor
    # endfor

    # ---------------------------------
    # Chutar um ponto extremo inicial,
    # ordenar tudo e depois conferir
    # se estah anti-horario
    # ---------------------------------

    ponto_aux = pontos_extremos[0]
    indice_B = 0
    indice_B_aux = 0
    while B_aux != []:
        encontrou_ponto = False
        segmento = B_aux[indice_B_aux]
        while segmento != []:
            for k in segmento:
                if(ponto_aux == k[0]):
                    encontrou_ponto = True
                    B[indice_B].append(k[0])
                    ponto_aux = k[1]
                    segmento.remove(k)
                elif(ponto_aux == k[1]):
                    encontrou_ponto = True
                    B[indice_B].append(k[1])
                    ponto_aux = k[0]
                    segmento.remove(k)
            # endfor
            if(not encontrou_ponto):
                if(indice_B_aux + 1 >= len(B_aux)):
                    indice_B_aux = 0
                    segmento = B_aux[0]
                else:
                    indice_B_aux += 1
                    segmento = B_aux[indice_B_aux]
        # endwhile
        B[indice_B].append(ponto_aux)
        B_aux.pop(indice_B_aux)
        indice_B_aux = 0
        indice_B += 1
    # endfor

    # Transforma a matriz de contorno em um vetor unidimensioanal
    pontos_contorno = [w for k in B for w in k[1:]]
    # print(pontos_contorno)

    # Inverte o contorno caso esteja no sentido horario
    if(geometry.poly_clockwise(x, y, pontos_contorno)):
        for k in B:
            k.reverse()
        B.reverse()

    return B

#--------- Escreve o arquivo para ser lido pelo cvfem ---------#


def escreve_arquivo_cvfem(arq, n, x, y, n_s, S, n_b, B):
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


#--------- Principal ---------#
arq = sys.argv[1]
[mesh_format,
 num_nos,
 x,
 y,
 num_contornos,
 eh_contorno,
 num_elem_total,
 num_elem_ponto,
 elementos_ponto,
 num_elem_line,
 elementos_line,
 num_elem_tri,
 elementos_tri
 ] = read_msh(arq)

[precisa,
 num_nos,
 x,
 y,
 num_contornos,
 eh_contorno,
 num_elem_total,
 num_elem_ponto,
 elementos_ponto,
 elementos_line,
 elementos_tri
 ] = remove_desconexos(num_nos, x, y, num_contornos, eh_contorno, num_elem_total, num_elem_ponto, elementos_ponto, num_elem_line, elementos_line, num_elem_tri, elementos_tri)

[S,
 S_temp,
 B
 ] = elemento_para_suporte(num_elem_tri, elementos_line, elementos_tri, x, y, num_nos, eh_contorno)

n_s = [len(i)-1 for i in S]
n_b = [len(i) for i in B]

aux = []
for i in range(len(S)):
    if(len(S[i]) < max(n_s)+1):
        aux = [0 for j in range(max(n_s)+1-len(S[i]))]
        S[i].extend(aux)
        aux = []

escreve_arquivo_cvfem('saida' + arq.split('.')
                      [0].capitalize() + '.dat', num_nos, x, y, n_s, S, n_b, B)

