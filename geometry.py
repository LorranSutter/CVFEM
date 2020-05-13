import sys
import math


def distance(x1, x2, y1, y2, z1=0, z2=0):
    """
    distance(x1,x2,y1,y2,z1=0,z2=0) -> float

    Set the distance between two points

    Params
    ----------
    x1 -> coordinate x of the first point
    x2 -> coordinate x of the second point
    y1 -> coordinate y of the first point
    y2 -> coordinate y the second point
    z1 -> coordinate z of the first point (optional)
    z2 -> coordinate z the second point (optional)

    Return
    -------
    float -> real number
    """
    return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)


def convex_hull_index(points, points_index):
    """
    convex_hull_index(points,points_index) -> list

    Find the points that determine the minimum box of a set of points.
    Returns a list of points counterclockwise,
    along with a list of their respective indices.

    Params
    ----------
    points       -> coordinate list of points
                    [[x0,y0],[x1,y1],...,[xn,yn]]
    points_index -> list of point indexes
                    [[idx0,idy1],[idx0,idy1],...,[idxn,idyn]]

    Return
    -------
    list_points       -> coordinate list of points
                         ordered counterclockwise
    list_points_index -> list of point indexes
                         ordered counterclockwise
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

    Find the points that determine the minimum box of a set of points.
    Returns a list of points counterclockwise

    Params
    -----------------------------------------
    points -> coordinate list of points
              [[x0,y0],[x1,y1],...,[xn,yn]]

    Return
    ----------------------------------------------
    list_points -> coordinate list of points
                   ordered counterclockwise
    """
    print('here')

    points = sorted(points)
    print(points)

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

    Determines whether the points of the given polygon are clockwise

    Params
    ----------
    x             -> list of coordinates x
    y             -> list of coordinates y
    poly_points   -> list of points
                     [p0,p1,p2,...,pn]

    Return
    -------
    bool -> Is it ordered or not clockwise
    """

    sum_output = 0
    for k in range(len(poly_points)-1):
        sum_output += (x[poly_points[k]] - x[poly_points[k-1]]) * \
            (y[poly_points[k+1]] - y[poly_points[k-1]])
        sum_output -= (x[poly_points[k+1]] - x[poly_points[k-1]]) * \
            (y[poly_points[k]] - y[poly_points[k-1]])

    if(sum_output == 0):
        print('\n It is not a polygon \n')
        sys.exit()

    return True if sum_output < 0 else False
