# Implementation: Aggelos Pikrakis, pikrakis@unipi.gr
# Educational material.
# Distributed under the GNU General Public License v2.0
###################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import other_functions
###################################################################################################

def edit_distance(y, x):
    # computes the standard edit (Levenstein) distance between two symbol strings
    lx = len(x)
    ly = len(y)
    s = float('inf') * np.ones((ly, lx))
    predi = np.ones((ly, lx)) * (-1)
    predj = np.ones((ly, lx)) * (-1)
    s[0, 0] = int(y[0] != x[0])
    # first column
    for i in range(1,ly):
        s[i, 0] = s[i-1, 0]+ 1
        predi[i, 0] = i - 1
        predj[i, 0] = 0

    # first row
    for j in range(1, lx):
        s[0, j] = s[0, j-1] + 1
        predi[0, j] = 0
        predj[0, j] = j - 1


    # all other cost grid elements
    for i in range(1,ly):
        for j in range(1,lx):
            q= [s[i-1, j-1]+int(y[i] != x[j]), s[i-1, j]+1, s[i, j-1] + 1]
            s[i, j] = np.amin(q)
            if np.argmin(q)==0:
                predi[i, j] = i - 1
                predj[i, j] = j - 1
            elif np.argmin(q)==1:
                predi[i, j] = i - 1
                predj[i, j] = j
            else:
                predi[i, j] = i
                predj[i, j] = j - 1

    return s[ly-1, lx-1], s, predi, predj
###################################################################################################


def dtw_distance_sc(y, x, gp_constraints):
    # computes the dynamic time warping cost between two sequences of numbers using the Sakoe-Chiba
    # local path constraints
    lx = len(x)
    ly = len(y)
    s = float('inf') * np.ones((ly, lx))
    predi = np.ones((ly, lx)) * (-1)
    predj = np.ones((ly, lx)) * (-1)
    s[0, 0] = abs(y[0] - x[0])

    # first column
    for i in range(1,ly):
        if abs(i-0)<=gp_constraints:
            s[i, 0] = s[i-1, 0] + abs(y[i] - x[0])
            predi[i, 0] = i-1
            predj[i, 0] = 0


    # first row
    for j in range(1, lx):
        if abs(0 - j) <= gp_constraints:
            s[0, j] = s[0, j-1] + abs(y[0] - x[j])
            predi[0, j] = 0
            predj[0, j] = j-1

    # all other cost grid elements
    for i in range(1,ly):
        for j in range(1,lx):
            if abs(i - j) <= gp_constraints:
                q = [s[i-1, j-1], s[i-1, j], s[i, j-1]]
                s[i, j] = np.amin(q) + abs(y[i] - x[j])
                if np.argmin(q)==0:
                    predi[i, j] = i - 1
                    predj[i, j] = j - 1
                elif np.argmin(q)==1:
                    predi[i, j] = i - 1
                    predj[i, j] = j
                else:
                    predi[i, j] = i
                    predj[i, j] = j - 1
    return s[ly-1, lx-1], s, predi, predj
###################################################################################################


def dtw_distance_ita(y, x, gp_constraints):
    # computes the dynamic time warping cost between two sequences of numbers using the Itakura
    # local path constraints
    lx = len(x)
    ly = len(y)
    s = float('inf')*np.ones((ly,lx))
    predi = np.ones((ly,lx))*(-1)
    predj = np.ones((ly, lx))*(-1)
    s[0, 0] = abs(y[0] - x[0])

    # first row
    if abs(0-1)<=gp_constraints:
        s[0, 1] = s[0, 0] + abs(y[0]-x[1])
        predi[0, 1] = 0
        predj[0, 1] = 0

    # all other cost grid elements
    for i in range(1,ly):
        for j in range(1,lx):
            if abs(i-j) <= gp_constraints:
                offset = 0
                if i > 1:
                    if j > 1 and predi[i, j-1] == i:
                        q = [s[i - 1, j - 1], s[i - 2, j - 1]]
                        offset = 1
                    else:
                        q = [s[i, j-1], s[i-1, j-1], s[i-2, j-1]]
                else:
                    if j > 1 and predi[i, j-1] == i:
                        q = [s[i - 1, j - 1]]
                        offset = 1
                    else:
                        q = [s[i, j - 1], s[i - 1, j - 1]]
                s[i, j] = np.amin(q) + abs(y[i] - x[j])
                predi[i, j] = i - offset - np.argmin(q)
                predj[i, j] = j - 1

    return s[ly-1, lx-1], s, predi, predj
###################################################################################################


def lcs(y, x):
    # computes the length of the longest common subsequence of y and x
    lx = len(x)
    ly = len(y)
    frwrdi = np.ones((ly, lx)) * (-1)
    frwrdj = np.ones((ly, lx)) * (-1)
    s = np.zeros((ly, lx))

    # last column
    for i in range(0, ly):
        s[i, lx-1] = int(y[i] == x[lx-1])

    # last row
    for j in range(0, lx):
        s[ly-1, j] = int(y[ly-1] == x[j])

    for i in range(ly-2, -1, -1):
        for j in range(lx-2, -1, -1):
            m1 = np.amax(s[i + 1,j + 1:lx])
            m2 = np.amax(s[i + 1:ly, j + 1])
            if np.amax([m1, m2]) == m1:
                w = np.argmax(s[i + 1, j + 1:lx]) + 1
                frwrdi[i, j] = i + 1
                frwrdj[i, j] = j + w
            else:
                w = np.argmax(s[i + 1:ly, j + 1]) + 1
                frwrdi[i, j] = i + w
                frwrdj[i, j] = j+1
            s[i, j] = int(y[i] == x[j]) + np.amax([m1,m2])
    return np.amax(s), s, frwrdi, frwrdj

###################################################################################################

def backtracking(s, predi, predj):

    # performs backtracking on a dissimilarity grid, starting from the node that corresponds to the last element of
    # each sequence

    v = s.shape
    xc = int(v[0]) - 1
    yc = int(v[1]) - 1
    bp=[]

    while (xc > -1) and (yc > -1):
        bp.append([xc, yc])
        tmpxc = int(predi[xc, yc])
        tmpyc = int(predj[xc, yc])
        xc = tmpxc
        yc = tmpyc
    return bp[::-1]
###################################################################################################


def forwardtracking(s, frwrdi, frwrdj):
    # performs forward tracking on a similarity grid, starting from the node that has accumulated the maximun similarity
    # value. This function can be used in cojunction with the lcs implementation

    xc, yc = other_functions.aux_fun(s, np.max)
    bp = []

    while (xc > -1) and (yc > -1):
        bp.append([xc, yc])
        tmpxc = int(frwrdi[xc, yc])
        tmpyc = int(frwrdj[xc, yc])
        xc = tmpxc
        yc = tmpyc
    return bp


def print_alignment(y, x, bp):
    # prints the alignment of y vs x, inserting gaps when needed
    L = len(bp)
    for k in range(L-1):
        print('%3d <----> %3d' % (y[bp[k][0]],x[bp[k][1]]))
        if bp[k+1][0] > bp[k][0]+1:
            for m in range(bp[k][0]+1, bp[k+1][0]):
                print('%3d <----> %3s'% (y[m], '-'))
        if bp[k+1][1] > bp[k][1]+1:
            for m in range(bp[k][1]+1, bp[k+1][1]):
                print('%3s <----> %3d' % ('-', x[m]))
    print('%3d <----> %3d' % (y[bp[-1][0]], x[bp[-1][1]]))

