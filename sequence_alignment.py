# Implementation: Aggelos Pikrakis, pikrakis@unipi.gr
# Educational material.
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


def dtw_distance_sc(y, x):
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
        s[i, 0] = s[i-1, 0] + abs(y[i] - x[0])
        predi[i, 0] = i-1
        predj[i, 0] = 0


    # first row
    for j in range(1, lx):
        s[0, j] = s[0, j-1] + abs(y[0] - x[j])
        predi[0, j] = 0
        predj[0, j] = j-1

    # all other cost grid elements
    for i in range(1,ly):
        for j in range(1,lx):
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


def dtw_distance_ita(y, x):
    # computes the dynamic time warping cost between two sequences of numbers using the Itakura
    # local path constraints
    lx = len(x)
    ly = len(y)
    s = float('inf')*np.ones((ly,lx))
    predi = np.ones((ly,lx))*(-1)
    predj = np.ones((ly, lx))*(-1)
    s[0, 0] = abs(y[0] - x[0])

    # first row
    s[0, 1] = s[0, 0] + abs(y[0]-x[1])
    predi[0, 1] = 0
    predj[0, 1] = 0

    # all other cost grid elements
    for i in range(1,ly):
        for j in range(1,lx):
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
    s = np.zeros((ly, lx))
    # last column
    for i in range(0, ly):
        s[i, lx-1] = int(y[i] == x[lx-1])

    # last row
    for j in range(0, lx):
        s[ly-1, j] = int(y[ly-1] == x[j])

    for i in range(ly-2,-1,-1):
        for j in range(lx-2,-1,-1):
            m1 = np.amax(s[i+1,j+1:lx])
            m2 = np.amax(s[i+1:ly,j+1])
            s[i, j] = int(y[i] == x[j]) + np.amax([m1,m2])
    return np.amax(s), s

###################################################################################################

