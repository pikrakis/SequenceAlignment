# Implementation: Aggelos Pikrakis, pikrakis@unipi.gr
# Educational material.
# Distributed under the GNU General Public License v2.0
###################################################################################################
import random
import numpy as np


def cr_random_list(minrange, maxrange, listlength):
    # creates a list of length listlength of random integers in [minrange,maxrange]
    a=[]
    for i in range(0,listlength):
        a.append(random.randint(minrange,maxrange))
    return a


def extract_best_path(s, predi, predj):
    # c = list(np.where(s == np.min(s)))
    # xc = c[0][0]
    # yc = c[1][0]

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


def aux_fun(s, f):
    c = list(np.where(s == f(s)))
    return c



