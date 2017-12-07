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


def aux_fun(s, f):
    c = list(np.where(s == f(s)))
    cx = c[0][-1]
    cy = c[1][-1]
    return cx, cy



