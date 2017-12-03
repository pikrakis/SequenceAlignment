import random

def cr_random_list(minrange, maxrange, listlength):
    # creates a list of length listlength of random integers in [minrange,maxrange]
    a=[]
    for i in range(0,listlength):
        a.append(random.randint(minrange,maxrange))
    return a

