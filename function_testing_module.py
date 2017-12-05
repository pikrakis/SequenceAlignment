# Implementation: Aggelos Pikrakis, pikrakis@unipi.gr
# Educational material.
###################################################################################################
# this is where function testing takes place


import sequence_alignment
import other_functions
import matplotlib.pyplot as plt
import numpy as np

test = [-1, 2, 2, 2, 1, 1, 1, 1, -1, 2]
ref = [-1, 2, 1, -1, 2]
c, s, predi, predj = sequence_alignment.edit_distance(test,ref)
print('\ndtw cost = ' + str(c) + '\n')
bp = other_functions.extract_best_path(s, predi, predj)
print('best path = ' + str(bp))
bp = np.array(bp)

plt.plot(bp[:, 1], bp[:, 0], 'ro')
plt.plot(bp[:, 1], bp[:, 0], linewidth=2.0)
plt.ylabel('test sequence')
plt.xlabel('reference sequence')
plt.title('alignment path')
plt.show()
