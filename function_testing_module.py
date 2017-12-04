# this is where function testing takes place

import sequence_alignment
import other_functions

test = [-1, 2, 1, -1, -1, 2]
ref = [-1, 2, 1, -1, 2]
c, s, predi, predj = sequence_alignment.dtw_distance_ita(test,ref)
print('\ndtw cost = ' + str(c) + '\n')
bp = other_functions.extract_best_path(s, predi, predj)
print('best path = ' + str(bp))
