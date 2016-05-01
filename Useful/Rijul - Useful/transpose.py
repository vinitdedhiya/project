''' This function returns the transpose of the 1D row matrix entered.'''

import numpy as np
def transpose(A):
    n = len(A)
    A_T = np.ones([n,1]);
    for i in range(n):
        A_T[[i],[0]] = A[i];

    return A_T