from control import matlab
import matplotlib.pyplot as plt
import numpy as np

A=np.array([[0, 1], 
            [-1,  -1.8]])
b=np.array([[0],
           [1]])
Uc=matlab.ctrb(A, b)
print(Uc)
rank=np.linalg.matrix_rank(Uc)
print(rank)