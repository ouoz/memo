from control import matlab
import matplotlib.pyplot as plt
import numpy as np

Gs = matlab.tf(0.5, [1, 0.5])
print(Gs)
T = 0.1
Gz = matlab.c2d(Gs, T, method = 'zoh') #'zoh'0次ホールド
print(Gz)
