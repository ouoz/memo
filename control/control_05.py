from control import matlab
import numpy as np
import matplotlib.pyplot as plt

P=np.array([[0,1],
            [-1,0]])
q=np.array([[1],
            [0]])
Ev=np.linalg.eigvals(P)
print('Ev=',Ev)

x0=np.array([[1],
             [0]])
N=10
x=x0
y=x.T
for i in range(N):
    x=P@x
    y=np.r_[y,x.T]
#plt.plot(y)
#print(y)

Pol=np.array([0.5, 0.1])
k=matlab.acker(P,q,Pol)
#print('PP k=', k)
Ev=np.linalg.eigvals(P-q@k)
#print('Ev=', Ev)

x=x0
y=x.T
for i in range(N):
    x=(P-q@k)@x
    y=np.r_[y,x.T]
#plt.plot(y)

Wx=np.identity(2)
Wu=1
H,Pol,k=matlab.dare(P,q,Wx,Wu)
print('LQ k=',k)

x=x0
y=x.T
for i in range(N):
    x=(P-q@k)@x
    y=np.r_[y,x.T]
plt.plot(y)