from control import matlab
import matplotlib.pyplot as plt
import numpy as np



'''''''''''''''
以下は仮の値
'''''''''''''''
J=0.000537
mp=0.0425
Lp=0.09485
myup=0.0000685
Lm=0.192
Jm=0.00418
myum=0.114
Km=1.111
g=9.8
h11=J
h12=mp*Lp*Lm
h22=Jm+mp*Lm*Lm
DLT=h11*h22-h12*h12
T=0.005
''''''''''''''''''

A=np.array([[0, 1, 0, 0], 
            [h22*mp*g*Lp/DLT, -h22*myup/DLT, 0, h12*myum/DLT],
            [0, 0, 0, 1],
            [-h12*mp*g*Lp/DLT, h12*myup/DLT, 0, -h11*myum/DLT]])
b=np.array([[0],
           [-h12*Km/DLT],
           [0],
           [h11*Km/DLT]])
C=np.array([[1, 0, 0, 0],
            [0, 0, 1, 0]])
d=np.array([[0],
            [0]])


#可制御性の検証
Uc=matlab.ctrb(A, b)
print(Uc)
rank=np.linalg.matrix_rank(Uc)
print(rank)
Ev=np.linalg.eigvals(A)#極の表示、システムが安定するかがわかる
print('si=',Ev)


#可観測性の検証
Uo=matlab.obsv(A, C)
print(Uo)  #出力に三点リーダーが出てくるのはいっぱいあるのを省略している
rank=np.linalg.matrix_rank(Uo)
print(rank)

#離散時間系でもやってみましょう！
print("ここからは離散時間系での実行")

S=matlab.ss(A,b,C,d)
P=matlab.c2d(S,T,method='zoh')    #continuous to discrete

#科制御性の検証
Uc=matlab.ctrb(P.A,P.B)
rank=np.linalg.matrix_rank(Uc)
print(rank)
Ev=np.linalg.eigvals(P.A)#極の表示、システムが安定するかがわかる
print('zi=',Ev)

#可観測性の検証
Uo=matlab.obsv(P.A,P.C)
rank=np.linalg.matrix_rank(Uo)
print(rank)



#各種シミュレーション
x0=np.array([[0.1],[0.1],[0.1],[0.1]])
N=200
t=np.arange(0,N*T+T,T)


x=x0
y=x.T
for i in range(N):
    x=P.A@x
    y=np.r_[y,x.T]
#plt.plot(t,y)

Pol=np.array([0.9,0.9,0.9,0.9])
k=matlab.acker(P.A,P.B,Pol)
print('PP K=',k)

x=x0
y=x0.T
for i in range(N):
    x=(P.A-P.B@k)@x
    y=np.r_[y,x.T]
#plt.plot(t,y)

Wx=np.identity(4)
Wu=1
#H,Pol,k=matlab.dare(P.A,P.B,Wx,Wu)
p=P.A
q=P.B
y=np.array([[0,0,0,0]])
H=Wx
for i in range(2000):
    H=p.T@H@p-p.T@H@q@q.T@H@p/(Wu+q.T@H@q)+Wx
    k=q.T@H@p/(Wu+q.T@H@q)
    y=np.r_[y,k]
print('LQ k=', k)
#plt.plot(y)

x=x0
y=x0.T
for i in range(N):
    x=(P.A-P.B@k)@x
    y=np.r_[y,x.T]
plt.plot(t,y)


      