from control import matlab		#制御モジュール
import numpy as np			#数値計算モジュール
import matplotlib.pyplot as plt

J=5.37e-4				#倒立振子パラメータ
mp=4.25e-2
Lp=9.485e-2
myup=6.85e-5
Lm=0.192
Jm=4.18e-3
myum=0.114
Km=1.111
g=9.8

h11=J					#システム行列パラメータ
h12=mp*Lp*Lm
h22=Jm+mp*Lm*Lm
DLT=h11*h22-h12*h12
T=5e-3					#サンプリング周期

A=np.array([[0,1,0,0],			#システム行列設定
            [h22*mp*g*Lp/DLT,-h22*myup/DLT,0,h12*myum/DLT],
            [0,0,0,1],
            [-h12*mp*g*Lp/DLT,h12*myup/DLT,0,-h11*myum/DLT]])
b=np.array([[0],[-h12*Km/DLT],[0],[h11*Km/DLT]])
C=np.array([[1,0,0,0],
            [0,0,1,0]])
d=np.array([[0],[0]])

Uc=matlab.ctrb(A,b)			#連続系可制御性行列作成
print(Uc)				#その表示
rank=np.linalg.matrix_rank(Uc)		#そのランク調査
print(rank)				#ランクが4なら可制御
Ev=np.linalg.eigvals(A)			#固有値=特性根計算
print(Ev)				#左半面なら安定

S=matlab.ss(A,b,C,d)			#システム行列の設定
P=matlab.c2d(S,T,method='zoh')		#離散系へ変換

Uc=matlab.ctrb(P.A,P.B)			#離散系での可制御性行列
print(Uc)
rank=np.linalg.matrix_rank(Uc)		#ランク調査
print(rank)				#ランクが4なら可制御
Ev=np.linalg.eigvals(P.A)		#固有値=特性根計算
print(Ev)				#単位円内なら安定
Uo=matlab.obsv(P.A,P.C)			#可観測性行列作成
print(Uo)
rank=np.linalg.matrix_rank(Uo)		#ランク調査
print(rank)				#ランクが4なら可観測

x0=np.array([[0.1],[0.1],[0.1],[.1]])
N=200
t=np.arange(0,N*T+T,T)

x=x0
y=x0.T
for i in range(N):
    x=P.A@x
    y=np.r_[y,x.T]
#plt.plot(t,y)

Pol=np.array([0.9, 0.9, 0.9, 0.9])
k=matlab.acker(P.A,P.B,Pol)
print('PP k=',k)

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
print('LQ k=',k)
plt.plot(y)

x=x0
y=x0.T
for i in range(N):
    x=(P.A-P.B@k)@x
    y=np.r_[y,x.T]
#plt.plot(t,y)
