from control import matlab
import numpy as np
import matplotlib.pyplot as plt

A=np.array([[0, 1], [-1, -1.8]])
b=np.array([[0], [1]])
c=np.array([1, 0])
d=0
T=0.3

S=matlab.ss(A, b, c, d)
#print(S)
Gs=matlab.ss2tf(S)
#print(Gs)
P=matlab.c2d(S, T, method='zoh')
#print(P)
Gz=matlab.ss2tf(P) #PはSを離散時間系に変換したやつ
#print(Gz) #p.13の問題6.の答え
#matlab.bode(S ,matlab.logspace(-2, 2)) #ボード線図を描こう！
yt, tt = matlab.step(S, np.arange(0, 10, 0.01)) #まずは連続時間系
#plt.plot(tt, yt) #連続時間系での出力をグラフにする
yk, tk = matlab.step(P, np.arange(0, 10, T)) #次は離散時間系
plt.plot(tk, yk, marker='o') #離散時間系のプロット、プロット点も描画する
