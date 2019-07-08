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

#可観測性の検証
Uo=matlab.obsv(P.A,P.C)
rank=np.linalg.matrix_rank(Uo)
print(rank)


'''
中島「磯野～！安定判別しようぜ～！」
カツオ「ええ～！？そんな面倒くさいことしたくないよ～」
中島「そこでコンピュータの力を借りるのさ、固有値問題なんてすぐ解けるよ」
'''
#というわけで安定判別します

#まずは連続時間系
Ev=np.linalg.eigvals(A)
print(Ev) #0を臨界点としているので、極が8.776875の場合は不安定

#次に離散時間系
Ev=np.linalg.eigvals(P.A)
print(Ev) #1を臨界点としているので、極が1.04486144の場合は不安定
