# -*- coding: utf-8 -*-
"""
Created on Sat May  1 11:52:55 2021
注射疫苗的状态
@author: pc
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
accessrate = 0.8
N = 100000          # 湖北省为6000 0000
E_0 = 0
I_0 = 1
R_0 = N*accessrate
D_0 = 0
S_0 = N- E_0 - I_0 - R_0
beta1 = 0.78735     # 真实数据拟合得出
beta2 = 0.15747
eps = 0.1
# r2 * beta2 = 2
sigma = 1/14         # 1/14, 潜伏期的倒数  潜伏期内发病的概率
gamma = 0.03        #    治愈的概率
gamma2 = 0.0015   #死亡率
r= 100         # 政府干预措施决定----每天接触的人

T = 300

#ode求解
INI = [S_0, E_0, I_0, R_0 ,D_0]
cnt=0
def SEIR(inivalue, t):
    X = inivalue
    Y = np.zeros(5)
    # S数量
    Y[0] = - (r * beta1 * X[0] * X[2]) / (N-X[4]) - (r * beta2 * X[0] * X[1]) / (N-X[4])
    # E数量
    Y[1] = (r * beta1 * X[0] * X[2]) / (N-X[4]) + (r * beta2 * X[0] * X[1]) / (N-X[4]) - sigma * X[1] - eps * X[1]
    # I数量
    Y[2] = sigma * X[1] - gamma * X[2] - gamma2 * X[2]
    # R数量
    Y[3] = gamma * X[2] + eps * X[1]
    # D
    Y[4] = gamma2 * X[2]
    return Y

T_range = np.arange(0, T+1)
Res = spi.odeint(SEIR, INI, T_range)
S_t = Res[:, 0]
E_t = Res[:, 1]
I_t = Res[:, 2]
R_t = Res[:, 3]
D_t = Res[:, 4]
# for i in range(T-1):
#     if I_t[i+1]-I_t[i]<0 and I_t[i+1]-I_t[i]>-10:
#         print(i)
print(max(D_t))
plt.plot(S_t, color='blue', label='Susceptibles')#, marker='.')
plt.plot(E_t, color='grey', label='Exposed')
plt.plot(I_t, color='red', label='Infected')
plt.plot(R_t, color='green', label='Removed')
plt.plot(D_t,color='yellow',label='Death')
plt.xlabel('Day')
plt.ylabel('Number')
plt.title('SEIR Model')
plt.legend()
plt.show()