# -*- coding: utf-8 -*-
"""
Created on Sat May  1 11:52:55 2021
未注射疫苗的状态
@author: pc
"""
#截至目前已确诊人数103626，累计死亡4857
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
N = 326766748          # 湖北省为6000 0000
E_0 = 140000
I_0 = 1067289
R_0 = 125949
D_0 = 68000
S_0 = N - E_0 - I_0 - R_0
beta1 = 0.78735     # 真实数据拟合得出---传染率
beta2 = 0.15747   #潜伏期的传染率
eps = 0.01  #无症状自愈率
# r2 * beta2 = 2
sigma = 1/14         # 1/14, 潜伏期的倒数  潜伏期内发病的概率
gamma = 0.03        #    治愈的概率
gamma2 = 0.005   #死亡率
r=  10     # 政府干预措施决定

T = 60

#ode求解
INI = [S_0, E_0, I_0, R_0 ,D_0]
def SEIR(inivalue, t):
    X = inivalue
    Y = np.zeros(5)
    # S数量
    Y[0] = - (r * beta1 * X[0] * 0.1*X[2]) / (N-X[4]-0.9*X[2]) - (r * beta2 * X[0] *X[1]) / (N-X[4]-0.9*X[2])
    # E数量
    Y[1] = (r * beta1 * X[0] * 0.1*X[2]) / (N-X[4]-0.9*X[2]) + (r * beta2 * X[0] *X[1]) / (N-X[4]-0.9*X[2])- sigma * X[1] - eps * X[1]
    # I数量
    Y[2] = sigma * X[1] - gamma * X[2] - gamma2 * X[2]
    # R数量
    Y[3] = gamma * X[2] + eps * X[1]
    #D
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
for i in range(T):
    if sigma*(E_t[i+1]-E_t[i])<0.5 and sigma*(E_t[i+1]-E_t[i])>-10:
        print(i)
        print(D_t[i])
        break
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