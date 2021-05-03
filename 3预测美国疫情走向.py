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
N = 330000000      #美国人口
E_0 = 700000
I_0 = 8923624
R_0 = 20003419
D_0 = 531682
A_0 = 40540474
S_0 = N - E_0 - I_0 - R_0
eps = 0  #无症状自愈率
# r2 * beta2 = 2
eta = 1*0.1 # 转化为I后没有隔离的比例
beta1 = 0.004 # 真实数据拟合得出---传染率
beta2 = 0.0006  # 潜伏期的传染率
sigma = 1/14         # 1/14, 潜伏期的倒数  潜伏期内发病的概率
gamma = 0.05        #    治愈的概率
gamma2 = 0.00035   #死亡率
zeta=1200000 #接种速率
r=  50     # 政府干预措施决定
T = 30
cnt=0
#ode求解
INI = [S_0, E_0, I_0, R_0 ,D_0,A_0]
def SEIR(inivalue, t):
    X = inivalue
    Y = np.zeros(6)
    global eta
    global beta1
    global beta2
    global zeta
    global cnt
    #S
    Y[0] = - (r * beta1 * X[0] * eta * X[2]) / (N - X[4] - (1 - eta) * X[2]) - (r * beta2 * X[0] * X[1]) / (
                    N - X[4] - (1 - eta) * X[2]) - zeta * X[0] / (X[0] + X[1])
        # E数量
    Y[1] = (r * beta1 * X[0] * eta * X[2]) / (N - X[4] - (1 - eta) * X[2]) + (r * beta2 * X[0]  * X[1]) / (
                    N - X[4] - (1 - eta) * X[2]) - sigma * X[1] - eps * X[1] - zeta * X[1] / (X[0] + X[1])
        # I数量
    Y[2] = sigma * X[1] - gamma * X[2] - gamma2 * X[2] + zeta * X[1] / (X[0] + X[1])
        # R数量
    Y[3] = gamma * X[2] + eps * X[1]
        # D
    Y[4] = gamma2 * X[2]
        # A 接种人数
    Y[5] = zeta * X[0] / (X[0] + X[1])  # 考虑到注射疫苗前需要核算检测
    # if sigma * X[1] - gamma * X[2] - gamma2 * X[2] + zeta * X[1] / (X[0] + X[1])>100:
    #     eta=0.5
    #     beta1 = 0.78735*0.05  # 真实数据拟合得出---传染率
    #     beta2 = 0.15747*0.05 # 潜伏期的传染率
    if X[0]<40000000:#老人与小孩不宜接种
        zeta=0
    return Y

T_range = np.arange(0, T+1)
Res = spi.odeint(SEIR, INI, T_range)
S_t = Res[:, 0]
E_t = Res[:, 1]
I_t = Res[:, 2]
R_t = Res[:, 3]
D_t = Res[:, 4]
A_t = Res[:, 5]
for i in range(T):
    if sigma*(E_t[i+1]-E_t[i])>100:
        print(i)
        print(D_t[i])
        break
print(max(R_t))
print(max(D_t))
print(max(A_t))
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