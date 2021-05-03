# -*- coding: utf-8 -*-
"""
Created on Sat May  1 11:52:55 2021
注射疫苗的状态
@author: pc
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import random
# accessrate = 0.3
# effecrate1 = 0.85#根据疫苗保护率80%转换得到
# effecrate2 = 0.85
arr = []
for accessrate in np.arange(0,1.1,0.1):
    for effecrate1 in  np.arange(0.1,1.1,0.1):
        effecrate2=effecrate1
        N = 100000          # 湖北省为6000 0000
        E_0 = 0
        I_0 = 1
        R_0 = 0
        A_0 = N*accessrate
        D_0 = 0
        I2_0=0
        S_0 = N- E_0 - I_0 - R_0
        beta1 = 0.78735     # 真实数据拟合得出---传染率
        beta2 = 0.15747   #潜伏期的传染率
        beta3 = 0.85  #---第二种病毒传染率
        eps = 0.1   #潜伏期内自愈率
        u=0.2 #从第一种转换为第二种病毒的概率
        # r2 * beta2 = 2
        sigma = 1/14         # 1/14, 潜伏期的倒数  潜伏期内发病的概率
        gamma = 0.03        #    治愈的概率
        gamma2 = 0.0015   #死亡率
        gamma3 = 0.003 #第二种病毒致死率
        r= 10         # 政府干预措施决定----每天接触的人
        eta=0.1#转化为I后没有隔离的比例
        T = 150

        #ode求解
        INI = [S_0, E_0, I_0, R_0 ,D_0,I2_0 ,A_0]
        # beta1 = 0.78735     # 真实数据拟合得出---传染率# beta2 = 0.15747   #潜伏期的传染率# beta3 = 0.85  #---第二种病毒传染率
        # eps = 0.1   #潜伏期内自愈率# u=0.2 #从第一种转换为第二种病毒的概率# # r2 * beta2 = 2
        # sigma = 1/14       潜伏期的倒数  潜伏期内发病的概率     # r= 1         # 政府干预措施决定----每天接触的人
        # gamma = 0.03  治愈的概率     #gamma2 = 0.0015   #死亡率      # gamma3 = 0.003 #第二种病毒致死率
        def SEIR(inivalue, t):
            X = inivalue
            Y = np.zeros(7)
            # S数量
            Y[0] = - (r * beta1 *eta*X[0] * X[2]) / (N-X[4]-(1-eta)*X[2]-(1-eta)*X[5]) - (r * beta2 *X[0] * X[1]) / (N-X[4]-(1-eta)*X[2]-(1-eta)*X[5])  - (r * beta3 *eta* X[0] * X[5]) /(N-X[4]-(1-eta)*X[2]-(1-eta)*X[5])
            # E数量
            Y[1] = (r * beta1 * X[0] *eta* X[2]) / (N-X[4]-(1-eta)*X[2]-(1-eta)*X[5])  + (r * beta2 * X[0] * X[1]) / (N-X[4]-(1-eta)*X[2]-(1-eta)*X[5]) + (r * beta3 *eta* X[0] * X[5]) / (N-X[4]-(1-eta)*X[2]-(1-eta)*X[5]) - sigma * X[1] - eps * X[1] + X[6]*r*eta*beta1*(X[2]*(1-effecrate1) + X[5]*eta*beta1*(1-effecrate2))/(N-X[4]-(1-eta)*X[2]-(1-eta)*X[5])
            # I数量
            Y[2] = sigma * X[1] - gamma * X[2] - gamma2 * X[2] - u * X[2]
            # R数量
            Y[3] = gamma * X[2] + eps * X[1] + gamma * X[5] - u*X[2]
            # D
            Y[4] = gamma2 * X[2] + gamma3 *X[5]
            # I2 感染第二种病毒的人
            Y[5] =u*X[2]-X[5]*(gamma+gamma3)
            # A 接种疫苗的人数
            Y[6]=-X[6]*r*(X[2]*eta*beta1*(1-effecrate1) + X[5]*eta*beta1*(1-effecrate2))/(N-X[4]-(1-eta)*X[2]-(1-eta)*X[5])
            return Y

        T_range = np.arange(0, T+1)
        Res = spi.odeint(SEIR, INI, T_range)
        S_t = Res[:, 0]
        E_t = Res[:, 1]
        I_t = Res[:, 2]
        R_t = Res[:, 3]
        D_t = Res[:, 4]
        I2_t = Res[:,5]
        A_t= Res[:,6]
        arr.append(max(D_t))
        #plt.plot(effecrate1,maxx,'-',color='blue')
        # for i in range(T-1):
        #     if I_t[i+1]-I_t[i]<0 and I_t[i+1]-I_t[i]>-10:
        #         print(i)
        # for i in range(T):
        #     if sigma*(E_t[i+1]-E_t[i])<0.5 and sigma*(E_t[i+1]-E_t[i])>-10:
        #         print(i)
        #         break
        # print(max(D_t))
        # plt.plot(S_t, color='blue', label='Susceptibles')#, marker='.')
        # plt.plot(E_t, color='grey', label='Exposed')
        # plt.plot(I_t, color='red', label='Infected')
        # plt.plot(R_t, color='green', label='Removed')
        # plt.plot(D_t,color='yellow',label='Death')
        # plt.plot(I2_t,color='black',label='Infected2')
        # plt.plot(A_t,color='pink',label='Accept')
        # plt.xlabel('Day')
        # plt.ylabel('Number')
        # plt.title('SEIR Model')
        # plt.legend()
effect=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
m=[]
print(len(arr))
for i in range(len(arr)):
    if i!=0 and i%10==0:
        plt.plot(effect,m,'-',color='blue',marker='.')
        x=int((1.1-(i/100))*10)/10
        plt.annotate(str(x),xy=(effect[0],m[0]))
        m.clear()
    m.append(arr[i])
plt.xlabel('effective')
plt.ylabel('InfectAll')
plt.show()


