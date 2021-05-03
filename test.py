import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
N = 60000000        # 湖北省为6000 0000
E_0 = 0
I_0 = 1
R_0 = 0
D_0=0
S_0 = N - E_0 - I_0 - R_0
beta1 = 0.78735     # 真实数据拟合得出
beta2 = 0.15747
# r2 * beta2 = 2
sigma = 0.1         # 1/14, 潜伏期的倒数
gamma = 0.1         # 1/7, 感染期的倒数
gamma2 = 0.0005
             # 政府干预措施决定
T = 60

#ode求解
INI = [S_0, E_0, I_0, R_0 , D_0]
def SEIR(inivalue, t):
    X = inivalue
    Y = np.zeros(5)
    if t<10:
        r=5
        # S数量
        Y[0] = - (r * beta1 * X[0] * X[2]) / N - (r * beta2 * X[0] * X[1]) / N
        # E数量
        Y[1] = (r * beta1 * X[0] * X[2]) / N + (r * beta2 * X[0] * X[1]) / N - sigma * X[1]
        # I数量
        Y[2] = sigma * X[1] - gamma * X[2]
        # R数量
        Y[3] = gamma * X[2]
        # D
        Y[4] = gamma2 * X[2]
    else :
        r=1
        # S数量
        Y[0] = - (r * beta1 * X[0] * X[2]) / N - (r * beta2 * X[0] * X[1]) / N
        # E数量
        Y[1] = (r * beta1 * X[0] * X[2]) / N + (r * beta2 * X[0] * X[1]) / N - sigma * X[1]
        # I数量
        Y[2] = sigma * X[1] - gamma * X[2]
        # R数量
        Y[3] = gamma * X[2]
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