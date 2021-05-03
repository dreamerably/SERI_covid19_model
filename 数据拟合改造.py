import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# 数据预处理常用库
from sklearn import preprocessing
import time
from datetime import datetime
from scipy import integrate, optimize
train = pd.read_csv(r"C:\Users\pc\Desktop\covid19\train.csv")
#train.Province_State.fillna("None", inplace=True)
population = float(46750238)
country_df = pd.DataFrame()
country_df['ConfirmedCases'] = train.loc[train['Province_State']=='Hubei'].ConfirmedCases.diff().fillna(0)
country_df = country_df[5:]
country_df['day_count'] = list(range(1,len(country_df)+1))

ydata = [i for i in country_df.ConfirmedCases]
xdata = country_df.day_count
ydata = np.array(ydata, dtype=float)
xdata = np.array(xdata, dtype=float)


N = population
inf0 = ydata[0]
sus0 = N - inf0
rec0 = 0.0

def sir_model(y, x, beta, gamma):
    sus = -beta * y[0] * y[1] / N
    rec = gamma * y[1]
    inf = -(sus + rec)
    return sus, inf, rec

def fit_odeint(x, beta, gamma):
    return integrate.odeint(sir_model, (sus0, inf0, rec0), x, args=(beta, gamma))[:,1]

popt, pcov = optimize.curve_fit(fit_odeint, xdata, ydata, maxfev=500000)
fitted = fit_odeint(xdata, *popt)

plt.plot(xdata, ydata, 'o')
plt.plot(xdata, fitted)
plt.title("SIR 模型根据西班牙感染者数量进行拟合")
plt.ylabel("感染人口数")
plt.xlabel("天数")
plt.show()
print("最佳参数: beta =", popt[0], " and gamma = ", popt[1])