#!/usr/bin/python
# -*- coding: UTF-8 -*-
"""
python解方程
"""

from scipy.optimize import fsolve

# r治愈率  d死亡率
def solve_function(unsolved_value):
    r, d, y = unsolved_value[0], unsolved_value[1], unsolved_value[2]
    return [
       r*(1-y**31)/(1-y)-0.78,
        d*(1-y**21)/(1-y)-0.018,
        r+d+y-1,
    ]


solved = fsolve(solve_function, [0,0,0])
print(solved)

print("Program done!")

"""
运行结果：
[-1.  3.  5.]
Program done!
"""