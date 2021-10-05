#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 22:33:39 2020

@author: zhangyuxin
"""

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['SimHei']
#plt.rcParams['axes.unicode_minus'] = False
col = ['c','darkgreen','brown']

#t1 = [21.52,105.81,40.12]
#t2 = [21.01,104.92,40.15]
#t0 = [22.30,110.32,43.06]

#t1 = [12.31,54.62,19.65]
#t2 = [11.72,53.82,19.80]
#t0 = [13.42,57.82,21.31]
#
t1 = [8.42,38.2,14.21]
t2 = [7.5,37.2,14.21]
t0 = [9.43,40.65,15.6]

title = ['CG','SP','BT']
index0 = np.array([0,0.6,1.2])
index1 = index0 + 0.1
index2 = index1 + 0.1

plt.bar(index0,t0,color=col[0],width=0.1,label='Compact')
plt.bar(index1,t1,color=col[1],width=0.1,label='Eagermap')
plt.bar(index2,t2,color=col[2],width=0.1,label='LMC')
plt.xticks(index0 + 0.1, title)
plt.legend()
plt.ylabel("average execution time/s")
plt.show()



