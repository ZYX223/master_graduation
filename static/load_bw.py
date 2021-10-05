#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 21:17:48 2021

@author: zhangyuxin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_txt(file_name):
    lines = open(file_name)
    x=[];y=[]
    for line in lines:
        if line[0] == '#':
            continue
        if line[0] == 'E':
            break
        tmp_list = line.split('\t')
        x.append(eval(tmp_list[2][1:-1]))
        y.append(eval(tmp_list[1]))
    return x,y


X,Y = load_txt("node0.txt")
M,N = load_txt("node1.txt")
P,Q = load_txt("all node.txt")

plt.title('loaded_latency')
plt.xlabel('bandwidth[mb/s]')
plt.ylabel('latency[ns]')
#plt.plot(X,Y)
#plt.plot(M,N)
plt.plot(P,Q)
plt.show()