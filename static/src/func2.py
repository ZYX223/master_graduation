#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 18:50:42 2020

@author: zhangyuxin
"""

import pandas as pd

data = pd.read_csv("1593007632_timestamped_results/cg.B.x-1706-as_matrix.csv",header=None)


thread_pair=[(0,6),(1,12),(2,10),(3,5),(4,14),(7,26),(8,11),(9,29),(13,23),(15,18),(16,31),(17,22),(19,25),(20,21),(24,28),(27,30)]

dt = {}
for item in thread_pair:
    dt[thread_pair.index(item)] = data[data.shape[0]-1-item[0]][item[1]]
    
dt = sorted(dt.items(),key=lambda x:x[1],reverse=True)

th = []
for i in dt:
    th.append(thread_pair[i[0]])
    
bind = {}

for i in range(len(th)):
    bind[i] = th[i]
bind_map = {}
bind = sorted(bind.items(),key=lambda x:x[1][0],reverse=False)
for item in bind:
    bind_map[item[1][0]] = item[0]*2
    bind_map[item[1][1]] = item[0]*2+1
bind_map = sorted(bind_map.items(),key=lambda x:x[0],reverse=False)
bind_list = []
for i in bind_map:
    bind_list.append(i[1])
pu_map = {0:0, 1:16, 2:1, 3:17, 4:2, 5:18, 6:3, 7:19, 8:4, 9:20, 10:5, 11:21, 12:6, 13:22, 14:7, 15: 23, 16:8, 17:24, 
          18:9, 19:25, 20:10, 21:26, 22:11, 23:27, 24:12, 25:28, 26:13, 27:29, 28:14, 29:30, 30:15, 31:31}
p = []
for i in bind_list:
    p.append(pu_map[i])
