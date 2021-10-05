#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 16:18:35 2020

@author: zhangyuxin
"""

import numpy as np
import pandas as pd

def generate_all():
    feq_list = [0]*num
    for time in time_list:
        tlist = list(data[data['timestamp']==time]['thread'])
        for tid in tlist:
            feq_list[tid]+=1
    for i in range(num):
        for j in range(num):
            if i != j:
                mat[num-1-i][j] += (feq_list[i]+feq_list[j])
    np.savetxt(file_name+"_time_mat.csv",mat,fmt='%d',delimiter=',')
    return

num=32
# sp.B.x-19277-time_line
# sp.B.x-30824-time_line
file_name="lu.B.x-14527-time_line.csv"
path = "timelinefile/"+file_name
# load file
title = ['timestamp','thread','accesslen']
data = pd.read_csv(path,header=None,names=title)

# produce data
time_slice=1000000
data['timestamp'] //= time_slice
data = data[data['timestamp']<int(4e12)]

# timeline
df = data.groupby(['timestamp']).count()
df=df.reset_index()
df['access_count'] = df['accesslen']
df = df.drop(['thread','accesslen'],axis=1)

# generate
mat = np.array(np.zeros((num,num)),dtype=int)
time_list = list(df['timestamp'])
generate_all()
#time_list = time_list[470:520]
#ind = 0
#for time in time_list:
#    tlist = list(data[data['timestamp']==time]['thread'])
#    ls = [0]*num
#    for tid in tlist:
#        ls[tid]+=1
#    for i in range(num):
#        for j in range(num):
#            if i != j:
#                mat[num-1-i][j] += (ls[i]+ls[j])
#    ind+=1
#    print("{}/{}".format(ind,len(time_list)))
#np.savetxt(file_name+"_time_mat.csv",mat,fmt='%d',delimiter=',')
