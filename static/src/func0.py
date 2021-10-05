#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:47:44 2020

@author: zhangyuxin
"""
from sklearn.cluster import KMeans
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys


title = ['timestamp','thread1','thread2']
data = pd.read_csv("cg.B.x-10539-timestamp.csv",header=None,names=title)


# 线程对提取
thread_pairs = [(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),
               (12,13),(14,15)]

# 遍历pandas 

df = data
for index, row in df.iterrows():
    if ((row['thread1'],row['thread2']) not in thread_pairs) and ((row['thread2'],row['thread1']) not in thread_pairs):
        data = data.drop(index=(data.loc[(data['timestamp']==row['timestamp'])].index))
# 聚类
t = np.array(data['timestamp']).reshape(-1,1)
kmeans = KMeans(n_clusters=7,max_iter=300,random_state=0).fit(t)
res = list(kmeans.labels_)

# 组（类）内总通信量
before = res[0]
group_comm = {}
group_comm[before] = 0
group_timestamp = {}
tmplist = []
for i in range(len(res)):
    if res[i] == before:
       group_comm[res[i]]+=1
       tmplist.append(int(t[i]))
    else:
        group_timestamp[before] = tmplist
        tmplist = []
        before = res[i]
        group_comm[before]=1
        tmplist.append(int(t[i]))

group_timestamp[before] = tmplist

# 组（类）内线程对通信量统计
tmp = 0
group_thread_pair ={}
pair_comm = {}
for it in group_timestamp.items():
    for time in it[1]:
        tmp = data.loc[data["timestamp"] == time]
        t1 = int(tmp['thread1']);t2 = int(tmp['thread2'])
        if (t1,t2) in pair_comm :
            pair_comm[(t1,t2)] += 1
        elif (t2,t1) in pair_comm :
            pair_comm[(t2,t1)] += 1
        else:
            pair_comm[(t1,t2)] = 1
    group_thread_pair[it[0]] = pair_comm
    pair_comm = {}

        

def define_k(t):
    SSE = []  # 存放每次结果的误差平方和
    for k in range(2,9):
        estimator = KMeans(n_clusters=k)  # 构造聚类器
        estimator.fit(t.reshape(-1,1))
        SSE.append(estimator.inertia_)
    X = range(2,9)
    plt.xlabel('k')
    plt.ylabel('SSE')
    plt.plot(X,SSE,'o-')
    plt.show()
    
mat = np.loadtxt("1592494356_timestamped_results/bt.B.x-21907-timestamp.csv",delimiter=',',dtype=int)

pu_map = {0:0, 1:16, 2:1, 3:17, 4:2, 5:18, 6:3, 7:19, 8:4, 9:20, 10:5, 11:21, 12:6, 13:22, 14:7, 15: 23, 16:8, 17:24, 
          18:9, 19:25, 20:10, 21:26, 22:11, 23:27, 24:12, 25:28, 26:13, 27:29, 28:14, 29:30, 30:15, 31:31}

l = [4,2,6,8,18,9,5,20,0,10,7,1,3,14,19,22,24,12,23,16,28,29,13,15,26,17,21,30,27,11,31,25]
#l = list(range(32))
p = []
for i in l:
    p.append(pu_map[i])

ave_time_cut = np.zeros((8,8),dtype=int)

for index in range(len(thread_pairs)):
    for j in range(index+1,len(thread_pairs)):
        p11=thread_pairs[index][0];p12=thread_pairs[index][1]
        p21=thread_pairs[j][0];p22=thread_pairs[j][1]
        d1 = data[(((data["thread1"] == p11) & (data["thread2"] == p12)) | ((data["thread1"] == p12) & (data["thread2"] == p11)))]
        d2 = data[(((data["thread1"] == p21) & (data["thread2"] == p22)) | ((data["thread1"] == p22) & (data["thread2"] == p21)))]
        
        if not (d1.shape[0] and d2.shape[0]):
            ave_time_cut[7-index][j] = sys.maxsize
            ave_time_cut[7-j][index] = sys.maxsize
            continue
        
        d = d1.append(d2)
        d.sort_values(by='timestamp',inplace=True,ascending=True)
        
        tmp = 0
        res = []
        for i in range(len(d)-1):
            row1 = d.iloc[i]
            row2 = d.iloc[i+1]
            t11 = int(row1['thread1']);t12 = int(row1['thread2'])
            t21 = int(row2['thread1']);t22 = int(row2['thread2'])
            # 线程对相同
            if ((t11,t12) == (t21,t22)) or ((t12,t11) == (t21,t22)):
                if tmp:
                    res.append(tmp);tmp=0
                else:
                    pass
            # 线程对不同
            else:
                if tmp:
                    tmp = min(tmp,row2['timestamp'] - row1['timestamp'])
                    res.append(tmp);tmp=0
                else:
                    tmp = row2['timestamp'] - row1['timestamp']
        if tmp: 
            res.append(tmp)
            tmp=0
        
        ave = int(sum(res)/len(res))
        ave_time_cut[7-index][j] = ave
        ave_time_cut[7-j][index] = ave
df = pd.DataFrame(ave_time_cut)
df.to_csv("res.csv",header=False,index=False)
    
