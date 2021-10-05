#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 15:08:11 2020

@author: zhangyuxin
"""
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import numpy as np

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

class Add():
    def __init__(self,a,b):
        self.a=a
        self.b=b
    def add(self):
        return self.a+self.b

titie = ['thread_id','time_stamp','address','mem_level','cost','acces_type']
data1 = pd.read_csv("text 3.dat",header=None,delimiter=' ',names=titie)

t = np.array(data1['time_stamp']).reshape(-1,1)
t = np.sort(t,axis=0)
kmeans = KMeans(n_clusters=3,max_iter=300,random_state=0).fit(t)
res = list(kmeans.labels_)

# 线程对提取
thread_pairs = [(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),
               (12,13),(14,15)]

df = data1
for index, row in df.iterrows():
    for tup in thread_pairs:
        if (row['thread_id'] == tup[0]) or (row['thread_id'] == tup[1]):
            #row['thread_id'] = thread_pairs.index(tup)
            df.loc[index,'thread_id'] = thread_pairs.index(tup)
            break
        
data1 = df
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


data1 = data1.sort_values(by='time_stamp',ascending=True)
#d = data1[508:1063]
d = data1
dc = {}

for index,row in d.iterrows():
    if row['thread_id'] in dc:
        dc[row['thread_id']]+=1
    else:
        dc[row['thread_id']]=1


dc = sorted(dc.items(),key=lambda x:x[1],reverse=True)

bind_list = []
for i in dc:
    bind_list.append(i[0]*2)
    bind_list.append(i[0]*2+1)
pu_map = {0:0, 1:16, 2:1, 3:17, 4:2, 5:18, 6:3, 7:19, 8:4, 9:20, 10:5, 11:21, 12:6, 13:22, 14:7, 15: 23, 16:8, 17:24, 
          18:9, 19:25, 20:10, 21:26, 22:11, 23:27, 24:12, 25:28, 26:13, 27:29, 28:14, 29:30, 30:15, 31:31}
l=[10,11,8,9,12,13,4,5,2,3,6,7,0,1,14,15]
p = []
for i in l:
    p.append(pu_map[i])
    
    
#ave_time_cut = np.zeros((8,8),dtype=int)
#
#for index in range(len(thread_pairs)):
#    for j in range(index+1,len(thread_pairs)):
#        p11=thread_pairs[index][0];p12=thread_pairs[index][1]
#        p21=thread_pairs[j][0];p22=thread_pairs[j][1]
#        d1 = data[(((data["thread1"] == p11) & (data["thread2"] == p12)) | ((data["thread1"] == p12) & (data["thread2"] == p11)))]
#        d2 = data[(((data["thread1"] == p21) & (data["thread2"] == p22)) | ((data["thread1"] == p22) & (data["thread2"] == p21)))]
#        
#        if not (d1.shape[0] and d2.shape[0]):
#            ave_time_cut[15-index][j] = sys.maxsize
#            ave_time_cut[15-j][index] = sys.maxsize
#            continue
#        
#        d = d1.append(d2)
#        d.sort_values(by='timestamp',inplace=True,ascending=True)
#        
#        tmp = 0
#        res = []
#        for i in range(len(d)-1):
#            row1 = d.iloc[i]
#            row2 = d.iloc[i+1]
#            t11 = int(row1['thread1']);t12 = int(row1['thread2'])
#            t21 = int(row2['thread1']);t22 = int(row2['thread2'])
#            # 线程对相同
#            if ((t11,t12) == (t21,t22)) or ((t12,t11) == (t21,t22)):
#                if tmp:
#                    res.append(tmp);tmp=0
#                else:
#                    pass
#            # 线程对不同
#            else:
#                if tmp:
#                    tmp = min(tmp,row2['timestamp'] - row1['timestamp'])
#                    res.append(tmp);tmp=0
#                else:
#                    tmp = row2['timestamp'] - row1['timestamp']
#        if tmp: 
#            res.append(tmp)
#            tmp=0
#        
#        ave = int(sum(res)/len(res))
#        ave_time_cut[15-index][j] = len(res)
#        ave_time_cut[15-j][index] = len(res)
#df = pd.DataFrame(ave_time_cut)
#df.to_csv("res.csv",header=False,index=False)

