#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 15:46:58 2020

@author: zhangyuxin
"""

import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import os

pit = np.loadtxt("bt.A.x.comm_time_count.csv",delimiter=',',dtype=int).T
time = pit[0].reshape(-1,1).copy()


# 线程对提取
thread_pairs = [(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),
               (12,13),(14,15),(16,18),(17,19),(20,21),
               (22,23),(24,25),(26,27),(28,29),(30,31)]

for root, dirs, files in os.walk('bt.A.x'):
    if ".DS_Store" in files:
        files.remove(".DS_Store")
    names = files

pair_dt = {}
time_pair_comm = {}
pair_comm = {}
for ti in list(time):
    tmp = np.loadtxt("bt.A.x/"+str(int(ti))+".csv",delimiter=',',dtype=int)
    if len(tmp.shape):
        for pair in thread_pairs:
            if pair[0]+1 <= tmp.shape[0] and pair[1]+1 <= tmp.shape[0]:
                if tmp[tmp.shape[0]-pair[0]-1][pair[1]]:
                    pair_comm[pair] = tmp[tmp.shape[0]-pair[0]-1][pair[1]]
                    if int(ti) in pair_dt:
                        pair_dt[int(ti)] += tmp[tmp.shape[0]-pair[0]-1][pair[1]]
                    else:
                        pair_dt[int(ti)] = tmp[tmp.shape[0]-pair[0]-1][pair[1]]
        
        time_pair_comm[int(ti)] = pair_comm;pair_comm = {}
pit[0] = np.arange(len(pit[0]))/len(pit[0]) * 100

p0 = pit[0].reshape(-1,1)
p2 = pit[1].reshape(-1,1)





# time 
tim = [];comm = []
for it in pair_dt.items():
    tim.append(it[0])
    comm.append(it[1])

#time = time[p2!=0].reshape(-1,1)
#c= p0[p2!=0]
#d = p2[p2!=0]
t = np.array(tim)
comm = np.array(comm)

#kmeans = XMeans(random_state = 1).fit(t)
kmeans = KMeans(n_clusters=3,max_iter=300,random_state=1).fit(t.reshape(-1,1),sample_weight=comm)
res = list(kmeans.labels_)

dc={}
# color bar
dic = {0:'g',1:'b',2:'r',3:'c',4:'m',5:'y',6:'k',7:'sandybrown',8:'pink',9:'teal',10:'orange',11:'gold',12:'grey',13:'honeydew',14:'lawngreen',15:'khaki'}
col=[]
for i in kmeans.labels_:
    col.append(dic[i])
    if i in dc:
        dc[i]+=1
    else:
        dc[i]=1


tim = np.array(tim)
tim = np.arange(len(tim))/len(tim) * 100
plt.bar(tim,comm,color=col,width=0.1)
plt.show()

pit[2] = -1;pit[1] = 0
for i,ic in enumerate(list(t)):
    pit[2][list(time).index(ic)] = kmeans.labels_[i]
    pit[1][list(time).index(ic)] = list(comm)[i]

col2 = [];count = 0
for i in list(pit[2]):
    if i == -1: 
        col2.append('lime')
    else:
        col2.append(col[count])
        count+=1


plt.bar(list(pit[0]),list(pit[1]),color=col2,width=0.5)

plt.show()



# 组内总通信量
key = kmeans.labels_[0]
group_comm = {};s = 0
group_count = {}
group_timelist = {}
tmplist = []
for i in range(len(comm)):
    if kmeans.labels_[i] in group_count:
        group_count[kmeans.labels_[i]]+=1
    else:
        group_count[kmeans.labels_[i]] =1
    if kmeans.labels_[i] == key:
        s += comm[i]
        tmplist.append(int(t[i]))
    else:
        group_comm[key] = s
        group_timelist[key] = tmplist
        tmplist = [];tmplist.append(int(t[i]))
        key = kmeans.labels_[i]
        s=0;s+=comm[i]
group_comm[key] = s
group_timelist[key] = tmplist

# 组内线程对通信量统计
tmp = 0
group_thread_pair ={}
pair_comm = {}
for item in group_timelist.items():
    for te in item[1]:
        for pair in thread_pairs:
            if pair in time_pair_comm[te].keys():
                if pair in pair_comm:
                    pair_comm[pair] += time_pair_comm[te][pair]
                else:
                    pair_comm[pair] = time_pair_comm[te][pair]
    
    group_thread_pair[item[0]] = pair_comm;pair_comm={}
            
#gruop_pairs = [[(30,31),(24,25),(14,15),(12,13),(10,11),(8,9),(6,7),(26,27)],[(0,1),(16,19),(22,23),(20,21),(17,18),(2,3),(4,5),(28,29)]]
#map()

def map_(gruop_pairs,mat):
    node = 2
    group_res = {0:[],1:[]}
    for i in range(node):
        chosen = [1] * len(gruop_pairs[i])
        for j in range(len(gruop_pairs[i])):
            wmax = -1
            if chosen[j]:
                for k in range(len(gruop_pairs[i])):
                    if j!=k and chosen[k]:
                        m1 = gruop_pairs[i][j][0]
                        m2 = gruop_pairs[i][j][1]
                        m3 = gruop_pairs[i][k][0]
                        m4 = gruop_pairs[i][k][1]
                        w = group_comm(m1,m2,m3,m4,mat)
                        if wmax < w:
                            wmax = w
                            tmp_k = k
                            tmp = [(m1,m2),(m3,m4)].copy()
                group_res[i].append(tmp)
                chosen[j] = 0;chosen[tmp_k] = 0
    return group_res


def group_comm(i,j,m,n,mat):
    return mat[mat.shape[0]-i][m] + mat[mat.shape[0]-i][n] + mat[mat.shape[0]-j][m] + mat[mat.shape[0]-j][n]


def define_k(t,comm):
    SSE = []  # 存放每次结果的误差平方和
    for k in range(2,9):
        estimator = KMeans(n_clusters=k)  # 构造聚类器
        estimator.fit(t.reshape(-1,1),sample_weight=comm)
        SSE.append(estimator.inertia_)
    X = range(2,9)
    plt.xlabel('k')
    plt.ylabel('SSE')
    plt.plot(X,SSE,'o-')
    plt.show()





