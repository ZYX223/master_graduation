#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 21:22:13 2020

@author: zhangyuxin
"""
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np


pit = np.loadtxt("sp.W.x.comm_time_count.csv",delimiter=',').T
time = pit[0].reshape(-1,1).copy()

pit[0] = np.arange(len(pit[0]))/len(pit[0]) * 100

p0 = pit[0].reshape(-1,1)
p2 = pit[2].reshape(-1,1)
time = time[p2!=0]
c= p0[p2!=0]
d = p2[p2!=0]
t = time-time[0]
kmeans = KMeans(n_clusters=16,max_iter=300,random_state=1).fit(t.reshape(-1,1))

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

pit[1] = -1
#plt.bar(c,d,color=col,width=0.3)
for i,ic in enumerate(c):
    pit[1][list(pit[0]).index(ic)] = kmeans.labels_[i]

col2 = [];count = 0
for i in list(pit[1]):
    if i == -1: 
        col2.append('lime')
    else:
        col2.append(col[count])
        count+=1


plt.bar(list(pit[0]),list(pit[2]),color=col2,width=0.01)

plt.show()

# 组内总通信量
key = kmeans.labels_[0]
group_comm = {};s = 0
group_count = {}
group_timelist = {}
tmplist = []
for i in range(len(d)):
    if kmeans.labels_[i] in group_count:
        group_count[kmeans.labels_[i]]+=1
    else:
        group_count[kmeans.labels_[i]] =1
    if kmeans.labels_[i] == key:
        s += d[i]
        tmplist.append(int(time[i]))
    else:
        group_comm[key] = s
        group_timelist[key] = tmplist
        tmplist = [];tmplist.append(int(time[i]))
        key = kmeans.labels_[i]
        s=0;s+=d[i]
group_comm[key] = s
group_timelist[key] = tmplist


"""
# 改写文件名（如果需要的话）
def rename_file():
    for root, dirs, files in os.walk('pit'):  
        root+='/'
        if ".DS_Store" in files:
            files.remove(".DS_Store")
        for file in files:
            dstfile = file.split('.')
            name = str(int(eval(dstfile[0])))+'.csv'
            try:
                os.rename(root+file,root+name)
            except Exception as e:
                print(e)
                print('rename file fail\r\n')

#rename_file()
# 组内线程通信情况
threads_num = 8
for root, dirs, files in os.walk('pit'):
    if ".DS_Store" in files:
        files.remove(".DS_Store")
    names = files
     
group_comm_mat = {}
data = np.zeros((threads_num,threads_num),dtype = int)
for item in group_timelist.items():
    for file in item[1]:
        try:
            tmp= np.loadtxt("pit/"+str(file)+".csv",delimiter=',',dtype=int)
            if tmp.shape[0]< data.shape[0]:
                data[-tmp.shape[0]:,:tmp.shape[0]]+=tmp
            else:
                data+=tmp
        except:
#            for t in range(1,6):
#                 if str(file+t)+".csv" in names:
#                     break
#            tmp= np.loadtxt("pit/"+str(file+t)+".csv",delimiter=',',dtype=int)
            if tmp.shape[0]< data.shape[0]:
                data[-tmp.shape[0]:,:tmp.shape[0]]+=tmp
            else:
                data+=tmp
    group_comm_mat[item[0]] = data
    data = np.zeros((threads_num,threads_num),dtype = int)
        
# 组内线程对通信量统计
group_thread_pair ={}
pair_count = {}
for item in group_comm_mat.items():
    mat = item[1]
    for i in range(threads_num):
        for j in range(i+1,threads_num):
            pair_count[(i,j)] = mat[threads_num-1-i][j]
    group_thread_pair[item[0]] = pair_count
    pair_count = {}
"""


