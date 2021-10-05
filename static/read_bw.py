#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 21:53:03 2021

@author: zhangyuxin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_bw(path):
    data = pd.read_csv(path)
    bw_0 = np.array(data['SKT0.7'][7:],dtype=float)
    ind = np.argwhere(bw_0 >=1000)
    start_time = ind[0][0];end_time = ind[-1][0]
    bw_0 = bw_0[start_time:end_time+1]
    bw_1 = np.array(data['SKT1.7'][7:],dtype=float)
    bw_1 = bw_1[start_time:end_time+1]
    #bw_total = bw_0+bw_1
#    plt.plot(range(len(bw_0)),bw_0,label='skt0')
#    plt.plot(range(len(bw_1)),bw_1,label='skt1')
#    print("skt diff:",np.sum(np.abs(bw_0[0:30]-bw_1[0:30])))
#    print("skt diff:",np.sum(np.abs(bw_0[30:]-bw_1[30:])))
    #plt.plot(range(len(bw_1)),bw_total)
    return bw_0,bw_1

def read_qpi():
    path="qpi.csv"
    data = pd.read_csv(path)
    lmb0 = np.array(data['Unnamed: 47'][1:],dtype=float)
    rmb0 = np.array(data['Unnamed: 48'][1:],dtype=float)
    lmb1 = np.array(data['Unnamed: 64'][1:],dtype=float)
    rmb1 = np.array(data['Unnamed: 65'][1:],dtype=float)
    read_0 = np.array(data['Unnamed: 49'][1:],dtype=float)
    write_0 = np.array(data['Unnamed: 50'][1:],dtype=float)
    read_1 = np.array(data['Unnamed: 66'][1:],dtype=float)
    write_1 = np.array(data['Unnamed: 67'][1:],dtype=float)
    qpi =  np.array(data['Unnamed: 22'][1:],dtype=float)
    qpi0 = np.array(data['SKT0trafficOut'][1:],dtype=float)
    qpi1 = np.array(data['SKT1trafficOut'][1:],dtype=float)
    llc_miss_lantcy_0 = np.array(data['LLCRDMISSLAT (ns)'][1:],dtype=float)
    llc_miss_lantcy_1 = np.array(data['Unnamed: 107'][1:],dtype=float)
    local_access = np.array(data['Unnamed: 14'][1:],dtype=int)
    
#    plt.plot(range(len(lmb0)),lmb0+rmb0,label='skt0_total')
#    plt.plot(range(len(lmb0)),lmb1+rmb1,label='skt1_total')
#    plt.plot(range(len(rmb0)),llc_miss_lantcy_0,label='skt0_lantcy')
#    plt.plot(range(len(rmb1)),llc_miss_lantcy_1,label='skt1_lantcy')
    s0 = 5;e0 = 61
    s1 = 66;e1 = 133
    s2 = 102;e2 = 140
    s3 = 150;e3 = 188
    
    std0 = np.std([read_0[s0:e0+1],read_1[s0:e0+1]],axis=0,ddof=1)
    std1 = np.std([read_0[s1:e1+1],read_1[s1:e1+1]],axis=0,ddof=1)
#    std2 = np.std([read_0[s2:e2+1],read_1[s2:e2+1]],axis=0,ddof=1)
#    std3 = np.std([read_0[s3:e3+1],read_1[s3:e3+1]],axis=0,ddof=1)
    
    sum_0 = np.sum(read_0[s0:e0+1] + read_1[s0:e0+1])
    sum_1 = np.sum(read_0[s1:e1+1] + read_1[s1:e1+1])
#    sum_2 = np.sum(read_0[s2:e2+1] + read_1[s2:e2+1])
#    sum_3 = np.sum(read_0[s3:e3+1] + read_1[s3:e3+1])
    
#    print("memory access sum:",sum_0,sum_1,sum_2,sum_3)
    ave_std0 = np.mean(std0)
    ave_std1 = np.mean(std1)
#    ave_std2 = np.mean(std2)
#    ave_std3 = np.mean(std3)
    
    print("mc imbalabce:",ave_std0,ave_std1)
    qpi_0 = np.sum(qpi[s0:e0+1])
    qpi_1 = np.sum(qpi[s1:e1+1])
#    qpi_2 = np.sum(qpi[s2:e2+1])
#    qpi_3 = np.sum(qpi[s3:e3+1])
    qpi_imporve1 = (qpi_0 - qpi_1) / qpi_0
#    qpi_imporve2 = (qpi_0 - qpi_2) / qpi_0
#    qpi_imporve3 = (qpi_0 - qpi_3) / qpi_0
    print("qpi volume:",qpi_0,qpi_1)
#    print("qpi improve:",qpi_imporve1,qpi_imporve2,qpi_imporve3)
#    
    plt.plot(range(len(lmb0)),read_0*1000,label='skt0_total')
    plt.plot(range(len(lmb0)),read_1*1000,label='skt1_total')
    local_0 = np.mean(local_access[s0:e0+1])
    local_1 = np.mean(local_access[s1:e1+1])
#    local_2 = np.mean(local_access[s2:e2+1])
#    local_3 = np.mean(local_access[s3:e3+1])

#    print("local memory access ratio:",local_0,local_1,local_2,local_3)
#    plt.plot(range(len(rmb1)),qpi0,label='qpi_skt0')
#    plt.plot(range(len(rmb1)),qpi1,label='qpi_skt1')
    l0 = llc_miss_lantcy_0[s0:e0+1]+llc_miss_lantcy_1[s0:e0+1]
    l1 = llc_miss_lantcy_0[s1:e1+1]+llc_miss_lantcy_1[s1:e1+1]
#    l2 = llc_miss_lantcy_0[s2:e2+1]+llc_miss_lantcy_1[s2:e2+1]
#    l3 = llc_miss_lantcy_0[s3:e3+1]+llc_miss_lantcy_1[s3:e3+1]

    ave_0 = np.mean(l0);ave_1 = np.mean(l1)
#    ave_2 = np.mean(l2);ave_3 = np.mean(l3)
    print("ave_latency:",ave_0 /2,ave_1 /2)
#    plt.plot(range(len(rmb1)),qpi)
#    
    
    plt.legend()
    return qpi0,qpi1

def read_latency():
    path = "latency.csv"
    title = ['skt0','skt1']
    data = pd.read_csv(path,header=None,names=title)
    s0 = 1;e0 = 70
    s1 = 75;e1 = 132
    s2 = 137;e2 = 194
#    s3 = 231;e3 = 298
    skt0_latency = np.array(data['skt0'],dtype=int)
    skt1_latency = np.array(data['skt1'],dtype=int)
    
    l0 = skt0_latency[s0:e0+1]+skt1_latency[s0:e0+1]
    l1 = skt0_latency[s1:e1+1]+skt1_latency[s1:e1+1]
    l2 = skt0_latency[s2:e2+1]+skt1_latency[s2:e2+1]
#    l3 = skt0_latency[s3:e3+1]+skt1_latency[s3:e3+1]
    ave_0 = np.mean(l0);ave_1 = np.mean(l1)
    ave_2 = np.mean(l2)
    
    print("ave_latency:",ave_0,ave_1,ave_2)
    return

def read_numa():
    path  = "numa.csv"
    data = pd.read_csv(path,error_bad_lines=False)
    start=170;end=2344
    loc_acc = np.array(data['Local DRAM accesses'][start:end])
    remote_acc = np.array(data['Remote DRAM accesses'][start:end])
    loc_acc = loc_acc[loc_acc != 'Local DRAM accesses'].astype('int')
    remote_acc = remote_acc[remote_acc != 'Remote DRAM accesses'].astype('int')
    acc = np.zeros((32,),dtype=int)
    i=0
    while i < len(loc_acc):
        tmp_loc_acc = loc_acc[i:i+32]
        tmp_rem_acc = remote_acc[i:i+32]
        acc+= (tmp_loc_acc+tmp_rem_acc)
        i+=33
    return acc

def plot_npb():
    df = sns.load_dataset('NPB',data_home='/Users/zhangyuxin/Desktop/static',cache=True)
    sns.set_style(style="whitegrid")
    # Execution time improvement(%) Latency(ns) Local_ratio(%) QPI improve(%) Imblance
    ax = sns.barplot(x="App",y="QPI improve(%)",hue="Function",data = df)
    fig = ax.get_figure()
    fig.savefig("NPB/QPI improve(%)", dpi = 300)

    
#    sns.barplot(x,y,palette='Blues_d')
    
def plot_rotor():
    df = sns.load_dataset('rotor_data',data_home='/Users/zhangyuxin/Desktop/static',cache=True)
    sns.set_style(style="whitegrid")
    # Exec_time(s) QPI(MB) num_threads Latency(ns) Local_ratio(%) num_threads
    ax = sns.pointplot(x="Number of threads",y="Imbalance",hue="Function",data = df)
    fig = ax.get_figure()
    fig.savefig("Imbalance", dpi = 300)

def plot_parsec():
    df = sns.load_dataset('PARSEC',data_home='/Users/zhangyuxin/Desktop/static',cache=True)
    sns.set_style(style="whitegrid")
    # Execution time improvement(%) QPI improve(%) Latency(ns) Imblance
    ax = sns.barplot(x="APP",y="Imblance",hue="Function",data = df)
    fig = ax.get_figure()
    fig.savefig("PARSEC/Imblance", dpi = 300)

def plot_overhead():
    df = sns.load_dataset('Overhead',data_home='/Users/zhangyuxin/Desktop/static',cache=True)
    sns.set_style(style="whitegrid")
    # Execution time improvement(%) QPI improve(%) Latency(ns) Imblance
    ax = sns.barplot(x="App",y="Extra overhead",hue="Function",data = df)
    fig = ax.get_figure()
    fig.savefig("overhead/Execution time(%)", dpi = 300)

def plot_kmaf():
    df = sns.load_dataset('kMAF',data_home='/Users/zhangyuxin/Desktop/static',cache=True)
    sns.set_style(style="whitegrid")
    # Execution time improvement(%) QPI improve(%) Latency(ns) Imblance
    ax = sns.barplot(x="App",y="Execution time improvement(%)",hue="Function",data = df)
    plt.xticks(rotation=15)
    fig = ax.get_figure()
    fig.savefig("kmaf/Execution time improvement(%)", dpi = 300)
    
def norm(l):
    l = l/np.max(l)*100
    return l

def normalization(data):
    _range = float(np.max(data) - np.min(data))
    return (data - np.min(data)) / _range


if __name__ == '__main__':
#    path = "file.csv"
#    bw_0,bw_1 = read_bw(path)
#    qpi0,qpi1 =read_qpi()
#    plot_parsec()
#    read_latency()
#    plot()
    plot_rotor()
#    df = sns.load_dataset('qpi',data_home='/Users/zhangyuxin/Desktop/static',cache=True)
#    acc = read_numa()
#    plot_kmaf()
#     plot_overhead()
    
    
    