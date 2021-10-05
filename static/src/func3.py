#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 23:20:32 2020

@author: zhangyuxin
"""

import pandas as pd

titie = ['thread_id','time_stamp','address','mem_level','cost','acces_type']
data = pd.read_csv("text.dat",header=None,delimiter=' ',names=titie)

data1 = pd.read_csv("text 3.dat",header=None,delimiter=' ',names=titie)
thread_list = list(data['thread_id'])
thread_list1 = list(data1['thread_id'])

data2 = pd.read_csv("callsite_dump_2.dat",header=None,delimiter=' ',names=titie)
data2 = data2[data2['mem_level']=="Local_RAM_Hit"]
#data2['time_stamp'] = data2['time_stamp'] // 100000000



#thread_list2 = list(data2['thread_id'])
#thread_pair=[(0,6),(1,12),(2,10),(3,5),(4,14),(7,26),(8,11),(9,29),(13,23),(15,18),(16,31),(17,22),(19,25),(20,21),(24,28),(27,30)]
#
#
#
#dt ={}
#for item in thread_pair:
#    for index, row in data2.iterrows():
#        if (item[0] == row['thread_id']) or (item[1] == row['thread_id']):
#            if thread_pair.index(item) in dt:
#                dt[thread_pair.index(item)] += row['cost']
#            else:
#                dt[thread_pair.index(item)] = row['cost']
#                
#dt1 ={}
#for item in thread_pair:
#    for index, row in data1.iterrows():
#        if (item[0] == row['thread_id']) or (item[1] == row['thread_id']):
#            if thread_pair.index(item) in dt1:
#                dt1[thread_pair.index(item)] += row['cost']
#            else:
#                dt1[thread_pair.index(item)] = row['cost']
    
# 统计线程对 对DRAM的访问情况
