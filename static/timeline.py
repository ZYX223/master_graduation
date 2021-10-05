#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 17:06:03 2020

@author: zhangyuxin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mapping import Map,Eagermap,greedy_lb_mapping
import time
import seaborn as sns

def standardize_obj(values, mean=None, std=None, excludes=None):
    """
    Standardize a
    Args:
        values (np.ndarray): 1-D `float32` array, the KPI observations.
        mean (float): If not :obj:`None`, will use this `mean` to standardize
            `values`. If :obj:`None`, `mean` will be computed from `values`.
            Note `mean` and `std` must be both :obj:`None` or not :obj:`None`.
            (default :obj:`None`)
        std (float): If not :obj:`None`, will use this `std` to standardize
            `values`. If :obj:`None`, `std` will be computed from `values`.
            Note `mean` and `std` must be both :obj:`None` or not :obj:`None`.
            (default :obj:`None`)
        excludes (np.ndarray): Optional, 1-D `int32` or `bool` array, the
            indicators of whether each point should be excluded for computing
            `mean` and `std`. Ignored if `mean` and `std` are not :obj:`None`.
            (default :obj:`None`)

    Returns:
        np.ndarray: The standardized `values`.
        float: The computed `mean` or the given `mean`.
        float: The computed `std` or the given `std`.
    """
    values = np.asarray(values, dtype=np.float32)
    if len(values.shape) != 1:
        raise ValueError('`values` must be a 1-D array')
    if (mean is None) != (std is None):
        raise ValueError('`mean` and `std` must be both None or not None')
    if excludes is not None:
        excludes = np.asarray(excludes, dtype=np.bool)
        if excludes.shape != values.shape:
            raise ValueError('The shape of `excludes` does not agree with '
                             'the shape of `values` ({} vs {})'.
                             format(excludes.shape, values.shape))

    if mean is None:
        if excludes is not None:
            val = values[np.logical_not(excludes)]
        else:
            val = values
        mean = val.mean()
        std = val.std()

    return (values - mean) / std, mean, std


def smoothing_extreme_values(values):
    """
    In general,the ratio of anomaly points in a time series is less than 5%[1].
    As such,simply remove the top5% data which deviate the most from the mean 
    value,and use linear interpolation to fill them.
    
    Args:
        values(np.ndarray) is a time series which has been preprosessed by linear 
        interpolation and standardization(to have zero mean and unit variance)
    
    Returns:
        np.ndarray: The smoothed `values`
    """
    
    values = np.asarray(values, np.float32)
    if len(values.shape) != 1:
        raise ValueError('`values` must be a 1-D array')
#    if (values.mean() != 0) or (values.std() != 1):
#        raise ValueError('`values` must be standardized to have zero mean and unit variance')
    
    #get the deviation of each point from zero mean
    values_deviation = np.abs(values)
    
    #the abnormal portion
    abnormal_portion = 0.05
    
    #replace the abnormal points with linear interpolation
    abnormal_max = np.max(values_deviation)
    abnormal_index = np.argwhere(values_deviation >= abnormal_max * (1-abnormal_portion))
    abnormal = abnormal_index.reshape(len(abnormal_index))
    normal_index = np.argwhere(values_deviation < abnormal_max * (1-abnormal_portion))
    normal = normal_index.reshape(len(normal_index))
    normal_values = values[normal]
    abnormal_values = np.interp(abnormal,normal,normal_values)
    values[abnormal] = abnormal_values
    
    
    return values


def extract_baseline(values,w):
    """
    A simple but effective method for removing noises if to apply moving 
    average with a small sliding window(`w`) on the KPI(`values`),separating 
    its curve into two parts:baseline and residuals.
    For a KPI,T,with a sliding window of length of `w`,stride = 1,for each 
    point x(t),the corresponding point on the baseline,denoted as x(t)*,is the 
    mean of vector (x(t-w+1),...,x(t)).
    Then the diffrence between x(t) and x(t)* is called a residuals.
    
    Args:
        values(np.ndarray): time series after preprocessing and smoothed
        
    Returns:
        tuple(np.ndarray,np.float32,np.float32):
            np.ndarray: the baseline of rawdata;
            np.float32: the mean of input values after moving average;
            np.float32: the std of input values after moving average.
        np.ndarray:the residuals between rawdata between baseline
        
        
    """
    #moving average to get the baseline
    baseline = np.convolve(values,np.ones((w,))/w,mode='valid')
    #get the residuals,the difference between raw series and baseline
    residuals = values[w-1:] - baseline
    
    return baseline,residuals

def slide_window(y):
    low_ind = int(len(y)*0.05)
    low_value = np.mean(np.sort(y)[:low_ind])
    ind = np.argwhere(y <=low_value)
    if len(ind) == 0:
        return np.array([])
    ls = [];ls.append(ind[0])
    cycle = 100
    for i,ele in enumerate(ind[1:]):
        if len(y)-1 - ele > cycle and ele-ls[-1] >= cycle:
            ls.append(ele)
    return np.array(ls)


def group_max(ls):
    ls = sorted(ls,reverse=True)
    feq_list = [0]*num
    cur_tids = set() #set                
    for tup in ls:
        if len(cur_tids) == num:
            break
        tar_group=tup[1]
        # 回填
        for tid in tar_group:
            if tid in cur_tids:
                continue
            feq_list[tid]+=1
        cur_tids = cur_tids | set(tar_group)
    return feq_list

def group_weight(ls):
    feq_list = np.zeros((num,1),dtype=float)
    feq_ls = []
    for tup in ls:
        tar_group=tup[1]
        weight = tup[0]
#        weight = 1
        tmp_feq = np.zeros((num,1),dtype=float)
        for tid in tar_group:
            tmp_feq[tid]+=1
#            feq_list[tid]+=1*weight
        tmp_feq = tmp_feq*weight
        feq_ls.append(tmp_feq)
        feq_list+=tmp_feq
    return feq_list,feq_ls

def generate_feqlist(ls):
    group_list = [];bw_list = []
    time_list = list(df['timestamp'])[start:end]
    pre_ind = 0
    for ind in ls:
        ind = int(ind)
        if ind == 0:
            continue
        tmp_list = time_list[pre_ind:ind+1]
        group = [];group_acc = 0
        for time in tmp_list:
            tlist = list(data[data['timestamp']==time]['thread'])
            group_acc += int(df[df['timestamp']==time]['access_count'])
            group+=tlist
        bw = group_acc / len(tmp_list)
        bw_list.append(bw)
        group_list.append(group)
        pre_ind = ind
    # last group
    tmp_list = time_list[ind:-1]
    group = [];group_acc = 0
    for time in tmp_list:
        tlist = list(data[data['timestamp']==time]['thread'])
        group_acc += int(df[df['timestamp']==time]['access_count'])
        group+=tlist
    bw = group_acc / len(tmp_list)
    bw_list.append(bw)
#    plt.plot(range(len(bw_list)),bw_list)
#    plt.show()
    group_list.append(group)
    tar_ = list(zip(bw_list,group_list))
    feq_list,feq_ls = group_weight(tar_)
    return feq_list,feq_ls

def get_cycles_to_ns_scale(tsc_frequency_khz):
    return int((1000000 * 2**10) / tsc_frequency_khz)

def cycles_to_nsec(cycles,tsc_frequency_khz):
    scale_factor = get_cycles_to_ns_scale(tsc_frequency_khz)
    return cycles * scale_factor *2**(-10)

def csv2mat(path):
    mat = pd.read_csv(path,header=None)
    
    mat = mat.iloc[::-1]
    mat = mat.reset_index()
    mat = mat.drop(['index'],axis=1)
    if len(mat)>32:
        mat = mat.drop([0])
        mat = mat.drop(columns = [0])
        mat = mat.reset_index()
        
        mat = mat.drop(['index'],axis=1)
        
    return np.array(mat)

def produce(mat,feq_list,mechine_config):
    map_ = Map(mat,feq_list,mechine_config)
    #map_.tid2node(16)
    map_.iter_anlyis()
    map_res = map_.map_res
    map_pos_ls = map_.map_pos_ls
    return map_res,map_pos_ls

def test_eagermap(mat,mechine_config):
    eagermap = Eagermap(mat,mechine_config['aritities'])
    eagermap.MapAlgorithm()
    map_res = eagermap.res
    return map_res

def lb_mapping(mat,feq_list,mechine_config):
    map_ = Map(mat,feq_list,mechine_config)
    map_pos = map_.lb_mapping()
    map_res = map_.pre_level
    eager_r,eager_acc,eager_pos = map_.eager_test()
    return map_res,map_pos,eager_pos

def write_pos(map_pos,path):
    f = open(path,'w+')
#    f.truncate()
    for pu in map_pos:
        if map_pos.index(pu) == len(map_pos)-1:
            f.write(str(pu))
            continue
        f.write(str(pu)+',')
    f.close()
    return

def deal_file(data):
    data = data[data['thread'] > 0]
    data['thread']-=1
    data = data.reset_index()
    return data

def deal_file2(data):
    data[data['thread'] > 0]-=2
    data = data.reset_index()
    return data

def lb_greedy(mat,feq_list,mechine_config):
    map_ = greedy_lb_mapping(mat,feq_list,mechine_config)
    map_pos = map_.MapAlgorithm()
    map_res = map_.cur_tids
    eager_ = Map(mat,feq_list,mechine_config)
    eager_r,eager_acc,eager_pos = eager_.eager_test()
    return map_res,map_pos,eager_pos


def svd():
    matrix = list()
    time_list = list(df['timestamp'])[start:end]
    for tsc in time_list:
        tlist = list(data[data['timestamp']==tsc]['thread'])
        cnt_list = [0]*num
        for tid in tlist:
            cnt_list[tid]+=1
        matrix.append(cnt_list)
    matrix= np.array(matrix)
    U,sigma,VT = np.linalg.svd(matrix)
    for i in range(U.shape[0]):
        flag = np.sign(U[:,i][abs(U[:,i])==abs(U[:,i]).max()][0])
        U[:,i] = flag*U[:,i]
        VT[i,:] = flag*VT[i,:]
    
    nsigma = (sigma/sum(sigma)).round(3)
    plt.plot(range(0,len(nsigma)+1),[1]+list(nsigma),'k-',range(0,len(nsigma)+1),[1]+list(nsigma),'k.')
    plt.title('Scree plot for OD flows',size = 15)
    plt.xlabel('Singular Values',size = 15)
    plt.ylabel('Magnitude',size = 15)
    plt.xlim(-2,30)
    plt.show()
    
    sns.set_style('whitegrid',{"xtick.major.size": 10 , "ytick.major.size": 10})
    a = 10
    b = 3
    fig, axs = plt.subplots(6, 1, figsize=(8, 30), sharey=True)
    
    for i in range(6):
        plt.sca(axs.reshape(1,6)[0][i])
        plt.plot(range(1,len(nsigma)+1),U[:,i],'k-',range(1,len(nsigma)+1),U[:,i],'k.')
    
        plt.title('Temporal flow '+str(i+1)+', normalised singular value:'+str(nsigma[i]))
        #plt.legend(frameon = 1)
        plt.ylim(-0.4,0.6)

    plt.show()
    
def norm(x,y):
    _range = max(x,y) - min(x,y)
    return (x- min(x,y)) / _range ,(y- min(x,y)) / _range

if __name__ == '__main__':
    num=32
    # load file
    title = ['timestamp','thread','accesslen']
    img_name = "a.out-18826"
    file_name = img_name+"-time_line.csv"
    path = "timelinefile/"+file_name
    data = pd.read_csv(path,header=None,names=title)
    data = deal_file(data)
#    data = deal_file2(data)
    # produce data
    cpu_feq=2100000                                          # cpu feq(khz)
    time_slice=1000000
    data['timestamp'] = data['timestamp'].astype("str")
    length = len(data['timestamp'][0])
    data = data[data['timestamp'].str.len() == length]
    data['timestamp'] = data['timestamp'].astype("float64")
    data['timestamp'] -= data['timestamp'][0]
    data = data.sort_values(by='timestamp')
    data['timestamp'] = cycles_to_nsec(data['timestamp'],cpu_feq)
    data['timestamp'] = data['timestamp'].astype("int")
    data['timestamp'] //= time_slice
    
    # timeline
    df = data.groupby(['timestamp']).count()
    df=df.reset_index()
    df['access_count'] = df['accesslen']
    df = df.drop(['thread','accesslen'],axis=1)
    # plot
    x = np.array(df['timestamp']).reshape(-1,1);x-=x[0]
    y = np.array(df['access_count'])
    fig = plt.figure();ax = fig.add_subplot(1,1,1)
    start = 0
    end = len(df)
    x=x[start:end]
    y=y[start:end]
    w = 1
    y = smoothing_extreme_values(y)
    tup,res = extract_baseline(y,w)
    x=x[w-1:]
    y=tup
    ind_list = slide_window(y)
    a=list(y)
    ax.set_title("rotor35-omp")
    ax.set_xlabel('timeline')
    ax.set_ylabel('access amount')
    ax.plot(x,y,linewidth=0.6)
    # plot vlines
    r = []
    for i in ind_list:
        r.append(x[i])
#    plt.vlines(r, 0, np.max(y), 'r', '--', label='example')
    #plt.show()
    plt.savefig(img_name+'.jpg', dpi=300)
    # comm mat and access list
    start_t = time.time()
    feq_list,feq_ls = generate_feqlist(ind_list)
    
#    feq_list = pd.read_csv("lu.B.x.32.6.mem_access.csv",header=None,dtype=int)
#    feq_list = feq_list = feq_list[1]
#    feq_list = feq_list[::-1]
#    feq_list = np.array(feq_list).reshape(-1,1)
    
    end_t = time.time()
    print('generate_feqlist time cost:',(end_t-start_t),'s')
    path = "/Users/zhangyuxin/Desktop/static/mat/fluidanimate-11521-as_matrix.csv"
    mat = csv2mat(path)
#    mat = np.zeros((num,num),dtype=int)
    pu_map = {0:0, 1:16, 2:1, 3:17, 4:2, 5:18, 6:3, 7:19, 8:4, 9:20, 10:5, 11:21, 12:6, 13:22, 14:7, 15: 23, 16:8, 17:24, 
          18:9, 19:25, 20:10, 21:26, 22:11, 23:27, 24:12, 25:28, 26:13, 27:29, 28:14, 29:30, 30:15, 31:31}
    mechine_config = {'cpus':32,'nodes':2,'pu_pos':pu_map,'aritities':[1,2,8,2]}
    #map_res,map_pos_ls = produce(mat,feq_list,mechine_config)
    #eager_res = test_eagermap(mat,mechine_config)
    #map_res,map_pos,eager_pos = lb_mapping(mat,feq_list,mechine_config)
    map_res,map_pos,eager_pos = lb_greedy(mat,feq_list,mechine_config)
    write_pos(eager_pos,"eager_pos.txt")
    write_pos(map_pos,"map_pos.txt")
    end_t2 = time.time()
    print('total time cost:',(end_t2-start_t),'s')