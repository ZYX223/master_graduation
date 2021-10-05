#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 15:32:10 2021

@author: zhangyuxin
"""
import numpy as np
import matplotlib.pyplot as plt
import copy

class Greedy_paired():
    def __init__(self,mat,pair_len):
        self.mat = mat
        self.global_mat = mat
        self.pair_len = pair_len
        if self.mat.shape[0] != self.mat.shape[1]:
            print("* The communication matrix has to be a square.")
            raise ValueError
        self.length = mat.shape[0]
        self.ind2tid = {x:[x] for x in range(self.length)}
        self.res = list()
        self.paired = list()
        
    def winner(self,owner,vec):
        max_=-1;win=0
        for i in range(len(vec)):
            if self.ind2tid[i] not in self.paired:
                if i != owner:
                    if vec[i] > max_:
                        max_=vec[i]
                        win=i
        return win
    
    def save_res(self,ls1,ls2):
        self.res.append(ls1+ls2)
        self.paired.append(ls1)
        self.paired.append(ls2)
        return
    
    def greedy_pair(self):
        iter_=2
        while iter_ <= self.pair_len:
        #for iter_ in range(2,self.pair_len+2,2):
            self.paired = list()
            self.res = list()
            for i in range(self.length):
                if self.ind2tid[i] not in self.paired:
                    # index
                    win = self.winner(i,self.mat[i])
                    self.save_res(self.ind2tid[i],self.ind2tid[win])
                    #self.save_res((i,win))
            self.length//=2
            self.ind2tid = {}
            #last not exec
            if iter_ == self.pair_len:
                return self.res
            self.generate_iter_mat()
            iter_*=2
        return
    
    def bcast_comm(self,pair,cur):
        sum_=0
        for tid in pair:
            sum_+=self.global_mat[tid][cur]
        return sum_
    
    def winner2(self,pair):
        max_=-1;max_tid=-1;
        for tid in range(self.length):
            if tid not in pair and tid not in self.paired:
                tmp = self.bcast_comm(pair,tid)
                if tmp > max_:
                    max_tid = tid
                    max_=tmp
        return max_tid
    
    def generate_one_level(self,start):
        cnt = 1;pair = []
        pair.append(start)
        while cnt < self.pair_len:
            tar_tid = self.winner2(pair)
            pair.append(tar_tid)
            cnt+=1
        return pair
    
    def save_res2(self,pair):
        self.res.append(pair)
        for tid in pair:
            self.paired.append(tid)
        return
    
    def greedy_pair2(self):
        for i in range(self.length):
            if i not in self.paired:
                pair = self.generate_one_level(i)
                self.save_res2(pair)
        return
    
    def generate_iter_mat(self):
        tmp_mat = np.zeros((self.length,self.length))
        for i in range(len(self.res)):
            for j in range(len(self.res)):
                if i!=j:
                    tmp_mat[i][j] = self.cal_tup(self.res[i],self.res[j])
            self.ind2tid[i] = self.res[i]
        self.mat = tmp_mat
        return
        
    def cal_tup(self,g1,g2):
        sum_=0
        for i in range(len(g1)):
            for j in range(len(g2)):
                sum_+=self.global_mat[g1[i]][g2[j]]
        return sum_
    
class Eagermap():
    def __init__(self,mat,aritities):
        self.mat = mat
        self.global_mat = mat
        self.aritities = aritities
        if self.mat.shape[0] != self.mat.shape[1]:
            print("* The communication matrix has to be a square.")
            raise ValueError
        self.length = mat.shape[0]
        self.res = list()
        self.paired = list()
        self.ind2tid = {x:[x] for x in range(self.length)}
    
    def bcast_comm(self,pair,cur):
        sum_=0
        for i in pair:
            for j in cur:
                sum_+=self.global_mat[i][j]
        return sum_
    
    def winner2(self,pair):
        max_=-1;max_tid=-1;
        for tid in range(self.mat.shape[0]):
            if self.ind2tid[tid] not in pair and self.ind2tid[tid] not in self.paired:
                sum_=[]
                for x in pair:
                    sum_+=x
                tmp = self.bcast_comm(sum_,self.ind2tid[tid])
                if tmp > max_:
                    max_tid = tid
                    max_=tmp
        return max_tid
    
    def save_res(self,pair):
        group=[]
        for x in pair:
            group+=x
            self.paired.append(x)
        self.res.append(group)
        return
    
    def generate_one_level(self,start,pair_len):
        cnt = 1;pair = []
        pair.append(self.ind2tid[start])
        #pair=self.ind2tid[start]
        while cnt < pair_len:
            tar_tid = self.winner2(pair)
            pair.append(self.ind2tid[tar_tid])
            cnt+=1
        return pair
    
    def cal_tup(self,g1,g2):
        sum_=0
        for i in range(len(g1)):
            for j in range(len(g2)):
                sum_+=self.global_mat[g1[i]][g2[j]]
        return sum_
    
    def generate_iter_mat(self):
        tmp_mat = np.zeros((self.length,self.length))
        for i in range(len(self.res)):
            for j in range(len(self.res)):
                if i!=j:
                    tmp_mat[i][j] = self.cal_tup(self.res[i],self.res[j])
            self.ind2tid[i] = self.res[i]
        self.mat = tmp_mat
        return
    
    def MapAlgorithm(self):
        for pair_len in self.aritities[::-1]:
            for i in range(self.mat.shape[0]):
                if self.ind2tid[i] not in self.paired:
                    pair = self.generate_one_level(i,pair_len)
                    self.save_res(pair)
            self.ind2tid = {}
            self.length//=pair_len
            self.generate_iter_mat()
            if pair_len == self.mat.shape[0]:
                return
            self.paired = list()
            self.res = list()
        return

class Map():
    def __init__(self,mcomm,acc_list,machcine_config):
        self.mcomm = mcomm
        self.acc_list = acc_list
        if self.mcomm.shape[0] != self.mcomm.shape[1]:
            print("* The communication matrix has to be a square.")
            raise ValueError
        self.length = mcomm.shape[0]
        self.cpus = machcine_config['cpus']
        self.nodes = machcine_config['nodes']
        self.machcine_config = machcine_config
        self.iter_ls = []
        self.map_res = {x:[] for x in range(self.nodes)}
        self.map_res_ls = []
        self.best_map_res = {}
        self.map_pos_ls = []
        self.eager_res = []
        self.pre_level = [[x] for x in range(self.length)]
        self.grouped = []
    
    def cal_iter(self):
        per_node = self.length // self.nodes
        cnt = int(np.log2(per_node))
        self.iter_ls = [2**x for x in range(1,cnt+1)]
        return
    
    def cal_group_acc(self,res):
        group_acc_dt = {}
        for group in res:
            sum_=0
            for tid in group:
                sum_+=self.acc_list[tid][0]
            group_acc_dt[sum_] = group
        return group_acc_dt
    
    def sctter2node(self,group_acc_dt):
        sorted_acc = sorted(group_acc_dt.items(), key = lambda i: i[0],reverse=True)
        if len(sorted_acc) == self.nodes:
            for i,group in enumerate(sorted_acc):
                self.map_res[i%self.nodes].append(sorted_acc[i][1])
            return
        i = 0;j = len(sorted_acc)-1
        cnt = 0
        while i<j:
            self.map_res[cnt%self.nodes].append(sorted_acc[i][1])
            self.map_res[cnt%self.nodes].append(sorted_acc[j][1])
            i+=1;j-=1
            cnt+=1
#        for i,group in enumerate(sorted_acc):
#            self.map_res[i%self.nodes].append(group[1])
        return
    
    def tid2node(self,pair_len):
        greedy_paired = Greedy_paired(self.mcomm,pair_len)
        #res = greedy_paired.greedy_pair()
        greedy_paired.greedy_pair2()
        group_acc_dt = self.cal_group_acc(greedy_paired.res)
        self.sctter2node(group_acc_dt)
        return
    
    def cal_balance(self):
        node_acc_list = []
        for i in range(self.nodes):
            sum_=0
            for pair in self.map_res[i]:
                for tid in pair:
                    sum_+=self.acc_list[tid][0]
            node_acc_list.append(sum_)
        print("node....",node_acc_list)
        acc_balance = np.var(node_acc_list)
        romote_comm = 0
        for i in range(self.nodes):
            for j in range(i+1,self.nodes):
                romote_comm += self.cal_remote(self.map_res[i],self.map_res[j])
        
        return romote_comm,acc_balance
    
    def cal_remote(self,n1,n2):
        sum_=0
        for pair1 in n1:
            for pair2 in n2:
                for tid1 in pair1:
                    for tid2 in pair2:
                        sum_+=self.mcomm[tid1][tid2]
        return sum_
    
    def norm_(self,ls):
        ls = np.array(ls)
        return 100*ls/np.max(ls)
    
    def anlyis(self,remote_list,acc_balance):
        min_ = 2**31 ;min_ind=0
        mapping_quality = []
        for i in range(len(remote_list)):
            tar = abs(remote_list[i] - acc_balance[i])
            if tar < min_:
                min_ = tar
                min_ind = self.iter_ls[i]
                self.best_map_res = self.map_res_ls[i]
            mapping_quality.append(tar)
        print("mapping_quality",mapping_quality)
        print("the best pair len is:",min_ind)
        print("best map res:",self.best_map_res)
        return 
    
    def iter_anlyis(self):
        self.cal_iter()
        remote_list = []
        acc_balance = []
        for pair_len in self.iter_ls:
            self.map_res = {x:[] for x in range(self.nodes)}
            self.tid2node(pair_len)
            r,acc=self.cal_balance()
            remote_list.append(r)
            acc_balance.append(acc)
            self.map_res_ls.append(self.map_res)
            self.map_pos_ls.append(self.map2core())
            print("pair length:",pair_len)
            print("remote....",r)
            print(self.map_res)
            print("\n")
        
        eager_r,eager_acc = self.eager_test()
        print("\n")
        remote_list.append(eager_r),acc_balance.append(eager_acc)
        self.iter_ls.append(self.iter_ls[-1]*2)
        self.anlyis(self.norm_(remote_list),self.norm_(acc_balance))
        plt.plot(self.iter_ls,self.norm_(remote_list))
        plt.plot(self.iter_ls,self.norm_(acc_balance))
        plt.legend(["reomte commuication","access var"])
        plt.show()
        self.iter_ls.pop()
    
    def eager_test(self):
        eagermap = Eagermap(self.mcomm,self.machcine_config['aritities'])
        eagermap.MapAlgorithm()
        eager_res = eagermap.res[0]
        length = len(eager_res) // self.nodes
        node_acc_list = [];node_tid_list = []
        cnt = 0
        for node in range(self.nodes):
            # per node
            sum_=0;tmp_ls=[]
            for i in range(cnt,cnt+length):
                sum_+=self.acc_list[eager_res[i]][0]
                tmp_ls.append(eager_res[i])
            node_acc_list.append(sum_)
            node_tid_list.append(tmp_ls)
            cnt+=length
#            # node 1
#            sum_=0
#            for i in range(ind,len(eager_res)):
#                sum_+=self.acc_list[eager_res[i]][0]
#            node_acc_list.append(sum_)
        print("eager_test...........................................")
        print("node....",node_acc_list)
        min_diff = np.std([node_acc_list[0],node_acc_list[1]],ddof=1)
        print(min_diff)
        acc_balance = np.var(node_acc_list)
        romote_comm = 0
        for i in range(self.nodes):
            for j in range(i+1,self.nodes):
                for tid1 in node_tid_list[i]:
                    for tid2 in node_tid_list[j]:
                        romote_comm+=self.mcomm[tid1][tid2]
        print("remote access",romote_comm)
        eager_pos = self.map2core2(node_tid_list)    
        print(node_tid_list)
        return romote_comm,acc_balance,eager_pos
                
    def map2core(self):
        map_pos = [0]*self.length
        pu_pos = self.machcine_config['pu_pos']
        pu=0
        for i in range(self.nodes):
            for pair in self.map_res[i]:
                for tid in pair:
                    map_pos[tid]=pu_pos[pu]
                    pu+=1
        return map_pos
    
    def cal_sum_(self,group,tmp_group):
        sum_=0
        for i in group:
            for j in tmp_group:
                sum_+=self.mcomm[i][j]
        return sum_
    
    def locality_rank(self,group,cur_level_group):
        rank_dt = []
        for tmp_group in cur_level_group:
            if tmp_group == group:
                continue
            sum_ = self.cal_sum_(group,tmp_group)
            
            rank_dt.append((sum_,[group,tmp_group]))
        
        rank_ls = sorted(rank_dt, key = lambda i: i[0],reverse=True)
        return rank_ls
    
    def justify_load(self,pair):
        g1 = pair[0]
        g2 = pair[1]
        g = g1+g2
        load_threshold = sum(self.acc_list) / self.length * len(g)
        vec = self.acc_list -np.mean(self.acc_list)
        vec = vec[vec>0]
        per_node_num = len(vec) // self.nodes
        if per_node_num <1: per_node_num=1
        cost_error = np.max(vec) * per_node_num
        total_load = 0
        for tid in g:
            total_load += self.acc_list[tid]
            
        if total_load <= load_threshold:
            return True
        elif total_load - load_threshold  <= cost_error:
            return True
        else:
            return False
    
    def generate_one_group(self,group,cur_level_group,max_tids):
        res = []
        res.append(group)
        rank_ls = self.locality_rank(group,cur_level_group)
        cur_level_group.remove(group)
        length_res = 0
        for pair in rank_ls:
            if self.justify_load(pair[1]):
                for g in pair[1]:
                    length_res+=len(g)
                if length_res <= max_tids:
                    cur_level_group.remove(pair[1][1])
                    res.append(pair[1][1])
                    if length_res == max_tids:
                        break
                else:
                    for g in pair[1]:
                        length_res-=len(g)
                    continue
            else:
                if length_res < max_tids:
                    continue
                elif length_res >= max_tids:
                    break
                    
        g_ls = []
        for item in res:
            self.grouped.append(item)
            g_ls+=item
        cur_level_group.append(g_ls)
        return cur_level_group

    def generate_one_level(self,max_tids):
        cur_level_group = copy.deepcopy(self.pre_level)
        for group in self.pre_level:
            if group in self.grouped:
                continue
            cur_level_group = self.generate_one_group(group,cur_level_group,max_tids)
            
        self.pre_level = copy.deepcopy(cur_level_group)
        return
    
    def lb_sanalyis(self):
        sor = self.acc_list
        i = 0
        node0 = [];node1=[]
        
        for i in range(len(self.pre_level)):
            if i==0:
                for tid in self.pre_level[i]:
                    node0.append(sor[tid])
            if i==1:
                for tid in self.pre_level[i]:
                    node1.append(sor[tid])
        print(sum(node1),sum(node0))
        min_diff = abs(sum(node1)-sum(node0))
        print(min_diff)
        romote_comm = 0
        nodes=2
        for i in range(nodes):
            for j in range(i+1,nodes):
                romote_comm += self.cal_remote2(self.pre_level[i],self.pre_level[j])
        print('remote access',romote_comm)
        return
    
    def cal_remote2(self,n1,n2):
        sum_=0
        for tid1 in n1:
            for tid2 in n2:
                sum_+=self.mcomm[tid1][tid2]
        return sum_
    
    def map2core2(self,map_res):
        map_pos = [0]*self.length
        pu_pos = self.machcine_config['pu_pos']
        pu=0
        for i in range(self.nodes):
            for tid in map_res[i]:
                map_pos[tid]=pu_pos[pu]
                pu+=1
        return map_pos
    
    def lb_mapping(self):
        per_node = self.length // self.nodes
        iter_num = int(np.log2(per_node))
        for level in range(iter_num):
            max_tids = 2**(level+1)
            self.generate_one_level(max_tids)
        self.lb_sanalyis()
        map_pos=self.map2core2(self.pre_level)
        return map_pos
    
class greedy_lb_mapping():
    def __init__(self,mat,acc_list,machcine_config):
        self.mat = mat
        self.length = mat.shape[0]
        self.acc_list = acc_list
        self.cpus = machcine_config['cpus']
        self.nodes = machcine_config['nodes']
        self.machcine_config = machcine_config
        self.cur_tids = [[x] for x in range(self.length)]
        self.overload_tids = []
    
    def bcast_comm(self,n1,n2):
        sum_=0
        for i in n1:
            for j in n2:
                sum_+=self.mat[i][j]
        return sum_
    
    def locality_rank(self,group):
        rank_dt = []
        for tmp_tid in self.cur_tids:
            if group == tmp_tid:
                continue
            sum_ = self.bcast_comm(group,tmp_tid)
            rank_dt.append((sum_,[group,tmp_tid]))
        rank_ls = sorted(rank_dt, key = lambda i: i[0],reverse=True)
        return rank_ls
    
    def cal_acc_list(self):
        vec = self.acc_list -np.mean(self.acc_list)
        self.sp_tids = np.argwhere(vec >0)
        vec = vec[vec>0]
        self.per_node_num = len(vec) // self.nodes
        if self.per_node_num <1: 
            self.per_node_num=1
        return
    
    def justify_load(self,pair):
        g2 = pair[1]
        if g2[0] in self.sp_tids:
            self.group_sp+=1
        if self.group_sp <= self.per_node_num:
            return True
        else:
            self.group_sp-=1
            return False
    
    def cal_confidence_interval(self,p,num):
        f = np.delete(self.acc_list,p)
        f = np.sort(f)
        f = f[::-1]
        max_ = np.sum(f[:num])/num
        min_ = np.sum(f[-num:])/num
        return [min_,max_]
    
    def justify_load2(self,pair):
        p1 = pair[0]
        p2 = pair[1]
        p = p1+p2
        total_load = 0
        for tid in p:
            total_load += self.acc_list[tid]
        last_value = np.sum(self.acc_list) / self.nodes - total_load
        last_num = self.length // self.nodes - len(p)
        if last_num == 0:
            if p2 in self.overload_tids:
                return False
            else:
                return True
        ave_last_value = last_value / last_num
        interval = self.cal_confidence_interval(p,last_num)
        if ave_last_value >= interval[0] and ave_last_value <= interval[1]:
            return True
        else:
            self.overload_tids.append(p2)
            return False
        
    def winner(self,group):
        rank_ls = self.locality_rank(group)
        for pair in rank_ls:
            if self.justify_load2(pair[1]):
                return pair[1][1]
            else:
                continue
        
        return rank_ls[0][1][1]
    
    def generate_one_group(self,group_len):
        cnt = 1
        group=self.cur_tids[0]
        while cnt < group_len:
            tar = self.winner(group)
            self.cur_tids.remove(group)
            self.cur_tids.remove(tar)
            group+=tar
            self.cur_tids.append(group)
            cnt+=1
        return 
            
    def cal_remote(self,n1,n2):
        sum_=0
        for tid1 in n1:
            for tid2 in n2:
                sum_+=self.mat[tid1][tid2]
        return sum_
    
    def lb_sanalyis(self):
        sor = self.acc_list
        i = 0
        node0 = [];node1=[]
        for i in range(len(self.cur_tids)):
            if i==0:
                for tid in self.cur_tids[i]:
                    node0.append(sor[tid])
            if i==1:
                for tid in self.cur_tids[i]:
                    node1.append(sor[tid])
        print(sum(node0),sum(node1))
        min_diff = np.std([sum(node0),sum(node1)],ddof=1)
        print(min_diff)
        romote_comm = 0
        nodes=2
        for i in range(nodes):
            for j in range(i+1,nodes):
                romote_comm += self.cal_remote(self.cur_tids[i],self.cur_tids[j])
        print('remote access',romote_comm)
        return
    
    def map2core(self,map_res):
        map_pos = [0]*self.length
        pu_pos = self.machcine_config['pu_pos']
        pu=0
        for i in range(self.nodes):
            for tid in map_res[i]:
                map_pos[tid]=pu_pos[pu]
                pu+=1
        return map_pos
       
    def MapAlgorithm(self):
        self.cal_acc_list()
        per_node_tids = self.length // self.nodes
        for node in range(self.nodes -1):
            self.overload_tids.clear()
            self.generate_one_group(per_node_tids)
        
        tmp = self.cur_tids[-(self.nodes-1):]
        last_node = self.cur_tids[:-(self.nodes-1)]
        ls = []
        for tid in last_node:
            ls+=tid
        self.cur_tids = tmp
        self.cur_tids.append(ls)
        self.lb_sanalyis()
        map_pos=self.map2core(self.cur_tids)
        return map_pos
               
    