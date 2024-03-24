
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:10:44 2021

@author: zyll1
"""
import numpy as np
import itertools, time
from scipy import stats
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
import copy
import multiprocessing
from multiprocessing import Pool
from scipy.stats import pearsonr
import pandas as pd
import random
import pingouin
import codecs
import scipy
from scipy import sparse
import scipy.sparse
import os, sys

def text_save(filename, data):
	file = open(filename,'w')
	for i in range(len(data)):
		s = str(data[i]).replace('[','').replace(']','')
		s = s.replace("'",'').replace(',','') +'\n'
		file.write(s)
	file.close()


def ind(gene,gene_name):
	xx=[gene_name.index(i) for i in gene]
	return xx


def text(file_name):
	with open (file_name,'r+') as f:
		line = f.readline()
		while line:
			yield line.split()
			line = f.readline()

def montage(file_name):
    
	a=[]
	for i in text(file_name):
		a.append(i)
	return a

def  lexicon1(file_name):
	data=montage(file_name)
	data_dic={}
	for i in range(len(data[0])):
		data_dic[data[0][i]]=[data[j][i] for j in range(1,len(data))]
	gene_name=list(data_dic.keys())
	return data_dic,gene_name


help_msg="""
usage: python CVP_inferring_causality.py -m=m -infile=infile_name -threshold=thres -outfile=outfile_name
Options and arguments:
-m: The parameter 'm' indicates the grouping for m-fold cross-validation. The default value is 3.
-threshold: The parameter 'threshold' is a correlation threshold produced by program 'pcc_threshold.py', the fault threshold value is 0.
-infile: The parameter 'infile' indicates the file name for input. The input file is a data matrix that each row represents a sample and each column represents a variable. The file 'data_example.txt' is an example file for input file.
-outfile: The parameter 'outfile' indicates the name of output file to store the causality with format '.txt'. The default value is 'CVP_result.txt'. The output file has three columns, the first column indicates the cause variable, the second column indicates the effect variable and the third column indicates the causal strength.
"""





data=pd.read_csv('LUAD_case_matrix.csv', sep = ',',header=None,index_col=0)
all_gene=list(data.index)



f=open('LUAD_network.txt')
node_center=[]
nodes_nei=[]
data_matrix=[]
while True:
    line=f.readline()
    if len(line)==0:
        break
    fields=line.split('\t')
    center_g=all_gene[int(fields[0])-1]
    nei_g=[]
    for g in range(1,len(fields)-1):
        nei=all_gene[int(fields[g])-1]
        nei_g.append(nei)
    nei_g.append(all_gene[int((fields[len(fields)-1])[0:-1])-1])
    node_center.append(center_g)
    data_matrix.append(data.iloc[int(fields[0])-1,])
    nodes_nei.append(nei_g)
f.close()


path0='./CVP-main/TCGA'
path0='.'+os.sep
path_name='data_LUAD.txt'
path=path0+path_name

data_new=np.array(data_matrix).T
data_new=np.log(data_new+1)
#pcc0_matrix=np.corrcoef(data_new.T)
#pcc_0=pcc0_matrix
#ptcc_tho=0.1
train=data_new[0:59,:]
test=data_new[59:559,:]


gene_name=node_center
nei_list=nodes_nei



def ttr(gene1,gene_1,gn,sample_test):
	if len(gene_1)>=1: 
		erro=[]
		k=0
		y_index=ind([gene1],gn)
		x_index=ind(gene_1,gn)
		reg=LinearRegression().fit(train[:,x_index],train[:,y_index])#通过训练集模拟回归模型
		erro_1=(reg.predict(sample_test[:,x_index])-sample_test[:,y_index])**2  #计算测试集的残差平方和
		#erro.append(erro_1)
		erro_avg=erro_1[0][0]
	else:
		erro_avg=0		
	return erro_avg

def regulon(gene1,genenei,sample_test):
	gene2=genenei
	end=[]
	end_1=[gene1]
	end_2=[0]
	if len(gene2)>=2:
		erro_1=ttr(gene1,gene2,gene_name,sample_test)
		##########
		gene_b=copy.deepcopy(gene2)
		gene2_len=list(range(len(gene2)))
		for j in gene2_len:
			gene_b.remove(gene2[j])
			erro_2=ttr(gene1,gene_b,gene_name,sample_test)
			erro_j=erro_1-erro_2
			if erro_j<0:
				gene_b.append(gene2[j])
				end_1.append(gene2[j])
				end_2.append(erro_j)
			else:
				erro_1=erro_2
		end.append([end_1,end_2])
	else:
		end.append([end_1,end_2])
	return end

def  pairs1(filename):
	data=montage(filename)
	gstr=[data[i][0:int(len(data[i])/2)] for i in range(len(data))]
	rdata=[data[i][int(len(data[i])/2)+1:len(data[i])] for i in range(len(data))]
	rdata= [list(map(float, rdata[i])) for i in range(len(rdata))]
	z=np.concatenate((rdata), axis=0)
	z_tho=0
	g=[]
	g0=[]
	for i in range(len(gstr)):
		gstr0=[gstr[i][0]]
		for j in range(len(rdata[i])):
			if abs(rdata[i][j])>z_tho:
				gstr0.append(gstr[i][j+1])
				gname=[gstr[i][j+1],gstr[i][0],-rdata[i][j]]
				g0.append(gname)
		if len(gstr0)>1:
			g.append(gstr0)
	return g0


import math
es=0.00000000001
def entropy(Plist,exp_list):
    if len(Plist)>=2:
        result=0
        for x in range(len(Plist)):
            result+=-(1/math.log(len(Plist)))*(Plist[x]*abs(exp_list[x]))*math.log((Plist[x]*(abs(exp_list[x])+es)))
    else:
        #print('Error loading data')
        result=0
    return result

def  re(sample_test,gene_name,predict_pairs1):
    #mn=min(abs(sample_test[0]));
    #mx=max(abs(sample_test[0]));
    #sample_test=(abs(sample_test)-mn)/(mx-mn)
    sample_test=sample_test
    sample_entropy=[]
    for i in range(len(gene_name)):
        g=gene_name[i]
        indeed_result=[]
        indeed_exp=[]
        outdeed_result=[]
        outdeed_exp=[]
        outdeed_ft=[]
        for j in range(len(predict_pairs1)):
           predict_list=predict_pairs1[j]
           if g in predict_list:
              
              g_index=predict_list.index(g)
              if g_index==0:
                  outdeed_result.append(predict_list[2])
                  outdeed_exp.append(sample_test[0][gene_name.index(predict_list[1])])

                  out=gene_name.index(predict_list[1])
                  out_ft=abs((sample_test[0][out]-np.mean(train[:,out]))/(np.std(train[:,out])+es))
                  outdeed_ft.append(out_ft)

              else:
                 indeed_result.append(predict_list[2])
                 indeed_exp.append(sample_test[0][gene_name.index(predict_list[0])])
        merge_list=outdeed_exp+indeed_exp
        if outdeed_result==[]:
             outdeed_local_entropy=0
             outdeed_ft=0
        else:
             out_list = [(x/sum(outdeed_result)) for x in outdeed_result]
             outdeed_exp = [(x/(sum(merge_list)+es)) for x in outdeed_exp]
             outdeed_local_entropy=entropy(out_list,outdeed_exp)
             outdeed_local_entropy=(len(outdeed_result)/(len(outdeed_result)+len(indeed_result)))*outdeed_local_entropy
        
        if indeed_result==[]:
           indeed_local_entropy=0
        else:
           in_list = [(x/sum(indeed_result)) for x in indeed_result]
           indeed_exp = [(x/(sum(merge_list)+es)) for x in indeed_exp]
           indeed_local_entropy=entropy(in_list,indeed_exp)
           indeed_local_entropy=(len(indeed_result)/(len(outdeed_result)+len(indeed_result)))*indeed_local_entropy
        
        sm=abs((sample_test[0][i]-np.mean(train[:,i]))/(np.std(train[:,i])+es))
        #local_entropy=(outdeed_local_entropy+indeed_local_entropy)*sm
        local_entropy=(outdeed_local_entropy*np.mean(outdeed_ft))+(indeed_local_entropy*sm)     
        sample_entropy.append(local_entropy)
    entropy_list=copy.deepcopy(sample_entropy)
    sample_entropy.sort(reverse=True)
    entropy_sum=np.mean(sample_entropy[0:int(0.05*len(gene_name))])
    
    return (entropy_list,entropy_sum)


if __name__=="__main__":
    all_sample_entropy=[]
    entropy_matrix=[]
    for ii in range(len(test)):
        sample_test=test[ii:ii+1,:]
        rt=[]
        k=0;
        for i in gene_name:
            gene_nei=nei_list[k]
            k=k+1
            #mygene = [i for i in gene_nei if i != 0]
            mygene = gene_nei
            regu=regulon(i,mygene,sample_test)
            rt.append(regu)
        text_save(path0+path_name[0][0:-4]+'_LUAD_pcc0_r.txt',rt)
        f = codecs.open(path0+path_name[0][0:-4]+'_LUAD_pcc0_r.txt', mode='r', encoding='utf-8')
        line = f.readline() 
        targets = []
        regulators=[]
        while line:
          a = line.split()
          b = a[0:1] 
          b1=a[0:int(len(a)/2)]
          targets.append(b[0])
          regulators.append(b1)
          line = f.readline()
        f.close()
        print(ii)
        filename=path0+path_name[0][0:-4]+'_LUAD_pcc0_r.txt'
        predict_pairs1=pairs1(path0+path_name[0][0:-4]+'_LUAD_pcc0_r.txt')
        (entropy_ls,single_sample_entropy)=re(sample_test,gene_name,predict_pairs1)
        entropy_matrix.append(entropy_ls)
        all_sample_entropy.append(single_sample_entropy)
        os.remove(path0+path_name[0][0:-4]+'_LUAD_pcc0_r.txt')



entropy_matrix_data=np.array(entropy_matrix).T
entropy_matrix_result=pd.DataFrame(data=entropy_matrix_data)
entropy_matrix_result.to_csv('LUAD_entropy_matrix.csv')
