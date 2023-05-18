import numpy as np
import math
from openpyxl import Workbook
import datetime
Association_md=np.loadtxt('miRNA-disease数字编号.txt',dtype=int)
Association_sm=np.loadtxt('SM-miRNA数字关联.txt',dtype=int)
A_md=np.zeros((541,383),dtype=int)
A_dm=np.zeros((383,541),dtype=int)
A_sm=np.zeros((831,541),dtype=int)
A_ms=np.zeros((541,831),dtype=int)
for i in Association_md:
    A_md[i[0]-1][i[1]-1]=1
    A_dm[i[1]-1][i[0]-1]=1
for i in Association_sm:
    A_sm[i[0]-1][i[1]-1]=1
    A_ms[i[1]-1][i[0]-1]=1
similarity_sm=np.loadtxt('SM similarity matrix.txt')
similarity_miRNA=np.loadtxt('miRNA similarity matrix.txt')
similarity1_disease=np.loadtxt('疾病语义类似性矩阵1.txt')
similarity2_disease=np.loadtxt('疾病语义类似性矩阵2.txt')
weight_disease=np.loadtxt('疾病语义类似性加权矩阵.txt',dtype=int)
SM_name=np.loadtxt('SM编号.txt',dtype=str)
miRNA_name=np.loadtxt('miRNA编号.txt',dtype=str)
w=0.1**0.5
#初始化
wb1 = Workbook()
ws1 = wb1.active
ws1.title = "global"
similarity_disease=0.5*(similarity1_disease+similarity2_disease)
F0=0
for i in range(383):
    for k in range(541):
        F0+=(A_dm[i][k])**2
γ_d=1/(F0/383)
for i in range(383):
    for j in range(383):
            if weight_disease[i][j] == 0:
                F1=(np.linalg.norm(A_dm[i]-A_dm[j]))**2
                KD_d=math.exp(-γ_d*F1)
                similarity_disease[i][j]=KD_d
#上面是计算高斯并整合

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

a=0
distance1=1
distance2=1
distance3=1
F_s=np.random.rand(831,100)
F_m=np.random.rand(541,100)
F_d=np.random.rand(383,100)
#初始化
while distance1>=10**-8 or distance2>=10**-8 or distance3>=10**-8:#迭代结束条件之一
    #SM层更新F_s
    X=A_sm.dot(F_m)+0.1*similarity_sm.dot(F_s)
    T_s=np.zeros((831,831))
    similarity_sm_sum=similarity_sm.sum(axis=1)
    for u in range(831):
        T_s[u][u]=similarity_sm_sum[u]
    Y=((1-w**2)*A_sm*(F_s.dot(F_m.T))+w**2*F_s.dot(F_m.T)).dot(F_m)+0.1*T_s.dot(F_s)+0.1*F_s
    updateF_s=F_s*np.sqrt(X/Y)#更新低秩矩阵
    #miRNA层
    X=A_md.dot(F_d)+A_ms.dot(updateF_s)+0.1*similarity_miRNA.dot(F_m)
    T_m=np.zeros((541,541))
    similarity_miRNA_sum=similarity_miRNA.sum(axis=1)
    for u in range(541):
        T_m[u][u]=similarity_miRNA_sum[u]
    Y1=((1-w**2)*A_ms*(F_m.dot(updateF_s.T))+w**2*F_m.dot(updateF_s.T)).dot(updateF_s)+0.1*T_m.dot(F_m)+0.1*F_m
    Y2=((1-w**2)*A_md*(F_m.dot(F_d.T))+w**2*F_m.dot(F_d.T)).dot(F_d)
    Y=Y1+Y2
    updateF_m=F_m*np.sqrt(X/Y)#更新低秩矩阵
    #disease层        
    X=A_dm.dot(updateF_m)+0.1*similarity_disease.dot(F_d)
    T_d=np.zeros((383,383))
    similarity_disease_sum=similarity_disease.sum(axis=1)
    for u in range(383):
        T_d[u][u]=similarity_disease_sum[u]
    Y=((1-w**2)*A_dm*(F_d.dot(updateF_m.T))+w**2*F_d.dot(updateF_m.T)).dot(updateF_m)+0.1*T_d.dot(F_d)+0.1*F_d
    updateF_d=F_d*np.sqrt(X/Y)#更新低秩矩阵
    distance1=np.linalg.norm(updateF_s-F_s)
    distance2=np.linalg.norm(updateF_m-F_m)
    distance3=np.linalg.norm(updateF_d-F_d)
    F_s=updateF_s.copy()
    F_m=updateF_m.copy()
    F_d=updateF_d.copy()
    a+=1
    if a==100: break#迭代次数达到100迭代结束
A_sm_new=F_s.dot(F_m.T)
SZ_all1=np.random.randint(1, 2, size=(831,541))
SY=np.nonzero(SZ_all1-A_sm)
scores=A_sm_new[SY]                           
scores_paixu_SY=np.argsort(-scores)
for i in scores_paixu_SY:
    ws1.append([SM_name[SY[0][i]][0],miRNA_name[SY[1][i]][0],scores[i]])

wb1.save('结果\所有候选样本的预测结果.xlsx')

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
             
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                