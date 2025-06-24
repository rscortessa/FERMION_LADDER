#!/usr/bin/env python
# coding: utf-8

# In[4]:


import netket as nk
from netket.operator.fermion import destroy as c
from netket.operator.fermion import create as cdag
from netket.operator.fermion import number as nc
from scipy.sparse.linalg import eigsh
import numpy as np
import sys
import os
import FERMION_APPROACH as FA
import Methods.class_WF as class_WF
from scipy.sparse import csr_matrix, save_npz, load_npz

# In[5]:


parameters_st=sys.argv
n_par=len(parameters_st)
parameters=[int(parameters_st[x]) for x in range(1,n_par)]
DG=0.01

N=parameters[0]
t0=1.0
t1=parameters[1]*DG
d=parameters[2]*DG
t2=1.0
N_sites=2*N
eps=10**(-16)
t2=1.0
v1s=[0.1*i for i in range(11)]
geo="G2"

if geo=="G1":
    Ham=FA.H_ladder_G1
elif geo=="G2":
    Ham=FA.H_ladder_G2
    
    
# In[ ]:


MASTER_DIR="N_"+str(parameters[0])+"T1_"+str(parameters[1])+"D_"+str(parameters[2])+"OBC"
try:
    os.mkdir(MASTER_DIR)
except:
    print("The folder",(MASTER_DIR,"already exists. All the files will be stored there"))
         


# In[ ]:


var=[N,parameters[1],parameters[2]]
name_var=[MASTER_DIR+"/OBSN_","T1_","D_"]
pubvar=class_WF.publisher(name_var,var,["T2","SL2","SL4","C","E"+geo])
pubvar.create()


# In[2]:

print("WORKING WITH:")
print("N=",N,"t_0=",t0,"t_1=",t1,"t_2",t2,"f=",d)
ii=0
for ii in range(len(v1s)):
    print("v1=",v1s[ii])
    H=Ham(N,t0,t1,t2,v1s[ii],v1s[ii],d)
    N_op=FA.N_fermion(N)
    N_op=N_op.to_sparse()
    E,sparse_eig_vec=FA.GS_WF(H,eps)
    aux=sparse_eig_vec[:,0].conjugate()
    n=aux.T@N_op@sparse_eig_vec[:,0]
    print(v1s[ii],n)
    save_npz(MASTER_DIR+"/"+"N_"+str(parameters[0])+"T1_"+str(parameters[1])+"V1_"+"D_"+str(parameters[2])+"v2_"+str(ii)+geo+".npz",sparse_eig_vec)
    pubvar.write([0.0,v1s[ii],0.0,n])
    ii+=1
pubvar.close()


# In[3]:





# In[134]:




