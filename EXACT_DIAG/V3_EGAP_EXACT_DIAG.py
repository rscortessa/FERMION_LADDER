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


# The parameters of the model are defined as it follows:
# The energy scale
DG=0.01
# The other parameters...
N=8
t0=1.0

tps=[i*20 for i in range(6)]
tms=[i*20 for i in range(6)]
vps=[i*40 for i in range(5)]
vms=[i*40 for i in range(5)]

geo="G2"

if geo=="G1":
    Ham=FA.H_ladder_G1
elif geo=="G2":
    Ham=FA.H_ladder_G2
    
    
# In[ ]:


MASTER_DIR="GAPS_N_"+str(N)+"D_"+str(d)
try:
    os.mkdir(MASTER_DIR)
except:
    print("The folder",(MASTER_DIR,"already exists. All the files will be stored there"))
         


# In[ ]:


var=[N,d]
name_var=[MASTER_DIR+"/N_","D_"]
pubvar=class_WF.publisher(name_var,var,["V1","V2","T1","T2","SGAP","E"+geo])
pubvar.create()

# In[2]:

N_fermions=[N-1,N,N+1]

print("WORKING WITH:")
for tp in tps: 
    for tm in tms:
        for vp in vps:
            for vm in vms:
                t1=int((tp+tm)/2.0)
                t2=int((tp-tm)/2.0)
                v1=int((vp+vm)/2.0)
                v2=int((vp-vm)/2.0)
                En=np.zeros(3)
                for n_f in range(len(N_fermions)):
                    print("t1=",t1*DG,"v1=",v1*DG,"v2=",v2*DG,"t2=",t2*DG,"N=",N_fermions[n_f])
                    H=Ham(N,t0,t1*DG,t2*DG,v1*DG,v2*DG,d*DG,N_fermions[n_f])
                    E,sparse_eig_vec=FA.GS_WF(H,eps)
                    En[n_f]=E
                    save_npz(MASTER_DIR+"/"+"N_"+str(N)+"V1_"+str(v1)+"V2_"+str(v2)+"T1_"+str(t1)+"T2_"+str(t2)+"N"+str(N_fermions[n_f])+"D"+str(d)+geo+".npz",sparse_eig_vec)
                    
                DE=2*En[1]-En[2]-En[0]
                pubvar.write([0.0,v1,0.0,v2,0.0,t1,0.0,t2,0.0,DE,0.0,E[0]])
pubvar.close()


# In[3]:





# In[134]:




