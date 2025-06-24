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
t1s=[10,100]
v1s=[50,200]
d=50
N_sites=2*N
# Fault tolerance
eps=10**(-16)

#Here we define the grid based on the parameters....
Nv2=10
v2s=[i/(Nv2*1.0) for i in range(Nv2+1)]
Nt2=10
t2s=[i/(Nt2*1.0) for i in range(Nt2+1)]

geo="G2"

if geo=="G1":
    Ham=FA.H_ladder_G1
elif geo=="G2":
    Ham=FA.H_ladder_G2
    
    
# In[ ]:


MASTER_DIR="RESULTSL"+str(N)+"D"+str(d)+"PBC"
try:
    os.mkdir(MASTER_DIR)
except:
    print("The folder",(MASTER_DIR,"already exists. All the files will be stored there"))
         


# In[ ]:


var=[N,str(d)]
name_var=[MASTER_DIR+"/N_","D_"]
pubvar=class_WF.publisher(name_var,var,["V1","V2","T1","T2","SL2","SL4","C","E"+geo])
pubvar.create()


# In[2]:

print("WORKING WITH:")
for t1 in t1s: 
    for v1 in v1s:
        for v2 in v2s:
            for t2 in t2s:
                print("t1=",t1*DG,v1*DG,v2*v1*DG,t2*t1*DG)
                H=Ham(N,t0,t1*DG,t2*t1*DG,v1*DG,v2*v1*DG,d*DG)
                E,sparse_eig_vec=FA.GS_WF(H,eps)
                save_npz(MASTER_DIR+"/"+"N_"+str(N)+"V1_"+str(v1)+"V2_"+str(v2)+"T1_"+str(t1)+"T2_"+str(t2)+"D"+str(d)+geo+".npz",sparse_eig_vec)
                print("ecco li")
                SL2=FA.S_part(sparse_eig_vec,N_sites,2,eps)
                print("uno ")
                SL4=FA.S_part(sparse_eig_vec,N_sites,4,eps)
                c_charge=SL2-SL4
                print("due ")
                pubvar.write([0.0,v1,0.0,v2,0.0,t1,0.0,t2,0.0,SL2,0.0,SL4,0.0,c_charge,0.0,E[0]])
pubvar.close()


# In[3]:





# In[134]:




