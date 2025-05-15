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


# In[5]:


parameters_st=sys.argv
n_par=len(parameters_st)
parameters=[int(parameters_st[x]) for x in range(1,n_par)]
DG=0.01

N=parameters[0]
t0=1.0
t1=parameters[1]*DG
v1=parameters[2]*DG
v2=v1
d=parameters[3]*DG
Nt2=10
N_sites=2*N
eps=10**(-16)
t2=[i/(Nt2*1.0) for i in range(Nt2+1)]


# In[ ]:


MASTER_DIR="N_"+str(parameters[0])+"T1_"+str(parameters[1])+"V1_"+str(parameters[2])+"D_"+str(parameters[3])
try:
    os.mkdir(MASTER_DIR)
except:
    print("The folder",(MASTER_DIR,"already exists. All the files will be stored there"))
         


# In[ ]:


var=[N,t1,v1,d]
name_var=[MASTER_DIR+"/N_","T1_","V1_","D_"]
pubvar=class_WF.publisher(name_var,var,["T2","C","E"])
pubvar.create()


# In[2]:

print("WORKING WITH:")
print("N",N,"t0",t0,"t1",t1,"t2","V1",v1,"V2",v2,"f",d)
for t_aux in t2:
    print("t2=",t_aux)
    H=FA.H_ladder(N,t0,t1,t_aux,v1,v2,d)
    E,sparse_eig_vec=FA.GS_WF(H,eps)
    c_charge=FA.S_part(sparse_eig_vec,N_sites,2,eps)-FA.S_part(sparse_eig_vec,N_sites,4,eps)
    pubvar.write([0.0,t_aux,0.0,c_charge,0.0,E[0]])
pubvar.close()


# In[3]:





# In[134]:




