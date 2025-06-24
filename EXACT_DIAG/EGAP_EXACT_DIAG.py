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
v1=parameters[2]*DG
v2=v1
d=parameters[3]*DG
Nt2=10
N_sites=2*N
eps=10**(-16)
t2=[i/(Nt2*1.0) for i in range(Nt2+1)]
geo="G2"

if geo=="G1":
    Ham=FA.H_ladder_G1
elif geo=="G2":
    Ham=FA.H_ladder_G2
    
    
# In[ ]:

PBC="PBC"
MASTER_DIR="GAPS_N_"+str(parameters[0])+"T1_"+str(parameters[1])+"V1_"+str(parameters[2])+"D_"+str(parameters[3])+PBC
try:
    os.mkdir(MASTER_DIR)
except:
    print("The folder",(MASTER_DIR,"already exists. All the files will be stored there"))
         


# In[ ]:


var=[N,parameters[1],parameters[2],parameters[3]]
name_var=[MASTER_DIR+"/N_","T1_","V1_","D_"]
pubvar=class_WF.publisher(name_var,var,["T2","DE"+geo])
pubvar.create()


# In[2]:

print("WORKING WITH:")
print("N=",N,"t_0=",t0,"t_1=",t1,"V_1=",v1,"V_2=",v2,"f=",d)

N_fermions=[N-1,N,N+1]

ii=0

for t_aux in t2:
    En=np.zeros(3)
    for n_f in range(len(N_fermions)):
        print("t2=",t_aux)
        H=Ham(N,t0,t1,t_aux,v1,v2,d,N_fermions[n_f])
        E,sparse_eig_vec=FA.GS_WF(H,eps)
        En[n_f]=E
        save_npz(MASTER_DIR+"/"+"N_"+str(parameters[0])+"T1_"+str(parameters[1])+"V1_"+str(parameters[2])+"D_"+str(parameters[3])+"T2_"+str(ii)+"N_F"+str(N_fermions[n_f])+geo+".npz",sparse_eig_vec)
        
    DE=2*En[1]-En[2]-En[0]
    pubvar.write([0.0,t_aux,0.0,DE])
    
    ii+=1
    
pubvar.close()


# In[3]:





# In[134]:




