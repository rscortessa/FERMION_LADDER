#!/usr/bin/env python
# coding: utf-8

# In[52]:


import netket as nk
from netket.operator.fermion import destroy as c
from netket.operator.fermion import create as cdag
from netket.operator.fermion import number as nc
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigvalsh
from scipy.sparse import csr_matrix
import numpy as np
from scipy.stats import log

# In[150]:


def S_part(sparse_eig_vec,N_sites,part,eps):
    
    x=2**(int(N_sites/part))
    rho=csr_matrix(np.zeros((x,x)),dtype=complex)
    N_iter=2**(int((part-1)*N_sites/part))

    for i in range(N_iter):
        aux=sparse_eig_vec[0+x*i:x*(i+1)]
        #print(aux.shape)
        aux_2=aux.getH()
        #print(aux_2.shape)
        B=aux@aux_2
        #print(B.shape)
        rho+=B
    
    
    #s_values,s_vec=eigsh(rho,k=x-2)
    s_values=eigvalsh(rho.toarray())
    renorm_s_values=s_values[np.abs(s_values)>eps]
    #print(s_values)
    #print(np.sum(s_values))
    
    SL_2=-1.0*renorm_s_values@np.log(renorm_s_values)
    

    return SL_2


def H_ladder_G1(N,t0,t1,t2,v1,v2,d):

    PBC=True
    if PBC==False:
        V1_edges=[(i,i+N,0) for i in range(N)]
        V2_edges=[(i+1,i+N,1) for i in range(N-1)]
        to_edges=[(i,i+1,2) for i in range(N-1)]+[(i,i+1,2) for i in range(N,2*N-1)]
        edges=V1_edges+V2_edges+to_edges
    else:
        V1_edges=[(i,i+N,0) for i in range(N)]
        V2_edges=[((i+1)%N,i+N,1) for i in range(N)]
        to_edges=[(i,(i+1)%N,2) for i in range(N)]+[(i+N,(i+1)%N+N,2) for i in range(N)]
        edges=V1_edges+V2_edges+to_edges

    #print("V1")
    #print(V1_edges)
    #print("V2")
    #print(V2_edges)
    #print("to")
    #print(to_edges)
        
    
    f=np.exp(1j*np.pi*d)
    
    g =nk.graph.Graph(edges=edges)
    N_sites=g.n_nodes
    colored_edges = list(zip(g.edges(), g.edge_colors))
    
    hi = nk.hilbert.SpinOrbitalFermions(N_sites, s=None)

    H=0.0
    for ((i,k),a) in colored_edges:
        if a==0:
            H+=-t1* (cdag(hi,i,dtype=complex)*c(hi,k)+cdag(hi,k)*c(hi,i))
            H+=v1* (nc(hi,i)*nc(hi,k))
        elif a==1:
            H+=-t2* (cdag(hi,i)*c(hi,k)+cdag(hi,k)*c(hi,i))
            H+=v2* (nc(hi,i)*nc(hi,k))
        elif a==2:
            H+=-t0* (cdag(hi,i)*c(hi,k)*f+cdag(hi,k)*c(hi,i)*np.conjugate(f))
    return H

def H_ladder_G2(N,t0,t1,t2,v1,v2,d,n_fermions=None):

    PBC=False
    if PBC==False:
        V1_edges=[(i,i+1,0) for i in range(0,2*N,2)]
        V2_edges=[(i+1,i+2,1) for i in range(0,2*(N-1),2)]
        to_edges=[(i,i+2,2) for i in range(0,2*(N-1),2)]+[(i,i+2,2) for i in range(1,2*N-1,2)]
        edges=V1_edges+V2_edges+to_edges
    else:
        V1_edges=[(i,i+1,0) for i in range(0,2*N,2)]
        V2_edges=[(i+1,(i+2)%(2*N),1) for i in range(0,2*N,2)]
        to_edges=[(i,(i+2)%(2*N),2) for i in range(0,2*N,2)]+[(i,(i+2)%(2*N),2) for i in range(1,2*N+1,2)]
        edges=V1_edges+V2_edges+to_edges

    #    print("V1")
    #    print(V1_edges)
    #    print("V2")
    #    print(V2_edges)
    #    print("to")
    #    print(to_edges)
        
    
    f=np.exp(1j*np.pi*d)
    
    g =nk.graph.Graph(edges=edges)
    N_sites=g.n_nodes
    colored_edges = list(zip(g.edges(), g.edge_colors))
    
    hi = nk.hilbert.SpinOrbitalFermions(N_sites, s=None,n_fermions=n_fermions)

    H=0.0
    for ((i,k),a) in colored_edges:
        if a==0:
            H+=-t1* (cdag(hi,i,dtype=complex)*c(hi,k)+cdag(hi,k)*c(hi,i))
            H+=v1* (nc(hi,i)*nc(hi,k))
        elif a==1:
            H+=-t2* (cdag(hi,i)*c(hi,k)+cdag(hi,k)*c(hi,i))
            H+=v2* (nc(hi,i)*nc(hi,k))
        elif a==2:
            H+=-t0* (cdag(hi,i)*c(hi,k)*f+cdag(hi,k)*c(hi,i)*np.conjugate(f))
    return H


 
def GS_WF(H,eps):
    #eig_vals,eig_vec=eigsh(H.to_sparse(),k=1,which="SA")
    
    eig_vals,eig_vec=eigsh(H.to_sparse())
    #print(H.to_sparse())
    #print(eig_vals)
    #print(len(eig_vec[0]))
    
    eig_vals,eig_vec=eigsh(H.to_sparse(),k=1,which="SA")
    print(eig_vals)
    
    sparse_eig_vec = csr_matrix(np.where(np.abs(eig_vec) >= eps, eig_vec, 0.0))
    
    return eig_vals,sparse_eig_vec

