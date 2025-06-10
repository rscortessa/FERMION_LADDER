using LinearAlgebra
using Statistics
using StatsBase
BLAS.set_num_threads(4)
using ITensors, ITensorMPS
include("fermions.jl")

using .Models
let
  # Define number of spins and create spin-one indices
  N = parse(Int,ARGS[1])
  NS = parse(Int,ARGS[2]) 
  V1s = [0.0,0.5,2.0]
  t1s = [0.1,1.0]
  t2s = [0.1*i for i in 1:10]
 
 
  filename_0 = "GAPM5L" * ARGS[1] * "NS" * ARGS[2] * "MPS" * ".txt" # Output file optimization
  filename_1 = "OPTM5L" * ARGS[1] * "NS" * ARGS[2] * "MPS" * ".txt" # Output file optimization
  

  println("(N,V1,V2,t0,t1,t2,NS,NR,NC)=(",ARGS,")")
  
  
  println("Good")
  # Run DMRG to find the ground state
  nsweeps = 10
  maxdim = [16,32,64,128,256,256,256,400,400,512,1024,1024,1024,1024,1024]
  cutoff = 1E-10

  #config=random_electron_state_strings(N,4)
  
  Numbers=[N-1,N,N+1]
  En=[0.0,0.0,0.0]
  sites = siteinds("Electron",N,conserve_nf=true)
  open(filename_1,"a") do file1
      write(file1,"sweep","\t","N","\t","t1","\t","V1","\t","t2","\t","E","\t","EVAR","\n")
  end

  open(filename_0,"a") do file
      write(file,"t1","\t","V1","\t","t2","\t","GAP","\n")
      
      for t1 in t1s
      	  for V1 in V1s
      	      for t2 in t2s
	      	  
	      	  os=ladder_H(N,V1,V1,1.0,t1,t2,false)
  		  
  		  # Initialize a random MPS
  		  H = MPO(os, sites)
		  
		  for m in 1:length(Numbers)
      		      config = fill("Emp", N)
		      if Numbers[m]<=N
		         config[1:Numbers[m]] .= "Up"  # Fill first Numbers[m] sites with electrons (adjust as needed)
		      else
			 config[1:Numbers[m]-Numbers[m]%N] .= "Up"
	   		 tt=Numbers[m]%N
	   		 config[1:tt] .= "UpDn"
		      end
		      psi0 = randomMPS(sites,config; linkdims=10)  # linkdims can be adjusted
	
		      for kk in 1:nsweeps
  
			  aux = maxdim[kk] * (kk <= length(maxdim) ? 1 : 0) + maxdim[end]*(kk > length(maxdim) ? 1 : 0)
			  n=1
			  energy,psi = dmrg(H,psi0; nsweeps=n , maxdim=[aux], cutoff)
  			  psi0=psi
  			  H2 = real(inner(H,psi,H,psi))
			  En[m]=energy
			  Evar=H2-energy^2
			  open(filename_1,"a") do file1
      			      write(file1,"$kk","\t","$m","\t","$t1","\t","$V1","\t","$t2","\t","$energy","\t","$Evar","\n")
			  end
		      end
		      a=totalqn(psi0) 
  	          end
  		  GAP=2*En[2]-En[3]-En[1]
		  println("$t1","\t","$V1","\t","$t2","\t","$GAP","\n")
		  write(file,"$t1","\t","$V1","\t","$t2","\t","$GAP","\n")
              end    	     
	  end
      end
      end
end
