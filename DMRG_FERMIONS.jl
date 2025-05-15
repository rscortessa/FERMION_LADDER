
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
  V1 = parse(Float64,ARGS[2])
  V2 = parse(Float64,ARGS[3])
  t0 = parse(Float64,ARGS[4])
  t1 = parse(Float64,ARGS[5])
  t2 = parse(Float64,ARGS[6])
  NS = parse(Int,ARGS[7])
  sites = siteinds("Fermion", 2*N)
  # Initialize a random MPS
  
  
  filename_0 = "DATAM5L" * ARGS[1] * "V1_" * ARGS[2]* "V2_" * ARGS[3] * "t1_" * ARGS[5] * "t2_" * ARGS[6] * "NS" * ARGS[7] * "MPS" * ".txt" # Output file optimization
  filename_1 = "OPTM5L" * ARGS[1] * "V1_" * ARGS[2]* "V2_" * ARGS[3] * "t1_" * ARGS[5] * "t2_" * ARGS[6] * "NS" * ARGS[7] * "MPS" * ".txt" # Output file optimization

  println("(N,V1,V2,t0,t1,t2,NS,NR,NC)=(",ARGS,")")
  os=ladder_H(N,V1,V2,t0,t1,t2)
  H = MPO(os, sites)
  println("Good")
  # Run DMRG to find the ground state
  nsweeps = 10
  maxdim = [16,32,64,128,256,256,256,400,400,512,1024,1024,1024,1024,1024]
  cutoff = 1E-10
  psi0 = random_mps(sites)
  
  for kk in 1:nsweeps
  
      	aux = maxdim[kk] * (kk <= length(maxdim) ? 1 : 0) + maxdim[end]*(kk > length(maxdim) ? 1 : 0)
	n=1
	energy, psi = dmrg(H, psi0; nsweeps=n , maxdim=[aux], cutoff)
  	psi0=psi
  	H2 = inner(H,psi,H,psi)
	
	open(filename_1,"a") do file
      	    write(file,"E","\t","E_var","\t","Bd","\n")
            write(file,"$energy","\t","$H2","\t","$aux","\n")
  	end
  end	
  	   

  # Call the function to sample and save to the file
  sample_mps_to_file(psi0, filename_0, NS)
  
  
  
end
