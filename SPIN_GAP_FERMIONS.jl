
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
  V1s = [0.0,0.5,2.0]
  t1s = [0.1,1.0]
  t2s = [0.1*i for i in 1:10]
  NS = parse(Int,ARGS[2])
  sites = siteinds("Electron",N)
  # Initialize a random MPS
  
  
  filename_0 = "SGAPM5L" * ARGS[1] *"NS" * ARGS[2] * "MPS" * ".txt" # Output file optimization
  filename_1 = "OPTM5L" * ARGS[1]  * "NS" * ARGS[2] * "MPS" * ".txt" # Output file optimization
  
  println("(N,V1,V2,t0,t1,t2,NS,NR,NC)=(",ARGS,")")
  
  sites = siteinds("Electron",N,conserve_nf=true)
  
  println("Good")
  
  # Run DMRG to find the ground state
  
  nsweeps = 10
  maxdim = [16,32,64,128,256,256,256,400,400,512,1024,1024,1024,1024,1024]
  cutoff = 1E-10
  
  # Generate Random MPS
  
  config = fill("Emp", N)
  config[1:N] .= "Up"  # Fill first Numbers[m] sites with electrons (adjust as needed)
  psi0 = randomMPS(sites,config; linkdims=10)  # linkdims can be adjusted
  psi1 = randomMPS(sites,config; linkdims=10)  # linkdims can be adjusted

  E0 = 0.0
  E1 = 0.0
  
  open(filename_0,"a") do file
      write(file,"t1","\t","V1","\t","t2","\t","E0","\t","E1","\t","SGAP","\n")
  end

  open(filename_1,"a") do file1
      write(file1,"sweep","\t","state","\t","t1","\t","V1","\t","t2","\t","E","\t","EVAR","\n")
  end


        
  for t1 in t1s
      for V1 in V1s
      	  for t2 in t2s
	  
	      os=ladder_H(N,V1,V1,1.0,t1,t2,false)
	      H = MPO(os, sites)
	  
	      #GROUND STATE DMRG
              for kk in 1:nsweeps
  
      	          aux = maxdim[kk] * (kk <= length(maxdim) ? 1 : 0) + maxdim[end]*(kk > length(maxdim) ? 1 : 0)
		  n=1
		  energy, psi = dmrg(H, psi0; nsweeps=n , maxdim=[aux], cutoff)
  		  psi0=psi
  		  H2 =real(inner(H,psi,H,psi))
		  E0=energy
		  Evar=H2-E0^2		  
		  open(filename_1,"a") do file1
     		       write(file1,"$kk","\t","0","\t","$t1","\t","$V1","\t","$t2","\t","$E0","\t","$Evar","\n")
                  end
	      end

	      #FIRST EXCITED STATE DMRG
              for kk in 1:nsweeps
  
      	          aux = maxdim[kk] * (kk <= length(maxdim) ? 1 : 0) + maxdim[end]*(kk > length(maxdim) ? 1 : 0)
		  n=1
		  energy, psi = dmrg(H,[psi0],psi1; nsweeps=n , maxdim=[aux], cutoff)
  		  psi1=psi
  		  H2 = real(inner(H,psi,H,psi))
		  E1=energy
		  Evar=H2-E1^2
		  open(filename_1,"a") do file1
     		       write(file1,"$kk","\t","1","\t","$t1","\t","$V1","\t","$t2","\t","$E1","\t","$Evar","\n")
                  end
	      end


	      
	      open(filename_0,"a") do file
	          EGAP=E1-E0
                  write(file,"$t1","\t","$V1","\t","$t2","\t","$E0","\t","$E1","\t","$EGAP","\n")
		  
  	      end
	      
  	     
  	  end
      end
  end
  
    
end
