
module Models

export sample_mps_to_file, ladder_H

using ITensors, ITensorMPS

function sample_mps_to_file(psi::MPS, filename::String, N::Int)
    # Open a file for writing
    open(filename, "w") do file
        for _ in 1:N
            # Sample a configuration from the MPS
            config = sample!(psi)
	    # config.-= 3/2
	    # config.*=2
	    config.=round.(Int,(config.-3/2)*2)
            # Convert the configuration into a string
            # The configuration `config` is a vector of spins at each site
            config_str = join(string.(config), " ")  # Convert each spin to a string
            println(file, config_str)  # Write to file
        end
    end
    println("Sampling complete. Configurations saved to $filename")
end

#MODELS:


function ladder_H(N::Int,V1::Float64,V2::Float64,t0::Float64,t1::Float64,t2::Float64)
  os = OpSum()
  f=exp(im*0.5*pi)
  fc=exp(-im*0.5*pi)
  for j = 1:N-1
  
     os += -t0*f,"Cdag",j,"C",j+1
     os += -t0*fc,"C",j+1,"Cdag",j

     os += -t0*f,"Cdag",N+j,"C",N+j+1
     os += -t0*fc,"C",N+j+1,"Cdag",N+j

     os += -t1,"Cdag",j,"C",j+N
     os += -t1,"Cdag",j+N,"C",j

     os += -t2,"Cdag",j+1,"C",j+N
     os += -t2,"Cdag",j+N,"C",j+1

     os += V1,"N",j,"N",j+N
     os += V2,"N",j+1,"N",j+N
     
  end

  os += -t1,"Cdag",N,"C",2*N
  os += -t1,"Cdag",2*N,"C",N

  os += V1,"N",N,"N",2*N
  


  return os
end



end

