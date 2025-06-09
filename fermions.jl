
module Models

export sample_mps_to_file, ladder_H, entanglement_entropy,entanglement_mps_to_file,random_electron_state_strings

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


function ladder_H(N::Int,V1::Float64,V2::Float64,t0::Float64,t1::Float64,t2::Float64,PBC::Bool)

  if PBC
     N_aux=N
  else
     N_aux=N-1
  end
  
  os = OpSum()
  f=exp(im*0.5*pi)
  fc=exp(-im*0.5*pi)
  
  for j = 1:N_aux
  
     os += -t0*f,"Adagup",j,"Aup",j%N+1
     os += -t0*fc,"Adagup",j%N+1,"Aup",j

     os += -t0*f,"Adagdn",j,"Adn",j%N+1
     os += -t0*fc,"Adagdn",j%N+1,"Adn",j

     os += -t2,"Adagup",j%N+1,"Adn",j
     os += -t2,"Adagdn",j,"Aup",j%N+1

     os +=V2,"Nup",j%N+1,"Ndn",j
     
  end

  for j = 1:N
  
     os += -t1,"Adagup",j,"Adn",j
     os += -t1,"Adagdn",j,"Aup",j
     os += V1,"Nupdn",j

  end

  return os
end



function entanglement_entropy(psi::MPS,b::Int)

  psi = orthogonalize(psi, b)
  U,S,V = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))
  SvN = 0.0
  
  for n = 1:dim(S, 1)

     p = S[n,n]^2
     SvN -= p * log(p)

  end
  return SvN
end


function entanglement_mps_to_file(psi::MPS, filename::String, N::Int)
    # Open a file for writing
    B=Int(N/2)
    # Write to a file
    open(filename, "w") do file
        for jj in 1:B
            # Entanglement Entropy
	    S=entanglement_entropy(psi,jj)
	    println(file,jj,"\t",S)  # Write to file
        end
    end
end


function random_electron_state_strings(m::Int, num_electrons::Int)
    @assert num_electrons ≤ m "Number of electrons can't exceed number of sites"

    # Start with all empty
    state_strings = fill("Emp", m)

    # Choose `num_electrons` distinct random positions
    occupied_sites =  StatsBase.sample(1:m, num_electrons; replace=false)

    # Assign "Up" or "Dn" randomly to those sites
    for site in occupied_sites
        state_strings[site] = rand(Bool) ? "Up" : "Dn"
    end

    return state_strings
end

end

