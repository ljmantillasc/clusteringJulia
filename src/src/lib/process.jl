module Process

function calculaU(data::Vector{Array{Float64,1}} , clusters::Int64)
  dim = size(data);
  aux = rand(1.0:100.0,dim[1],clusters)
  for j=1:dim[1]
    suma = sum(aux[j,:])
    for i=1:clusters
      aux[j,i] = aux[j,i]./suma
    end
  end
  return [aux[j,:] for j=1:dim[1]]
end

function calculaU(data::Vector{Array{Float64,2}} , clusters::Int64)
  dim = size(data);
  aux = rand(1.0:100.0,dim[1],clusters)
  for j=1:dim[1]
    suma = sum(aux[j,:])
    for i=1:clusters
      aux[j,i] = aux[j,i]./suma
    end
  end
  return [aux[j,:] for j=1:dim[1]]
end

function generaU(data::Vector{Array{Float64,1}}, clusters::Int64)
  dim = size(data)
  u = calculaU(data,clusters)
  return u
end

function generaU(data::Vector{Array{Float64, 2}}, clusters::Int64)
  dim = size(data)
  u = calculaU(data,clusters)
  return u
end

function generaC(data::Vector{Array{Float64,1}}, clusters::Int64)
  dim =size(data[1])
  c = [zeros(dim[1]) for i=1:clusters ]
  return c
end

function generaC(data::Vector{Array{Float64,2}}, clusters::Int64)
  dim =size(data[1])
  c = [zeros(dim) for i=1:clusters ]
  return c
end

function generaE(data::Array{Array{Float64,1}})
  dim = size(data)
  e = zeros(dim[1])
  return e
end

function generaMatrizGC(tamX::Int64, tamY::Int64)
    return trunc.(Int64,10 * rand(tamX,tamY))
end

end
