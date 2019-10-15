function powBaseExponente(base, exponente)
  return base ^ exponente
end

function powVectorExponente(v, e)
  aux = zeros(size(v)[1])
  a = 1
  for i in v
    aux[a] = i ^ e
    a += 1
  end
  return aux
end

function powArrayExponente(v, e)
  return v.^e
end

function multNumeros(a, b)
  return a*b
end

function multVectorXNumero(v, e)
  return [v...]*e
end

function multArrayXNumero(v, e)
  return v*e
end

function multMatrixXNumero(v, e)
  return v*e
end

function multMatrixXNumero(m::Array{Float64,2}, e)
  return m.*e
end

function sumArray(v)
  aux = sum(v)
end

function distEuclidea(v1, v2)
  return sqrt( sum( powArrayExponente( v1-v2, [2.0] ) ) )
end

function norma(v1, v2)
  return sqrt( sum( powArrayExponente( v1-v2, [2.0] ) ) )
end
