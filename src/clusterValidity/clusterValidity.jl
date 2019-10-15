module clusterValidity
  using LinearAlgebra

  function powArrayExponente(v, e)
    return v.^e
  end

  function multVectorXNumero(v, e)
    return [v...]*e
  end

  function distEuclidea(v1, v2)
    return sqrt( sum( powArrayExponente( v1-v2, [2.0] ) ) )
  end

  function distMatrix(data, matc, matu, dim, clusters)
    mat = zeros(dim[1], clusters)
    for j=1:dim[1]
      for i=1:clusters
        mat[j,i] = distEuclidea(data[j],matc[i])*(matu[j][i])^2
      end
    end
    return sum(mat)
  end

  function partitionCoefficient(matu, dim)
    return sum( sum([matu[j].^2 for j=1:dim[1]]) )/dim[1]
  end

  function partitionCoefficientM(matu, dim, clusters)
    return 1-(clusters/(clusters-1))*(1-partitionCoefficient(matu, dim))
  end

  function partitionEntropy(matu, dim)
    dimmatu = size(matu[1])
    return -1.0*sum( sum([matu[j][i] * log( 2 , matu[j][i]) for j=1:dim[1],i=1:dimmatu[1] ]) )/dim[1]
  end

  function fukuyamaSugeno(data, matc, matu, dim, clusters)
    mat1 = zeros(dim[1], clusters)
    mat2 = zeros(dim[1], clusters)
    for j=1:dim[1]
      for i=1:clusters
        mat1[j,i] = distEuclidea(data[j],matc[i])*(matu[j][i])^2
      end
    end
    centroida = sum(matc)/clusters
    for j=1:dim[1]
      for i=1:clusters
        mat2[j,i] = distEuclidea(matc[i],centroida)*(matu[j][i])^2
      end
    end
    return sum(mat1) - sum(mat2)
  end

  function fuzzyHypervolume(data, matc, matu, m, dim, clusters)
    return sum([sqrt(det(sum([(matu[j][i]^2)*( (data[j]-matc[i])*(data[j]-matc[i])' ) for j=1:dim[1]])/sum([matu[j][i] for j=1:dim[1]]))) for i=1:clusters])
  end

  function xieBeni(data, matc, matu, dim, m, clusters)
    aux  = [distEuclidea(data[j], matc[i]) * (matu[j][i]^m) for j=1:dim[1], i=1:clusters]
    aux2 = [distEuclidea(matc[j], matc[i]) for j=1:clusters, i=1:clusters]
    for j=1:clusters
      aux2[j,j] = maximum(aux2)
    end
    return sum(aux)/(dim[1]*(minimum(aux2)) )
  end

  function testAllIndexA(data, matu, matc, clusters, m, dim)
    return [partitionCoefficient(matu, dim),
            partitionCoefficientM(matu, dim, clusters),
            partitionEntropy(matu, dim),
            fukuyamaSugeno(data, matc, matu, dim, clusters),
            fuzzyHypervolume(data, matc, matu, m, dim, clusters),
            xieBeni(data, matc, matu, dim, m, clusters)]
  end
end
