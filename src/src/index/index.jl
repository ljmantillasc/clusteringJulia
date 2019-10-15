module Index
  #Blue Green Red Red_Edge Near_IR
  function NDVI(auxdata)
    dim = size(auxdata)
    aux = [(auxdata[j,5]-auxdata[j,3])/(auxdata[j,5]+auxdata[j,3]) for j=1:dim[1]]
    return aux
  end

  function NDWI(auxdata)
    dim = size(auxdata)
    aux =  [(auxdata[j,2]-auxdata[j,5])/(auxdata[j,2]+auxdata[j,5]) for j=1:dim[1]]
  end

  function SR(auxdata)
    dim = size(auxdata)
    aux =  [(auxdata[j,5])/(auxdata[j,3]) for j=1:dim[1]]
    return aux
  end

  function NHFD(auxdata)
    dim = size(auxdata)
    aux =  [(auxdata[j,4]-auxdata[j,1])/(auxdata[j,4]+auxdata[j,1]) for j=1:dim[1]]
    return aux
  end

  function SAVI(auxdata)
    dim = size(auxdata)
    aux =  [(auxdata[j,5]-auxdata[j,3])/(auxdata[j,4]+auxdata[j,1]+0.5)*(1.5) for j=1:dim[1]]
    return aux
  end

end
