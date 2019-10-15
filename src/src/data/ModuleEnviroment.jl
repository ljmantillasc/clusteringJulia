module ModuleEnviroment
  using Core
  using Base
  using PyCall
  using CSV
  @pyimport lib

  function openDataImage(path)
    #lib.openImage("C:/imagen1.tif")
    lib.openImage(path)
    size  = lib._getSize()
    dim   = size[1].*size[2]
    band  = lib._getBandCount()
    matrix = zeros(dim, band)
    #aux   = lib._getVectorBand(1)
    for i = 1: band
      aux = lib._getVectorBand(i)
      for j=1:dim
        matrix[j,i] = aux[j]
      end
    end
    return matrix
  end

  function openDataCustomImage(x1,y1,x2,y2,path)
    #tener en cuenta x1 = 0     x2 = 500 => 0,1,2,3,4,...,499
    lib.openImage(path)
    dim   = (x2-x1)*(y2-y1)
    band  = lib._getBandCount()
    matrix = zeros(dim, band)
    for i = 1: band
      aux = lib._getVectorCustomSize(x1,y1,x2,y2,i)
      for j=1:dim
        matrix[j,i] = aux[j]
      end
    end
    return matrix
  end

  function openDataCustomVectorImage(x1,y1,x2,y2,path)
    #tener en cuenta x1 = 0     x2 = 500 => 0,1,2,3,4,...,499
    lib.openImage(path)
    dim   = (x2-x1)*(y2-y1)
    band  = lib._getBandCount()
    matrix = zeros(dim, band)
    for i = 1: band
      aux = lib._getCustomVector(x1,y1,x2,y2,i)
      for j=1:dim
        matrix[j,i] = aux[j]
      end
    end
    return matrix
  end

  function generateImage(w, h, data, name)
    lib.createImage(w, h, data, name)
  end

  function loadMatrixCSV(path, numberArchives, rows, columns)
    #genera un arreglo de matrices, para ser trabajadas como elementos
    #path = "data"
    data = zeros(numberArchives, rows, columns)
    for i=1:numberArchives
      aux = CSV.read("$path/$i.csv" ;  delim=',')
      data[i,:,:] = [aux[i,j] for i=1:49,j=1:512]
    end
    return data
  end

end
