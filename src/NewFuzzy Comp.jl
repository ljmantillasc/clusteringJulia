path = "C:/Users/luism/OneDrive/Escritorio/parallel"
#path = pwd()
#using Distributed
include(path*"/src/lib/BasicOperation.jl")
include(path*"/src/lib/matrixOperation.jl")
include(path*"/src/parallel/ParallelOperation.jl")
include(path*"/src/secuencial/SecuencialOperation.jl")
include(path*"/src/data/ModuleEnviroment.jl")
include(path*"/src/lib/process.jl")
include(path*"/src/clusterValidity/clusterValidity.jl")

using HDF5
using JLD2
using DataFrames
using LinearAlgebra
using CSV

em = ModuleEnviroment
so = SecuencialOperation
ps = Process
cv = clusterValidity

TAM      = 10000
data     = Any
auxdata  = Any
clusters = 3
matu     = Any
mate     = Any
matc     = Any
m        = 2
dim      = Any
tamMatrix = 700


matuant  = Any
matIndexfcc = Any
matIndexfccp = Any
matIndexfcm = Any
matIndexnewPolar = Any

function prepareAndLoadData(string, cls)
  global clusters
  global matu
  global mate
  global matc
  global data
  global dim
  global auxdata
  clusters = cls
  println("hasta aqui")
  @load path*"/data/$string" auxdata
  dim = size(auxdata)
  println(dim)
  data    = [auxdata[i,:] for i=1:dim[1]]
  matu = ps.generaU(data, clusters)
  mate = ps.generaE(data)
  matc = ps.generaC(data, clusters)
end

function distMatrixAux(matrix)
  mat = zeros(dim[1], clusters)
  for j=1:dim[1]
    for i=1:clusters
        mat[j,i] = distEuclidea(data[j],matc[i])*m*sum(sum(matu[1,:].^(m-1)))
    end
  end
  return mat
end

function writeDataCSV(name, matrix)
  dt = DataFrame(matrix)
  #writetable("out/"*string(name)*".csv", dt, separator = ';', header = true)
  CSV.write("out/"*string(name), dt)
end

function distMatrix()
  global data
  global matc
  global m
  global matu
  global dim
  mat = zeros(dim[1], clusters)
  for j=1:dim[1]
    for i=1:clusters
      mat[j,i] = distEuclidea(data[j],matc[i])
    end
  end
  return [mat[j,i] for j=1:dim[1], i=1:clusters]
end

function seleccionaOrderedClusters(name)
  global dim
  global matu
  global tamMatrix
  aux = [sum(findall(matu[j] .== maximum(matu[j]))) for j=1:dim[1]]
  ModuleEnviroment.generateImage(tamMatrix,tamMatrix,aux,name)
  return aux
end

function calculaU()
  global matu
  mat = distMatrix()
  aux = zeros(dim[1], clusters)
  for k=1:dim[1]
    for i=1:clusters
      aux[k,i] = (1+sum([(mat[k,i]/mat[k,j])^(2/(m-1)) for j=1:clusters]))^(-1)
    end
  end
  for j=1:dim[1]
    matu[j] = aux[j,:]/sum(aux[j,:])
  end
end

function calculaC()
  global matc
  global data
  global matu
  global dim
  global m
  aux = so.exeSecOperationTwoARG(powVectorExponente,matu,[m for i =1:dim[1]])
  suma = sum( aux )
  for j=1:dim[1]
    aux[j] = aux[j]./suma
  end
  for i=1:clusters
    matc[i] = sum(so.exeSecOperationTwoARG(multVectorXNumero, data, [aux[j][i] for j=1:dim[1]] ))
  end
end

function calculaCParallel()
  global matc
  global data
  global matu
  global dim
  global m
  aux = zeros(dim[1], clusters)
  for j=1:dim[1]
      aux[j,:] = powVectorExponente(matu[j], m)
  end
  suma = sum([aux[i,:] for i=1:dim[1]])
  for j=1:dim[1]
    aux[j,:] = aux[j,:]./suma
  end
  mat = zeros(dim[1], clusters, dim[2])
  #@async begin
    for i=1:clusters
      for j=1:dim[1]
        mat[j,i,:] = multVectorXNumero(data[j], aux[j,i])
      end
      matc[i] = sum([ mat[j,i,:] for j=1:dim[1]])
    end
  #end
end

function fcmAlgorithm(iter, error, name)
  matuant = deepcopy(matu)
  for i=1:iter
    calculaCParallel()
    calculaU()
    #aux = maximum([maximum(abs.(matuant[j]-matu[j])) for j=1:dim[1]])
    aux = sum( sum( [ abs.(matuant[j]-matu[j]) ./ matuant[j]  for j=1:dim[1] ] ) / (dim[1] * clusters) )
    println(aux)
    if aux < error && i > 8
      break
    else
      matuant = deepcopy(matu)
    end
  end
  seleccionaOrderedClusters("$name fcm")
end

function calculaAngulos()
  mat = zeros(dim[1],clusters)
  for j=1:dim[1]
    aux = sqrt(sum(data[j].^2))
    for i=1:clusters
      aux2 = sqrt(sum(matc[i].^2))
      mat[j,i] = acos( sum(data[j].*matc[i]) / (aux * aux2) )
    end
  end
  aux = mat .- minimum(mat)
  aux = (aux/maximum(aux))*360
  e = 2.71828
  mat1 = zeros(dim[1],clusters)
  for j=1:dim[1]
    for i=1:clusters
      mat1[j,i] = 2*e^( (aux[j,i]/360)*2pi )
    end
  end
  aux = mat1 .- minimum(mat1)
  aux = (aux / maximum(aux))*256

  return aux
end

function calculaNewU()
  global matu
  mat = distMatrix()
  matri = 1 ./(calculaAngulos().+1)
  aux = zeros(dim[1], clusters)
  for k=1:dim[1]
    for i=1:clusters
      aux[k,i] = 1/(1+(matri[k,i]*mat[k,i]*mat[k,i]))
    end
  end
  for j=1:dim[1]
    matu[j] = aux[j,:]/sum(aux[j,:])
  end
end

function fcmAlgorithmNewPolar(iter, error, name)
  matuant = deepcopy(matu)
  for i=1:iter
    calculaCParallel()
    calculaNewU()
    #aux = maximum([maximum(abs.(matuant[j]-matu[j])) for j=1:dim[1]])
    aux = sum( sum( [ abs.(matuant[j]-matu[j]) ./ matuant[j]  for j=1:dim[1] ] ) / (dim[1] * clusters) )
    println(aux)
    if aux < error
      break
    else
      matuant = deepcopy(matu)
    end
  end
  seleccionaOrderedClusters("$name polar")
end

function pruebasEspecificasE(cls)
    global matIndexnfcc
    global matIndexfcm
    global matIndexnewPolar
    global matu
    global matc
    global clusters
    global m
    matIndexnfcc = zeros(cls,6)
    matIndexfcc = zeros(cls,6)
    matIndexnewPolar = zeros(cls,6)
    #prepareAndLoadData(3200, 2900, 3500, 3200, "C:/imagen.tif")
    #save("imagen1 auxdata.jld","auxdata",auxdata)
    #for ml=20:36
    for ml in [112]
      #for l=2:cls
      for l in [3,4,5,6,7]
          #prepareAndLoadData("dataimage$ml.jld2",l)
          matuCopy = deepcopy(matu)
          #fccAlgorithm(70,0.02,ml*100+l)
          #matIndexnfcc[l,:] = cv.testAllIndexA(data, matu, matc, clusters, m, dim)

          prepareAndLoadData("dataimage$ml.jld2",l)
          matuCopy = deepcopy(matu)
          #matu = deepcopy(matuCopy)
          fcmAlgorithm(70,0.001,ml*100+l)
          matIndexfcc[l,:] = cv.testAllIndexA(data, matu, matc, clusters, m, dim)

          prepareAndLoadData("dataimage$ml.jld2",l)
          matu = deepcopy(matuCopy)
          fcmAlgorithmNewPolar(70,0.001,ml*100+l)
          matIndexnewPolar[l,:] = cv.testAllIndexA(data, matu, matc, clusters, m, dim)

          #writeDataCSV("nfcc$ml.csv",matIndexnfcc )
          writeDataCSV("fcm$ml.csv",matIndexfcc )
          writeDataCSV("polar$ml.csv",matIndexnewPolar )
      end
    end
end

function generaDatos()
  count1 = 1
  x = 2000
  y = 2000
  for m=1:10
    for l=1:10
      prepareAndLoadData(x,y,x+700,y+700,"c:/imagen.tif")
      #println("$x $y")
      y = y + 100
      @save "dataimage$count1.jld2" auxdata
      count1 = count1 + 1
    end
    x = x + 100
    y = 2000
  end
end

function prepareAndLoadData(x1, y1, x2, y2, path)
  global clusters
  global matu
  global mate
  global matc
  global data
  global dim
  global auxdata
  auxdata = em.openDataCustomImage(x1,y1,x2,y2,path)
  #dim = size(auxdata)
  #data    = [auxdata[i,:] for i=1:dim[1]]
  #matu = ps.generaU(data, clusters)
  #mate = ps.generaE(data)
  #matc = ps.generaC(data, clusters)
end

function cargatabla(nombre)
  final = zeros(7,6)
  for k=2:14
    aux2 = CSV.read("out/indices/$nombre$k.csv" ;  delim=',')
    aux = [ aux2[i,j] for i=1:7,j=1:6]
    final = final + aux
  end
  return final/13
end

function generaIMG()
  for i=23:23
    for j=2:8
      cadena = pwd()*"\\$i"*"0"*"$j fcm.jld2"
      @load cadena aux
      ModuleEnviroment.generateImage(tamMatrix,tamMatrix,aux,"$i"*"0"*"$j fcm")
      @load pwd()*"\\$i"*"0"*"$j polar.jld2" aux
      ModuleEnviroment.generateImage(tamMatrix,tamMatrix,aux,"$i"*"0"*"$j polar")
    end
  end
end
