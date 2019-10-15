#C:\Users\luism\OneDrive\Escritorio\parallel
path = "C:/Users/luism/OneDrive/Escritorio/parallel"
using Distributed
include(path*"/src/lib/BasicOperation.jl")
include(path*"/src/lib/matrixOperation.jl")
include(path*"/src/parallel/ParallelOperation.jl")
include(path*"/src/secuencial/SecuencialOperation.jl")
include(path*"/src/data/ModuleEnviroment.jl")
include(path*"/src/lib/process.jl")
#include(path*"/src/clusterValidity/clusterValidity.jl")


#using ModuleEnviroment
#using ParallelOperation
#using SecuencialOperation
#using Process
#using clusterValidity
using HDF5
using JLD2
using DataFrames
using CSV

clusters = 3

function loadMatrixData()
    global dim
    global matrixData
    matrixData = @load "image20x20.jld" matrixData
    #matrixData = matrixData[1:6,1:6]
    dim = size(matrixData)
end

function prepareAndLoadData16bands(x1, y1, x2, y2, path)
  global clusters
  global matu
  global mate
  global matc
  global data
  global dim
  global auxdata
  auxdata = ModuleEnviroment.openDataCustomVectorImage(x1,y1,x2,y2,path)
  dim = size(auxdata)
  data = [auxdata[i,:] for i=1:dim[1]]
  matu = Process.generaU(data, clusters)
  mate = Process.generaE(data)
  matc = Process.generaC(data, clusters)
end

function prepareAndLoadData(x1, y1, x2, y2, path)
  global clusters
  global matu
  global mate
  global matc
  global data
  global dim
  global auxdata
  auxdata = ModuleEnviroment.openDataCustomImage(x1,y1,x2,y2,path)
  dim = size(auxdata)
end

function normalizaColummns(matrix)
    #matrix = matrix.*matrix
    #matrix = matrix'
    dim = size(matrix)
    for j=1:dim[1]
        suma = sum(matrix[j,:])
        for i=1:dim[2]
            matrix[j,i] = matrix[j,i]/suma
        end
    end
    #matrix = matrix'
    #dim = size(matrix)
    #for j=1:dim[1]
    #    suma = sum(matrix[j,:])
    #    for i=1:dim[2]
    #        matrix[j,i] = matrix[j,i]/suma
    #    end
    #end
end


function omp(A, v, limit)
    dim = size(A)
    omega = []
    x = []
    r = deepcopy(v)
    #process
    x = []
    xj = []
    indexMax = 0
    for i=1: limit
        xj = [sum(A[j,:].*r) for j=1:dim[1]]
        indexMax = indiceMaxVector(abs.(xj))
        maxVector = xj[indexMax[1]]
        append!(omega, indexMax)
        append!(x, maxVector)
        r = r - (maxVector*A[indexMax[1],:])
    end
    return x, omega
end

function ompAnter(A, v, limit)
    dim = size(A)
    omega = []
    x = []
    r = deepcopy(v)
    #process
    x = []
    xj = []
    for i=1: limit
        xj = [sum(A[j,:].*r) for j=1:dim[1]]
        indexMax = indiceMaxVector(abs.(xj))
        maxVector = xj[indexMax[1]]
        append!(omega, indexMax)
        append!(x, maxVector)
        #modificacion orthogonal
        dimOmega = size(omega)
        auxDistOmega = [distEuclidea( (A[omega[j],:].*x[j]), r)^2 for j=1:dimOmega[1]]
        minDist = minimum(auxDistOmega)
        indexMin = indiceMinVector(auxDistOmega)
        r = r - A[omega[indexMin],:]*x[indexMin]
    end
    return x, omega
end

function indiceMaxVector(vector)
    maxv = maximum(vector)
    index = findall(a->a==maxv, vector)
    return index
end

function indiceMinVector(vector)
    minv = minimum(vector)
    index = findall(a->a==minv, vector)[1]
    return index
end


function generaPrueba(A, limit)
    dim = size(A)
    aux = zeros(dim[1])
    for i = 1: limit
        o = omp(A, A[100*i,:], 70)
        for i in o aux[i] =1 end
        ModuleEnviroment.generateImage( 700, 700,trunc.(Int64,aux),"nuevo $i ")
    end
end


function algClustering(A, c, limit)
    dim = size(A)
    mat = zeros(dim[1],c)
    centroides = zeros(c, dim[2])
    matSparce = zeros(dim[1], c)
    lista = [i for i=1:c]
    for i=1:dim[1]
        index = rand(lista)
        mat[i, index] = 1.0
    end
    #calcular centroides
    for i=1:c
        aux = A.*mat[:,i] / sum(mat[:,i])
        centroides[i,:] = sum([aux[j ,:] for j=1:dim[1]])
    end
    # extraer vector esparza
    for l = 1: limit
        for j=1:dim[1]
            x, omega = omp(centroides, A[j,:], c)
            matSparce[j,:] = omega
        end
        for j=1:dim[1]
            mat[j,:] = zeros(c)
            mat[j,trunc.(Int64,matSparce[j,1])] = 1
        end
    end
    aux = zeros(dim[1])
    for j=1:dim[1]
        aux[j] = indiceMaxVector(mat[j,:])[1]
    end
    ModuleEnviroment.generateImage( 700, 700,trunc.(Int64, aux),"nuevo A $l ")
    #generar segmentación
    return mat, matSparce, centroides, aux
end

function algClusteringM1(A, c, limit)
    dim = size(A)
    mat = zeros(dim[1],c)
    centroides = zeros(c, dim[2])
    matSparce = zeros(dim[1], c)
    lista = [i for i=1:c]
    for i=1:dim[1]
        index = rand(lista)
        mat[i, index] = 1.0
    end
    #calcular centroides
    for i=1:c
        aux = A.*mat[:,i] / sum(mat[:,i])
        centroides[i,:] = sum([aux[j ,:] for j=1:dim[1]])
    end
    # extraer vector esparza
    for l = 1: limit
        for j=1:dim[1]
            x, omega = omp(centroides, A[j,:], c)
            matSparce[j,:] = omega
        end
        for j=1:dim[1]
            mat[j,:] = zeros(c)
            indice = trunc.(Int64, sum(matSparce[j,:])/c)
            mat[j, indice] = 1
        end
    end
    aux = zeros(dim[1])
    for j=1:dim[1]
        aux[j] = indiceMaxVector(mat[j,:])[1]
    end
    ModuleEnviroment.generateImage( 700, 700,trunc.(Int64, aux),"nuevo Aux ")
    #generar segmentación
    return mat, matSparce, centroides, aux
end

function kernelClusteringOriginal(A, c, limit)
    dim = size(A)
    mat = zeros(dim[1],c)
    centroides = zeros(c, dim[2])
    matSparce = zeros(dim[1], c)
    lista = [i for i=1:c]
    for i=1:dim[1]
        index = rand(lista)
        mat[i, index] = 1.0
    end
    #calcular centroides
    for i=1:c
        aux = A.*mat[:,i] / sum(mat[:,i])
        centroides[i,:] = sum([aux[j ,:] for j=1:dim[1]])
    end
    println("$centroides")
    # extraer vector esparza
    for l = 1: limit
        for j=1:dim[1]
            x, omega = omp(centroides, A[j,:], c)
            matSparce[j,:] = omega
        end
        for j=1:dim[1]
            mat[j,:] = zeros(c)
            indice = trunc.(Int64, sum(matSparce[j,:])/c)
            mat[j, indice] = 1
        end
    end
    return mat
end

function kernelClustering(A, c, limit)
    dim = size(A)
    mat = zeros(dim[1],c)
    centroides = zeros(c, dim[2])
    matSparce = zeros(dim[1], c)
    lista = [i for i=1:c]
    for i=1:dim[1]
        index = i%c+1
        mat[i, index] = 1.0
    end
    #calcular centroides
    for i=1:c
        aux = A.*mat[:,i] / sum(mat[:,i])
        centroides[i,:] = sum([aux[j ,:] for j=1:dim[1]])
    end
    println("$centroides")
    # extraer vector esparza
    for l = 1: limit
        for j=1:dim[1]
            x, omega = omp(centroides, A[j,:], c)
            matSparce[j,:] = omega
        end
        for j=1:dim[1]
            mat[j,:] = zeros(c)
            indice = trunc.(Int64, sum(matSparce[j,:])/c)
            mat[j, indice] = 1
        end
    end
    return mat
end


function seleccionaMaximoElementos(matrix)
    dim = size(matrix)
    aux = zeros(dim[2])
    for i = 1: dim[2]
        aux[i] = sum(matrix[:,i])
    end
    return indiceMaxVector(aux)[1]
end

function extraeSubMatriz(vectorPertenenciaDura , A)
    dim = size(A)
    numElementos = trunc.(Int64, sum(vectorPertenenciaDura))
    subMatrix = zeros(numElementos, dim[2])
    indice = 1
    for i=1:dim[1]
        if vectorPertenenciaDura[i] == 1
            subMatrix[indice,:] = A[i,:]
            indice = indice + 1
        end
    end
    return subMatrix
end

function separarClasificacion(vectorMatClass , matrixClasificacion)
    dim = size(vectorMatClass)
    matAux = zeros(dim[1],2)
    indice = 1
    for i=1:dim[1]
        if vectorMatClass[i] == 1
            if matrixClasificacion[indice, 1] == 1
                matAux[i,1] = 1.0
            else
                matAux[i,2] = 1.0
            end
            indice = indice + 1
        end
    end
    return matAux
end

function fastSparceClustering(A, c)
    dim = size(A)
    matrixData = deepcopy(A)
    matClas = zeros(dim[1],c)
    indiceMayor = 1
    matClas[:,1] = ones(dim[1])
    for j=1:c-1
        mat = kernelClustering(matrixData, 2, 2)
        auxMat = separarClasificacion(matClas[:,indiceMayor] , mat)
        matClas[:,indiceMayor] = auxMat[:,1]
        matClas[:,j+1] = auxMat[:,2]
        auxdt = sum(matClas)
        indiceMayor = seleccionaMaximoElementos(matClas)
        matrixData = extraeSubMatriz(matClas[:,indiceMayor] , A)
    end
    aux = zeros(dim[1])
    for j=1:dim[1]
        aux[j] = indiceMaxVector(matClas[j,:])[1]
    end
    ModuleEnviroment.generateImage( 837, 837,trunc.(Int64, aux),"nuevo salida $c")
    return matClas, aux
end

function writeDataCSV(name, matrix)
  dt = DataFrame(matrix)
  #writetable("out/"*string(name)*".csv", dt, separator = ';', header = true)
  CSV.write("out/"*string(name)*".csv", dt, separator = ';')
end

#
function loadCSV()
    aux = CSV.read("D:/data.csv", delim=';')
    auxdata = [aux[i,j] for i=1:5000,j=1:2]
end


function creaDiccionario(mat)
    dim = size(mat)
    promedio = sum([mat[j,:] for j=1:dim[1]]) / dim[1]
    matrixProyection = [sum(mat[j,:].*promedio) for j=1:dim[1]]
    puntoCorte = sum(matrixProyection)/dim[1]
    subgrupos = zeros(dim[1],2)
    for j=1:dim[1]
        if matrixProyection[j] > puntoCorte
            subgrupos[j,1] = 1
        else
            subgrupos[j,2] = 1
        end
    end
    centroides = zeros(2, dim[2])
    #calcular centroides
    for i=1:2
        aux = mat.*subgrupos[:,i] / sum(subgrupos[:,i])
        centroides[i,:] = sum([aux[j ,:] for j=1:dim[1]])
    end
end
