module matrixOperation


function matrixToVector(matrixA)
    dim = size(matrixA)
    vector = Vector{Int64}(dim[1]*dim[2])
    for j=1:dim[1]
        vector[(j-1)*dim[2]+1: j*dim[2]] = matrixA[j,:]
    end
    return vector
end

#calculate de row average
function rowSubsetAverage(matrixData, column)
    dim = size(matrixData)
    return sum(matrixData[:,column])/dim[1]
end

#calculate the column average
function columnSubsetAverage(matrixData, row)
    dim = size(matrixData)
    return sum(matrixData[row,:])/dim[2]
end

#calculate matriz average
function submatrixAverage(matrixData)
    dim = size(matrixData)
    return sum(sum(matrixData))/(dim[1]*dim[2])
end

function submatrixAverage_VectorExclude(matrixData, vector_I, vector_J)
    dim = size(matrixData)
    aux = matrixData.*vector_I'.*vector_J
    return sum(sum(aux))/(dim[1]*dim[2])
end

function rowMatrixAllZeros(matrixData, row)
    matrixData[row,:] = 0
end

function columnMatrixAllZeros(matrixData, column)
    matrixData[:,column] = 0
end


end
