module ParallelOperation

function exeParlOperationTwoARG(f::Function,args,arg)
  return pmap(f,args,arg)
end

function exeParlOperationOneARG(f::Function,args)
  return pmap(f,args)
end

end
