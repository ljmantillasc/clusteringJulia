module SecuencialOperation

function exeSecOperationTwoARG(f::Function,args,arg)
  tam = size(args)
  return [f(args[i],arg[i]) for i=1:tam[1]]
end

function exeSecOperationOneARG(f::Function,args)
  tam = size(args)
  return [f(args[i]) for i=1:tam[1]]
end

end
