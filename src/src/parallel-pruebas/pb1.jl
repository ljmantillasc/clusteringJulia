
c1 = Channel(32)
c2 = Channel(32)

function foo()
    i = 1
    while i>300
        data = take!(c1)
        result = data + rand(1)
        put!(c2, result)
    end
end

for _ in 1:3
    @async foo()
end
