using KrylovKit,LinearAlgebra

include("../src/iED.jl")



#= N = 100
H = randn(N,N) |> x -> x .+ x'
β = 10

#= @time Es,Us,ϵ,iter = KrylovTruncate(β,H)
Us[1] =#
@time eigen(H)
@time calcEg(H) =#

function binary_search_known(arr, v, k, target)
    if target == v
        return k
    elseif target < v
        left = 1
        right = k - 1
    else
        left = k + 1
        right = length(arr)
    end

    while left <= right
        mid = div((left + right),2) # 使用整数除法
        if arr[mid] == target
            return mid
        elseif arr[mid] < target
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return -1  # 表示未找到
end

a = [1,3,4,6,7,8,19,37,45,59,156]

binary_search_known(a,19,7,19)
