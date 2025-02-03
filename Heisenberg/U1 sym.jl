
include("utils.jl")

function Ham_U1(N = 8)
    state = U1Symm(N, 0)
    H = zeros(2^N,2^N)
    for (ia,a) in enumerate(state)
        for i in 1:N
            j = mod(i,N) + 1
            if bit(a,i) == bit(a,j)
                H[ia,ia] += 1/4
            else
                H[ia,ia] += - 1/4
                b = flip(a,[i,j])
                ib = findfirst(x -> x == b, state)
                H[ia,ib] += 1/2
            end
        end
    end

    return H
end

N = 8
H = Ham_U1(N)
Eig = eigen(H)
E = Eig.values 
U = Eig.vectors 

