

include("../Heisenberg/utils.jl")

function getnd(N,state)
    state = bitstring(state)[end-N+1:end]
    nd = let 
        cnt = 0
        for i in 1:div(N,2)
            reverse(state)[2*i-1:2i] == "11" && (cnt += 1)
        end
        cnt
    end
    return nd
end

Nx,Ny = 2,2
N = 2*Nx*Ny

getnd(N,)
