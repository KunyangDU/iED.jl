include("../../src/iED.jl")
include("../utils.jl")

function Ham_U1(N)
    @show N
    state = U1Symm(N, 0)
    H = zeros(length(state),length(state))
    for (ia,a) in enumerate(state)
        for i in 1:N-1
            j = mod(i,N) + 1
            # j = i + 1
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
lsN = [14,]
foldername = "Heisenberg/data"
lsE = zeros(length(lsN))
for (i,N) in enumerate(lsN)
    H = Ham_U1(N)
    E,~ = calcEg(H)
    lsE[i] = E / N - 1/4
    @show E / N - 1/4
end
@save "$(foldername)/lsN_obc.jld2" lsN
@save "$(foldername)/lsE_obc.jld2" lsE
