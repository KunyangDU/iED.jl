using CairoMakie,LaTeXStrings,KrylovKit

include("utils.jl")
include("../src/iED.jl")


function main_U1pLatt(N = 4,m=0,lsk=0:N-1;Jxy=1,Jzz=1)
    state = U1Symm(N,m)
    Hs = Vector(undef,length(lsk))
    Hs .= nothing
    stateinfos = Vector(undef,length(lsk))
    stateinfos .= nothing
    ks = Vector(undef,length(lsk))
    ks .= nothing
    for (ik,k) in enumerate(lsk)
        pstate,Rs = LattSymm(N, k, state)
        isempty(pstate) && continue
        L = length(pstate)
        H = zeros(L,L) * 1im
        for (ia,a) in enumerate(pstate)
            for i in 1:N
                j = mod(i,N) + 1
                if bit(a,i) == bit(a,j)
                    H[ia,ia] += 1/4 * Jzz
                else
                    H[ia,ia] += - 1/4 * Jzz
                    b = flip(a,[i,j])
                    b,l = LattGauge(b, N)
                    ib = bifind(pstate,a,ia,b)
                    if !isnothing(ib)
                        H[ia,ib] += (Jxy) * 1/2 * sqrt(Rs[ia]/Rs[ib]) * exp(1im*k*l * 2pi / N)
                    end
                end
            end
        end
        Hs[ik] = hermitianize(H)
        stateinfos[ik] = Dict("pstate" => pstate,"Rs" => Rs)
        ks[ik] = k
    end
    return Hs,stateinfos,ks
end

function getS2(N,m,k;Jzz = 1,Jxy = 1)
    pstate,Rs = LattSymm(N, k, U1Symm(N,m))
    isempty(pstate) && return nothing
    L = length(pstate)
    H = zeros(L,L) * 1im
    for (ia,a) in enumerate(pstate)
        for i in 1:N
            for j in 1:N 
                i == j && continue
                if bit(a,i) == bit(a,j)
                    H[ia,ia] += 1/4 * Jzz
                else
                    H[ia,ia] += - 1/4 * Jzz
                    b = flip(a,[i,j])
                    b,l = LattGauge(b, N)
                    ib = bifind(pstate,a,ia,b)
                    if !isnothing(ib)
                        H[ia,ib] += (Jxy) * (1/2) * sqrt(Rs[ia]/Rs[ib]) * exp(1im*k*l * 2pi / N)
                    end
                end
            end
        end
    end
    return hermitianize(H) + 3*N/4 * diagm(ones(L))
end

N = 22
for N in [2,6,12,16,18,20]
params = (Jzz =1,Jxy = 1)

lsm = -N:N
Hdata = Dict()

@time for m in [0]
    #@show m
    #@time "main_U1pLatt" 
    Hs,stateinfos,ks = main_U1pLatt(N,m,[0,div(N,2)];params...)
    for i in eachindex(ks)
        (isnothing(Hs[i]) || isnothing(stateinfos[i])) && continue
        H = Hs[i]
        stateinfo = stateinfos[i]
        k = ks[i]
        pstate = stateinfo["pstate"]
        E,U = calcEg(H)
        Hdata[(m,k)] = Dict(
            "H" => H,
            "E" => E,
            "U" => U,
            "stateinfo" => stateinfo,
            #"S2" => getS2(N,m,k),
        )
    end
    E0 = minimum([Hdata[(m,k)]["E"] for k in [0,div(N,2)]]) / N - 1/4
    @show N,E0
end
end








