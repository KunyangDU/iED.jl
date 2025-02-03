
include("utils.jl")

function main(N = 4,m=0,lsk=0:N-1;Jxy=1,Jzz=1)
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
                    ib = findfirst(x -> x == b, pstate)
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

N = 16
Hs,stateinfos,ks = main(N,0,[0,])
H = Hs[1]

eigvals(H)

