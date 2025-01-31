

include("utils.jl")

function main_U1(N = 4)
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

    return - eigvals(H)[1] / N
end

function main_U1pLatt(N = 4)
    lsk = (-N/2+1):N/2
    state = U1Symm(N,0)
    E0s = Vector(undef,N)
    for (ik,k) in enumerate(lsk)
        pstate,Rs = LattSymm(N, k, state)
        isempty(pstate) && continue
        L = length(pstate)
        H = zeros(L,L) * 1im
        for (ia,a) in enumerate(pstate)
            for i in 1:N
                j = mod(i,N) + 1
                if bit(a,i) == bit(a,j)
                    H[ia,ia] += 1/4
                else
                    H[ia,ia] += - 1/4
                    b = flip(a,[i,j])
                    b,l = gauge(b, N)
                    ib = findfirst(x -> x == b, pstate)
                    if !isnothing(ib)
                        H[ia,ib] += 1/2 * sqrt(Rs[ia]/Rs[ib]) * exp(1im*k*l * 2pi / N)
                    end
                end
            end
        end
        H = hermitianize(H)
        E0s[ik] = - eigvals(H)[1] / N
    end
    E0s
end

N = 16
#@time E1 = main_U1(N)
#@time E2 = main_U1pLatt(N)
#@show maximum(E2) ≈ E1;
#E2

#= 
reversion symmetry、spin inversion symmetry的哈密顿量
热容、磁化率的热力学量
=#




