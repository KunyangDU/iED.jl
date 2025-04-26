using CairoMakie,LaTeXStrings

include("../utils.jl")

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
                    ib = findfirst(x -> x == b, pstate)
                    if !isnothing(ib)
                        H[ia,ib] += (Jxy) * (1/2) * sqrt(Rs[ia]/Rs[ib]) * exp(1im*k*l * 2pi / N)
                    end
                end
            end
        end
    end
    return hermitianize(H) + 3*N/4 * diagm(ones(L))
end

time0 = time()

lsN = [12]
params = (Jzz =1,Jxy = 1)
lsβ = 1 ./ range(0.025,2,100)

figsize = (width = 200,height = 200)

fig = Figure()
axc = Axis(fig[1,1];figsize...,
xlabel = L"T\cdot J",ylabel=L"C")
axχ = Axis(fig[1,2];figsize...,
xlabel = L"T\cdot J",ylabel=L"\chi")

ylims!(axc,-0.02,0.5)
ylims!(axχ,-0.006,0.165)

ChighT = lines!(axc,1 ./ lsβ, 3/13*lsβ.^2,color=:red,label="high-T",linewidth = 2)
lines!(axχ,1 ./ lsβ, 1/4*lsβ,color=:red,label=L"high-T",linewidth = 2)

text!(axc,1.3,0.2;text=L"3/13T",color=:red,fontsize = 16)
text!(axχ,1.2,0.13;text=L"1/4T",color=:red,fontsize = 16)

label = [L"N=4",L"N=8",L"N=12",L"N=16"]

for (iN,N) in enumerate(lsN)
    @show N
    lsm = -N:N
    Hdata = Dict()

    for m in lsm
        Hs,stateinfos,ks = main_U1pLatt(N,m;params...)
        for (i,k) in enumerate(ks)
            (isnothing(Hs[i]) || isnothing(stateinfos[i])) && continue
            H = Hs[i]
            stateinfo = stateinfos[i]
            # k = ks[i]
            pstate = stateinfo["pstate"]
            F = eigen(H)
            E = real.(F.values)
            U = F.vectors
            Hdata[(m,k)] = Dict(
                "H" => H,
                "E" => E,
                "U" => U,
                "stateinfo" => stateinfo,
                "S2" => getS2(N,m,k),
            )
        end
        
    end

    Fs = zeros(length(lsβ))
    Cs = zeros(length(lsβ))
    χs = zeros(length(lsβ))
    for (iβ,β) in enumerate(lsβ)
        Obs = Dict()

        for (m,k) in keys(Hdata)
            H = Hdata[(m,k)]["H"]
            E = Hdata[(m,k)]["E"]
            U = Hdata[(m,k)]["U"]
            S2 = Hdata[(m,k)]["S2"]

            obs = Dict()

            obs["expE"] = @. exp(-β*E)
            obs["H"] = E
            obs["H2"] = E .^ 2
            obs["m"] = m*ones(length(Hdata[(m,k)]["stateinfo"]["pstate"]))
            obs["m2"] = diag(real.(U'*S2*U)) / 3

            Obs[(m,k)] = obs
        end

        F,C,χ = let 
            Z = sum([sum(Obs[mk]["expE"]) for mk in keys(Obs)])
            obs = zeros(4)
            for (m,k) in keys(Obs)
                obs += map(x -> sum(Obs[(m,k)]["expE"] .* Obs[(m,k)][x]) / Z,["H2","H","m2","m"])
            end
            C = @. (obs[1] - obs[2] ^2) * β ^2 
            χ = @. (obs[3] - obs[4] ^2) * β
            
            - log(Z) / β, C, χ
        end
        Fs[iβ] = F
        Cs[iβ] = C
        χs[iβ] = χ
    end

    lines!(axc,1 ./ lsβ,Cs / N,label = label[iN],linewidth = 2)
    lines!(axχ,1 ./ lsβ,χs / N,label = label[iN],linewidth = 2)
    Es = [minimum(Hdata[key]["E"]) for key in keys(Hdata)]
    @show minimum(Es) / N - 1/4
end

period = time() - time0
Legend(fig[1,3],axχ)

Label(fig[1,1:2][1, 1, TopLeft()], "Heisenberg chain, U(1)⊗ℤ,  $(round(period;digits=2))s",
        fontsize = 18,
        font = :bold,
        padding = (0, -450, 15, 0),
        halign = :center
)

resize_to_layout!(fig)
display(fig)

save("Heisenberg/figures/heisenberg.pdf",fig)
