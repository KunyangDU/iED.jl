include("../../src/iED.jl")
include("../utils.jl")

function main_U1(N,m;Jz=1,Jxy=1)
    # @show N
    state = U1Symm(N, m)
    H = zeros(length(state),length(state))
    for (ia,a) in enumerate(state)
        for i in 1:N-1
            j = mod(i,N) + 1
            # j = i + 1
            if bit(a,i) == bit(a,j)
                H[ia,ia] += Jz/4
            else
                H[ia,ia] += - Jz/4
                b = flip(a,[i,j])
                ib = findfirst(x -> x == b, state)
                H[ia,ib] += Jxy/2
            end
        end
    end
    return H,state
end
function S2_U1(N,m)
    # @show N
    state = U1Symm(N, m)
    S = zeros(length(state),length(state))
    for (ia,a) in enumerate(state)
        for i in 1:N-1,j in 1:N-1
            i == j && continue
            # for j = mod(i,N) + 1
            # j = i + 1
            if bit(a,i) == bit(a,j)
                S[ia,ia] += 1/4
            else
                S[ia,ia] += - 1/4
                b = flip(a,[i,j])
                ib = findfirst(x -> x == b, state)
                S[ia,ib] += 1/2
            end
        end
    end
    return hermitianize(S) + 3*N/4 * diagm(ones(length(state)))
end

dataname = "Heisenberg/obc/data"
time0 = time()

lsN = [14]
params = (Jz =1,Jxy = 1)
lsβ = vcat(2. .^ (-5:1:-1), 1:10)
figsize = (width = 200,height = 200)

for (iN,N) in enumerate(lsN)
    @show N
    lsm = -N:N
    Hdata = Dict()

    for m in lsm
        H,state = main_U1(N,m;params...)
        F = eigen(H)
        E = real.(F.values)
        U = F.vectors
        Hdata[(m,)] = Dict(
            "H" => H,
            "E" => E,
            "U" => U,
            "state" => state,
            "S2" => S2_U1(N,m),
        )
        
    end

    Fs = zeros(length(lsβ))
    Cs = zeros(length(lsβ))
    χs = zeros(length(lsβ))
    for (iβ,β) in enumerate(lsβ)
        Obs = Dict()

        for (m,) in keys(Hdata)
            H = Hdata[(m,)]["H"]
            E = Hdata[(m,)]["E"]
            U = Hdata[(m,)]["U"]
            S2 = Hdata[(m,)]["S2"]

            obs = Dict()

            obs["expE"] = @. exp(-β*E)
            obs["H"] = E
            obs["H2"] = E .^ 2
            obs["m"] = m*ones(length(Hdata[(m,)]["state"]))
            obs["m2"] = diag(real.(U'*S2*U)) / 3

            Obs[(m,)] = obs
        end

        F,C,χ = let 
            Z = sum([sum(Obs[key]["expE"]) for key in keys(Obs)])
            obs = zeros(4)
            for (m,) in keys(Obs)
                obs += map(x -> sum(Obs[(m,)]["expE"] .* Obs[(m,)][x]) / Z,["H2","H","m2","m"])
            end
            C = @. (obs[1] - obs[2] ^2) * β ^2 
            χ = @. (obs[3] - obs[4] ^2) * β
            
            - log(Z) / β, C, χ
        end
        Fs[iβ] = F
        Cs[iβ] = C
        χs[iβ] = χ
    end

    fig = Figure()
    axc = Axis(fig[1,1];figsize...,
    xscale = log10,
    xlabel = L"T\cdot J",ylabel=L"C")
    axχ = Axis(fig[1,2];figsize...,
    xscale = log10,
    xlabel = L"T\cdot J",ylabel=L"\chi")

    scatterlines!(axc,1 ./ lsβ,Cs / N,linewidth = 2)
    scatterlines!(axχ,1 ./ lsβ,χs / N,linewidth = 2)

    period = time() - time0

    Label(fig[1,1:2][1, 1, TopLeft()], "Heisenberg chain 1x$(N), U(1),  $(round(period;digits=2))s",
            fontsize = 18,
            font = :bold,
            padding = (0, -450, 15, 0),
            halign = :center
    )

    resize_to_layout!(fig)
    display(fig)

    save("Heisenberg/figures/heisenberg_obc_1x$(N)_$(params).pdf",fig)

    data = Dict(
        "F" => Fs,
        "C" => Cs,
        "χ" => χs,
        "Hdata" => Hdata
    )
    @save "$(dataname)/data_obc_1x$(N)_$(params).jld2" data
end


