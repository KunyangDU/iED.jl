using LinearAlgebra,CairoMakie,JLD2

include("../Heisenberg/utils.jl")
include("utils.jl")
include("../src/iED.jl")

function getspinocc(N_orb::Int,ket::Int)
    maskup = convert(Int64,sum(2. .^ ((1:2:N_orb-1)) ))
    maskdn = convert(Int64,sum(2. .^ ((0:2:N_orb-2)) ))
    up = (ket & maskup) >> 1
    dn = ket & maskdn
    return up,dn
end

function getmag(N_orb::Int,ket::Int)
    up,dn = getspinocc(N_orb,ket)
    return (count1s(up) - count1s(dn))/2
end
function getnd(N_orb::Int,ket::Int)
    up,dn = getspinocc(N_orb,ket)
    return count1s(up & dn)
end

function ϵ(k)
    return -2sum(cos.(k))
end

function getk(L::Int;condition = :obc)
    if condition == :obc
        return @. pi * (1:L) / (L+1)
    elseif condition == :pbc
        return @. 2pi * (1:L) / L
    end
end

function getk(Lx::Int,Ly::Int)
    if Ly == 1
        lsk = getk(Lx;condition = :obc)
    else
        lskx = getk(Lx;condition = :obc)
        lsky = getk(Ly;condition = :pbc)
        lsk = [[kx,ky] for kx in lskx,ky in lsky][:]
    end
    return lsk
end

function ue(β::Number,Lx::Int,Ly::Int)
    lsk = getk(Lx,Ly)
    lsum = @.  ϵ(lsk) / (1 + exp( β * ϵ(lsk)))
    return sum(lsum) / Lx / Ly
end

function fe(β::Number,Lx::Int,Ly::Int)
    lsk = getk(Lx,Ly)
    return - sum(@. log(1+exp(-β*(ϵ(lsk))))) / β / Lx / Ly
end

function ce(β::Number,Lx::Int,Ly::Int)
    lsk = getk(Lx,Ly)
    return β^2/2 * sum(@. ϵ(lsk)^2/(1 + cosh(β * ϵ(lsk)))) / Lx / Ly
end



function main(Nx,Ny,μ,U,basis = 0:(2^N_orb - 1))
    orbitals = [(i, j, σ) for j in 1:Ny for i in 1:Nx for σ in 1:2]
    N_orb = length(orbitals)
    H_single = zeros(N_orb, N_orb)
    for (m, (i, j, σ)) in enumerate(orbitals)
        neighbors = []
        i < Nx && push!(neighbors, (i+1, j, σ))
        i > 1 && push!(neighbors, (i-1, j, σ))
        j < Ny && push!(neighbors, (i, j+1, σ))
        j > 1 && push!(neighbors, (i, j-1, σ))
        if Ny != 1
            j == Ny &&  push!(neighbors, (i, 1, σ))
            j == 1 && push!(neighbors, (i, Ny, σ))
        end
        
        for (ni, nj, nσ) in neighbors
            n = findfirst(orb -> orb == (ni, nj, nσ), orbitals)
            isnothing(n) && continue
            H_single[m, n] = -1.0
            H_single[n, m] = -1.0
        end
    end

    N_states = length(basis)

    H = zeros(Float64, N_states, N_states)
    for (iket,ket) in enumerate(basis)
        H[iket,iket] += -(μ + U/2) * count1s(ket)
        H[iket,iket] += U*getnd(N_orb,ket)
        for n in 1:N_orb, m in 1:N_orb
            (H_single[m, n] == 0) && continue

            (ket & (1 << (n-1)) == 0) && continue
            ket1 = ket ⊻ (1 << (n-1))
            sign1 = (-1)^count1s(ket & ((1 << (n-1)) - 1))

            (ket1 & (1 << (m-1)) != 0) && continue
            ket2 = ket1 | (1 << (m-1))
            sign2 = (-1)^count1s(ket1 & ((1 << (m-1)) - 1))

            ket2_idx = bifind(basis,ket,iket,ket2)
            H[ket2_idx, iket] += H_single[m, n] * sign1 * sign2
        end
    end

    F = eigen(H)
    eigenvectors = F.vectors
    eigenvalues = F.values
    return eigenvalues,eigenvectors
end

U = 0
Nx,Ny = 2,4
N = Nx*Ny
N_orb = 2N

lsβ = vcat(2. .^ (-15:1:-1), 1:10)
lsT = 1 ./ lsβ

data = Dict()

for μ in 0
    for N0 in N_orb:-1:0,M0 in -div(N,2):1/2:div(N,2)
        Neff = min(N0,N_orb-N0)
        !(M0 in -Neff/2:Neff/2) && continue
        basis = let 
            tmp = []
            for s in 0:2^N_orb-1
                count1s(s) == N0 && getmag(N_orb,s) == M0 && push!(tmp,  s)
            end
            tmp
        end
        isempty(basis) && println("---------") 
        E,V = main(Nx,Ny,μ,U,basis)

        tmpdata = Dict(
            "E" => E,
            "V" => V,
            "N" => N0,
            "M" => M0
        )

        data[(μ,N0,M0)] = tmpdata
    end

    Fs = zeros(length(lsβ))
    Us = zeros(length(lsβ))
    Ns = zeros(length(lsβ))
    Ms = zeros(length(lsβ))
    Ces = zeros(length(lsβ))
    χs = zeros(length(lsβ))

    for (iβ,β) in enumerate(lsβ)
        Z = sum(vcat([exp.(- β * data[key]["E"]) for key in keys(data)]...))
        Uo = sum(vcat([exp.(- β * data[key]["E"]) .* data[key]["E"] for key in keys(data)]...)) / Z
        Ntmp = sum(vcat([exp.(- β * data[key]["E"]) .* data[key]["N"] for key in keys(data)]...)) / Z
        Uo2 = sum(vcat([exp.(- β * data[key]["E"]) .* (data[key]["E"] .^ 2) for key in keys(data)]...)) / Z

        Fs[iβ] = -log(Z) / β / N
        Us[iβ] = Uo / N
        Ns[iβ] = Ntmp / N
        Ces[iβ] = (Uo2 - Uo^2) * β ^ 2 / N
    end

    figsize = (width = 300,
    height = 200,)

    fig = Figure()
    axf = Axis(fig[1,1];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    axU = Axis(fig[1,2];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    axN = Axis(fig[2,1];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    ylims!(axN,0.9,1.1)

    axCe = Axis(fig[2,2];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    if U == 0
        lines!(axf,lsT,2*fe.(lsβ,Nx,Ny);color = :red)
        lines!(axU,lsT,2*ue.(lsβ,Nx,Ny);color = :red)
        lines!(axCe,lsT,2*ce.(lsβ,Nx,Ny);color = :red)
        scatter!(axf,lsT,Fs)
        scatter!(axU,lsT,Us)
        scatter!(axN,lsT,Ns)
        scatter!(axCe,lsT,Ces)
    else
        scatterlines!(axf,lsT,Fs)
        scatterlines!(axU,lsT,Us)
        scatterlines!(axN,lsT,Ns)
        scatterlines!(axCe,lsT,Ces)
    end

    

    resize_to_layout!(fig)
    display(fig)

    save("FermiHubbard/figures/spinful_hubbard_$(Nx)x$(Ny)_U=$(U).pdf",fig)
    save("FermiHubbard/figures/spinful_hubbard_$(Nx)x$(Ny)_U=$(U).png",fig)

    @save "FermiHubbard/data/data_$(Nx)x$(Ny)_U=$(U).jld2" data
    @save "FermiHubbard/data/Fs_$(Nx)x$(Ny)_U=$(U).jld2" Fs
    @save "FermiHubbard/data/Us_$(Nx)x$(Ny)_U=$(U).jld2" Us
    @save "FermiHubbard/data/Ns_$(Nx)x$(Ny)_U=$(U).jld2" Ns
    @save "FermiHubbard/data/Ces_$(Nx)x$(Ny)_U=$(U).jld2" Ces
end



