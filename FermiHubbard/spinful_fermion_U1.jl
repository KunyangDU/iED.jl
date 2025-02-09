using LinearAlgebra,CairoMakie

include("../Heisenberg/utils.jl")
include("utils.jl")

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

function main(Nx,Ny,μ,U,basis = 0:(2^N_orb - 1))
    orbitals = [(i, j, σ) for j in 1:Ny for i in 1:Nx for σ in 1:2]
    # 生成轨道列表 (i, j, spin)
    N_orb = length(orbitals)
    # 构建单粒子哈密顿量（最近邻跃迁）
    H_single = zeros(N_orb, N_orb)
    for (m, (i, j, σ)) in enumerate(orbitals)
        # 寻找最近邻
        neighbors = []
        i < Nx && push!(neighbors, (i+1, j, σ))
        i > 1 && push!(neighbors, (i-1, j, σ))
        j < Ny && push!(neighbors, (i, j+1, σ))
        j > 1 && push!(neighbors, (i, j-1, σ))
        
        # 设置跃迁矩阵元
        for (ni, nj, nσ) in neighbors
            n = findfirst(orb -> orb == (ni, nj, nσ), orbitals)
            isnothing(n) && continue
            H_single[m, n] = -1.0
            H_single[n, m] = -1.0  # 保证厄米性
        end
    end

    # 生成所有基矢（二进制位表示占据情况）
    N_states = length(basis)

    # 构建多体哈密顿量
    H = zeros(Float64, N_states, N_states)
    for ket in basis
        ket_idx = findfirst(x -> x == ket,basis)
        H[ket_idx,ket_idx] += -(μ + U/2) * count_ones(ket)
        maskup = convert(Int64,sum(2. .^ ((1:2:N_orb-1)) ))
        maskdn = convert(Int64,sum(2. .^ ((0:2:N_orb-2)) ))
        up = (ket & maskup) >> 1
        dn = ket & maskdn
        H[ket_idx,ket_idx] += U*count_ones(dn & up)
        for n in 1:N_orb, m in 1:N_orb
            (H_single[m, n] == 0) && continue

            # 湮灭算符c_n作用
            (ket & (1 << (n-1)) == 0) && continue
            ket1 = ket ⊻ (1 << (n-1))
            sign1 = (-1)^count_ones(ket & ((1 << (n-1)) - 1))

            # 产生算符c^†_m作用
            (ket1 & (1 << (m-1)) != 0) && continue
            ket2 = ket1 | (1 << (m-1))
            sign2 = (-1)^count_ones(ket1 & ((1 << (m-1)) - 1))

            # 更新矩阵元
            ket2_idx = findfirst(x -> x == ket2,basis)
            H[ket2_idx, ket_idx] += H_single[m, n] * sign1 * sign2
        end
    end

    # 对角化并输出基态能量
    F = eigen(H)
    eigenvectors = F.vectors
    eigenvalues = F.values
    return eigenvalues,eigenvectors
end


U = 8
Nx,Ny = 2,3
N = Nx*Ny
N_orb = 2N
#lsμ = (U/2 + 2) .* range(-1,1,3*(U+4) + 1)

lsβ = 10 .^ range(log10.([5e-3,20])...,40)
lsT = 1 ./ lsβ

data = Dict()

for μ in 0

    for N0 in 0:2N
        basis = let 
            tmp = []
            for s in 0:2^N_orb-1
                count_ones(s) == N0 && push!(tmp,  s)
            end
            tmp
        end
        E,V = main(Nx,Ny,μ,U,basis)

        tmpdata = Dict(
            "E" => E,
            "V" => V,
            "N" => N0
        )

        data[(μ,N0)] = tmpdata

    end

    Fs = zeros(length(lsβ))
    Us = zeros(length(lsβ))
    Ns = zeros(length(lsβ))
    Ces = zeros(length(lsβ))

    for (iβ,β) in enumerate(lsβ)
        Z = sum(vcat([exp.(- β * data[key]["E"]) for key in keys(data)]...))
        U = sum(vcat([exp.(- β * data[key]["E"]) .* data[key]["E"] for key in keys(data)]...)) / Z
        Ntmp = sum(vcat([exp.(- β * data[key]["E"]) .* data[key]["N"] for key in keys(data)]...)) / Z
        U2 = sum(vcat([exp.(- β * data[key]["E"]) .* (data[key]["E"] .^ 2) for key in keys(data)]...)) / Z

        Fs[iβ] = -log(Z) / β / N
        Us[iβ] = U / N
        Ns[iβ] = Ntmp / N
        Ces[iβ] = (U2 - U^2) * β ^ 2 / N
    end

    figsize = (width = 300,
    height = 200,)

    fig = Figure()
    axf = Axis(fig[1,1];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    scatterlines!(axf,lsT,Fs)

    axU = Axis(fig[1,2];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    scatterlines!(axU,lsT,Us)

    axN = Axis(fig[2,1];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    ylims!(axN,0.9,1.1)

    scatterlines!(axN,lsT,Ns)

    axCe = Axis(fig[2,2];
    xscale = log10,
    figsize...,
    xlabel = L"T")

    scatterlines!(axCe,lsT,Ces)

    resize_to_layout!(fig)
    display(fig)


end



