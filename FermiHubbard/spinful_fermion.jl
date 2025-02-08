using LinearAlgebra,CairoMakie

include("../Heisenberg/utils.jl")
include("utils.jl")
function main(Nx,Ny,μ,U)


    # 生成轨道列表 (i, j, spin)
    orbitals = [(i, j, σ) for j in 1:Ny for i in 1:Nx for σ in 1:2]
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
    basis = 0:(2^N_orb - 1)
    N_states = length(basis)

    # 构建多体哈密顿量
    H = zeros(Float64, N_states, N_states)
    for ket in basis
        ket_idx = Int(ket) + 1
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
            ket2_idx = Int(ket2) + 1
            H[ket2_idx, ket_idx] += H_single[m, n] * sign1 * sign2
        end
    end

    # 对角化并输出基态能量
    F = eigen(H)
    eigenvectors = F.vectors
    eigenvalues = F.values
    return eigenvalues[1],eigenvectors[:,1]
end


U = 8
Nx,Ny = 4,1
N = Nx*Ny
lsμ = (U/2 + 2) .* range(-1,1,3*(U+4) + 1)

Nop = getNop(2*N)

ns = []
Es = []
@time for μ in lsμ
    E,V = main(Nx,Ny,μ,U)
    n = V' * Nop * V / N
    push!(ns,n)
    push!(Es,E)
end

fig = Figure()
ax = Axis(fig[1,1];
title = "Fermi Hubbard model",
width = 400,
height = 200,
ylabel = L"\mu",
xlabel = L"n")

lines!(ax,[0,2],U/2*[1,1];color = :grey,linestyle = :dash)
lines!(ax,[0,2],-U/2*[1,1];color = :grey,linestyle = :dash)
scatterlines!(ax,ns,lsμ;label = "DS-R1 code")
axislegend(ax,position = :lt)
resize_to_layout!(fig)
display(fig)

save("FermiHubbard/figures/spinful_fermion.png",fig)

