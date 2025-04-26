using LinearAlgebra
include("model.jl")
include("../src/iED.jl")
# # 辅助函数：快速二进制1的计数
# count1s(x::Int) = sum(ones(Int, sizeof(x)*8) do n, i
#     n + ((x >> (i-1)) & 1)
# end)

# 辅助函数：在排序数组中二分查找
function bifind(arr, val, start=1, stop=length(arr))
    while start ≤ stop
        mid = (start + stop) ÷ 2
        if arr[mid] < val
            start = mid + 1
        else
            stop = mid - 1
        end
    end
    return start ≤ length(arr) && arr[start] == val ? start : nothing
end

function main(Nx, Ny, μ, U, basis=0:(2^(Nx*Ny)-1))
    # 生成轨道列表 (i, j)
    orbitals = [(i, j) for j in 1:Ny for i in 1:Nx]
    N_orb = length(orbitals)
    
    # 构建单粒子哈密顿量（最近邻跃迁）
    H_single = zeros(N_orb, N_orb)
    for (m, (i, j)) in enumerate(orbitals)
        # 寻找最近邻
        neighbors = []
        i < Nx && push!(neighbors, (i+1, j))
        i > 1 && push!(neighbors, (i-1, j))
        j < Ny && push!(neighbors, (i, j+1))
        j > 1 && push!(neighbors, (i, j-1))
        
        # 周期性边界条件
        if Ny != 1
            j == Ny && push!(neighbors, (i, 1))
            j == 1 && push!(neighbors, (i, Ny))
        end
        
        # 设置跃迁矩阵元
        for (ni, nj) in neighbors
            n = findfirst(orb -> orb == (ni, nj), orbitals)
            isnothing(n) && continue
            H_single[m, n] = -1.0
            H_single[n, m] = -1.0  # 保证厄米性
        end
    end

    # 生成所有基矢（二进制位表示占据情况）
    N_states = length(basis)

    # 构建多体哈密顿量
    H = zeros(Float64, N_states, N_states)
    for (iket, ket) in enumerate(basis)
        # 化学势项
        H[iket, iket] += -μ * count1s(ket)
        
        # 相互作用项（若需要）
        if U != 0
            H[iket, iket] += U * count1s(ket & (ket >> 1)) # 示例相互作用形式
        end
        
        # 跃迁项
        for n in 1:N_orb, m in 1:N_orb
            (H_single[m, n] == 0) && continue

            # 湮灭算符c_n作用
            (ket & (1 << (n-1)) == 0) && continue
            ket1 = ket ⊻ (1 << (n-1))
            sign1 = (-1)^count1s(ket & ((1 << (n-1)) - 1))

            # 产生算符c^†_m作用
            (ket1 & (1 << (m-1)) != 0) && continue
            ket2 = ket1 | (1 << (m-1))
            sign2 = (-1)^count1s(ket1 & ((1 << (m-1)) - 1))

            # 更新矩阵元
            ket2_idx = bifind(basis, ket2)
            isnothing(ket2_idx) && continue
            H[ket2_idx, iket] += H_single[m, n] * sign1 * sign2
        end
    end

    # 对角化
    F = eigen(Symmetric(H)) # 利用对称性加速
    return F.values, F.vectors
end

# 参数设置
Nx, Ny = 2, 2
N_orb = Nx * Ny
μ = 0.0  # 化学势
U = 0.0  # 相互作用强度（自由费米子时设为0）

# 按粒子数分块对角化
data = Dict()
for N0 in 0:N_orb
    # 生成固定粒子数子空间
    basis = filter(s -> count1s(s) == N0,0:(2^N_orb-1))
    isempty(basis) && continue
    
    # 对角化子空间哈密顿量
    E, V = main(Nx, Ny, μ, U, basis)
    
    # 存储结果
    data[N0] = Dict(
        "E" => E,
        "V" => V,
        "N" => N0
    )
end

# 热力学量计算（示例）
β = 1.0  # 逆温度
Z = sum(exp.(-β * data[N0]["E"]) for N0 in keys(data))
E_avg = sum(exp(-β * data[N0]["E"])' * data[N0]["E"] for N0 in keys(data)) / Z
println("平均能量：", E_avg)