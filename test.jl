using LinearAlgebra
using Combinatorics

# 系统参数
const Lx = 3  # 晶格大小
const Ly = 3
const N = Lx * Ly  # 总格点数
const t = 1.0  # 跳跃强度
const U = 0.0  # 相互作用强度
const N_up = 2  # 上自旋粒子数
const N_down = 2  # 下自旋粒子数

# 生成所有满足 N_up 和 N_down 的实空间基矢
function generate_states()
    states = []
    # 每个格点用两位表示: (上自旋, 下自旋)
    for bits in 0:(2^(2*N) - 1)
        up = 0
        down = 0
        for i in 0:N-1
            up += (bits >> (2*i)) & 1
            down += (bits >> (2*i + 1)) & 1
        end
        if up == N_up && down == N_down
            push!(states, bits)
        end
    end
    return states
end

states = generate_states()

# 平移操作
function translate(state, dx, dy)
    translated = 0
    for i in 0:Lx-1
        for j in 0:Ly-1
            old_pos = i * Ly + j
            ni = (i + dx) % Lx
            nj = (j + dy) % Ly
            new_pos = ni * Ly + nj
            # 提取原位置的两比特
            bits = (state >> (2*old_pos)) & 0b11
            translated |= bits << (2*new_pos)
        end
    end
    return translated
end

# 生成平移群
translation_group = [(dx, dy) for dx in 0:Lx-1, dy in 0:Ly-1]

# 寻找轨道
visited = Set()
orbits = []
for s in states
    if !(s in visited)
        orbit = []
        for (dx, dy) in translation_group
            ts = translate(s, dx, dy)
            if !(ts in orbit)
                push!(orbit, ts)
            end
        end
        for ts in orbit
            push!(visited, ts)
        end
        push!(orbits, orbit)
    end
end

# 构建动量基矢
k_points = [(2π*n/Lx, 2π*m/Ly) for n in 0:Lx-1, m in 0:Ly-1]
k_bases = Dict()

for orbit in orbits
    for (kx, ky) in k_points
        coeffs = []
        for (dx, dy) in translation_group
            ts = translate(orbit[1], dx, dy)
            phase = exp(-im * (kx*dx + ky*dy))
            push!(coeffs, phase)
        end
        # 归一化
        norm = sqrt(length(coeffs))
        k_state = round(Int64,sum(coeffs) / norm)  # 假设实空间基矢正交
        if abs(k_state) > 1e-6
            if !haskey(k_bases, (kx, ky))
                k_bases[(kx, ky)] = []
            end
            push!(k_bases[(kx, ky)], k_state)
        end
    end
end

# 构建哈密顿量矩阵
H_k = Dict()
for (k, basis) in k_bases
    @show basis
    dim = length(basis)
    H = zeros(ComplexF64, dim, dim)
    # 库仑项 (对角)
    for i in 1:dim
        state = basis[i]
        coulomb = 0
        for site in 0:N-1
            up = (state >> (2*site)) & 1
            down = (state >> (2*site + 1)) & 1
            coulomb += U * up * down
        end
        H[i, i] = coulomb
    end
    # 跳跃项 (简化为动量空间对角项)
    H += -t * 2*(cos(k[1]) + cos(k[2])) * Matrix(I, dim, dim)
    H_k[k] = H
end

# 对角化每个 k 点的哈密顿量
eigenvalues = Dict()
for (k, H) in H_k
    vals = eigvals(H)
    eigenvalues[k] = vals
    println("k = $k: $vals")
end