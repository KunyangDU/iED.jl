
using JLD2,LinearAlgebra

include("model.jl")
include("../src/iED.jl")
# 晶格参数设置
Lx = 4   # x方向格点数
Ly = 1   # y方向格点数
N = Lx * Ly
t = 1.0  # 跃迁振幅
L = Lx*Ly
# 生成邻接表（包含所有最近邻有向边）
function generate_neighbors()
    neighbors = Tuple{Int,Int}[]
    for site in 1:N
        x = (site-1) % Lx + 1
        y = (site-1) ÷ Lx + 1
        
        # 右邻
        if x < Lx
            push!(neighbors, (site, site+1))
        end
        # 左邻
        if x > 1
            push!(neighbors, (site, site-1))
        end
        # 上邻
        if y < Ly
            push!(neighbors, (site, site+Lx))
        end
        # 下邻
        if y > 1
            push!(neighbors, (site, site-Lx))
        end

        if Ly != 1
            if y == 1
                push!(neighbors, (site, site+Lx*(Ly-1)))
            end 
            if y == Ly
                push!(neighbors, (site, site-Lx*(Ly-1)))
            end
        end
        
    end
    return neighbors
end

# 生成所有邻接关系
neighbors = generate_neighbors()

# 将整数状态转换为有序占据位列表
function state_to_occupation(s::Int)
    occ = Int[]
    for i in 1:N
        if (s & (1 << (i-1))) != 0
            push!(occ, i)
        end
    end
    return sort!(occ)
end

# 将占据列表转换为整数状态
function occupation_to_state(occ::Vector{Int})
    s = 0
    for i in occ
        s |= 1 << (i-1)
    end
    return s
end

# 计算费米子算符作用后的符号
function compute_sign(original::Vector{Int}, new::Vector{Int})
    # 通过有序列表比较确定置换次数
    sign = 1
    i = 1
    j = 1
    while i <= length(original) && j <= length(new)
        if original[i] == new[j]
            i += 1
            j += 1
        else
            sign *= -1
            j += 1
        end
    end
    return sign
end

# 构建完整哈密顿量
function build_full_hamiltonian(states = 0:2^N-1)
    dim = length(states)
    H = zeros(Float64, dim, dim)
    
    for s in states
        occ = state_to_occupation(s)
        
        # 遍历所有邻接对
        for (i,j) in neighbors
            # 湮灭j产生i的过程
            if (s & (1 << (j-1))) != 0  # 检查j是否被占据
                # 创建中间态（湮灭j）
                temp_occ = filter(x -> x != j, occ)
                
                if !(i in temp_occ)  # 检查i是否未被占据
                    # 插入i并保持有序
                    pos = searchsortedfirst(temp_occ, i)
                    new_occ = insert!(copy(temp_occ), pos, i)
                    
                    # 计算符号
                    @show findfirst(==(j), occ),occ,j
                    sign = (-1)^(findfirst(==(j), occ)-1) * (-1)^(pos-1)
                    
                    # 转换到新状态
                    s_prime = occupation_to_state(new_occ)
                    
                    # 更新矩阵元
                    H[s_prime+1, s+1] += t * sign
                end
            end
        end
    end
    
    # 确保厄米性
    for i in 1:dim, j in i+1:dim
        H[i,j] = H[j,i] = (H[i,j] + H[j,i])/2
    end
    
    return H
end

N = div(Lx*Ly,2)
fulls = 0:2^(Lx*Ly)-1
states = Int64[]
for s in fulls
    if count1s(s) == N
        push!(states,s)
    end
end

# 生成并显示2x2晶格哈密顿量
H = build_full_hamiltonian(states)

E,~ = calcEg(H)
E/L - ue(100,Lx,Ly)

