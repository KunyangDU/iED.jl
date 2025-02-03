using LinearAlgebra

# code by DeepSeek

# 定义泡利矩阵 (自旋1/2算符)
Sx = [0 1; 1 0]
Sy = [0 -im; im 0]
Sz = [1 0; 0 -1]
S0 = diagm(ones(2))

# 计算多个量子算符的张量积
function tensor_product(ops)
    result = ops[1]
    for op in ops[2:end]
        result = kron(result, op)
    end
    return result
end

# 构建一维海森堡模型哈密顿量（开放边界条件）
function construct_hamiltonian(N)
    dim = 2^N
    H = zeros(ComplexF64, dim, dim)
    
    # 遍历所有最近邻对
    for i in 1:N-1
        # 生成所有位置的单位矩阵列表
        ops_left = [S0 for _ in 1:i-1]
        ops_right = [S0 for _ in 1:N-1-i]
        
        # 构造相互作用项 S·S = Sx⊗Sx + Sy⊗Sy + Sz⊗Sz
        for S in [Sx, Sy, Sz]
            term = tensor_product(vcat(ops_left, [S, S], ops_right))
            H += 0.25 * term  # 0.25来自自旋算符的1/2系数
        end
    end

    for S in [Sx, Sy, Sz]
        term = tensor_product(vcat([S,],repeat([S0,], N-2),[S,]))
        H += 0.25 * term  # 0.25来自自旋算符的1/2系数
    end
    
    return (H)  # 海森堡模型哈密顿量为实对称矩阵
end

# 系统参数
N = 12  # 自旋链长度

# 构建哈密顿量并对角化
H = construct_hamiltonian(N)
eigenvalues = eigvals(Hermitian(H))

println("系统尺寸: $N 个自旋")
eigenvalues[1] / N 

