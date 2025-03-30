using Arpack
using LinearAlgebra
using SparseArrays
using Plots

struct Block
    length::Int
    basis_size::Int
    operator_dict::Dict{Symbol,AbstractMatrix{Float64}}
end

struct EnlargedBlock
    length::Int
    basis_size::Int
    operator_dict::Dict{Symbol,AbstractMatrix{Float64}}
end


isvalid(block::Union{Block,EnlargedBlock}) =
    all(op -> size(op) == (block.basis_size, block.basis_size), values(block.operator_dict))

model_d = 2  

Sz1 = [0.5 0.0; 0.0 -0.5]  
Sp1 = [0.0 1.0; 0.0 0.0]  

H1 = [0.0 0.0; 0.0 0.0] 

function H2(Sz1, Sp1, Sz2, Sp2)  
    J = 1.0
    Jz = 1.0
    return (J / 2) * (kron(Sp1, Sp2') + kron(Sp1', Sp2)) + Jz * kron(Sz1, Sz2)
end


initial_block = Block(1, model_d, Dict{Symbol,AbstractMatrix{Float64}}(
    :H => H1,
    :conn_Sz => Sz1,
    :conn_Sp => Sp1,
))

function enlarge_block(block::Block)
 
    mblock = block.basis_size
    o = block.operator_dict


    I1 = sparse(1.0I, model_d, model_d)
    I_block = sparse(1.0I, mblock, mblock)
    enlarged_operator_dict = Dict{Symbol,AbstractMatrix{Float64}}(
        :H => kron(o[:H], I1) + kron(I_block, H1) + H2(o[:conn_Sz], o[:conn_Sp], Sz1, Sp1),
        :conn_Sz => kron(I_block, Sz1),
        :conn_Sp => kron(I_block, Sp1),
    )

    return EnlargedBlock(block.length + 1,
                         block.basis_size * model_d,
                         enlarged_operator_dict)
end

function rotate_and_truncate(operator, transformation_matrix)

    return transformation_matrix' * (operator * transformation_matrix)
end

function single_dmrg_step(sys::Block, env::Block, m::Int)


    @assert isvalid(sys)
    @assert isvalid(env)


    sys_enl = enlarge_block(sys)
    if sys === env  
        env_enl = sys_enl
    else
        env_enl = enlarge_block(env)
    end

    @assert isvalid(sys_enl)
    @assert isvalid(env_enl)

    # Construct the full superblock Hamiltonian.
    m_sys_enl = sys_enl.basis_size
    m_env_enl = env_enl.basis_size
    sys_enl_op = sys_enl.operator_dict
    env_enl_op = env_enl.operator_dict
    I_sys_enl = sparse(1.0I, m_sys_enl, m_sys_enl)
    I_env_enl = sparse(1.0I, m_env_enl, m_env_enl)
    superblock_hamiltonian = kron(sys_enl_op[:H], I_env_enl) + kron(I_sys_enl, env_enl_op[:H]) +
                             H2(sys_enl_op[:conn_Sz], sys_enl_op[:conn_Sp], env_enl_op[:conn_Sz], env_enl_op[:conn_Sp])

  
    superblock_hamiltonian = (superblock_hamiltonian + superblock_hamiltonian') / 2
    (energy,), psi0 = eigs(superblock_hamiltonian, nev=1, which=:SR)


    psi0 = transpose(reshape(psi0, (env_enl.basis_size, sys_enl.basis_size)))
    rho = Hermitian(psi0 * psi0')

   
    fact = eigen(rho)
    evals, evecs = fact.values, fact.vectors
    permutation = sortperm(evals, rev=true)

    my_m = min(length(evals), m)
    indices = permutation[1:my_m]
    transformation_matrix = evecs[:, indices]

    truncation_error = 1 - sum(evals[indices])
    println("truncation error: ", truncation_error)


    new_operator_dict = Dict{Symbol,AbstractMatrix{Float64}}()
    for (name, op) in sys_enl.operator_dict
        new_operator_dict[name] = rotate_and_truncate(op, transformation_matrix)
    end

    newblock = Block(sys_enl.length, my_m, new_operator_dict)

    return newblock, energy
end

function infinite_system_algorithm(L::Int, m::Int)
    block = initial_block
    while 2 * block.length < L
        println("L = ", block.length * 2 + 2)
        block, energy = single_dmrg_step(block, block, m)
        println("E/L = ", energy / (block.length * 2))
    end
end

infinite_system_algorithm(100, 20)


energies = []  # 用于存储每次迭代的能量
function infinite_system_algorithm(L::Int, m::Int)
    block = initial_block
    while 2 * block.length < L
        block, energy = single_dmrg_step(block, block, m)
        push!(energies, energy / (block.length * 2))
    end
    plot(energies, label="Energy vs Iteration", xlabel="Iteration", ylabel="Energy")
end

