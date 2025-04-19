using ITensors, ITensorMPS
using LinearAlgebra


function calculate_heisenberg_model(N=20)
    sites = siteinds("S=1/2", N)
    ampo = AutoMPO()
    for j in 1:N-1
        ampo += 0.5, "S+", j, "S-", j+1
        ampo += 0.5, "S-", j, "S+", j+1
        ampo += "Sz", j, "Sz", j+1
    end
    H = MPO(ampo, sites)
    psi0 = randomMPS(sites, 10)
    sweeps = Sweeps(5)
    setmaxdim!(sweeps, 10, 20, 100, 200, 400)
    setcutoff!(sweeps, 1E-10)
    energy, psi = dmrg(H, psi0, sweeps)
    entropies = [my_entanglement_entropy(psi, j) for j in 1:N-1]
    return energy, entropies
end

function my_entanglement_entropy(psi::MPS, bond::Int)
    N = length(psi)
    @assert 1 <= bond < N "Bond index out of range"

  
    psi_copy = copy(psi)
    orthogonalize!(psi_copy, bond)

    A = psi_copy[bond]
    site_index = siteind(psi_copy, bond)
    right_link = linkind(psi_copy, bond)

    
    U, S, V = svd(A, (site_index, right_link))

 
    singular_values = diag(S)
    singular_values ./= norm(singular_values)


    entropy = -sum(s -> s^2 * log2(s^2), singular_values)

    return isnan(entropy) ? 0.0 : entropy
end

function main()
    energy, entropies = calculate_heisenberg_model()
    println("基态能量: ", energy)
    println("\n纠缠熵分布 (使用自定义函数计算):")
    for (j, S) in enumerate(entropies)
        println("Bond $j: S = ", round(S; digits=4))
    end
end

main()
