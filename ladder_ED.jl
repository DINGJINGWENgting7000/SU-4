using LinearAlgebra
using SparseArrays
using KrylovKit

L = 4                # Number of sites (ladder length)
dim = 6              # Local Hilbert space dimension (6*6)      # Total Hilbert space dimension
deltaU = 0.5           # Interaction strength (deltaU)
J=1.0             # Coupling strength (J)
# Define the SU(4) generators and operators
Sd = sparse([2 0 0 0 0 0;
             0 0 0 0 0 0;
             0 0 0 0 0 0;
             0 0 0 0 0 0;
             0 0 0 0 0 0;
             0 0 0 0 0 -2])

H1 = sparse([0 0 0 0 0 0;
             0 1/2 0 0 0 0;
             0 0 1/2 0 0 0;
             0 0 0 -1/2 0 0;
             0 0 0 0 -1/2 0;
             0 0 0 0 0 0])

H2 = sparse([2/3 0 0 0 0 0;
             0 -1/3 0 0 0 0;
             0 0 1/3 0 0 0;
             0 0 0 -1/3 0 0;
             0 0 0 0 1/3 0;
             0 0 0 0 0 -2/3])

H3 = sparse([1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 0 -1 0 0 0;
             0 0 0 1 0 0;
             0 0 0 0 -1 0;
             0 0 0 0 0 -1])

Ep = [
    sparse([0 0 0 0 0 0;
            0 0 0 1 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0]),

    sparse([0 0 0 -1 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 1;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0]),

    sparse([0 1 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 1;
            0 0 0 0 0 0]),

    sparse([0 0 0 0 -1 0;
            0 0 0 0 0 -1;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0]),

    sparse([0 0 1 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 -1;
            0 0 0 0 0 0;
            0 0 0 0 0 0]),

    sparse([0 0 0 0 0 0;
            0 0 1 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0])
]

Em = [e' for e in Ep]

# Build the Hamiltonian components for the ladder system
function build_hamiltonian(L=1.0, deltaU=0.5)
    local Coupling = spzeros(dim^2, dim^2) 
    for i in 1:6
        Coupling += 0.5 * (kron(Ep[i], Em[i]) +  kron(Em[i], Ep[i]))
    end
    Coupling += kron(H1, H1) + kron(H2, H2) + kron(H3, H3)
    return Coupling
end


# Build the full Hamiltonian for the ladder system
function build_full_hamiltonian(L, deltaU)
    H = spzeros(dim^(L), dim^(L))  
    local Coupling = build_hamiltonian(L, deltaU)  
    for i in 1:(L - 1)  
        pre = kron(sparse(I, dim^(i-1), dim^(i-1)), Coupling)
        post = kron(pre, sparse(I, dim^(L-i-1), dim^(L-i-1)))
        H += post
    end
    for i in 1:L
        H += deltaU * kron( kron(sparse(I, dim^(i-1), dim^(i-1)), Sd), sparse(I, dim^(L-i), dim^(L-i)) )
    end
    return H
end

# Diagonalize the Hamiltonian to find the eigenvalues
function diagonalize(H, n=5)
    dimH = size(H, 1)
    vals, _, _ = eigsolve(x -> H * x, dimH, n, :SR, tol=1e-6)
    return vals
end

# Build and diagonalize the Hamiltonian for a given system size L and deltaU
H = build_full_hamiltonian(L, deltaU)
vals = diagonalize(H, 5)

# Print the results
println("Ground state energy (L = $L): ", vals[1])
println("Lowest 5 eigenvalues: ")
println(vals)