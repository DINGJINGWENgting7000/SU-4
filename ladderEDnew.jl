using LinearAlgebra
using SparseArrays
using KrylovKit

L = 12           # Number of sites (ladder length)
Dim_site = 6              # Local Hilbert space dimension (6*6)
J = 1.0             # Coupling strength (J)

# SU(4) generators and operators
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
Em = [sparse(Matrix(e')) for e in Ep]  # Convert adjoints to sparse matrices

#one-site operators
function embed_one_site(op::SparseMatrixCSC, i::Int, L::Int)
    left = sparse(I, Dim_site^(i - 1), Dim_site^(i - 1))
    right = sparse(I, Dim_site^(L - i), Dim_site^(L - i))
    return kron(kron(left, op), right)
end

# two-site operator
function embed_two_site(op1::SparseMatrixCSC, op2::SparseMatrixCSC, i::Int, L::Int)
    left = sparse(I, Dim_site^(i - 1), Dim_site^(i - 1))
    right = sparse(I, Dim_site^(L - i - 1), Dim_site^(L - i - 1))
    return kron(kron(kron(left, op1), op2), right)
end

# 
function build_hamiltonian(L, deltaU)
    H = spzeros(Dim_site^L, Dim_site^L)

    # coupling terms
    for i in 1:(L-1)
        Hij = spzeros(Dim_site^L, Dim_site^L)
        for a in 1:6
            opL = Em[a]
            opR = Ep[a]
            Hij += J/2*embed_two_site(opL, opR, i, L)   # T^a_i T^a_{i+1}
            Hij += J/2*embed_two_site(opR, opL, i, L)   # T^a_{i+1} T^a_i
        end
        for Hgen in (H1, H2, H3)
            Hij += J*embed_two_site(Hgen, Hgen, i, L)   # H_i H_{i+1}
        end
        H += Hij  
    end

    # on-site terms
    for i in 1:L
        H += deltaU * embed_one_site(Sd, i, L)  
    end

    return H
end

# diagonalization function using KrylovKit
function diagonalize(H, n=5)
    dimH = size(H, 1)
    vals, _, _ = eigsolve(x -> H * x, dimH, n, :SR, tol=1e-6)
    return vals
end

# Print the results
results = []
open("scan_deltaU_ED.txt", "w") do file
    println(file, "\n\n=======================")
    println(file, "Exact Diagonalization of SU(4) ladder Model")
    println(file, "=======================")
    println(file, "L = $L, J = $J")
    println(file, "=======================")
    deltaU_list = -2.4:0.2:2.4
    deltaU_values = collect(deltaU_list)
    for deltaU_val in deltaU_values
        H = build_hamiltonian(L, deltaU_val)
        vals = diagonalize(H, 5)
        push!(results, (deltaU_val, vals[1]))
        println(file, "deltaU = $deltaU_val, Ground state energy = $(vals[1])")
    end
    println(file, " ============= DONE =============")
end

# Read the results from the file and plot
deltaUs = Float64[]
energies = Float64[]
open("scan_deltaU_ED.txt", "r") do file
    seen_deltaUs = Set{Float64}()  # To track unique deltaU values
    for line in eachline(file)
        if occursin("deltaU", line)
            match_result = match(r"deltaU\s*=\s*([-0-9.]+),\s*Ground state energy\s*=\s*([-0-9.]+)", line)
            if match_result !== nothing
                deltaU_val = parse(Float64, match_result.captures[1])
                energy_val = parse(Float64, match_result.captures[2])
                if !(deltaU_val in seen_deltaUs)  
                    push!(deltaUs, deltaU_val)
                    push!(energies, energy_val)
                    push!(seen_deltaUs, deltaU_val)
                end
            end
        end
    end
end
plot(deltaUs, energies;
     xlabel = "deltaU",
     ylabel = "Ground State Energy",
     title = "Ground State Energy vs deltaU(ED)",
     lw = 2,
     marker = :circle,
     legend = false,
     grid = true)

# Save the plot, overwriting the existing file
savefig("E0-deltaU(ED).png")
println("Plot energy-deltaU.png")