using ITensors, ITensorMPS, Serialization
using LinearAlgebra
using Plots
using Base: match  
using Statistics  # Import Statistics for the mean function

Dim_site = 6  
L = 20  # Length of the ladder
J = 1.0       
deltaU_list = -2.0:1.0:2.0
deltaU_values = collect(deltaU_list)

# Define the operator space for each site
function ITensors.space(::SiteType"IdxSU4_6d")
  return Dim_site
end
function ITensors.siteind(::SiteType"IdxSU4_6d", i::Int; kwargs...)
  return Index(Dim_site; tags="Site,IdxSU4_6d", kwargs...)
end

sites = siteinds("IdxSU4_6d", L)

## Define the corresponding Operator of other generators. There should be 15 in total.
function ITensors.op(::OpName"Sd", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>1, s=>1] = 2.0+0im
    Op[s'=>6, s=>6] = -2.0+0im
    return Op
end
function ITensors.op(::OpName"H1", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>2, s=>2] = 1/2.0+0im
    Op[s'=>3, s=>3] = 1/2.0+0im
    Op[s'=>4, s=>4] = -1/2.0+0im
    Op[s'=>5, s=>5] = -1/2.0+0im
    return Op
  end
  
  function ITensors.op(::OpName"H2", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>1, s=>1] = 2/3.0+0im
    Op[s'=>2, s=>2] = -1/3.0+0im
    Op[s'=>3, s=>3] = 1/3.0+0im
    Op[s'=>4, s=>4] = -1/3.0+0im
    Op[s'=>5, s=>5] = 1/3.0+0im
    Op[s'=>6, s=>6] = -2/3.0+0im
    return Op
  end
  
  function ITensors.op(::OpName"H3", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>1, s=>1] = 1.0+0im
    Op[s'=>2, s=>2] = 1.0+0im
    Op[s'=>3, s=>3] = -1.0+0im
    Op[s'=>4, s=>4] = 1.0+0im
    Op[s'=>5, s=>5] = -1.0+0im
    Op[s'=>6, s=>6] = -1.0+0im
    return Op
  end
  
  function ITensors.op(::OpName"Ep1", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>2, s=>4] = 1.0+0im
    Op[s'=>3, s=>5] = 1.0+0im
    return Op
  end
  
  function ITensors.op(::OpName"Ep2", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>1, s=>4] = -1.0+0im
    Op[s'=>3, s=>6] = 1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Ep3", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>1, s=>2] = 1.0+0im
    Op[s'=>5, s=>6] = 1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Ep4", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>1, s=>5] = -1.0+0im
    Op[s'=>2, s=>6] = -1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Ep5", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>1, s=>3] = 1.0+0im
    Op[s'=>4, s=>6] = -1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Ep6", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>2, s=>3] = 1.0+0im
    Op[s'=>4, s=>5] = 1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Em1", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>4, s=>2] = 1.0+0im
    Op[s'=>5, s=>3] = 1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Em2", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>4, s=>1] = -1.0+0im
    Op[s'=>6, s=>3] = 1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Em3", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>2, s=>1] = 1.0+0im
    Op[s'=>6, s=>5] = 1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Em4", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>5, s=>1] = -1.0+0im
    Op[s'=>6, s=>2] = -1.0+0im
    return Op
  end
  
  function ITensors.op( ::OpName"Em5", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>3, s=>1] = 1.0+0im
    Op[s'=>6, s=>4] = -1.0+0im
    return Op
  end
  
  function ITensors.op(::OpName"Em6", ::SiteType"IdxSU4_6d", s::Index)
    Op = emptyITensor(ComplexF64, s, s')
    Op[s'=>3, s=>2] = 1.0+0im
    Op[s'=>5, s=>4] = 1.0+0im
    return Op
  end


  #Dict communication for SU(4) generators
function generate_su4_generators(s::Index)
    stype = SiteType("IdxSU4_6d")
    gens = Dict{String, ITensor}()

    # Cartan generators
    for name in ["H1", "H2", "H3"]
        gens[name] = op(OpName(name), stype, s)  # Convert String to OpName
    end

    # Raising operators
    for i in 1:6
        name = "Ep$(i)"
        gens[name] = op(OpName(name), stype, s)  # Convert String to OpName
    end

    # Lowering operators
    for i in 1:6
        name = "Em$(i)"
        gens[name] = op(OpName(name), stype, s)  # Convert String to OpName
    end

    # Optional
    gens["Sd"] = op(OpName("Sd"), stype, s)  # Convert String to OpName
    return gens
end

function commutator_norm(A::ITensor, B::ITensor)
  comm = commutator(A, B)
  return norm(comm)
end
function check_all_commutators()
  s = siteind("IdxSU4_6d", 1)
  gens = generate_su4_generators(s)
  for (nameA, A) in gens
      for (nameB, B) in gens
          if nameA != nameB
              nrm = commutator_norm(A, B)
              println("Commutator of $nameA and $nameB: norm = ", round(nrm, digits=6))
          end
      end
  end
end

# Define the Hamiltonian for the ladder system
function su4_ladder_hamiltonian(J::Real, deltaU::Real, L::Int)
  ampo = AutoMPO()
  # Diagonal generators (Cartan elements)
  diag_gens = ["H1", "H2", "H3"]
  # Off-diagonal generator pairs
  offdiag_pairs = [
    ("Ep1", "Em1"),
    ("Ep2", "Em2"),
    ("Ep3", "Em3"),
    ("Ep4", "Em4"),
    ("Ep5", "Em5"),
    ("Ep6", "Em6"),
    ("Em1", "Ep1"),
    ("Em2", "Ep2"),
    ("Em3", "Ep3"),
    ("Em4", "Ep4"),
    ("Em5", "Ep5"),
    ("Em6", "Ep6"),
  ]

  # Nearest-neighbor interaction terms
  for i in 1:(L - 1)
    for gen in diag_gens
      add!(ampo, J, gen, i, gen, i+1)
    end
    for (gen1, gen2) in offdiag_pairs
      add!(ampo, J/2, gen1, i, gen2, i+1)
    end
  end
  # Onsite term: deltaU * Sd
  for i in 1:L
    add!(ampo, deltaU, "Sd", i)
  end
  return MPO(ampo, sites)
end

# Define the result structure
struct DMRGResult
  deltaU::Float64
  energy::Float64
  psi::MPS
end
results = DMRGResult[]

# Perform DMRG calculation and save results to a file
open("scandeltaU_dmrg.txt", "w") do file
    println(file, "\n\n=======================")
    println(file, "DMRG of SU(4) ladder Model")
    println(file, "=======================")
    println(file, "L = $L, J = $J")
    println(file, "=======================")
    
    for deltaU_val in deltaU_values
        H = su4_ladder_hamiltonian(J, deltaU_val, L)  
        sweeps = Sweeps(5) 
        maxdim!(sweeps, 50, 100, 150, 200, 200)
        cutoff!(sweeps, 1e-14)
        psi0 = randomMPS(sites)  
        energy, psi = dmrg(H, psi0, sweeps) 
        println(file, "deltaU = $deltaU_val, Ground state energy = $energy")
        filename = "gs_mps_deltaU_$(deltaU_val).itensors"
        open(filename, "w") do io
            Serialization.serialize(io, psi)
        end
        push!(results, DMRGResult(deltaU_val, energy, psi))  # Store result
    end
    println(file, " ============= DONE =============")
end

# Read the results from the file and plot E0-deltaU(DMRG)
deltaUs = Float64[]
energies = Float64[]
open("scandeltaU_dmrg.txt", "r") do file
    seen_deltaUs = Set{Float64}()
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
     title = "Ground State Energy vs deltaU(DMRG)",
     lw = 2,
     marker = :circle,
     legend = false,
     grid = true)

# Save the plot, overwriting the existing file
savefig("E0-deltaU(DMRG).png")
println("Plot saved as E0-deltaU(DMRG).png")

# H3  
function calculate_corr_vs_r(ψ::MPS, sites, max_r::Int=6)
  L = length(sites)
  corr = zeros(Float64, max_r)
  stype = SiteType("IdxSU4_6d")  # Define the SiteType
  center = div(L, 2) 
  indices = [center - 1, center, center + 1]  
  for r in 1:max_r
    vals = Float64[]
    for i in indices
      if i + r <= L  
        ampo = AutoMPO()
        add!(ampo, 1.0, "H3", i, "H3", i + r)  # Use AutoMPO to define the operator
        mpo = MPO(ampo, sites)  # Convert AutoMPO to MPO
        v = inner(ψ, Apply(mpo, ψ))  # Use Apply to ensure indices match
        push!(vals, real(v))
      end
    end
    corr[r] = mean(vals)  # Use mean from Statistics
  end
  return corr
end

# Plotting correlation vs distance  
max_r = div(L, 2)  
plt = plot(xlabel="Distance r", ylabel="⟨H3_i H3_{i+r}⟩",
           title="Correlation vs Distance for Various δU", lw=2, marker=:circle, grid=true)

for r in results
  corr = calculate_corr_vs_r(r.psi, sites, max_r)
  plot!(plt, 1:max_r, corr, label="δU=$(r.deltaU)", linestyle=:auto, marker=:auto) 
end

savefig("correlation_vs_distance.png")
println("Plots saved.")