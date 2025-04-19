using ITensors, ITensorMPS
using LinearAlgebra
const Dim_site = 6  # Hilbert Dimension of each site
const N = 4         # Number of rungs
const deltaU = 0.5  # Interaction parameter
const J = 1.0       # exchange energy parameter

# Define the operator space for each site
function ITensors.space(::SiteType"IdxSU4_6d")
  return Dim_site
end
function ITensors.siteind(::SiteType"IdxSU4_6d", i::Int; kwargs...)
  return Index(Dim_site; tags="Site,IdxSU4_6d", kwargs...)
end
# Create the sites for the ladder system
sites = siteinds("IdxSU4_6d", N)

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
        gens[name] = op(name, stype, s)
    end

    # Raising operators
    for i in 1:6
        name = "Ep$(i)"
        gens[name] = op(name, stype, s)
    end

    # Lowering operators
    for i in 1:6
        name = "Em$(i)"
        gens[name] = op(name, stype, s)
    end

    # Optional
    gens["Sd"] = op("Sd", stype, s)
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
function su4_ladder_hamiltonian(J::Real, deltaU::Real, N::Int)
  ampo = AutoMPO()

  # Diagonal generators (Cartan elements)
  diag_gens = ["H1", "H2", "H3"]

  # Off-diagonal generators: raising and lowering pairs
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

  # Loop over each site
  for i in 1:N
    # Diagonal terms S_i^2
    for gen in diag_gens
      add!(ampo, J, gen, i, gen, i)
    end

    # Off-diagonal terms S^{ab}_i S^{ba}_i
    for (gen1, gen2) in offdiag_pairs
      add!(ampo, J, gen1, i, gen2, i)
    end

    # Interaction term (second part of the Hamiltonian)
       add!(ampo, deltaU, "Sd", i, "Sd", i)
      end
  return MPO(ampo, sites)
end 

open("output_ladder.txt", "w") do file
  deltaU_values = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
  for deltaU in deltaU_values
      println(file, "Running for deltaU = ", deltaU)
      H = su4_ladder_hamiltonian(J, deltaU, N)
      sweeps = Sweeps(10)
      maxdim!(sweeps, 10, 20, 100, 100, 200)
      cutoff!(sweeps, 1e-15)
      psi0 = randomMPS(sites)
      energy, psi = dmrg(H, psi0, sweeps)
      println(file, "Ground state energy = ", deltaU, ": ", energy)
  end
end
