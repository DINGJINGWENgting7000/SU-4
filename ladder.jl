using ITensors, ITensorMPS


const Dim_site = 6  # Dimension of each site
const L = 4         # Number of ladder sites
const deltaU = 0.5  # Interaction parameter
const J = 1.0       # Coupling constant
# Create the sites for the ladder system
sites = siteinds("IdxSU4_6d", 2L)

# Define the operator space for each site
ITensors.space(::SiteType"IdxSU4_6d") = Dim_site

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

  function build_full_hamiltonian(sites, deltaU)
    ampo = AutoMPO()


    for i in 1:2:(2L-2)
        for j in 1:6
            add!(ampo, 0.5, "Ep$j", i, "Em$j", i+2)
            add!(ampo, 0.5, "Em$j", i, "Ep$j", i+2)
        end
    end

    for i in 1:L
        si = 2*i - 1  
        sj = 2*i      
        for j in 1:6
            add!(ampo, 0.5, "Ep$j", si, "Em$j", sj)
            add!(ampo, 0.5, "Em$j", si, "Ep$j", sj)
        end
    end

    #
    for i in 1:L
        add!(ampo, deltaU, "Sd", i, "Sd", i)
    end

    return MPO(ampo, sites)
end
function build_full_hamiltonian(sites, deltaU, J)
    ampo = AutoMPO()
    L = div(length(sites), 2)

   #= Inter-chain interactions
    #for i in 1:2:(2L - 2)
       # for j in 1:6
       #     add!(ampo, 0.5 * J, "Ep$j", i, "Em$j", i+2)
            add!(ampo, 0.5 * J, "Em$j", i, "Ep$j", i+2)
        end
        add!(ampo, J, "H1", i, "H1", i+2)
        add!(ampo, J, "H2", i, "H2", i+2)
        add!(ampo, J, "H3", i, "H3", i+2)
    end
=#
   ## Intra-chain interactions
    for i in 1:L
        si = 2*i - 1
        sj = 2*i
        for j in 1:6
            add!(ampo, 0.5 * J, "Ep$j", si, "Em$j", sj)
            add!(ampo, 0.5 * J, "Em$j", si, "Ep$j", sj)
        end
        add!(ampo, J, "H1", si, "H1", sj)
        add!(ampo, J, "H2", si, "H2", sj)
        add!(ampo, J, "H3", si, "H3", sj)
    end

    # ΔU 
    for i in 1:2L
        add!(ampo, deltaU, "Sd", i)
    end
    return MPO(ampo, sites)
end


# Build the full Hamiltonian
H = build_full_hamiltonian(sites, deltaU)

# Create a random MPS initial state
psi0 = randomMPS(sites)

# Set up DMRG sweeps
sweeps = Sweeps(6)
setmaxdim!(sweeps, 20, 50, 100, 200)
setcutoff!(sweeps, 1E-10)

# Perform DMRG calculation
energy, psi = dmrg(H, psi0, sweeps)

# Output the result
println("Ground state energy: ", energy)