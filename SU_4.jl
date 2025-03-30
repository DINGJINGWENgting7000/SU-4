using ITensors, ITensorMPS

global Dim_site = 6  ## Define the dimension of Hilbert space in a site.
global N = 4  ## Define the number of sites in the chain

function ITensors.space(::SiteType"IdxSO6")
    return Dim_site
end

## Define the corresponding Operator of other generators. There should be 15 in total.
function ITensors.op(::OpName"H1", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>2, s=>2] = 1/2.0+0im
  Op[s'=>3, s=>3] = 1/2.0+0im
  Op[s'=>4, s=>4] = -1/2.0+0im
  Op[s'=>5, s=>5] = -1/2.0+0im
  return Op
end

function ITensors.op(::OpName"H2", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>1, s=>1] = 2/3.0+0im
  Op[s'=>2, s=>2] = -1/3.0+0im
  Op[s'=>3, s=>3] = 1/3.0+0im
  Op[s'=>4, s=>4] = -1/3.0+0im
  Op[s'=>5, s=>5] = 1/3.0+0im
  Op[s'=>6, s=>6] = -2/3.0+0im
  return Op
end

function ITensors.op(::OpName"H3", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>1, s=>1] = 1.0+0im
  Op[s'=>2, s=>2] = 1.0+0im
  Op[s'=>3, s=>3] = -1.0+0im
  Op[s'=>4, s=>4] = 1.0+0im
  Op[s'=>5, s=>5] = -1.0+0im
  Op[s'=>6, s=>6] = -1.0+0im
  return Op
end

function ITensors.op(::OpName"Ep1", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>2, s=>4] = 1.0+0im
  Op[s'=>3, s=>5] = 1.0+0im
  return Op
end

function ITensors.op(::OpName"Ep2", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>1, s=>4] = -1.0+0im
  Op[s'=>3, s=>6] = 1.0+0im
  return Op
end

function ITensors.op( ::OpName"Ep3", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>1, s=>2] = 1.0+0im
  Op[s'=>5, s=>6] = 1.0+0im
  return Op
end

function ITensors.op( ::OpName"Ep4", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>1, s=>5] = -1.0+0im
  Op[s'=>2, s=>6] = -1.0+0im
  return Op
end

function ITensors.op( ::OpName"Ep5", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>1, s=>3] = 1.0+0im
  Op[s'=>4, s=>6] = -1.0+0im
  return Op
end

function ITensors.op( ::OpName"Ep6", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>2, s=>3] = 1.0+0im
  Op[s'=>4, s=>5] = 1.0+0im
  return Op
end

function ITensors.op( ::OpName"Em1", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>4, s=>2] = 1.0+0im
  Op[s'=>5, s=>3] = 1.0+0im
  return Op
end

function ITensors.op( ::OpName"Em2", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>4, s=>1] = -1.0+0im
  Op[s'=>6, s=>3] = 1.0+0im
  return Op
end

function ITensors.op( ::OpName"Em3", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>2, s=>1] = 1.0+0im
  Op[s'=>6, s=>5] = 1.0+0im
  return Op
end

function ITensors.op( ::OpName"Em4", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>5, s=>1] = -1.0+0im
  Op[s'=>6, s=>2] = -1.0+0im
  return Op
end

function ITensors.op( ::OpName"Em5", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>3, s=>1] = 1.0+0im
  Op[s'=>6, s=>4] = -1.0+0im
  return Op
end

function ITensors.op(::OpName"Em6", ::SiteType"IdxSO6", s::Index)
  Op = emptyITensor(ComplexF64, s, s')
  Op[s'=>3, s=>2] = 1.0+0im
  Op[s'=>5, s=>4] = 1.0+0im
  return Op
end

# Define the sites
sites = siteinds("IdxSO6", N)

#  AutoMPO 
function build_h_1d_chain(sites; J1=1.0, J2=1.0, J3=1.0, J4=1.0, J5=1.0, J6=1.0, J7=1.0, J8=1.0, J9=1.0)
    N = length(sites)  # 系统大小
    ampo = AutoMPO()
    for i in 1:N-1  
        # H
        ampo += J1, "H1", i, "H1", i+1
        ampo += J1, "H2", i, "H2", i+1
        ampo += J1, "H3", i, "H3", i+1
        # E
        ampo += 0.5*J4, "Ep1", i, "Em1", i+1
        ampo += 0.5*J4, "Em1", i, "Ep1", i+1
        ampo += 0.5*J5, "Ep2", i, "Em2", i+1
        ampo += 0.5*J5, "Em2", i, "Ep2", i+1
        ampo += 0.5*J6, "Ep3", i, "Em3", i+1
        ampo += 0.5*J6, "Em3", i, "Ep3", i+1
        ampo += 0.5*J7, "Ep4", i, "Em4", i+1
        ampo += 0.5*J7, "Em4", i, "Ep4", i+1
        ampo += 0.5*J8, "Ep5", i, "Em5", i+1
        ampo += 0.5*J8, "Em5", i, "Ep5", i+1
        ampo += 0.5*J9, "Ep6", i, "Em6", i+1
        ampo += 0.5*J9, "Em6", i, "Ep6", i+1
    end
    return MPO(ampo, sites)
end

# Hamiltonian
H = build_h_1d_chain(sites; J1=1.0, J2=1.0, J3=1.0, J4=1.0, J5=1.0, J6=1.0, J7=1.0, J8=1.0, J9=1.0)


# DMRG 
psi0 = randomMPS(sites)
sweeps = Sweeps(5)
setmaxdim!(sweeps, 10, 20, 50)
setcutoff!(sweeps, 1E-10)
energy, psi = dmrg(H, psi0, sweeps)
println("基态能量: ", energy)


