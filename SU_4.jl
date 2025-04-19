using ITensors, ITensorMPS

global Dim_site = 6  ## Define the dimension of Hilbert space in a site.
global Lx = 2  ## Define the length of the chain
global Ly = 3 ## Define the width of the chain
global n_in_cell = 2 ## Define the number of sublattices in a unit cell
function ITensors.space(::SiteType"IdxSU4_6d")
    return Dim_site
end

## Define the corresponding Operator of other generators. There should be 15 in total.
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
 #Define the sites
 N = Lx * Ly * n_in_cell
 sites = siteinds("IdxSU4_6d", N)

function pos_to_idx(cell_pos,sublat_idx,Lx,Ly,n_in_cell)
  # Convert the cell position and sublattice index to a linear index
  pbc_pos = [mod1(cell_pos[1],Lx), mod1(cell_pos[2],Ly)]
  return (pbc_pos[1]-1)*Ly*n_in_cell + (pbc_pos[2]-1)*n_in_cell + sublat_idx
end

# Define the triangular lattice in open boundary condition
Lattice_triangle = Vector{Tuple{Int, Int}}()
for x in 1:Lx, y in 1:Ly, sublat_idx in 1:n_in_cell
    i  = pos_to_idx([x,y], sublat_idx, Lx, Ly, n_in_cell)

    ## right 
    if x < Lx
      j  = pos_to_idx([x+1,y], sublat_idx, Lx, Ly, n_in_cell)
      push!(Lattice_triangle, (i,j))
    end

    ## up
    j  = pos_to_idx([x,y+1], sublat_idx, Lx, Ly, n_in_cell)
    push!(Lattice_triangle, (i,j))

    ## up left
    if x >1 
        j  = pos_to_idx([x-1,y+1], sublat_idx, Lx, Ly, n_in_cell)
        push!(Lattice_triangle, (i,j))
    end
end
#  AutoMPO 
function build_h_2d(sites; Lattice, Lx, Ly, n_in_cell, J1=1.0, J2=1.0, J3=1.0, J4=1.0, J5=1.0, J6=1.0, J7=1.0, J8=1.0, J9=1.0)
    @assert N == Lx * Ly * n_in_cell

    ampo = AutoMPO()
    for (i,j) in Lattice  
        # H
        ampo += J1, "H1", i, "H1", j
        ampo += J2, "H2", i, "H2", j
        ampo += J3, "H3", i, "H3", j
        # E
        ampo += 0.5*J4, "Ep1", i, "Em1", j
        ampo += 0.5*J4, "Em1", i, "Ep1", j
        ampo += 0.5*J5, "Ep2", i, "Em2", j
        ampo += 0.5*J5, "Em2", i, "Ep2", j
        ampo += 0.5*J6, "Ep3", i, "Em3", j
        ampo += 0.5*J6, "Em3", i, "Ep3", j
        ampo += 0.5*J7, "Ep4", i, "Em4", j
        ampo += 0.5*J7, "Em4", i, "Ep4", j
        ampo += 0.5*J8, "Ep5", i, "Em5", j
        ampo += 0.5*J8, "Em5", i, "Ep5", j
        ampo += 0.5*J9, "Ep6", i, "Em6", j
        ampo += 0.5*J9, "Em6", i, "Ep6", j
    end
    return MPO(ampo, sites)
end

# Hamiltonian
H = build_h_2d(sites; Lattice=Lattice_triangle, Lx=Lx, Ly=Ly, n_in_cell=2,J1=1.0, J2=1.0, J3=1.0, J4=1.0, J5=1.0, J6=1.0, J7=1.0, J8=1.0, J9=1.0)

# DMRG 
psi0 = randomMPS(sites)
sweeps = Sweeps(10)
setmaxdim!(sweeps, 10, 20, 50)
setcutoff!(sweeps, 1E-15)
energy, psi = dmrg(H, psi0, sweeps)
println("基态能量: ", energy)

