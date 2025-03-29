using ITensors, ITensorMPS


global Dim_site = 6  ## Define the dimension of Hilbert space in a site.
Nx = 2        # Length of x direction of lattice
Ny = 3        # Length of y direction of lattice
N = Nx*Ny

function ITensors.space(::SiteType"IdxSO6")
    return Dim_site
end

## Define the operator
function ITensors.op!(Op::ITensor,::OpName"I",::SiteType"IdxSO6",s::Index)
    Op[s'=>1,s=>1] = 1
    Op[s'=>2,s=>2] = 1
    Op[s'=>3,s=>3] = 1
    Op[s'=>4,s=>4] = 1
    Op[s'=>5,s=>5] = 1
    Op[s'=>6,s=>6] = 1
end


## Define the corresponding Operator of other generators. There should be 15 in total.





#@show IdxSO6

s = Index(6,"IdxSO6")
@show Op = op("I",s)

# Define the sites
sites = siteinds("IdxSO6",N)


## Construct the MPO
# 使用AutoMPO自动构建MPO
ampo = AutoMPO()
for i in 1:N-1
  for a in 1:15  # 遍历所有15个生成元
    ampo += (J, "T$a", i, "T$a", i+1)
  end
end

# 将AutoMPO转换为MPO
H = MPO(ampo, sites)

# 4. 验证关键属性
# ------------------------------
@show maxlinkdim(H)  # 检查MPO的键维度
@show inner(H[1], H[1])  # 验证算符归一化

