# types.jl — 连续介质力学核心类型
# 参考: Landau Vol.7, Timoshenko, Fung

using LinearAlgebra

# ============================================================
# StressTensor — 应力张量
# ============================================================
struct StressTensor
    xx::Float64; yy::Float64; zz::Float64; xy::Float64; xz::Float64; yz::Float64
end
stress_matrix(s::StressTensor)=[s.xx s.xy s.xz; s.xy s.yy s.yz; s.xz s.yz s.zz]
StressTensor(M::AbstractMatrix)=StressTensor(M[1,1],M[2,2],M[3,3],M[1,2],M[1,3],M[2,3])

"静水应力: p = tr(σ)/3"
hydrostatic_pressure(s::StressTensor) = (s.xx+s.yy+s.zz)/3

"偏应力张量: sᵢⱼ = σᵢⱼ - pδᵢⱼ"
function deviatoric_stress(s::StressTensor)
    p=hydrostatic_pressure(s)
    StressTensor(s.xx-p, s.yy-p, s.zz-p, s.xy, s.xz, s.yz)
end

"von Mises from deviatoric: σ_vm=√(3/2 s:s)"
function von_mises_deviatoric(s::StressTensor)
    sdev=deviatoric_stress(s); sqrt(1.5*(sdev.xx^2+sdev.yy^2+sdev.zz^2+2(sdev.xy^2+sdev.xz^2+sdev.yz^2)))
end

"应力不变量"
function stress_invariants(sigma::StressTensor)
    s=principal_stresses(sigma)
    I1=s[1]+s[2]+s[3]
    J2=((s[1]-s[2])^2+(s[2]-s[3])^2+(s[3]-s[1])^2)/6
    J3=((2*s[1]-s[2]-s[3])*(2*s[2]-s[3]-s[1])*(2*s[3]-s[1]-s[2]))/27
    (I1=I1,J2=J2,J3=J3)
end

"主应力 (特征值, 降序)"
principal_stresses(sigma::StressTensor)=sort(eigvals(stress_matrix(sigma)),rev=true)

# ============================================================
# StrainTensor — 应变张量
# ============================================================
struct StrainTensor
    xx::Float64; yy::Float64; zz::Float64; xy::Float64; xz::Float64; yz::Float64
end
strain_matrix(eps::StrainTensor)=[eps.xx eps.xy eps.xz; eps.xy eps.yy eps.yz; eps.xz eps.yz eps.zz]
StrainTensor(M::AbstractMatrix)=StrainTensor(M[1,1],M[2,2],M[3,3],M[1,2],M[1,3],M[2,3])

"小变形应变张量: ε=½(∇u+∇uᵀ)"
strain_tensor(du::AbstractMatrix)=0.5*(du+transpose(du))

"体积应变: θ = tr(ε) = ε_xx+ε_yy+ε_zz"
dilatation(eps::StrainTensor)=eps.xx+eps.yy+eps.zz

"偏应变: eᵢⱼ = εᵢⱼ - θ/3 δᵢⱼ"
function deviatoric_strain(eps::StrainTensor)
    th=dilatation(eps)/3; StrainTensor(eps.xx-th,eps.yy-th,eps.zz-th,eps.xy,eps.xz,eps.yz)
end

"主应变"
principal_strains(eps::StrainTensor)=sort(eigvals(strain_matrix(eps)),rev=true)

"等效应变 (von Mises): ε_eq = √(2/3 e:e)"
function equivalent_strain(eps::StrainTensor)
    ed=deviatoric_strain(eps); sqrt(2/3*(ed.xx^2+ed.yy^2+ed.zz^2+2(ed.xy^2+ed.xz^2+ed.yz^2)))
end

# ============================================================
# ElasticMaterial
# ============================================================
struct ElasticMaterial
    E::Float64; nu::Float64; rho::Float64; G::Float64; K::Float64; lam::Float64; mu::Float64
end
function ElasticMaterial(E::Float64, nu::Float64, rho=0.0)
    G=E/(2*(1+nu)); K=E/(3*(1-2nu)); lam=E*nu/((1+nu)*(1-2nu)); mu=G
    ElasticMaterial(E,nu,rho,G,K,lam,mu)
end

# ============================================================
# DeformationGradient — 变形梯度 F_ij = partial x_i / partial X_j
# 极分解: F = R·U = V·R (R旋转, U/V伸缩)
# Green-Lagrange: E = 1/2(F^T F - I)
# Almansi: e = 1/2(I - F^{-T} F^{-1})
# ============================================================
struct DeformationGradient
    F11::Float64; F12::Float64; F13::Float64
    F21::Float64; F22::Float64; F23::Float64
    F31::Float64; F32::Float64; F33::Float64
end
DeformationGradient() = DeformationGradient(1,0,0,0,1,0,0,0,1)
function defgrad_to_matrix(F::DeformationGradient)
    [F.F11 F.F12 F.F13; F.F21 F.F22 F.F23; F.F31 F.F32 F.F33]
end
defgrad_det(F::DeformationGradient) = det(defgrad_to_matrix(F))

"右Cauchy-Green张量: C = F^T F"
function right_cauchy_green(F::DeformationGradient)
    m = defgrad_to_matrix(F); m'*m
end

"左Cauchy-Green(Finger)张量: b = F F^T"
function left_cauchy_green(F::DeformationGradient)
    m = defgrad_to_matrix(F); m*m'
end

"Green-Lagrange应变: E = 1/2(C - I)"
function green_lagrange_strain(F::DeformationGradient)
    C = right_cauchy_green(F); 0.5*(C - I)
end

"Almansi-Euler应变: e = 1/2(I - b^{-1})"
function almansi_strain(F::DeformationGradient)
    b = left_cauchy_green(F); 0.5*(I - inv(b))
end

"极分解: F = R·U (Newton迭代, 参考Higham 1986)"
function polar_decomposition(F::DeformationGradient)
    Fm = defgrad_to_matrix(F); C = Fm'*Fm
    Uk = Matrix{Float64}(I,3,3)
    for _ in 1:15; Uk = 0.5*(Uk + C/Uk); end
    R = Fm/Uk; R, Uk
end

"体积比: J = det(F), J>0 (不可穿透)"
volume_ratio(F::DeformationGradient) = defgrad_det(F)

"等容部分: F_bar = J^{-1/3} F (det F_bar = 1)"
function isochoric_part(F::DeformationGradient)
    J = defgrad_det(F); s = J^(-1/3)
    DeformationGradient(s*F.F11,s*F.F12,s*F.F13,s*F.F21,s*F.F22,s*F.F23,s*F.F31,s*F.F32,s*F.F33)
end

"Hencky(对数)应变: ln(V), V来自F=VR极分解"
function hencky_strain(F::DeformationGradient)
    _, V = polar_decomposition(F)
    # 近似: log(V) = logm(V) via eigen decomposition
    evals = eigvals(V); evecs = eigvecs(V)
    logV = evecs * Diagonal(log.(evals)) * evecs'
    StrainTensor(logV)
end

# ============================================================
# RotationTensor — 有限旋转张量
# 通过轴角表示: R = I + sin(θ)[a]× + (1-cos(θ))[a]×²
# ============================================================
struct RotationTensor
    angle::Float64; axis::Vector{Float64}  # rad, unit vector
end
function rotation_matrix(R::RotationTensor)
    a = R.axis; th = R.angle
    ax = [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
    I + sin(th)*ax + (1-cos(th))*ax*ax
end
"旋转张量合成 (Rodrigues参数)"
function compose_rotations(R1::RotationTensor, R2::RotationTensor)
    R1m = rotation_matrix(R1); R2m = rotation_matrix(R2)
    Rm = R1m * R2m
    th = acos((tr(Rm)-1)/2)
    ax = [R1m[3,2]-R1m[2,3], R1m[1,3]-R1m[3,1], R1m[2,1]-R1m[1,2]]/(2*sin(th))
    RotationTensor(th, ax/norm(ax))
end

# ============================================================
# MaterialSymmetry — 弹性对称性群
# ============================================================
@enum MaterialSymmetry begin
    SYM_ISOTROPIC
    SYM_TRANSVERSELY_ISO
    SYM_ORTHOTROPIC
    SYM_MONOCLINIC
    SYM_TRICLINIC
end
"独立弹性常数数量"
n_independent_constants(::Type{Val{:isotropic}}) = 2
n_independent_constants(::Type{Val{:trans_iso}}) = 5
n_independent_constants(::Type{Val{:orthotropic}}) = 9
n_independent_constants(::Type{Val{:monoclinic}}) = 13
n_independent_constants(::Type{Val{:triclinic}}) = 21

# ============================================================
# RateOfDeformation — 变形率张量 D = sym(L), L = dot(F)F^{-1}
# 旋率张量: W = skew(L)
# ============================================================
"变形率: D = 1/2(L + L^T), L = velocity_gradient"
function rate_of_deformation(L::AbstractMatrix)
    0.5*(L + L')
end
"旋率: W = 1/2(L - L^T)"
function spin_tensor(L::AbstractMatrix)
    0.5*(L - L')
end

# ============================================================
# 材料客观性 (Objectivity) 检查
# ============================================================
"检查张量在刚体旋转下是否客观: Q T Q^T (适用于Cauchy应力)"
function objective_transform(T::AbstractMatrix, Q::AbstractMatrix)
    Q * T * Q'
end

export StressTensor, StrainTensor, ElasticMaterial, DeformationGradient, RotationTensor
export stress_matrix, strain_matrix, strain_tensor
export hydrostatic_pressure, deviatoric_stress, von_mises_deviatoric, stress_invariants, principal_stresses
export dilatation, deviatoric_strain, principal_strains, equivalent_strain
export defgrad_to_matrix, defgrad_det, right_cauchy_green, left_cauchy_green
export green_lagrange_strain, almansi_strain, polar_decomposition, volume_ratio
export isochoric_part, hencky_strain
export rotation_matrix, compose_rotations
export MaterialSymmetry, n_independent_constants
export rate_of_deformation, spin_tensor, objective_transform
