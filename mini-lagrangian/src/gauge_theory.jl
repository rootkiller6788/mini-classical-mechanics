# gauge_theory.jl — 规范对称性与拉格朗日力学
# 参考：Goldstein Ch.1, Peskin & Schroeder Ch.15
# L2 Core Concepts: gauge symmetry
# L8 Advanced Topics: gauge invariance in classical mechanics

"""
规范变换: A → A + ∇χ, φ → φ - ∂χ/∂t

电磁场在规范变换下不变:
  E = -∇φ - ∂A/∂t → E' = E
  B = ∇ × A → B' = B

拉格朗日量 L = ½mv² - q(φ - v·A) 在规范变换下改变一个全时间导数:
  L → L + q dχ/dt
因此运动方程不变 (作用量 S → S + q[χ(t₂) - χ(t₁)], δS 不变).
"""

# 规范变换
function gauge_transform_A(A::Function, chi::Function, r::Vector{Float64}, t::Float64; h=1e-6)
    grad_chi = zeros(3)
    for i in 1:3
        rp = copy(r); rp[i] += h
        rm = copy(r); rm[i] -= h
        grad_chi[i] = (chi(rp, t) - chi(rm, t))/(2h)
    end
    A_orig = A(r, t)
    return A_orig + grad_chi
end

function gauge_transform_phi(phi::Function, chi::Function, r::Vector{Float64}, t::Float64; h=1e-6)
    dchi_dt = (chi(r, t+h) - chi(r, t-h))/(2h)
    phi_orig = phi(r, t)
    return phi_orig - dchi_dt
end

"""验证电磁拉格朗日量在规范变换下的不变性: L' = L + q dχ/dt"""
function verify_gauge_invariance(q::Float64, chi::Function, r::Vector{Float64}, v::Vector{Float64}, t::Float64, phi::Function, A::Function; h=1e-6)
    m = 1.0
    phi_old = phi(r, t)
    A_old = A(r, t)
    L_old = 0.5*m*dot(v,v) - q*phi_old + q*dot(v, A_old)
    A_new = gauge_transform_A(A, chi, r, t; h=h)
    phi_new = gauge_transform_phi(phi, chi, r, t; h=h)
    L_new = 0.5*m*dot(v,v) - q*phi_new + q*dot(v, A_new)
    dchi_dt = (chi(r, t+h) - chi(r, t-h))/(2h)
    delta_L_theory = q * dchi_dt
    delta_L_actual = L_new - L_old
    return L_old, L_new, delta_L_actual, delta_L_theory
end

"""Chern-Simons 荷 (2+1D): j^0 = (1/2π) ε^{0μν} F_{μν}"""
function chern_simons_charge(Bz::Function, area_bounds::Vector{Float64}; n_grid=50)
    x_min, x_max, y_min, y_max = area_bounds
    dx = (x_max-x_min)/n_grid; dy = (y_max-y_min)/n_grid
    flux = 0.0
    for i in 1:n_grid, j in 1:n_grid
        x = x_min + (i-0.5)*dx; y = y_min + (j-0.5)*dy
        flux += Bz(x, y) * dx * dy
    end
    return flux
end

"""磁场从矢势计算: B = ∇ × A"""
function magnetic_field_from_A(A::Function, r::Vector{Float64}, t::Float64; h=1e-6)
    x,y,z = r[1], r[2], r[3]
    dAz_dy = (A([x,y+h,z], t)[3] - A([x,y-h,z], t)[3])/(2h)
    dAy_dz = (A([x,y,z+h], t)[2] - A([x,y,z-h], t)[2])/(2h)
    Bx = dAz_dy - dAy_dz
    dAx_dz = (A([x,y,z+h], t)[1] - A([x,y,z-h], t)[1])/(2h)
    dAz_dx = (A([x+h,y,z], t)[3] - A([x-h,y,z], t)[3])/(2h)
    By = dAx_dz - dAz_dx
    dAy_dx = (A([x+h,y,z], t)[2] - A([x-h,y,z], t)[2])/(2h)
    dAx_dy = (A([x,y+h,z], t)[1] - A([x,y-h,z], t)[1])/(2h)
    Bz = dAy_dx - dAx_dy
    return [Bx, By, Bz]
end

"""电场从势计算: E = -∇φ - ∂A/∂t"""
function electric_field_from_potentials(phi::Function, A::Function, r::Vector{Float64}, t::Float64; h=1e-6)
    x,y,z = r[1], r[2], r[3]
    dphi_dx = (phi([x+h,y,z], t) - phi([x-h,y,z], t))/(2h)
    dphi_dy = (phi([x,y+h,z], t) - phi([x,y-h,z], t))/(2h)
    dphi_dz = (phi([x,y,z+h], t) - phi([x,y,z-h], t))/(2h)
    Ex = -dphi_dx; Ey = -dphi_dy; Ez = -dphi_dz
    dA_dt = (A(r, t+h) - A(r, t-h))/(2h)
    return [Ex - dA_dt[1], Ey - dA_dt[2], Ez - dA_dt[3]]
end

"""Maxwell 场强张量 F_{μν} = ∂_μ A_ν - ∂_ν A_μ"""
function field_strength_tensor(A::Function, r::Vector{Float64}, t::Float64; h=1e-6)
    F = zeros(4, 4)
    E = electric_field_from_potentials(phi -> 0.0, A, r, t; h=h)
    B = magnetic_field_from_A(A, r, t; h=h)
    for i in 1:3
        F[i+1, 1] = E[i]; F[1, i+1] = -E[i]
    end
    F[2,3] = B[3]; F[3,2] = -B[3]
    F[3,4] = B[1]; F[4,3] = -B[1]
    F[4,2] = B[2]; F[2,4] = -B[2]
    return F
end

"""拉格朗日量中加全导数项: L -> L + dF/dt 不改变运动方程"""
function total_derivative_addition(L::Function, F::Function, q::Vector{Float64}, qdot::Vector{Float64}, t::Float64; h=1e-6)
    L_val = L(q, qdot)
    dF_dq = zeros(length(q))
    for i in 1:length(q)
        qp = copy(q); qp[i] += h; qm = copy(q); qm[i] -= h
        dF_dq[i] = (F(qp, t) - F(qm, t))/(2h)
    end
    dF_dt = (F(q, t+h) - F(q, t-h))/(2h)
    dF_dt_total = dot(dF_dq, qdot) + dF_dt
    return L_val + dF_dt_total
end

"""Wess-Zumino 项: L_WZ = θ qdot (拓扑项, 平庸运动方程)"""
function wess_zumino_lagrangian(theta::Float64, qdot::Float64)
    return theta * qdot
end

"""Dirac 磁单极子矢势 (奇异弦沿负 z 轴)"""
function dirac_monopole_A(g::Float64, r::Float64, theta::Float64, phi::Float64)
    sθ = sin(theta)
    if sθ < 1e-12
        return [0.0, 0.0, 0.0]
    end
    A_phi = g * (1 - cos(theta)) / (r * sθ)
    Ax = -A_phi * sin(phi)
    Ay = A_phi * cos(phi)
    Az = 0.0
    return [Ax, Ay, Az]
end

"""Berry 联络: A_n(R) = i ⟨n|∇_R n⟩"""
function berry_connection(state_func::Function, R::Vector{Float64}; h=1e-6)
    n = length(R)
    A = zeros(n)
    for i in 1:n
        Rp = copy(R); Rp[i] += h; Rm = copy(R); Rm[i] -= h
        psi_p = state_func(Rp); psi_m = state_func(Rm)
        psi_R = state_func(R)
        A[i] = imag(dot(psi_R, (psi_p - psi_m)/(2h)))
    end
    return A
end

"""Berry 曲率: Ω_{ij} = ∂_i A_j - ∂_j A_i"""
function berry_curvature(connection::Function, R::Vector{Float64}; h=1e-6)
    n = length(R)
    Omega = zeros(n, n)
    for i in 1:n, j in 1:n
        Rp = copy(R); Rp[j] += h; Rm = copy(R); Rm[j] -= h
        A_p = connection(Rp); A_m = connection(Rm)
        dAi_dRj = (A_p[i] - A_m[i])/(2h)
        # dAj/dRi
        Rpi = copy(R); Rpi[i] += h; Rmi = copy(R); Rmi[i] -= h
        Api = connection(Rpi); Ami = connection(Rmi)
        dAj_dRi = (Api[j] - Ami[j])/(2h)
        Omega[i,j] = dAi_dRj - dAj_dRi
    end
    return Omega
end

"""Aharonov-Bohm 相位: Δφ = (q/ħ) ∮ A·dl = 2π Φ/Φ₀"""
function aharonov_bohm_phase(q::Float64, A::Function, path::Vector{Vector{Float64}}; hbar=1.0)
    phase = 0.0
    for i in 1:length(path)-1
        r_mid = (path[i] + path[i+1]) / 2
        dl = path[i+1] - path[i]
        phase += q/hbar * dot(A(r_mid, 0.0), dl)
    end
    return phase
end
