# generalized.jl — 广义坐标与坐标变换
# 参考：Goldstein Ch.1-2

struct GeneralizedCoords
    n::Int
    values::Vector{Float64}
end
GeneralizedCoords(v::Vector{Float64}) = GeneralizedCoords(length(v), v)

struct GeneralizedVelocities
    n::Int
    values::Vector{Float64}
end
GeneralizedVelocities(v::Vector{Float64}) = GeneralizedVelocities(length(v), v)

struct GeneralizedState
    t::Float64
    q::GeneralizedCoords
    qdot::GeneralizedVelocities
end

"""极坐标 (r,theta) -> 笛卡尔"""
function polar_to_cartesian(r::Float64, theta::Float64)
    return r*cos(theta), r*sin(theta)
end

function polar_jacobian(r::Float64, theta::Float64)
    return [cos(theta) -r*sin(theta); sin(theta) r*cos(theta)]
end

"""球坐标 (r,theta,phi) -> 笛卡尔"""
function spherical_to_cartesian(r::Float64, theta::Float64, phi::Float64)
    return r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)
end

function spherical_jacobian(r::Float64, theta::Float64, phi::Float64)
    st, ct = sin(theta), cos(theta)
    sp, cp = sin(phi), cos(phi)
    return [st*cp r*ct*cp -r*st*sp; st*sp r*ct*sp r*st*cp; ct -r*st 0.0]
end

"""质量矩阵 M(q): T = 0.5*qdot'*M*qdot"""
function mass_matrix(masses::Vector{Float64}, jacobian::Matrix{Float64})
    n_dof = size(jacobian, 2)
    n_particles = length(masses)
    M = zeros(n_dof, n_dof)
    for i in 1:n_dof, j in 1:n_dof
        for k in 1:n_particles
            M[i,j] += masses[k] * (
                jacobian[3*k-2,i]*jacobian[3*k-2,j] +
                jacobian[3*k-1,i]*jacobian[3*k-1,j] +
                jacobian[3*k,  i]*jacobian[3*k,  j])
        end
    end
    return M
end

"""广义动量 p = M*qdot"""
generalized_momenta(M::Matrix{Float64}, qdot::Vector{Float64}) = M * qdot

## ============================================================
## 更多坐标系统
## ============================================================

"""柱坐标 (rho, phi, z) → 笛卡尔"""
function cylindrical_to_cartesian(rho::Float64, phi::Float64, z::Float64)
    return rho*cos(phi), rho*sin(phi), z
end

function cylindrical_jacobian(rho::Float64, phi::Float64)
    cp, sp = cos(phi), sin(phi)
    return [cp -rho*sp 0.0; sp rho*cp 0.0; 0.0 0.0 1.0]
end

"""椭圆坐标 (mu, nu): x=c*cosh(mu)*cos(nu), y=c*sinh(mu)*sin(nu)"""
function elliptic_to_cartesian(c::Float64, mu::Float64, nu::Float64)
    return c*cosh(mu)*cos(nu), c*sinh(mu)*sin(nu)
end

function elliptic_jacobian(c::Float64, mu::Float64, nu::Float64)
    return [c*sinh(mu)*cos(nu) -c*cosh(mu)*sin(nu);
            c*cosh(mu)*sin(nu)  c*sinh(mu)*cos(nu)]
end

"""双摆笛卡尔坐标（从广义坐标 θ1,θ2）"""
function double_pendulum_cartesian(L1::Float64, L2::Float64, q::Vector{Float64})
    th1, th2 = q[1], q[2]
    x1 = L1*sin(th1);               y1 = -L1*cos(th1)
    x2 = L1*sin(th1)+L2*sin(th2);   y2 = -L1*cos(th1)-L2*cos(th2)
    return [x1,y1], [x2,y2]
end

function double_pendulum_jacobian(L1::Float64, L2::Float64, q::Vector{Float64})
    th1, th2 = q[1], q[2]
    # 4行2列: [∂x1/∂θ1 ∂x1/∂θ2; ∂y1/∂θ1 ∂y1/∂θ2; ∂x2/∂θ1 ∂x2/∂θ2; ∂y2/∂θ1 ∂y2/∂θ2]
    return [L1*cos(th1) 0.0;
            L1*sin(th1) 0.0;
            L1*cos(th1) L2*cos(th2);
            L1*sin(th1) L2*sin(th2)]
end

## ============================================================
## Christoffel 符号（质量矩阵的度量联络）
## ============================================================

"""
Christoffel 符号 Γ^i_{jk} = 0.5 Σ_l M^{il}(∂M_{lj}/∂q_k + ∂M_{lk}/∂q_j - ∂M_{jk}/∂q_l)
用于 EL 方程的几何表述: qddot^i + Γ^i_{jk} qdot^j qdot^k = (M^{-1} ∂L/∂q)^i
"""
function christoffel_symbols(M_func::Function, q::Vector{Float64}; eps_val=1e-6)
    n = length(q)
    M = M_func(q)
    Minv = inv(M)
    Gamma = zeros(n, n, n)
    for i in 1:n, j in 1:n, k in 1:n
        # ∂M_{lj}/∂q_k
        dM_dq = zeros(n, n, n)
        for a in 1:n
            qp = copy(q); qp[a] += eps_val; qm = copy(q); qm[a] -= eps_val
            Mp = M_func(qp); Mm = M_func(qm)
            for l in 1:n
                dM_dq[l,:,a] = (Mp[l,:] - Mm[l,:]) / (2*eps_val)
            end
        end
        for l in 1:n
            Gamma[i,j,k] += 0.5 * Minv[i,l] * (dM_dq[l,j,k] + dM_dq[l,k,j] - dM_dq[j,k,l])
        end
    end
    return Gamma
end

"""
测地线方程（无势时的 EL 方程）: qddot^i + Γ^i_{jk} qdot^j qdot^k = 0
"""
function geodesic_acceleration(M_func::Function, q::Vector{Float64}, qdot::Vector{Float64})
    n = length(q)
    M = M_func(q)
    Gamma = christoffel_symbols(M_func, q)
    qddot = zeros(n)
    for i in 1:n
        for j in 1:n, k in 1:n
            qddot[i] -= Gamma[i,j,k] * qdot[j] * qdot[k]
        end
    end
    Minv = inv(M)
    return Minv \ qddot  # 实际上 qddot 已经包含了 M^{-1}
end

## ============================================================
## 常用拉格朗日量的动能部分
## ============================================================

"""双摆质量矩阵 M(θ1,θ2)"""
function double_pendulum_mass_matrix(m1::Float64, m2::Float64, L1::Float64, L2::Float64, q::Vector{Float64})
    th1, th2 = q[1], q[2]; cd = cos(th1 - th2)
    M11 = (m1+m2)*L1^2
    M12 = m2*L1*L2*cd
    M22 = m2*L2^2
    return [M11 M12; M12 M22]
end

"""球摆质量矩阵（广义坐标 θ, φ）"""
function spherical_pendulum_mass_matrix(m::Float64, L::Float64, q::Vector{Float64})
    th = q[1]
    return [m*L^2 0.0; 0.0 m*L^2*sin(th)^2]
end

"""中心力问题质量矩阵（极坐标 r, θ）"""
function central_force_mass_matrix(m::Float64, q::Vector{Float64})
    return [m 0.0; 0.0 m*q[1]^2]
end

"""带电粒子在磁场中的质量矩阵 + 规范势贡献
L = 0.5*m*v² + q*v·A → ∂L/∂v = m*v + q*A
质量矩阵部分始终是 m*I"""
function em_mass_matrix(m::Float64, n::Int, q_coords::Vector{Float64})
    return m * Matrix{Float64}(I, n, n)
end
