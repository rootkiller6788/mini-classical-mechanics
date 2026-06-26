# hamilton_jacobi.jl — Hamilton-Jacobi 方程
# 参考：Goldstein Ch.10, Arnold Ch.9

"""
Hamilton-Jacobi 方程: H(q, ∂S/∂q, t) + ∂S/∂t = 0
S(q, α, t) 是 Hamilton 主函数

如果 H 不显含 t: S(q, α, t) = W(q, α) - Et
→ 时间无关 H-J 方程: H(q, ∂W/∂q) = E

分离变量法适用于可积系统。
"""

"""
一维 H-J 求解器
W(q, E) = ∫ p(q, E) dq, 其中 p(q, E) 从 H(q, p) = E 解出
"""
function hamilton_jacobi_1d(H_func::Function, q_range::Tuple{Float64,Float64}, E::Float64, n_points::Int=1000)
    q_min, q_max = q_range
    dq = (q_max - q_min) / (n_points - 1)
    W_vals = Float64[]
    q_vals = Float64[]
    W = 0.0
    for i in 1:n_points
        q = q_min + (i-1)*dq
        # 从 H(q, p) = E 解出 p（取正根）
        # 用二分法解 p
        p = solve_p_from_H(H_func, q, E)
        push!(q_vals, q); push!(W_vals, W)
        W += p * dq
    end
    return q_vals, W_vals
end

"""
二分法从 H(q, p) = E 解 p
"""
function solve_p_from_H(H_func::Function, q::Float64, E::Float64; p_max=100.0, tol=1e-10)
    # 求 H(q, 0) vs E 判断是否有解
    H0 = H_func([q], [0.0])    
    if H0 > E
        return 0.0  # 经典禁区
    end
    lo, hi = 0.0, p_max
    for _ in 1:60
        mid = (lo + hi) / 2
        H_val = H_func([q], [mid])
        if H_val < E; lo = mid; else; hi = mid; end
        if hi - lo < tol; break; end
    end
    return (lo + hi) / 2
end

"""
谐振子 H-J 解: H = p²/(2m) + ½mω²q² = E
p = √(2mE - m²ω²q²)
W(q, E) = ∫ √(2mE - m²ω²q²) dq
→ 解析: W = (E/ω)*arcsin(q/a) + (q/2)*√(a²-q²), a=√(2E/(mω²))
"""
function harmonic_oscillator_HJ(m::Float64, omega::Float64, E::Float64, q_vals::Vector{Float64})
    a = sqrt(2*E/(m*omega^2))
    W_vals = Float64[]
    for q in q_vals
        if abs(q) > a
            push!(W_vals, NaN); continue
        end
        W = (E/omega)*asin(q/a) + (q/2)*sqrt(a^2 - q^2)
        push!(W_vals, W)
    end
    return W_vals
end

"""
作用量-角变量（从 H-J 角度）
J = ∂W/∂(action_param)
对于一维周期运动:
J(E) = (1/2π) ∮ p dq = (1/π) ∫_{q_min}^{q_max} √(2m(E-U(q))) dq
"""
function action_from_HJ(U::Function, m::Float64, E::Float64, q_range::Tuple{Float64,Float64}; n_points=2000)
    q_min, q_max = q_range
    dq = (q_max - q_min) / (n_points - 1)
    J = 0.0
    for i in 1:n_points
        q = q_min + (i-1)*dq
        U_val = U(q)
        if U_val < E
            p = sqrt(2*m*(E - U_val))
            w = (i==1 || i==n_points) ? 0.5 : 1.0
            J += w * p * dq
        end
    end
    return J / π
end

"""
H-J 方程的可分离性检测
对于 H = Σ H_i(q_i, p_i)，如果每个 H_i 只依赖一对 (q_i, p_i)，
则 H-J 方程可分离: S = Σ S_i(q_i)
"""
function is_separable(H_func::Function, n::Int)
    # 试探性检测：检查 ∂²H/∂q_i∂p_j = 0 for i≠j
    q_test = ones(n); p_test = ones(n)
    eps_val = 1e-5
    for i in 1:n, j in 1:n
        if i == j; continue; end
        q_ip = copy(q_test); q_ip[i] += eps_val
        p_jp = copy(p_test); p_jp[j] += eps_val
        H_pp = H_func(q_ip, p_jp)
        H_p = H_func(q_ip, p_test)
        H__p = H_func(q_test, p_jp)
        H__ = H_func(q_test, p_test)
        cross = (H_pp - H_p - H__p + H__) / eps_val^2
        if abs(cross) > 1e-6
            return false
        end
    end
    return true
end
