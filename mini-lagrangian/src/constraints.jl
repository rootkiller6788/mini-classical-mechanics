# constraints.jl — 拉格朗日乘子法与约束系统
# 参考：Goldstein Ch.2, MIT 8.012 Lecture 21

"""
拉格朗日乘子法处理完整约束 f(q,t) = 0
修正拉格朗日量: L' = L + λ*f(q,t)
扩展EL方程: d/dt(∂L/∂qdot) - ∂L/∂q = λ * ∂f/∂q + f(q,t) = 0

非完整约束 f(q, qdot, t) = 0 同理，∂f/∂q 换成 ∂f/∂qdot
"""

struct Constraint
    name::String
    func::Function          # f(q, t) = 0 或 f(q, qdot, t) = 0
    grad_q::Function        # ∂f/∂q
    grad_qdot::Function     # ∂f/∂qdot（完整约束时为 nothing）
    is_holonomic::Bool
end

"""
完整约束: f(q, t) = 0
例：质点被约束在半径为 R 的球面上: r² - R² = 0
"""
function holonomic_constraint(name::String, f::Function, grad_f::Function)
    return Constraint(name, f, grad_f, (q,qdot,t) -> zeros(length(q)), true)
end

"""
非完整约束: f(q, qdot, t) = 0
例：冰刀/滚轮约束: dx*sin(θ) - dy*cos(θ) = 0
"""
function nonholonomic_constraint(name::String, f::Function, grad_q::Function, grad_qdot::Function)
    return Constraint(name, f, grad_q, grad_qdot, false)
end

"""
含拉格朗日乘子的扩展 ODE 系统
状态: [q; qdot; λ]
EL: d/dt(∂L/∂qdot) - ∂L/∂q = Σ λ_k * ∂f_k/∂q
约束: f_k(q, qdot, t) = 0
"""
struct ConstrainedLagrangianSystem
    base_sys::Any             # EulerLagrangeSystem
    constraints::Vector{Constraint}
    n_constraints::Int
end

"""
将约束拉格朗日系统转换为一阶 ODE
"""
function constrained_el_to_ode(sys::ConstrainedLagrangianSystem)
    base = sys.base_sys
    m = sys.n_constraints
    n = base.n_dof

    function f(t, y)
        q = y[1:n]
        qdot = y[n+1:2n]
        lambda = y[2n+1:2n+m]

        # 基础 EL 力
        dL_dq = base.grad_L_q(q, qdot)
        dL_dqdot = base.grad_L_qdot(q, qdot)
        M = base.mass_matrix(q)

        # Mdot 数值估算
        Mdot = zeros(n, n)
        eps_val = 1e-6
        for i in 1:n, j in 1:n
            for k in 1:n
                qp = copy(q); qp[k] += eps_val
                qm = copy(q); qm[k] -= eps_val
                Mdot[i,j] += (base.mass_matrix(qp)[i,j] - base.mass_matrix(qm)[i,j]) / (2*eps_val) * qdot[k]
            end
        end

        # 约束力: Σ λ_k * grad f_k
        constraint_force = zeros(n)
        for k in 1:m
            cons = sys.constraints[k]
            grad = cons.is_holonomic ? cons.grad_q(q, 0.0) : cons.grad_q(q, qdot, t)
            constraint_force += lambda[k] * grad
        end

        rhs = dL_dq + constraint_force - Mdot * qdot
        qddot = M \ rhs

        # 约束方程的时间导数（指数稳定化）
        # d/dt(f_k) + α*f_k = 0 → 稳定到约束面
        alpha = 10.0  # 稳定化参数
        constraint_rhs = zeros(m)
        for k in 1:m
            cons = sys.constraints[k]
            if cons.is_holonomic
                # df/dt = grad_f · qdot → 加稳定项
                grad = cons.grad_q(q, t)
                constraint_rhs[k] = -dot(grad, qdot) - alpha * cons.func(q, t)
            else
                grad = cons.grad_q(q, qdot, t)
                grad_d = cons.grad_qdot(q, qdot, t)
                constraint_rhs[k] = -dot(grad, qdot) - alpha * cons.func(q, qdot, t)
            end
        end

        return [qdot; qddot; constraint_rhs]
    end
    return f
end

"""
常见约束系统

1. 质点约束在球面上: x² + y² + z² - R² = 0
"""
function spherical_constraint(R::Float64)
    function f(q, t) 
        return q[1]^2 + q[2]^2 + q[3]^2 - R^2
    end
    function grad_f(q, t)
        return [2q[1], 2q[2], 2q[3]]
    end
    return holonomic_constraint("sphere_R$(R)", f, grad_f)
end

"""
2. 质点约束在曲线上: y = sin(x)
"""
function curve_constraint(curve_func::Function, curve_grad::Function)
    function f(q, t)
        return q[2] - curve_func(q[1])
    end
    function grad_f(q, t)
        return [-curve_grad(q[1]), 1.0]
    end
    return holonomic_constraint("curve", f, grad_f)
end

"""
3. 两个质点通过刚性杆连接: |r1 - r2|² - L² = 0
"""
function rigid_rod_constraint(L::Float64)
    function f(q, t)
        dx = q[1] - q[4]; dy = q[2] - q[5]; dz = q[3] - q[6]
        return dx^2 + dy^2 + dz^2 - L^2
    end
    function grad_f(q, t)
        dx = q[1] - q[4]; dy = q[2] - q[5]; dz = q[3] - q[6]
        return [2dx, 2dy, 2dz, -2dx, -2dy, -2dz]
    end
    return holonomic_constraint("rigid_rod_L$(L)", f, grad_f)
end

"""
4. 纯滚动约束（非完整）: dx - R*cos(θ)*dφ = 0, dy - R*sin(θ)*dφ = 0
圆盘在平面上无滑滚动
"""
function rolling_disk_constraints(R::Float64)
    # q = [x, y, θ, φ]  (位置x, 位置y, 朝向角, 自转角)
    function f1(q, qdot, t)
        return qdot[1] - R * cos(q[3]) * qdot[4]
    end
    function f2(q, qdot, t)
        return qdot[2] - R * sin(q[3]) * qdot[4]
    end
    function grad_q1(q, qdot, t)
        return [0.0, 0.0, R*sin(q[3])*qdot[4], 0.0]
    end
    function grad_qdot1(q, qdot, t)
        return [1.0, 0.0, 0.0, -R*cos(q[3])]
    end
    function grad_q2(q, qdot, t)
        return [0.0, 0.0, -R*cos(q[3])*qdot[4], 0.0]
    end
    function grad_qdot2(q, qdot, t)
        return [0.0, 1.0, 0.0, -R*sin(q[3])]
    end
    return [
        nonholonomic_constraint("rolling_x", f1, grad_q1, grad_qdot1),
        nonholonomic_constraint("rolling_y", f2, grad_q2, grad_qdot2),
    ]
end
