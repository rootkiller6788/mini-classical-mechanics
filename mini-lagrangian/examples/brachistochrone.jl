#!/usr/bin/env julia
# examples/brachistochrone.jl — 最速降线问题
# 变分法的经典：在重力场中找两点间的最速滑行路径

include("../src/Lagrangian.jl")
using .Lagrangian

function main()
    println("="^60)
    println("  Brachistochrone Problem — Calculus of Variations")
    println("="^60)

    g=9.81

    # 从 (0,0) 到 (x1,y1) 的最速降线
    # 时间: T = ∫ ds/v = ∫ √(1+y'²)/√(2gy) dx
    # Lagrangian: L(y, y') = √((1+y'²)/y)
    # 不显含 x → Beltrami identity: L - y'*∂L/∂y' = const

    # 解析解: 摆线
    # x(θ) = R(θ - sinθ), y(θ) = R(1 - cosθ)
    # 其中 R 由终点条件确定

    function cycloid_point(R, theta)
        x = R*(theta - sin(theta))
        y = R*(1 - cos(theta))
        return x, y
    end

    # 给定终点 (x1, y1)，求解 R 和 θ_end
    # 需要解: x1/R = θ1 - sin(θ1), y1/R = 1 - cos(θ1)
    # → (1-cosθ1)/(θ1-sinθ1) = y1/x1
    function find_cycloid_params(x1, y1; tol=1e-10, max_iter=100)
        ratio = y1 / x1
        # 二分法找 θ1
        lo, hi = 1e-6, 2π-1e-6
        for _ in 1:max_iter
            mid = (lo + hi) / 2
            val = (1 - cos(mid)) / (mid - sin(mid))
            if val < ratio; lo = mid; else; hi = mid; end
            if hi - lo < tol; break; end
        end
        theta1 = (lo + hi) / 2
        R = y1 / (1 - cos(theta1))
        return R, theta1
    end

    # 几个不同的终点
    test_points = [(1.0, 0.5), (2.0, 1.0), (1.0, 1.0), (0.5, 0.8)]
    
    println("Brachistochrone from (0,0) to various endpoints:\n")
    for (x1, y1) in test_points
        R, theta1 = find_cycloid_params(x1, y1)
        travel_time = theta1 * sqrt(R / g)  # T = θ1*√(R/g)
        straight_time = sqrt(2*(x1^2 + y1^2)/(g*y1/(sqrt(x1^2+y1^2))))
        
        println("  To ($x1, $y1):")
        println("    Cycloid: R=$(round(R,digits=4)), θ_end=$(round(rad2deg(theta1),digits=1))°")
        println("    Travel time: $(round(travel_time,digits=4)) s")
        println("    Straight line time: $(round(straight_time,digits=4)) s")
        println("    Speedup: $(round(straight_time/travel_time*100-100,digits=1))% faster")
    end

    # 用离散变分法近似验证
    println("\n--- Discrete variational check ---")
    x1, y1 = 1.0, 0.5
    R, theta1 = find_cycloid_params(x1, y1)
    
    n_segments = 50
    dt_theta = theta1 / n_segments
    T_disc = 0.0
    for i in 1:n_segments
        th0 = (i-1)*dt_theta; th1 = i*dt_theta
        x0,y0 = cycloid_point(R, th0); x1c,y1c = cycloid_point(R, th1)
        ds = sqrt((x1c-x0)^2 + (y1c-y0)^2)
        v_mid = sqrt(2*g*(y0+y1c)/2)
        T_disc += ds / v_mid
    end
    println("  Continuous theory: $(round(theta1*sqrt(R/g),digits=6)) s")
    println("  Discrete (n=50):    $(round(T_disc,digits=6)) s")

    println("\n✅ brachistochrone.jl done.")
end
main()
