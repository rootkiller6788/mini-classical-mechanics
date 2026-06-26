# functional.jl — 泛函分析基础: 泛函求值, Gateaux/Frechet导数, Beltrami恒等式
# 参考: Gelfand & Fomin "Calculus of Variations", Courant & Hilbert

using LinearAlgebra

# ============================================================
# 泛函求值
# ============================================================
"数值求泛函值: J[y] = ∫_a^b F(x,y(x),y''(x)) dx (梯形法则)"
function evaluate_functional(F::Function, y::Function, yp::Function, a::Float64, b::Float64; n=1000)
    dx=(b-a)/(n-1); J=0.0
    for i in 1:n; x=a+(i-1)*dx
        w=(i==1||i==n)?0.5:1.0; J+=w*F(x,y(x),yp(x))*dx
    end; J
end

"泛函范数: ‖y‖ = max|y| + max|y''| + √(∫y²dx)"
function functional_norm(y::Function, yp::Function, a::Float64, b::Float64; n=500)
    dx=(b-a)/(n-1); norm_int=0.0; max_y=0.0; max_yp=0.0
    for i in 1:n; x=a+(i-1)*dx
        yv,ypv=abs(y(x)),abs(yp(x))
        max_y=max(max_y,yv); max_yp=max(max_yp,ypv)
        w=(i==1||i==n)?0.5:1.0; norm_int+=w*yv^2*dx
    end; max_y+max_yp+sqrt(norm_int)
end

"泛函在两点间的增量 ΔJ = J[y+η] - J[y]"
function functional_increment(F::Function, y::Function, yp::Function, eta::Function, etap::Function, a::Float64, b::Float64; n=500)
    yt(x)=y(x)+eta(x); ypt(x)=yp(x)+etap(x)
    evaluate_functional(F,yt,ypt,a,b;n=n)-evaluate_functional(F,y,yp,a,b;n=n)
end

# ============================================================
# 导数: Gateaux, Frechet, 泛函导数
# ============================================================
"Gateaux 导数: δJ[y;η] = lim_{ε→0} (J[y+εη]-J[y])/ε"
function gateaux_derivative(F::Function, y::Function, yp::Function, eta::Function, etap::Function, a::Float64, b::Float64; eps_val=1e-6, n=1000)
    yf(x)=y(x)+eps_val*eta(x); ypf(x)=yp(x)+eps_val*etap(x)
    (evaluate_functional(F,yf,ypf,a,b;n=n)-evaluate_functional(F,y,yp,a,b;n=n))/eps_val
end

"泛函导数 (Euler-Lagrange 算子): δF/δy = ∂F/∂y - d/dx(∂F/∂y'')"
function functional_derivative(F::Function, dFdy::Function, dFdyp::Function, y::Function, yp::Function, a::Float64, b::Float64; n=1000)
    dx=(b-a)/(n-1); result=zeros(n); xs=range(a,b,length=n)
    for i in 1:n; x=xs[i]; yv=y(x); ypv=yp(x)
        if i==1; d=(dFdyp(x+dx,y(x+dx),yp(x+dx))-dFdyp(x,yv,ypv))/dx
        elseif i==n; d=(dFdyp(x,yv,ypv)-dFdyp(x-dx,y(x-dx),yp(x-dx)))/dx
        else; d=(dFdyp(x+dx,y(x+dx),yp(x+dx))-dFdyp(x-dx,y(x-dx),yp(x-dx)))/(2dx); end
        result[i]=dFdy(x,yv,ypv)-d
    end; xs, result
end

"Beltrami 恒等式: F - y''∂F/∂y'' = const (当 F 不显含 x)"
function beltrami_constant(F::Function, dFdyp::Function, y::Function, yp::Function, x::Float64)
    F(x,y(x),yp(x))-yp(x)*dFdyp(x,y(x),yp(x))
end

"验证 Beltrami 恒等式沿极值曲线是否守恒"
function verify_beltrami(F::Function, dFdyp::Function, y::Function, yp::Function, a::Float64, b::Float64; n=200)
    xs=range(a,b,length=n); vals=[beltrami_constant(F,dFdyp,y,yp,x) for x in xs]
    xs, vals, maximum(vals)-minimum(vals)
end

# ============================================================
# 弱极值与强极值
# ============================================================
"弱极值: C¹ 邻近函数中最小; 强极值: C⁰ 邻近中最小"
function weak_extremum_test(J::Function, y0::Function, perturbations::Vector{Function}, a::Float64, b::Float64; n=500)
    J0=J(y0,a,b;n=n); results=Bool[]
    for eta in perturbations
        yt(x)=y0(x)+eta(x)
        push!(results, J(yt,a,b;n=n) >= J0)
    end; all(results)
end

"比较两条曲线的泛函值"
function compare_functionals(F::Function, y1::Function, yp1::Function, y2::Function, yp2::Function, a::Float64, b::Float64; n=500)
    J1=evaluate_functional(F,y1,yp1,a,b;n=n)
    J2=evaluate_functional(F,y2,yp2,a,b;n=n)
    J1, J2, J2-J1, (J2-J1)/max(abs(J1),1e-15)
end

"一阶变分: δJ = ∫(∂F/∂y·η + ∂F/∂y''·η'')dx"
function first_variation_integral(dFdy::Function, dFdyp::Function, y::Function, yp::Function, eta::Function, etap::Function, a::Float64, b::Float64; n=500)
    dx=(b-a)/n; delta=0.0
    for i in 1:n; x=a+(i-0.5)*dx
        delta+=(dFdy(x,y(x),yp(x))*eta(x)+dFdyp(x,y(x),yp(x))*etap(x))*dx
    end; delta
end

"二阶变分公式: δ²J = ½∫[P η''² + Q η²]dx 其中 P=∂²F/∂y''², Q=∂²F/∂y²-d/dx(∂²F/∂y∂y'')"
function second_variation_formula(d2Fdyp2::Function, d2Fdy2::Function, eta::Function, etap::Function, a::Float64, b::Float64; n=500)
    dx=(b-a)/n; delta2=0.0
    for i in 1:n; x=a+(i-0.5)*dx
        delta2+=0.5*(d2Fdyp2(x)*etap(x)^2+d2Fdy2(x)*eta(x)^2)*dx
    end; delta2
end

# ============================================================
# 泛函的凸性
# ============================================================
"测试泛函是否凸: J[θy₁+(1-θ)y₂] ≤ θJ[y₁]+(1-θ)J[y₂]"
function is_convex_functional(F::Function, y1::Function, yp1::Function, y2::Function, yp2::Function, a::Float64, b::Float64; n=500, n_theta=5)
    J1=evaluate_functional(F,y1,yp1,a,b;n=n); J2=evaluate_functional(F,y2,yp2,a,b;n=n)
    convex=true
    for i in 1:n_theta
        theta=i/(n_theta+1)
        yt(x)=theta*y1(x)+(1-theta)*y2(x); ypt(x)=theta*yp1(x)+(1-theta)*yp2(x)
        Jt=evaluate_functional(F,yt,ypt,a,b;n=n)
        J_interp=theta*J1+(1-theta)*J2
        Jt > J_interp+1e-10 && (convex=false; break)
    end; convex
end

# ============================================================
# 泛函的 Frechet 导数 (强导数)
# ============================================================
"Frechet 导数: 对于线性泛函的范数逼近"
function frechet_derivative_norm(F::Function, y::Function, yp::Function, a::Float64, b::Float64; n=500, eps_vals=[1e-2,1e-3,1e-4])
    results=Float64[]
    for eps_val in eps_vals
        delta_max=0.0
        for _ in 1:10  # 随机扰动
            eta_func(x)=eps_val*randn()
            etap_func(x)=eps_val*randn()
            dJ=gateaux_derivative(F,y,yp,eta_func,etap_func,a,b;eps_val=eps_val,n=n)
            delta_max=max(delta_max,abs(dJ))
        end; push!(results,delta_max)
    end; eps_vals, results
end

export evaluate_functional, functional_norm, functional_increment
export gateaux_derivative, functional_derivative, beltrami_constant, verify_beltrami
export weak_extremum_test, compare_functionals
export first_variation_integral, second_variation_formula
export is_convex_functional, frechet_derivative_norm
