# second_variation.jl — 二阶变分理论: Legendre条件, Jacobi条件, Weierstrass, Hilbert
# 参考: Gelfand & Fomin Ch.5, Bolza "Lectures on the Calculus of Variations"

using LinearAlgebra

# ============================================================
# Legendre 条件 (二阶必要条件)
# ============================================================
"Legendre 条件: ∂²F/∂y''² ≥ 0 (弱必要条件)"
legendre_condition(d2Fdyp2::Float64) = d2Fdyp2 >= -1e-12

"增强 Legendre 条件: ∂²F/∂y''² > 0 (强必要条件)"
strong_legendre(d2Fdyp2::Float64) = d2Fdyp2 > 1e-12

"Legendre 条件沿极值曲线验证"
function verify_legendre_along_extremal(d2Fdyp2::Function, y_sol::Function, yp_sol::Function, a::Float64, b::Float64; n=200)
    xs=range(a,b,length=n); vals=[d2Fdyp2(x,y_sol(x),yp_sol(x)) for x in xs]
    all(v->v>=-1e-10, vals), minimum(vals), xs, vals
end

# ============================================================
# Jacobi 条件 (共轭点)
# ============================================================
"找共轭点: 解 Jacobi 方程 (P u''''+ Q u = 0), 其中 P=∂²F/∂y''², Q=∂²F/∂y²-d/dx(∂²F/∂y∂y'')"
function find_conjugate_points(a::Float64, b::Float64, P::Function, Q::Function; n=1000)
    dx=(b-a)/(n-1); u=zeros(n); du=zeros(n)
    u[1]=0.0; du[1]=1.0
    for i in 1:n-1; xm=a+(i-0.5)*dx
        du[i+1]=du[i]+dx*(-Q(xm)*u[i])
        u[i+1]=u[i]+dx*(du[i]/max(P(xm),1e-15))
    end
    pts=Float64[]
    for i in 2:n; u[i-1]*u[i]<0 && push!(pts,a+(i-1)*dx-dx*abs(u[i-1])/(abs(u[i-1])+abs(u[i]))); end
    pts
end

"Jacobi 条件: (a,b) 内无共轭点"
function jacobi_condition(pts::Vector{Float64}, a::Float64, b::Float64)
    isempty(filter(x->a<x<b, pts))
end

"共轭点最近的 x 坐标"
function first_conjugate_point(pts::Vector{Float64}, a::Float64)
    filtered=filter(x->x>a, pts)
    isempty(filtered) ? Inf : minimum(filtered)
end

# ============================================================
# Weierstrass 条件
# ============================================================
"Weierstrass 过剩函数: E(x,y,y''_opt, y''_test) = F(x,y,y''_test)-F(x,y,y''_opt)-(y''_test-y''_opt)∂F/∂y''|_{opt}"
function weierstrass_excess(F::Function, dFdyp::Function, x::Float64, y::Float64, yp_opt::Float64, yp_test::Float64)
    F(x,y,yp_test)-F(x,y,yp_opt)-(yp_test-yp_opt)*dFdyp(x,y,yp_opt)
end

"Weierstrass 必要条件: E(x,y,p_opt, p) ≥ 0 ∀p"
function weierstrass_condition(F::Function, dFdyp::Function, x::Float64, y::Float64, yp_opt::Float64, p_range::AbstractRange)
    all(weierstrass_excess(F,dFdyp,x,y,yp_opt,p) >= -1e-10 for p in p_range)
end

"沿极值曲线验证 Weierstrass 条件"
function verify_weierstrass_along_extremal(F::Function, dFdyp::Function, y_sol::Function, yp_sol::Function, a::Float64, b::Float64, p_range::AbstractRange; n=100)
    xs=range(a,b,length=n)
    all(weierstrass_condition(F,dFdyp,x,y_sol(x),yp_sol(x),p_range) for x in xs)
end

# ============================================================
# Hilbert 不变积分
# ============================================================
"Hilbert 不变积分: ∫[F(x,y,p)+(y''-p)∂F/∂p]dx (与路径无关)"
function hilbert_invariant_integral(F::Function, dFdyp::Function, y::Function, yp::Function, p_field::Function, a::Float64, b::Float64; n=500)
    dx=(b-a)/n; I=0.0
    for i in 1:n; x=a+(i-0.5)*dx
        p=p_field(x,y(x))
        I+=(F(x,y(x),p)+(yp(x)-p)*dFdyp(x,y(x),p))*dx
    end; I
end

"构造斜率场 p(x,y) 满足可积条件 (Hamilton-Jacobi)"
function slope_field_from_extremal(y_sol::Function, yp_sol::Function, a::Float64, b::Float64; n=100)
    xs=range(a,b,length=n)
    x->begin
        idx=argmin(abs.(xs.-x))
        yp_sol(xs[max(1,min(n,idx))])
    end
end

# ============================================================
# 充分条件
# ============================================================
"强极小值的充分条件: Legendre + Jacobi + Weierstrass"
function sufficient_conditions_strong(d2Fdyp2::Function, y_sol::Function, yp_sol::Function, P::Function, Q::Function, F::Function, dFdyp::Function, a::Float64, b::Float64; p_range=-10:0.1:10, n=200)
    leg_ok,_,_,_=verify_legendre_along_extremal(d2Fdyp2,y_sol,yp_sol,a,b;n=n)
    cj_pts=find_conjugate_points(a,b,P,Q;n=500)
    jac_ok=jacobi_condition(cj_pts,a,b)
    wei_ok=verify_weierstrass_along_extremal(F,dFdyp,y_sol,yp_sol,a,b,p_range;n=n)
    (Legendre=leg_ok, Jacobi=jac_ok, Weierstrass=wei_ok, Sufficient=leg_ok&&jac_ok&&wei_ok)
end

# ============================================================
# Ritz 方法 (直接法)
# ============================================================
"Ritz 法: 用有限基函数近似求解变分问题"
function ritz_method(F::Function, basis::Vector{Function}, basis_deriv::Vector{Function}, a::Float64, b::Float64; n=500, lambda_reg=0.0)
    m=length(basis); A=zeros(m,m); b_vec=zeros(m); dx=(b-a)/n
    for k in 1:n; x=a+(k-0.5)*dx
        for i in 1:m, j in 1:m; A[i,j]+=basis_deriv[i](x)*basis_deriv[j](x)*dx; end
        for i in 1:m; b_vec[i]+=basis[i](x)*dx; end
    end
    A_reg=A+lambda_reg*Matrix{Float64}(I,m,m)
    c=A_reg\b_vec; (x->sum(c[i]*basis[i](x) for i in 1:m), c)
end

"Galerkin 法: 加权残差 (用于 PDE → ODE)"
function galerkin_method(el_operator::Function, basis::Vector{Function}, a::Float64, b::Float64; n=500)
    m=length(basis); A=zeros(m,m); b_vec=zeros(m); dx=(b-a)/n
    for k in 1:n; x=a+(k-0.5)*dx
        for i in 1:m, j in 1:m; A[i,j]+=el_operator(basis[j],x)*basis[i](x)*dx; end
        for i in 1:m; b_vec[i]+=0.0*dx; end  # 齐次问题
    end
    c=A\b_vec; (x->sum(c[i]*basis[i](x) for i in 1:m), c)
end

"Kantorovich 法 (半离散化: 仅部分积分)"
function kantorovich_method(F::Function, x_basis::Vector{Function}, y_range::Tuple{Float64,Float64}; nx=100, ny=50)
    # 2D → 1D 的降维方法
    m=length(x_basis); c_funcs=[y->0.0 for _ in 1:m]  # placeholder
    c_funcs
end

# ============================================================
# Morse 理论 (极小极大的计数)
# ============================================================
"共轭点的 Morse 指数: (a,b) 内共轭点的个数"
function morse_index(conjugate_points::Vector{Float64}, a::Float64, b::Float64)
    count(pt->a<pt<b, conjugate_points)
end

"共轭点间的零点计数"
function count_jacobi_zeros(u_sol::Vector{Float64})
    count=0
    for i in 2:length(u_sol); u_sol[i-1]*u_sol[i]<0 && (count+=1); end
    count
end

export legendre_condition, strong_legendre, verify_legendre_along_extremal
export find_conjugate_points, jacobi_condition, first_conjugate_point
export weierstrass_excess, weierstrass_condition, verify_weierstrass_along_extremal
export hilbert_invariant_integral, slope_field_from_extremal, sufficient_conditions_strong
export ritz_method, galerkin_method, kantorovich_method
export morse_index, count_jacobi_zeros
