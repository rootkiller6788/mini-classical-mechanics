# euler_lagrange.jl — 欧拉-拉格朗日方程: BVP求解, 多变量, 高阶导数, 自然边界
# 参考: Gelfand & Fomin Ch.1-3, Goldstein Ch.2

using LinearAlgebra

# ============================================================
# 1D Euler-Lagrange BVP 求解
# ============================================================
"打靶法求解 1D Euler-Lagrange 两点边值问题 (y(a)=ya, y(b)=yb)"
function solve_el_bvp(F::Function, dFdy::Function, dFdyp::Function, d2Fdyp2::Function, a::Float64, b::Float64, ya::Float64, yb::Float64; n_grid=500, tol=1e-8)
    function shoot(yp0::Float64)
        dx=(b-a)/(n_grid-1); y=zeros(n_grid); yp=zeros(n_grid)
        y[1]=ya; yp[1]=yp0
        for i in 1:n_grid-1; x=a+(i-1)*dx
            d2=d2Fdyp2(x,y[i],yp[i]); d2=abs(d2)<1e-15 ? 1.0 : d2
            yp[i+1]=yp[i]+dx*(dFdy(x,y[i],yp[i])-dFdyp(x,y[i],yp[i]))/d2
            y[i+1]=y[i]+dx*yp[i]
        end; y[end]-yb
    end
    lo,hi=-20.0,20.0; flo,fhi=shoot(lo),shoot(hi)
    for _ in 1:60
        abs(flo)<tol && return lo; abs(fhi)<tol && return hi
        mid=(lo+hi)/2; fmid=shoot(mid)
        flo*fmid<0 ? (hi=mid;fhi=fmid) : (lo=mid;flo=fmid)
    end; (lo+hi)/2
end

# ============================================================
# 边界条件
# ============================================================
"自然边界条件: ∂F/∂y'' = 0 at x=a,b (当端点自由)"
natural_boundary_condition(F::Function, dFdyp::Function, x::Float64, y_val::Float64, yp_val::Float64) = dFdyp(x,y_val,yp_val)

"横截条件 (端点沿曲线 φ 变动): F+(y''_curve-y'')∂F/∂y'' = 0"
function transversality_condition(F::Function, dFdyp::Function, phi::Function, phip::Function, x::Float64, y_val::Float64, yp_val::Float64)
    F(x,y_val,yp_val)+(phip(x)-yp_val)*dFdyp(x,y_val,yp_val)
end

"循环坐标第一积分: ∂F/∂y'' = const (当 F 不显含 y)"
first_integral_cyclic(dFdyp::Function, x::Float64, y::Float64, yp::Float64) = dFdyp(x,y,yp)

"边界项: [η ∂F/∂y'']_a^b (部分积分留出的边界贡献)"
function boundary_term(dFdyp::Function, eta::Function, a::Float64, b::Float64, y_sol::Function, yp_sol::Function)
    eta(b)*dFdyp(b,y_sol(b),yp_sol(b))-eta(a)*dFdyp(a,y_sol(a),yp_sol(a))
end

# ============================================================
# 多变量 Euler-Lagrange 方程
# ============================================================
"多变量 EL (N个未知函数 y₁,...,y_N): ∂F/∂y_i - d/dx(∂F/∂y_i'') = 0, i=1..N"
function multi_variable_el_residual(F::Function, dFdy::Vector{Function}, dFdyp::Vector{Function}, ys::Vector{Function}, yps::Vector{Function}, x::Float64; dx=1e-4)
    N=length(ys); residuals=zeros(N)
    for i in 1:N
        d_term=(dFdyp[i](x+dx,[y(x+dx) for y in ys]...,[yp(x+dx) for yp in yps]...)-dFdyp[i](x-dx,[y(x-dx) for y in ys]...,[yp(x-dx) for yp in yps]...))/(2dx)
        residuals[i]=dFdy[i](x,[y(x) for y in ys]...,[yp(x) for yp in yps]...)-d_term
    end; residuals
end

# ============================================================
# 高阶导数 Euler-Lagrange 方程
# ============================================================
"高阶导数 EL: F(x,y,y'',y'''',...,y^{(m)})
    ∂F/∂y - d/dx(∂F/∂y'') + d²/dx²(∂F/∂y'''') - ... + (-1)^m d^m/dx^m(∂F/∂y^{(m)}) = 0"
function higher_order_el_operator(F::Function, dF_ders::Vector{Function}, y::Function, y_ders::Vector{Function}, x::Float64; dx=1e-4)
    m=length(dF_ders); result=0.0
    for k in 0:m-1
        term(xv)=dF_ders[k+1](xv,[f(xv) for f in y_ders]...)
        dk_term=numerical_derivative_k(term,x,k,dx)
        result+=(-1)^k*dk_term
    end; result
end

function numerical_derivative_k(f::Function, x::Float64, k::Int, dx::Float64)
    k==0 && return f(x)
    (numerical_derivative_k(f,x+dx,k-1,dx)-numerical_derivative_k(f,x-dx,k-1,dx))/(2dx)
end

# ============================================================
# 参数化 Euler-Lagrange
# ============================================================
"参数化曲线: 最小化 ∫F(t,x,ẋ)dt, 其中 x=(x₁,...,x_n)"
function parametric_el_residual(F::Function, dFdx::Function, dFdxdot::Function, x::Function, xdot::Function, t::Float64; dt=1e-4)
    d_term=(dFdxdot(t+dt,x(t+dt),xdot(t+dt))-dFdxdot(t-dt,x(t-dt),xdot(t-dt)))/(2dt)
    dFdx(t,x(t),xdot(t))-d_term
end

"最短路径 (测地线): F=√(Σ ẋ_i²)"
function geodesic_integrand(xdot::Vector{Float64})
    sqrt(sum(xd->xd^2, xdot))
end

"极小曲面 (面积泛函): F(u,∇u)=√(1+|∇u|²)"
function minimal_surface_integrand(ux::Float64, uy::Float64)
    sqrt(1+ux^2+uy^2)
end

# ============================================================
# 场论 Euler-Lagrange
# ============================================================
"场论 EL: ∂F/∂u - ∂/∂x(∂F/∂u_x) - ∂/∂y(∂F/∂u_y) = 0"
function field_el_residual(F::Function, dFdu::Function, dFdux::Function, dFduy::Function, u::Function, ux::Function, uy::Function, x::Float64, y::Float64; dxy=1e-3)
    term_x=(dFdux(x+dxy,y,ux(x+dxy,y),uy(x+dxy,y))-dFdux(x-dxy,y,ux(x-dxy,y),uy(x-dxy,y)))/(2dxy)
    term_y=(dFduy(x,y+dxy,ux(x,y+dxy),uy(x,y+dxy))-dFduy(x,y-dxy,ux(x,y-dxy),uy(x,y-dxy)))/(2dxy)
    dFdu(x,y,ux(x,y),uy(x,y))-term_x-term_y
end

"波动方程: F=½(u_t²-c²u_x²)"
function wave_equation_residual(c::Float64, u::Function, ut::Function, ux::Function, x::Float64, t::Float64; dt=1e-3, dx=1e-3)
    utt=(ut(x,t+dt)-ut(x,t-dt))/(2dt)
    uxx=(ux(x+dx,t)-ux(x-dx,t))/(2dx)
    utt-c^2*uxx
end

# ============================================================
# 断折极值线 (Weierstrass-Erdmann 角条件)
# ============================================================
"角条件: ∂F/∂y'' 和 F-y''∂F/∂y'' 在角点连续"
function corner_condition(F::Function, dFdyp::Function, x_c::Float64, yl::Float64, ypl::Float64, yr::Float64, ypr::Float64)
    p_cont=abs(dFdyp(x_c,yl,ypl)-dFdyp(x_c,yr,ypr))
    H_cont=abs((F(x_c,yl,ypl)-ypl*dFdyp(x_c,yl,ypl))-(F(x_c,yr,ypr)-ypr*dFdyp(x_c,yr,ypr)))
    p_cont<1e-10 && H_cont<1e-10
end

export solve_el_bvp, natural_boundary_condition, transversality_condition
export first_integral_cyclic, boundary_term
export multi_variable_el_residual, higher_order_el_operator
export parametric_el_residual, geodesic_integrand, minimal_surface_integrand
export field_el_residual, wave_equation_residual, corner_condition
