# waves.jl — 弹性波: P波/S波/Rayleigh/Love, 1D/2D波动方程, Snell定律
# 参考: Achenbach "Wave Propagation in Elastic Solids", Landau Vol.7

using LinearAlgebra

# ============================================================
# 体波波速
# ============================================================
p_wave_speed(lam::Float64, mu::Float64, rho::Float64) = sqrt((lam+2mu)/rho)
s_wave_speed(mu::Float64, rho::Float64) = sqrt(mu/rho)
p_wave_speed_E(E::Float64, nu::Float64, rho::Float64) = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2nu)))
s_wave_speed_E(E::Float64, nu::Float64, rho::Float64) = sqrt(E/(2*rho*(1+nu)))

# ============================================================
# 表面波
# ============================================================
"Rayleigh 表面波速 (Bergmann近似)"
function rayleigh_wave_speed(lam::Float64, mu::Float64, rho::Float64)
    cs = s_wave_speed(mu, rho)
    nu = lam/(2*(lam+mu))
    (0.862+1.14*nu)/(1+nu) * cs
end

"Rayleigh 波速 — 精确解 (Viktorov 1967 拟合)"
function rayleigh_wave_speed_exact(lam::Float64, mu::Float64, rho::Float64)
    nu = lam/(2*(lam+mu)); cs = s_wave_speed(mu, rho)
    (0.874 + 0.196*nu - 0.043*nu^2 - 0.055*nu^3) * cs
end

"Love 波速范围 (上层/下层)"
love_wave_speed_range(mu1::Float64, rho1::Float64, mu2::Float64, rho2::Float64) = s_wave_speed(mu1,rho1), s_wave_speed(mu2,rho2)

"Love 波色散关系: tan(kH√(c²/c₁²-1)) = (μ₂/μ₁)√(1-c²/c₂²)/√(c²/c₁²-1)"
function love_wave_dispersion(k::Float64, H::Float64, mu1::Float64, mu2::Float64, rho1::Float64, rho2::Float64, c_guess::Float64)
    c1 = s_wave_speed(mu1, rho1); c2 = s_wave_speed(mu2, rho2)
    !(c1 < c_guess < c2) && return NaN  # c 必须在c1和c2之间
    q1 = k*sqrt(abs(c_guess^2/c1^2 - 1))
    q2 = k*sqrt(abs(1 - c_guess^2/c2^2))
    return tan(q1*H) - (mu2/mu1)*q2/q1
end

# ============================================================
# 1D 波动方程 FDTD
# ============================================================
"1D 波动方程: ∂²u/∂t² = c² ∂²u/∂x² (FDTD)"
function wave_1d_solve(c::Float64, L::Float64, T::Float64, u0::Function, v0::Function; nx=200, CFL=0.5)
    dx=L/(nx-1); dt=CFL*dx/c; nt=Int(ceil(T/dt)); u=zeros(nx,nt+1)
    xs=range(0,L,length=nx)
    for i in 1:nx; u[i,1]=u0(xs[i]); end
    r=(c*dt/dx)^2
    for i in 2:nx-1; u[i,2]=u[i,1]+dt*v0(xs[i])+0.5*r*(u[i+1,1]-2u[i,1]+u[i-1,1]); end
    for n in 2:nt; for i in 2:nx-1; u[i,n+1]=2u[i,n]-u[i,n-1]+r*(u[i+1,n]-2u[i,n]+u[i-1,n]); end; end
    xs, range(0,dt*nt,length=nt+1), u
end

"1D 波动方程 — 达朗贝尔解: u(x,t)=F(x-ct)+G(x+ct)"
function d_alembert_solution(F::Function, G::Function, x::Float64, t::Float64, c::Float64)
    F(x-c*t) + G(x+c*t)
end

# ============================================================
# 反射/透射
# ============================================================
"声阻抗: Z = ρ c"
acoustic_impedance(rho::Float64, c::Float64) = rho*c

"反射系数 (垂直入射, 介质1→介质2)"
reflection_coefficient(Z1::Float64, Z2::Float64) = (Z2-Z1)/(Z2+Z1)

"透射系数 (垂直入射)"
transmission_coefficient(Z1::Float64, Z2::Float64) = 2*Z2/(Z2+Z1)

"Snell 定律 (折射): sinθ₂ = (c₂/c₁)sinθ₁"
function snell_angle(theta1::Float64, c1::Float64, c2::Float64)
    sin_th2 = c2/c1*sin(theta1); abs(sin_th2)>1.0 && return NaN; asin(sin_th2)
end

"全反射临界角: θ_c = arcsin(c₁/c₂)"
critical_angle(c1::Float64, c2::Float64) = c2>c1 ? asin(c1/c2) : NaN

"P波→S波转换 (自由表面反射)"
function free_surface_reflection_coefficients(cp::Float64, cs::Float64, theta_p::Float64)
    sp,cp_val=sin(theta_p),cos(theta_p)
    sin_ts=cs/cp*sp; sin_ts>1 && (sin_ts=1.0)
    cos_ts=sqrt(1-sin_ts^2)
    A=(cs^2*sin(2*theta_p)*sin_ts - cp^2*cos_ts^2)/(cs^2*sin(2*theta_p)*sin_ts + cp^2*cos_ts^2)
    B=2*cp*cs*sin(2*theta_p)*cos_ts/(cs^2*sin(2*theta_p)*sin_ts + cp^2*cos_ts^2)
    return (Rpp=A, Rps=B)
end

# ============================================================
# 色散关系
# ============================================================
"均匀介质色散: ω = c k, v_phase = c, v_group = c"
dispersion_uniform(k::Float64, c::Float64) = c*k, c

"群速度: v_group = dω/dk (中心差分)"
function group_velocity(omega::Function, k::Float64, dk=1e-4)
    (omega(k+dk) - omega(k-dk))/(2*dk)
end

"Timoshenko 梁色散: ω(k)"
function timoshenko_beam_dispersion(k::Float64, E::Float64, G::Float64, rho::Float64, nu::Float64, h::Float64)
    I=h^3/12; A=h; kappa=5/6
    a0=rho*A; a1=-rho*I*(1+E/(kappa*G)); a2=E*I
    c0=a2*k^4; c1=a1*k^2; c2=a0
    sqrt(max((-c1-sqrt(c1^2-4c2*c0))/(2c2),0.0))
end

export p_wave_speed, s_wave_speed, p_wave_speed_E, s_wave_speed_E
export rayleigh_wave_speed, rayleigh_wave_speed_exact
export love_wave_speed_range, love_wave_dispersion
export wave_1d_solve, d_alembert_solution
export acoustic_impedance, reflection_coefficient, transmission_coefficient
export snell_angle, critical_angle, free_surface_reflection_coefficients
export dispersion_uniform, group_velocity, timoshenko_beam_dispersion
