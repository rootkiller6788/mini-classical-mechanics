# maps.jl — 离散时间映射: Logistic, Hénon, Standard, Circle, Ikeda, Tent, Sinai
# 参考: Strogatz Ch.10, May 1976, Chirikov 1979

using LinearAlgebra

"Logistic map: x_{n+1} = r x_n (1-x_n)"
logistic_map(x::Float64, r::Float64) = r*x*(1-x)

"Cubic map: x_{n+1} = r x_n - x_n^3"
cubic_map(x::Float64, r::Float64) = r*x - x^3

"Sine map: x_{n+1} = r sin(π x_n)"
sine_map(x::Float64, r::Float64) = r*sin(pi*x)

"Hénon map: x_{n+1}=1-a x_n^2+y_n, y_{n+1}=b x_n"
henon_map(x::Float64, y::Float64, a=1.4, b=0.3) = 1-a*x^2+y, b*x

"Lozi map: x_{n+1}=1-a|x_n|+y_n, y_{n+1}=b x_n"
lozi_map(x::Float64, y::Float64, a=1.7, b=0.5) = 1-a*abs(x)+y, b*x

"Standard map (Chirikov-Taylor): p_{n+1}=p_n+K sin θ_n, θ_{n+1}=θ_n+p_{n+1} (mod 2π)"
function standard_map(theta::Float64, p::Float64, K=0.9716)
    p_new=mod(p+K*sin(theta)+pi,2pi)-pi; mod(theta+p_new,2pi), p_new
end

"Circle map: θ_{n+1}=θ_n+Ω-(K/2π)sin(2πθ_n) (mod 1)"
circle_map(theta::Float64, Omega=0.606661, K=0.8) = mod(theta+Omega-K/(2pi)*sin(2pi*theta), 1.0)

"Ikeda map: z_{n+1}=A+B z_n exp(i(κ-p/(1+|z|^2)))"
function ikeda_map(x::Float64, y::Float64, A=1.0, B=0.9, kappa=0.4, p=6.0)
    r2=x^2+y^2; phi=kappa-p/(1+r2); A+B*(x*cos(phi)-y*sin(phi)), B*(x*sin(phi)+y*cos(phi))
end

"Tent map: x_{n+1}=μ min(x_n, 1-x_n)"
tent_map(x::Float64, mu=2.0) = mu*min(x, 1-x)

"Gingerbreadman map: x_{n+1}=1-y_n+|x_n|, y_{n+1}=x_n"
gingerbreadman_map(x::Float64, y::Float64) = 1-y+abs(x), x

"Sinai map: x_{n+1}=x_n+y_n+δ cos(2π y_n) (mod1), y_{n+1}=x_n+2y_n (mod1)"
function sinai_map(x::Float64, y::Float64, delta=0.1)
    mod(x+y+delta*cos(2pi*y), 1.0), mod(x+2y, 1.0)
end

"Arnold cat map: x_{n+1}=2x_n+y_n (mod1), y_{n+1}=x_n+y_n (mod1)"
arnold_cat_map(x::Float64, y::Float64) = mod(2x+y,1.0), mod(x+y,1.0)

"Baker's map: x_{n+1}=2x_n (mod1), y_{n+1}=y_n/2 + (1/2) floor(2x_n)"
function bakers_map(x::Float64, y::Float64)
    sx = 2x > 1 ? 1.0 : 0.0
    mod(2x, 1.0), 0.5*y + 0.5*sx
end

"Gauss map: x_{n+1}=1/x_n mod 1 (x≠0)"
gauss_map(x::Float64) = 1/x - floor(1/x)

"Lorenz map (Poincaré 截面极大值映射)"
lorenz_map(x::Float64, a=2.05, b=-1.0) = a + b*abs(x)^0.5

"迭代映射 (含瞬态消除)"
function iterate_map(f::Function, x0, n_iter::Int, n_transient=500)
    x=x0
    for _ in 1:n_transient; x=f(x); end
    result=[x]
    for _ in 1:n_iter; x=f(x); push!(result,x); end
    result
end

export logistic_map, cubic_map, sine_map
export henon_map, lozi_map, standard_map, circle_map, ikeda_map
export tent_map, gingerbreadman_map, sinai_map, arnold_cat_map, bakers_map
export gauss_map, lorenz_map, iterate_map
