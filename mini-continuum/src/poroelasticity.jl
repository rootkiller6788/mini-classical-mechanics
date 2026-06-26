# poroelasticity.jl — 多孔弹性介质: Biot理论, Terzaghi, 扩散-变形耦合
# 参考: Wang "Theory of Linear Poroelasticity", Biot 1941

using LinearAlgebra

"Biot 有效应力: σ_eff = σ + α p I"
function biot_effective_stress(sigma::StressTensor, p::Float64, alpha::Float64)
    StressTensor(sigma.xx+alpha*p, sigma.yy+alpha*p, sigma.zz+alpha*p, sigma.xy, sigma.xz, sigma.yz)
end

"Biot 系数: α = 1 - K_d/K_s"
biot_coefficient(K_drained::Float64, K_solid::Float64) = 1.0 - K_drained/K_solid

"储水系数: S = φ/K_f + (α-φ)/K_s"
function storage_coefficient(phi::Float64, K_fluid::Float64, alpha::Float64, K_solid::Float64)
    phi/K_fluid + (alpha-phi)/K_solid
end

"Terzaghi 1D 固结: ∂p/∂t = c_v ∂²p/∂z²"
function terzaghi_consolidation_1d(cv::Float64, H::Float64, t_end::Float64, p0::Float64; nz=100)
    dz=H/(nz-1); dt_fac=0.4; dt=dt_fac*dz^2/cv; nt=Int(ceil(t_end/dt))
    p=zeros(nz,nt+1); zs=range(0,H,length=nz)
    p[2:end-1,1].=p0; r=cv*dt/dz^2
    for n in 1:nt
        for i in 2:nz-1
            p[i,n+1]=p[i,n]+r*(p[i+1,n]-2*p[i,n]+p[i-1,n])
        end; p[1,n+1]=0.0; p[nz,n+1]=0.0
    end; zs, range(0,dt*nt,length=nt+1), p
end

"固结系数: c_v = k/(μ_f S)"
consolidation_coefficient(k::Float64, mu_f::Float64, S::Float64) = k/(mu_f*S)

"渗透率 (Darcy): q = -k/μ_f ∇p"
darcy_velocity(k::Float64, mu_f::Float64, grad_p::Float64) = -k/mu_f*grad_p

"Biot-Gassmann 饱和岩石体积模量: K_sat = K_d + α²/( (α-φ)/K_s + φ/K_f )"
function biot_gassmann(K_dry::Float64, K_solid::Float64, K_fluid::Float64, phi::Float64)
    alpha = biot_coefficient(K_dry, K_solid)
    K_dry + alpha^2/( (alpha-phi)/K_solid + phi/K_fluid )
end

"Skempton B系数: Δp = -B Δσ/3"
function skempton_B(K_dry::Float64, K_solid::Float64, K_fluid::Float64, phi::Float64)
    alpha = biot_coefficient(K_dry, K_solid)
    1/( 1 + phi*(1/K_fluid - 1/K_solid)/(1/K_dry - 1/K_solid) )
end

"多孔弹性波速 (快P波)"
function biot_fast_p_wave(K_sat::Float64, G::Float64, rho::Float64)
    sqrt((K_sat + 4*G/3)/rho)
end

export biot_effective_stress, biot_coefficient, storage_coefficient
export terzaghi_consolidation_1d, consolidation_coefficient, darcy_velocity
export biot_gassmann, skempton_B, biot_fast_p_wave
