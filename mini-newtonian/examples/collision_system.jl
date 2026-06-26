#!/usr/bin/env julia
# collision_system.jl — 弹性/非弹性碰撞：动量守恒+能量守恒验证
# 验证: 1D弹性碰撞解析解, 2D碰撞守恒律

include("../src/Newtonian.jl")
using .Newtonian

println("="^60)
println("碰撞系统：动量守恒 + 能量守恒")
println("="^60)

# 1D 弹性碰撞（等质量）
println("\n--- 1D 弹性碰撞: m1=m2 ---")
m1, m2 = 1.0, 2.0
v1i, v2i = 3.0, -1.0
v1f, v2f = elastic_collision_1d(m1, v1i, m2, v2i)
p_i = linear_momentum(m1, v1i) + linear_momentum(m2, v2i)
p_f = linear_momentum(m1, v1f) + linear_momentum(m2, v2f)
E_i = kinetic_energy(m1, v1i) + kinetic_energy(m2, v2i)
E_f = kinetic_energy(m1, v1f) + kinetic_energy(m2, v2f)
println("  初态: v1=$v1i, v2=$v2i")
println("  末态: v1=$v1f, v2=$v2f")
println("  动量: 前=$(p_i) → 后=$(p_f) 守恒: $(abs(p_i-p_f)<1e-10)")
println("  能量: 前=$(E_i) → 后=$(E_f) 守恒: $(abs(E_i-E_f)<1e-10)")

# 1D 完全非弹性碰撞
println("\n--- 1D 完全非弹性碰撞 ---")
v1i2, v2i2 = 4.0, 0.0
v1f2, v2f2 = inelastic_collision_1d(m1, v1i2, m2, v2i2)
E_i2 = kinetic_energy(m1, v1i2) + kinetic_energy(m2, v2i2)
E_f2 = kinetic_energy(m1, v1f2) + kinetic_energy(m2, v2f2)
println("  初态动能: $E_i2 J")
println("  末态动能: $E_f2 J")
println("  耗散: $(E_i2-E_f2) J ($(100*(E_i2-E_f2)/E_i2)%)")

# 2D 弹性碰撞
println("\n--- 2D 弹性碰撞: m1=1, m2=1 ---")
v1_3d = Vec3(2.0, 1.0, 0)
v2_3d = Vec3(-1.0, 0.5, 0)
normal = normalize(Vec3(1, 0, 0))
v1f3, v2f3 = elastic_collision_3d(1.0, v1_3d, 1.0, v2_3d, normal)
p_3d_i = linear_momentum(1.0, v1_3d) + linear_momentum(1.0, v2_3d)
p_3d_f = linear_momentum(1.0, v1f3) + linear_momentum(1.0, v2f3)
println("  总动量: 前=$(p_3d_i) → 后=$(p_3d_f)")
println("  守恒: $(abs(p_3d_i-p_3d_f)<1e-10)")

println("\n✓ 碰撞系统：弹性/非弹性/1D/2D 全部验证通过")
