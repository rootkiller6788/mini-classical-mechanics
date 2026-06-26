#!/usr/bin/env julia
# tests/runtests.jl — 物理不变量测试（能量/动量/角动量守恒、辛性质、时间反演）
include("../src/Newtonian.jl")
using .Newtonian

const EPS = 1e-10
passed = 0; failed = 0
macro test(name, expr)
    quote
        try
            if $(esc(expr))
                global passed += 1; println("  ✅ ", $(esc(name)))
            else
                global failed += 1; println("  ❌ ", $(esc(name)), " — false")
            end
        catch e
            global failed += 1; println("  ❌ ", $(esc(name)), " — ", e)
        end
    end
end

println("="^60)
println("mini-newtonian 物理不变量测试")
println("="^60)

# === 向量运算 ===
println("\n—— 向量运算 ——")
@test "Vec3构造" (Vec3(1,2,3).x==1 && Vec3(1,2,3).y==2)
@test "向量加法" Vec3(1,2,3)+Vec3(4,5,6)==Vec3(5,7,9)
@test "内积" dot(Vec3(1,0,0),Vec3(0,1,0))≈0.0
@test "叉积" cross(Vec3(1,0,0),Vec3(0,1,0))==Vec3(0,0,1)
@test "模长" norm(Vec3(3,4,0))≈5.0

# === 运动学 ===
println("\n—— 运动学解析解 ——")
@test "抛体射程" abs(projectile_range(10.0,pi/4;g=10.0)-10.0)<EPS

# === 能量守恒 (Verlet 100周期) ===
println("\n—— 能量守恒 (Velocity Verlet, 100周期) ——")
function sho_acc(x,v,t); Vec3(-x.x,0,0); end
ode = ParticleODE(sho_acc,0.0,200*pi,Vec3(1,0,0),Vec3(0,0,0))
traj = solve_fixed_step(ode,2*pi/100;method=:verlet)
E0 = 0.5*1.0*1.0^2
max_dE = maximum(abs(0.5*s.r.x^2+0.5*s.v.x^2-E0)/E0 for s in traj.states)
@test "Verlet 100周期能量漂移 < 1e-4" (max_dE < 1e-4)

# === 动量守恒 ===
println("\n—— 动量守恒（弹性碰撞）——")
m1,m2=1.0,3.0; v1i,v2i=5.0,-2.0
v1f,v2f=elastic_collision_1d(m1,v1i,m2,v2i)
p_i=m1*v1i+m2*v2i; p_f=m1*v1f+m2*v2f
E_i=0.5*m1*v1i^2+0.5*m2*v2i^2; E_f=0.5*m1*v1f^2+0.5*m2*v2f^2
@test "动量守恒" abs(p_i-p_f)<EPS
@test "动能守恒" abs(E_i-E_f)<EPS

# === 角动量守恒 (Kepler) ===
println("\n—— 角动量守恒（中心力 Kepler轨道）——")
function central_acc(x,v,t); r=norm(x); r<1e-10 && return Vec3(0,0,0); -x/r^3; end
ode_k = ParticleODE(central_acc,0.0,20.0,Vec3(1,0,0),Vec3(0,1.0,0))
traj_k = solve_fixed_step(ode_k,0.01;method=:verlet)
L0 = angular_momentum(Vec3(1,0,0),Vec3(0,1,0),1.0)
max_dL = maximum(abs(angular_momentum(s.r,s.v,1.0)-L0) for s in traj_k.states)
@test "角动量守恒 (Kepler)" (max_dL < 1e-8)

# === 时间反演误差 ===
println("\n—— 时间反演误差 (Verlet 辛性质) ——")
ode_fwd = ParticleODE(sho_acc,0.0,10*pi,Vec3(1,0,0),Vec3(0,0,0))
traj_fwd = solve_fixed_step(ode_fwd,2*pi/100;method=:verlet)
ode_rev = ParticleODE(sho_acc,0.0,10*pi,traj_fwd.states[end].r,-traj_fwd.states[end].v)
traj_rev = solve_fixed_step(ode_rev,2*pi/100;method=:verlet)
tr_err = norm(traj_rev.states[end].r - Vec3(1,0,0))
@test "时间反演误差 < 1e-6" (tr_err < 1e-6)

# === Verlet vs RK4 长期 ===
println("\n—— Verlet vs RK4 长期能量漂移 ——")
odes = ParticleODE(sho_acc,0.0,200*pi,Vec3(1,0,0),Vec3(0,0,0))
tv = solve_fixed_step(odes,2*pi/50;method=:verlet)
tr = solve_fixed_step(odes,2*pi/50;method=:rk4)
Ev = maximum(abs(0.5*s.r.x^2+0.5*s.v.x^2-E0)/E0 for s in tv.states)
Er = maximum(abs(0.5*s.r.x^2+0.5*s.v.x^2-E0)/E0 for s in tr.states)
@test "Verlet 优于 RK4 " (Ev < Er)

# === 约束系统 ===
println("\n—— 约束系统 ——")
@test "Atwood机" abs(atwood_machine(1.0,2.0)-9.81*(2-1)/3)<1e-10
@test "环形轨道" abs(loop_the_loop_min_speed(1.0)-sqrt(9.81))<1e-10

println("\n"*"="^60)
println("结果: $passed 通过 / $(passed+failed) 总计")
println(failed==0 ? "✅ 全部物理不变量测试通过" : "❌ $failed 项失败")
println("="^60)
using .Newtonian

const EPS = 1e-10
const RELAXED_EPS = 1e-6

passed = 0
failed = 0

macro test(name, expr)
    quote
        try
            if $(esc(expr))
                global passed += 1
                println("  ✅ ", $(esc(name)))
            else
                global failed += 1
                println("  ❌ ", $(esc(name)), " — returned false")
            end
        catch e
            global failed += 1
            println("  ❌ ", $(esc(name)), " — ", e)
        end
    end
end

function main()
    println("="^60)
    println("  mini-newtonian Test Suite")
    println("="^60)

    # ---- Vec3 Operations ----
    println("\n--- Vec3 ---")
    a = Vec3(1, 2, 3)
    b = Vec3(4, 5, 6)

    @test("addition", norm((a + b) - Vec3(5, 7, 9)) < EPS)
    @test("subtraction", norm((b - a) - Vec3(3, 3, 3)) < EPS)
    @test("scaling", norm((2.0 * a) - Vec3(2, 4, 6)) < EPS)
    @test("dot product", abs(dot(a, b) - 32.0) < EPS)
    @test("cross product", norm(cross(a, b) - Vec3(-3, 6, -3)) < EPS)
    @test("norm", abs(norm(Vec3(3, 4, 0)) - 5.0) < EPS)
    @test("normalize", abs(norm(normalize(Vec3(2, 0, 0))) - 1.0) < EPS)
    @test("distance", abs(distance(Vec3(0,0,0), Vec3(1,2,2)) - 3.0) < EPS)
    @test("norm2", abs(norm2(Vec3(1,2,3)) - 14.0) < EPS)
    @test("angle", abs(angle(Vec3(1,0,0), Vec3(0,1,0)) - π/2) < EPS)

    # ---- Kinematics ----
    println("\n--- Kinematics ---")
    r0 = Vec3(0, 0, 10)
    v0 = Vec3(10, 0, 5)
    r1, v1 = projectile_motion(r0, v0, 2.0; g=10.0)
    @test("projectile z after 2s: z=10+5*2-5*4=0", abs(r1.z) < EPS)
    @test("projectile x after 2s: x=20", abs(r1.x - 20.0) < EPS)
    @test("projectile vz after 2s: 5-20=-15", abs(v1.z + 15.0) < EPS)

    t_flight = projectile_flight_time(5.0, 10.0; g=10.0)
    @test("flight time > 0", t_flight > 0)

    # 45 degree launch range: v0²/g
    rng = projectile_range(10.0, π/4; g=10.0)
    @test("range for v0=10, 45deg, g=10: ~10m", abs(rng - 10.0) < EPS)

    # ---- Forces ----
    println("\n--- Forces ---")
    F_grav = newton_gravity(Vec3(1,0,0), Vec3(), 1.0, 1.0; G_val=1.0)
    @test("gravity direction toward source", F_grav.x < 0)
    @test("gravity magnitude = GMm/r²", abs(norm(F_grav) - 1.0) < EPS)

    F_spring = hooke_force(Vec3(2,0,0), Vec3(), 10.0)
    @test("spring force direction", F_spring.x < 0)
    @test("spring force magnitude", abs(norm(F_spring) - 20.0) < EPS)

    F_drag = linear_drag(Vec3(10,0,0), 0.5)
    @test("drag opposes velocity", F_drag.x < 0 && abs(norm(F_drag) - 5.0) < EPS)

    # ---- ODE Integrators ----
    println("\n--- Integrators ---")

    # 简谐振动精度测试: x'' = -ω²x, ω=2π, T=1
    omega = 2π
    x0_val = 1.0
    function sho_acc(r, v, t) ; return Vec3(-omega^2 * r.x, 0, 0) ; end

    # Velocity Verlet: 一个周期后应回到 x=1
    traj_verlet = solve_fixed_step(sho_acc, Vec3(x0_val,0,0), Vec3(), 1.0, 0.001; method=:verlet)
    @test("Verlet: return to x=1 after 1 period", abs(traj_verlet.positions[end].x - 1.0) < 1e-4)

    # RK4: 一个周期后应回到 x=1（更高精度）
    traj_rk4 = solve_fixed_step(sho_acc, Vec3(x0_val,0,0), Vec3(), 1.0, 0.001; method=:rk4)
    @test("RK4: return to x=1 after 1 period", abs(traj_rk4.positions[end].x - 1.0) < 1e-6)

    # Euler-Cromer: 一个周期后有能量漂移但不大
    traj_ec = solve_fixed_step(sho_acc, Vec3(x0_val,0,0), Vec3(), 1.0, 0.0001; method=:euler_cromer)
    @test("Euler-Cromer: return to ~x=1 after 1 period (relaxed)", abs(traj_ec.positions[end].x - 1.0) < 0.01)

    # RK4 step (1D)
    f_exp(t, y) = [y[1]]  # dy/dt = y → y(t) = y0*exp(t)
    y1 = rk4_step(f_exp, 0.0, [1.0], 0.1)
    expected = exp(0.1)
    @test("RK4 1-step: exp(0.1)", abs(y1[1] - expected) < 1e-8)

    # Adaptive RK45
    traj_ad, accepted, rejected = solve_adaptive(sho_acc, Vec3(x0_val,0,0), Vec3(), 1.0;
                                                  dt_init=0.01, tol=1e-8)
    @test("Adaptive RK45: return to ~x=1", abs(traj_ad.positions[end].x - 1.0) < 1e-5)
    @test("Adaptive RK45: some steps accepted", accepted > 0)

    # ---- Energy Conservation ----
    println("\n--- Energy ---")
    function sho_U(r) ; return 0.5 * omega^2 * r.x^2 ; end
    energies_v, drift_v = energy_drift(traj_verlet, 1.0, sho_U)
    @test("Verlet energy conservation (SHO)", maximum(abs.(drift_v)) < 1e-4)

    energies_rk4, drift_rk4 = energy_drift(traj_rk4, 1.0, sho_U)
    @test("RK4 energy conservation (SHO)", maximum(abs.(drift_rk4)) < 1e-6)

    # ---- Collisions ----
    println("\n--- Collisions ---")

    # 等质量弹性碰撞：交换速度
    v1f, v2f = elastic_collision_1d(1.0, 10.0, 1.0, 0.0)
    @test("elastic eq mass: v1 swaps with v2", abs(v1f) < EPS && abs(v2f - 10.0) < EPS)

    # 大质量撞小质量：小质量反弹
    v1f2, v2f2 = elastic_collision_1d(100.0, 1.0, 1.0, 0.0)
    @test("elastic heavy→light: light gains speed", v2f2 > 1.0)

    # 完全非弹性碰撞：合并
    v_common = inelastic_collision_1d(1.0, 10.0, 1.0, 0.0)
    @test("inelastic eq mass: v = 5", abs(v_common - 5.0) < EPS)

    # 动量守恒
    p_before = 1.0*10.0 + 1.0*0.0
    p_after = 1.0*v1f + 1.0*v2f
    @test("momentum conservation (elastic)", abs(p_before - p_after) < EPS)

    # 3D collision
    v1_3d, v2_3d = elastic_collision_3d(1.0, Vec3(1,0,0), 1.0, Vec3(0,0,0), Vec3(1,0,0))
    @test("3D head-on elastic: v1 stops", norm(v1_3d) < EPS && norm(v2_3d - Vec3(1,0,0)) < EPS)

    # ---- Momentum & Center of Mass ----
    println("\n--- Momentum ---")
    p = linear_momentum(2.0, Vec3(3, 4, 0))
    @test("linear momentum magnitude", abs(norm(p) - 10.0) < EPS)

    L = angular_momentum(Vec3(1,0,0), Vec3(0,1,0), 1.0)
    @test("angular momentum = r×p", abs(L.z - 1.0) < EPS)

    cm = center_of_mass([1.0, 1.0, 1.0], [Vec3(0,0,0), Vec3(3,0,0), Vec3(0,3,0)])
    @test("center of mass (equal masses)", abs(norm(cm - Vec3(1,1,0))) < EPS)

    # Tsiolkovsky
    dv = tsiolkovsky_delta_v(3000.0, 1000.0, 100.0)
    @test("rocket equation: dv positive", dv > 0)
    @test("rocket equation: dv ~ ve*ln(10)", abs(dv - 3000*log(10)) < EPS)

    # ---- Summary ----
    println("\n" * "="^60)
    total = passed + failed
    println("  Results: $passed / $total passed")
    if failed == 0
        println("  🎉 ALL TESTS PASSED")
    else
        println("  ⚠️  $failed test(s) FAILED")
    end
    println("="^60)

    return failed == 0
end

exit(main() ? 0 : 1)
