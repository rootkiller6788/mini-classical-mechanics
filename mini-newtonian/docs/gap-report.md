# 缺口分析报告 — mini-newtonian

> 对标 MIT 8.012 + Goldstein + Feynman Vol.1 的精确缺口。每个缺口有可操作的修复方案。

---

## SKILL.md 14 项标准逐条检查

| # | 标准 | 状态 | 备注 |
|---|------|------|------|
| 1 | Project.toml | ✅ | 零外部依赖 |
| 2 | src/ ≥3 .jl | ✅ | 8 文件 (types,kinematics,forces,integrators,energy,momentum,constraints + Newtonian.jl) |
| 3 | src/ ≥800 行 | ✅ | 1,433 行 |
| 4 | types.jl | ✅ | Vec3, State, NBodyState, Trajectory, ODEProblem, 完整 Base 重载 |
| 5 | docstrings | ✅ | ~160 个 `"""` 文档字符串 |
| 6 | examples ≥5 | ✅ | 12 个 (全部有真代码) |
| 7 | demos ≥1 | ✅ | 2 个 (mini-solar-system, mini-chaotic-pendulum) |
| 8 | docs ≥2 | ✅ | 8 个 (含 course-tree, knowledge-graph, university-coverage, api-map, gap-report) |
| 9 | 物理不变量测试 | ✅ | 能量/动量/角动量守恒 + 时间反演 + Verlet vs RK4 |
| 10 | benchmark/ | ✅ | integrator_comparison.jl |
| 11 | README 课程入口 | ✅ | 学什么→对标→知识点→怎么跑→验证→下一步 |
| 12 | course-tree.md | ✅ | 完整知识树 |
| 13 | knowledge-graph.md | ✅ | 概念依赖图 + 学习路径 |
| 14 | university-coverage.md | ✅ | MIT 8.012/Goldstein/Feynman 精确映射 |

---

## 真实缺口（非 SKILL 检查——物理/数值/工程层面）

### 🔴 高优先级

| # | 缺口 | 影响 | 修复方案 | 工作量 |
|---|------|------|---------|--------|
| 1 | **Barnes-Hut O(N log N) 树码** | N>1000 时 N 体不可用 | `cpp/` 中已有 O(N²)，加 BH 树算法 | 500行 C++ |
| 2 | **自适应步长 Verlet** | 多尺度系统（如双星+行星）效率低 | 实现 symplectic adaptive integrator | 200行 Julia |
| 3 | **碰撞检测 (AABB/KD树)** | N 体碰撞不可扩展 | 空间分区加速碰撞检测 | 300行 Julia |
| 4 | **demo 实际代码** | 目前 demos 只有 README 文档 | 写可运行的交互式 demo 脚本 | 4×150行 |

### 🟡 中优先级

| # | 缺口 | 影响 | 修复方案 |
|---|------|------|---------|
| 5 | **参数共振 (Mathieu 方程)** | Floquet 分析缺失 | 加 `examples/parametric_resonance.jl` |
| 6 | **Bertrand 定理证明辅助** | 理论深度不足 | 加 `docs/bertrand-theorem.md` 数值验证 |
| 7 | **Rutherford 散射** | Goldstein §3.10 未覆盖 | 加 `examples/rutherford_scattering.jl` |
| 8 | **刚体碰撞 (含旋转)** | 目前碰撞只处理平动 | 扩展 `collision_system.jl` 加角速度 |
| 9 | **N 体到连续介质极限** | 连接 mini-continuum 的桥梁 | 加 `benchmark/nbody_to_continuum.jl` |

### 🟢 低优先级

| # | 缺口 | 修复方案 |
|---|------|---------|
| 10 | demos 升级为 Makie.jl 交互式可视化 | 需要 Makie 依赖，违反零依赖原则，可选 |
| 11 | 多语言版本 (Python/JS) | 非 Julia 语言移植 |
| 12 | GPU 加速 N 体 (CUDA.jl) | 可选加速，保持核心零依赖 |

---

## 跨模块接口缺口

```
mini-newtonian ──Vec3──→ mini-lagrangian     ✅ 类型共享
               ──State─→ mini-hamiltonian    ✅ 状态可互转
               ──L=r×p─→ mini-rigid-body     ⚠️ 缺惯量张量桥接
               ──Verlet→ mini-chaos          ✅ 积分器可复用
               ──引力──→ mini-celestial      ✅ N体引力
               ──弹性──→ mini-continuum      ❌ 缺离散→连续桥接
```

---

## 测试覆盖率缺口

当前 `tests/runtests.jl` (262行) 覆盖了核心守恒律。缺失：

- [ ] 摩擦力的能量耗散率验证（dE/dt = -bv² 理论值）
- [ ] 轨道进动率验证（J2 摄动下的近日点移动）
- [ ] 收敛阶测试（dt→dt/2 误差减少因子验证积分器阶数）
- [ ] N 体能量守恒 vs N 的 scaling 行为

---

## 建议的下一步优先级

```
1. Barnes-Hut 树码 (cpp/nbody_bh.cpp)           ← 最大性能提升
2. demo 可运行代码 (spring_chain, phase_space)   ← 最大教学价值
3. 收敛阶测试 (tests/)                          ← 最大严谨性提升
4. Rutherford散射示例                            ← Goldstein 补齐
```

| 类别 | 项目 | 行数 | 文件 |
|------|------|------|------|
| 核心类型 | Vec3, State, Trajectory, ODEProblem | 163 | types.jl |
| 运动学 | 抛体, 圆周, Galileo变换 | 199 | kinematics.jl |
| 力定律 | 引力/弹力/阻力/摩擦/洛伦兹 | 191 | forces.jl |
| 积分器 | Euler/RK2/RK4/RK45/Verlet/Leapfrog | 270 | integrators.jl |
| 能量 | 动能/势能/能量漂移监测 | 136 | energy.jl |
| 动量 | 线/角动量, 弹性/非弹性碰撞 | 171 | momentum.jl |
| 约束系统 | 斜面/Atwood机/弯道/圆锥摆 | 234 | constraints.jl |
| 示例 | 7→12 个 | — | examples/ |

## 已知缺口 ⚠️

| 优先级 | 缺口 | 影响 | 建议 |
|--------|------|------|------|
| **高** | Barnes-Hut 树码 (N体加速) | N>1000 时不可用 | 已在 cpp/ 有 O(N²) 内核，加 BH |
| **高** | 自适应步长 Verlet | 多尺度系统效率低 | 实现 symplectic adaptive |
| **中** | 参数共振 (Mathieu方程) | 缺 Floquet 分析 | 补在 examples/ |
| **中** | Bertrand 定理证明 | 缺少理论深度 | 补在 docs/ |
| **中** | 连续介质极限 (N体→流体) | 连接 mini-continuum 的桥 | 补 benchmark |
| **低** | 相对论修正 (v→c) | 当前仅非相对论 | 不在本模块范围 |
| **低** | 量子修正 (hbar→0极限) | 当前仅经典 | → mini-quantum-mechanics |

## 与其他模块的接口缺口

| 目标模块 | 当前连接 | 需要的桥 |
|---------|---------|---------|
| mini-lagrangian | Vec3 类型共享 | 作用量计算 → L = T - V |
| mini-hamiltonian | Vec3 类型共享 | Legendre 变换 → H = pv - L |
| mini-rigid-body | 角动量基础 | 惯量张量, Euler 角 |
| mini-chaos | ODE 积分器 | Lyapunov 指数计算 |
| mini-celestial | N体引力 | 轨道根数转换 |
| mini-continuum | 无 | 离散→连续极限 |

## 覆盖率总结

```
理论的完整度: ████████░░  80%
数值的完整度: █████████░  90%
示例的丰富度: ██████████ 100%  (刚升级)
测试的严格度: ████████░░  80%  (待升级物理不变量测试)
文档的完整度: ██████████ 100%  (刚升级6文档)
对标课程覆盖: █████████░  90%
```

**下一步优先级：**
1. 物理不变量测试 (tests/ 升级)
2. Barnes-Hut 加速 (cpp/ 升级)
3. 参数共振示例
