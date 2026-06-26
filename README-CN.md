# Mini Classical Mechanics（迷你经典力学）

一套**从零开始、零依赖的 C 语言实现**，涵盖大学层次的经典力学。每个子模块对标 MIT（及其他顶尖大学）课程，将 Goldstein、Landau、Marion-Thornton 等教材中的公式转化为可运行的 C 仿真代码。

## 子模块

| 子模块 | 主题 | 参考课程 |
|------|------|---------|
| [mini-newtonian](mini-newtonian/) | 牛顿定律、质点动力学、简谐振子、有心力、非惯性系 | MIT 8.01, Stanford PHYSICS 41 |
| [mini-lagrangian](mini-lagrangian/) | 拉格朗日表述、广义坐标、变分法、欧拉-拉格朗日方程、约束系统 | MIT 8.09, Landau Vol.1 |
| [mini-hamiltonian](mini-hamiltonian/) | 哈密顿力学、相空间、正则变换、泊松括号、可积性、作用-角变量 | MIT 8.09, Goldstein Ch.8-10 |
| [mini-rigid-body](mini-rigid-body/) | 刚体运动学、欧拉角、惯性张量、欧拉方程、无力矩运动、重陀螺 | MIT 2.003, Goldstein Ch.4-5 |
| [mini-continuum](mini-continuum/) | 连续介质力学、应力应变张量、本构关系、弹性力学、连续梁 | MIT 2.071, Landau Vol.7 |
| [mini-celestial-mechanics](mini-celestial-mechanics/) | 开普勒问题、二体轨道动力学、摄动理论、限制性三体问题、拉格朗日点 | MIT 8.288, Murray-Dermott |
| [mini-chaos](mini-chaos/) | 非线性动力学、相图、分岔、李雅普诺夫指数、庞加莱截面、KAM 理论 | MIT 18.S197, Strogatz |
| [mini-variational-principle](mini-variational-principle/) | 最小作用量原理、哈密顿原理、Noether 定理、对称性与守恒律 | MIT 8.09, Landau Vol.1 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个子模块自带 `Makefile`、`include/`、`src/`、`tests/`、`examples/`、`docs/`
- **理论到代码的映射** — 每个模块的 `docs/` 目录包含公式推导和课程对标
- **物理仿真库** — ODE 积分器（RK4、Verlet、辛积分）、相空间可视化、轨道传播器

## 构建方式

每个子模块相互独立。进入子模块目录后运行：

```bash
cd mini-newtonian
make all    # 构建库和示例
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
0. mini-classical-mechanics/
├── mini-newtonian/               # 牛顿定律、质点动力学
├── mini-lagrangian/              # 拉格朗日表述、变分法
├── mini-hamiltonian/             # 哈密顿力学、相空间
├── mini-rigid-body/              # 刚体运动学与动力学
├── mini-continuum/               # 连续介质力学
├── mini-celestial-mechanics/     # 天体力学、轨道动力学
├── mini-chaos/                   # 非线性动力学与混沌
├── mini-variational-principle/   # 最小作用量与 Noether 定理

```

## 许可证

MIT
