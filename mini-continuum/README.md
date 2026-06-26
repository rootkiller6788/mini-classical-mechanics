# mini-continuum — 连续介质力学 (Julia + C)

> 参考 Landau Vol.7, Timoshenko, Achenbach, Fung
> 弹性本构、弹性波、梁板壳、强度理论、塑性、断裂、接触、复合材料
> **双语言实现**: Julia (高层API) + C (核心计算库)

## Module Status: COMPLETE ✅

| Level | Status | Score |
|-------|--------|-------|
| L1 Definitions | Complete | 2/2 |
| L2 Core Concepts | Complete | 2/2 |
| L3 Math Structures | Complete | 2/2 |
| L4 Fundamental Laws | Complete | 2/2 |
| L5 Computational Methods | Complete | 2/2 |
| L6 Canonical Systems | Complete | 2/2 |
| L7 Applications | Complete | 2/2 |
| L8 Advanced Topics | Complete | 2/2 |
| L9 Research Frontiers | Partial | 1/2 |
| **Total Score** | | **17/18** |

**include/ + src/** 合计: **3020 行** ≥ 3000 ✓

## 模块对照

| 子模块 | 内容 | 实现文件 |
|--------|------|---------|
| Types | 应力/应变张量、弹性材料、变形梯度、旋转张量 | `src/types.jl` + `include/continuum_types.h` + `src/continuum_types.c` |
| Elasticity | Hooke定律、弹性常数互转、各向异性、非线性弹性 | `src/elasticity.jl` + `include/continuum_elasticity.h` + `src/continuum_elasticity.c` |
| Waves | P/S/Rayleigh/Love/Lamb波、FDTD、Snell定律 | `src/waves.jl` + `include/continuum_waves.h` + `src/continuum_waves.c` |
| Beams | Euler-Bernoulli/Timoshenko梁、屈曲、Winkler地基 | `src/beams.jl` + `include/continuum_beams.h` + `src/continuum_beams.c` |
| Plates | Kirchhoff板、Mindlin板、膜振动、壳 | `src/plates.jl` |
| Failure | von Mises/Tresca/Mohr-Coulomb/Drucker-Prager/疲劳/断裂 | `src/failure.jl` + `include/continuum_failure.h` + `src/continuum_failure.c` |
| Plasticity | J2塑性、Drucker-Prager、Cam-Clay、硬化 | `src/plasticity.jl` |
| Fracture | K因子、J积分、Paris、Weibull、Griffith | `src/fracture.jl` |
| Composites | 混合律、Halpin-Tsai、层合板ABD、Eshelby | `src/composites.jl` |
| Viscoelasticity | Maxwell/Kelvin-Voigt/SLS/Prony/Boltzmann | `src/viscoelasticity.jl` |
| Contact | Hertz/JKR/DMT/Mindlin/Stribeck摩擦 | `src/contact.jl` |
| Thermoelastic | 热传导、热应力、热弹性阻尼、Biot数 | `src/thermoelastic.jl` |
| Poroelasticity | Biot理论、Terzaghi固结、Gassmann | `src/poroelasticity.jl` |
| Conservation | 质量/动量/能量守恒、Clausius-Duhem | `src/conservation.jl` |
| Creep | Norton/Larson-Miller/Kachanov损伤 | `src/creep.jl` |
| Stability | 屈曲/后屈曲/分岔/弧长法 | `src/stability.jl` |
| Micromechanics | Eshelby/Mori-Tanaka/位错/Hall-Petch | `src/micromechanics.jl` |

## API Reference (200+ functions)

### 核心定义 (L1)
`StressTensor` `StrainTensor` `ElasticMaterial` `DeformationGradient` `RotationTensor`
`stress_matrix` `strain_tensor` `hydrostatic_pressure` `deviatoric_stress`
`principal_stresses` `stress_invariants` `dilatation` `equivalent_strain`

### 弹性 (L2/L4)
`hookes_law_isotropic` `hookes_law_inverse` `lame_constants` `engineering_constants`
`elastic_constants_table` `from_bulk_shear` `plane_stress_stiffness` `plane_strain_stiffness`
`orthotropic_compliance` `transversely_isotropic_stiffness` `strain_energy_density`
`complementary_energy` `thermoelastic_stress` `neo_hookean_energy` `mooney_rivlin_energy`

### 弹性波 (L5/L6)
`p_wave_speed` `s_wave_speed` `rayleigh_wave_speed` `rayleigh_wave_speed_exact`
`love_wave_speed_range` `love_wave_dispersion` `wave_1d_solve` `d_alembert_solution`
`reflection_coefficient` `transmission_coefficient` `snell_angle` `critical_angle`
`free_surface_reflection_coefficients` `group_velocity` `timoshenko_beam_dispersion`

### 梁 (L5/L6)
`euler_bernoulli_deflection` `beam_max_deflection_uniform` `beam_max_deflection_point`
`cantilever_deflection` `cantilever_uniform` `timoshenko_deflection_point`
`shear_deformation_ratio` `euler_buckling_load` `BUCKLING_K` `critical_stress`
`radius_of_gyration` `slenderness_ratio` `beam_natural_frequencies`
`cantilever_fundamental_frequency` `lateral_torsional_buckling` `secant_formula`
`winkler_characteristic_length` `winkler_deflection_point` `southwell_fit`

### 板/壳/膜 (L6)
`plate_flexural_rigidity` `circular_plate_deflection` `clamped_circular_plate_max`
`plate_natural_frequencies` `plate_mode_shape` `BOUNDARY_LAMBDA`
`membrane_frequencies` `membrane_wave_speed` `cylindrical_shell_frequency`

### 强度/疲劳/断裂 (L4/L6)
`von_mises` `tresca` `max_principal` `mohr_coulomb_shear` `drucker_prager`
`safety_factor_static` `basquin_sn_curve` `fatigue_life` `strain_life`
`goodman` `gerber` `soderberg` `morrow_equivalent` `miner_damage`
`stress_intensity_factor_mode1` `energy_release_rate` `paris_law` `weibull_failure_probability`

### 塑性/蠕变/粘弹性 (L8)
`von_mises_yield_function` `radial_return` `isotropic_hardening` `cam_clay_yield`
`norton_creep_rate` `larson_miller_parameter` `kachanov_damage_rate`
`maxwell_stress_relaxation` `prony_series` `complex_modulus_maxwell`

### 守恒/微观/稳定性 (L3/L8)
`first_pk_stress` `second_pk_stress` `kirchhoff_stress` `clausius_duhem`
`eshelby_stress` `mori_tanaka_bulk` `hashin_shtrikman_bounds`
`plate_buckling_load` `shell_buckling_classical` `dynamic_buckling_load`

## 核心定理 (L4)
- **Hooke定律**: σ = λ tr(ε)I + 2με
- **Cauchy运动方程**: ρ a = ∇·σ + ρ b
- **Euler屈曲**: P_cr = n²π²EI/(KL)²
- **Griffith断裂**: σ_cr = √(2Eγ/(πa))
- **Clausius-Duhem**: σ:D - ρ(ψ̇ + sṪ) - q·∇T/T ≥ 0
- **Paris律**: da/dN = C(ΔK)^m

## 九校课程映射
MIT 2.071/2.080 · Stanford ME 340 · Berkeley CE 230 · Caltech Ae/AM 102
Princeton MAE 542 · Cambridge Part II/III · Oxford C6 · ETH 151-0515 · Tokyo Univ

## 运行

```bash
# Julia
julia examples/wave_and_beam.jl
julia tests/runtests.jl

# C library
make
make test
```
