# nbody_wrapper.jl — Julia → C++ N体加速度内核调用封装
# 用法:
#   include("cpp/nbody_wrapper.jl")
#   acc = nbody_cpp_accel(masses, positions; G_val=G, softening=0.0)

module NBodyCPP

using ..Newtonian: Vec3

# 尝试加载 C++ 共享库
const lib_paths = [
    joinpath(@__DIR__, "build", "libnbody.so"),
    joinpath(@__DIR__, "build", "libnbody.dylib"),
    joinpath(@__DIR__, "build", "libnbody.dll"),
]

const lib_loaded = Ref{Bool}(false)

function __init__()
    for p in lib_paths
        if isfile(p)
            try
                # Julia 会自动处理 ccall
                lib_loaded[] = true
                return
            catch
            end
        end
    end
    @warn "C++ nbody library not found. Using pure Julia fallback."
end

"""
C++ 加速 N 体加速度计算
如果 C++ 库可用则调用，否则退回纯 Julia 实现
"""
function nbody_accel(masses::Vector{Float64}, positions::Vector{Vec3};
                     G_val::Float64=Newtonian.G, softening::Float64=0.0)
    n = length(masses)
    n == length(positions) || error("masses and positions must have same length")

    acc = [Vec3() for _ in 1:n]

    if lib_loaded[]
        # 准备 C 数组
        m_arr = masses
        x_arr = [p.x for p in positions]
        y_arr = [p.y for p in positions]
        z_arr = [p.z for p in positions]
        ax_arr = zeros(Float64, n)
        ay_arr = zeros(Float64, n)
        az_arr = zeros(Float64, n)

        ccall((:nbody_accel_direct, lib_paths[1]), Cvoid,
              (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
               Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64, Float64),
              n, m_arr, x_arr, y_arr, z_arr, ax_arr, ay_arr, az_arr,
              G_val, softening)

        for i in 1:n
            acc[i] = Vec3(ax_arr[i], ay_arr[i], az_arr[i])
        end
    else
        # 纯 Julia fallback (对称性优化)
        for i in 1:n
            for j in i+1:n
                dr = positions[i] - positions[j]
                dist2 = Newtonian.norm2(dr) + softening^2
                dist = sqrt(dist2)
                dist3 = dist2 * dist
                factor = G_val / dist3

                acc[i] = acc[i] - masses[j] * factor * dr
                acc[j] = acc[j] + masses[i] * factor * dr
            end
        end
    end

    return acc
end

export nbody_accel

end # module NBodyCPP
