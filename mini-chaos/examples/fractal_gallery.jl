#!/usr/bin/env julia
include("../src/Chaos.jl"); using .Chaos

function main()
    println("="^60)
    println("  Fractal Gallery")
    println("="^60)

    println("\n--- Mandelbrot Set ---")
    for (cx, cy) in [(-0.5, 0.0), (0.25, 0.5), (-1.0, 0.0), (-0.75, 0.1)]
        iter = mandelbrot_iter(cx, cy; max_iter=50)
        println("  M($cx, $cy) → $iter iterations")
    end

    println("\n--- Julia Sets ---")
    for (cx, cy) in [(-0.7, 0.27), (0.355, 0.355), (-0.4, 0.6)]
        println("  J(c=$cx+$(cy)i):")
        for (zx, zy) in [(0.5, 0.0), (-0.5, 0.5), (0.0, -0.5)]
            iter = julia_iter(zx, zy, cx, cy; max_iter=50)
            print("    z₀=$zx+$(zy)i → $iter  ")
        end
        println()
    end

    println("\n--- Newton Fractal (z³=1) ---")
    for (x, y) in [(0.5, 0.2), (-0.4, 0.7), (-0.4, -0.7)]
        root = newton_fractal(x, y; max_iter=30)
        names = ["none", "1", "e^{2πi/3}", "e^{4πi/3}"]
        println("  z₀=$x+$(y)i → root $(names[root+1])")
    end

    # Sierpinski
    pts = sierpinski_triangle(1000)
    println("\n--- Sierpinski Triangle ---")
    println("  1000 points generated")
    println("  Theoretical dimension: D = log(3)/log(2) = $(round(log(3)/log(2),4))")

    println("\n✅ fractal_gallery.jl done.")
end
main()
