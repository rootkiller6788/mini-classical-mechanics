# fractals.jl — 分形: Mandelbrot, Julia, Newton分形, Sierpinski, 盒计数维数
# 参考: Peitgen & Richter, Falconer "Fractal Geometry"

"----------- Mandelbrot / Julia -----------"
"Mandelbrot 集迭代数 (c固定, z₀=0)"
function mandelbrot_iter(cx::Float64, cy::Float64; max_iter=100, escape=2.0)
    x=y=0.0
    for n in 1:max_iter; xt=x*x-y*y+cx; y=2*x*y+cy; x=xt; x*x+y*y>escape^2 && return n; end
    max_iter
end

"Julia 集迭代数 (c固定, z₀可变)"
function julia_iter(zx::Float64, zy::Float64, cx::Float64, cy::Float64; max_iter=100, escape=2.0)
    x,y=zx,zy
    for n in 1:max_iter; xt=x*x-y*y+cx; y=2*x*y+cy; x=xt; x*x+y*y>escape^2 && return n; end
    max_iter
end

"Mandelbrot 集图像 (矩阵)"
function mandelbrot_set(x_range, y_range; nx=400, ny=400, max_iter=100)
    xs=range(x_range[1],x_range[2],length=nx)
    ys=range(y_range[1],y_range[2],length=ny)
    img=zeros(Int,ny,nx)
    for (i,x) in enumerate(xs), (j,y) in enumerate(ys); img[j,i]=mandelbrot_iter(x,y;max_iter=max_iter); end
    xs,ys,img
end

"Julia 集图像"
function julia_set(x_range, y_range, cx, cy; nx=400, ny=400, max_iter=100)
    xs=range(x_range[1],x_range[2],length=nx)
    ys=range(y_range[1],y_range[2],length=ny)
    img=zeros(Int,ny,nx)
    for (i,x) in enumerate(xs), (j,y) in enumerate(ys); img[j,i]=julia_iter(x,y,cx,cy;max_iter=max_iter); end
    xs,ys,img
end

"距离估计法 Mandelbrot (平滑着色)"
function mandelbrot_distance(cx::Float64, cy::Float64; max_iter=200, escape=2.0)
    x=y=dx=dy=0.0
    for n in 1:max_iter
        xt=x*x-y*y+cx; yt=2*x*y+cy
        dxt=2*(x*dx-y*dy)+1; dyt=2*(x*dy+y*dx)
        x,y=xt,yt; dx,dy=dxt,dyt
        r=sqrt(x*x+y*y)
        r>escape && return 0.5*log(r)*r/sqrt(dx*dx+dy*dy)
    end
    return 0.0
end

"----------- Newton 分形 -----------"
"Newton分形: z_{n+1}=z_n-(z³-1)/(3z²), 按收敛到哪个根上色"
function newton_fractal(x::Float64, y::Float64; max_iter=50, tol=1e-6)
    z=complex(x,y); roots=[complex(1,0),complex(-0.5,sqrt(3)/2),complex(-0.5,-sqrt(3)/2)]
    for n in 1:max_iter; z2=z*z; z=z-(z2*z-1)/(3*z2)
        for (k,root) in enumerate(roots); abs(z-root)<tol && return k; end
    end
    0
end

"Newton分形图像"
function newton_fractal_image(x_range, y_range; nx=400, ny=400, max_iter=50)
    xs=range(x_range[1],x_range[2],length=nx)
    ys=range(y_range[1],y_range[2],length=ny)
    img=zeros(Int,ny,nx)
    for (i,x) in enumerate(xs), (j,y) in enumerate(ys); img[j,i]=newton_fractal(x,y;max_iter=max_iter); end
    xs,ys,img
end

"----------- 经典分形 -----------"
"Sierpinski 三角 (混沌游戏法)"
function sierpinski_triangle(n_points::Int)
    vertices=[[0.0,0.0],[1.0,0.0],[0.5,sqrt(3)/2]]
    p=[0.5,0.5]; points=[copy(p)]
    for _ in 1:n_points
        v=vertices[rand(1:3)]; p=[(p[1]+v[1])/2,(p[2]+v[2])/2]; push!(points,copy(p))
    end
    points
end

"Koch 雪花 (L-system 近似)"
function koch_curve(n_iter::Int)
    segments=[[(0.0,0.0),(1.0,0.0)]]
    for _ in 1:n_iter
        new_segs=[]
        for ((x1,y1),(x2,y2)) in segments
            dx,dy=x2-x1,y2-y1
            a=(x1+dx/3,y1+dy/3); b=(x1+dx/2-dy*sqrt(3)/6,y1+dy/2+dx*sqrt(3)/6)
            c=(x1+2dx/3,y1+2dy/3)
            push!(new_segs,((x1,y1),a),(a,b),(b,c),(c,(x2,y2)))
        end
        segments=new_segs
    end
    [s[1] for s in segments]; push!([s[1] for s in segments], segments[end][2])
end

"----------- 分形维数 -----------"
"盒计数维数"
function box_counting_dimension(points::Vector{Vector{Float64}}, epsilons::Vector{Float64})
    Ns=Float64[]
    for eps_val in epsilons
        boxes=Set{Tuple{Int,Int}}()
        for p in points; push!(boxes,(Int(floor(p[1]/eps_val)),Int(floor(p[2]/eps_val)))); end
        push!(Ns,Float64(length(boxes)))
    end
    x_vals=log.(1 ./ epsilons); y_vals=log.(Ns)
    n=length(x_vals)
    (n*sum(x_vals.*y_vals)-sum(x_vals)*sum(y_vals))/(n*sum(x_vals.^2)-sum(x_vals)^2)
end

"单一 eps 维数估计"
function box_counting_dimension_single(points::Vector{Vector{Float64}}; n_eps=20)
    xmin=minimum(p[1] for p in points); xmax=maximum(p[1] for p in points)
    ymin=minimum(p[2] for p in points); ymax=maximum(p[2] for p in points)
    L=max(xmax-xmin,ymax-ymin)
    epsilons=[L/2^i for i in 2:n_eps+1 if L/2^i > 1e-10]
    length(epsilons) < 3 && return 0.0
    box_counting_dimension(points, epsilons)
end

"Koch 曲线的理论维数: D = log(4)/log(3)"
koch_dimension() = log(4)/log(3)

"Sierpinski 三角理论维数: D = log(3)/log(2)"
sierpinski_dimension() = log(3)/log(2)

export mandelbrot_iter, julia_iter, mandelbrot_set, julia_set, mandelbrot_distance
export newton_fractal, newton_fractal_image
export sierpinski_triangle, koch_curve
export box_counting_dimension, box_counting_dimension_single, koch_dimension, sierpinski_dimension
