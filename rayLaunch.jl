include("rkMethods.jl");
include("sonarDiffEq.jl")
include("plot.jl")

using GLMakie;
using DifferentialEquations;

function LaunchSingleGenRayWithODE(
    r₀ :: AbstractFloat,
    z₀ :: AbstractFloat,
    θ₀ :: AbstractFloat,
    tspan :: Tuple{<:AbstractFloat,<:AbstractFloat},
    f :: Function,
    df₁ :: Function, 
    df₂ :: Function  
)
    ray = InitGenericRay(r₀, z₀, θ₀);
    u₀ = [ray.r, ray.z, ray.ξ, ray.ζ];

    p = [f,df₁,df₂, 6400.0, 0.0];
    
    prob = ODEProblem(ΔGenericRayAttribute!,u₀,tspan,p);
    sol = solve(prob,Tsit5(),reltol = 1e-8, abstol = 1e-8);

    ss = sol.t;
    rs = map(x -> x[1], sol.u);
    zs = map(x -> x[2], sol.u);
    ξs = map(x -> x[3], sol.u);
    ζs = map(x -> x[4], sol.u);
    return (ss,rs,zs,ξs,ζs)
end

function LaunchSingleGenRay(
    r₀ :: AbstractFloat,
    z₀ :: AbstractFloat,
    θ₀ :: AbstractFloat,
    tspan :: Tuple{<:AbstractFloat,<:AbstractFloat},
    stepSize :: AbstractFloat,
    maxStep :: Integer,
    rkData :: RKData, 
    f :: Function,
    df₁ :: Function, 
    df₂ :: Function  
)
end

function LaunchGenRaysWithODE(
    angleRange ::AbstractFloat,
    angleOffset :: AbstractFloat,
    r₀ :: AbstractFloat, z₀ :: AbstractFloat,
    tspan :: Tuple{<:AbstractFloat, <:AbstractFloat},
    numberOfRays :: Integer, 
    f :: Function,
    df₁ :: Function, 
    df₂ :: Function 
)
    angleStep = angleRange / numberOfRays;

    results = []

    for i = 1:numberOfRays
        
        θ₀ = angleStep * i + angleOffset;

        res = LaunchSingleGenRayWithODE(
            r₀, z₀, θ₀, tspan, f, df₁, df₂
        );
        results = [ results ; res];
    end
    return results;
end

function LaunchGenRays(
    angleRange ::AbstractFloat,
    angleOffset :: AbstractFloat,
    r₀ :: AbstractFloat, z₀ :: AbstractFloat,
    numberOfRays :: Integer, 
    stepSize :: AbstractFloat, maxStep :: Integer,
    rkData :: RKData)
    angleStep = angleRange / numberOfRays
    
    res = [];

    for i = 0:numberOfRays
        (arclengths,rays) = LaunchSingleGenRay(r₀,z₀,angleStep * (i) + angleOffset,stepSize,maxStep,rkData) 

        res = [res; rays];
    end
    
    rs = map(x -> x.r,res);
    zs = map(x -> x.z,res);
    plot(rs,zs);
end

