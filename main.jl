# using DifferentialEquations;
import Pkg;
Pkg.add("GLMakie")
using GLMakie;
using LinearAlgebra;

# struct RayData 
#     position :: Vector{Float32}
#     angle    :: Float32     
#     speed    :: Float32
# end

# struct RayChangeData 
#     ΔDistance :: Vector{Float32}
#     ΔAngle    :: Float32 
#     ΔSpeed    :: Float32
# end

# struct LayerData
#     temperature :: Float32 
#     salinity    :: Float32
#     depth       :: Float32

# end

"""
    Data for storing a ray, the elapsed time between each step, and the maximum number of steps. 
"""
struct SimulationData
    ray :: RayData

    stepSize :: Float32 #Time between each step 
    maxStep  :: Int32
end    


"""
    Data for holding different values to be visualized using plots.
"""
struct PlotData 
    speedVals           :: Vector{Float32}
    positionVals        :: Vector{Vector{Float32}}
    angleVals           :: Vector{Float32}
    horizontalPropVals  :: Vector{Float32}
    verticalPropVals    :: Vector{Float32}
    timeVals            :: Vector{Float32}
end

function PreparePlotData(rays :: Vector{RayData},timeVals :: Vector{Float32}) :: PlotData
    speedVals    = map(x -> (x.speed), rays);
    positionVals = map(x -> (x.position), rays);
    angleVals    = map(x -> (x.angle), rays);
    
    changeAmount = size(positionVals,1) - 1;
    
    horizontalPropVals = Vector{Float32}(undef, changeAmount);
    verticalPropVals   = Vector{Float32}(undef, changeAmount);
    for i = 1:(changeAmount-1)
        horizontalPropVals[i] = positionVals[i+1][1] - positionVals[i][1];
        verticalPropVals[i]   = positionVals[i+1][2] - positionVals[i][2];
    end 

    return PlotData(speedVals,positionVals,angleVals,horizontalPropVals,verticalPropVals,timeVals);
end

function PlotData(plotDatas)
    

    xPlots = [];
    yPlots = [];

    fig = Figure();

    titles = [
        ("Speed Relative to Elapsed Time"),
        ("Speed Relative to Depth")
    ] 

    labels = [
        ("Elapsed Time (s)","Speed (m/s)"),
        ("Depth (m)","Speed (m/s)")
    ];


    Axis(fig[1,1], title = "Ray Paths")

    for i in eachindex(labels)
        for j in eachindex(plotDatas)
            Axis(fig[2,1][j, i], title = titles[i],
                xlabel = labels[i][1],
                ylabel = labels[i][2]
            );
        end
    end

    
    for i in eachindex(plotDatas)
        plotData = plotDatas[i];

        lines!(
            fig[2,1][i,2],
            plotData.timeVals,
            plotData.speedVals,
            );
        lines!(
            fig[2,1][i,1],
            map(x -> x[2],plotData.positionVals),
            plotData.speedVals,
        );
        
        append!(xPlots,map(x -> x[1],plotData.positionVals));
        append!(yPlots,map(x -> x[2],plotData.positionVals)); 
    end
    
    scatter!(
        fig[1,1],
        xPlots,
        yPlots,
    );
    return fig;
end


# function EulerMethod(data :: SimulationData)
#     ray = data.ray;
    
#     rays = Vector{RayData}(undef,data.maxStep);
#     timeVals = Vector{Float32}(undef,data.maxStep);

#     for i = 1:data.maxStep 
#         rays[i] = ray;
#         timeVals[i] = data.stepSize * i;

#         ray = ApplyPropagationChange(ray,ΔPropagation(ray,data.stepSize));      
#     end 
#     return (timeVals,rays);
# end

# function HeunMethod(data :: SimulationData)
#     ray = data.ray;

#     rays = Vector{RayData}(undef, data.maxStep);
#     timeVals = Vector{Float32}(undef, data.maxStep);

#     for i = 1:data.maxStep
#         rays[i] = ray;
#         timeVals[i] = data.stepSize * i;
        
#         intermidiateRay = ApplyPropagationChange(ray,ΔPropagation(ray, data.stepSize));
        
#         Δray₀ = ΔPropagation(ray,data.stepSize / 2);
#         Δray₁ = ΔPropagation(intermidiateRay,data.stepSize / 2); # ??? tᵢ₊₁ = tᵢ + h ???

#         ray = ApplyPropagationChange(ApplyPropagationChange(ray,Δray₀),Δray₁);
#     end
#     return (timeVals,rays)
# end


# function RK45(data :: SimulationData)
#     ray = data.ray;

#     rays = Vector{RayData}(undef, data.maxStep);
#     timeVals = Vector{Float32}(undef, data.maxStep);

#     coMatrix = [
#         0.0 0.0 0.0 0.0; 
#         1.0/2.0 1/2.0 0.0 0.0;
#         1.0/2.0 0.0 1.0/2.0 0.0;
#         1.0 0.0 0.0 1.0
#     ];

#     weightVec = [ 1/6 1/3 1/3 1/6];

#     for i = 1:data.maxStep
#         rays[i] = ray;
#         timeVals[i] = data.stepSize * i;

#         ray = GeneralRKMethod(data,coMatrix,weightVec);
#     end

#     return (timeVals,rays);
# end

# function GeneralRKMethod(data :: SimulationData, coMatrix :: Matrix{Float32}, weightVec :: Vector{Float32}) 
#     cache = ComputeRKIntermidiateValues(data, coMatrix);
        
#     s = length(weightVec);

#     accΔRay = AccPropagationChange(cache,weightVec);

#     return ApplyPropagationChange(data.ray,accΔRay,1.0);
# end

# function ComputeRKIntermidiateValues(data :: SimulationData , coMatrix :: Matrix{Float32})
    
#     (m,n) = (size(coMatrix,1),size(coMatrix,2));
#     if(n != m)
#         error("Bad matrix dimensions!");
#     end

#     cache = Vector{Float32}(undef,m);

#     cache[1] = ΔPropagation(data.ray,data.stepSize * coMatrix[1][1]);


#     for i = 2:m

#         for j = 2:n
#             cache[i] = ΔPropagation(
#                 ApplyPropagationChange(
#                     data.ray,
#                     cache[j],
#                     coMatrix[i][j]),
#                 data.stepSize * coMatrix[i][1]
#                 );
#         end

#     end

#     return cache;
# end

# function AccPropagationChange(Δrays :: Vector{RayChangeData}, weightVec :: Vector{Float32})
#     accΔProp  = [0.0,0.0,0.0];
#     accΔAngle = 0.0;
#     accΔSpeed = 0.0;
    
#     for i = 1:length(Δrays)
#         accΔProp  += Δrays[i].ΔDistance * weightVec;
#         accΔAngle += Δrays[i] * weightVec;
#         accΔSpeed += Δrays[i] * weightVec;
#     end

#     return RayChangeData(accΔProp,accΔAngle,accΔSpeed);
# end

# function ApplyPropagationChange(ray₀ :: RayData, Δray :: RayChangeData, w :: Float32)::RayData
#     ray₁ = RayData(
#         ray₀.position + Δray.ΔDistance * w,
#         ray₀.angle + Δray.ΔAngle * w,
#         ray₀.speed + Δray.ΔSpeed * w
#         );
# end

# function ΔPropagation(ray::RayData, step::Float32)::RayChangeData
#     ξ = cos(ray.angle) / ray.speed;

#     dirVec  = [cos(ray.angle),sin(ray.angle), 0.0];
#     propVec = ray.speed * step .* dirVec;

#     depth₀ = ray.position[2];
#     speed₀  = ray.speed;

#     # This feels weird . . . 
#     # Normally you would get the new depth information from raytracing, and intersection certain layers of water
#     # So please review this
#     tempPos = ray.position + propVec;
#     depth₁  = tempPos[2];

#     Δy  = (ray.position + propVec)[2] - depth₀  

#     layerData = LayerData(2.0,34.7,depth₁); # TODO: Make this not constant

#     speed₁ = CoppensFormula(layerData);
#     ΔSpeed = speed₁ - speed₀

#     v₀ = ΔSpeed / Δy;

#     α = ξ^2 * speed₀^2;
#     β = ξ^2 * speed₁^2;

#     if (β >= 1)
#         Δx = 2 / (ξ * v₀) * sqrt(1-α)
#     else
#         Δx = 1/(ξ * v₀) * ( sqrt(1-α) - sqrt(1-β) )
#     end

#     ΔDistance = Vector([Δx,Δy,0.0])

#     localθ₀ = atan((tempPos[2] - ray.position[2])/ (tempPos[1] - ray.position[1]));
#     localθ₁ = atan(Δy / Δx);

#     Δθ = localθ₁ - localθ₀;

#     ΔRay = RayChangeData(ΔDistance,Δθ,ΔSpeed);
#     print("
#         $α\t 
#         $β\t
#         $tempPos\t
#         $Δθ\t 
#         \n  
#     ")

#     return ΔRay;
# end

# function ConstantForwardPropFunc(ray :: RayData, step :: Float32)::RayData
#     propVec = ray.speed * step .* [cos(ray.angle),sin(ray.angle),0.0];

#     pos₁ = ray.position + propVec;


#     RayData(pos₁, ray.angle,ray.speed);
# end



# function CoppensFormula(layerData :: LayerData)
#     d = layerData.depth / 1000.0;
#     s = layerData.salinity / 1000.0;
#     t = layerData.temperature / 10;

#     c₀ = CoppensFormulaBase(LayerData(layerData.temperature,layerData.salinity,0.0));

#     c₀ + (16.23 + 0.253 * t) * d + (0.213 - 0.1 * t) * d^2 + (0.016 + 0.0002 * (s - 35.0)) * (s - 35.0) * t * d; 
# end

# function CoppensFormulaBase(layerData :: LayerData)
#     d = layerData.depth / 1000.0     # Depth in kilometers
#     s = layerData.salinity / 1000.0  # Salinity in parts per thousand
#     t = layerData.temperature / 10;  # temperature

#     1449.05 + 45.7 * t - 5.21 * t^2 + 0.23 * t^3 + (1.333 - 0.126 * t + 0.009 * t^2) * (s - 35.0)
# end


function ShootRays(
    numberOfRays :: Int,
    timeInterval :: Float64,
    maxSteps :: Int,
    initPos :: Vector{Float64} = [0.0,550.0,0.0]
    )
    
    if (numberOfRays == 0)
        error("Number of rays must be at least 1")
    end

    interval = pi / numberOfRays; 


    startingLayer = LayerData(2.0,34.7,250.0);
    startingSpeed = CoppensFormula(startingLayer);

    xPosGlob = [];
    yPosGlob = [];

    rays = [];
    plotDatas = [];

    for i = 0:(numberOfRays-1)

        offset = pi/2.0
        angle = i * interval - offset;
        
        angle = abs(angle) < 0.00001 ? 0.00001 : angle

        print("{$angle}\n");
        ray = RayData(initPos,angle,startingSpeed);
        simData = SimulationData(ray,timeInterval,maxSteps);
        
        (timeVals,rays) = EulerMethod(simData);
        
        push!(plotDatas,PreparePlotData(rays,timeVals));
    end

    print("Simulation finished...")
    
    PlotData(plotDatas);
end

