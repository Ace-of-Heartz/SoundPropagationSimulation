
"""
    Data for storing a ray's attributes
"""
struct RayData 
    position :: Vector{Float32}
    angle    :: Float32     
    speed    :: Float32
end


"""
    Data for storing the change in a ray's attributes
"""
struct RayChangeData 
    ΔDistance :: Vector{Float32}
    ΔAngle    :: Float32 
    ΔSpeed    :: Float32
end


"""
    Data for calculating the speed of sound in a given layer of sea water. 
"""
struct LayerData
    temperature :: Float32 
    salinity    :: Float32
    depth       :: Float32

end



"""
    Compute the speed of sound in sea water using Coppen's formula
"""
function CoppensFormula(layerData :: LayerData)
    d = layerData.depth / 1000.0;
    s = layerData.salinity / 1000.0;
    t = layerData.temperature / 10;

    c₀ = CoppensFormulaBase(LayerData(layerData.temperature,layerData.salinity,0.0));

    c₀ + (16.23 + 0.253 * t) * d + (0.213 - 0.1 * t) * d^2 + (0.016 + 0.0002 * (s - 35.0)) * (s - 35.0) * t * d; 
end

"""
    Compute the speed of sound in sea water using Coppen's formula at 0 depth
"""
function CoppensFormulaBase(layerData :: LayerData)
    d = layerData.depth / 1000.0     # Depth in kilometers
    s = layerData.salinity / 1000.0  # Salinity in parts per thousand
    t = layerData.temperature / 10;  # temperature

    1449.05 + 45.7 * t - 5.21 * t^2 + 0.23 * t^3 + (1.333 - 0.126 * t + 0.009 * t^2) * (s - 35.0)
end


"""
    Accumulate a vector of change in ray attributes with weights. 
"""
function AccPropagationChange(Δrays :: Vector{RayChangeData}, weightVec :: Vector{Float32}) :: RayChangeData
    accΔProp  = [0.0,0.0,0.0];
    accΔAngle = 0.0;
    accΔSpeed = 0.0;
    
    for i = 1:length(Δrays)
        accΔProp  += Δrays[i].ΔDistance * weightVec;
        accΔAngle += Δrays[i] * weightVec;
        accΔSpeed += Δrays[i] * weightVec;
    end

    return RayChangeData(accΔProp,accΔAngle,accΔSpeed);
end

"""
    Applies the calculated changes in ray's attribute to a given ray. 
"""
function ApplyPropagationChange(ray₀ :: RayData, Δray :: RayChangeData, w :: Float32)::RayData
    ray₁ = RayData(
        ray₀.position + Δray.ΔDistance * w,
        ray₀.angle + Δray.ΔAngle * w,
        ray₀.speed + Δray.ΔSpeed * w
        );
end


"""
    Calculates the changes in ray's attributes given the elapsed time between steps. 
"""
function ΔPropagation(ray::RayData, step::Float32)::RayChangeData
    ξ = cos(ray.angle) / ray.speed;

    dirVec  = [cos(ray.angle),sin(ray.angle), 0.0];
    propVec = ray.speed * step .* dirVec;

    depth₀ = ray.position[2];
    speed₀  = ray.speed;

    # This feels weird . . . 
    # Normally you would get the new depth information from raytracing, and intersection certain layers of water
    # So please review this
    tempPos = ray.position + propVec;
    depth₁  = tempPos[2];

    Δy  = (ray.position + propVec)[2] - depth₀  

    layerData = LayerData(2.0,34.7,depth₁); # TODO: Make this not constant

    speed₁ = CoppensFormula(layerData);
    ΔSpeed = speed₁ - speed₀

    v₀ = ΔSpeed / Δy;

    α = ξ^2 * speed₀^2;
    β = ξ^2 * speed₁^2;

    if (β >= 1)
        Δx = 2 / (ξ * v₀) * sqrt(1-α)
    else
        Δx = 1/(ξ * v₀) * ( sqrt(1-α) - sqrt(1-β) )
    end

    ΔDistance = Vector([Δx,Δy,0.0])

    localθ₀ = atan((tempPos[2] - ray.position[2])/ (tempPos[1] - ray.position[1]));
    localθ₁ = atan(Δy / Δx);

    Δθ = localθ₁ - localθ₀;

    ΔRay = RayChangeData(ΔDistance,Δθ,ΔSpeed);
    print("
        $α\t 
        $β\t
        $tempPos\t
        $Δθ\t 
        \n  
    ")

    return ΔRay;
end