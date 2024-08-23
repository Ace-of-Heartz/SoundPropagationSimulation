include("sonarForms.jl");
#

"""
    Data for storing a ray, the elapsed time between each step, and the maximum number of steps. 
"""
struct SimulationData
    ray :: RayData

    stepSize :: Float32 #Time between each step 
    maxStep  :: Int32
end    

"""
    ConstantForwardPropFunc(ray :: RayData, step :: Float32) :: RayData

    Function used for debugging simple raymarching
"""
function ConstantForwardPropFunc(ray :: RayData, step :: Float32) :: RayData
    propVec = ray.speed * step .* [cos(ray.angle),sin(ray.angle),0.0];

    pos₁ = ray.position + propVec;


    RayData(pos₁, ray.angle,ray.speed);
end

function March(data :: SimulationData)
    ray = data.ray;
    
    rays = Vector{RayData}(undef,data.maxStep);
    timeVals = Vector{Float32}(undef,data.maxStep);

    for i = 1:data.maxStep 
        rays[i] = ray;
        timeVals[i] = data.stepSize * i;

        ray = ConstantForwardPropFunc(ray,data.stepSize);      
    end 
    return (timeVals,rays);
end


"""
    EulerMethod(data :: SimulationData)

    Implementation of Euler's method for calculating the propagation of sound in sea water.
    
"""
function EulerMethod(data :: SimulationData)
    ray = data.ray;
    
    rays = Vector{RayData}(undef,data.maxStep);
    timeVals = Vector{Float32}(undef,data.maxStep);

    for i = 1:data.maxStep 
        rays[i] = ray;
        timeVals[i] = data.stepSize * i;

        ray = ApplyPropagationChange(ray,ΔPropagation(ray,data.stepSize),Float32(1.0));      
    end 
    return (timeVals,rays);
end


"""
    HeunMethod(data :: SimulationData)

    Implementation of Heun's method for calculating the propagation of sound in sea water.

"""
function HeunMethod(data :: SimulationData)
    ray = data.ray;

    rays = Vector{RayData}(undef, data.maxStep);
    timeVals = Vector{Float32}(undef, data.maxStep);

    for i = 1:data.maxStep
        rays[i] = ray;
        timeVals[i] = data.stepSize * i;
        
        intermidiateRay = ApplyPropagationChange(ray,ΔPropagation(ray, data.stepSize),1.0);
        
        Δray₀ = ΔPropagation(ray,data.stepSize / 2);
        Δray₁ = ΔPropagation(intermidiateRay,data.stepSize / 2); # ??? tᵢ₊₁ = tᵢ + h ???

        ray = ApplyPropagationChange(ApplyPropagationChange(ray,Δray₀,1.0),Δray₁,1.0);
    end
    return (timeVals,rays)
end

"""
    RK45(data :: SimulationData)

    Implementation of Runge Kutta 4(5) for calculating the propagation of sound in sea water.
"""
function RK45(data :: SimulationData)
    ray = data.ray;

    rays = Vector{RayData}(undef, data.maxStep);
    timeVals = Vector{Float32}(undef, data.maxStep);

    coMatrix :: Matrix{Float32} = [
        0.0 0.0 0.0 0.0; 
        1.0/2.0 1.0/2.0 0.0 0.0;
        1.0/2.0 0.0 1.0/2.0 0.0;
        1.0 0.0 0.0 1.0
    ];

    weightVec :: Vector{Float32} = [ 1.0/6.0 , 1.0/3.0 , 1.0/3.0 , 1.0/6.0];

    for i = 1:data.maxStep
        rays[i] = ray;
        timeVals[i] = data.stepSize * i;

        ray = GeneralRKMethod(data,coMatrix,weightVec);
    end

    return (timeVals,rays);
end

"""
    GeneralRKMethod(data :: SimulationData, coMatrix :: Matrix{Float32}, weightVec :: Vector{Float32})

    Calculates the next ray's attributes using a general coefficient matrix 
    and a vector of weights defined by various Runge Kutta methods.
"""
function GeneralRKMethod(data :: SimulationData, coMatrix :: Matrix{Float32}, weightVec :: Vector{Float32}) 
    cache = ComputeRKIntermidiateValues(data, coMatrix);
        
    s = length(weightVec);

    accΔRay = AccPropagationChange(cache,weightVec);

    return ApplyPropagationChange(data.ray,accΔRay,Float32(1.0));
end

"""
    ComputeRKIntermidiateValues(data :: SimulationData , coMatrix :: Matrix{Float32})

    Calculates the intermidiate changes of the ray attributes for further computations.
"""
function ComputeRKIntermidiateValues(data :: SimulationData , coMatrix :: Matrix{Float32})
    
    (m,n) = (size(coMatrix,1),size(coMatrix,2));
    if(n != m)
        error("Bad matrix dimensions!");
    end

    cache = Vector{RayChangeData}(undef,m);

    cache[1] = ΔPropagation(data.ray,data.stepSize + data.stepSize * coMatrix[1, 1]);

    for i = 2:m

        for j = 2:i
            cache[i] = ΔPropagation(
                ApplyPropagationChange(
                    data.ray,
                    cache[j-1],
                    coMatrix[i, j]),
                data.stepSize * coMatrix[i, 1]
                );
        end

    end

    return cache;
end