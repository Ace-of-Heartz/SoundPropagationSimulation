include("rkData.jl");
include("soundFormulas.jl")
include("utils.jl");


"""
    Data for storing a ray, the elapsed time between each step, and the maximum number of steps. 
"""
struct SimulationData
    ray :: GenericRayData

    stepSize :: Float32 #Time between each step 
    maxStep  :: Int32
end    

"""
    RK45(data :: SimulationData)

    Implementation of Runge Kutta 4(5) for calculating the propagation of sound in sea water.
"""
function RKMethod(data :: SimulationData, rkData :: RKData)
    ray = data.ray;

    rays = Vector{GenericRayData}(undef, data.maxStep);
    arclength = Float32(0.0);
    arclengths = []

    for i = 1:data.maxStep
        arclengths = [arclengths ; arclength];
        rays[i] = ray;

        arclength += data.stepSize; 
        ray = GeneralRKMethod(ray,data.stepSize,rkData);
        
    end

    return (arclengths,rays);
end 

# TODO: Test this out
function RKMethod(
    diffEq :: Function,
    u₀,
    stepSize,
    maxStep,
    p ,
    rkData :: RKData
)
    u = u₀
    du = Vector{<:AbstractFloat}(0.0,4);
    
    s = 0.0;
    
    arclengths = [];
    rayAttributes = [];

    for i = 0:maxStep 
        arclengths = [arclengths ; s];
        rayAttributes = [rayAttributes; u];

        u = ComputeRKNextValue(diffEq,u,p,rkData,stepSize);
        s += stepSize;
    end

    return (arclengths,rays);
end




"""
    GeneralRKMethod(data :: SimulationData, coMatrix :: Matrix{Float32}, weightVec :: Vector{Float32})

    Calculates the next ray's attributes using a general coefficient matrix 
    and a vector of weights defined by various Runge Kutta methods.
"""
function GeneralRKMethod(ray :: GenericRayData,stepSize :: Float32, rkData :: RKData) 
    intValues = ComputeRKIntermidiateValues(ray,stepSize, rkData.comatrix);
    ΔRay = SumRayAttributeChanges(intValues,rkData.weightVector);

    return ApplyRayAttributesChanges(ray,ΔRay);
end

"""
    ComputeRKIntermidiateValues(data :: SimulationData , coMatrix :: Matrix{Float32})

    Calculates the intermidiate changes of the ray attributes for further computations.
"""
function ComputeRKIntermidiateValues(ray :: GenericRayData, stepSize :: Float32 , coMatrix :: Matrix{Float32})
    
    (m,n) = (size(coMatrix,1),size(coMatrix,2));
    if(n != m)
        error("Bad matrix dimensions!");
    end

    intValues = Vector{GenericRayDataChange}(undef,m);

    intValues[1] = ΔRayAttributes(stepSize + stepSize * coMatrix[1,1],ray) #This is probably not right, recalc this please

    for i = 2:m

        for j = 2:i
            intValues[i] = ΔRayAttributes(
                stepSize + stepSize * coMatrix[i, 1],
                ApplyRayAttributesChanges(
                    ray,
                    intValues[j-1],
                    coMatrix[i, j])
                );
        end
    end

    return intValues;
end

function ComputeRKNextValue(diffEq :: Function, u, p, rkData ::RKData, stepSize :: AbstractFloat)
    (m,n) = (size(coMatrix,1),size(coMatrix,2));
    if(m != n-1)
        error("Bad matrix dimensions!");
    end

    ks = Vector{GenericRayDataChange}(undef,m);

    du = diffEq(u,p);
    ks[1] = du;

    # Because our differential equation won't depend on arclength, we can forgo the entire section about using different stepsizes
    for i = 2:m
        kᵢ = u;
        for j = 1:(i-1)
            kᵢ += ks[j] * rkData.comatrix[i,j]
        end
        ks[i] = kᵢ;
        du = diffEq(u,p,t) .* stepSize; 
    end

    uᵢ₊₁ = u;
    for i in eachindex(rkData.weightVector)
        uᵢ₊₁ += ks[i] * rkData.weightVector[i];
    end

    return uᵢ₊₁;
end