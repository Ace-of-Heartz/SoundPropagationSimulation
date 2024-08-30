include("rayData.jl");

function ApplyRayAttributesChanges(ray :: GenericRayData, Δray ::GenericRayDataChange) :: GenericRayData
    return GenericRayData(
        ray.r + Δray.r,
        ray.z + Δray.z,
        ray.ξ + Δray.ξ,
        ray.ζ + Δray.ζ,
        ray.c + Δray.c 
    )
end

function ApplyRayAttributesChanges(ray :: GenericRayData, Δray ::GenericRayDataChange, w :: Float32) :: GenericRayData
    return GenericRayData(
        ray.r + Δray.r * w,
        ray.z + Δray.z * w,
        ray.ξ + Δray.ξ * w,
        ray.ζ + Δray.ζ * w,
        ray.c + Δray.c * w 
    )
end

function SumRayAttributeChanges(Δrays :: Vector{GenericRayDataChange}, ws :: Vector{Float32}) :: GenericRayDataChange 

    Δr :: Float32 = 0.0
    Δz :: Float32 = 0.0
    Δc :: Float32 = 0.0
    Δξ :: Float32 = 0.0
    Δζ :: Float32 = 0.0
    for ray in Δrays
        Δr = ray.r;
        Δz = ray.z;
        Δc = ray.c;
        Δξ = ray.ξ;
        Δζ = ray.ζ;
    end
    ΔsumRay = GenericRayDataChange(Δr,Δz,Δξ,Δζ,Δc); 


    return ΔsumRay;
end

function CylinderToCartesian(pos :: Tuple{<:AbstractFloat, <: AbstractFloat}) :: Tuple{<:AbstractFloat, <:AbstractFloat, <:AbstractFloat}
    (r, u) = pos;
    return (cos(r),sin(r),u); 
end

