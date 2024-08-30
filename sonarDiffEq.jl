include("soundFormulas.jl")


"""
    ΔGenericRayAttribute(du, u, p, t)

    Function for calculating the properties of sound rays using differential equations.
    Modifies the content of "du" by storing the computed data in it.
"""
function ΔGenericRayAttribute!(du, u, p, t)
    r, z, ξ, ζ = u; 
    f, ∂₁c, ∂₂c, D, NA = p; 
    
    l = LayerData(2.0,34.7,z);

    c = f(l);

    α = -1.0/(c^2);     
    
    du[1] = c * ξ;
    du[2] = c * ζ;
    du[3] = α * ∂₁c(l);
    du[4] = α * ∂₂c(l);

    zᵢ₊₁ = z + du[2];

    if (zᵢ₊₁ > D)
        h = (D - z) / (c * ζ) 
        t = p[5] + h;

        du[2] = h * c * ζ;
        du[4] = 0.0;
        u[4] *= -1;
    end

    if (zᵢ₊₁ < 0.0)
        h = (0 - z) / (c * ζ) 
        t = p[5] + h;

        du[2] = h * c * ζ;
        du[4] = 0.0;
        u[4] *= -1;
    end
    p[5] = t;
end

"""
    ΔGenericRayAttribute!(u,p,t)

    Returns "du" after computing it's content.
"""
function ΔGenericRayAttribute!(u, p)
    r, z, ξ, ζ = u; 
    f, ∂₁c, ∂₂c, D, NA = p; 
    
    l = LayerData(2.0,34.7,z);

    c = f(l);

    α = -1.0/(c^2);     
    
    du = Vector{Float64}(0.0,4);
    du[1] = c * ξ;
    du[2] = c * ζ;
    du[3] = α * ∂₁c(l);
    du[4] = α * ∂₂c(l);

    zᵢ₊₁ = z + du[2];

    if (zᵢ₊₁ > D)
        h = (D - z) / (c * ζ) 
        t = p[5] + h;

        du[2] = h * c * ζ;
        du[4] = 0.0;
        u[4] *= -1;
    end

    if (zᵢ₊₁ < 0.0)
        h = (0 - z) / (c * ζ) 
        t = p[5] + h;

        du[2] = h * c * ζ;
        du[4] = 0.0;
        u[4] *= -1;
    end
    p[5] = t;

    return du;
end

"""
    InitGenericRay(r₀ :: Float32, z₀ :: Float32, θ₀ :: Float32)

    Initializes a GenericRayData structure from a starting position and initial takeoff angle.
"""
function InitGenericRay(r₀ :: AbstractFloat, z₀ :: AbstractFloat, θ₀ :: AbstractFloat)
    l = LayerData(2.0,34.7,z₀)
    c₀ = CoppensFormula(l);

    ξ = cos(θ₀) / c₀; 
    ζ = sin(θ₀) / c₀;
    
    return GenericRayData(r₀,z₀,ξ,ζ,c₀);
end