include("soundData.jl");

"""
    Compute the speed of sound in sea water using Coppen's formula
"""
function CoppensFormula(layerData :: LayerData)
    d = layerData.depth / 1000.0;
    s = layerData.salinity / 1000.0;
    t = layerData.temperature / 10;

    c₀ = CoppensFormulaBase(LayerData(layerData.temperature,layerData.salinity,0.0));

    return c₀ + (16.23 + 0.253 * t) * d + (0.213 - 0.1 * t) * d^2 + (0.016 + 0.0002 * (s - 35.0)) * (s - 35.0) * t * d; 
end

function ∂CoppensOnR(layerData :: LayerData)
    return 0; 
end

function ∂CoppensOnZ(layerData :: LayerData)
    return
end

"""
    Compute the speed of sound in sea water using Coppen's formula at 0 depth
"""
function CoppensFormulaBase(layerData :: LayerData)
    d = layerData.depth / 1000.0     # Depth in kilometers
    s = layerData.salinity / 1000.0  # Salinity in parts per thousand
    t = layerData.temperature / 10.0;  # temperature

    1449.05 + 45.7 * t - 5.21 * t^2 + 0.23 * t^3 + (1.333 - 0.126 * t + 0.009 * t^2) * (s - 35.0)
end

function CompactMackenzieFormula(layerData :: LayerData)
    d = layerData.depth;
    s = layerData.salinity;
    t = layerData.temperature; 
    return 1448.96 + 4.591 * t - 0.05304 * t^2 + 0.0002374 * t^3 + 0.016 * d;
end

function ∂CompMacOnR(layerData ::LayerData)
    return 0;
end

function ∂CompMacOnZ(layerData ::LayerData)
    return 0.016;
end 

function MackenzieFormula(layerData :: LayerData)
    d = layerData.depth;
    s = layerData.salinity;
    t = layerData.temperature; 
    return 1448.96 + 4.591 * t - 0.05304 * t^2 + 0.0002374 * t^3 + 0.016 * d + (1.34 - 0.01025 * t) * (s - 35) + 1.675 * 10^(-7) * d^2 - 7.139 * 10^(-13) * t * d^3;
end

function ∂MacOnR(layerData :: LayerData)
    return 0;
end

function ∂MacOnZ(layerData :: LayerData)
    d = layerData.depth;
    s = layerData.salinity;
    t = layerData.temperature; 

    res = 0.016 + 3.350 * 10^(-7) * d - 7.139 * 10^(-13) * t * 3 * d^2; 
    return res
end
