# using DifferentialEquations;
import Pkg;
Pkg.add("GLMakie")
using GLMakie;
using LinearAlgebra;

include("./plot.jl");
include("./sonarForms.jl");
include("./rkMethods.jl");




"""
    ShootRays(
    numberOfRays :: Int,
    timeInterval :: Float64,
    maxSteps :: Int,
    initPos :: Vector{Float64} = [0.0,550.0,0.0]
    )

    Shootrays in a directions with uniform distribution 

"""
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

    rays = [];
    plotDatas = [];

    for i = 0:(numberOfRays-1)

        offset = pi;
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

function ShootRaysStraight(    numberOfRays :: Int,
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

    rays = [];
    plotDatas = [];

    for i = 0:(numberOfRays-1)

        offset = pi/2.0
        angle = i * interval - offset;
        
        angle = abs(angle) < 0.00001 ? 0.00001 : angle

        print("{$angle}\n");
        ray = RayData(initPos,angle,startingSpeed);
        simData = SimulationData(ray,timeInterval,maxSteps);
        
        (timeVals,rays) = RK45(simData);
        
        push!(plotDatas,PreparePlotData(rays,timeVals));
    end

    print("Simulation finished...")
    
    PlotData(plotDatas);
end

