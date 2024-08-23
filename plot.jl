include("sonarForms.jl");

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


"""
    PreparePlotData(rays :: Vector{RayData},timeVals :: Vector{Float32}) :: PlotData
    speedVals    = map(x -> (x.speed), rays);

    Prepares the interesting data for plotting
"""
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



"""
    PlotData(plotDatas)

    Plots the interesting data
"""
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


    Axis(fig[1,1], title = "Ray Paths", yreversed = true);

    # for i in eachindex(labels)
    #     for j in eachindex(plotDatas)
    #         Axis(fig[2,1][j, i], title = titles[i],
    #             xlabel = labels[i][1],
    #             ylabel = labels[i][2]
    #         );
    #     end
    # end

    
    for i in eachindex(plotDatas)
        plotData = plotDatas[i];
        
        # lines!(
        #     fig[2,1][i,2],
        #     plotData.timeVals,
        #     plotData.speedVals,
        #     );
        # lines!(
        #     fig[2,1][i,1],
        #     map(x -> x[2],plotData.positionVals),
        #     plotData.speedVals,
        # );
        
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