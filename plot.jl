
using GLMakie

function PlotSoundSpeedProfile(
    ds :: AbstractArray{<:AbstractFloat}, 
    ts :: AbstractArray{<:AbstractFloat}, 
    ss :: AbstractArray{<:AbstractFloat},
    speedFunc :: Function
    )

    if (length(ds) != length(ts) || length(ds) != length(ss))
        error("Number of values for depth, temperature and salinity must be equal.");
    end

    cs = map((d,t,s) -> speedFunc(d,t,s),ds,ts,ss)
    ps = map((d,s,c) -> (d,s,c), ds,ss,cs);

    fig = Figure()
    ax = Axis3(fig[1,1],xlabel = "Depth", ylabel = "Salinity",zlabel = "Velocity")
    lines!(
        ax,
        ps[2:length(ps)],
        color = ts[2:length(ps)],
        linewidth = 5.0,
        joinstyle = :round
    )    
    display(fig)
end

function PlotSoundSpeedProfile(
    ds :: LinRange{<:AbstractFloat,<:Integer}, 
    ts :: LinRange{<:AbstractFloat,<:Integer}, 
    ss :: LinRange{<:AbstractFloat,<:Integer},
    speedFunc :: Function
)
    dArray = collect(ds);
    tArray = collect(ts);
    sArray = collect(ss);

    PlotSoundSpeedProfile(dArray,tArray,sArray,speedFunc);
end

function PlotRayAttributes((ss,rs,zs,ξs,ζs) :: Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat},Vector{<:AbstractFloat},Vector{<:AbstractFloat},Vector{<:AbstractFloat}}) 
    fig = Figure();

    Axis(fig[1,1], title = "Ray ξ Attribute", xlabel = "Arclength (m)", ylabel = "ξ");
    Axis(fig[1,2], title = "Ray ζ Attribute", xlabel = "Arclength (m)", ylabel = "ζ");
    Axis(fig[2,1], title = "Ray Propagation", xlabel = "Horizontal Position (m)", ylabel = "Depth (m)",yreversed = true);
    
    PlotRayLines((ss,rs,zs,ξs,ζs),fig);
    display(fig)
end

function PlotRayLines(
    (ss,rs,zs,ξs,ζs) :: Tuple{Vector{<:AbstractFloat},Vector{<:AbstractFloat},Vector{<:AbstractFloat},Vector{<:AbstractFloat},Vector{<:AbstractFloat}},
    fig :: Figure
)
    lines!(fig[1,1],ss, ξs);
    lines!(fig[1,2],ss, ζs);
    lines!(fig[2,1],rs,zs);
end

function PlotMultipleRayAttributes(res)
    fig = Figure();
    
    Axis(fig[1,1], title = "Ray ξ Attribute", xlabel = "Arclength (m)", ylabel = "ξ");
    Axis(fig[1,2], title = "Ray ζ Attribute", xlabel = "Arclength (m)", ylabel = "ζ");
    Axis(fig[2,1], title = "Ray Propagation", xlabel = "Horizontal Position (m)", ylabel = "Depth (m)", yreversed = true);
    
    for r in res
        PlotRayLines(r,fig);
    end 
    display(fig)
end


function PlotBathymetry()
