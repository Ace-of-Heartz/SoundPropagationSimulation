struct RKData
    comatrix :: Matrix{<:AbstractFloat}
    weightVector :: Vector{<:AbstractFloat}
end


function GetRK45Codata() :: RKData
    comatrix :: Matrix{<:AbstractFloat} = [
        0.0 0.0 0.0 0.0; 
        1.0/2.0 1.0/2.0 0.0 0.0;
        1.0/2.0 0.0 1.0/2.0 0.0;
        1.0 0.0 0.0 1.0
    ];

    weightVector :: Vector{<:AbstractFloat} = [ 1.0/6.0 , 1.0/3.0 , 1.0/3.0 , 1.0/6.0];


    return RKData(comatrix,weightVector);
end 

