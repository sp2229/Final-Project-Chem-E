mutable struct MSULatticeModel

    # data -
    connectivity::Array{Int64,2}
    data::Array{Float64,1}

    # constructor -
    MSULatticeModel() = new();
end