import Base.+

function +(buffer::Array{String,1}, line::String)
    push!(buffer, line)
end

function read_reaction_file(path_to_file::String)::Array{String,1}

    # initialize -
    vff_file_buffer = String[]
    vff_reaction_array = Array{String,1}()

    # Read in the file -
    open("$(path_to_file)", "r") do file
        for line in eachline(file)
            +(vff_file_buffer,line)
        end
    end

    # process -
    for reaction_line ∈ vff_file_buffer
        
        # skip comments and empty lines -
        if (occursin("//", reaction_line) == false && 
            isempty(reaction_line) == false)
        
            # grab -
            push!(vff_reaction_array,reaction_line)
        end
    end

    # return -
    return vff_reaction_array
end