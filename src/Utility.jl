function generate_reaction_string(reaction_string::String)

    # split the reaction into its components -
	component_array = split(reaction_string,',');

    # reaction name -
	reaction_name = String.(component_array[1]);
	reactant_phrase = String.(component_array[2]);
	product_phrase = String.(component_array[3]);

    # reaction string 
    return "$(reactant_phrase) -> $(product_phrase)"
end

function generate_reaction_strings(reactions::Array{String,1}; expand::Bool = false)::Array{String,1}

    # initialize -
    reaction_string_array = Array{String,1}()

    # should we expand the reversible reactions?
    reactions_to_process = reactions;
    if (expand == true)
        reactions_to_process = expand_reversible_reactions(reactions);
    end

    # first: let's discover the species list -
	for reaction_string ∈ reactions_to_process
        
        # rstring -
        rstring = generate_reaction_string(reaction_string)
        push!(reaction_string_array, rstring)
    end

    # return -
    return reaction_string_array
end