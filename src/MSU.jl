function build(model::Type{MSULatticeModel}; ṅₒ::Float64, L::Int64, u::Float64, d::Float64)::MSULatticeModel

    # build an empty lattice model -
    lattice_model = MSULatticeModel();
    
    # alias some stuff -
    number_of_levels = L;

    # compute connectivity - 
    number_items_per_level = [i for i = 1:number_of_levels]
    tmp_array = Array{Int64,1}()
    theta = 0
    for value in number_items_per_level
        for _ = 1:value
            push!(tmp_array, theta)
        end
        theta = theta + 1
    end

    N = sum(number_items_per_level[1:(number_of_levels-1)])
    connectivity_index_array = Array{Int64,2}(undef, N, 3)
    for row_index = 1:N

        # index_array[row_index,1] = tmp_array[row_index]
        connectivity_index_array[row_index, 1] = row_index
        connectivity_index_array[row_index, 2] = row_index + 1 + tmp_array[row_index]
        connectivity_index_array[row_index, 3] = row_index + 2 + tmp_array[row_index]
    end
    lattice_model.connectivity = connectivity_index_array

    # init the tree -
    total_number_of_lattice_nodes = connectivity_index_array[end, end]
    number_of_nodes_to_evaluate = connectivity_index_array[end, 1]
    tree_value_array = Array{Float64,1}(undef, total_number_of_lattice_nodes) # nodes x 3 = col1: underlying price, col2: intrinsic value, col3: option price

    # First: let's compute the underlying price on the lattice -
    tree_value_array[1, 1] = ṅₒ
    for node_index ∈ 1:number_of_nodes_to_evaluate

        # get index -
        parent_node_index = connectivity_index_array[node_index, 1]
        up_node_index = connectivity_index_array[node_index, 2]
        down_node_index = connectivity_index_array[node_index, 3]

        # compute prices -
        parent_price = tree_value_array[parent_node_index, 1]
        up_price = parent_price * u
        down_price = parent_price * d

        # store prices -
        tree_value_array[up_node_index, 1] = up_price
        tree_value_array[down_node_index, 1] = down_price
    end

    # Second: let's add the data to the lattice model -
    lattice_model.data = tree_value_array

    # return the model -
    return lattice_model
end

function build_children_dictionary(nodes::Dict{Int64,Array{Int64,1}})::Dict{Int64,Array{Int64,1}}

    # initialize -
    children_index_dictionary = Dict{Int64, Array{Int64,1}}();
    
    # build the kids for the root node -
    children_index_dictionary[1] = nodes[1];

    # get the keys for nodes dictionary -
    number_of_tree_levels = length(keys(nodes));
    for l ∈ 1:(number_of_tree_levels-1)
        
        # get the nodes on this level -
        level_node_array = nodes[l];

        # what is the next level index?
        next_level_index = l+1;
        for i ∈ 1:next_level_index

            node_index = level_node_array[i];

            # compute the descendants -
            # initialize new darray -
            darray = Array{Int64,1}(undef,2)
            darray[1] = descendant(node=node_index, nextlevel=next_level_index, branch=0);
            darray[2] = descendant(node=node_index, nextlevel=next_level_index, branch=1);
            children_index_dictionary[node_index] = darray;
        end
    end

    # return -
    return children_index_dictionary;
end

function build_nodes_dictionary(levels::Int64)::Dict{Int64,Array{Int64,1}}

    # initialize -
    index_dict = Dict{Int64, Array{Int64,1}}()

    counter = 0
    for l = 0:levels
        
        # create index set for this level -
        index_array = Array{Int64,1}()
        for _ = 1:(l+1)
            counter = counter + 1
            push!(index_array,counter)
        end

        index_dict[l] = index_array
    end

    # return -
    return index_dict
end

function descendant(; node::Int64=0, nextlevel::Int64=1, branch::Int64=0)
    return (node+nextlevel+branch)
end

function build_probability_dictionary(model::MSULatticeModel, levels::Int64)::Dict{Int64, Array{Float64,1}}

    # initialize -
    probability_dict = Dict{Int64, Array{Float64,1}}()
    p = model.p

    for l = 0:levels
        
        # initialize -
        probability_array = Array{Float64,1}()
        
        # generate k range -
        karray = range(0,step=1,stop=l) |> collect

        for k ∈ karray
            tmp = binomial(l,k)*p^(k)*(1-p)^(l-k)
            push!(probability_array,tmp)
        end

        # grab the array - note: we have to reverse (d move is first, we need the other way arround)
        probability_dict[l] = reverse(probability_array)
    end

    # return -
    return probability_dict
end
