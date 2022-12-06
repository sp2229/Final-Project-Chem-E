function expa(stoichiometric_matrix::Array{Float64,2})

    # initialize -
	(ℳ,ℛ) = size(stoichiometric_matrix)
    T = [Matrix{Float64}(I,ℛ,ℛ) transpose(stoichiometric_matrix)];
	
    # outer loop -
    for col_index = 1:ℳ

        # adjust the col index (from the org matrix -)
	    local_col_index = col_index + ℛ;
        
        # get the current RHS -
        RHS = T[:,(ℛ + 1):end] # this gets the cols belonging to the stoichiometric matrix -

        # Pick a col -
	    TMP_COL = RHS[:,col_index]; # this is grabbing a metabolite -

        # how many zeros, + and - in this col -
        idx_zero_array = findall(x->x==0.0,TMP_COL)
        idx_positive_array = findall(x->x>0.0,TMP_COL)
        idx_negative_array = findall(x->x<0.0,TMP_COL)

        # make new T -
        TNEW = T[idx_zero_array,:]

        # zero out non-zero elements -
        TMP_ARR = Array{Array{Float64,1},1}()
        for idx_positive ∈ idx_positive_array
            for idx_negative ∈ idx_negative_array
            
                # setup α and β -
			    α = abs(T[idx_positive,local_col_index]);
			    β = abs(T[idx_negative,local_col_index]); 

                # compute a tmp row -
                TMP_ROW = α*T[idx_negative,:] .+ β*T[idx_positive,:];

                # grab the tmp row -
                push!(TMP_ARR,TMP_ROW)
            end
        end
    
		# ok, so we have a new set of possible pathways, check for independent pathways -
        P = transpose(hcat(TMP_ARR...))		
		if (isempty(P) == false)
			# Update T -
			T = vcat(TNEW,P);
		end
    end # main -

	# ok, check for indepedent pathways -
	rows_to_keep_set = compute_convex_independent_rows(T)

	# grab: which rows should we keep?
    T_IND = T[(rows_to_keep_set |> collect),:]
	
    # return - 
    return T_IND

end # function -

function compute_convex_independent_rows(P)::Set{Int64}

	# ok, so we have a new set of possible pathways, check for independent pathways -
	(NRP,_) = size(P)
	IPA = Matrix{Float64}(I,NRP,NRP)
	for i ∈ 1:NRP

		# get rᵢ -
		rᵢ = P[i,:]

		# compare rᵢ against rⱼ i neq j
		for j ∈ 1:NRP
			if (i != j)

				# get rⱼ -
				rⱼ = P[j,:]
				
				# check: is convex independent?
				if (is_convex_independent(rᵢ,rⱼ) == true)
					IPA[i,j] = 1.0					
				else
					IPA[i,j] = 0.0	
				end
			end
		end # inner 
	end # outer

	# rows we should keep -
	rows_to_keep_set = Set{Int64}()
	for i ∈ 1:NRP
		sᵢ = sum(IPA[i,:])
		if (sᵢ == NRP)
			push!(rows_to_keep_set,i)
		end
	end

	# return -
	return rows_to_keep_set
end

function is_convex_independent(rᵢ,rⱼ)::Bool
	idx_zero_rᵢ = findall(x->x==0.0,rᵢ)
	idx_zero_rⱼ = findall(x->x==0.0,rⱼ)	
	return issubset(idx_zero_rᵢ,idx_zero_rⱼ) ? false : true
end