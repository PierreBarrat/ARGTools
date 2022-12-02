
write(filename::AbstractString, arg::ARG) = write(filename, extended_newick(arg))

function extended_newick(arg::ARG)
    hybrids = Dict()
    nwk = extended_newick(arg.root[end], nothing, hybrids)
    return nwk*";"
end

function color_as_label(color::Array{Bool, 1})
    color_int = Array{Int, 1}(undef, 0)
    for (i, c) in enumerate(color)
        if c
            push!(color_int, i-1)
        end
    end
    name = "&segments={$(color_int)}"
    name = replace(name," " => "" )
    name = replace(name,"[" => "" )
    name = replace(name,"]" => "" )
    name = "["*name*"]"
    return name
end

function extended_newick(a_r::ARGNode, a_r_anc::Union{Nothing, ARGNode}, hybrids::Dict)
	nwk = ""
	# Dealing with children of a_r
	if (!a_r.isleaf && !(length(a_r.anc)>1)) || ((length(a_r.anc)>1) && !haskey(hybrids, a_r.label))
		nwk *= "("
		for c in a_r.children
			nwk *= extended_newick(c, a_r, hybrids)
			nwk *= ","
		end
		nwk = nwk[1:end-1] # Removing trailing ','
		nwk *= ")"
	end
	# Adding a_r to hybrid if necessary
	if (length(a_r.anc)>1) && !haskey(hybrids, a_r.label)
		i = length(hybrids) +1
		hybrids[a_r.label] = "#H$(i)"
	end
	# Writing a_r itself
	a_r_label = if (length(a_r.anc)>1)
		a_r.label * hybrids[a_r.label]
	else
		a_r.label
	end
	nwk *= a_r_label
    color = a_r.color
    for (i, a) in enumerate(a_r.anc)
        if a == a_r_anc
            color = a_r.anccolor[i]
        end
    end
	nwk *= color_as_label(color)

	return nwk
end