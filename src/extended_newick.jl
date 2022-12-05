write(filename::AbstractString, arg::ARG; pruned_singletons=true) = write(filename, extended_newick(arg; pruned_singletons))

isshared(n) = (sum(n.color) == length(n.color))
"""
	clade_depth(node::TreeNode)

Topologic distance from `node` to leaves.
"""
function clade_depth(node::TreeNode)
	d = 0
	_node = node
	while !_node.isleaf
		_node = _node.children[1]
		d += _node.tau
	end
	return d
end

function extended_newick(arg::ARG; pruned_singletons=true)
    hybrids = Dict()
    if pruned_singletons
        nwk = extended_newick_pruned(arg)
    else
        nwk = extended_newick(arg.root[end], nothing, hybrids)*";"
    end
    return nwk
end

function color_as_label(color::Array{Bool, 1})
    color_int = Array{Int, 1}(undef, 0)
    for (i, c) in enumerate(color)
        if c
            push!(color_int, i-1)
        end
    end
    name = "&segments={$(color_int)}"
    name = replace(name," " => "", "[" => "", "]" => "")
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

function extended_newick_pruned(arg::ARG)
	hybrids = Dict()

	if arg.root[1] == arg.root[2]
		# (i)
        str = extended_newick_pruned(arg.root[1], nothing, hybrids)
		return str * ";"
	else
		if !isshared(arg.root[1]) && !isshared(arg.root[2])
			#(ii)
			str1 = extended_newick_pruned(arg.root[1], nothing, hybrids)
			str2 = extended_newick_pruned(arg.root[2], nothing, hybrids)
			return "($(str1),$(str2))GlobalRoot:0.;"
		elseif isshared(arg.root[1])
			#(iii)
			str = extended_newick_pruned(arg.root[2], nothing, hybrids)
			return str * ";"
		elseif isshared(arg.root[2])
			str = extended_newick_pruned(arg.root[1], nothing, hybrids)
			return str * ";"
		end
	end
end

function extended_newick_pruned(a_r::ARGNode, a_r_anc::Union{Nothing, ARGNode}, hybrids::Dict)
	nwk = ""
	is_hybrid = (length(filter(x -> !isnothing(x), a_r.anc))>1)
	if !a_r.isleaf && (!is_hybrid  || (is_hybrid && !haskey(hybrids, a_r.label)))
		if (is_hybrid && !haskey(hybrids, a_r.label))
			nwk *= "("
		end
		nwk *= "("
		for c in a_r.children
			nwk *= extended_newick_pruned(c, a_r, hybrids)
			nwk *= ","
		end
		nwk = nwk[1:end-1] # Removing trailing ','
		nwk *= ")"
	end
	# Adding a_r to hybrid if necessary
	if is_hybrid
        if !haskey(hybrids, a_r.label)
            i = length(hybrids) +1
            color = Array{Bool, 1}(undef, 2)
            for c in 1:length(a_r.anccolor[end])
                color[c] = any(map(l->l[c],a_r.anccolor))
            end
            color_label = color_as_label(color)
            hybrids[a_r.label] = "hybrid_node_$(i)#H$(i)"
			if a_r.isleaf
				nwk *= "("
			end
            a_r_label = "$(a_r.label)$(color_label))" *hybrids[a_r.label]
        else
            a_r_label = hybrids[a_r.label]
        end
	else
        # Writing a_r itself
        a_r_label = a_r.label
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