write(filename::AbstractString, arg::ARG; pruned_singletons=true) = write(filename, extended_newick(arg; pruned_singletons))

isshared(n) = (sum(n.color) == length(n.color))
"""
	clade_depth(node::TreeNode)

Topologic distance from `node` to leaves.
"""
function clade_depth(node::ARGNode)
	d = 0
	_node = node
	while !_node.isleaf
		_node = _node.children[1]
		d += _node.tau[1]
	end
	return d
end

function clade_depth_no_tau(node::ARGNode)
	d = 0
	if !node.isleaf
		d = max([clade_depth_no_tau(c) for c in node.children]) + 1
	else
		d = 1
	end
	return d
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

function extended_newick(arg::ARG; pruned_singletons=true)
	hybrids = Dict()
	strs = Array{String,1}(undef, 0)
	roots = Array{ARGNode}(undef, 0)
	for i in 1:length(arg.root)
		try
			r = arg.root[i]
			if !isnothing(r)
				push!(roots, r)
			end
		catch e
			println("There is sth wrong with the ARG")
		end
	end
	roots = unique(roots)
	#roots = unique(filter(x -> !isnothing(x), arg.root))
	if ismissing(roots[1].tau[1])
		sort!(roots, by=clade_depth_no_tau, rev=true)
	else
		sort!(roots, by=clade_depth, rev=true)
	end
	color_seen = zeros(length(roots[1].color))

	for r in roots
		seen = false
		seen_new = color_seen
		for (i,(c,c_seen)) in enumerate(zip(r.color, color_seen))
			if c && c==c_seen
				seen = true
			end
			if c
				seen_new[i] = 1
			end
		end
		if !seen
			if pruned_singletons
				str = extended_newick_pruned(r, nothing, hybrids)
			else
				str = extended_newick(r, nothing, hybrids)
			end
        	push!(strs, str)
			color_seen = seen_new
		end
	end
	nwk = ""
	if length(strs) >1
		nwk *= "("
		for str in strs
			nwk *= str
			nwk *= ","
		end
		nwk = nwk[1:end-1] # Removing trailing ','
		nwk *= ")GlobalRoot:0."
	else
		nwk= strs[1]
	end
	return nwk * ";"
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
            color = Array{Bool, 1}(undef, length(a_r.color))
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