# function read_arg(str::AbstractString)

# end

function parse_extended_newick(str::AbstractString; degree=2)
	# S = (str[end] == ';') ? str[1:end-1] : str
	# strict=false below to avoid warnings for nodes with one child only.
	tree = parse_newick_string(str; force_new_labels=true)
	node_data = Dict(n => parse_node_label(n) for n in nodes(tree))
	arg = arg_from_ext_tree(tree, node_data)
	return arg
end


function arg_from_ext_tree(tree, node_data; degree=2)
	arg = ARG(;degree)
	add_node_to_arg!(arg, tree.root, nothing, node_data, 1; degree)
	return arg
end

function add_node_to_arg!(arg, n::TreeNode, a, node_data, i; degree=2, hybrids=[])
	# Color of node
	color = zeros(Bool, degree)
	for c in node_data[n][2]["segments"]
		color[c+1] = true
	end

	# If n is not a reassorted node
	an = if isempty(node_data[n][2]["reassortment"])
		label = string(node_data[n][1])
		an = ARGNode(;
			degree,
			anc = [a],
			anccolor = [color],
			color = color,
			label = isempty(label) ? "ARG_NODE_$i" : label,
			tau = Union{Missing, Float64}[n.tau],
			isroot = n.isroot ? color : zeros(Bool, degree),
			isleaf = n.isleaf
		)
		arg.nodes[an.label] = an
		an.isleaf && (arg.leaves[an.label] = an)
		if n.isroot
			for i in findall(color)
				arg.root[i] =  an
			end
		end

		an

	# If n is a reassorted node that is not already in the ARG
	# Then we add it as a normal ARG node, but using the reassortment name as a label
	elseif !haskey(arg.nodes, string(split(string(node_data[n][1]),"#")[1]))
		#label = string(node_data[n][2]["reassortment"]) #this is where the name is removed
		label = string(split(string(node_data[n][1]),"#")[1])
		push!(hybrids, label)
		#label = string(node_data[n][1])
		an = ARGNode(;
			degree,
			anc = [a],
			anccolor = [color],
			color = color,
			label = isempty(label) ? "ARG_NODE_$i" : label,
			tau = Union{Missing, Float64}[n.tau],
			isroot = n.isroot ? color : zeros(Bool, degree),
			isleaf = n.isleaf
		)
		arg.nodes[an.label] = an
		an.isleaf && (arg.leaves[an.label] = an)
		if n.isroot
			for i in findall(color)
				arg.root[i] =  an
			end
		end

		an

	# Else, add n to an already existing reassorted node
	# In this case, color represents the color to the ancestor.
	# --> it must be added to `an.color` and pushed to `an.anccolor`
	else
		label = string(split(string(node_data[n][1]),"#")[1])
		an = arg.nodes[label]
		push!(an.anc, a)
		push!(an.anccolor, color)
		push!(an.tau, n.tau)
		an.color += color
		n.isroot && (an.isroot += color)

		# If `n` is not a leaf, then `an` should not be one
		# Only one of the set of reassorted nodes will be potentially non-leaf
		if !n.isleaf
			an.isleaf = false
			delete!(arg.leaves, an.label)
		end

		# If n is a root and an was not one (we added color to an.isroot already)
		if n.isroot && an.isroot .== color
			for i in findall(color)
				arg.root[i] =  an
			end
		end
		an
	end

	##check if a is infact a root -> i.e. an internal root
	if !isnothing(a) && !all(a.isroot) && a.label ∉ hybrids
		##this happens if the a has lost a color that was seen in n
		for i in 1:length(color)
			if color[i] && !a.color[i]
				##this is a missing color -> root of this subtree
				arg.root[i] =  a
				a.isroot[i] = true
			end
		end
	end

	# Adding children
	for (k,c) in enumerate(n.child)
		cn = add_node_to_arg!(arg, c, an, node_data, i+k; degree, hybrids)
		if !in(cn, an.children)
			push!(an.children, cn)
		end
	end

	return an
end


parse_node_label(n::TreeNode, delim="__") = parse_node_label(n.label, delim)
function parse_node_label(n::AbstractString, delim="__")
	# Remove randstr that node2tree added
	m_delim = match(Regex(delim), n)
	label = if !isnothing(m_delim)
		n[1:(m_delim.offset-1)]
	else
		n
	end

	# Get annotation at the end of the label
	m_ann = match(r"\[&(.)*\]", label)

	label, annotation = if !isnothing(m_ann)
		(label[1:(m_ann.offset-1)], label[m_ann.offset:end])
	else
		(label, "")
	end
	# Parse annotation into a dict
	annotation = parse_annotation(annotation)

	# Parse the label for '#i' indicating a reassorted node
	m_rea = match(r"#[RH][0-9]*", label)
	if !isnothing(m_rea)
		haskey(annotation, "reassortment") && error("Node $n has 'reassortment' in its annotation.")
		annotation["reassortment"] = m_rea.match
		# label = label[1:(m_rea.offset-1)]
	else
		annotation["reassortment"] = ""
	end
	return label, annotation
end

function parse_annotation(str::AbstractString)
	# str with format "[&k1=v1,k2=v2,...]"
	# !!!!! For now, format "[&key={i1,i2,...}]"
	if length(str) < 2 || str[1:2] != "[&" || str[end] != ']'
		return Dict()
	end
	S = str[3:end-1]
	keys, vals = get_key_vals(S)
	if length(keys) != length(vals)
		error("Error in parsing annotation $(str): different number of keys and values.")
	end
	out = Dict()
	for (key, val) in zip(keys, vals)
		pval = if key == "segments"
			valstr = split(val[2:end-1], ',') # Skipping the '{' and '}'
			if valstr == SubString{String}[""]
				[]
			else
				[parse(Int, s) for s in valstr]
			end
		else
			val
		end
		out[key] = pval
	end

	return out
end

# Find "=" signs that are not between curly brackets
function get_key_vals(S)
	bra = 0
	i0 = 1
	keys = []
	vals = []
	for (i,x) in enumerate(S)
		if x == '{'
			bra += 1
		elseif x == '}'
			bra -= 1
		end
		if bra < 0
			error("Inconsistent curly brackets $S")
		end
		#
		if bra == 0 && x == '='
			push!(keys, S[i0:i-1])
			i0 = i+1
		end
		#
		if bra == 0 && (x == ',' || i == length(S))
			x == ',' ? push!(vals, S[i0:i-1]) : push!(vals, S[i0:i])
			i0 = i+1
		end
	end
	return keys, vals
end

write(filename::AbstractString, arg::ARG; pruned_singletons=true, tau=false) = write(filename, extended_newick(arg; pruned_singletons, tau))

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

function extended_newick(arg::ARG; pruned_singletons=true, tau=false)
	hybrids = Dict()
	strs = Array{String,1}(undef, 0)
	roots = Array{ARGNode}(undef, 0)
	for i in 1:length(arg.root)
		r = arg.root[i]
		if !isnothing(r)
			push!(roots, r)
		end
	end
	roots = unique(roots)
	
	#sort roots, the oldest root is the most exterior root and needs to root the ARG
	if roots[1].isleaf || ismissing(roots[1].children[1].tau[1])
		sort!(roots, by=clade_depth_no_tau, rev=true)
	else
		sort!(roots, by=clade_depth, rev=true)
	end
	## keep track of colors that have been added -> if roots[2] has a color that was 
	## seen in roots[1] it means that roots[2] was added to the nwk with the subtree of roots[1]
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
			str = extended_newick(r, nothing, hybrids; pruned_singletons, tau)
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
		nwk *= ")GlobalRoot"*color_as_label(color_seen)*":0."
	else
		nwk= strs[1]
	end
	return nwk * ";"
end

function extended_newick(a_r::ARGNode, a_r_anc::Union{Nothing, ARGNode}, hybrids::Dict; pruned_singletons=true, tau=false)
	nwk = ""
	is_hybrid = (length(filter(x -> !isnothing(x), a_r.anc))>1)
	if !a_r.isleaf && (!is_hybrid  || (is_hybrid && !haskey(hybrids, a_r.label)))
		if pruned_singletons && (is_hybrid && !haskey(hybrids, a_r.label))
			nwk *= "("
		end
		nwk *= "("
		for c in a_r.children
			nwk *= extended_newick(c, a_r, hybrids; pruned_singletons, tau)
			nwk *= ","
		end
		nwk = nwk[1:end-1] # Removing trailing ','
		nwk *= ")"
	end
	if is_hybrid
		#add to hybrid, and name if not yet seen
        if !haskey(hybrids, a_r.label)
			i = length(hybrids) +1
			if pruned_singletons
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
				hybrids[a_r.label] = "#H$(i)"
				a_r_label =  a_r.label * hybrids[a_r.label]
			end
        else
			if pruned_singletons
            	a_r_label = hybrids[a_r.label]
			else
				a_r_label = a_r.label * hybrids[a_r.label]
			end
        end
	else
        # Writing a_r itself
        a_r_label = a_r.label
	end

	#write node information
	nwk *= a_r_label
    color = a_r.color
    for (i, a) in enumerate(a_r.anc)
        if a == a_r_anc
            color = a_r.anccolor[i]
        end
    end
	nwk *= color_as_label(color)
	if tau && !ismissing(a_r.tau[1])
		nwk *= ":"*string(a_r.tau[1])
	end

	return nwk
end


