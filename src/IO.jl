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

function add_node_to_arg!(arg, n::TreeNode, a, node_data, i; degree=2)
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
	elseif !haskey(arg.nodes, node_data[n][2]["reassortment"])
		label = string(node_data[n][2]["reassortment"])
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
		label = string(node_data[n][2]["reassortment"])
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

	# Adding children
	for (k,c) in enumerate(n.child)
		cn = add_node_to_arg!(arg, c, an, node_data, i+k; degree)
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
			[parse(Int, s) for s in valstr]
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


function write_extnewick(arg::ARG)
	g = add_internal_above_reassorted_leaves(arg)
	reassorted_nodes = get_reassorted_nodes(arg)
end

# Array of nodes with more than one parent
function get_reassorted_nodes(arg)
	reassorted_nodes = ARGNode[]
	for (label, n) in arg.nodes
		if length(n.anc) > 1
			push!(reassorted_nodes, n.label)
		end
	end

	return reassorted_nodes
end

# Add singletons internal nodes above reassorted leaves
# Useful for the convention that one segment goes to reassorted nodes with leaves, and the
# other to the ones without leaves
function add_internal_above_reassorted_leaves(arg::ARG)

end
