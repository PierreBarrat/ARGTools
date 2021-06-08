"""
	MCCs_from_arg(arg)

For each arg leave, find the set of leaves that are connected to it by full branches.
"""
function MCCs_from_arg(arg)
	mccs = Dict{String, Array{String,1}}()
	visited = Dict(x => false for x in values(arg.nodes))
	for (s1,n1) in arg.leaves
		if !visited[n1]
			mccs[s1] = String[s1]
			visited[n1] = true
			for (s2,n2) in arg.leaves
				if !visited[n2] && is_linked_pair(n1,n2)
					push!(mccs[s1], s2)
					visited[n2] = true
				end
			end
		end
	end
	# return sort(mccs, by=x->length(x[2]))
	return sort([sort(m) for m in values(mccs)], lt=clt)
end

function MCCs_pairs_from_arg(arg)
	MCCs = Dict()
	for i in 1:arg.degree, j in (i+1):arg.degree
		MCCs[i,j] = _MCCs_from_arg(arg, i, j)
		MCCs[j,i] = MCCs[i,j]
	end
	return MCCs
end
function _MCCs_from_arg(arg, i::Vararg{Int})
	mccs = Array{Array{String,1}}(undef, 0)
	visited = Dict(x=>false for x in keys(arg.leaves))
	for (s1, n1) in Iterators.filter(x->!visited[x[1]], arg.leaves)
		a = n1
		k = find_common_branch(a, i...)
		while !isnothing(k)
			a = a.anc[k]
			k = find_common_branch(a, i...)
		end
		# n1 is in the MCC with root `a`
		push!(mccs, mcc_from_root(a, i...))
		for x in mccs[end]
			visited[x] = true
		end
	end
	return sort(mccs; lt=clt)
end

function mcc_from_root(r, i...)
	mcc = String[]
	mcc_from_root!(mcc, r, i...)
	return sort(mcc)
end
function mcc_from_root!(mcc, r, i...)
	for c in r.children
		if is_common_branch(r, c, i...)
			mcc_from_root!(mcc, c, i...)
		end
	end
	r.isleaf && push!(mcc, r.label)
	return nothing
end

"""
	is_common_branch(n[, i::Integer, j::Vararg{Integer}])

Check if one of the branches above `n` is common to colors `[i, j...]`, or to all trees
  if `i` and `j` are not provided.
"""
function is_common_branch(n)
	return mapreduce(x->!x, *, n.isroot; init=true) &&
		n.degree == length(n.color) &&
		length(n.anc) == 1
end
function is_common_branch(n, i::Integer, j::Vararg{Integer})
	if n.isroot[i] || !mapreduce(x->!n.isroot[x], *, j; init=true)
		return false
	else
		for ab in n.anccolor
			if ab[i] && mapreduce(x->ab[x], *, j; init=true)
				return true
			end
		end
		return false
	end
end
"""
	is_common_branch(a, n, i::Integer, j::Vararg{Integer})

Check if the branch going from `n` to `a` is common to colors `[i, j...]`.
	Return `nothing` if `a` is not a direct ancestor of `n`.
"""
function is_common_branch(a, n, i::Integer, j::Vararg{Integer})
	k = find_ancestor_index(a, n)
	if isnothing(k)
		return false
	else
		return n.anccolor[k][i] && mapreduce(x->n.anccolor[k][x], *, j; init=true)
	end
end
"""
	find_common_branch(n, i::Integer, j::Vararg{Integer})

If `n` has an ancestor branch common to colors `i` and `j`, return its index.
  Otherwise, return nothing.
"""
function find_common_branch(n, i::Integer, j::Vararg{Integer})
	if n.isroot[i] || !mapreduce(x->!n.isroot[x], *, j; init=true)
		return nothing
	else
		for (k,ab) in enumerate(n.anccolor)
			if ab[i] && mapreduce(x->ab[x], *, j; init=true)
				return k
			end
		end
		return nothing
	end
end


function clt(x,y)
    if length(x) < length(y)
        return true
    elseif length(x) > length(y)
        return false
    else
        return x[1] < y[1]
    end
end


"""
	coherent_subtrees(r::ARGNode)
"""
function coherent_subtrees(r::ARGNode)
	for c in r.color
		if !c
			@error "Node should exist for all segments."
		end
	end
	#
	if r.isleaf
		out = [r.label]
	else
		out = String[]
		for c in r.children
			if length(c.anc) == 1 && sum(c.color) == length(c.color)
				append!(out, coherent_subtrees(c))
			end
		end
	end
	return out
end

function coherent_subtrees(arg)
	out = Array{String,1}[]
	for n in values(arg.nodes)
		if sum(n.color) == length(n.color) && length(n.anc) > 1
			sb = coherent_subtrees(n)
			!isempty(sb) && push!(out, coherent_subtrees(n))
		end
	end
	return out
end

"""
Goal: Can I link the two nodes with fully common branches?
Method:
- go up from `n1` until a non-common branch is found. Store nodes found.
- go up from `n2` until a non-common branch is found (false) or a previously stored node is found (true).
"""
function is_linked_pair(n1::ARGNode, n2::ARGNode)
	if n1 == n2
		return true
	end
	d = Dict{String,Bool}()
	a = n1
	while is_common_branch(a)
		a = a.anc[1]
		d[a.label] = true
	end
	#
	a = n2
	while is_common_branch(a)
		a = a.anc[1]
		if get(d, a.label, false)
			return true
		end
	end
	return false
end
