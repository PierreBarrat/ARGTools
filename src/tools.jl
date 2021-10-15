function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Array{Bool,1}, tau=missing)
	for c in findall(color)
		regraft!(n, oldanc, newanc, c, tau)
	end
end
"""
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Array{Bool,1}, tau=missing)
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, tau=missing)

Regraft node `n` from `oldanc` to `newanc` for color `color`. 
## Note
- `oldanc` has to be an ancestor of `n` for color `color`
- `newanc` and `n` have to be of color `color`
## Warning  
Does not handle `Nothing`. Regrafting a root node **or** to `Nothing` will fail. 
"""
function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, tau=missing; w=false)
	i_old = findfirst(x->x==oldanc, n.anc)
	i_child = findfirst(x->x==n, oldanc.children)
	if isnothing(i_old)
		error("Attempting to regraft node from uncorrect ancestor: $(oldanc.label) not ancestor of $(n.label)")
	elseif isnothing(i_child)
		error("Attempting to regraft node from uncorrect ancestor: $(n.label) not child of $(oldanc.label)")
	elseif length(oldanc.children) == 1
		w && @warn "Removing unique child from internal node ($(oldanc.label), $(n.label))"
	end
	if !n.anccolor[i_old][color]
		error("Attempting to regraft node from uncorrect ancestor: branch from $(n.label) to $(oldanc.label) is not of color $color.")
	end
	if !newanc.color[color]
		error("Attempting to regraft node to ancestor of incorrect color.")
	elseif !n.color[color]
		error("Attempting to regraft node for the wrong color.")
	end
	#
	regraft!(n, oldanc, newanc, color, i_old, tau)
end
"""
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, i_old::Int64)

Core function for regrafting. Does not handle errors.
"""
function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, i_old::Int64, tau=missing)
	# Changes for n and newanc
	graft!(n, newanc, color, tau) 
	# Changes for n and oldanc
	cut_branch!(n, i_old, color)
end
"""
	graft!(n::ARGNode, a::ARGNode, clr::Int64)

Graft `n` onto `a` for color `clr`. 
"""
function graft!(n::ARGNode{T}, a::ARGNode{T}, clr::Int64, tau=missing) where T
	# Changes for `a`
	if !isnothing(a) && !is_ancestor(a, n)[1]
		push!(a.children, n)
	end
	# Changes for `n`
	ia = findfirst(x->x==a, n.anc) # Was `a` already an ancestor of `n` for another color? 
	if isnothing(ia) # We have to create a new ancestor for `n`
		push!(n.anc, a)
		push!(n.anccolor, _color(clr, length(n.color)))
		push!(n.tau, tau)
		push!(n.data, T())
	else # Just set the branch of the right color
		n.anccolor[ia][clr] = true
	end
end

"""
	prune!(n::ARGNode; nochildren=true)

Prune `n`. Return `n`. If `nochildren`, `n` should not have children when pruned. 
"""
function prune!(arg::ARG, n::ARGNode)
	if length(n.children)>0
		@error "Trying to prune node `n.label` with children."
	end
	while !isempty(n.anc)
		for c in findall(n.anccolor[1])
			cut_branch!(n, 1, c)
		end
	end
	delete!(arg.nodes, n.label)
	n.isleaf && delete!(arg.leaves, n.label)
	return n
end

"""
Cut branch from `n` to `n.anc[i]` for color `clr`. Color of `n.anc[i]` is changed if necessary. 1
# Note
`n` is not regrafted and may not have any ancestor left after this.   
`n.anc` may have no child after this. 
"""
function cut_branch!(n::ARGNode, i::Int64, clr::Int64)
	# Checks
	length(n.anc) < i && @error "Cannot delete anc #$i from $(n.label): too few ancestors"
	# Delete branches
	a = n.anc[i]
	n.anccolor[i][clr] = false
	if !(|)(n.anccolor[i]...)
		deleteat!(n.anc, i)
		deleteat!(n.tau, i)
		deleteat!(n.anccolor, i)
		deleteat!(n.data, i)
		if !isnothing(a)
			ic = findfirst(x->x.label==n.label, a.children)
			deleteat!(a.children, ic)
		end
	else 
	end
	# Correct colors if necessary
	correct_color!(a)
end
function cut_branch!(n::ARGNode, idx::Array{Int64}, clr::Int64)
	length(n.anc) < max(idx...) && @error "Cannot delete anc #$(max(idx...)) from $(n.label): too few ancestors"
	# 
	todel = Int64[]
	alist = []
	for i in idx
		a = n.anc[i]
		push!(alist, a)
		n.anccolor[i][clr] = false
		if !(|)(n.anccolor[i]...)
			push!(todel, i)
			if !isnothing(a)
				ic = findfirst(x->x.label==n.label, a.children)
				deleteat!(a.children, ic)
			end
		end
	end
	deleteat!(n.anc, todel)
	deleteat!(n.anccolor, todel)
	for a in alist
		correct_color!(a)
	end
end

###
###


function correct_color!(n::ARGNode)
	if !n.isleaf
		x = zeros(Bool, length(n.color))
		for c in n.children
			i = findfirst(x->x==n, c.anc)
			x .= (|).(x, c.anccolor[i])
		end
		n.color = x
	end
end
correct_color!(n::Nothing) = nothing

"""
	has_singletons(arg)

Singletons are nodes `n` such that `length(n.children) == 1`.
"""
function has_singletons(arg)
	for n in values(arg.nodes)
		if length(n.children) == 1
			return true
		end
	end
	return false
end
"""
	prune_singletons!(arg::ARG; v=false)

Singletons are nodes `n` such that `length(n.children) == 1`.  
"""
function prune_singletons!(arg::ARG; v=false, Nit=1e3)
	npruned = 0
	nit = 0
	while has_singletons(arg) && nit < Nit
		for n in values(arg.nodes)
			if length(n.children) == 1
				# Graft the child onto ancestors
				v && println("Prune $(n.label).\n Ancestors $([x.label for x in n.anc]).\n Child $(n.children[1].label).")
				for (i,a) in enumerate(n.anc)
					ic = findfirst(==(n), n.children[1].anc) # Index of n in n.children[1].anc
					clr = n.anccolor[i]
					regraft!(n.children[1], n, a, clr, n.children[1].tau[ic] + n.tau[i])
				end
				prune!(arg, n)
				npruned += 1
			end
		end
		if length(arg.nodes) == 1
			@warn "ARG with one node only."
			break
		end
		nit += 1
	end
	# check_arg(arg)
	v && println("Pruned $npruned nodes")
	return nothing
end

"""
	has_lone_nodes(arg::ARG)

A lone node `n` is such that `!n.isleaf && length(n.children)==0`. This can arise when simulating an ARG.
"""
function has_lone_nodes(arg::ARG)
	for n in values(arg.nodes)
		if !n.isleaf && length(n.children) == 0
			return true
		end
	end
	return false
end
"""
	prune_lone_nodes!(arg::ARG; v=false)

A lone node `n` is such that `!n.isleaf && length(n.children)==0`.
This can happen when simulating the ARG.
*Note to self*: how can this happen?!
"""
function prune_lone_nodes!(arg::ARG; v=false, Nit=1e3)
	npruned = 1
	nit = 0
	while has_lone_nodes(arg) && nit < Nit
		for n in values(arg.nodes)
			if !n.isleaf && length(n.children) == 0
				prune!(arg, n)
				npruned += 1
			end
		end
		if length(arg.nodes) == 1
			@warn "ARG with one node only."
			break
		end
		nit += 1
	end
	v && println("Pruned $npruned nodes")
	return nothing
end


"""
	find_ancestor_index(a, c)

Find index if `a` in the ancestors of `c`.
"""
find_ancestor_index(a, c) = findfirst(x -> x==a, c.anc)
"""
	is_ancestor(a::ARGNode, c::ARGNode, color::Vararg{Int64})

Check if `a` is an ancestor of `c` for colors in `color`. Also check whether `c` is in `a.children`. 		
"""
function is_ancestor(a::ARGNode, c::ARGNode, color::Vararg{Integer})

	i_anc = findfirst(x->!isnothing(x) && x.label==a.label, c.anc)
	flag = !isnothing(i_anc) && !isnothing(findfirst(x->x.label==c.label, a.children)) 
	if flag
		for clr in color
			flag *= c.anccolor[i_anc][clr]
		end
	end
	return flag
end
"""
	is_ancestor(a::ARGNode, c::ARGNode)

Check if `a` is an ancestor of `c`. Return a `Bool` as well as an array of colors for which `a` is an ancestor of `c`. 		
"""
function is_ancestor(a::ARGNode, c::ARGNode)
	# println()
	# println(a.label)
	# println(c.label)
	out = zeros(Bool, length(c.color))
	flag = false
	for clr in findall(c.color)
		if is_ancestor(a,c,clr)
			try
				out[clr] = true
			catch
				println(a.label)
				println(c.label)
			end
			flag = true
		end
	end
	return flag, out
end
is_ancestor(a::Nothing, c::ARGNode) = ((|)(c.isroot...), c.isroot)
is_ancestor(a::Nothing, c::Nothing) = (false, zeros(Bool, 0))
is_ancestor(a::ARGNode, c::Nothing) = (false, zeros(Bool, 0))

"""
	get_children_index(a::ARGNode, clr::Int64)

Get index of children of `a` for color `clr`. 
"""
function get_children_index(n::ARGNode, clr::Int)
	idx = Int64[]
	if !n.color[clr]
		error("$(n.label) is not of color $clr")
	end
	for (ic,c) in enumerate(n.children)
		# Does `c` has ancestor `n` for color `clr`?
		# Going through all the ancestors of `c` to check
		for (ia,(a, ac)) in enumerate(zip(c.anc, c.anccolor))
			if a == n && ac[clr]
				push!(idx, ic)
				break
			end
		end
	end
	return idx
end
"""
	get_children(a::ARGNode, clr::Int64)

Get children of `a` for color `clr`. 
"""
get_children(a::ARGNode, clr::Int64) = a.children[get_children_index(a, clr)]
child(a::ARGNode, clr::Int) = a.children[get_children_index(a, clr)]


function get_ancestor_index(a::ARGNode, clr::Int)
	for (i, ac) in enumerate(a.anccolor)
		if ac[clr]
			return i
		end
	end
	return nothing
end

ancestor(a::ARGNode, clr::Int) = a.anc[get_ancestor_index(a, clr)]
