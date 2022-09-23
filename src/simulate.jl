module SimulateARG

using ARGTools, TreeTools
using Random
using Distributions

export simulate, SimState, SimParam

let verbose::Bool = false, vverbose::Bool = false
	global v() = verbose
	global vv() = vverbose
	global set_verbose(v) = (verbose = v)
	global set_vverbose(v) = (vverbose = v)
end

let disc_t::Int, exact_t::Float64
	global reset_discrete_t() = (disc_t = 1)
	global inc_discrete_t() = (disc_t += 1)
	global get_discrete_t() = (disc_t)
	global reset_exact_t() = (exact_t=0.)
	global inc_exact_t(t) = (exact_t += t)
	global get_exact_t() = (exact_t)
end
# const global K = 2 # Two trees

# What matters is n/(N*r) where n is the size of the current ancestry
# If 2/(Nr) is small, it's likely that a split happens before a coalescence before the last two ancestors coalesce. 
struct SimParam
	N::Int # Global pop. size
	r::Float64 # Per individual recombination rate
	ρ::Float64 # Population reassortment rate (i.e. N*r)
	n0::Int # Initial sample size
	Tmax::Int # Maximal number of generations
	K::Int # Degree of the ARG
	nmax::Int # Maximal sample size
	s::Float64 # Rate at which leaves are added to the population
end
mutable struct SimState
	arg::ARG
	node_times::Dict{String, Float64} # Exact time of each node
	found_color_root::Array{Bool,1} # Roots found for different colors. 
	pop_per_color::Array{Int,1} # Number of current ancestors that have a given color
	eligible_for_reassortment::Array{String,1} # Nodes of the arg eligible for reassortment
	eligible_for_coalescence::Array{String,1} # 
end
mutable struct SimSnapshot # Used to record history of the simulation
	N::Real # Population size
	nc::Int # Number of nodes eligible for coalescence
	nr::Int # Number of nodes eligible for reassortment
	pop_per_color::Array{Int,1} # Number of nodes per color
	event::Union{Symbol, Nothing} # :coa or :split 
	discrete_t::Int # Number of the event
	exact_t::Float64 # Time of the event
end

"""
	simulate(param::SimParam)
"""
function simulate(
	param::SimParam;
	verbose=false, 
	vverbose = false,
	prune_singletons=true,
	popvar=t->1, 
	output_history = false,
	simtype = :yule
)
	set_verbose(verbose)
	set_vverbose(vverbose)
	v() && println("Simulating an ARG (simtype $simtype)")
	# 
	simstate = initiate(param)
	reset_discrete_t()
	reset_exact_t()
	simhistory = SimSnapshot[SimSnapshot(
		param.N,
		length(simstate.eligible_for_coalescence),
		length(simstate.eligible_for_reassortment), 
		simstate.pop_per_color,
		nothing,
		get_discrete_t(),
		get_exact_t()
	)]
	n_tot = param.n0

	sim = true
	while sim
		N = param.N * popvar(get_exact_t())
		s = if n_tot >= param.nmax
			0
		elseif true in simstate.found_color_root
			0
		else
			param.s
		end

		τ, etype = choose_event(
			param.r,
			N,
			length(simstate.eligible_for_coalescence), 
			length(simstate.eligible_for_reassortment),
			s,
			simtype,
		)
		vv() && println()
		v() && println("Event : $etype - Time $τ")
		if etype == :coa
			coalescence!(simstate, τ)
		elseif etype == :split
			split!(simstate, τ)
		elseif etype == :add_leaf
			add_leaf!(simstate, τ, param.K)
			n_tot += 1
		end
		inc_discrete_t()
		inc_exact_t(τ)
		vv() && println("Eligible for coalescence: $(length(simstate.eligible_for_coalescence)) nodes - $(simstate.eligible_for_coalescence)")
		vv() && println("Eligible for reassortment: $(length(simstate.eligible_for_reassortment)) nodes - $(simstate.eligible_for_reassortment)")
		# Storing Snapshot
		push!(simhistory, SimSnapshot(N, length(simstate.eligible_for_coalescence), 
			length(simstate.eligible_for_reassortment), 
			simstate.pop_per_color,
			etype,
			get_discrete_t(),
			get_exact_t()))
		# Stop ? 
		if halt_condition(simstate, param.Tmax)
			sim = false
		elseif get_discrete_t() > param.Tmax # No more recomb to finish simulation
			v() && println("Stopping: Maximum number of iterations reached ($(param.Tmax)).")
			ρ = 0
		end
	end
	set_roots_ancestry!(simstate.arg)
	prune_singletons && ARGTools.prune_singletons!(simstate.arg)
	ARGTools.check_arg(simstate.arg)
	if output_history
		return simstate.arg, simhistory
	else
		return simstate.arg
	end
end
"""
	simulate(
		N,r,n0;
		K=2,
		s = 0,
		nmax = n0,
		Tmax=1e6,
		verbose=false,
		vverbose=false,
		popvar=t->1,
		prune_singletons=true,
		output_history=false,
		simtype=:yule,
	)
"""
function simulate(
	N,r,n0;
	K=2, 
	s = 0,
	nmax = n0,
	Tmax=1e6,
	verbose=false, 
	vverbose=false, 
	popvar=t->1, 
	prune_singletons=true,
	output_history=false,
	simtype=:yule,
)

	simulate(
		SimParam(N, r, N*r, n0, Tmax, K, nmax, s);
		verbose, vverbose, popvar, prune_singletons, output_history, simtype,
	)
end


"""
	initiate(param::SimParam)

Create `param.n0` `ARGNode` structures with uninitialized parents. 
"""
function initiate(param::SimParam)
	arg = ARG(degree=param.K)
	for i in 1:param.n0
		an = ARGNode(
			degree=param.K,
			anc = Vector{Any}(undef, 0),
			label="$(i)_0",
			isroot = zeros(Bool, param.K),
			isleaf = true
		)
		arg.nodes[an.label] = an
		arg.leaves[an.label] = an
	end

	return SimState(
		arg,
		Dict(x=>0. for x in keys(arg.nodes)),
		zeros(Bool,param.K),
		param.n0*ones(Int,param.K),
		collect(keys(arg.nodes)),
		collect(keys(arg.nodes))
	)
end

"""
	halt_condition(simstate::SimState, Tmax::Int)
"""
function halt_condition(simstate::SimState, Tmax::Int)
	if (&)(simstate.found_color_root...) 
		v() && println("Stopping: Roots for all colors have been found.")
		return true
	end
	return false
end

"""
	choose_event(r, N, n, nr, s, simtype)

Choose type of the next event. Return the time to the time to next event as well as its type `:coa` or `:split`.  
r, n and nr are resp. the population reassortment rate, the number of nodes available for coalescence and the number of nodes available for reassortment. 
"""
function choose_event(r, N, n, nr, s, simtype)
	iTr = r*nr
	if simtype == :kingman 
		iTc = n*(n-1) /2. /N
	elseif simtype == :yule
		iTc = (n-1)/2. /N
	elseif simtype == :flu
		iTc = n^0.2 * (n-1) /2. /N
	end
	iTs = s

	tc = rand(Exponential(1/iTc))
	tr = rand(Exponential(1/iTr))
	ts = rand(Exponential(1/iTs))

	if tc <= min(tr, ts)
		return tc, :coa
	elseif tr <= min(tc, ts)
		return tr, :split
	elseif ts <= min(tc, tr)
		return ts, :add_leaf
	end
end

function add_leaf!(simstate::SimState, t, degree)
	τ = get_exact_t() + t

	an = ARGNode(
		label = "$(length(simstate.arg.leaves)+1)_0",
		anc = Vector{Any}(undef, 0),
		isroot = zeros(Bool, degree),
		isleaf = true,
		degree = degree,
		tau = Vector{Float64}(undef, 0)
	)
	simstate.arg.nodes[an.label] = an
	simstate.arg.leaves[an.label] = an

	simstate.node_times[an.label] = τ
	push!(simstate.eligible_for_coalescence, an.label)
	push!(simstate.eligible_for_reassortment, an.label)
	for i in 1:length(simstate.pop_per_color)
		simstate.pop_per_color[i] += 1
	end

	return an
end

"""
	coalescence!(simstate::SimState, t)

Find two nodes to coalesce in `simstate.arg`. 
"""
function coalescence!(simstate::SimState, t)
	# Chosing two nodes
	if length(simstate.eligible_for_coalescence) == 1
		@error "Can't coalesce singleton population"
	end
	i = rand(1:length(simstate.eligible_for_coalescence))
	j = rand(1:length(simstate.eligible_for_coalescence))
	while j==i 
		j = rand(1:length(simstate.eligible_for_coalescence))
	end
	n1 = simstate.arg.nodes[simstate.eligible_for_coalescence[i]]
	n2 = simstate.arg.nodes[simstate.eligible_for_coalescence[j]]
	t1 = t + get_exact_t() - simstate.node_times[n1.label]
	t2 = t + get_exact_t() - simstate.node_times[n2.label]

	# Coalesce
	new_node = coalescence!(simstate.arg, n1, n2, t1, t2, simstate)
	for (i,c) in enumerate(n1.color .& n2.color)
		simstate.pop_per_color[i] -= c # If n1 and n2 have color c, remove one 
		if simstate.pop_per_color[i] == 1 && !simstate.found_color_root[i]
			# new_node is the root for color i
			simstate.found_color_root[i] = true
			new_node.isroot[i] = true
			v() && println("Found root for color $i.")
			simstate.arg.root[i] = new_node
		end
	end

	# Fix simstate
	# `simstate.eligible_for_coalescence` 
	deleteat!(simstate.eligible_for_coalescence, (min(i,j), max(i,j)))
	push!(simstate.eligible_for_coalescence, new_node.label)
	# simstate.eligible_for_reassortment
	# Note: if `new_node` is a root, I remove it from eligible_for_reassortment. This greatly simplifies the code... 
	idx = findall(x->in(x, (n1.label, n2.label)), simstate.eligible_for_reassortment)
	deleteat!(simstate.eligible_for_reassortment, idx)
	if new_node.degree > 1 && !(|)(new_node.isroot...)
		push!(simstate.eligible_for_reassortment, new_node.label)
	end

	# Node times
	simstate.node_times[new_node.label] = simstate.node_times[n1.label] + n1.tau[1]
	delete!(simstate.node_times, n1.label)
	delete!(simstate.node_times, n2.label)

	# `pop_per_color`
	# If new_node is the root for some colors and does not have any other color, remove it from eligible_for_coalescence
	if sum(new_node.isroot) == new_node.degree
		idx = findall(x->x==new_node.label, simstate.eligible_for_coalescence)
		deleteat!(simstate.eligible_for_coalescence, idx)
	end

	# Add new node to ARG
	simstate.arg.nodes[new_node.label] = new_node
end
"""	
	coalescence!(arg::ARG, n1::ARGNode, n2::ARGNode, t1, t2)

Coalesce `n1` and `n2` into a single `ARGNode`, and adds it to `arg`. 
"""
function coalescence!(arg::ARG, n1::ARGNode, n2::ARGNode, t1, t2, simstate)
	vv() && println("Attempting to coalesce $(n1.label) and $(n2.label)")
	vv() && println("Respective branch lenghts: $t1 - $t2")
	# Parent node
	new_label = "internal_$(get_discrete_t())"
	new_color = convert(Array{Bool,1}, n1.color .| n2.color)
	new_color .*= (!).(simstate.found_color_root) # Colors of roots don't propagate up
	new_node = ARGNode(children = [n1,n2],
		anc = Array{Union{ARGNode{TreeTools.MiscData},Nothing}}(nothing, 0),
		color = new_color,
		degree = sum(new_color),
		label=new_label,
		isroot=zeros(Bool, length(new_color)),
		isleaf=false)
	vv() && println("New node $(new_label) with color $new_color")
	# Children nodes
	for (n,t) in zip((n1,n2), (t1,t2))
		push!(n.anc, new_node)	# Pushing in case n is the root node of some color. In this case it already has an ancestor (nothing)
		push!(n.tau, t)
		push!(n.data, TreeTools.MiscData())
		push!(n.anccolor, copy(n.color) .* new_node.color)
	end
	#
	return new_node
end

"""
"""
function split!(simstate::SimState, t)
	# Choose a node
	n = simstate.arg.nodes[sample(simstate.eligible_for_reassortment)]

	# Choose a color split
	# Use ghost segments for other colors so that reassortment rate stays constant during the simulation
	# Split all segments, including ghosts
	colors = shuffle(1:simstate.arg.degree)
	bar = rand(1:(simstate.arg.degree-1))
	cs1, cs2 = sort!(colors[1:bar]), sort!(colors[(bar+1):end])

	# Map this to colors in `n`
	c1 = cs1[findall(c -> n.color[c], cs1)]
	c2 = cs2[findall(c -> n.color[c], cs2)]

	# println("color split: ", cs1, " / ", cs2)
	# println("node colors: ", n.color)
	# println("final split: ", c1, " / ", c2)
	# println()
	if isempty(c1) || isempty(c2)
		return n
	else
		return split!(simstate, n, c1, c2, t)
	end
end

function split!(simstate::SimState, n, c1, c2, t)
	@assert !isempty(c1) && !isempty(c2)
	# Split it backwards
	a1, a2 = split!(n, c1, c2, t + get_exact_t() - simstate.node_times[n.label])

	# node times
	for a in (a1,a2)
		simstate.node_times[a.label] = t + get_exact_t()
	end
	delete!(simstate.node_times, n.label)

	# update eligible_for_reassortment
	deleteat!(simstate.eligible_for_reassortment, findfirst(x->x==n.label, simstate.eligible_for_reassortment))
	for a in (a1,a2)
		a.degree > 1 && push!(simstate.eligible_for_reassortment, a.label)
	end

	# update eligible for coalescence
	deleteat!(simstate.eligible_for_coalescence, findfirst(x->x==n.label, simstate.eligible_for_coalescence))
	push!(simstate.eligible_for_coalescence, a1.label, a2.label)

	# Add ancestors to arg
	simstate.arg.nodes[a1.label] = a1
	simstate.arg.nodes[a2.label] = a2

	return a1, a2
end
function split!(n::ARGNode, c1::Array{Int,1}, c2::Array{Int,1}, t)
	vv() && println("Attempting to split node $(n.label) with color $(n.color).")
	vv() && println("Colors $c1 going left and $c2 going right.")
	# Parent nodes
	new_label1 = "internal1_$(get_discrete_t())"
	new_label2 = "internal2_$(get_discrete_t())"
	new_clr1 = ARGTools._color(c1, length(n.color))
	new_clr2 = ARGTools._color(c2, length(n.color))
	a1 = ARGNode(children = [n],
		anc = Array{Union{ARGNode{TreeTools.MiscData},Nothing}}(nothing, 0),
		color = new_clr1,
		degree = sum(new_clr1),
		label = new_label1,
		isroot=zeros(Bool, length(new_clr1)),
		isleaf = false)
	a2 = ARGNode(children = [n],
		anc = Array{Union{ARGNode{TreeTools.MiscData},Nothing}}(nothing, 0),
		color = new_clr2,
		degree = sum(new_clr2),
		label = new_label2,
		isroot=zeros(Bool, length(new_clr2)),
		isleaf = false)
	# Child
	push!(n.anc, a1)
	push!(n.anccolor, copy(new_clr1))
	push!(n.tau, t)
	push!(n.data, TreeTools.MiscData())

	push!(n.anc, a2)
	push!(n.anccolor, copy(new_clr2))
	push!(n.tau, t)
	push!(n.data, TreeTools.MiscData())
	#
	return a1, a2
end

"""
	set_roots_ancestry!(arg::ARG)

For nodes that are roots for a color `c`, removes all ancestry corresponding to `c` and adds `nothing` as an ancestor for `c`.
"""
function set_roots_ancestry!(arg::ARG)
	for (c,ar) in enumerate(arg.root)
		# Cut branches for indices `todel` and color `c`
		todel = Int[]
		for (i,clr) in enumerate(ar.anccolor)
			clr[c] && push!(todel, i)
		end
		!isempty(todel) && ARGTools.cut_branch!(ar, todel, c)
		# Adding `nothing` as ancestor
		push!(ar.anc, nothing)
		push!(ar.anccolor, ARGTools._color(c, arg.degree))
		push!(ar.tau, missing)
		push!(ar.data, TreeTools.MiscData())
	end
	ARGTools.prune_lone_nodes!(arg)
end


end
