using Pkg
Pkg.status
using TreeKnit
using TreeKnit.SRG
using TreeTools
using ARGTools

# Utility function
function get_r(ρ, n, N, simtype::Symbol)
    if simtype == :kingman
        return ρ * (n/2) / N
    elseif simtype == :yule
        return ρ / N
    elseif simtype == :flu
    	return ρ * (n/2)^0.2 / N
    else
        @error "Unrecognized `simtype`."
    end
end

# Parameters of the ARG simulation
N = 10_000 # pop size
n = 6 # Number of lineages
ρ = 1 # Reassortment rate scaled to coalescence rate
simtype = :kingman
r = get_r(ρ, n, N, simtype) # Absolute reassortment rate

# Simulating the ARG - 2 segments
prune_singletons=false
arg = ARGTools.SimulateARG.simulate(N, r, n; simtype, prune_singletons);
#ARGTools.write(outfolder * "argtools_arg.nwk", arg)

# The trees for the 2 segments
t1, t2 = ARGTools.trees_from_ARG(arg);

# The real MCCs
rMCCs = ARGTools.MCCs_from_arg(arg)

# Writing results
outfolder = "ARG_n$(n)_rho$(ρ)_simtype_$(simtype)/"
mkpath(outfolder)
write_newick(outfolder * "tree1.nwk", t1)
write_newick(outfolder * "tree2.nwk", t2)
write_mccs(outfolder * "real_MCCs.dat", rMCCs)

# Write ARG as extended newick
# Currently ARGTools does not have the functionality to write extended newick
# this is implemented in TreeKnit -- but only for 2 trees!!!
arg_TK, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, rMCCs)
TreeKnit.write(outfolder * "arg.nwk", arg_TK)


