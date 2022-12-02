# ARGTools

Julia package for backward time simulation of the reassortment and coalescence process of $k \geq 1$-genetic segments. ARGTools ancestral reassortment graphs (ARGs), individual phylogenetic trees for each segment and maximally compatible clades (MCCs). This package is primarily used for generation data [TreeKnit](https://pierrebarrat.github.io/TreeKnit.jl). For more information about ARGs and MCCs see TreeKnit's [documentation](https://pierrebarrat.github.io/TreeKnit.jl/) and the following paper:

> TreeKnit: Inferring Ancestral Reassortment Graphs of influenza viruses   
> Pierre Barrat-Charlaix, Timothy G. Vaughan, Richard A. Neher
> bioRxiv 2021.12.20.473456; doi: https://doi.org/10.1101/2021.12.20.473456

Please cite this article if you use TreeKnit! 

## Simulations

ARGTools simulates a backward time coalescence and reassortment process of populations. Branch lengths correspond to generations. Coalescence is a stochastic process with rate parameter 
$$\gamma_c = \left(\frac{n_c}{2}\right)^{\alpha} \frac{(n_c-1)}{2N},$$
$N$ is the population size and $n_c$ is the number of nodes available for coalescence, $\alpha$ is a parameter that determines tree structure. In the standard Kingman model $\alpha = 1$, in ARGTools this corresponds to `simtype = :kingman`. As TreeKnit is primarily used for influenza data ARGTools also allows the user the option to simulate trees that are more flu-like with `simtype = :flu`, this method sets $\alpha = 0.2$. When combined with reassortment this model leads to the unfortunate side-effect that trees with higher reassortment rates also have longer branches. Reassortment is also modeled as stochastic process with rate parameter $$\gamma_r = r n_r.$$

```
# Parameters of the ARG simulation
N = 10_000 # pop size
n = 6 # Number of lineages at time t=0 (number of lineages decreases with time)
simtype = :kingman
r = 0.2 # Absolute reassortment rate

# Simulating the ARG - 2 segments
arg = ARGTools.SimulateARG.simulate(N, r, n; simtype);

# The trees for the 2 segments
t1, t2 = ARGTools.trees_from_ARG(arg; node_data = TreeTools.MiscData);

# The real MCCs
rMCCs = ARGTools.MCCs_from_arg(arg)
```

To perform simulations ARGTools defines two (not disjoint) sets of nodes; nodes available for reassortment ($n_r$) and nodes available for coalescence ($n_c$). By default all nodes are available for coalescence and the simulation is over when there is only one node left. A node is only available for reassortment when it has more than 1 segment. At starting time, generation $t=0$ all nodes have $k>1$-segments. In the code segments are referred to as colors. 

When reassortment occurs at generation $t$ a random node $n$ is chosen from the set of nodes available for reassortment ($n_r$) and is replaced by two ancestral nodes in the next generation $t+1$. This replacement occurs in both the set $n_r$ and $n_c$. Node $n$ has at least $2$ segments, these are randomly split amongst the two ancestral nodes. If an ancestral node has under two segments it can no longer undergo reassortment and will be removed from the $n_r$ set. 

When coalescence occurs at generation $t$ two random nodes are selected from the set of nodes available for coalescence ($n_c$) and are replaced by a new node with the union of their segments in the next generation $t+1$. The simulation is finished when there is only one node in the set $n_c$ remaining. 

The time to the next coalescent or recombination event is determined by sampling from the exponential distribution $$t \sim Exp(\gamma_c + \gamma_r).$$ Whether the next event is a reassortment or coalescence event is determined by random sampling, if $$rand((0,1)) \leq \frac{\gamma_c}{\gamma_c + \gamma_r}$$ a coalescence has occurred, otherwise a recombination event has occurred.  

## Extended Newick 

ARGs are often visualized in extended newick format [Cardona, G., Rosselló, F. & Valiente, G. Extended Newick: it is time for a standard representation of phylogenetic networks. BMC Bioinformatics 9, 532 (2008). https://doi.org/10.1186/1471-2105-9-532](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-532). Currently, ARGTools does not allow exporting to extended newick format, however the `ARG_simulation_nwk.jl` gives an example how TreeKnit can be used to export ARGs in this format. 