using ARGTools
using Test
using TreeTools


println("### Extended newick")

file = "$(dirname(pathof(ARGTools)))/..//tests/trees/small_arg_test.nwk"
s = read(file, String)       
arg = ARGTools.parse_extended_newick(s)
label = ARGTools.parse_node_label("ARGNode_3NIZmScq#H1[&segments={0}]")
@testset "parse label" begin
    @test label[1] == "ARGNode_3NIZmScq#H1"
end
ext_newick = ARGTools.extended_newick(arg; pruned_singletons=false)
@testset "extended_newick - simple" begin
    @test s == ext_newick
end

## test simulations with pruned_singletons
N = 10_000 # pop size
n = 4 # Number of lineages
simtype = :kingman
r = (n/2)^0.2 / N # reassortment rate (here same as coalescence rate)
K = 2
# Simulating the ARG - 2 segments 
arg = ARGTools.SimulateARG.simulate(N, r, n; K, simtype, seed=12);
ext_newick = ARGTools.extended_newick(arg; pruned_singletons=true)
true_ext_nwk = "((((2_0[&segments={0,1}],1_0[&segments={0,1}])internal_4[&segments={0,1}])hybrid_node_1#H1[&segments={0}],"*
"(3_0[&segments={0,1}],4_0[&segments={0,1}])internal_1[&segments={0,1}])internal_10[&segments={1}],hybrid_node_1#H1[&segments={1}])internal_11[&segments={1}];"
@testset "extended_newick - simulated ARGs, K=2" begin
    @test ext_newick == true_ext_nwk
end

K = 3
# Simulating the ARG - 3 segments 
arg = ARGTools.SimulateARG.simulate(N, r, n; K, simtype, seed=12);
ext_newick = ARGTools.extended_newick(arg; pruned_singletons=true)
true_ext_nwk_3 = "(((((((2_0[&segments={0,1,2}])hybrid_node_1#H1[&segments={0}],(1_0[&segments={0,1,2}])hybrid_node_2#H2[&segments={2}])internal_9[&segments={0,2}])hybrid_node_3#H3[&segments={2}],"*
"((3_0[&segments={0,1,2}],4_0[&segments={0,1,2}])internal_1[&segments={0,1,2}])hybrid_node_4#H4[&segments={0,1}])internal_13[&segments={0,1,2}],"*
"(((hybrid_node_1#H1[&segments={1,2}],hybrid_node_2#H2[&segments={0,1}])internal_6[&segments={0,1,2}],hybrid_node_4#H4[&segments={2}])internal_8[&segments={0,1,2}])hybrid_node_5#H5[&segments={0}])internal_14[&segments={0,1,2}],"*
"hybrid_node_3#H3[&segments={0}])internal_15[&segments={1,2}],hybrid_node_5#H5[&segments={1,2}])internal_17[&segments={2}];"
@testset "extended_newick - simulated ARGs, K=3" begin
    @test ext_newick == true_ext_nwk_3
end

file = "$(dirname(pathof(ARGTools)))/..//tests/trees/arg_different_roots.nwk"
s = read(file, String)         
arg = ARGTools.parse_extended_newick(s)
ext_newick = ARGTools.extended_newick(arg; pruned_singletons=false)
@testset "extended_newick - different roots" begin
    @test s == ext_newick
end