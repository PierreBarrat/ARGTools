using ARGTools
using Test
using TreeTools


println("### Extended newick")

file = "$(dirname(pathof(ARGTools)))/..//tests/trees/small_arg_test.nwk"
s = read(file, String)       
arg = ARGTools.parse_extended_newick(s)
label = ARGTools.parse_node_label("ARGNode_3NIZmScq#H1[&segments={0}]")
ext_newick = ARGTools.extended_newick(arg; pruned_singletons=false)
@testset "extended_newick" begin
    @test s == ext_newick
end

file = "$(dirname(pathof(ARGTools)))/..//tests/trees/arg_different_roots.nwk"
s = read(file, String)         
arg = ARGTools.parse_extended_newick(s)
label = ARGTools.parse_node_label("ARGNode_3NIZmScq#H1[&segments={0}]")
ext_newick = ARGTools.extended_newick(arg; pruned_singletons=false)
@testset "extended_newick" begin
    @test s == ext_newick
end