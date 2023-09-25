require(ape)
require(Rcpp)

generateTest <- function() {
  ntip<-25
  #test_tree <- rtree(ntip, br = sample(0:10,100, replace=T), digits = 0, equiprob = TRUE)
  test_tree <- rtree(ntip, br = sample.int(10,100, replace = "TRUE"), digits = 0, equiprob = TRUE)
  #test_tree <- read.tree(text= "((1:10,2:7)15:6,((3:2,(4:7,5:9)18:1)17:6,((6:6,(7:8,8:5)21:9)20:8,(9:1,(10:7,(11:8,(12:2,13:1)25:4)24:10)23:2)22:10)19:8)16:7)14;")
  test_tree$tip.label <- seq(1, ntip)
  test_tree$node.label <- seq(ntip+1, ntip + test_tree$Nnode)
  
  newick_string <- write.tree(test_tree, append = FALSE, digits = 10, tree.names = FALSE)
  #newick_string <- "(1:5,((2:8,((3:6,4:8)20:6,((5:6,6:5)22:7,(7:9,8:9)23:3)21:3)19:10)18:2,((9:10,10:1)25:4,(11:8,((12:3,13:6)28:2,(14:9,15:3)29:7)27:3)26:9)24:9)17:2)16;"
  backwardParameterAlgorithm_aic(newick_string, 3)
  
  #plot.phylo(test_tree, use.edge.length = TRUE, show.node.label = TRUE)
  #edgelabels(test_tree$edge.length, bg = "yellow", col="black")
}

sourceCpp('~/Stats Research/R Files/pthreadsRcpp.cpp')
generateTest()

