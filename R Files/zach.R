require(ape)
require(Rcpp)
require(doParallel)

#generateTest <- function() {
#cl <- makeCluster(detectCores())
#registerDoParallel(cl)
x <- c()
yFor <- c()
yBack <- c()
yBest <- c()
yAct <- c()
forSlopes <- c()
backSlopes <- c()
actSlopes <- c()
highestNTips <- c()
#for (j in 3:10) {
#print(j)
for (ntip in 13:13) {
  #print("Tips")
  #print(ntip)
  for (i in 1:1) {
    #print(i)
    #test_tree <- rtree(ntip, br = sample(0:10,100, replace=T), digits = 0, equiprob = TRUE)
    ntip <- 13
    test_tree <- rtree(ntip, br = sample.int(10,100, replace = "TRUE"), digits = 0, equiprob = TRUE)
    test_tree$tip.label <- seq(1, ntip)
    test_tree$node.label <- seq(ntip+1, ntip + test_tree$Nnode)
    
    newick_string <- write.tree(test_tree, append = FALSE, digits = 10, tree.names = FALSE)
    print(newick_string)
    
    testVec <- mainFunc(newick_string)
    
    # plot.phylo(test_tree, use.edge.length = FALSE, show.node.label = TRUE)
    # edgelabels(test_tree$edge.length, bg = "yellow", col="black")
    x <- c(x, ntip)
    yFor <- c(yFor, testVec[1])
    yBack <- c(yBack, testVec[2])
    yAct <- c(yAct, testVec[3])
    if (testVec[1] < testVec[2]) {
      yBest <- c(yBest, testVec[1])
    }
    else {
      yBest <- c(yBest, testVec[2])
    }
  }
}
#tempForSlope <- lm(yFor~x)$coefficients[[2]]
#forSlopes <- c(forSlopes, tempForSlope)

#tempBackSlope <- lm(yBack~x)$coefficients[[2]]
#backSlopes <- c(backSlopes, tempBackSlope)

#tempActSlope <- lm(yAct~x)$coefficients[[2]]
#actSlopes <- c(actSlopes, tempActSlope)

#highestNTips <- c(highestNTips, j)
#}
#plot(forSlopes ~ highestNTips, col = "blue", xlim = c(1,11), ylim = c(8.5,9.7), main = "Slopes vs. Max Number of Tips in Loop", xlab = "Max Number of Tips in Loop", ylab = "Slope of Regression")
#points(backSlopes ~ highestNTips, col = "red")
#legend(x="topleft",legend = c("Forward", "Backward", "Brute Force"), fill = c("blue", "red", "green"))
#points(actSlopes ~ highestNTips, col = "green")

#Group of code for finding difference between forward/backward and brute
#plot(yBack - yAct~ x, col = "orange", main = "Difference in Fast and Brute Force Algorithms", xlab = "Number of Tips", ylab = "Backwards - Brute Force", xlim = c(3,13), ylim = c(-1, 5))
#points(yFor - yAct~ x, col = "green", pch = 2)
#legend(x = "topleft", legend = c("Backward - Brute", "Forward - Brute"), fill = c("orange", "green"))
#abline(lm(yBack - yAct~ x), col = "orange")
#abline(lm(yFor-yAct~x), col = "green")
#print("Backward - Brute")
#print(summary(yBack-yAct))
#print(summary(lm(yBack-yAct~x)))
#print("Forward - Brute")
#print(summary(yFor-yAct))
#print(summary(lm(yFor-yAct~x)))

#Group of code for plotting the difference in best fast algo vs brute force
#plot(yBest - yAct ~ x, col = "purple", main = "Difference in Best Fast Algorithm\nand Brute Force Algorithm", xlab = "Number of Tips", ylab = "AIC Difference (Best Fast Algo - Brute)")
#abline(lm(yBest-yAct~x), col = "purple")
#print("AIC Differences")
#print(summary(yBest-yAct))
#print("Regression Output")
#print(summary(lm(yBest-yAct~x)))
#}

sourceCpp('~/Stats Research/R Files/zach_cpp.cpp')
start <- Sys.time()
generateTest()
end <- Sys.time()
print(end - start)

#y <- c(1,2,3,4)
#x <- c(1,2,3,4)
#testVar <- lm(y ~ x)
#newVar <- testVar$coefficients
#print(newVar[[2]])