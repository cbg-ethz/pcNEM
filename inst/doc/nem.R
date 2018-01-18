### R code from vignette source 'nem.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: no.nonsense
###################################################
rm(list=ls())


###################################################
### code chunk number 2: nem.Rnw:94-96
###################################################
library(nem)
data("BoutrosRNAi2002")


###################################################
### code chunk number 3: nem.Rnw:121-122
###################################################
res.disc <- nem.discretize(BoutrosRNAiExpression,neg.control=1:4,pos.control=5:8,cutoff=.7)


###################################################
### code chunk number 4: nem.Rnw:125-133
###################################################
disc <- cbind(res.disc$neg,res.disc$pos,res.disc$dat)
e.2fold <- BoutrosRNAiExpression[res.disc$sel,]

#--- hierarchisch clustern
library(e1071)
hc    <- hclust(as.dist(hamming.distance(disc[,9:16])))
e.2fold <- e.2fold[hc$order, ]
disc  <- disc [hc$order, ]


###################################################
### code chunk number 5: data_cont
###################################################
#--- CONTINUOUS DATA
#pdf("data_cont.pdf",width=5,height=13)
par(las=2,mgp=c(5.5,1,0),mar=c(6.7,7,4,1),cex.lab=1.7,cex.main=2)
image(x   = 1:ncol(e.2fold),
y   = 1:nrow(e.2fold),
z   = scale(t(e.2fold)),
main= "Original data",
xlab= "Experiments",
xaxt= "n",
ylab= "E-genes",
yaxt= "n",
col = gray(seq(0,.95,length=10))
)
abline(v=c(4,8,10,12,14)+.5)
axis(1,1:ncol(e.2fold),colnames(e.2fold))
axis(2,1:nrow(e.2fold),rownames(e.2fold))
#dev.off()


###################################################
### code chunk number 6: data_disc
###################################################
#--- DISCRETE DATA
#pdf("data_disc.pdf",width=5,height=13)
par(las=2,mgp=c(5.5,1,0),mar=c(6.7,7,4,1),cex.lab=1.7,cex.main=2)
image(x   = 1:ncol(disc),
      z   = t(disc),
      main= "Discretized data",
      xlab= "Experiments",
      xaxt= "n",
      ylab= "",
      yaxt= "n",
      col = gray(seq(.95,0,length=10))
      )
abline(v=c(4,8,10,12,14)+.5)
axis(1,1:ncol(e.2fold),colnames(e.2fold))
#dev.off()


###################################################
### code chunk number 7: nem.Rnw:190-191
###################################################
res.disc$para


###################################################
### code chunk number 8: nem.Rnw:224-227
###################################################
hyper = set.default.parameters(unique(colnames(res.disc$dat)), para=res.disc$para)
result <- nem(res.disc$dat,inference="search",control=hyper, verbose=FALSE)
result


###################################################
### code chunk number 9: nem.Rnw:233-236
###################################################
plot.nem(result,what="graph")
plot.nem(result,what="mLL")
plot.nem(result,what="pos")


###################################################
### code chunk number 10: nem.Rnw:239-254
###################################################
Sgenes <- unique(colnames(res.disc$dat))
models <- enumerate.models(Sgenes)
best5 <- -sort(-unique(result$mLL))[1:5]
col<-c("yellow","yellow","green","blue")
names(col) = Sgenes
library(Rgraphviz)
for (i in 1:5) {
   graph <- as(models[[which(result$mLL == best5[i])[1]]]-diag(4),"graphNEL")
   pdf(file=paste("topo",i,".pdf",sep=""))
   par(cex.main=5)
   plot(graph,
        nodeAttrs=list(fillcolor=col),
        main=paste("-",i, "-"))        
   dev.off()
   }


###################################################
### code chunk number 11: scores1
###################################################
plot.nem(result,what="mLL")


###################################################
### code chunk number 12: pos1
###################################################
plot.nem(result,what="pos")


###################################################
### code chunk number 13: nem.Rnw:294-298
###################################################
hyper$type="FULLmLL"
hyper$hyperpara=c(1,9,9,1)
result2 <- nem(res.disc$dat,inference="search",control=hyper,verbose=FALSE)
result2


###################################################
### code chunk number 14: nem.Rnw:301-311
###################################################
best5 <- -sort(-unique(result2$mLL))[1:5]
for (i in 1:5) {
   graph <- as(models[[which(result2$mLL == best5[i])[1]]]-diag(4),"graphNEL")
   pdf(file=paste("topo2",i,".pdf",sep=""))
   par(cex.main=5)
   plot(graph,
        nodeAttrs=list(fillcolor=col),
        main=paste("-",i, "-"))
   dev.off()
   }


###################################################
### code chunk number 15: scores2
###################################################
plot.nem(result2,what="mLL")


###################################################
### code chunk number 16: pos2
###################################################
plot.nem(result2,what="pos")


###################################################
### code chunk number 17: nem.Rnw:352-354
###################################################
resultPairs <- nem(res.disc$dat,inference="pairwise",control=hyper, verbose=FALSE)
resultPairs


###################################################
### code chunk number 18: graph3
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultPairs$graph)
plot.nem(resultPairs, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 19: nem.Rnw:371-373
###################################################
resultTriples <- nem(res.disc$dat,inference="triples",control=hyper, verbose=FALSE)
resultTriples


###################################################
### code chunk number 20: graph4
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultTriples$graph)
plot.nem(resultTriples, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 21: nem.Rnw:388-390
###################################################
resultGreedy <- nem(res.disc$dat,inference="nem.greedy",control=hyper, verbose=FALSE)
resultGreedy


###################################################
### code chunk number 22: graph44
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultGreedy$graph)
plot.nem(resultGreedy, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 23: nem.Rnw:429-431 (eval = FALSE)
###################################################
## resultMN <- nem(res.disc$dat,inference="ModuleNetwork",control=hyper, verbose=FALSE) # this will do exactly the same as the exhaustive search
## resultMN.orig <- nem(res.disc$dat,inference="ModuleNetwork.orig",control=hyper, verbose=FALSE) # this will do exactly the same as the exhaustive search


###################################################
### code chunk number 24: nem.Rnw:439-441
###################################################
resGreedyMAP <- nem(BoutrosRNAiLods, inference="nem.greedyMAP", control=hyper, verbose=FALSE)
resGreedyMAP


###################################################
### code chunk number 25: graph55
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resGreedyMAP$graph)
plot.nem(resGreedyMAP, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 26: nem.Rnw:466-470
###################################################
hyper$Pm = diag(4)
hyper$lambda = 10
resultRegularization <- nem(res.disc$dat, inference="search", control=hyper, verbose=FALSE)
resultRegularization


###################################################
### code chunk number 27: graph6
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultRegularization$graph)
plot.nem(resultRegularization, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 28: nem.Rnw:496-497
###################################################
resultModsel <- nemModelSelection(c(0.01,1,100),res.disc$dat, control=hyper,verbose=FALSE)


###################################################
### code chunk number 29: graph7
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultModsel$graph)
plot.nem(resultModsel, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 30: nem.Rnw:513-517
###################################################
hyper$lambda=0	# set regularization parameter to 0
hyper$Pm	# this is our prior graph structure
resultBayes <- nem(res.disc$dat, control=hyper, verbose=FALSE) # now we use Bayesian model selection to incorporate the prior information
resultBayes


###################################################
### code chunk number 31: graph77
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultBayes$graph)
plot.nem(resultBayes, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 32: nem.Rnw:537-543 (eval = FALSE)
###################################################
## logdensities = getDensityMatrix(myPValueMatrix,dirname="DiagnosticPlots")
## nem(logdensities[res.disc$sel,],type="CONTmLLBayes",inference="search")
## nem(logdensities[res.disc$sel,],type="CONTmLLMAP",inference="search")
## 
## preprocessed <- nem.cont.preprocess(BoutrosRNAiExpression,neg.control=1:4,pos.control=5:8)
## nem(preprocessed$prob.influenced,type="CONTmLL",inference="search")


###################################################
### code chunk number 33: nem.Rnw:557-564 (eval = FALSE)
###################################################
## mydat = filterEGenes(Porig, logdensities) # a-priori filtering of E-genes
## 
## hyper$selEGenes = TRUE
## hyper$type = "CONTmLLBayes"
## resAuto = nem(mydat,control=hyper) # use filtered data to estimate a network; perform automated subset selection of E-genes with tuning of the parameter delta
## 
## most.relevant = getRelevantEGenes(as(resAuto$graph, "matrix"), mydat)$selected # returns all E-genes with partial log-likelihood > 0.


###################################################
### code chunk number 34: nem.Rnw:578-582 (eval = FALSE)
###################################################
## significance=nem.calcSignificance(disc$dat,res) # assess statistical significance
## bootres=nem.bootstrap(res.disc$dat) # bootstrapping on E-genes
## jackres=nem.jackknife(res.disc$dat) # jackknife on S-genes
## consens=nem.consensus(res.disc$dat) # bootstrap & jackknife on S-genes


###################################################
### code chunk number 35: nem.Rnw:589-591
###################################################
hyper$mode="binary_Bayesian"
resultBN = nem(res.disc$dat, inference="BN.greedy", control=hyper)


###################################################
### code chunk number 36: graph8
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultBN$graph)
plot.nem(resultBN, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 37: nem.Rnw:605-610
###################################################
hyper$mcmc.nburnin = 10 # much to few in practice
hyper$mcmc.nsamples = 100  
hyper$type = "CONTmLLDens"
hyper$Pm = NULL
resulteminem = nem(BoutrosRNAiLods, inference="mc.eminem", control=hyper)


###################################################
### code chunk number 38: graph8
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resultBN$graph)
plot.nem(resulteminem, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 39: nem.Rnw:625-630 (eval = FALSE)
###################################################
## data("Ivanova2006RNAiTimeSeries")
## dim(dat)
## control = set.default.parameters(dimnames(dat)[[3]], para=c(0.1,0.1))
## net = nem(dat, inference="dynoNEM", control=control)
## plot.nem(net, SCC=FALSE, plot.probs=TRUE)


###################################################
### code chunk number 40: nem.Rnw:648-651 (eval = FALSE)
###################################################
## data(SahinRNAi2008)
## control = set.default.parameters(setdiff(colnames(dat.normalized),"time"), map=map.int2node, type="depn",debug=FALSE) # set mapping of interventions to perturbed nodes
## net = nem(dat.normalized, control=control) # greedy hillclimber to find most probable network structure


###################################################
### code chunk number 41: nem.Rnw:659-661
###################################################
resEdgeInf = infer.edge.type(result, BoutrosRNAiLogFC)
plot.nem(resEdgeInf, SCC=FALSE)


###################################################
### code chunk number 42: graph100
###################################################
col<-c("yellow","yellow","green","blue")
names(col) = nodes(resEdgeInf$graph)
plot.nem(resEdgeInf, nodeAttrs=list(fillcolor=col), SCC=FALSE)


###################################################
### code chunk number 43: nem.Rnw:675-676
###################################################
plotEffects(res.disc$dat,result)


###################################################
### code chunk number 44: plot_effects
###################################################
plotEffects(res.disc$dat,result)


###################################################
### code chunk number 45: nem.Rnw:696-698
###################################################
result.scc <- SCCgraph(result$graph,name=TRUE)
plot(result.scc$graph)


###################################################
### code chunk number 46: scc
###################################################
# col2<-c("yellow","blue","green")
# names(col2) = nodes(result.scc$graph)
plot(result.scc$graph)#,nodeAttrs=list(fillcolor=col2))


###################################################
### code chunk number 47: nem.Rnw:728-729
###################################################
toLatex(sessionInfo())


