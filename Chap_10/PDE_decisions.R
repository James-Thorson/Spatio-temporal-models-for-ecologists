
setwd( R'(C:\Users\James.Thorson\Desktop\Git\Spatio-temporal-models-for-ecologists\Chap_10)' )

#
library(diagram)

# Visualize graph
box.labels <- c( "Drift: \n Vector fields for advection",
            "Taxis: \n Covariates for preference function",
            "Diffusion: \n Covariates for diffusion rates",
            "Diffusion-taxis-drift \n partial differential equation",
            "Continuous space, \n continuous time",
            "Discretized space, \n continuous time",
            "tracks",
            "tracks",
            "densities",
            "Multivariate state-space model",
            "Matrix exponential in \n Hidden Markov model",
            "Matrix exponential in \n species distribution model" )
names(box.labels) = letters[seq_len(length(box.labels))]
box.size = ifelse( names(box.labels)=="d", 0.3, 0.15 )
box.type = c( "square",
              "square",
              "square",
              "circle",
              "square",
              "square",
              "square",
              "square",
              "square",
              "circle",
              "circle",
              "circle" )
pos = matrix( c(
  1, 5,
  2, 5,
  3, 5,
  2, 4,
  1, 3,
  3, 3,
  1, 2,
  2, 2,
  3, 2,
  1, 1,
  2, 1,
  3, 1
), byrow=TRUE, ncol=2)
pos[,1] = (pos[,1]-0.5) / 3
pos[,2] = (pos[,2]-0.5) / 5

# Construct M
M1 = array(0, dim=c(length(box.labels),length(box.labels)), dimnames=list(names(box.labels),names(box.labels)))
M1['a', 'd'] = ""
  M1['b', 'd'] = ""
  M1['c', 'd'] = ""
M1['d', 'e'] = ""
  M1['d', 'f'] = ""
M1['e', 'g'] = ""
M1['f', 'h'] = ""
  M1['f', 'h'] = ""
  M1['f', 'i'] = ""
M1['g', 'j'] = ""
  M1['h', 'k'] = ""
  M1['i', 'l'] = ""

# Plot
png( "Decision_tree_for_movement.png", width=6.5, height=4, res=200, un="in")
  par(  mar=c(0,2,2,0) )
  plotmat( t(M1), pos=pos, name = box.labels, box.type = box.type,
           curve = 0, lwd = 3, arr.pos=0.7, box.lwd = 1, box.size=box.size,
           box.prop = 0.2, box.cex=0.7, main="Decision tree for movement models", shadow.size=0 )
  #mtext( x=0, y=(5-0.5)/5, las=3, labels=c("Mechanisms") )
  mtext( side=2, at=(5:1-0.5)/5, las=3, text=c("Mechanisms", "Model", "Discretization", "Goal", "Method"), cex=0.9 )
dev.off()


