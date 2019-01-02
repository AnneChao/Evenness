#*******************************************************************************
#********************************************************************************
# 
#
## R scipts "Evenness" for Chao and Ricotta (2019) paper. 
## This R code is for computing Figure 3, Figure 4 and Figure 5 of Chao and Ricotta (2019) paper.
# NOTE: The packages "ggplot2", "dplyr", "ade4", "reshape2", "ggpubr", "phytools", "ape" must be 
# installed and loaded before running the scripts. 
# 
#
#
# The following R scripts include two parts:  
# (1). Script for computing the profiles for six classes of evenness measures (Figure 3 in Chao and Ricotta's paper).
# (2). Script for computing the contribution of each species/node to taxonomic dissimarity and/or phylogenetic
#      dissimarity measures (Figures 4 and 5 in Chao and Ricotta's paper)
#
# 
#
#*******************************************************************************
#*******************************************************************************


####################################################################################
#
# (1). Computing the profiles for six classes of evenness measures (Figure 3)
#
####################################################################################
library(dplyr)
library(ade4)
library(phytools)
library(ggplot2)
library(ape)
library(reshape2)
library(ggpubr)

qD <- function(p,q){
  p <- p[p>0]
  if(q!=1){
    (sum(p^q))^(1/(1-q))
  }else{
    exp(-sum(p*log(p)))
  }
}

#' @param x is an observed species-by-assemblage frequency matrix. 
#' @param q.order the setting is 0 to 2 divided with 0.05 intervals.
#' @return is the values of all evenness index(6 types).


new_fun <- function(x,q.order){
  FUN <- qD
  n <- sum(x)
  p <- x/n
  q_profile_evenness <- function(q){
    qDest <- FUN(p,q)
    #S <- sum(x>0)
    S <- sum(x>0)
    E1 <- ifelse(q!=1, (1-qDest^(1-q))/(1-S^(1-q)), log(qDest)/log(S))
    E2 <- ifelse(q!=1, (1-qDest^(q-1))/(1-S^(q-1)), log(qDest)/log(S))
    E3 <- (qDest-1)/(S-1)
    E4 <- (1-1/qDest)/(1-1/S)
    E5 <- log(qDest)/log(S)
    if(q==0){
      p <- p[p>0]
      nu <- abs(p - (1/S))
      nu <- nu[nu > 0]
      sub1 <- (sum(log(abs(nu)))/sum(nu>0)-(log(1-1/S)+(1-S)*log(S))/S)
      E6 <- 1-exp(sub1)
    }else{
      p <- p[p>0]
      E6 <- 1-(sum(abs(p-1/S)^q)/((1-1/S)^q+(S-1)*S^(-q)))^(1/q)
    }
    
    #E6 <- ifelse(q=1, 1-sum(abs(p-1/S)^(1-q))/(abs(1-1/S)^(1-q)+)
    return(c(E1,E2,E3,E4,E5,E6))
  }
  out <- as.matrix(t(sapply(q.order, q_profile_evenness)))
  colnames(out) <- c("E1", "E2", "E3", "E4", "E5", "E6")
  out
}



tax_q_profile <- function(x, name1){
  x <- as.data.frame(x)
  x$q.order <-as.character(q)
  x1 <- melt(x, id.vars = c("q.order"))
  names(x1) <- c("q", "Habitat", "evenness")
  ggplot(x1, aes(q, evenness))+
    geom_line(aes(color = Habitat, group = Habitat, linetype = Habitat), size = 1.1)+
    scale_linetype_manual(values=c("dashed", "1111", "solid"))+
    theme_bw()+
    theme(axis.title = element_text(size = 14), 
          axis.text = element_text(size = 12), 
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          legend.position = "bottom")+
    scale_x_discrete(breaks=seq(0, 2, 0.5))+
    theme(legend.key.width = unit(2,"cm"))+
    theme(plot.title = element_text(size=20, face="bold.italic",hjust = 0.5))+
    ggtitle(name1)+
    xlab("Diversity order q")
  #ylim(c(0, 1))
}
##########caculate evenness##########

data <- read.table("Alpine data.txt")
q <- seq(0, 2, 0.05)
name1 <- c("E1", "E2", "E3", "E4", "E5", "E6")
Ricotta_ind_evenness1 <- lapply(list(data[, 1], data[, 2], data[, 3]), function(x) new_fun(x = x, q.order = q))
Ricotta_ind_evenness2 <- array(unlist(Ricotta_ind_evenness1), c(length(q), 6, 3))
Ricotta_ind_evenness3 <- aperm(Ricotta_ind_evenness2, c(1, 3, 2))
dimnames(Ricotta_ind_evenness3)[[1]] <- paste0("q = ", q)
dimnames(Ricotta_ind_evenness3)[[2]] <- colnames(data)
dimnames(Ricotta_ind_evenness3)[[3]] <- name1

##########plot by value
ggarrange(tax_q_profile(Ricotta_ind_evenness3[, , 1], name1[1])+ylim(c(0.4, 1)),
          tax_q_profile(Ricotta_ind_evenness3[, , 2], name1[2])+ylim(c(0.4, 1)),
          tax_q_profile(Ricotta_ind_evenness3[, , 3], name1[3])+ylim(c(0.4, 1)),
          tax_q_profile(Ricotta_ind_evenness3[, , 4], name1[4])+ylim(c(0.4, 1)),
          tax_q_profile(Ricotta_ind_evenness3[, , 5], name1[5])+ylim(c(0.4, 1)),
          tax_q_profile(Ricotta_ind_evenness3[, , 6], name1[6])+ylim(c(0.4, 1)),
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom")








####################################################################################
#
# (2). Computing the contribution of each species/node (Figures 4 and 5)
#
####################################################################################

#'@param type is tax. or phy.
#'@param type2 is "species" or "k"."species" means each species how much contribution for dissimilarity (Jaccard-type dissimilarity or Sorensen-type dissimilarity).
#' "k" means each location how much contribution for dissimilarity (Jaccard-type dissimilarity or Sorensen-type dissimilarity).
#' @return dissimilarity (Jaccard-type dissimilarity or Sorensen-type dissimilarity)

dis1 <- function(x, q, type = "tax", type2 = "species", tree = NULL){
  if(type2 == "species"){
    FUN <- rowSums
  }else{
    FUN <- colSums
  }
  if(type == "tax"){
    x <- as.matrix(x)
    x <- x[rowSums(x)>0, ]
    N <- ncol(x)
    zbar <- rowSums(x)/N
    x1 <- x[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    if(q==0){
      UqN <- FUN(x==0)/((N-1)*(sum(rowSums(x)>0)))
      CqN <- FUN(x==0)/((N-1)*(sum(apply(x, 2, function(i){sum(i>0)}))))
    }else if(q==2){
      UqN <- FUN((x1-zbar1)^2)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1-zbar1)^2)/((1-N^(1-q))*sum(x1^q))
    }else if(q!=1){
      UqN <- FUN((x1)^q-(zbar1)^q)/((N^q-N)*sum(zbar1^q))
      CqN <- FUN((x1)^q-(zbar1)^q)/((1-N^(1-q))*sum(x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(x1*log(x2), na.rm = T)/((sum(x)*log(N)))
      CqN <- UqN
    }
  }else{
    Li <- c(tree1$leaves, tree1$nodes)
    cumtree = function(a, tree){
      a <- a[names(tree$leaves)]
      for(i in 1:length(tree$parts)){
        a[1+length(a)] <- sum(a[tree$parts[[i]]])
        names(a)[length(a)] <- names(tree$parts)[i]
      }
      a
    }
    ai <- apply(x, 2, cumtree, tree1)
    wt <- apply(ai, 1, function(x1)(sum(x1))^q/sum(Li*rowSums(ai, na.rm = T)^q))
    N <- ncol(ai)
    zbar <- rowSums(ai)/N
    x1 <- ai[zbar>0, ]
    zbar1 <- zbar[zbar>0]
    Li <- Li[zbar>0]
    T1 <- sum(rowSums(x1)*Li)
    if(q==0){
      if(type2 == "species"){
        rn <- nrow(x1)
        UqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))})/((N-1)*sum(Li)) 
        CqN <- sapply(1:rn, function(i){(Li[i]*sum(x1[i, ]==0))/((N-1)*sum(Li*rowSums(x1!=0)))})
      }else{
        UqN <- apply(x1, 2, function(x){sum(Li[x==0])})/((N-1)*sum(Li)) 
        CqN <- apply(x1, 2, function(x){sum(Li[x==0])/((N-1)*sum(Li*colSums(x1!=0)))})
      }
      
    }else if(q==2){
      UqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1-zbar1)^2), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else if(q!=1){
      UqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((N^q-N)*sum(Li*zbar1^q))
      CqN <- FUN(Li*((x1)^q-(zbar1)^q), na.rm = T)/((1-N^(1-q))*sum(Li*x1^q))
    }else{
      x2 <- x1/zbar1
      UqN <- FUN(Li*x1*log(x2), na.rm = T)/(T1*log(N))
      CqN <- UqN
    }
  }
  
  # c(sum(UqN), sum(CqN))
  rbind(UqN, CqN)
}
#' plot each species how much contribution for dissimilarity (Jaccard-type dissimilarity or Sorensen-type dissimilarity).
#' @param  data is the value of output.
#' @param title_name is the title name of plot. 

draw_dis_spe <- function(data, title_name, type = "tax"){
  colnames(data) <- c("q = 0", "q = 1", "q = 2")
  data <- melt(data)
  g <- ggplot(data, aes(x = as.factor(Var1), y = value, fill = Var2))+
    geom_col(width = 0.2)+
    facet_grid(Var2~., scales = "free_y")+
    theme_bw()+
    # ylim(c(0, max(data[, 3])))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3), 
          axis.title = element_text(size = 14), 
          plot.title = element_text(hjust = 0.5))+
    guides(fill=FALSE)+
    ggtitle(title_name)
  
  if(type == "tax"){
    g <- g +
      xlab("Species")+
      ylab("Species contribution")
  }else{
    g <-  g +
      xlab("Species/node")+
      ylab("Species/node contribution")
  }
  return(g)
}

######plot tree####
plot.phylog <- function (x, y = NULL,
                         f.phylog = 0.5, cleaves = 1, cnodes = 0,
                         labels.leaves = names(x$leaves), clabel.leaves = 1,
                         labels.nodes = names(x$nodes), clabel.nodes = 0,
                         sub = "", csub = 1.25, possub = "bottomleft", draw.box = FALSE, ...)
{
  if (!inherits(x, "phylog")) 
    stop("Non convenient data")
  leaves.number <- length(x$leaves)
  leaves.names <- names(x$leaves)
  nodes.number <- length(x$nodes)
  nodes.names <- names(x$nodes)
  if (length(labels.leaves) != leaves.number) labels.leaves <- names(x$leaves)
  if (length(labels.nodes) != nodes.number) labels.nodes <- names(x$nodes)
  leaves.car <- gsub("[_]"," ",labels.leaves)
  nodes.car <- gsub("[_]"," ",labels.nodes)
  mar.old <- par("mar")
  on.exit(par(mar=mar.old))
  
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  
  if (f.phylog < 0.05) f.phylog <- 0.05 
  if (f.phylog > 0.95) f.phylog <- 0.95 
  
  maxx <- max(x$droot)
  plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
               yaxt = "n", xlim = c(-maxx*0.15, maxx/f.phylog), ylim = c(-0.05, 1), xaxs = "i", 
               yaxs = "i", frame.plot = FALSE)
  
  x.leaves <- x$droot[leaves.names]
  x.nodes <- x$droot[nodes.names]
  if (is.null(y)) y <- (leaves.number:1)/(leaves.number + 1)
  else y <- (leaves.number+1-y)/(leaves.number+1)
  names(y) <- leaves.names
  xcar <- maxx*1.05
  xx <- c(x.leaves, x.nodes)
  
  if (clabel.leaves > 0) {
    for (i in 1:leaves.number) {
      text(xcar, y[i], leaves.car[i], adj = 0, cex = par("cex") * 
             clabel.leaves)
      #segments(xcar, y[i], xx[i], y[i], col = grey(0.7))
    }
  }
  yleaves <- y[1:leaves.number]
  xleaves <- xx[1:leaves.number]
  if (cleaves > 0) {
    for (i in 1:leaves.number) {
      points(xx[i], y[i], pch = 21, bg=1, cex = par("cex") * cleaves)
    }
  }
  yn <- rep(0, nodes.number)
  names(yn) <- nodes.names
  y <- c(y, yn)
  for (i in 1:length(x$parts)) {
    w <- x$parts[[i]]
    but <- names(x$parts)[i]
    y[but] <- mean(y[w])
    b <- range(y[w])
    segments(xx[but], b[1], xx[but], b[2])
    x1 <- xx[w]
    y1 <- y[w]
    x2 <- rep(xx[but], length(w))
    segments(x1, y1, x2, y1)
  }
  if (cnodes > 0) {
    for (i in nodes.names) {
      points(xx[i], y[i], pch = 21, bg="white", cex = cnodes)
      
    }
  }
  #if (clabel.nodes > 0) {
  text(xx[names(x.nodes)], y[names(x.nodes)], nodes.car, 
       clabel.nodes, col = 2 , pos = 4)
  #}
  x <- (x.leaves - par("usr")[1])/(par("usr")[2]-par("usr")[1])
  y <- y[leaves.names]
  xbase <- (xcar - par("usr")[1])/(par("usr")[2]-par("usr")[1])
  if (csub>0) scatterutil.sub(sub, csub=csub, possub=possub)
  if (draw.box) box()
  if (cleaves > 0) points(xleaves, yleaves, pch = 21, bg=1, cex = par("cex") * cleaves)
  
  return(invisible(list(xy=data.frame(x=x, y=y), xbase= xbase, cleaves=cleaves)))
}
######arrange data######
data1 <- read.table("Alpine data.txt")
tree <- read.table("Alpine phylo_tree.txt", header = F)[1,1]
tree1 <- newick2phylog(tree)
plot.phylog(tree1,
            draw.box = TRUE, labels.nodes = names(tree1$nodes), clabel.leaves = 1, clabel.nodes = 1)


######caculate for tax. dissimarity and then plot.#####
t01 <- t(dis1(data1, 0, type = "tax", type2 = "species"))
t11 <- t(dis1(data1, 1, type = "tax", type2 = "species"))
t21 <- t(dis1(data1, 2, type = "tax", type2 = "species"))

tax_UqN_r <- cbind(t01[, 1], t11[, 1], t21[, 1])
tax_CqN_r <- cbind(t01[, 2], t11[, 2], t21[, 2])
draw_dis_spe(tax_UqN_r, "Jaccard-type taxonomic dissimilarity")
draw_dis_spe(tax_CqN_r, "Sorensen-type taxonomic dissimilarity")

######caculate for phy. dissimarity and plot#####

p01 <- t(dis1(data1, 0, type = "phy", type2 = "species", tree = tree1))
p11 <- t(dis1(data1, 1, type = "phy", type2 = "species", tree = tree1))
p21 <- t(dis1(data1, 2, type = "phy", type2 = "species", tree = tree1))

phy_UqN_r <- cbind(p01[, 1], p11[, 1], p21[, 1])
phy_CqN_r <- cbind(p01[, 2], p11[, 2], p21[, 2])

draw_dis_spe(phy_UqN_r, "Jaccard-type phylogenetic dissimilarity", type = "phy")
draw_dis_spe(phy_CqN_r, "Sorensen-type phylogenetic dissimilarity", type = "phy")

##################################################


