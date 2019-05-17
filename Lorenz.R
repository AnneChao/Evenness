#*******************************************************************************
#********************************************************************************
# 
#
## R scipts "Lorenz" for Chao and Ricotta (2019) evenness paper. 
## The R code is for plotting Figure 6 of Chao and Ricotta's evenness paper.
# NOTE: The packages "ggplot2" must be installed and loaded before running the scripts. 
# 
#
#*******************************************************************************
#*******************************************************************************






#====Function for Lorenz Curve====
#'Run the two functions:"lcs_table" and "lcs_plot" first
#'lcs_table is use to compute the values of each curve
#' @param ps is an observed species-by-assemblage frequency matrix. 
#' @return is the values of each curve.

require(ggplot2)
lcs_table = function(ps){
  lc = function(p,j){
    p = p[p>0]
    s = length(p)
    p = p/sum(p)
    p = sort(p)
    out_x = sapply(1:s, function(i){i/s})
    out_y = sapply(1:s, function(i){sum(p[1:i])})
    out = data.frame(x = c(0,out_x),y = c(0,out_y),source = names(ps)[j])
  }
  curves = lapply(1:length(ps),function(j){lc(ps[[j]],j)})
  curves = do.call(rbind,curves)
}

#'lcs_plot is used to plot each curve
#' @param curves is output of function lcs_table.
#' @return is the plot for Lorenz curve.

lcs_plot = function(curves){
  tmp = data.frame(x=c(0,1),y = c(0,1), source = "Completely even")
  curves = rbind(tmp,curves)
  curves$source = factor(curves$source, levels = unique(curves$source))
  g = ggplot(data = curves,aes(x = x,y = y,color = source, linetype = source))
  g = g + geom_line(size = 1.2)+
    geom_point(size = 3.5)+
    theme_bw()+
    theme(legend.position="bottom")+
    labs(col = "")+
    guides(linetype = FALSE)+
    ggtitle("Lorenz curve")+
    xlab("Cumulative proportion of species")+
    ylab("Cumulative proportion of abundance")+
    theme(plot.title = element_text(size=20, face="bold.italic",hjust = 0.5))
  g
}


#====Example====
#the input must be a list with each element's name
com_1 = c(10, 2)
com_2 = c(10, 10, 2, 1, 1)
ps = list(com_1=com_1, com_2 = com_2)
result_table = lcs_table(ps)
lcs_plot(result_table)
