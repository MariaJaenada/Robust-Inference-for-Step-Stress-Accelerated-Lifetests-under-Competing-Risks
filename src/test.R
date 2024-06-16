rm(list=ls())

N = 2*180
stressvec = c(35,45)

ITvec=c(15,25,35,45,55,65,75)
tauvec = c(45,75)
C=2

x0 = 30
t0 = 50
alpha = 0.05
alpha0= 0.05

theta = c(5,-0.02,6.2,-0.04)
theta = -c(log(0.004),0.02, log(0.003),0.04)

#theta.cont.vec = list(theta, theta.cont, theta.cont2, theta.cont3)
cont.vec = c(0,0.05,0.1,0.15,0.2,0.25,0.3) #0.4,0.6)
beta.list = c(0.2,0.4,0.6,0.8,1)

source("G:/Mi unidad/UCM/ARTICULOS/One shot devices/SSALT Competing risks/codes/dpd_estimation_competing.R")
B=500

plot.error <- function(res.matrix, epsilons, xlabel, ylabel="MSE", tl = "", beta.list = c(0.2,0.4,0.6,0.8,1)){
  
  
  x11()
  colores= c("#0D0D0D","#3B3B3B", "#1874CD",  "#87CEFA", "#008B00", "#FF7F00", "#EE7942", "#CD2990", "#EE7AE9", "#CD3333")
  formas=c(17,3,15,5,2,7,1,8,9,10,11,12)
  
  plot(epsilons, res.matrix[,1],col=colores[1],"o",lty=1, pch=formas[1],lwd=1.6, cex=1.05, ylim=c(min(res.matrix),max(res.matrix)), 
       xlab = xlabel, ylab = ylabel)
  for (j in 2:(length(beta.list)+1)){
    lines(epsilons ,res.matrix[,j],col=colores[j],"o",lty=j, pch=formas[j],lwd=1.6,cex=1.05)
  }
  legend("topleft", title= expression(beta), inset=c(0,0.0), c("MLE",as.character(c( 0.2, 0.4, 0.6, 0.8,1))), col = colores, lty = c(1,4:8),pch= formas, cex =0.75 )
  title(main = tl)
}

plot.mean <- function(res, epsilons, xlabel, ylabel="MSE", tl = "", beta.list = c(0.2,0.4,0.6,0.8,1)){
  
  
  x11()
  colores= c("#0D0D0D","#3B3B3B", "#1874CD",  "#87CEFA", "#008B00", "#FF7F00", "#EE7942", "#CD2990", "#EE7AE9", "#CD3333")
  formas=c(17,3,15,5,2,7,1,8,9,10,11,12)
  
  plot(epsilons, res[[1]][,1],col=colores[1],"o",lty=1, pch=formas[1],lwd=1.6, cex=1.05,
       xlab = xlabel, ylab = ylabel, ylim = c(min(unlist(res)), max(unlist(res))))
  for (j in 2:(length(beta.list)+1)){
    lines(epsilons ,res[[j]][,1],col=colores[j],"o",lty=j, pch=formas[j],lwd=1.6,cex=1.05)
  }
  legend("topleft", title= expression(beta), inset=c(0,0.0), c("MLE",as.character(c( 0.2, 0.4, 0.6, 0.8,1)), "Opt"), col = colores, lty = c(1,4:8),pch= formas, cex =0.75 )
  title(main = tl)
}

plot.sesitivity <- function(theta, tauvec, stressvec, ITvec, C){
 
  betas = seq(0,1,by=0.01)
  sensitivity.matrix = matrix(NA, nrow = length(betas), ncol = 2*C)
  
  for(b in 1:length(betas)){
    sensitivity.matrix[b,] = sensitivity.points(theta, tauvec, stressvec, ITvec, C, betas[b])
  }
  tl = c(expression(widehat(a)["01"]^beta), expression(widehat(a)["11"]^beta), expression(widehat(a)["02"]^beta), expression(widehat(a)["12"]^beta))
 
   #plot
  #colores= c("#0D0D0D", "#1874CD",  "#87CEFA", "#008B00", "#FF7F00", "#EE7942", "#CD2990", "#EE7AE9", "#CD3333")
  colores <- c( "#E69F00", "#009E73",  "#0072B2",  "#CC79A7")
  for (k in 1:(2*C)){
    #plot(betas, sensitivity.matrix[,k], col=colores[k], type = "l", lty=1,lwd=1.6, xlab = expression(beta), ylab = "IF")
    x11()
    ggplot(data=data.frame("beta" = betas, "IF" = sensitivity.matrix[,k]), aes(x=beta, y=IF))+
    geom_line(size=1, color= colores[k])+
    labs(title=" ",x=expression(beta), y = "sefl-standarized sensitivity")+
    scale_color_brewer(palette="Blues")+
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background=element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray"))
    }
  
}

# theta = -c(log(0.004),0.04, log(0.003),0.03)
#plot.sesitivity(theta, tauvec, stressvec, ITvec, C)


# for (s in 3:15){
#   step.cont <<- s
#   res.matrix = simulate(theta, cont.vec,  N, tauvec, stressvec, ITvec, C=2, B=10)
#   plot.error(res.matrix, cont.vec, xlabel="cont", ylabel="MSE", tl = toString(s))
#   
# }


res.matrix = simulate(theta, cont.vec,  N, tauvec, stressvec, ITvec, C=2, B=B, beta.list = beta.list)
plot.error(res.matrix, cont.vec, xlabel="cont", ylabel="MSE", beta.list = beta.list)

#MEAN LIFETIME
res.asymptotic = simulate.charasteristic(theta,  cont.vec,  N, tauvec, stressvec, ITvec,
                              x0=30, t0=50, alpha = 0.05,
                              C=2, B=B, beta.list = c(0.2,0.4,0.6,0.8,1),
                              charasteristic="mean.lifetime", alpha0=0.05)

plot.mean(res.asymptotic$direct, cont.vec, xlabel="e", ylabel="MSE", tl = "", beta.list = c(0.2,0.4,0.6,0.8,1))
saveRDS(res.asymptotic, "asymptotic_meanlifetime.rds")

# readRDS("asymptotic_trial_meanlifetime.rds")
res.bootstrap = simulate.bootstrap.charasteristic(theta,  cont.vec,  N, tauvec, stressvec, ITvec, x0=30, t0=50, alpha = 0.05,
                                                   C=2, R=1000, B=100, beta.list = c(0.2,0.4,0.6,0.8,1),
                                                   charasteristic="mean.lifetime", alpha0=0.05)
saveRDS(res.bootstrap, "bootstrap_meanlifetime.rds")
# readRDS(res.bootstrap, "bootstrap_trial_meanlifetime.rds")

##create tables summary
summary=list()
nbetas = length(beta.list)+1
for (cont in 1:length(cont.vec)){
  summary[[cont]] = matrix(NA, nrow = nbetas, ncol = 6)
  for (beta in 1:nbetas){
    summary[[cont]][beta,1] = res.asymptotic$direct[[beta]][cont,5]*100
    summary[[cont]][beta,2] = res.asymptotic$direct[[beta]][cont,6]
    
    summary[[cont]][beta,3] = res.asymptotic$transformed[[beta]][cont, 5]*100
    summary[[cont]][beta,4] = res.asymptotic$transformed[[beta]][cont, 6]
    
    summary[[cont]][beta,5] = res.bootstrap[[beta]][cont, 5]*100
    summary[[cont]][beta,6] = res.bootstrap[[beta]][cont, 6]
  }
}

library(xtable)
  for (cont in 1:length(cont.vec)){print(xtable(summary[[cont]]))}


#RELIABILITY
res.asymptotic = simulate.charasteristic(theta,  cont.vec,  N, tauvec, stressvec, ITvec,
                                         x0=30, t0=50, alpha = 0.05,
                                         C=2, B=B, beta.list = c(0.2,0.4,0.6,0.8,1),
                                         charasteristic="reliability", alpha0=0.05)

plot.mean(res.asymptotic$direct, cont.vec, xlabel="e", ylabel="MSE", tl = "", beta.list = c(0.2,0.4,0.6,0.8,1))

saveRDS(res.asymptotic, "asymptotic_reliability.rds")
# res.asymptotic = readRDS("results/asymptotic_trial_reliability.rds")

res.bootstrap = simulate.bootstrap.charasteristic(theta,  cont.vec,  N, tauvec, stressvec, ITvec, x0=30, t0=50, alpha0 = 0.5,
                                                  C=2, R=1000, B=100, beta.list = c(0.2,0.4,0.6,0.8,1),
                                                  charasteristic="reliability", alpha=0.05)
saveRDS(res.bootstrap, "bootstrap_reliability.rds")
#res.bootstrap = readRDS("results/bootstrap_trial_reliability.rds")

##create tables summary
summary=list()
nbetas = length(beta.list)+1
for (cont in 1:length(cont.vec)){
  summary[[cont]] = matrix(NA, nrow = nbetas, ncol = 6)
  for (beta in 1:nbetas){
    summary[[cont]][beta,1:2] = res.asymptotic$direct[[beta]][cont,5:6]
    summary[[cont]][beta,3:4] = res.asymptotic$transformed[[beta]][cont, 5:6]
    summary[[cont]][beta,5:6] = res.bootstrap[[beta]][cont, 5:6]
  }
}

library(xtable)
for (cont in 1:length(cont.vec)){print(xtable(summary[[cont]]))}



#MEDIAN
res.asymptotic = simulate.charasteristic(theta,  cont.vec,  N, tauvec, stressvec, ITvec,
                                         x0=30, t0=50, alpha0 = 0.5,
                                         C=2, B=B, beta.list = c(0.2,0.4,0.6,0.8,1),
                                         charasteristic="quantile", alpha=0.05)

plot.mean(res.asymptotic$direct, cont.vec, xlabel="e", ylabel="MSE", tl = "", beta.list = c(0.2,0.4,0.6,0.8,1))

saveRDS(res.asymptotic, "asymptotic_quantile.rds")
# res.asymptotic = readRDS("results/asymptotic_trial_quantile.rds")

res.bootstrap = simulate.bootstrap.charasteristic(theta,  cont.vec,  N, tauvec, stressvec, ITvec, x0=30, t0=50, alpha0 = 0.5,
                                                  C=2, R=1000, B=100, beta.list = c(0.2,0.4,0.6,0.8,1),
                                                  charasteristic="quantile", alpha=0.05)
saveRDS(res.bootstrap, "bootstrap_quantile.rds")
# res.bootstrap = readRDS("results/bootstrap_trial_quantile.rds")

##create tables summary
summary=list()
nbetas = length(beta.list)+1
for (cont in 1:length(cont.vec)){
  summary[[cont]] = matrix(NA, nrow = nbetas, ncol = 6)
  for (beta in 1:nbetas){
    summary[[cont]][beta,1:2] = res.asymptotic$direct[[beta]][cont,5:6]
    summary[[cont]][beta,3:4] = res.asymptotic$transformed[[beta]][cont, 5:6]
    summary[[cont]][beta,5:6] = res.bootstrap[[beta]][cont, 5:6]
  }
}

library(xtable)
for (cont in 1:length(cont.vec)){print(xtable(summary[[cont]]))}