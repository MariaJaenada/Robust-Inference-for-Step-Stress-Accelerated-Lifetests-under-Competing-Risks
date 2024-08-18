rm(list=ls())

source("dpd_estimation_competing.R")

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
cont = 0.2
beta.list = c(0.2,0.4,0.6,0.8,1)
initial = c(7, -0.01, 4, -0.01)
 
n = simulate.sample(theta, cont, N, tauvec, stressvec, ITvec,  seed = 2024, C=C)
n

estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial = initial, C=C, beta.list = beta.list)
estimators

