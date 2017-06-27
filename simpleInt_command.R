#dat1 to dat8 provided in the ec016.RData are 8 sequences shown in the paper of 
#Haonan. The unit is second.
#D is the length of domain of these sequences, u_max is the length of domain we want
#to estimate the covariance function, r_est are the grids for calculation.
D=5453
u_max=2
r_est=1:3000/3000*u_max

#pair the data and calculate the distance between them
distance18=pair_dat_cross(dat1,dat8,u_max)
distance28=pair_dat_cross(dat2,dat8,u_max)
distance38=pair_dat_cross(dat3,dat8,u_max)
distance48=pair_dat_cross(dat4,dat8,u_max)
distance58=pair_dat_cross(dat5,dat8,u_max)
distance68=pair_dat_cross(dat6,dat8,u_max)
distance78=pair_dat_cross(dat7,dat8,u_max)
distance88=pair_dat_modi(dat8,u_max)

#here the h is selected randomly. A CV function is given in the file of function,
#but right now it is either too time-consuming or not reliable. The calculation of 
#aa18-aa88 is quite time-consuming(about 20mins for each of them), and the result 
#has already been provided in the data set
h=0.01
# aa18=isotro_corr_est(distance18,r_est,h,u_max,D)
# aa28=isotro_corr_est(distance28,r_est,h,u_max,D)
# aa38=isotro_corr_est(distance38,r_est,h,u_max,D)
# aa48=isotro_corr_est(distance48,r_est,h,u_max,D)
# aa58=isotro_corr_est(distance58,r_est,h,u_max,D)
# aa68=isotro_corr_est(distance68,r_est,h,u_max,D)
# aa78=isotro_corr_est(distance78,r_est,h,u_max,D)
# aa88=isotro_corr_est(distance88,r_est,h,u_max,D)

#A quick view of the second order intensity
plot(aa18-length(dat1)*1.0*length(dat8)/D/D~r_est,cex=0.1)
plot(aa28-length(dat2)*1.0*length(dat8)/D/D~r_est,cex=0.1)
plot(aa38-length(dat3)*1.0*length(dat8)/D/D~r_est,cex=0.1)
plot(aa48-length(dat4)*1.0*length(dat8)/D/D~r_est,cex=0.1)
plot(aa58-length(dat5)*1.0*length(dat8)/D/D~r_est,cex=0.1)
plot(aa68-length(dat6)*1.0*length(dat8)/D/D~r_est,cex=0.1)
plot(aa78-length(dat7)*1.0*length(dat8)/D/D~r_est,cex=0.1)
plot(aa88-length(dat8)*1.0*length(dat8)/D/D~r_est,cex=0.1)

#Calculate the matrix
res_vec=aa88[1:2400]
X=matrix(0,2400,600)
for(i in 1:2400){
  X[i,]=i-(-600):(-1)
}
X18=matrix(aa18[X],2400,600)
X28=matrix(aa28[X],2400,600)
X38=matrix(aa38[X],2400,600)
X48=matrix(aa48[X],2400,600)
X58=matrix(aa58[X],2400,600)
X68=matrix(aa68[X],2400,600)
X78=matrix(aa78[X],2400,600)
X88=matrix(aa88[X],2400,600)



library(MASS)  # Package needed to generate correlated precictors
library(glmnet)


#code with spline
library("splines")
nknots=30
ndf=nknots+4
base=bs(x=(1:600/600*0.4)[1:600],knots=c(0:nknots)/nknots*0.4, degree = 3)
X_tilde_18=X18%*%base
X_tilde_28=X28%*%base
X_tilde_38=X38%*%base
X_tilde_48=X48%*%base
X_tilde_58=X58%*%base
X_tilde_68=X68%*%base
X_tilde_78=X78%*%base
X_tilde_88=X88%*%base

#spline with 7 input and itself

X=cbind(X_tilde_18,X_tilde_28,X_tilde_38,X_tilde_48,X_tilde_58,X_tilde_68,X_tilde_78,X_tilde_88)
fit.lasso.spline=glmnet(X, res_vec, family="gaussian", alpha=1,intercept=T,lambda=0.1)
coef.lasso.spline=coef(fit.lasso.spline)[-c(1)]
coef.lasso.1=coef.lasso.spline[1:ndf]
coef.lasso.2=coef.lasso.spline[(ndf+1):(2*ndf)]
coef.lasso.3=coef.lasso.spline[(2*ndf+1):(3*ndf)]
coef.lasso.4=coef.lasso.spline[(3*ndf+1):(4*ndf)]
coef.lasso.5=coef.lasso.spline[(4*ndf+1):(5*ndf)]
coef.lasso.6=coef.lasso.spline[(5*ndf+1):(6*ndf)]
coef.lasso.7=coef.lasso.spline[(6*ndf+1):(7*ndf)]
coef.lasso.8=coef.lasso.spline[(7*ndf+1):(8*ndf)]
#quick view of different K_i function. Remark: some of them are very close to 0, like coef.lasso.4
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.1)),cex=0.1)
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.2)),cex=0.1)
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.3)),cex=0.1)
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.4)),cex=0.1)
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.5)),cex=0.1)
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.6)),cex=0.1)
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.7)),cex=0.1)
plot(cbind((600:1/600*0.4)[1:600],base%*%as.vector(coef.lasso.8)),cex=0.1)
abline(h=0)
plot(cbind(1,X)%*%coef(fit.lasso.spline),cex=0.1)
points(res_vec,cex=0.1,col=2)
