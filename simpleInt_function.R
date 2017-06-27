library("locpol")
##calculate the k_th moments of standard Epanechnikov kernel on the intervel (s-h,s+h)\cap [0,u_max](for isotropic situation)
Int_general=function(s,h,k,u_max){
  if(s<=u_max-h&s>=h) return(ifelse(k%%2==0,3*h^(k+1)/(k+1)/(k+3),0))
  if(s<h) return(3/4*(2*h^(k+1)/(k+1)/(k+3)+(-s)^(k+3)/(k+3)/h^2-(-s)^(k+1)/(k+1)))
  s=u_max-s
  if(k%%2==0) return(3/4*(2*h^(k+1)/(k+1)/(k+3)+(-s)^(k+3)/(k+3)/h^2-(-s)^(k+1)/(k+1)))
  return(-3/4*(2*h^(k+1)/(k+1)/(k+3)+(-s)^(k+3)/(k+3)/h^2-(-s)^(k+1)/(k+1)))
}

#Estimation of second-order intensity for isotropic situation(For estimate at a specific position est_r)
isotro_corr_est_CV=function(distance,est_r,h,u_max,D){
  s=est_r
  I1=Int_general(s,h,0,u_max)
  I2=Int_general(s,h,1,u_max)
  I3=I2
  I4=Int_general(s,h,2,u_max)
  
  temp=distance-s
  temp_flag=temp<h&temp>-h
  temp=temp[temp_flag]
  temp1=EpaK(temp/h)
  A1=temp1/(D-distance[temp_flag])
  A2=sum(A1*temp)
  A1=sum(A1)
  
  temp=matrix(c(I1,I2,I3,I4),2,2)
  return(solve(temp,c(A1,A2))[1])
}

#Estimation of second-order intensity for isotropic situation(For general usage)
isotro_corr_est=function(distance,est_r,h,u_max,D){
  
  result=NULL
  for(i in 1:length(est_r)){
    print(i)
    result[i]=isotro_corr_est_CV(distance,est_r[i],h,u_max,D)
  }
  return(result)
}

#pair data point and calculate the distance
pair_dat_cross=function(dat1,dat2,u_max){
  flag1=rep(0,length(dat1))
  temp_record=sum(dat2<(dat1[1]-u_max))
  flag1[1]=temp_record
  for(i in 2:length(dat1)){
    temp=T
    temp1=temp_record+5
    while(temp){
      if(dat2[temp1]<(dat1[i]-u_max))
        temp1=temp1+5
      else
        temp=F
    }
    temp_record=sum(dat2[(flag1[i-1]+1):temp1]<(dat1[i]-u_max))+temp_record
    flag1[i]=temp_record
  }
  
  flag2=rep(0,length(dat1))
  temp_record=sum(dat2<(dat1[1]))
  flag2[1]=temp_record
  for(i in 2:length(dat1)){
    temp=T
    temp1=temp_record+5
    while(temp){
      if(temp1>=length(dat2)){
        temp1=length(dat2)
        break
      }
      if(dat2[temp1]<(dat1[i]))
        temp1=temp1+5
      else
        break
    }
    temp_record=sum(dat2[(flag2[i-1]+1):temp1]<(dat1[i]))+temp_record
    if(temp_record==length(dat2)){
      flag2[i:length(flag2)]=length(dat2)
      break
    }
    flag2[i]=temp_record
  }
  
  record=flag2-flag1
  result=rep(0,sum(record))
  flag=1
  for(i in 2:length(dat1)){
    if(record[i]==0) next
    result[flag:(flag+record[i]-1)]=dat1[i]-dat2[(flag1[i]+1):(flag2[i])]
    flag=flag+record[i]
  }
  return(result)
} 

pair_dat_modi=function(dat,u_max){
  flag=1
  record=rep(0,length(dat))
  for(i in 2:length(dat)){
    record[i]=sum((dat[i]-dat[flag:(i-1)])<u_max)
    flag=flag+i-flag-record[i]
  }
  result=rep(0,sum(record))
  flag=1
  for(i in 2:length(dat)){
    if(record[i]==0) next
    result[flag:(flag+record[i]-1)]=dat[i]-dat[(i-record[i]):(i-1)]
    flag=flag+record[i]
  }
  return(result)
} 

#calculation of the cross-validation

CV_LS_2_cross=function(dat1,dat2,h,est_r,test_sample=1000){
  paired=pair_dat_cross(dat1,dat2)
  distance=abs(paired[,1]-paired[,2])
  distance=distance[distance<max(est_r)]
  sample_leave=sample(1:length(distance),size=test_sample)
  LS=NULL
  for(i in 1:test_sample){
    dat_test=distance[sample_leave[i]]
    dat_train=distance[-sample_leave[i]]
    est_test=isotro_corr_est_CV(dat_train,dat_test,h)
    LS[i]=((est_test))
  }
  all=isotro_corr_est(distance,est_r,h=h)
  return(2*mean(LS)*length(distance)-1/length(r_est)*sum(all^2)*max(r_est))
}

