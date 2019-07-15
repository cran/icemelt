

.packageName <- "icemelt"

###############   IM function implements the imputation method   #################

im= function(datamat,wmat,rfix,gridlen,ntimp,nximp)
{

  n= nrow(datamat)  
  m= ntimp             
  result= NULL

  delta_temp= datamat[,3]
  k= sum(delta_temp)
  
  left= datamat[,1]
  right1= datamat[,2]
  z2= datamat[,4]
  
  rcpos= which(delta_temp==0)
  notrcpos= (1:n)[-rcpos]
 
  liequalri= which(right1[notrcpos]==left[notrcpos])
  liminusrismall= which(right1[notrcpos]-left[notrcpos]<= gridlen)
  toosmall= notrcpos[liminusrismall]
  ltoo= length(toosmall)


  tmat= NULL
  numgrid= NULL
  csg= 0
  cumsumgrid= NULL
  for(i5 in 1:k)
  {
    i6= notrcpos[i5]

    newseq= seq(left[i6], right1[i6], gridlen)
      
    len_newseq= length(newseq)
    csg= csg+len_newseq
    tmat= c(tmat, newseq)
    numgrid= c(numgrid, len_newseq)
    cumsumgrid= c(cumsumgrid, csg)
  }

  maxgrid= max(numgrid)
  totngrid= cumsumgrid[k]

  lessthan1= liminusrismall[-match(liequalri,liminusrismall)]
  tmat[cumsumgrid[lessthan1]]= 0.5*(left[notrcpos[lessthan1]]+right1[notrcpos[lessthan1]])

  ttilmat= tmat 
  wbar=apply(wmat, 1, mean)
    
  r= 10000
  isd= sample(1:90000,1)  


####################### Turnbull codes ###########################

  r6_left= round(left,3)
  r6_right= round(right1,3)
  turn_tau= sort(unique(c(r6_left,r6_right)))
  len_tt= length(turn_tau)
  survobj= survfit(Surv(turn_tau[1:(len_tt-1)],rep(1,len_tt-1))~1)
  init_s= c(1, survobj$surv)
  turn_p= -diff(init_s)

  turn_tau12= cbind(turn_tau[-len_tt],turn_tau[-1])
  alpha_fn= function(xx,infi,supr) ifelse(xx[1]>=infi & xx[2]<=supr, 1, 0)
  alpha_mat= apply(turn_tau12,1,alpha_fn,infi=left,supr=right1)
  id_zero= which(apply(alpha_mat==0,1,all))
  if(length(id_zero)>0) alpha_mat= alpha_mat[-id_zero,]
  
  new_s= init_s
  
  turn_err= 1
  while(turn_err>0.01)
  {
    old_s= new_s
    turn_p= -diff(old_s)
    dj_temp= t(t(alpha_mat)*turn_p)/as.vector(alpha_mat%*%turn_p)
    turn_dj= apply(dj_temp,2,sum)
    turn_yj= rev(cumsum(rev(turn_dj)))
    new_s= c(1,cumprod(1-turn_dj/turn_yj))
    new_s[which(is.na(new_s)==T)]= 0
    turn_err= sum(abs((new_s[-len_tt]-old_s[-len_tt])))
  }  
  
  surv_left= new_s[match(r6_left,turn_tau)]
  surv_right= new_s[match(r6_right,turn_tau)]


  surv_left= log(surv_left)
  surv_right= log(1+surv_right)

  lrcpos= length(rcpos)
  slsl= surv_left*surv_left
  srsr= surv_right*surv_right

  des_mat=cbind(1, z2, surv_left, surv_right, z2*surv_left, z2*surv_right, surv_left*surv_right)
  
  out100=lm(wbar~ des_mat-1)
  my_gamma_init= as.vector(out100$coef)
  my_gamma_init= c(my_gamma_init)


  n_el=length(my_gamma_init)
  nn= n_el*(n_el+1)/2
  nn_el= n*n_el

  
  ncolp=ncol(des_mat)   
  storage.mode(des_mat)<-"double"
  capr=20000
  store_sigma2u=as.double(rep(0, capr))
  store_myx=matrix(0, nrow=1000, ncol=n)
  storage.mode(store_myx)<-"double"


  a_sigmau= a_sigmax=1
  b_sigmau= b_sigmax=1

  theta_rep= nximp
  mtr= m*theta_rep

  aic_vec= rep(0,4)
  aic_x= NULL
  sumsqwdiff= sum((wmat-wbar)^2)
  aic_part1= sumsqwdiff


  for(j1 in 1:4)
  {
    ncomp=j1
    store_sigma2x=matrix(0, nrow=capr, ncol=ncomp)
    store_gamma=matrix(0, nrow=capr, ncol=(ncomp*ncolp))
    store_vec=matrix(0, nrow=capr, ncol=ncomp)
    store_pi=matrix(0, nrow=capr, ncol=ncomp)

    storage.mode(store_gamma)<-"double"
    storage.mode(store_sigma2x)<-"double"
    storage.mode(store_vec)<-"double"
    storage.mode(store_pi)<-"double"

    prmn=matrix(rnorm(ncomp*ncolp), nrow=ncomp, ncol=ncolp)
    prmn[, 1]=quantile(wbar,   prob=1/(ncomp+1)+(0:(ncomp-1))/(ncomp+1))
    prvar=matrix(5, nrow=ncomp, ncol=ncolp)
    storage.mode(prmn)<-"double"
    storage.mode(prvar)<-"double"

    outxk=.Fortran("mh_parkcmp", 
    des_mat, isd=as.integer(runif(1, 10, 100000)), mm=as.integer(ncol(wmat)), 
    n=as.integer(n), ncolp=as.integer(ncolp), ncomp=as.integer(ncomp), 
    ndim=as.integer(ncomp*ncolp),  capr=as.integer(capr), wbar=as.double(wbar),
    sumsqwdiff, a_sigmau=as.double(a_sigmau), a_sigmax=as.double(a_sigmax),
    b_sigmau=as.double(b_sigmau), b_sigmax=as.double(b_sigmax), prmn,
    prvar, output1=store_gamma, output2=store_myx, output3=store_pi, 
    output4=store_sigma2u, output5=store_sigma2x,  output6=store_vec)


    if(ncomp==1)
    {
      post_gamma= apply(outxk$output1[(0.5*capr):capr,],2,mean)
      post_s2x= mean(outxk$output5[(0.5*capr):capr,])
      post_pi= mean(outxk$output3[(0.5*capr):capr,])
      post_s2u= mean(outxk$output4[(0.5*capr):capr])
    }
    else
    {
      post_gamma= apply(outxk$output1[(0.5*capr):capr,],2,mean)
      post_s2x= apply(outxk$output5[(0.5*capr):capr,],2,mean)
      post_pi= apply(outxk$output3[(0.5*capr):capr,],2,mean)
      post_s2u= mean(outxk$output4[(0.5*capr):capr])
    }
    
    aic_part2= 0
    for(j2 in 1:ncomp)
    {
      aic_part3= -0.5*(log(post_s2x[j2]+0.5*post_s2u)-log(0.5*post_s2x[j2]*post_s2u))+
      log(post_pi[j2])-0.5*log(post_s2x[j2])
      aic_part4= (wbar-des_mat%*%post_gamma[((j2-1)*ncolp+1):(j2*ncolp)])^2
      aic_part2= aic_part2 + exp(-0.5*aic_part4/(post_s2x[j2]+0.5*post_s2u)+ aic_part3)
    }
    lik_w= -0.5*aic_part1/post_s2u - n*log(2*pi*post_s2u) +sum(log(aic_part2))
    aic_vec[j1]= 2*(ncolp*ncomp+ 2*ncomp+ ncomp-1 + 1)- 2*lik_w
    aic_x=cbind(aic_x, t(outxk$output2[seq(1000,1000-(theta_rep-1)*50,-50),]))

  }

  aic_sel= which(aic_vec==min(aic_vec))
  final_x= aic_x[,(((aic_sel-1)*theta_rep+1):(aic_sel*theta_rep))]

######################   X imputation step complete   ###########################

  tvec= sort(unique(c(tmat,left)))
  ltvec= length(tvec)
  
  orderforhold= match(ttilmat,tvec)
  
  beta_old= c(-1.5,1.5)
  hallleft= log(left)
  hallold= log(tmat)

  beta_neww= rep(0,2)
  mysd= rep(0,2)

  storage.mode(final_x)= "double"
  
  allbeta_hope= matrix(0,2,mtr)
  storage.mode(allbeta_hope)= "double"

  if(rfix==0)
  {
    out_im= .Fortran("im_r0", beta_neww= as.double(beta_neww), as.double(beta_old), 
                       as.integer(cumsumgrid), as.integer(delta_temp), final_x, as.double(hallleft),
                       as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                       as.integer(lrcpos), as.integer(ltvec), as.integer(m), as.integer(maxgrid), 
                       as.integer(mtr), mysd= as.double(mysd), as.integer(n), as.integer(notrcpos), 
                       as.integer(numgrid), as.integer(orderforhold), as.integer(rcpos), 
                       as.double(rfix), as.integer(theta_rep), as.integer(totngrid), 
                       as.double(ttilmat), as.double(tvec), as.double(z2))
  }
  else
  {
    out_im= .Fortran("im_genr", beta_neww= as.double(beta_neww), as.double(beta_old),
                  as.integer(cumsumgrid), as.integer(delta_temp), final_x, as.double(hallleft),
                  as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                  as.integer(lrcpos), as.integer(ltvec), as.integer(m), as.integer(maxgrid),
                  as.integer(mtr), mysd= as.double(mysd), as.integer(n),
                  as.integer(notrcpos), as.integer(numgrid), as.integer(orderforhold),
                  as.integer(rcpos), as.double(rfix), as.integer(theta_rep),
                  as.integer(totngrid), as.double(ttilmat), as.double(tvec), as.double(z2))
  }
  
return_list= list(beta1.est= out_im$beta_neww[1], beta2.est= out_im$beta_neww[2], beta1.sd= out_im$mysd[1], beta2.sd= out_im$mysd[2])

}


#########################################################################################
####################   NV function implements the naive method   ########################

nv= function(datamat,wmat,rfix,gridlen,ntimp)
{

  n= nrow(datamat)
  m= ntimp
  result= NULL

  delta_temp= datamat[,3]
  k= sum(delta_temp)
  
  left= datamat[,1]
  right1= datamat[,2]
  z2= datamat[,4]
  
  rcpos= which(delta_temp==0)
  notrcpos= (1:n)[-rcpos]
  lrcpos= length(rcpos)
 
  liequalri= which(right1[notrcpos]==left[notrcpos])
  liminusrismall= which(right1[notrcpos]-left[notrcpos]<= gridlen)
  toosmall= notrcpos[liminusrismall]
  ltoo= length(toosmall)

  tmat= NULL
  numgrid= NULL
  csg= 0
  cumsumgrid= NULL
  for(i5 in 1:k)
  {
    i6= notrcpos[i5]

    newseq= seq(left[i6], right1[i6], gridlen)
      
    len_newseq= length(newseq)
    csg= csg+len_newseq
    tmat= c(tmat, newseq)
    numgrid= c(numgrid, len_newseq)
    cumsumgrid= c(cumsumgrid, csg)
  }

  maxgrid= max(numgrid)
  totngrid= cumsumgrid[k]

  lessthan1= liminusrismall[-match(liequalri,liminusrismall)]
  tmat[cumsumgrid[lessthan1]]= 0.5*(left[notrcpos[lessthan1]]+right1[notrcpos[lessthan1]])

  ttilmat= tmat

  if(ncol(wmat)==1)
    wbar= wmat[,1]
  else
    wbar=apply(wmat, 1, mean)

  
  r= 10000
  isd= sample(1:90000,1)

  tvec= sort(unique(c(tmat,left)))
  ltvec= length(tvec)
  orderforhold= match(ttilmat,tvec)
  
  beta_old= c(-1.5,1.5)
  hallleft= log(left)
  hallold= log(tmat)

  beta_neww= rep(0,2)
  mysd= rep(0,2)

  final_x= wbar

  if(rfix==0)
  {
    out_nv= .Fortran("nv_r0", beta_neww= as.double(beta_neww), as.double(beta_old), 
                    as.integer(cumsumgrid), as.integer(delta_temp), as.double(final_x), 
                    as.double(hallleft), as.double(hallold), as.integer(isd), as.integer(k), 
                    as.double(left), as.integer(lrcpos), as.integer(ltvec), as.integer(m),
                    as.integer(maxgrid), mysd= as.double(mysd), as.integer(n), 
                    as.integer(notrcpos), as.integer(numgrid), as.integer(orderforhold), 
                    as.integer(rcpos), as.double(rfix), as.integer(totngrid), as.double(ttilmat),
                    as.double(tvec), as.double(z2))
  }
  else
  {
  out_nv= .Fortran("nv_genr", beta_neww= as.double(beta_neww), as.double(beta_old), 
                  as.integer(cumsumgrid), as.integer(delta_temp), as.double(final_x),
                  as.double(hallleft), as.double(hallold), as.integer(isd),
                  as.integer(k), as.double(left), as.integer(lrcpos),
                  as.integer(ltvec), as.integer(m), as.integer(maxgrid),
                  mysd= as.double(mysd), as.integer(n), as.integer(notrcpos),
                  as.integer(numgrid), as.integer(orderforhold), as.integer(rcpos),
                  as.double(rfix), as.integer(totngrid), as.double(ttilmat),
                  as.double(tvec), as.double(z2))
  
  }

return_list= list(beta1.est= out_nv$beta_neww[1], beta2.est= out_nv$beta_neww[2], beta1.sd= out_nv$mysd[1], beta2.sd= out_nv$mysd[2])

}


#########################################################################################
###########   RC function implements the regression calibration method   ################


rc= function(datamat,wmat,rfix,gridlen,ntimp)
{
  n= nrow(datamat)
  m= ntimp
  result= NULL

  delta_temp= datamat[,3]
  k= sum(delta_temp)
  nrep= ncol(wmat)

  left= datamat[,1]
  right1= datamat[,2]
  z2= datamat[,4]
  
  rcpos= which(delta_temp==0)
  notrcpos= (1:n)[-rcpos]
  lrcpos= length(rcpos)
 
  liequalri= which(right1[notrcpos]==left[notrcpos])
  liminusrismall= which(right1[notrcpos]-left[notrcpos]<= gridlen)
  toosmall= notrcpos[liminusrismall]
  ltoo= length(toosmall)


  tmat= NULL
  numgrid= NULL
  csg= 0
  cumsumgrid= NULL
  for(i5 in 1:k)
  {
    i6= notrcpos[i5]

    newseq= seq(left[i6], right1[i6], gridlen)
      
    len_newseq= length(newseq)
    csg= csg+len_newseq
    tmat= c(tmat, newseq)
    numgrid= c(numgrid, len_newseq)
    cumsumgrid= c(cumsumgrid, csg)
  }

  maxgrid= max(numgrid)
  totngrid= cumsumgrid[k]

  lessthan1= liminusrismall[-match(liequalri,liminusrismall)]
  tmat[cumsumgrid[lessthan1]]= 0.5*(left[notrcpos[lessthan1]]+right1[notrcpos[lessthan1]])

  ttilmat= tmat
  wbar=apply(wmat, 1, mean)

  r= 10000
  isd= sample(1:90000,1)


##############   regression calibration estimate of X    #################

  estsigmau2= sum((wmat-wbar)^2)/(n*(nrep-1))
  wibarminuswbar= wbar-mean(wbar)

  nu= (n-1)*nrep
  estsigmax2= (nrep*sum(wibarminuswbar^2) - (n-1)*estsigmau2)/nu
  estsigmaxz= nrep*sum(wibarminuswbar*(z2-mean(z2)))/nu
  estsigmaz2= var(z2)
  
  rcmat= matrix(c(estsigmax2+estsigmau2/nrep, estsigmaxz, estsigmaxz, estsigmaz2),2,2)

  xrc= as.vector(mean(wbar)+ c(estsigmax2, estsigmaxz)%*% solve(rcmat)%*% rbind(wibarminuswbar, z2-mean(z2)))

# Bootstrap to calculate sigma_alpha
  maxbt= 200   # number of bootstrap samples
  bootout= matrix(0,maxbt,3)
  for(bt in 1:maxbt)
  {
    boot= sample(1:n,n,replace=T)
    z2boot= z2[boot]
    wmatboot= wmat[boot,]
    wbarboot=apply(wmatboot, 1, mean)

    estsigmau2= sum((wmatboot-wbarboot)^2)/(n*(nrep-1))
    wibarminuswbar= wbarboot-mean(wbarboot)

    nu= (n-1)*nrep
    estsigmax2= (nrep*sum(wibarminuswbar^2) - (n-1)*estsigmau2)/nu
    estsigmaxz= nrep*sum(wibarminuswbar*(z2boot-mean(z2boot)))/nu
    estsigmaz2= var(z2boot)

    mydet= estsigmaz2*(estsigmax2 + estsigmau2/nrep) - estsigmaxz*estsigmaxz
    alpha1= (estsigmax2*estsigmaz2 - estsigmaxz*estsigmaxz)/mydet
    alpha2= estsigmaxz*estsigmau2/(nrep*mydet)
    
    alpha0= mean(wbarboot)*(1-alpha1)-alpha2*mean(z2boot)
    
    bootout[bt,]= c(alpha0,alpha1,alpha2)
  }

  sigma_alpha= cov(bootout)
  tvec= sort(unique(c(tmat,left)))
  ltvec= length(tvec)
  
  orderforhold= match(ttilmat,tvec)
  
  beta_old= c(-1.5,1.5)
  hallleft= log(left)
  hallold= log(tmat)

  beta_neww= rep(0,2)
  mysd= rep(0,2)

  final_x= xrc
  storage.mode(sigma_alpha)= "double"

  if(rfix==0)
  {
    out_rc= .Fortran("rc_r0", beta_neww= as.double(beta_neww), as.double(beta_old), as.integer(cumsumgrid),
                       as.integer(delta_temp), as.double(final_x), as.double(hallleft),
                       as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                       as.integer(lrcpos), as.integer(ltvec), as.integer(m),
                       as.integer(maxgrid), mysd= as.double(mysd),
                       as.integer(n), as.integer(notrcpos), as.integer(numgrid), 
                       as.integer(orderforhold), as.integer(rcpos), as.double(rfix),
                       sigma_alpha, as.integer(totngrid), as.double(ttilmat),
                       as.double(tvec), as.double(wbar), as.double(z2))
  }
  else
  {
  out_rc= .Fortran("rc_genr", beta_neww= as.double(beta_neww), as.double(beta_old), as.integer(cumsumgrid),
                       as.integer(delta_temp), as.double(final_x), as.double(hallleft),
                       as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                       as.integer(lrcpos), as.integer(ltvec),
                       as.integer(m), as.integer(maxgrid), mysd= as.double(mysd),
                       as.integer(n), as.integer(notrcpos), as.integer(numgrid), 
                       as.integer(orderforhold), as.integer(rcpos), as.double(rfix),
                       sigma_alpha, as.integer(totngrid), as.double(ttilmat),
                       as.double(tvec), as.double(wbar), as.double(z2))
  
  }

return_list= list(beta1.est= out_rc$beta_neww[1], beta2.est= out_rc$beta_neww[2], beta1.sd= out_rc$mysd[1], beta2.sd= out_rc$mysd[2])

}

##############################################################################
###############     IM method when wmat has 1 column     #####################


imw1= function(datamat,wmat,rfix,gridlen,ntimp,nximp,sigma2u)
{
  n= nrow(datamat)  
  m= ntimp             
  result= NULL

  delta_temp= datamat[,3]
  k= sum(delta_temp)
  
  left= datamat[,1]
  right1= datamat[,2]
  z2= datamat[,4]
  
  rcpos= which(delta_temp==0)
  notrcpos= (1:n)[-rcpos]
 
  liequalri= which(right1[notrcpos]==left[notrcpos])
  liminusrismall= which(right1[notrcpos]-left[notrcpos]<= gridlen)
  toosmall= notrcpos[liminusrismall]
  ltoo= length(toosmall)


  tmat= NULL
  numgrid= NULL
  csg= 0
  cumsumgrid= NULL
  for(i5 in 1:k)
  {
    i6= notrcpos[i5]

    newseq= seq(left[i6], right1[i6], gridlen)
      
    len_newseq= length(newseq)
    csg= csg+len_newseq
    tmat= c(tmat, newseq)
    numgrid= c(numgrid, len_newseq)
    cumsumgrid= c(cumsumgrid, csg)
  }

  maxgrid= max(numgrid)
  totngrid= cumsumgrid[k]

  lessthan1= liminusrismall[-match(liequalri,liminusrismall)]
  tmat[cumsumgrid[lessthan1]]= 0.5*(left[notrcpos[lessthan1]]+right1[notrcpos[lessthan1]])

  ttilmat= tmat
  wbar= wmat[,1]
    
  r= 10000
  isd= sample(1:90000,1)  


####################### Turnbull codes ###########################

  r6_left= round(left,3)
  r6_right= round(right1,3)
  turn_tau= sort(unique(c(r6_left,r6_right)))
  len_tt= length(turn_tau)
  survobj= survfit(Surv(turn_tau[1:(len_tt-1)],rep(1,len_tt-1))~1)
  init_s= c(1, survobj$surv)
  turn_p= -diff(init_s)

  turn_tau12= cbind(turn_tau[-len_tt],turn_tau[-1])
  alpha_fn= function(xx,infi,supr) ifelse(xx[1]>=infi & xx[2]<=supr, 1, 0)
  alpha_mat= apply(turn_tau12,1,alpha_fn,infi=left,supr=right1)
  id_zero= which(apply(alpha_mat==0,1,all))
  if(length(id_zero)>0) alpha_mat= alpha_mat[-id_zero,]
  
  new_s= init_s
  
  turn_err= 1
  while(turn_err>0.01)
  {
    old_s= new_s
    turn_p= -diff(old_s)
    dj_temp= t(t(alpha_mat)*turn_p)/as.vector(alpha_mat%*%turn_p)
    turn_dj= apply(dj_temp,2,sum)
    turn_yj= rev(cumsum(rev(turn_dj)))
    new_s= c(1,cumprod(1-turn_dj/turn_yj))
    new_s[which(is.na(new_s)==T)]= 0
    turn_err= sum(abs((new_s[-len_tt]-old_s[-len_tt])))      #/old_s[-len_tt]))
  }  
  
  surv_left= new_s[match(r6_left,turn_tau)]
  surv_right= new_s[match(r6_right,turn_tau)]


  surv_left= log(surv_left)
  surv_right= log(1+surv_right)

  lrcpos= length(rcpos)
  slsl= surv_left*surv_left
  srsr= surv_right*surv_right

  des_mat=cbind(1, z2, surv_left, surv_right, z2*surv_left, z2*surv_right, surv_left*surv_right)
  
  out100=lm(wbar~ des_mat-1)
  my_gamma_init= as.vector(out100$coef)
  my_gamma_init= c(my_gamma_init)


  n_el=length(my_gamma_init)
  nn= n_el*(n_el+1)/2
  nn_el= n*n_el

  
  ncolp=ncol(des_mat)   
  storage.mode(des_mat)<-"double"
  capr=20000
  store_myx=matrix(0, nrow=1000, ncol=n)
  storage.mode(store_myx)<-"double"


  a_sigmax=1
  b_sigmax=1

  theta_rep= nximp
  mtr= m*theta_rep

  aic_vec= rep(0,4)
  aic_x= NULL
  aic_part1= 0

  post_s2u= sigma2u


  for(j1 in 1:4)
  {

    ncomp=j1
    store_sigma2x=matrix(0, nrow=capr, ncol=ncomp)
    store_gamma=matrix(0, nrow=capr, ncol=(ncomp*ncolp))
    store_vec=matrix(0, nrow=capr, ncol=ncomp)
    store_pi=matrix(0, nrow=capr, ncol=ncomp)

    storage.mode(store_gamma)<-"double"
    storage.mode(store_sigma2x)<-"double"
    storage.mode(store_vec)<-"double"
    storage.mode(store_pi)<-"double"

    prmn=matrix(rnorm(ncomp*ncolp), nrow=ncomp, ncol=ncolp)
    prmn[, 1]=quantile(wbar,   prob=1/(ncomp+1)+(0:(ncomp-1))/(ncomp+1))
    prvar=matrix(5, nrow=ncomp, ncol=ncolp)
    storage.mode(prmn)<-"double"
    storage.mode(prvar)<-"double"

    outxk=.Fortran("mh_parkcmpw1", 
    des_mat, isd=as.integer(runif(1, 10, 100000)), mm=as.integer(ncol(wmat)), 
    n=as.integer(n), ncolp=as.integer(ncolp), ncomp=as.integer(ncomp), 
    ndim=as.integer(ncomp*ncolp),  capr=as.integer(capr), wbar=as.double(wbar),
    a_sigmax=as.double(a_sigmax), b_sigmax=as.double(b_sigmax), prmn,
    prvar, output1=store_gamma, output2=store_myx, output3=store_pi, 
    sigma2u= as.double(sigma2u), output5=store_sigma2x, output6=store_vec)


    if(ncomp==1)
    {
      post_gamma= apply(outxk$output1[(0.5*capr):capr,],2,mean)
      post_s2x= mean(outxk$output5[(0.5*capr):capr,])
      post_pi= mean(outxk$output3[(0.5*capr):capr,])
    }
    else
    {
      post_gamma= apply(outxk$output1[(0.5*capr):capr,],2,mean)
      post_s2x= apply(outxk$output5[(0.5*capr):capr,],2,mean)
      post_pi= apply(outxk$output3[(0.5*capr):capr,],2,mean)
    }
        
    aic_part2= 0
    for(j2 in 1:ncomp)
    {
      aic_part3= -0.5*(log(post_s2x[j2]+0.5*post_s2u)-log(0.5*post_s2x[j2]*post_s2u))+
      log(post_pi[j2])-0.5*log(post_s2x[j2])
      aic_part4= (wbar-des_mat%*%post_gamma[((j2-1)*ncolp+1):(j2*ncolp)])^2
      aic_part2= aic_part2 + exp(-0.5*aic_part4/(post_s2x[j2]+0.5*post_s2u)+ aic_part3)
    }
    lik_w= -0.5*aic_part1/post_s2u - n*log(2*pi*post_s2u) +sum(log(aic_part2))
    aic_vec[j1]= 2*(ncolp*ncomp+ 2*ncomp+ ncomp-1 + 1)- 2*lik_w
    aic_x=cbind(aic_x, t(outxk$output2[seq(1000,1000-(theta_rep-1)*50,-50),]))

  }

  aic_sel= which(aic_vec==min(aic_vec))
  final_x= aic_x[,(((aic_sel-1)*theta_rep+1):(aic_sel*theta_rep))]

######################   X imputation step complete   ###########################

  tvec= sort(unique(c(tmat,left)))
  ltvec= length(tvec)
  
  orderforhold= match(ttilmat,tvec)
  
  beta_old= c(-1.5,1.5)
  hallleft= log(left)
  hallold= log(tmat)

  beta_neww= rep(0,2)
  mysd= rep(0,2)

  
  storage.mode(final_x)= "double"
  
  allbeta_hope= matrix(0,2,mtr)
  storage.mode(allbeta_hope)= "double"

  if(rfix==0)
  {
    out_im= .Fortran("im_r0", beta_neww= as.double(beta_neww), as.double(beta_old), 
                       as.integer(cumsumgrid), as.integer(delta_temp), final_x, as.double(hallleft),
                       as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                       as.integer(lrcpos), as.integer(ltvec), as.integer(m), as.integer(maxgrid), 
                       as.integer(mtr), mysd= as.double(mysd), as.integer(n), as.integer(notrcpos), 
                       as.integer(numgrid), as.integer(orderforhold), as.integer(rcpos), 
                       as.double(rfix), as.integer(theta_rep), as.integer(totngrid), 
                       as.double(ttilmat), as.double(tvec), as.double(z2))
  }
  else
  {
    out_im= .Fortran("im_genr", beta_neww= as.double(beta_neww), as.double(beta_old),
                  as.integer(cumsumgrid), as.integer(delta_temp), final_x, as.double(hallleft),
                  as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                  as.integer(lrcpos), as.integer(ltvec), as.integer(m), as.integer(maxgrid),
                  as.integer(mtr), mysd= as.double(mysd), as.integer(n),
                  as.integer(notrcpos), as.integer(numgrid), as.integer(orderforhold),
                  as.integer(rcpos), as.double(rfix), as.integer(theta_rep),
                  as.integer(totngrid), as.double(ttilmat), as.double(tvec), as.double(z2))
  }

return_list= list(beta1.est= out_im$beta_neww[1], beta2.est= out_im$beta_neww[2], beta1.sd= out_im$mysd[1], beta2.sd= out_im$mysd[2])

}


#########################################################################################
###########   RC method when wmat has 1 column  ################


rcw1= function(datamat,wmat,rfix,gridlen,ntimp,sigma2u)
{
  n= nrow(datamat)
  m= ntimp
  result= NULL

  delta_temp= datamat[,3]
  k= sum(delta_temp)
  nrep= ncol(wmat)
  
  left= datamat[,1]
  right1= datamat[,2]
  z2= datamat[,4]
  
  rcpos= which(delta_temp==0)
  notrcpos= (1:n)[-rcpos]
  lrcpos= length(rcpos)
 
  liequalri= which(right1[notrcpos]==left[notrcpos])
  liminusrismall= which(right1[notrcpos]-left[notrcpos]<= gridlen)
  toosmall= notrcpos[liminusrismall]
  ltoo= length(toosmall)


  tmat= NULL
  numgrid= NULL
  csg= 0
  cumsumgrid= NULL
  for(i5 in 1:k)
  {
    i6= notrcpos[i5]

    newseq= seq(left[i6], right1[i6], gridlen)
      
    len_newseq= length(newseq)
    csg= csg+len_newseq
    tmat= c(tmat, newseq)
    numgrid= c(numgrid, len_newseq)
    cumsumgrid= c(cumsumgrid, csg)
  }

  maxgrid= max(numgrid)
  totngrid= cumsumgrid[k]

  lessthan1= liminusrismall[-match(liequalri,liminusrismall)]
  tmat[cumsumgrid[lessthan1]]= 0.5*(left[notrcpos[lessthan1]]+right1[notrcpos[lessthan1]])

  ttilmat= tmat
  wbar= wmat[,1]

  r= 10000
  isd= sample(1:90000,1)


##############   regression calibration estimate of X    #################

  estsigmau2= sigma2u
  wibarminuswbar= wbar-mean(wbar)

  nu= (n-1)*nrep
  estsigmax2= (nrep*sum(wibarminuswbar^2) - (n-1)*estsigmau2)/nu
  estsigmaxz= nrep*sum(wibarminuswbar*(z2-mean(z2)))/nu
  estsigmaz2= var(z2)
  
  rcmat= matrix(c(estsigmax2+estsigmau2/nrep, estsigmaxz, estsigmaxz, estsigmaz2),2,2)

  xrc= as.vector(mean(wbar)+ c(estsigmax2, estsigmaxz)%*% solve(rcmat)%*% rbind(wibarminuswbar, z2-mean(z2)))

# Bootstrap to calculate sigma_alpha
  maxbt= 200   # number of bootstrap samples
  bootout= matrix(0,maxbt,3)
  for(bt in 1:maxbt)
  {
    boot= sample(1:n,n,replace=T)
    z2boot= z2[boot]
    wmatboot= wmat[boot,]
    wbarboot= wmatboot
    wibarminuswbar= wbarboot-mean(wbarboot)

    nu= (n-1)*nrep
    estsigmax2= (nrep*sum(wibarminuswbar^2) - (n-1)*estsigmau2)/nu
    estsigmaxz= nrep*sum(wibarminuswbar*(z2boot-mean(z2boot)))/nu
    estsigmaz2= var(z2boot)

    mydet= estsigmaz2*(estsigmax2 + estsigmau2/nrep) - estsigmaxz*estsigmaxz
    alpha1= (estsigmax2*estsigmaz2 - estsigmaxz*estsigmaxz)/mydet
    alpha2= estsigmaxz*estsigmau2/(nrep*mydet)
    
    alpha0= mean(wbarboot)*(1-alpha1)-alpha2*mean(z2boot)
    
    bootout[bt,]= c(alpha0,alpha1,alpha2)
  }

  sigma_alpha= cov(bootout)
  tvec= sort(unique(c(tmat,left)))
  ltvec= length(tvec)
  
  orderforhold= match(ttilmat,tvec)
  
  beta_old= c(-1.5,1.5)
  hallleft= log(left)
  hallold= log(tmat)

  beta_neww= rep(0,2)
  mysd= rep(0,2)

  final_x= xrc
  storage.mode(sigma_alpha)= "double"
  allbeta_hope= matrix(0,2,m)
  storage.mode(allbeta_hope)= "double"


  if(rfix==0)
  {
    out_rc= .Fortran("rc_r0", beta_neww= as.double(beta_neww), as.double(beta_old), as.integer(cumsumgrid),
                       as.integer(delta_temp), as.double(final_x), as.double(hallleft),
                       as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                       as.integer(lrcpos), as.integer(ltvec), as.integer(m),
                       as.integer(maxgrid), mysd= as.double(mysd),
                       as.integer(n), as.integer(notrcpos), as.integer(numgrid), 
                       as.integer(orderforhold), as.integer(rcpos), as.double(rfix),
                       sigma_alpha, as.integer(totngrid), as.double(ttilmat),
                       as.double(tvec), as.double(wbar), as.double(z2))
  }
  else
  {
  out_rc= .Fortran("rc_genr", beta_neww= as.double(beta_neww), as.double(beta_old), as.integer(cumsumgrid),
                       as.integer(delta_temp), as.double(final_x), as.double(hallleft),
                       as.double(hallold), as.integer(isd), as.integer(k), as.double(left),
                       as.integer(lrcpos), as.integer(ltvec),
                       as.integer(m), as.integer(maxgrid), mysd= as.double(mysd),
                       as.integer(n), as.integer(notrcpos), as.integer(numgrid), 
                       as.integer(orderforhold), as.integer(rcpos), as.double(rfix),
                       sigma_alpha, as.integer(totngrid), as.double(ttilmat),
                       as.double(tvec), as.double(wbar), as.double(z2))
  
  }

return_list= list(beta1.est= out_rc$beta_neww[1], beta2.est= out_rc$beta_neww[2], beta1.sd= out_rc$mysd[1], beta2.sd= out_rc$mysd[2])

}


