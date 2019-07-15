!==========================================================================
!          PART 1: SUBROUTINES USED IN WHILE LOOP
!==========================================================================


! IMPUTATION OF t
! output: x_ik

subroutine imp_tr0(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                 lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                 x_ik, rcpos, rfix, totngrid, ttilmat, z2)
implicit none

integer(4):: i5, i6, i8, i10, i11, i12, k, m, n, lrcpos, sample_wr

integer(4):: cumsumgrid(k), maxgrid, myindx, notrcpos(k),&
            numgrid(k), outval, outx, rcpos(lrcpos), totngrid

real(8):: beta_old(2), bz(n), bzplusleft(k), bzplusleft_temp(n),&
         den, hallold(totngrid), hallleft(n), left(n), newz1(n),&
         num, outch1, outch2, outch3, pvec(maxgrid), ratio(maxgrid),&
         rfix, ttilmat(totngrid), x_ik(n,m), z2(n)

integer(4) :: is1, is2    

  common /unif_seeds/ is1, is2
  save /unif_seeds/


    x_ik(:,:)= 0.d0

    bz= beta_old(1)*newz1+beta_old(2)*z2

    bzplusleft(1)= bz(notrcpos(1))+ hallold(1)

    do i10= 2, k
      myindx= notrcpos(i10)
      bzplusleft(i10)= bz(myindx)+ hallold(cumsumgrid(i10-1)+1)
    end do

    myindx= notrcpos(1)

    outch1= exp(bzplusleft(1))
    outch2= exp(bz(myindx)+ hallold(cumsumgrid(1)))
    den= exp(-outch1)- exp(-outch2)

    ratio= 0.d0
    pvec= 0.d0
    
    outch3= exp(bz(myindx)+ hallold(1))
    num= exp(-outch1)- exp(-outch3)
    ratio(1)= num/den
    pvec(1)= ratio(1)

    do i6= 2, numgrid(1)
      outch3= exp(bz(myindx)+ hallold(i6))
      num= exp(-outch1)- exp(-outch3)
      ratio(i6)= num/den
      pvec(i6)= ratio(i6)-ratio(i6-1)
    end do

    if(den<0.000001d0) then
      pvec(1:numgrid(1))= 0.5d0
    end if


    do i8= 1, m
      outx= sample_wr(maxgrid, pvec)
      x_ik(myindx, i8)= ttilmat(outx)
!     print*, "outx= ", outx
    end do

  do i5= 2, k
      myindx= notrcpos(i5)

    outch1= exp(bzplusleft(i5))
    outch2= exp(bz(myindx)+ hallold(cumsumgrid(i5)))
    den= exp(-outch1)- exp(-outch2)

    ratio= 0.d0
    pvec= 0.d0
    
    outch3= exp(bz(myindx)+ hallold(cumsumgrid(i5-1)+1))
    num= exp(-outch1)- exp(-outch3)
    ratio(1)= num/den
    pvec(1)= ratio(1)

    do i6= 2, numgrid(i5)
      outch3= exp(bz(myindx)+ hallold(cumsumgrid(i5-1)+i6))
      num= exp(-outch1)- exp(-outch3)
      ratio(i6)= num/den
      pvec(i6)= ratio(i6)-ratio(i6-1)
    end do

    
    if(den<0.000001d0) then
      pvec(1:numgrid(i5))= 0.5d0
    end if

    do i8= 1, m
      outx= sample_wr(maxgrid, pvec)
      x_ik(myindx, i8)= ttilmat(cumsumgrid(i5-1)+outx)
!     print*, "outx= ", outx
    end do

  end do

  do i11= 1, lrcpos
    x_ik(rcpos(i11),:)= left(rcpos(i11))
  end do

! print*, "finished"

return
end subroutine


!============================================================================
!                PART 2: THE BETA AND H ESTIMATION
!                        AT EACH STEP OF ITERATION
!============================================================================

! BETA ESTIMATION
! output: beta_old (takes the same value as beta_new at the end)

subroutine betaestr0(beta_old, delta, hall, k, n, r, z)
implicit none

integer(4) :: k, n
integer(4) :: delta(n), i6

real(8) :: adj_df_beta(2,2), beta_old(2), beta_new(2), bz(n),&
          bzplush(n), det_df_beta, df_beta(2,2),&
          err, f_beta(2), hall(n), inv_df_beta(2,2),&
          outh(n), r, restofbetafn(n), z(n,2), zzprime(2,2)



    err= 1.d0
    do while(err> 0.01d0)
  
    bz= z(:, 1)*beta_old(1)+z(:, 2)*beta_old(2)
    bzplush= bz+hall
    
    outh= exp(bzplush)
    restofbetafn= delta-outh    ! for r=0, haz_fn = chaz_fn
    f_beta=(/sum(restofbetafn*z(:, 1)), sum(restofbetafn*z(:, 2))/)
!!  print*, f_beta

    df_beta= 0.d0
    do i6= 1, n
      zzprime(1,1)= z(i6, 1)*z(i6, 1)
      zzprime(1,2)= z(i6, 1)*z(i6, 2) 
      zzprime(2,1)= zzprime(1,2)
      zzprime(2,2)= z(i6, 2)*z(i6, 2)
      df_beta= df_beta - zzprime*outh(i6)
    end do
!!  print*, dbetafn
  
    det_df_beta= 1/(df_beta(1,1)*df_beta(2,2)-df_beta(1,2)*df_beta(2,1))
    adj_df_beta(1,:)= (/df_beta(2,2), -df_beta(1,2)/)
    adj_df_beta(2,:)= (/-df_beta(2,1), df_beta(1,1)/)
    inv_df_beta= adj_df_beta*det_df_beta
!!  print*, inv_df_beta

    beta_new(1)= beta_old(1)-inv_df_beta(1,1)*f_beta(1)-inv_df_beta(1,2)*f_beta(2)
    beta_new(2)= beta_old(2)-inv_df_beta(2,1)*f_beta(1)-inv_df_beta(2,2)*f_beta(2) 
    err= sum(abs((beta_new-beta_old)/beta_old))
    beta_old= beta_new
!   print*, err
!   print*, beta_old
  end do


return
end subroutine




!=============================================================
!          IMPUTATION METHOD for r=0
!============================================================

subroutine im_r0(beta_neww, beta_old, cumsumgrid,&
                         delta_temp, final_x, hallleft,&
                         hallold, isd, k, left,&
                         lrcpos, ltvec, m, maxgrid, mtr,&
                         mysd, n, notrcpos, numgrid,&
                         orderforhold, rcpos, rfix,&
                         theta_rep, totngrid, ttilmat,&
                         tvec, z2)
implicit none


integer(4):: avhid, ii, j, j2, k, k1, kk, kminus1, kminus2,&
            lrcpos, ltvec, m, maxgrid, mtr,&
            n, t1, t2, theta_rep, totngrid, tt,&
            sm, dellater, dellater1, dellater2
            
integer(4):: bfnindx1, bfnindx2, cumsumgrid(k), delta_temp(n),&
            left_tvec_match(n), notrcpos(k), numgrid(k),&
            orderforhold(totngrid), rcpos(lrcpos)

real(8):: beta_difff, beta_neww(2), beta_old(2), final_x(n, theta_rep),&
         left(n), mysd(2), rfix, ttilmat(totngrid), tvec(ltvec), z2(n)

real(8):: beta_older(2), beta_old_mean(2), hallleft(n),&
         hallold(totngrid), init_beta(2), newhavg(ltvec),&
         newz1(n), x_ik(n,m), z(n,2)


integer(4):: cc, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10,&
            allind(k,m), delone(k), delta(n), inddel(k), ind_del(k),&
            this, x_ik_sort_indx(n)

real(8):: allbeta(2,m), alldel(n,m), allh(k,m), allhall(n,m), allx_ymat(n,k,m),&
         allz(n,2,m), beta_old_sum(2), beta_old2(2), bz(n), hall(n), hall_old(n), hnew(k),&
         init_beta2(2), newhsum(ltvec), outvec(ltvec), theta_difff, theta_old(2),&
         x_ik_sort(n), x_t1t2(k), x_ymat(n,k), z_temp2(n,2)

real(8):: adj_siglstar(4), ahvar_mean(2), allbeta_hope(2,mtr),&
         allvar_hope(2,mtr), beta_diff(2,m), beta_hope(2,mtr),&
         beta_hope_beta_hope_t(2), bfnvec(n), bfnvec1(n), det_siglstar,&
         dhazy(n), dhazyhdiff(n), firstel1(n), firstel2(n), fourthel1(n),&
         fourthel2(n), hazy(n), hazyhdiff(n), hdiff(k-1), integrand

real(8):: mn, mybzplush(n), mybzvec(n), myfinalbeta(2), myfinalh(k),&
         myfinalhall(n), outdh(n), outh(n), outh1(n), outh2(n), outh3(n),&
         oyb(n), ratio, secondel1(n), secondel2(n), sig1(4),&
         sig2(4), sig_hope(2), sum_int, theta_diff, thirdel2(n),&
         varbeta(2), ymat(n,k), zbarvec(2), zminuszbar(n,2),&
         zzprime1(4), zzprime2(4)

real(8):: siglstar(4), sigustar(4), inv_siglstar(4)

real(8):: err1, err2, fxl1, fxl, fxmid1, fxmid, fxu1,&
         fxu, term1, xl, xmid, xu 

integer(4) :: is1, is2, isd, seed1, seed2        

    common /unif_seeds/ is1, is2
    save /unif_seeds/
    seed1=10005
    seed2=isd      
    CALL set_uniform(seed1, seed2)
  
    
    init_beta= beta_old
    beta_neww= 0.d0
    allbeta_hope= 0.d0
    allvar_hope= 0.d0
    hall_old= 0.d0


    do tt=1, theta_rep

      beta_old= init_beta
      newz1= final_x(:,tt)
      beta_difff= 1.0d0

      do while(beta_difff > 0.01d0)
      
        beta_older= beta_old

        CALL imp_tr0(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                   lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                   x_ik, rcpos, rfix, totngrid, ttilmat, z2)

!!      output: x_ik


!=============  BETA AND H ESTIMATION STEP ======================


        z_temp2(:,1)= newz1       
        z_temp2(:,2)= z2
  
        beta_old_sum= 0.d0
        newhsum= 0.d0
        init_beta2= beta_old


        do i6= 1, m
 
          beta_old2= init_beta2
          hall= 0.d0
  
          x_ik_sort= x_ik(:,i6)
          CALL bubble(x_ik_sort, n)           ! equivalent of sort()$x in R
          CALL my_match(x_ik_sort_indx, n, n, x_ik_sort, x_ik(:,i6))
                                              ! equivalent of sort()$ix in R
    
          z= z_temp2(x_ik_sort_indx,:)
          delta= delta_temp(x_ik_sort_indx)

          allz(:,:,i6)= z
          alldel(:,i6)= delta


          cc=1
          do i8= 1, n
            if(delta(i8)==1.d0) then
              x_t1t2(cc)= x_ik_sort(i8)
              inddel(cc)= i8
              cc= cc+1
            end if
          end do
    
          x_ymat=0.d0    
          do i7= 1, k
            where(x_ik_sort>= x_t1t2(i7)) x_ymat(:,i7)= 1.d0
          end do
    
          allind(:, i6)= inddel

      allx_ymat(:,:,i6)= x_ymat

      theta_difff= 1.d0
      do while(theta_difff> 0.01d0)
        theta_old= beta_old2

!       CALL hest(beta_old2, delta, hall, hnew, k, n, rfix, x_ymat, z)

!===============  H ESTIMATION STEP  ==========================


      xl= -20.d0
      xu= 20.d0
      xmid= 0.d0
      bz= z(:, 1)*beta_old2(1)+z(:, 2)*beta_old2(2)
  
  
      err1= 1
      do while(err1> 0.01d0)
 
        fxl1= sum(x_ymat(:,1)*exp(bz+xl))-1
        fxu1= sum(x_ymat(:,1)*exp(bz+xu))-1

        if (fxl1*fxu1.lt.0.d0) then
          xmid= (xl+xu)/2.d0
          fxmid1= sum(x_ymat(:,1)*exp(bz+xmid))-1
          if(fxmid1.lt.0.d0) then
            xl= xmid
          else
            xu= xmid
          end if

        else
!          print*, 'try different extremes 1'
          exit

        end if
        err1= abs(fxmid1)

      end do
    
      hnew(1)= xmid
  
      do i1= 2, k
        xl= -15.d0
        xu= 15.d0
        xmid= 0.d0

        err2= 1.d0
        term1=1+sum(x_ymat(:,i1)*exp(bz+hnew(i1-1)))

        do while(err2> 0.01d0)

          fxl= sum(x_ymat(:,i1)*exp(bz+xl))-term1
          fxu= sum(x_ymat(:,i1)*exp(bz+xu))-term1

          if (fxl*fxu.lt.0d0) then
            xmid= (xl+xu)/2.d0
            fxmid= sum(x_ymat(:,i1)*exp(bz+xmid))-term1
            if(fxmid.lt.0d0) then
              xl= xmid
            else
              xu= xmid
            end if

          else
  !          print*, 'try different extremes 2'
            exit
          end if
          err2= abs(fxmid)

        end do
        hnew(i1)= xmid
      end do
  
      hall= -50.d0

      i5= 0
      do i4= 1, n
        if(delta(i4).eq.1) then
          i5= i5+1
          delone(i5)= i4
        end if
      end do


      i3= 0
      do i2= minval(delone), n
        if(delta(i2).eq.1) then
          i3= i3+1
          hall(i2)= hnew(i3)
      
        else
          hall(i2)= hnew(i3)
        end if
      end do

!===============   H ESTIMATION ENDS HERE  =====================


      CALL betaestr0(beta_old2, delta, hall, k, n, rfix, z)

      theta_difff= sum(abs((beta_old2-theta_old)/theta_old))

    end do   !  end of while loop for theta_difff
    
    allbeta(:,i6)= beta_old2
    allhall(:,i6)= hall
    allh(:,i6)= hnew

    CALL newhallr0(hall, ltvec, n, outvec, tvec, x_ik_sort)

    newhsum= newhsum + outvec
    beta_old_sum= beta_old_sum + beta_old2

  end do    ! end of for loop i6= 1, m

 
  beta_old_mean= beta_old_sum/m
  newhavg= newhsum/m

!============  BETA AND H ESTIMATION ENDS HERE  =====================


        beta_old= beta_old_mean

        CALL my_match(left_tvec_match, n, ltvec, left, tvec)
!!      output: left_tvec_match

        hallleft= newhavg(left_tvec_match)
        hallold= newhavg(orderforhold)
 
        beta_difff= sum(abs((beta_old-beta_older)/beta_older))
        
      end do   ! end of while loop for beta convergence

      
!     print*, 'final_beta_theta= ', beta_old
!     beta_neww= beta_neww+beta_old
      
      
!!========================================
!!            Sigma Computation
!!========================================

     do kk= 1, m
        ymat= allx_ymat(:,:,kk)
        myfinalbeta= allbeta(:,kk)
        myfinalh= allh(:,kk)
        myfinalhall= allhall(:,kk)
        z= allz(:,:,kk)
            
        hdiff= myfinalh(2:k)-myfinalh(1:(k-1))
        mybzvec= z(:,1)*myfinalbeta(1)+z(:,2)*myfinalbeta(2)

        outh1= exp(mybzvec+myfinalhall)

        sig1=0.d0
        sig2=0.d0

        firstel1= 0.d0
        secondel1= 0.d0
        fourthel1= 0.d0
        
        firstel2= 0.d0
        secondel2= 0.d0
        thirdel2= 0.d0
        fourthel2= 0.d0

        do t1= 2, k
          bfnvec= exp(myfinalh(t1)-myfinalhall)
          mybzplush= mybzvec+myfinalh(t1)
         
          outh2= exp(mybzplush)
          hazy= outh2*ymat(:,t1)

          dhazy= hazy
     
          oyb= outh1*ymat(:,t1)*bfnvec
          zbarvec(1)= sum(z(:,1)*oyb)/sum(hazy)
          zbarvec(2)= sum(z(:,2)*oyb)/sum(hazy)

          zminuszbar(:,1)= z(:,1)-zbarvec(1)
          zminuszbar(:,2)= z(:,2)-zbarvec(2)

          hazyhdiff= hazy*hdiff(t1-1)
          dhazyhdiff= dhazy*hdiff(t1-1)

          firstel1= firstel1+zminuszbar(:,1)*zminuszbar(:,1)*hazyhdiff
          secondel1= secondel1+zminuszbar(:,1)*zminuszbar(:,2)*hazyhdiff
          fourthel1= fourthel1+zminuszbar(:,2)*zminuszbar(:,2)*hazyhdiff

          firstel2= firstel2+zminuszbar(:,1)*z(:,1)*dhazyhdiff
          secondel2= secondel2+zminuszbar(:,2)*z(:,1)*dhazyhdiff
          thirdel2= thirdel2+zminuszbar(:,1)*z(:,2)*dhazyhdiff
          fourthel2= fourthel2+zminuszbar(:,2)*z(:,2)*dhazyhdiff

        end do  ! end of for loop t1= 2, k
   
        sig1= (/sum(firstel1), sum(secondel1), sum(secondel1), sum(fourthel1)/)
        sig2= (/sum(firstel2), sum(secondel2), sum(thirdel2), sum(fourthel2)/)

        sigustar= sig1/n
        siglstar= sig2/n


        det_siglstar= 1.d0/(siglstar(1)*siglstar(4)-siglstar(2)*siglstar(3))
        adj_siglstar= (/siglstar(4), -siglstar(2), -siglstar(3), siglstar(1)/)
        inv_siglstar= adj_siglstar*det_siglstar
        
        varbeta(1)= inv_siglstar(1)*inv_siglstar(1)*sigustar(1) +&
                      inv_siglstar(1)*inv_siglstar(3)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(3)*inv_siglstar(3)*sigustar(4)
        varbeta(2)= inv_siglstar(2)*inv_siglstar(2)*sigustar(1) +&
                      inv_siglstar(2)*inv_siglstar(4)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(4)*inv_siglstar(4)*sigustar(4)


        avhid= m*(tt-1)+kk
        allvar_hope(1,avhid)= varbeta(1)/n
        allvar_hope(2,avhid)= varbeta(2)/n

        allbeta_hope(:,avhid)= allbeta(:,kk)

      end do      ! end of do loop for kk= 1, m

    end do        ! end of do loop for theta_rep

    beta_neww(1)= sum(allbeta_hope(1,:))/mtr
    beta_neww(2)= sum(allbeta_hope(2,:))/mtr


    ahvar_mean(1)= sum(allvar_hope(1,:))/mtr
    ahvar_mean(2)= sum(allvar_hope(2,:))/mtr

    beta_hope= allbeta_hope - spread(beta_neww,2,mtr)
    beta_hope_beta_hope_t(1)= sum(beta_hope(1,:)*beta_hope(1,:))
    beta_hope_beta_hope_t(2)= sum(beta_hope(2,:)*beta_hope(2,:))

    sig_hope= ahvar_mean + (mtr+1)*beta_hope_beta_hope_t/(mtr*(mtr-1))
    mysd= SQRT(sig_hope)


!print*, 'final ans= ', (/beta_neww,mysd/)
!print*, (/beta_neww,mysd/)

return
end subroutine






!========================================================================
! THIS FUNCTION IS USED IN THE BETA AND H ESTIMATION PART
! this function extends H estimates over points where H was not estimated
! input: a vector of SORTED t values where H has been estimated (xvec)
!        a vector of the corresponding H estimates (Hvec)
!        a vector over which the H estimates are to be expanded using the continuity property (tvec)
!        ltvec= length(tvec)
!        n= length(xvec)
! output: a vector of H estimates at all t-values included in tvec (outvec)

subroutine newhallr0(hvec, ltvec, n, outvec, tvec, xvec)
implicit none

integer(4):: i, ltvec, n
integer(4):: myord(n)
real(8):: hvec(n), minout, outvec(ltvec), tvec(ltvec), xvec(n)

  outvec(:)= 0.d0

  CALL my_match(myord, n, ltvec, xvec, tvec)
  do i= 1, n-1
    outvec(myord(i):(myord(i+1)-1))= hvec(i)
  end do

  outvec(myord(n):ltvec)= maxval(hvec)

  minout= minval(pack(outvec, outvec/=0.d0))
  where(outvec==0.d0)  outvec= minout-abs(minout)
return
end subroutine







!==========================================================================
!          PART 3: SUBROUTINES USED IN WHILE LOOP general r case
!==========================================================================

! IMPUTATION OF t
! output: x_ik

subroutine imp_t(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                 lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                 x_ik, rcpos, rfix, totngrid, ttilmat, z2)
implicit none

integer(4):: i5, i6, i8, i10, i11, i12, k, m, n, lrcpos, sample_wr

integer(4):: cumsumgrid(k), maxgrid, myindx, notrcpos(k),&
            numgrid(k), outval, outx, rcpos(lrcpos), totngrid

real(8):: beta_old(2), bz(n), bzplusleft(k), bzplusleft_temp(n),&
         den, hallold(totngrid), hallleft(n), left(n), newz1(n),&
         num, outch1, outch2, outch3, pvec(maxgrid), ratio(maxgrid),&
         rfix, ttilmat(totngrid), x_ik(n,m), z2(n)

integer(4) :: is1, is2    

  common /unif_seeds/ is1, is2
  save /unif_seeds/

!   print*, "loc 1"

    x_ik(:,:)= 0.d0

    bz= beta_old(1)*newz1+beta_old(2)*z2

    bzplusleft(1)= bz(notrcpos(1))+ hallold(1)

    do i10= 2, k
      myindx= notrcpos(i10)
      bzplusleft(i10)= bz(myindx)+ hallold(cumsumgrid(i10-1)+1)
    end do

!   print*, "loc 2"

    myindx= notrcpos(1)

    outch1= log(1+rfix*exp(bzplusleft(1)))/rfix
    outch2= log(1+rfix*exp(bz(myindx)+ hallold(cumsumgrid(1))))/rfix
    den= exp(-outch1)- exp(-outch2)

    ratio= 0.d0
    pvec= 0.d0
    
    outch3= log(1+rfix*exp(bz(myindx)+ hallold(1)))/rfix
    num= exp(-outch1)- exp(-outch3)
    ratio(1)= num/den
    pvec(1)= ratio(1)

!   print*, "loc 3"

    do i6= 2, numgrid(1)
      outch3= log(1+rfix*exp(bz(myindx)+ hallold(i6)))/rfix
      num= exp(-outch1)- exp(-outch3)
      ratio(i6)= num/den
      pvec(i6)= ratio(i6)-ratio(i6-1)
    end do

!   print*, "loc 4"

    if(den<0.000001d0) then
      pvec(1:numgrid(1))= 0.5d0
    end if


    do i8= 1, m
      outx= sample_wr(maxgrid, pvec)
      x_ik(myindx, i8)= ttilmat(outx)
!     print*, "outx= ", outx
    end do


  do i5= 2, k
      myindx= notrcpos(i5)

    outch1= log(1+rfix*exp(bzplusleft(i5)))/rfix
    outch2= log(1+rfix*exp(bz(myindx)+ hallold(cumsumgrid(i5))))/rfix
    den= exp(-outch1)- exp(-outch2)

    ratio= 0.d0
    pvec= 0.d0
    
    outch3= log(1+rfix*exp(bz(myindx)+ hallold(cumsumgrid(i5-1)+1)))/rfix
    num= exp(-outch1)- exp(-outch3)
    ratio(1)= num/den
    pvec(1)= ratio(1)

    do i6= 2, numgrid(i5)
      outch3= log(1+rfix*exp(bz(myindx)+ hallold(cumsumgrid(i5-1)+i6)))/rfix
      num= exp(-outch1)- exp(-outch3)
      ratio(i6)= num/den
      pvec(i6)= ratio(i6)-ratio(i6-1)
    end do

    
    if(den<0.000001d0) then
      pvec(1:numgrid(i5))= 0.5d0
    end if

    do i8= 1, m
      outx= sample_wr(maxgrid, pvec)
      x_ik(myindx, i8)= ttilmat(cumsumgrid(i5-1)+outx)
!     print*, "outx= ", outx
    end do

  end do

  do i11= 1, lrcpos
    x_ik(rcpos(i11),:)= left(rcpos(i11))
  end do

! print*, "finished"

return
end subroutine


!============================================================================
!                PART 4: THE BETA AND H ESTIMATION for general r case
!============================================================================

! BETA ESTIMATION
! output: beta_old (takes the same value as beta_new at the end)

subroutine betaest(beta_old, delta, hall, k, n, r, z)
implicit none

integer(4) :: k, n
integer(4) :: delta(n), i6

real(8) :: adj_df_beta(2,2), beta_old(2), beta_new(2), bz(n),&
          bzplush(n), det_df_beta, df_beta(2,2),&
          err, f_beta(2), hall(n), inv_df_beta(2,2),&
          outch(n), outh(n), r, restofbetafn(n), z(n,2), zzprime(2,2)



    err= 1.d0
    do while(err> 0.01d0)
  
    bz= z(:, 1)*beta_old(1)+z(:, 2)*beta_old(2)
    bzplush= bz+hall
    
    outh= exp(bzplush)/(1+r*exp(bzplush))
    outch= log(1+r*exp(bzplush))/r
    
    restofbetafn= delta-outch    ! for r=0, haz_fn = chaz_fn
    f_beta=(/sum(restofbetafn*z(:, 1)), sum(restofbetafn*z(:, 2))/)
!!  print*, f_beta

    df_beta= 0.d0
    do i6= 1, n
      zzprime(1,1)= z(i6, 1)*z(i6, 1)
      zzprime(1,2)= z(i6, 1)*z(i6, 2) 
      zzprime(2,1)= zzprime(1,2)
      zzprime(2,2)= z(i6, 2)*z(i6, 2)
      df_beta= df_beta - zzprime*outh(i6)
    end do
!!  print*, dbetafn
  
    det_df_beta= 1/(df_beta(1,1)*df_beta(2,2)-df_beta(1,2)*df_beta(2,1))
    adj_df_beta(1,:)= (/df_beta(2,2), -df_beta(1,2)/)
    adj_df_beta(2,:)= (/-df_beta(2,1), df_beta(1,1)/)
    inv_df_beta= adj_df_beta*det_df_beta
!!  print*, inv_df_beta

    beta_new(1)= beta_old(1)-inv_df_beta(1,1)*f_beta(1)-inv_df_beta(1,2)*f_beta(2)
    beta_new(2)= beta_old(2)-inv_df_beta(2,1)*f_beta(1)-inv_df_beta(2,2)*f_beta(2) 
    err= sum(abs((beta_new-beta_old)/beta_old))
    beta_old= beta_new
!   print*, err
!   print*, beta_old
  end do


return
end subroutine




!=============================================================
!            IMPUTATION METHOD for general r case
!=============================================================

subroutine im_genr(beta_neww, beta_old, cumsumgrid,&
                       delta_temp, final_x, hallleft,&
                       hallold, isd, k, left,&
                       lrcpos, ltvec, m, maxgrid, mtr,&
                       mysd, n, notrcpos, numgrid,&
                       orderforhold, rcpos, rfix,&
                       theta_rep, totngrid, ttilmat,&
                       tvec, z2)
implicit none


integer(4):: avhid, ii, j, j2, k, k1, kk, kminus1, kminus2,&
            lrcpos, ltvec, m, maxgrid, mtr,&
            n, t1, t2, theta_rep, totngrid, tt,&
            sm, dellater, dellater1, dellater2
            
integer(4):: bfnindx1, bfnindx2, cumsumgrid(k), delta_temp(n),&
            left_tvec_match(n), notrcpos(k), numgrid(k),&
            orderforhold(totngrid), rcpos(lrcpos)

real(8):: beta_difff, beta_neww(2), beta_old(2), final_x(n, theta_rep),&
         left(n), mysd(2), rfix, ttilmat(totngrid), tvec(ltvec), z2(n)

real(8):: beta_older(2), beta_old_mean(2), hallleft(n),&
         hallold(totngrid), init_beta(2), newhavg(ltvec),&
         newz1(n), x_ik(n,m), z(n,2)


integer(4):: cc, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10,&
            allind(k,m), delone(k), delta(n), inddel(k), ind_del(k),&
            this, x_ik_sort_indx(n)

real(8):: allbeta(2,m), alldel(n,m), allh(k,m), allhall(n,m), allx_ymat(n,k,m),&
         allz(n,2,m), beta_old_sum(2), beta_old2(2), bz(n), hall(n), hnew(k),&
         init_beta2(2), newhsum(ltvec), outvec(ltvec), theta_difff, theta_old(2),&
         x_ik_sort(n), x_t1t2(k), x_ymat(n,k), z_temp2(n,2)

real(8):: adj_siglstar(4), ahvar_mean(2), allbeta_hope(2,mtr),&
         allvar_hope(2,mtr), beta_diff(2,m), beta_hope(2,mtr),&
         beta_hope_beta_hope_t(2), bfnvec(n), bfnvec1(n), det_siglstar,&
         dhazy(n), dhazyhdiff(n), firstel1(n), firstel2(n), fourthel1(n),&
         fourthel2(n), hazy(n), hazyhdiff(n), hdiff(k-1), integrand

real(8):: mn, mybzplushmat(n,k), mybzvec(n), myfinalbeta(2), myfinalh(k),&
         myfinalhall(n), outdh(n), outh(n), outh1(n), outh2(n), outh3(n),&
         oyb(n), ratio, secondel1(n), secondel2(n), sig1(4),&
         sig2(4), sig_hope(2), sum_int, theta_diff, thirdel2(n),&
         varbeta(2), ymat(n,k), zbarvec(2), zminuszbar(n,2),&
         zzprime1(4), zzprime2(4)

real(8):: siglstar(4), sigustar(4), inv_siglstar(4)

real(8):: err1, err2, fxl1, fxl, fxmid1, fxmid, fxu1,&
         fxu, term1, xl, xmid, xu 

integer(4) :: is1, is2, isd, seed1, seed2        

    common /unif_seeds/ is1, is2
    save /unif_seeds/
    seed1=10005
    seed2=isd      
    CALL set_uniform(seed1, seed2)
  
  
    init_beta= beta_old
    beta_neww= 0.d0
    allbeta_hope= 0.d0
    allvar_hope= 0.d0


    do tt=1, theta_rep

      beta_old= init_beta
      newz1= final_x(:,tt)
      beta_difff= 1.0d0

      do while(beta_difff > 0.01d0)
      
        beta_older= beta_old

        CALL imp_t(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                   lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                   x_ik, rcpos, rfix, totngrid, ttilmat, z2)

!!      output: x_ik


!=============  BETA AND H ESTIMATION STEP ======================


        z_temp2(:,1)= newz1       
        z_temp2(:,2)= z2
  
        beta_old_sum= 0.d0
        newhsum= 0.d0
        init_beta2= beta_old


        do i6= 1, m
 
          beta_old2= init_beta2
          hall= 0.d0
  
          x_ik_sort= x_ik(:,i6)
          CALL bubble(x_ik_sort, n)           ! equivalent of sort()$x in R
          CALL my_match(x_ik_sort_indx, n, n, x_ik_sort, x_ik(:,i6))
                                              ! equivalent of sort()$ix in R
    
          z= z_temp2(x_ik_sort_indx,:)
          delta= delta_temp(x_ik_sort_indx)

          allz(:,:,i6)= z
          alldel(:,i6)= delta

!          print*, "step 1 done"

          cc=1
          do i8= 1, n
            if(delta(i8)==1.d0) then
              x_t1t2(cc)= x_ik_sort(i8)
              inddel(cc)= i8
              cc= cc+1
            end if
          end do
    
          x_ymat=0.d0    
          do i7= 1, k
            where(x_ik_sort>= x_t1t2(i7)) x_ymat(:,i7)= 1.d0
          end do
    
          allind(:, i6)= inddel
          allx_ymat(:,:,i6)= x_ymat

          theta_difff= 1.d0
          do while(theta_difff> 0.01d0)
            theta_old= beta_old2

!           CALL hest(beta_old2, delta, hall, hnew, k, n, rfix, x_ymat, z)

!===============  H ESTIMATION STEP  ==========================


      xl= -20.d0
      xu= 20.d0
      xmid= 0.d0
      bz= z(:, 1)*beta_old2(1)+z(:, 2)*beta_old2(2)
  
  
      err1= 1
      do while(err1> 0.01d0)
 
        fxl1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xl))/rfix)-1
        fxu1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xu))/rfix)-1

        if (fxl1*fxu1.lt.0.d0) then
          xmid= (xl+xu)/2.d0
          fxmid1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xmid))/rfix)-1
          if(fxmid1.lt.0.d0) then
            xl= xmid
          else
            xu= xmid
          end if

        else
!          print*, 'try different extremes 1'
          exit

        end if
        err1= abs(fxmid1)

      end do
    
      hnew(1)= xmid
  
      do i1= 2, k
        xl= -15.d0
        xu= 15.d0
        xmid= 0.d0

        err2= 1.d0
        term1=1+sum(x_ymat(:,i1)*log(1+rfix*exp(bz+hnew(i1-1)))/rfix)

        do while(err2> 0.01d0)

          fxl= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xl))/rfix)-term1
          fxu= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xu))/rfix)-term1

          if (fxl*fxu.lt.0d0) then
            xmid= (xl+xu)/2.d0
            fxmid= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xmid))/rfix)-term1
            if(fxmid.lt.0d0) then
              xl= xmid
            else
              xu= xmid
            end if

          else
!            print*, 'try different extremes 2'
            exit
          end if
          err2= abs(fxmid)

        end do
        hnew(i1)= xmid
      end do
  
      hall= -50.d0

      i5= 0
      do i4= 1, n
        if(delta(i4).eq.1) then
          i5= i5+1
          delone(i5)= i4
        end if
      end do


      i3= 0
      do i2= minval(delone), n
        if(delta(i2).eq.1) then
          i3= i3+1
          hall(i2)= hnew(i3)
      
        else
          hall(i2)= hnew(i3)
        end if
      end do

!===============   H ESTIMATION ENDS HERE  =====================


      CALL betaest(beta_old2, delta, hall, k, n, rfix, z)
    
      theta_difff= sum(abs((beta_old2-theta_old)/theta_old))

    end do   !  end of while loop for theta_difff
    
    allbeta(:,i6)= beta_old2
    allhall(:,i6)= hall
    allh(:,i6)= hnew

    CALL newhall(hall, ltvec, n, outvec, tvec, x_ik_sort)

    newhsum= newhsum + outvec
    beta_old_sum= beta_old_sum + beta_old2

  end do    ! end of for loop i6= 1, m

 
  beta_old_mean= beta_old_sum/m
  newhavg= newhsum/m

!============  BETA AND H ESTIMATION ENDS HERE  ==================


        beta_old= beta_old_mean

        CALL my_match(left_tvec_match, n, ltvec, left, tvec)
!!      output: left_tvec_match

        hallleft= newhavg(left_tvec_match)
        hallold= newhavg(orderforhold)
 
        beta_difff= sum(abs((beta_old-beta_older)/beta_older))

      end do   ! end of while loop for beta convergence

      
      
!!========================================
!!            Sigma Computation
!!========================================

     kminus1= k-1
     kminus2= k-2

     do kk= 1, m
        ymat= allx_ymat(:,:,kk)
        myfinalbeta= allbeta(:,kk)
        myfinalh= allh(:,kk)
        myfinalhall= allhall(:,kk)
        ind_del= allind(:,kk)
        z= allz(:,:,kk)

        hdiff= myfinalh(2:k)-myfinalh(1:(k-1))
        mybzvec= z(:,1)*myfinalbeta(1)+z(:,2)*myfinalbeta(2)

        outh1= exp(mybzvec+myfinalhall)/(1+rfix*exp(mybzvec+myfinalhall))

        do k1= 1, k
          mybzplushmat(:,k1)= exp(mybzvec+myfinalh(k1))
        end do


        sig1=0.d0
        sig2=0.d0

        firstel1= 0.d0
        secondel1= 0.d0
        fourthel1= 0.d0
        
        firstel2= 0.d0
        secondel2= 0.d0
        thirdel2= 0.d0
        fourthel2= 0.d0

        do t1= 2, k

          sum_int= 0.d0
          bfnvec= 1.d0

          do t2= t1, kminus2
          
            outh2= mybzplushmat(:,(t2+1))/(1+rfix*mybzplushmat(:,(t2+1)))
            hazy= outh2*ymat(:,(t2+1))
            outh3= mybzplushmat(:,(t2+1))/((1+rfix*mybzplushmat(:,(t2+1)))*(1+rfix*mybzplushmat(:,(t2+1))))
            dhazy= outh3*ymat(:,(t2+1))

            integrand= sum(dhazy)*hdiff(t2)/sum(hazy)
            sum_int= sum_int + integrand

            bfnindx1= ind_del(t2+1)
            bfnindx2= ind_del(t2+2)-1
            bfnvec(bfnindx1:bfnindx2)= exp(-sum_int)

          end do

          bfnindx1= ind_del(k)

          outh2= mybzplushmat(:,k)/(1+rfix*mybzplushmat(:,k))
          hazy= outh2*ymat(:,k)
          outh3= mybzplushmat(:,k)/((1+rfix*mybzplushmat(:,k))*(1+rfix*mybzplushmat(:,k)))
          dhazy= outh3*ymat(:,k)

          integrand= sum(dhazy)*hdiff(kminus1)/sum(hazy)
          sum_int= sum_int + integrand

          bfnvec(bfnindx1:n)= exp(-sum_int)

          if(t1== kminus1) then
            bfnindx1= ind_del(kminus1)
            bfnindx2= ind_del(k)-1
            bfnvec(bfnindx1:bfnindx2)= 1.d0
 
             
            outh2= mybzplushmat(:,k)/(1+rfix*mybzplushmat(:,k))
            hazy= outh2*ymat(:,k)
            outh3= mybzplushmat(:,k)/((1+rfix*mybzplushmat(:,k))*(1+rfix*mybzplushmat(:,k)))
            dhazy= outh3*ymat(:,k)
            
            integrand= sum(dhazy)*hdiff(kminus1)/sum(hazy)
            sum_int= sum_int + integrand

            bfnvec(ind_del(k):n)= exp(-sum_int)
                     
          else if(t1==k) then
            bfnvec(ind_del(k):n)= 1.d0
            
          end if

         
          outh2= mybzplushmat(:,t1)/(1+rfix*mybzplushmat(:,t1))
          hazy= outh2*ymat(:,t1)

          outh3= mybzplushmat(:,t1)/((1+rfix*mybzplushmat(:,t1))*(1+rfix*mybzplushmat(:,t1)))
          dhazy= outh3*ymat(:,t1)
      
          oyb= outh1*ymat(:,t1)*bfnvec
          zbarvec(1)= sum(z(:,1)*oyb)/sum(hazy)
          zbarvec(2)= sum(z(:,2)*oyb)/sum(hazy)

          zminuszbar(:,1)= z(:,1)-zbarvec(1)
          zminuszbar(:,2)= z(:,2)-zbarvec(2)

          hazyhdiff= hazy*hdiff(t1-1)
          dhazyhdiff= dhazy*hdiff(t1-1)

          firstel1= firstel1+zminuszbar(:,1)*zminuszbar(:,1)*hazyhdiff
          secondel1= secondel1+zminuszbar(:,1)*zminuszbar(:,2)*hazyhdiff
          fourthel1= fourthel1+zminuszbar(:,2)*zminuszbar(:,2)*hazyhdiff

          firstel2= firstel2+zminuszbar(:,1)*z(:,1)*dhazyhdiff
          secondel2= secondel2+zminuszbar(:,2)*z(:,1)*dhazyhdiff
          thirdel2= thirdel2+zminuszbar(:,1)*z(:,2)*dhazyhdiff
          fourthel2= fourthel2+zminuszbar(:,2)*z(:,2)*dhazyhdiff

        end do    ! end of for loop t1= 2, k
   
        sig1= (/sum(firstel1), sum(secondel1), sum(secondel1), sum(fourthel1)/)
        sig2= (/sum(firstel2), sum(secondel2), sum(thirdel2), sum(fourthel2)/)


        sigustar= sig1/n
        siglstar= sig2/n


        det_siglstar= 1.d0/(siglstar(1)*siglstar(4)-siglstar(2)*siglstar(3))
        adj_siglstar= (/siglstar(4), -siglstar(2), -siglstar(3), siglstar(1)/)
        inv_siglstar= adj_siglstar*det_siglstar
        
        varbeta(1)= inv_siglstar(1)*inv_siglstar(1)*sigustar(1) +&
                      inv_siglstar(1)*inv_siglstar(3)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(3)*inv_siglstar(3)*sigustar(4)
        varbeta(2)= inv_siglstar(2)*inv_siglstar(2)*sigustar(1) +&
                      inv_siglstar(2)*inv_siglstar(4)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(4)*inv_siglstar(4)*sigustar(4)


        avhid= m*(tt-1)+kk
        allvar_hope(1,avhid)= varbeta(1)/n
        allvar_hope(2,avhid)= varbeta(2)/n

        allbeta_hope(:,avhid)= allbeta(:,kk)

      end do      ! end of do loop for kk= 1, m

    end do        ! end of do loop for theta_rep

!    print*, "allbeta_hope(1,:)= ", allbeta_hope(1,:)
!    print*, "allbeta_hope(2,:)= ", allbeta_hope(2,:)

    beta_neww(1)= sum(allbeta_hope(1,:))/mtr
    beta_neww(2)= sum(allbeta_hope(2,:))/mtr


    ahvar_mean(1)= sum(allvar_hope(1,:))/mtr
    ahvar_mean(2)= sum(allvar_hope(2,:))/mtr

    beta_hope= allbeta_hope - spread(beta_neww,2,mtr)
    beta_hope_beta_hope_t(1)= sum(beta_hope(1,:)*beta_hope(1,:))
    beta_hope_beta_hope_t(2)= sum(beta_hope(2,:)*beta_hope(2,:))

    sig_hope= ahvar_mean + (mtr+1)*beta_hope_beta_hope_t/(mtr*(mtr-1))
    mysd= SQRT(sig_hope)


!print*, 'final ans= ', (/beta_neww,mysd/)! 
!print*, (/beta_neww, mysd/)

return
end subroutine






!========================================================================
! THIS FUNCTION IS USED IN THE BETA AND H ESTIMATION PART
! this function extends H estimates over points where H was not estimated
! input: a vector of SORTED t values where H has been estimated (xvec)
!        a vector of the corresponding H estimates (Hvec)
!        a vector over which the H estimates are to be expanded using the continuity property (tvec)
!        ltvec= length(tvec)
!        n= length(xvec)
! output: a vector of H estimates at all t-values included in tvec (outvec)

subroutine newhall(hvec, ltvec, n, outvec, tvec, xvec)
implicit none

integer(4):: i, ltvec, n
integer(4):: myord(n)
real(8):: hvec(n), minout, outvec(ltvec), tvec(ltvec), xvec(n)

  outvec(:)= 0.d0

  CALL my_match(myord, n, ltvec, xvec, tvec)
  do i= 1, n-1
    outvec(myord(i):(myord(i+1)-1))= hvec(i)
  end do

  outvec(myord(n):ltvec)= maxval(hvec)

  minout= minval(pack(outvec, outvec/=0.d0))
  where(outvec==0.d0)  outvec= minout-abs(minout)
return
end subroutine



!=============================================================
! MATCH FUNCTION
! vec1 is a vector whose position in vec2 is returned by outvec

subroutine my_match(outvec, s1, s2, vec1, vec2)
implicit none

integer(4):: i1, i2, s1, s2
integer(4):: outvec(s1)
real(8):: vec1(s1), vec2(s2)

  outvec(:)= 0.d0
  do i2= 1, s1
    do i1= 1, s2
      if(vec1(i2)==vec2(i1)) then
        outvec(i2)= i1
        exit
      end if
    end do
  end do
  
return
end subroutine



! MATCH FUNCTION USED IN IMPUTATION OF t
! val is a scalar value whose position in vec2 is returned by outval

subroutine my_match_val(outval, s, val, vec2)
implicit none

integer(4):: i1, s
integer(4):: outval
real(8):: val, vec2(s)

  do i1= 1, s
      if(val==vec2(i1)) then
        outval= i1
        exit
      end if
  end do
return
end subroutine


!==========================================
!RETURNS P, Q IN ASCENDING ORDER
! used in bubble sorting

subroutine order(p,q)
implicit none

real(8):: p, q, temp

  if (p>q) then
    temp=p
    p=q
    q=temp
  end if

return
end subroutine


!=========================================
!BUBBLE SORT

subroutine bubble(vec, n)
implicit none

integer(4):: i, j, n
real(8):: vec(n), temp

  do i=1, n
    do j=n, i+1, -1
!      CALL order(vec(j-1), vec(j))

      if (vec(j-1)>vec(j)) then
        temp=vec(j-1)
        vec(j-1)=vec(j)
        vec(j)=temp
      end if

    end do
  end do
  
return
end subroutine


!==============================================
! WITH REPLACEMENT SAMPLING
! input: length of the weight vector (s), weight vector (weight)


function sample_wr(s, weight)

integer(4):: j1, j2, s, sample_wr
real(8):: cumsum(s), weight(s), xx

integer(4):: is1, is2
double precision uniform

common /unif_seeds/ is1, is2
save /unif_seeds/

!print*, isd

  do j1= 1, s
    cumsum(j1)= sum(weight(1:j1))/sum(weight)
  end do


  xx= uniform()
    do j2= 1, s
      if(xx .lt. cumsum(j2)) then
        sample_wr= j2
        exit
      end if
    end do
  
end function



!==============================================================
!                     NAIVE METHOD for r=0
!==============================================================


subroutine nv_r0(beta_neww, beta_old, cumsumgrid,&
                       delta_temp, final_x, hallleft,&
                       hallold, isd, k, left,&
                       lrcpos, ltvec,&
                       m, maxgrid, mysd, n,&
                       notrcpos, numgrid,&
                       orderforhold, rcpos, rfix,&
                       totngrid, ttilmat, tvec, z2)
implicit none


integer(4):: avhid, ii, j, j2, k, k1, kk, kminus1, kminus2,&
            lrcpos, ltvec, m, maxgrid, n, t1, t2, totngrid, tt,&
            sm, dellater, dellater1, dellater2
            
integer(4):: bfnindx1, bfnindx2, cumsumgrid(k), delta_temp(n),&
            left_tvec_match(n), notrcpos(k), numgrid(k),&
            orderforhold(totngrid), rcpos(lrcpos)

real(8):: beta_difff, beta_neww(2), beta_old(2), final_x(n),&
         left(n), mysd(2), rfix, ttilmat(totngrid),&
         tvec(ltvec), z2(n)

real(8):: beta_older(2), beta_old_mean(2), hallleft(n),&
         hallold(totngrid), init_beta(2), newhavg(ltvec),&
         newz1(n), x_ik(n,m), z(n,2)


integer(4):: cc, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10,&
            allind(k,m), delone(k), delta(n), inddel(k), ind_del(k),&
            this, x_ik_sort_indx(n)

real(8):: allbeta(2,m), alldel(n,m), allh(k,m), allhall(n,m), allx_ymat(n,k,m),&
         allz(n,2,m), beta_old_sum(2), beta_old2(2), bz(n), hall(n), hnew(k),&
         init_beta2(2), newhsum(ltvec), outvec(ltvec), theta_difff, theta_old(2),&
         x_ik_sort(n), x_t1t2(k), x_ymat(n,k), z_temp2(n,2)

real(8):: adj_siglstar(4), ahvar_mean(2), allbeta_hope(2,m),&
         allvar_hope(2,m), beta_diff(2,m), beta_hope(2,m),&
         beta_hope_beta_hope_t(2), bfnvec(n), bfnvec1(n), det_siglstar,&
         dhazy(n), dhazyhdiff(n), firstel1(n), firstel2(n), fourthel1(n),&
         fourthel2(n), hazy(n), hazyhdiff(n), hdiff(k-1), integrand

real(8):: mn, mybzplush(n), mybzvec(n), myfinalbeta(2), myfinalh(k),&
         myfinalhall(n), outdh(n), outh(n), outh1(n), outh2(n), outh3(n),&
         oyb(n), ratio, secondel1(n), secondel2(n), sig1(4),&
         sig2(4), sig_hope(2), sum_int, theta_diff, thirdel2(n),&
         varbeta(2), ymat(n,k), zbarvec(2), zminuszbar(n,2),&
         zzprime1(4), zzprime2(4)

real(8):: siglstar(4), sigustar(4), inv_siglstar(4)

real(8):: err1, err2, fxl1, fxl, fxmid1, fxmid, fxu1,&
         fxu, term1, xl, xmid, xu 

integer(4) :: is1, is2, isd, seed1, seed2        

    common /unif_seeds/ is1, is2
    save /unif_seeds/
    seed1=10005
    seed2=isd      
    CALL set_uniform(seed1, seed2)
  
  
    init_beta= beta_old
    beta_neww= 0.d0
    allbeta_hope= 0.d0
    allvar_hope= 0.d0

!   print*, "entering theta loop"



      beta_old= init_beta
      newz1= final_x
      beta_difff= 1.0d0

      do while(beta_difff > 0.01d0)
      
        beta_older= beta_old
!       print*, "ok till here"

        CALL imp_tr0(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                   lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                   x_ik, rcpos, rfix, totngrid, ttilmat, z2)
                   
!!      output: x_ik


!=============  BETA AND H ESTIMATION STEP ======================


        z_temp2(:,1)= newz1       
        z_temp2(:,2)= z2
  
        beta_old_sum= 0.d0
        newhsum= 0.d0
        init_beta2= beta_old


        do i6= 1, m
 
          beta_old2= init_beta2
          hall= 0.d0
  
          x_ik_sort= x_ik(:,i6)
          CALL bubble(x_ik_sort, n)           ! equivalent of sort()$x in R
          CALL my_match(x_ik_sort_indx, n, n, x_ik_sort, x_ik(:,i6))
                                              ! equivalent of sort()$ix in R
    
          z= z_temp2(x_ik_sort_indx,:)
          delta= delta_temp(x_ik_sort_indx)

          allz(:,:,i6)= z
          alldel(:,i6)= delta

!          print*, "step 1 done"

          cc=1
          do i8= 1, n
            if(delta(i8)==1.d0) then
              x_t1t2(cc)= x_ik_sort(i8)
              inddel(cc)= i8
              cc= cc+1
            end if
          end do
    
          x_ymat=0.d0    
          do i7= 1, k
            where(x_ik_sort>= x_t1t2(i7)) x_ymat(:,i7)= 1.d0
          end do
    
          allind(:, i6)= inddel

      allx_ymat(:,:,i6)= x_ymat

      theta_difff= 1.d0
      do while(theta_difff> 0.01d0)
        theta_old= beta_old2

!       CALL hest(beta_old2, delta, hall, hnew, k, n, rfix, x_ymat, z)

!===============  H ESTIMATION STEP  ==========================


      xl= -20.d0
      xu= 20.d0
      xmid= 0.d0
      bz= z(:, 1)*beta_old2(1)+z(:, 2)*beta_old2(2)
  
  
      err1= 1
      do while(err1> 0.01d0)
 
        fxl1= sum(x_ymat(:,1)*exp(bz+xl))-1
        fxu1= sum(x_ymat(:,1)*exp(bz+xu))-1

        if (fxl1*fxu1.lt.0.d0) then
          xmid= (xl+xu)/2.d0
          fxmid1= sum(x_ymat(:,1)*exp(bz+xmid))-1
          if(fxmid1.lt.0.d0) then
            xl= xmid
          else
            xu= xmid
          end if

        else
 !         print*, 'try different extremes 1'
          exit

        end if
        err1= abs(fxmid1)

      end do
    
      hnew(1)= xmid
  
      do i1= 2, k
        xl= -15.d0
        xu= 15.d0
        xmid= 0.d0

        err2= 1.d0
        term1=1+sum(x_ymat(:,i1)*exp(bz+hnew(i1-1)))

        do while(err2> 0.01d0)

          fxl= sum(x_ymat(:,i1)*exp(bz+xl))-term1
          fxu= sum(x_ymat(:,i1)*exp(bz+xu))-term1

          if (fxl*fxu.lt.0d0) then
            xmid= (xl+xu)/2.d0
            fxmid= sum(x_ymat(:,i1)*exp(bz+xmid))-term1
            if(fxmid.lt.0d0) then
              xl= xmid
            else
              xu= xmid
            end if

          else
  !          print*, 'try different extremes 2'
            exit
          end if
          err2= abs(fxmid)

        end do
        hnew(i1)= xmid
      end do
  
      hall= -50.d0

      i5= 0
      do i4= 1, n
        if(delta(i4).eq.1) then
          i5= i5+1
          delone(i5)= i4
        end if
      end do


      i3= 0
      do i2= minval(delone), n
        if(delta(i2).eq.1) then
          i3= i3+1
          hall(i2)= hnew(i3)
      
        else
          hall(i2)= hnew(i3)
        end if
      end do

!===============   H ESTIMATION ENDS HERE  =====================


      CALL betaestr0(beta_old2, delta, hall, k, n, rfix, z)
    
      theta_difff= sum(abs((beta_old2-theta_old)/theta_old))

!      print*, 'beta_old2= ', beta_old2
      
    end do   !  end of while loop for theta_difff
    
    allbeta(:,i6)= beta_old2
    allhall(:,i6)= hall
    allh(:,i6)= hnew

    CALL newhallr0(hall, ltvec, n, outvec, tvec, x_ik_sort)

    newhsum= newhsum + outvec
    beta_old_sum= beta_old_sum + beta_old2

!    print*, "step 4 done" 

  end do    ! end of for loop i6= 1, m

 
  beta_old_mean= beta_old_sum/m
  newhavg= newhsum/m

!============  BETA AND H ESTIMATION ENDS HERE  =====================

!       print*, "betah_est is fine"

        beta_old= beta_old_mean

        CALL my_match(left_tvec_match, n, ltvec, left, tvec)
!!      output: left_tvec_match

        hallleft= newhavg(left_tvec_match)
        hallold= newhavg(orderforhold)
 
        beta_difff= sum(abs((beta_old-beta_older)/beta_older))

      end do   ! end of while loop for beta convergence

      
!     print*, 'final_beta_theta= ', beta_old
!     beta_neww= beta_neww+beta_old

     
      
      
      
!!========================================
!!            Sigma Computation
!!========================================

     do kk= 1, m
        ymat= allx_ymat(:,:,kk)
        myfinalbeta= allbeta(:,kk)
        myfinalh= allh(:,kk)
        myfinalhall= allhall(:,kk)
        z= allz(:,:,kk)
            
        hdiff= myfinalh(2:k)-myfinalh(1:(k-1))
        mybzvec= z(:,1)*myfinalbeta(1)+z(:,2)*myfinalbeta(2)

        outh1= exp(mybzvec+myfinalhall)

        sig1=0.d0
        sig2=0.d0

        firstel1= 0.d0
        secondel1= 0.d0
        fourthel1= 0.d0
        
        firstel2= 0.d0
        secondel2= 0.d0
        thirdel2= 0.d0
        fourthel2= 0.d0

        do t1= 2, k
          bfnvec= exp(myfinalh(t1)-myfinalhall)
          mybzplush= mybzvec+myfinalh(t1)
         
          outh2= exp(mybzplush)
          hazy= outh2*ymat(:,t1)

          dhazy= hazy
     
          oyb= outh1*ymat(:,t1)*bfnvec
          zbarvec(1)= sum(z(:,1)*oyb)/sum(hazy)
          zbarvec(2)= sum(z(:,2)*oyb)/sum(hazy)

          zminuszbar(:,1)= z(:,1)-zbarvec(1)
          zminuszbar(:,2)= z(:,2)-zbarvec(2)

          hazyhdiff= hazy*hdiff(t1-1)
          dhazyhdiff= dhazy*hdiff(t1-1)

          firstel1= firstel1+zminuszbar(:,1)*zminuszbar(:,1)*hazyhdiff
          secondel1= secondel1+zminuszbar(:,1)*zminuszbar(:,2)*hazyhdiff
          fourthel1= fourthel1+zminuszbar(:,2)*zminuszbar(:,2)*hazyhdiff

          firstel2= firstel2+zminuszbar(:,1)*z(:,1)*dhazyhdiff
          secondel2= secondel2+zminuszbar(:,2)*z(:,1)*dhazyhdiff
          thirdel2= thirdel2+zminuszbar(:,1)*z(:,2)*dhazyhdiff
          fourthel2= fourthel2+zminuszbar(:,2)*z(:,2)*dhazyhdiff

        end do  ! end of for loop t1= 2, k
   
        sig1= (/sum(firstel1), sum(secondel1), sum(secondel1), sum(fourthel1)/)
        sig2= (/sum(firstel2), sum(secondel2), sum(thirdel2), sum(fourthel2)/)

        sigustar= sig1/n
        siglstar= sig2/n


        det_siglstar= 1.d0/(siglstar(1)*siglstar(4)-siglstar(2)*siglstar(3))
        adj_siglstar= (/siglstar(4), -siglstar(2), -siglstar(3), siglstar(1)/)
        inv_siglstar= adj_siglstar*det_siglstar
        
        varbeta(1)= inv_siglstar(1)*inv_siglstar(1)*sigustar(1) +&
                      inv_siglstar(1)*inv_siglstar(3)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(3)*inv_siglstar(3)*sigustar(4)
        varbeta(2)= inv_siglstar(2)*inv_siglstar(2)*sigustar(1) +&
                      inv_siglstar(2)*inv_siglstar(4)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(4)*inv_siglstar(4)*sigustar(4)


        avhid= kk
        allvar_hope(1,avhid)= varbeta(1)/n
        allvar_hope(2,avhid)= varbeta(2)/n

        allbeta_hope(:,avhid)= allbeta(:,kk)

      end do      ! end of do loop for kk= 1, m


!    print*, "allbeta_hope(1,:)= ", allbeta_hope(1,:)
!    print*, "allbeta_hope(2,:)= ", allbeta_hope(2,:)

    beta_neww(1)= sum(allbeta_hope(1,:))/m
    beta_neww(2)= sum(allbeta_hope(2,:))/m


    ahvar_mean(1)= sum(allvar_hope(1,:))/m
    ahvar_mean(2)= sum(allvar_hope(2,:))/m

    beta_hope= allbeta_hope - spread(beta_neww,2,m)
    beta_hope_beta_hope_t(1)= sum(beta_hope(1,:)*beta_hope(1,:))
    beta_hope_beta_hope_t(2)= sum(beta_hope(2,:)*beta_hope(2,:))

    sig_hope= ahvar_mean + (m+1)*beta_hope_beta_hope_t/(m*(m-1))
    mysd= SQRT(sig_hope)


!print*, 'final ans= ', (/beta_neww,mysd/)
!print*, (/beta_neww,mysd/)

return
end subroutine






!#################################################################################
!                           NAIVE METHOD general r case
!#################################################################################

subroutine nv_genr(beta_neww, beta_old, cumsumgrid,&
                       delta_temp, final_x, hallleft,&
                       hallold, isd, k, left,&
                       lrcpos, ltvec, m,&
                       maxgrid, mysd, n,&
                       notrcpos, numgrid,&
                       orderforhold, rcpos, rfix,&
                       totngrid, ttilmat, tvec, z2)
implicit none


integer(4):: avhid, ii, j, j2, k, k1, kk, kminus1, kminus2,&
            lrcpos, ltvec, m, maxgrid, n, t1, t2, totngrid,&
            sm, dellater, dellater1, dellater2
            
integer(4):: bfnindx1, bfnindx2, cumsumgrid(k), delta_temp(n),&
            left_tvec_match(n), notrcpos(k), numgrid(k),&
            orderforhold(totngrid), rcpos(lrcpos)

real(8):: beta_difff, beta_neww(2), beta_old(2), final_x(n),&
         left(n), mysd(2), rfix, ttilmat(totngrid),&
         tvec(ltvec), z2(n)

real(8):: beta_older(2), beta_old_mean(2), hallleft(n),&
         hallold(totngrid), init_beta(2), newhavg(ltvec),&
         newz1(n), x_ik(n,m), z(n,2)


integer(4):: cc, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10,&
            allind(k,m), delone(k), delta(n), inddel(k), ind_del(k),&
            this, x_ik_sort_indx(n)

real(8):: allbeta(2,m), alldel(n,m), allh(k,m), allhall(n,m), allx_ymat(n,k,m),&
         allz(n,2,m), beta_old_sum(2), beta_old2(2), bz(n), hall(n), hnew(k),&
         init_beta2(2), newhsum(ltvec), outvec(ltvec), theta_difff, theta_old(2),&
         x_ik_sort(n), x_t1t2(k), x_ymat(n,k), z_temp2(n,2)

real(8):: adj_siglstar(4), ahvar_mean(2), allbeta_hope(2,m),&
         allvar_hope(2,m), beta_diff(2,m), beta_hope(2,m),&
         beta_hope_beta_hope_t(2), bfnvec(n), bfnvec1(n), det_siglstar,&
         dhazy(n), dhazyhdiff(n), firstel1(n), firstel2(n), fourthel1(n),&
         fourthel2(n), hazy(n), hazyhdiff(n), hdiff(k-1), integrand

real(8):: mn, mybzplushmat(n,k), mybzvec(n), myfinalbeta(2), myfinalh(k),&
         myfinalhall(n), outdh(n), outh(n), outh1(n), outh2(n), outh3(n),&
         oyb(n), ratio, secondel1(n), secondel2(n), sig1(4),&
         sig2(4), sig_hope(2), sum_int, theta_diff, thirdel2(n),&
         varbeta(2), ymat(n,k), zbarvec(2), zminuszbar(n,2),&
         zzprime1(4), zzprime2(4)

real(8):: siglstar(4), sigustar(4), inv_siglstar(4)

real(8):: err1, err2, fxl1, fxl, fxmid1, fxmid, fxu1,&
         fxu, term1, xl, xmid, xu 

integer(4) :: is1, is2, isd, seed1, seed2        

    common /unif_seeds/ is1, is2
    save /unif_seeds/
    seed1=10005
    seed2=isd      
    CALL set_uniform(seed1, seed2)
  
  
    init_beta= beta_old
    beta_neww= 0.d0
    allbeta_hope= 0.d0
    allvar_hope= 0.d0

!   print*, "entering theta loop"


      beta_old= init_beta
      newz1= final_x
      beta_difff= 1.0d0

      do while(beta_difff > 0.01d0)
      
        beta_older= beta_old
        
        CALL imp_t(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                   lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                   x_ik, rcpos, rfix, totngrid, ttilmat, z2)


!!      output: x_ik


!=============  BETA AND H ESTIMATION STEP ======================


        z_temp2(:,1)= newz1       
        z_temp2(:,2)= z2
  
        beta_old_sum= 0.d0
        newhsum= 0.d0
        init_beta2= beta_old


        do i6= 1, m

          beta_old2= init_beta2
          hall= 0.d0
  
          x_ik_sort= x_ik(:,i6)
          CALL bubble(x_ik_sort, n)           ! equivalent of sort()$x in R
          CALL my_match(x_ik_sort_indx, n, n, x_ik_sort, x_ik(:,i6))
                                              ! equivalent of sort()$ix in R
    
          z= z_temp2(x_ik_sort_indx,:)
          delta= delta_temp(x_ik_sort_indx)

          allz(:,:,i6)= z
          alldel(:,i6)= delta

!          print*, "step 1 done"

          cc=1
          do i8= 1, n
            if(delta(i8)==1.d0) then
              x_t1t2(cc)= x_ik_sort(i8)
              inddel(cc)= i8
              cc= cc+1
            end if
          end do
    
          x_ymat=0.d0    
          do i7= 1, k
            where(x_ik_sort>= x_t1t2(i7)) x_ymat(:,i7)= 1.d0
          end do
    
          allind(:, i6)= inddel

      allx_ymat(:,:,i6)= x_ymat

!      print*, "step 3 done"

      theta_difff= 1.d0
      do while(theta_difff> 0.01d0)
        theta_old= beta_old2

!       CALL hest(beta_old2, delta, hall, hnew, k, n, rfix, x_ymat, z)

!===============  H ESTIMATION STEP  ==========================


      xl= -20.d0
      xu= 20.d0
      xmid= 0.d0
      bz= z(:, 1)*beta_old2(1)+z(:, 2)*beta_old2(2)
  
  
      err1= 1
      do while(err1> 0.01d0)
 
        fxl1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xl))/rfix)-1
        fxu1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xu))/rfix)-1

        if (fxl1*fxu1.lt.0.d0) then
          xmid= (xl+xu)/2.d0
          fxmid1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xmid))/rfix)-1
          if(fxmid1.lt.0.d0) then
            xl= xmid
          else
            xu= xmid
          end if

        else
 !         print*, 'try different extremes 1'
          exit

        end if
        err1= abs(fxmid1)

      end do
    
      hnew(1)= xmid
  
      do i1= 2, k
        xl= -15.d0
        xu= 15.d0
        xmid= 0.d0

        err2= 1.d0
        term1=1+sum(x_ymat(:,i1)*log(1+rfix*exp(bz+hnew(i1-1)))/rfix)

        do while(err2> 0.01d0)

          fxl= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xl))/rfix)-term1
          fxu= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xu))/rfix)-term1

          if (fxl*fxu.lt.0d0) then
            xmid= (xl+xu)/2.d0
            fxmid= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xmid))/rfix)-term1
            if(fxmid.lt.0d0) then
              xl= xmid
            else
              xu= xmid
            end if

          else
  !          print*, 'try different extremes 2'
            exit
          end if
          err2= abs(fxmid)

        end do
        hnew(i1)= xmid
      end do
  
      hall= -50.d0

      i5= 0
      do i4= 1, n
        if(delta(i4).eq.1) then
          i5= i5+1
          delone(i5)= i4
        end if
      end do


      i3= 0
      do i2= minval(delone), n
        if(delta(i2).eq.1) then
          i3= i3+1
          hall(i2)= hnew(i3)
      
        else
          hall(i2)= hnew(i3)
        end if
      end do

!===============   H ESTIMATION ENDS HERE  =====================


      CALL betaest(beta_old2, delta, hall, k, n, rfix, z)
    
      theta_difff= sum(abs((beta_old2-theta_old)/theta_old))

!     print*, 'beta_old2= ', beta_old2

     
    end do   !  end of while loop for theta_difff
    
    allbeta(:,i6)= beta_old2
    allhall(:,i6)= hall
    allh(:,i6)= hnew

    CALL newhall(hall, ltvec, n, outvec, tvec, x_ik_sort)

    newhsum= newhsum + outvec
    beta_old_sum= beta_old_sum + beta_old2

!   print*, "step 4 done" 

  end do    ! end of for loop i6= 1, m

 
  beta_old_mean= beta_old_sum/m
  newhavg= newhsum/m

!============  BETA AND H ESTIMATION ENDS HERE  =====================

!       print*, "betah_est is fine"

        beta_old= beta_old_mean

        CALL my_match(left_tvec_match, n, ltvec, left, tvec)
!!      output: left_tvec_match

        hallleft= newhavg(left_tvec_match)
        hallold= newhavg(orderforhold)
 
        beta_difff= sum(abs((beta_old-beta_older)/beta_older))

      end do   ! end of while loop for beta convergence

      
!     print*, 'final_beta_theta= ', beta_old
!     beta_neww= beta_neww+beta_old

     
      
      
      
!!========================================
!!            Sigma Computation
!!========================================

     kminus1= k-1
     kminus2= k-2

     do kk= 1, m
        ymat= allx_ymat(:,:,kk)
        myfinalbeta= allbeta(:,kk)
        myfinalh= allh(:,kk)
        myfinalhall= allhall(:,kk)
        ind_del= allind(:,kk)
        z= allz(:,:,kk)

        hdiff= myfinalh(2:k)-myfinalh(1:(k-1))
        mybzvec= z(:,1)*myfinalbeta(1)+z(:,2)*myfinalbeta(2)

        outh1= exp(mybzvec+myfinalhall)/(1+rfix*exp(mybzvec+myfinalhall))

        do k1= 1, k
          mybzplushmat(:,k1)= exp(mybzvec+myfinalh(k1))
        end do


        sig1=0.d0
        sig2=0.d0

        firstel1= 0.d0
        secondel1= 0.d0
        fourthel1= 0.d0
        
        firstel2= 0.d0
        secondel2= 0.d0
        thirdel2= 0.d0
        fourthel2= 0.d0

        do t1= 2, k

          sum_int= 0.d0
          bfnvec= 1.d0

          do t2= t1, kminus2
          
            outh2= mybzplushmat(:,(t2+1))/(1+rfix*mybzplushmat(:,(t2+1)))
            hazy= outh2*ymat(:,(t2+1))
            outh3= mybzplushmat(:,(t2+1))/((1+rfix*mybzplushmat(:,(t2+1)))*(1+rfix*mybzplushmat(:,(t2+1))))
            dhazy= outh3*ymat(:,(t2+1))

            integrand= sum(dhazy)*hdiff(t2)/sum(hazy)
            sum_int= sum_int + integrand

            bfnindx1= ind_del(t2+1)
            bfnindx2= ind_del(t2+2)-1
            bfnvec(bfnindx1:bfnindx2)= exp(-sum_int)

          end do

          bfnindx1= ind_del(k)

          outh2= mybzplushmat(:,k)/(1+rfix*mybzplushmat(:,k))
          hazy= outh2*ymat(:,k)
          outh3= mybzplushmat(:,k)/((1+rfix*mybzplushmat(:,k))*(1+rfix*mybzplushmat(:,k)))
          dhazy= outh3*ymat(:,k)

          integrand= sum(dhazy)*hdiff(kminus1)/sum(hazy)
          sum_int= sum_int + integrand

          bfnvec(bfnindx1:n)= exp(-sum_int)

          if(t1== kminus1) then
            bfnindx1= ind_del(kminus1)
            bfnindx2= ind_del(k)-1
            bfnvec(bfnindx1:bfnindx2)= 1.d0
 
             
            outh2= mybzplushmat(:,k)/(1+rfix*mybzplushmat(:,k))
            hazy= outh2*ymat(:,k)
            outh3= mybzplushmat(:,k)/((1+rfix*mybzplushmat(:,k))*(1+rfix*mybzplushmat(:,k)))
            dhazy= outh3*ymat(:,k)
            
            integrand= sum(dhazy)*hdiff(kminus1)/sum(hazy)
            sum_int= sum_int + integrand

            bfnvec(ind_del(k):n)= exp(-sum_int)
          
          else if(t1==k) then

            bfnvec(ind_del(k):n)= 1.d0

          end if

         
          outh2= mybzplushmat(:,t1)/(1+rfix*mybzplushmat(:,t1))
          hazy= outh2*ymat(:,t1)

          outh3= mybzplushmat(:,t1)/((1+rfix*mybzplushmat(:,t1))*(1+rfix*mybzplushmat(:,t1)))
          dhazy= outh3*ymat(:,t1)
      
          oyb= outh1*ymat(:,t1)*bfnvec
          zbarvec(1)= sum(z(:,1)*oyb)/sum(hazy)
          zbarvec(2)= sum(z(:,2)*oyb)/sum(hazy)

          zminuszbar(:,1)= z(:,1)-zbarvec(1)
          zminuszbar(:,2)= z(:,2)-zbarvec(2)

          hazyhdiff= hazy*hdiff(t1-1)
          dhazyhdiff= dhazy*hdiff(t1-1)

          firstel1= firstel1+zminuszbar(:,1)*zminuszbar(:,1)*hazyhdiff
          secondel1= secondel1+zminuszbar(:,1)*zminuszbar(:,2)*hazyhdiff
          fourthel1= fourthel1+zminuszbar(:,2)*zminuszbar(:,2)*hazyhdiff

          firstel2= firstel2+zminuszbar(:,1)*z(:,1)*dhazyhdiff
          secondel2= secondel2+zminuszbar(:,2)*z(:,1)*dhazyhdiff
          thirdel2= thirdel2+zminuszbar(:,1)*z(:,2)*dhazyhdiff
          fourthel2= fourthel2+zminuszbar(:,2)*z(:,2)*dhazyhdiff

        end do    ! end of for loop t1= 2, k
   
        sig1= (/sum(firstel1), sum(secondel1), sum(secondel1), sum(fourthel1)/)
        sig2= (/sum(firstel2), sum(secondel2), sum(thirdel2), sum(fourthel2)/)


        sigustar= sig1/n
        siglstar= sig2/n


        det_siglstar= 1.d0/(siglstar(1)*siglstar(4)-siglstar(2)*siglstar(3))
        adj_siglstar= (/siglstar(4), -siglstar(2), -siglstar(3), siglstar(1)/)
        inv_siglstar= adj_siglstar*det_siglstar
        
        varbeta(1)= inv_siglstar(1)*inv_siglstar(1)*sigustar(1) +&
                      inv_siglstar(1)*inv_siglstar(3)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(3)*inv_siglstar(3)*sigustar(4)
        varbeta(2)= inv_siglstar(2)*inv_siglstar(2)*sigustar(1) +&
                      inv_siglstar(2)*inv_siglstar(4)*(sigustar(2)+sigustar(3)) +&
                      inv_siglstar(4)*inv_siglstar(4)*sigustar(4)


        avhid= kk
        allvar_hope(1,avhid)= varbeta(1)/n
        allvar_hope(2,avhid)= varbeta(2)/n

        allbeta_hope(:,avhid)= allbeta(:,kk)

      end do      ! end of do loop for kk= 1, m

    beta_neww(1)= sum(allbeta_hope(1,:))/m
    beta_neww(2)= sum(allbeta_hope(2,:))/m


    ahvar_mean(1)= sum(allvar_hope(1,:))/m
    ahvar_mean(2)= sum(allvar_hope(2,:))/m

    beta_hope= allbeta_hope - spread(beta_neww,2,m)
    beta_hope_beta_hope_t(1)= sum(beta_hope(1,:)*beta_hope(1,:))
    beta_hope_beta_hope_t(2)= sum(beta_hope(2,:)*beta_hope(2,:))

    sig_hope= ahvar_mean + (m+1)*beta_hope_beta_hope_t/(m*(m-1))
    mysd= SQRT(sig_hope)


!print*, 'final ans= ', (/beta_neww,mysd/)
!print*, (/beta_neww, mysd/)

return
end subroutine




!==========================================================================
!                    REGRESSION CALIBRATION with r=0
!==========================================================================

subroutine rc_r0(beta_neww, beta_old, cumsumgrid,&
                       delta_temp, final_x, hallleft,&
                       hallold, isd, k, left,&
                       lrcpos, ltvec,&
                       m, maxgrid, mysd, n,&
                       notrcpos, numgrid,&
                       orderforhold, rcpos, rfix,&
                       sigma_alpha, totngrid,&
                       ttilmat, tvec, wbar, z2)
implicit none


integer(4):: avhid, ii, j, j2, k, k1, kk, kminus1, kminus2,&
            lrcpos, ltvec, m, maxgrid, n, t1, t2, t3, t4,&
            totngrid, tt, sm, dellater, dellater1, dellater2
            
integer(4):: bfnindx1, bfnindx2, cumsumgrid(k), delta_temp(n),&
            left_tvec_match(n), notrcpos(k), numgrid(k),&
            orderforhold(totngrid), rcpos(lrcpos)

real(8):: beta_difff, beta_neww(2), beta_old(2), final_x(n),&
         left(n), mysd(2), rfix, ttilmat(totngrid),&
         tvec(ltvec), wbar(n), z2(n)

real(8):: beta_older(2), beta_old_mean(2), hallleft(n),&
         hallold(totngrid), init_beta(2), newhavg(ltvec),&
         newz1(n), x_ik(n,m), z(n,2)


integer(4):: cc, cc1, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10,&
            allind(k,m), delone(k), delta(n), inddel(k), ind_del(k),&
            t1minus1, this, x_ik_sort_indx(n)

real(8):: allbeta(2,m), alldel(n,m), allh(k,m), allhall(n,m), allx_ymat(n,k,m),&
         allz(n,2,m), beta_old_sum(2), beta_old2(2), bz(n), hall(n), hnew(k),&
         init_beta2(2), newhsum(ltvec), outvec(ltvec), theta_difff, theta_old(2),&
         x_ik_sort(n), x_t1t2(k), x_ymat(n,k), z_temp2(n,2)

real(8):: adj_siglstar(4), ahvar_mean(2), allbeta_hope(2,m),&
         allvar_hope(2,m), beta_diff(2,m), beta_hope(2,m),&
         beta_hope_beta_hope_t(2), bfnvec(n), det_siglstar,&
         dhazy(n), dhazyhdiff(n), firstel1(n), firstel2(n), fourthel1(n),&
         fourthel2(n), hazy(n), hazyhdiff(n), hdiff(k-1), integrand

real(8):: mn, mybzplush(n), mybzvec(n), myfinalbeta(2), myfinalh(k),&
         myfinalhall(n), outdh(n), outh(n), outh1(n), outh2(n), outh3(n),&
         oyb(n), ratio, secondel1(n), secondel2(n), sig1(4),&
         sig2(4), sig_hope(2), sum_int, theta_diff, thirdel2(n),&
         varbeta(2), ymat(n,k), zbarvec(2), zminuszbar(n,2),&
         zzprime1(4), zzprime2(4)

real(8):: siglstar(4), sigustar(4), sigustar_old(4), inv_siglstar(4)

real(8):: err1, err2, fxl1, fxl, fxmid1, fxmid, fxu1,&
         fxu, term1, xl, xmid, xu 

real(8):: cumsum_hdiff, delpart(n,3), dhazy_term1,&
         dhazy_term2, dhazy_term3, dhdalpha1(k-1),&
         dhdalpha2(k-1), dhdalpha3(k-1), dhdalpha_all1(n),&
         dhdalpha_all2(n), dhdalpha_all3(n), dhdalpha_den,&
         dhdalpha_numvec1, dhdalpha_numvec2, dhdalpha_numvec3,&
         dxdalpha(n,3), hazyvec, mydel(n), psi(2,3), psi_siga_psi(2,2),&
         psi_tr(3,2), sigma_alpha(3,3), zstar(2,n)


integer(4) :: is1, is2, isd, seed1, seed2        

    common /unif_seeds/ is1, is2
    save /unif_seeds/
    seed1=10005
    seed2=isd      
    CALL set_uniform(seed1, seed2)
  
  
    init_beta= beta_old
    beta_neww= 0.d0
    allbeta_hope= 0.d0
    allvar_hope= 0.d0



      beta_old= init_beta
      newz1= final_x
      beta_difff= 1.0d0
      
      dxdalpha(:,1)= 1.d0
      dxdalpha(:,2)= wbar
      dxdalpha(:,3)= z2
      
      zstar(1,:)= z2
      zstar(2,:)= final_x


      do while(beta_difff > 0.01d0)
      
        beta_older= beta_old

        CALL imp_tr0(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                   lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                   x_ik, rcpos, rfix, totngrid, ttilmat, z2)

!!      output: x_ik

!=============  BETA AND H ESTIMATION STEP ======================


        z_temp2(:,1)= newz1       
        z_temp2(:,2)= z2
  
        beta_old_sum= 0.d0
        newhsum= 0.d0
        init_beta2= beta_old


        do i6= 1, m

          beta_old2= init_beta2
          hall= 0.d0
  
          x_ik_sort= x_ik(:,i6)
          CALL bubble(x_ik_sort, n)           ! equivalent of sort()$x in R
          CALL my_match(x_ik_sort_indx, n, n, x_ik_sort, x_ik(:,i6))
                                              ! equivalent of sort()$ix in R
    
          z= z_temp2(x_ik_sort_indx,:)
          delta= delta_temp(x_ik_sort_indx)

          allz(:,:,i6)= z
          alldel(:,i6)= delta

!          print*, "step 1 done"

          cc=1
          do i8= 1, n
            if(delta(i8)==1.d0) then
              x_t1t2(cc)= x_ik_sort(i8)
              inddel(cc)= i8
              cc= cc+1
            end if
          end do
    
          x_ymat=0.d0    
          do i7= 1, k
            where(x_ik_sort>= x_t1t2(i7)) x_ymat(:,i7)= 1.d0
          end do

          allind(:, i6)= inddel
          allx_ymat(:,:,i6)= x_ymat

          theta_difff= 1.d0
          do while(theta_difff> 0.01d0)
            theta_old= beta_old2

!       CALL hest(beta_old2, delta, hall, hnew, k, n, rfix, x_ymat, z)

!===============  H ESTIMATION STEP  ==========================


      xl= -20.d0
      xu= 20.d0
      xmid= 0.d0
      bz= z(:, 1)*beta_old2(1)+z(:, 2)*beta_old2(2)
  
  
      err1= 1
      do while(err1> 0.01d0)
 
        fxl1= sum(x_ymat(:,1)*exp(bz+xl))-1
        fxu1= sum(x_ymat(:,1)*exp(bz+xu))-1

        if (fxl1*fxu1.lt.0.d0) then
          xmid= (xl+xu)/2.d0
          fxmid1= sum(x_ymat(:,1)*exp(bz+xmid))-1
          if(fxmid1.lt.0.d0) then
            xl= xmid
          else
            xu= xmid
          end if

        else
!          print*, 'try different extremes 1'
          exit

        end if
        err1= abs(fxmid1)

      end do
    
      hnew(1)= xmid
  
      do i1= 2, k
        xl= -15.d0
        xu= 15.d0
        xmid= 0.d0

        err2= 1.d0
        term1=1+sum(x_ymat(:,i1)*exp(bz+hnew(i1-1)))

        do while(err2> 0.01d0)

          fxl= sum(x_ymat(:,i1)*exp(bz+xl))-term1
          fxu= sum(x_ymat(:,i1)*exp(bz+xu))-term1

          if (fxl*fxu.lt.0d0) then
            xmid= (xl+xu)/2.d0
            fxmid= sum(x_ymat(:,i1)*exp(bz+xmid))-term1
            if(fxmid.lt.0d0) then
              xl= xmid
            else
              xu= xmid
            end if

          else
 !           print*, 'try different extremes 2'
            exit
          end if
          err2= abs(fxmid)

        end do
        hnew(i1)= xmid
      end do
  
      hall= -50.d0

      i5= 0
      do i4= 1, n
        if(delta(i4).eq.1) then
          i5= i5+1
          delone(i5)= i4
        end if
      end do


      i3= 0
      do i2= minval(delone), n
        if(delta(i2).eq.1) then
          i3= i3+1
          hall(i2)= hnew(i3)
      
        else
          hall(i2)= hnew(i3)
        end if
      end do

!===============   H ESTIMATION ENDS HERE  =====================


      CALL betaestr0(beta_old2, delta, hall, k, n, rfix, z)
    
      theta_difff= sum(abs((beta_old2-theta_old)/theta_old))
      
    end do   !  end of while loop for theta_difff
    
    allbeta(:,i6)= beta_old2
    allhall(:,i6)= hall
    allh(:,i6)= hnew

    CALL newhallr0(hall, ltvec, n, outvec, tvec, x_ik_sort)

    newhsum= newhsum + outvec
    beta_old_sum= beta_old_sum + beta_old2

  end do    ! end of for loop i6= 1, m

 
  beta_old_mean= beta_old_sum/m
  newhavg= newhsum/m

!============  BETA AND H ESTIMATION ENDS HERE  =====================

!       print*, "betah_est is fine"

        beta_old= beta_old_mean

        CALL my_match(left_tvec_match, n, ltvec, left, tvec)
!!      output: left_tvec_match

        hallleft= newhavg(left_tvec_match)
        hallold= newhavg(orderforhold)
 
        beta_difff= sum(abs((beta_old-beta_older)/beta_older))

      end do   ! end of while loop for beta convergence

      
!     print*, 'final_beta_theta= ', beta_old
!     beta_neww= beta_neww+beta_old

     
        
!!========================================
!!            Sigma Computation
!!========================================

     kminus1= k-1
     kminus2= k-2

     do kk= 1, m
        ymat= allx_ymat(:,:,kk)
        myfinalbeta= allbeta(:,kk)
        myfinalh= allh(:,kk)
        myfinalhall= allhall(:,kk)
        mydel= alldel(:,kk)
        z= allz(:,:,kk)
            
        hdiff= myfinalh(2:k)-myfinalh(1:(k-1))
        mybzvec= z(:,1)*myfinalbeta(1)+z(:,2)*myfinalbeta(2)

        outh1= exp(mybzvec+myfinalhall)
        
        
        sig1=0.d0
        sig2=0.d0

        firstel1= 0.d0
        secondel1= 0.d0
        fourthel1= 0.d0
        
        firstel2= 0.d0
        secondel2= 0.d0
        thirdel2= 0.d0
        fourthel2= 0.d0

        cumsum_hdiff= 0.d0
        dhdalpha_numvec1= 0.d0
        dhdalpha_numvec2= 0.d0
        dhdalpha_numvec3= 0.d0
        
        do t1= 2, k
          t1minus1= t1-1
          bfnvec= exp(myfinalh(t1)-myfinalhall)
          mybzplush= mybzvec+myfinalh(t1)
         
          outh2= exp(mybzplush)
          hazy= outh2*ymat(:,t1)
          hazyvec= sum(hazy)

          dhazy= hazy
          dhazy_term1= myfinalbeta(1)*sum(dhazy*dxdalpha(:,1)) 
          dhazy_term2= myfinalbeta(1)*sum(dhazy*dxdalpha(:,2)) 
          dhazy_term3= myfinalbeta(1)*sum(dhazy*dxdalpha(:,3)) 
     
          oyb= outh1*ymat(:,t1)*bfnvec
          zbarvec(1)= sum(z(:,1)*oyb)/sum(hazy)
          zbarvec(2)= sum(z(:,2)*oyb)/sum(hazy)

          zminuszbar(:,1)= z(:,1)-zbarvec(1)
          zminuszbar(:,2)= z(:,2)-zbarvec(2)

          hazyhdiff= hazy*hdiff(t1-1)
          dhazyhdiff= dhazy*hdiff(t1-1)

          firstel1= firstel1+zminuszbar(:,1)*zminuszbar(:,1)*hazyhdiff
          secondel1= secondel1+zminuszbar(:,1)*zminuszbar(:,2)*hazyhdiff
          fourthel1= fourthel1+zminuszbar(:,2)*zminuszbar(:,2)*hazyhdiff

          firstel2= firstel2+zminuszbar(:,1)*z(:,1)*dhazyhdiff
          secondel2= secondel2+zminuszbar(:,2)*z(:,1)*dhazyhdiff
          thirdel2= thirdel2+zminuszbar(:,1)*z(:,2)*dhazyhdiff
          fourthel2= fourthel2+zminuszbar(:,2)*z(:,2)*dhazyhdiff


          cumsum_hdiff= cumsum_hdiff+hdiff(t1minus1)
          dhdalpha_den= exp(cumsum_hdiff)

          dhdalpha_numvec1= dhdalpha_numvec1+ exp(cumsum_hdiff)*&
                            hdiff(t1minus1)*dhazy_term1/hazyvec
          dhdalpha_numvec2= dhdalpha_numvec2+ exp(cumsum_hdiff)*&
                            hdiff(t1minus1)*dhazy_term2/hazyvec
          dhdalpha_numvec3= dhdalpha_numvec3+ exp(cumsum_hdiff)*&
                            hdiff(t1minus1)*dhazy_term3/hazyvec

          dhdalpha1(t1minus1)= -dhdalpha_numvec1/dhdalpha_den
          dhdalpha2(t1minus1)= -dhdalpha_numvec2/dhdalpha_den
          dhdalpha3(t1minus1)= -dhdalpha_numvec3/dhdalpha_den

        end do  ! end of for loop t1= 2, k

        sig1= (/sum(firstel1), sum(secondel1), sum(secondel1), sum(fourthel1)/)
        sig2= (/sum(firstel2), sum(secondel2), sum(thirdel2), sum(fourthel2)/)

        sigustar_old= sig1/n
        siglstar= sig2/n

        dhdalpha_all1=0.d0
        dhdalpha_all2=0.d0
        dhdalpha_all3=0.d0

        cc1= 1 
        do t3= delone(2), n
          if (mydel(t3) .eq. 1) then
            cc1= cc1+1
          end if

          dhdalpha_all1(t3)= dhdalpha1(cc1-1)
          dhdalpha_all2(t3)= dhdalpha2(cc1-1)
          dhdalpha_all3(t3)= dhdalpha3(cc1-1)
        end do
                        
        delpart(:,1)= dhdalpha_all1 + myfinalbeta(1)*dxdalpha(:,1)
        delpart(:,2)= dhdalpha_all2 + myfinalbeta(1)*dxdalpha(:,2)
        delpart(:,3)= dhdalpha_all3 + myfinalbeta(1)*dxdalpha(:,3)

        psi_siga_psi= 0.d0
        do t4= 1, n
          psi(1,:)= outh1(t4)*zstar(1,t4)*delpart(t4,:)
          psi(2,:)= outh1(t4)*zstar(2,t4)*delpart(t4,:)
          psi_tr(:,1)= psi(1,:)
          psi_tr(:,2)= psi(2,:)
          psi_siga_psi= psi_siga_psi + MATMUL(MATMUL(psi,sigma_alpha),psi_tr)
        end do
        
        sigustar(1)= sigustar_old(1)+ psi_siga_psi(1,1)/n
        sigustar(2)= sigustar_old(2)+ psi_siga_psi(2,1)/n
        sigustar(3)= sigustar_old(3)+ psi_siga_psi(1,2)/n
        sigustar(4)= sigustar_old(4)+ psi_siga_psi(2,2)/n

        det_siglstar= 1.d0/(siglstar(1)*siglstar(4)-siglstar(2)*siglstar(3))
        adj_siglstar= (/siglstar(4), -siglstar(2), -siglstar(3), siglstar(1)/)
        inv_siglstar= adj_siglstar*det_siglstar
        
        varbeta(1)= inv_siglstar(1)*inv_siglstar(1)*sigustar(1) +&
                    inv_siglstar(1)*inv_siglstar(3)*(sigustar(2)+sigustar(3)) +&
                    inv_siglstar(3)*inv_siglstar(3)*sigustar(4)
        varbeta(2)= inv_siglstar(2)*inv_siglstar(2)*sigustar(1) +&
                    inv_siglstar(2)*inv_siglstar(4)*(sigustar(2)+sigustar(3)) +&
                    inv_siglstar(4)*inv_siglstar(4)*sigustar(4)

        avhid= kk
        allvar_hope(1,avhid)= varbeta(1)/n
        allvar_hope(2,avhid)= varbeta(2)/n

        allbeta_hope(:,avhid)= allbeta(:,kk)
        
      end do      ! end of do loop for kk= 1, m
      

    beta_neww(1)= sum(allbeta_hope(1,:))/m
    beta_neww(2)= sum(allbeta_hope(2,:))/m


    ahvar_mean(1)= sum(allvar_hope(1,:))/m
    ahvar_mean(2)= sum(allvar_hope(2,:))/m

    beta_hope= allbeta_hope - spread(beta_neww,2,m)
    beta_hope_beta_hope_t(1)= sum(beta_hope(1,:)*beta_hope(1,:))
    beta_hope_beta_hope_t(2)= sum(beta_hope(2,:)*beta_hope(2,:))

    sig_hope= ahvar_mean + (m+1)*beta_hope_beta_hope_t/(m*(m-1))
    mysd= SQRT(sig_hope)


!print*, 'final ans= ', (/beta_neww,mysd/)
!print*, (/beta_neww,mysd/)

return
end subroutine


!==========================================================================
!                 REGRESSION CALIBRATION for general r case
!==========================================================================


subroutine rc_genr(beta_neww, beta_old, cumsumgrid,&
                       delta_temp, final_x, hallleft,&
                       hallold, isd, k, left,&
                       lrcpos, ltvec,&
                       m, maxgrid, mysd, n,&
                       notrcpos, numgrid,&
                       orderforhold, rcpos, rfix,&
                       sigma_alpha, totngrid,&
                       ttilmat, tvec, wbar, z2)
implicit none


integer(4):: avhid, ii, j, j2, k, k1, kk, kminus1, kminus2,&
            lrcpos, ltvec, m, maxgrid,&
            n, t1, t2, t3, t4, totngrid, tt,&
            sm, dellater, dellater1, dellater2
            
integer(4):: bfnindx1, bfnindx2, cumsumgrid(k), delta_temp(n),&
            left_tvec_match(n), notrcpos(k), numgrid(k),&
            orderforhold(totngrid), rcpos(lrcpos)

real(8):: beta_difff, beta_neww(2), beta_old(2), final_x(n),&
         left(n), mysd(2), rfix, ttilmat(totngrid),&
         tvec(ltvec), wbar(n), z2(n)

real(8):: beta_older(2), beta_old_mean(2), hallleft(n),&
         hallold(totngrid), init_beta(2), newhavg(ltvec),&
         newz1(n), x_ik(n,m), z(n,2)


integer(4):: cc, cc1, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10,&
            allind(k,m), delone(k), delta(n), inddel(k), ind_del(k),&
            t1minus1, this, x_ik_sort_indx(n)

real(8):: allbeta(2,m), alldel(n,m), allh(k,m), allhall(n,m), allx_ymat(n,k,m),&
         allz(n,2,m), beta_old_sum(2), beta_old2(2), bz(n), hall(n), hnew(k),&
         init_beta2(2), newhsum(ltvec), outvec(ltvec), theta_difff, theta_old(2),&
         x_ik_sort(n), x_t1t2(k), x_ymat(n,k), z_temp2(n,2)

real(8):: adj_siglstar(4), ahvar_mean(2), allbeta_hope(2,m),&
         allvar_hope(2,m), beta_diff(2,m), beta_hope(2,m),&
         beta_hope_beta_hope_t(2), bfnvec(n), det_siglstar,&
         dhazy(n), dhazyhdiff(n), firstel1(n), firstel2(n), fourthel1(n),&
         fourthel2(n), hazy(n), hazyhdiff(n), hdiff(k-1), integrand

real(8):: mn, mybzplushmat(n,k), mybzvec(n), myfinalbeta(2), myfinalh(k),&
         myfinalhall(n), outdh(n), outh(n), outh1(n), outh2(n), outh3(n),&
         oyb(n), ratio, secondel1(n), secondel2(n), sig1(4),&
         sig2(4), sig_hope(2), sum_int, theta_diff, thirdel2(n),&
         varbeta(2), ymat(n,k), zbarvec(2), zminuszbar(n,2),&
         zzprime1(4), zzprime2(4)

real(8):: siglstar(4), sigustar(4), sigustar_old(4), inv_siglstar(4)

real(8):: err1, err2, fxl1, fxl, fxmid1, fxmid, fxu1,&
         fxu, term1, xl, xmid, xu 

real(8):: cumsum_hdiff, delpart(n,3), dhazy_term1,&
         dhazy_term2, dhazy_term3, dhdalpha1(k-1),&
         dhdalpha2(k-1), dhdalpha3(k-1), dhdalpha_all1(n),&
         dhdalpha_all2(n), dhdalpha_all3(n), dhdalpha_den,&
         dhdalpha_numvec1, dhdalpha_numvec2, dhdalpha_numvec3,&
         dxdalpha(n,3), hazyvec, mydel(n), psi(2,3), psi_siga_psi(2,2),&
         psi_tr(3,2), sigma_alpha(3,3), zstar(2,n)


integer(4) :: is1, is2, isd, seed1, seed2        

    common /unif_seeds/ is1, is2
    save /unif_seeds/
    seed1=10005
    seed2=isd      
    CALL set_uniform(seed1, seed2)
  
  
    init_beta= beta_old
    beta_neww= 0.d0
    allbeta_hope= 0.d0
    allvar_hope= 0.d0

      beta_old= init_beta
      newz1= final_x
      beta_difff= 1.0d0
      
      dxdalpha(:,1)= 1.d0
      dxdalpha(:,2)= wbar
      dxdalpha(:,3)= z2
      
      zstar(1,:)= z2
      zstar(2,:)= final_x


      do while(beta_difff > 0.01d0)
      
        beta_older= beta_old
!       print*, "ok till here"

        CALL imp_t(beta_old, cumsumgrid, hallold, hallleft, k, left,&
                   lrcpos, m, maxgrid, n, newz1, notrcpos, numgrid,&
                   x_ik, rcpos, rfix, totngrid, ttilmat, z2)

!!      output: x_ik


!=============  BETA AND H ESTIMATION STEP ======================


        z_temp2(:,1)= newz1       
        z_temp2(:,2)= z2
  
        beta_old_sum= 0.d0
        newhsum= 0.d0
        init_beta2= beta_old


        do i6= 1, m


          beta_old2= init_beta2
          hall= 0.d0
  
          x_ik_sort= x_ik(:,i6)
          CALL bubble(x_ik_sort, n)           ! equivalent of sort()$x in R
          CALL my_match(x_ik_sort_indx, n, n, x_ik_sort, x_ik(:,i6))
                                              ! equivalent of sort()$ix in R
    
          z= z_temp2(x_ik_sort_indx,:)
          delta= delta_temp(x_ik_sort_indx)

          allz(:,:,i6)= z
          alldel(:,i6)= delta


          cc=1
          do i8= 1, n
            if(delta(i8)==1.d0) then
              x_t1t2(cc)= x_ik_sort(i8)
              inddel(cc)= i8
              cc= cc+1
            end if
          end do
    
          x_ymat=0.d0    
          do i7= 1, k
            where(x_ik_sort>= x_t1t2(i7)) x_ymat(:,i7)= 1.d0
          end do

          allind(:, i6)= inddel
          allx_ymat(:,:,i6)= x_ymat

!      print*, "step 3 done"

      theta_difff= 1.d0
      do while(theta_difff> 0.01d0)
        theta_old= beta_old2

!       CALL hest(beta_old2, delta, hall, hnew, k, n, rfix, x_ymat, z)

!===============  H ESTIMATION STEP  ==========================


      xl= -20.d0
      xu= 20.d0
      xmid= 0.d0
      bz= z(:, 1)*beta_old2(1)+z(:, 2)*beta_old2(2)
  
  
      err1= 1
      do while(err1> 0.01d0)
 
        fxl1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xl))/rfix)-1
        fxu1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xu))/rfix)-1

        if (fxl1*fxu1.lt.0.d0) then
          xmid= (xl+xu)/2.d0
          fxmid1= sum(x_ymat(:,1)*log(1+rfix*exp(bz+xmid))/rfix)-1
          if(fxmid1.lt.0.d0) then
            xl= xmid
          else
            xu= xmid
          end if

        else
!          print*, 'try different extremes 1'
          exit

        end if
        err1= abs(fxmid1)

      end do
    
      hnew(1)= xmid
  
      do i1= 2, k
        xl= -15.d0
        xu= 15.d0
        xmid= 0.d0

        err2= 1.d0
        term1=1+sum(x_ymat(:,i1)*log(1+rfix*exp(bz+hnew(i1-1)))/rfix)

        do while(err2> 0.01d0)

          fxl= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xl))/rfix)-term1
          fxu= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xu))/rfix)-term1

          if (fxl*fxu.lt.0d0) then
            xmid= (xl+xu)/2.d0
            fxmid= sum(x_ymat(:,i1)*log(1+rfix*exp(bz+xmid))/rfix)-term1
            if(fxmid.lt.0d0) then
              xl= xmid
            else
              xu= xmid
            end if

          else
 !           print*, 'try different extremes 2'
            exit
          end if
          err2= abs(fxmid)

        end do
        hnew(i1)= xmid
      end do
  
      hall= -50.d0

      i5= 0
      do i4= 1, n
        if(delta(i4).eq.1) then
          i5= i5+1
          delone(i5)= i4
        end if
      end do


      i3= 0
      do i2= minval(delone), n
        if(delta(i2).eq.1) then
          i3= i3+1
          hall(i2)= hnew(i3)
      
        else
          hall(i2)= hnew(i3)
        end if
      end do

!===============   H ESTIMATION ENDS HERE  =====================


      CALL betaest(beta_old2, delta, hall, k, n, rfix, z)
    
      theta_difff= sum(abs((beta_old2-theta_old)/theta_old))

    end do   !  end of while loop for theta_difff
    
    allbeta(:,i6)= beta_old2
    allhall(:,i6)= hall
    allh(:,i6)= hnew

    CALL newhall(hall, ltvec, n, outvec, tvec, x_ik_sort)

    newhsum= newhsum + outvec
    beta_old_sum= beta_old_sum + beta_old2

  end do    ! end of for loop i6= 1, m

 
  beta_old_mean= beta_old_sum/m
  newhavg= newhsum/m

!============  BETA AND H ESTIMATION ENDS HERE  =====================

!       print*, "betah_est is fine"

        beta_old= beta_old_mean

        CALL my_match(left_tvec_match, n, ltvec, left, tvec)
!!      output: left_tvec_match

        hallleft= newhavg(left_tvec_match)
        hallold= newhavg(orderforhold)

 
        beta_difff= sum(abs((beta_old-beta_older)/beta_older))

      end do   ! end of while loop for beta convergence

      
!     print*, 'final_beta_theta= ', beta_old
!     beta_neww= beta_neww+beta_old

     
      
      
      
!!========================================
!!            Sigma Computation
!!========================================

     kminus1= k-1
     kminus2= k-2

     do kk= 1, m
        ymat= allx_ymat(:,:,kk)
        myfinalbeta= allbeta(:,kk)
        myfinalh= allh(:,kk)
        myfinalhall= allhall(:,kk)
        ind_del= allind(:,kk)
        mydel= alldel(:,kk)
        z= allz(:,:,kk)
            
        hdiff= myfinalh(2:k)-myfinalh(1:(k-1))
        mybzvec= z(:,1)*myfinalbeta(1)+z(:,2)*myfinalbeta(2)

        outh1= exp(mybzvec+myfinalhall)/(1+rfix*exp(mybzvec+myfinalhall))

        do k1= 1, k
          mybzplushmat(:,k1)= exp(mybzvec+myfinalh(k1))
        end do


        sig1=0.d0
        sig2=0.d0

        firstel1= 0.d0
        secondel1= 0.d0
        fourthel1= 0.d0
        
        firstel2= 0.d0
        secondel2= 0.d0
        thirdel2= 0.d0
        fourthel2= 0.d0

        cumsum_hdiff= 0.d0
        dhdalpha_numvec1= 0.d0
        dhdalpha_numvec2= 0.d0
        dhdalpha_numvec3= 0.d0
 
        
        do t1= 2, k

          t1minus1= t1-1

          sum_int= 0.d0
          bfnvec= 1.d0

          do t2= t1, kminus2
          
            outh2= mybzplushmat(:,(t2+1))/(1+rfix*mybzplushmat(:,(t2+1)))
            hazy= outh2*ymat(:,(t2+1))
            outh3= mybzplushmat(:,(t2+1))/((1+rfix*&
                   mybzplushmat(:,(t2+1)))*(1+rfix*mybzplushmat(:,(t2+1))))
            dhazy= outh3*ymat(:,(t2+1))

            integrand= sum(dhazy)*hdiff(t2)/sum(hazy)
            sum_int= sum_int + integrand

            bfnindx1= ind_del(t2+1)
            bfnindx2= ind_del(t2+2)-1
            bfnvec(bfnindx1:bfnindx2)= exp(-sum_int)

          end do

          bfnindx1= ind_del(k)

          outh2= mybzplushmat(:,k)/(1+rfix*mybzplushmat(:,k))
          hazy= outh2*ymat(:,k)
          outh3= mybzplushmat(:,k)/((1+rfix*mybzplushmat(:,k))*&
                 (1+rfix*mybzplushmat(:,k)))
          dhazy= outh3*ymat(:,k)

          integrand= sum(dhazy)*hdiff(kminus1)/sum(hazy)
          sum_int= sum_int + integrand

          bfnvec(bfnindx1:n)= exp(-sum_int)

          if(t1== kminus1) then
            bfnindx1= ind_del(kminus1)
            bfnindx2= ind_del(k)-1
            bfnvec(bfnindx1:bfnindx2)= 1.d0
 
             
            outh2= mybzplushmat(:,k)/(1+rfix*mybzplushmat(:,k))
            hazy= outh2*ymat(:,k)
            outh3= mybzplushmat(:,k)/((1+rfix*mybzplushmat(:,k))*&
                   (1+rfix*mybzplushmat(:,k)))
            dhazy= outh3*ymat(:,k)
            
            integrand= sum(dhazy)*hdiff(kminus1)/sum(hazy)
            sum_int= sum_int + integrand

            bfnvec(ind_del(k):n)= exp(-sum_int)
          
          else if(t1==k) then

            bfnvec(ind_del(k):n)= 1.d0

          end if

         
          outh2= mybzplushmat(:,t1)/(1+rfix*mybzplushmat(:,t1))
          hazy= outh2*ymat(:,t1)
          hazyvec= sum(hazy)

          outh3= mybzplushmat(:,t1)/((1+rfix*mybzplushmat(:,t1))*&
                 (1+rfix*mybzplushmat(:,t1)))
          dhazy= outh3*ymat(:,t1)

          dhazy_term1= myfinalbeta(1)*sum(dhazy*dxdalpha(:,1)) 
          dhazy_term2= myfinalbeta(1)*sum(dhazy*dxdalpha(:,2)) 
          dhazy_term3= myfinalbeta(1)*sum(dhazy*dxdalpha(:,3)) 
     
          oyb= outh1*ymat(:,t1)*bfnvec
          zbarvec(1)= sum(z(:,1)*oyb)/hazyvec
          zbarvec(2)= sum(z(:,2)*oyb)/hazyvec

          zminuszbar(:,1)= z(:,1)-zbarvec(1)
          zminuszbar(:,2)= z(:,2)-zbarvec(2)

          hazyhdiff= hazy*hdiff(t1-1)
          dhazyhdiff= dhazy*hdiff(t1-1)

          firstel1= firstel1+zminuszbar(:,1)*zminuszbar(:,1)*hazyhdiff
          secondel1= secondel1+zminuszbar(:,1)*zminuszbar(:,2)*hazyhdiff
          fourthel1= fourthel1+zminuszbar(:,2)*zminuszbar(:,2)*hazyhdiff

          firstel2= firstel2+zminuszbar(:,1)*z(:,1)*dhazyhdiff
          secondel2= secondel2+zminuszbar(:,2)*z(:,1)*dhazyhdiff
          thirdel2= thirdel2+zminuszbar(:,1)*z(:,2)*dhazyhdiff
          fourthel2= fourthel2+zminuszbar(:,2)*z(:,2)*dhazyhdiff


          cumsum_hdiff= cumsum_hdiff+sum(dhazy)*hdiff(t1minus1)/hazyvec
          dhdalpha_den= exp(cumsum_hdiff)

          dhdalpha_numvec1= dhdalpha_numvec1+ exp(cumsum_hdiff)*&
                            hdiff(t1minus1)*dhazy_term1/hazyvec
          dhdalpha_numvec2= dhdalpha_numvec2+ exp(cumsum_hdiff)*&
                            hdiff(t1minus1)*dhazy_term2/hazyvec
          dhdalpha_numvec3= dhdalpha_numvec3+ exp(cumsum_hdiff)*&
                            hdiff(t1minus1)*dhazy_term3/hazyvec

          dhdalpha1(t1minus1)= -dhdalpha_numvec1/dhdalpha_den
          dhdalpha2(t1minus1)= -dhdalpha_numvec2/dhdalpha_den
          dhdalpha3(t1minus1)= -dhdalpha_numvec3/dhdalpha_den

        end do  ! end of for loop t1= 2, k


        sig1= (/sum(firstel1), sum(secondel1), sum(secondel1), sum(fourthel1)/)
        sig2= (/sum(firstel2), sum(secondel2), sum(thirdel2), sum(fourthel2)/)

        sigustar_old= sig1/n
        siglstar= sig2/n



        cc1= 1
        dhdalpha_all1=0.d0
        dhdalpha_all2=0.d0
        dhdalpha_all3=0.d0

        do t3= delone(2), n
          if (mydel(t3) .eq. 1) then
            cc1= cc1+1
          end if

          dhdalpha_all1(t3)= dhdalpha1(cc1-1)
          dhdalpha_all2(t3)= dhdalpha2(cc1-1)
          dhdalpha_all3(t3)= dhdalpha3(cc1-1)
        end do
                        
        delpart(:,1)= dhdalpha_all1 + myfinalbeta(1)*dxdalpha(:,1)
        delpart(:,2)= dhdalpha_all2 + myfinalbeta(1)*dxdalpha(:,2)
        delpart(:,3)= dhdalpha_all3 + myfinalbeta(1)*dxdalpha(:,3)

        psi_siga_psi= 0.d0
        do t4= 1, n
          psi(1,:)= outh1(t4)*zstar(1,t4)*delpart(t4,:)
          psi(2,:)= outh1(t4)*zstar(2,t4)*delpart(t4,:)
          psi_tr(:,1)= psi(1,:)
          psi_tr(:,2)= psi(2,:)
          psi_siga_psi= psi_siga_psi + MATMUL(MATMUL(psi,sigma_alpha),psi_tr)
        end do

       
        sigustar(1)= sigustar_old(1)+ psi_siga_psi(1,1)/n
        sigustar(2)= sigustar_old(2)+ psi_siga_psi(2,1)/n
        sigustar(3)= sigustar_old(3)+ psi_siga_psi(1,2)/n
        sigustar(4)= sigustar_old(4)+ psi_siga_psi(2,2)/n

        det_siglstar= 1.d0/(siglstar(1)*siglstar(4)-siglstar(2)*siglstar(3))
        adj_siglstar= (/siglstar(4), -siglstar(2), -siglstar(3), siglstar(1)/)
        inv_siglstar= adj_siglstar*det_siglstar
        
        varbeta(1)= inv_siglstar(1)*inv_siglstar(1)*sigustar(1) +&
                    inv_siglstar(1)*inv_siglstar(3)*(sigustar(2)+sigustar(3)) +&
                    inv_siglstar(3)*inv_siglstar(3)*sigustar(4)
        varbeta(2)= inv_siglstar(2)*inv_siglstar(2)*sigustar(1) +&
                    inv_siglstar(2)*inv_siglstar(4)*(sigustar(2)+sigustar(3)) +&
                    inv_siglstar(4)*inv_siglstar(4)*sigustar(4)

        avhid= kk
        allvar_hope(1,avhid)= varbeta(1)/n
        allvar_hope(2,avhid)= varbeta(2)/n

        allbeta_hope(:,avhid)= allbeta(:,kk)
        
      end do      ! end of do loop for kk= 1, m
      
    beta_neww(1)= sum(allbeta_hope(1,:))/m
    beta_neww(2)= sum(allbeta_hope(2,:))/m

    ahvar_mean(1)= sum(allvar_hope(1,:))/m
    ahvar_mean(2)= sum(allvar_hope(2,:))/m

    beta_hope= allbeta_hope - spread(beta_neww,2,m)
    beta_hope_beta_hope_t(1)= sum(beta_hope(1,:)*beta_hope(1,:))
    beta_hope_beta_hope_t(2)= sum(beta_hope(2,:)*beta_hope(2,:))

    sig_hope= ahvar_mean + (m+1)*beta_hope_beta_hope_t/(m*(m-1))
    mysd= SQRT(sig_hope)


!print*, 'final ans= ', (/beta_neww,mysd/)
!print*, (/beta_neww,mysd/)

return
end subroutine



 

!===============================================================
!                 RANDOM SAMPLING CODES
!===============================================================
! Set seeds for the uniform random number generator.
!
!       Generate uniformly distributed random numbers using the 32-bit
!       generator from figure 3 of:
!       L'Ecuyer, P. Efficient and portable combined random number
!       generators, C.A.C.M., vol. 31, 742-749 & 774-?, June 1988.
!
!       The cycle length is claimed to be 2.30584E+18
!
!       Seeds can be set by calling the routine set_uniform
!
!       It is assumed that the Fortran compiler supports long variable
!       names, and integer(4).
!


        SUBROUTINE set_uniform(iseed1, iseed2)
        IMPLICIT NONE
        INTEGER(4) is1, is2, iseed1, iseed2
        common /unif_seeds/ is1, is2
        save /unif_seeds/

        is1 = iseed1
        is2 = iseed2
        return
        end


!================================================
! Normal random number generator
! Here v is variance
 
        DOUBLE PRECISION FUNCTION ran_nor(m,v)
        IMPLICIT NONE
        DOUBLE PRECISION m, v, uniform, pi,z,huge
        INTEGER(4) is1, is2
        common /unif_seeds/ is1, is2
        save /unif_seeds/
 
        pi = 3.1415926536d0
       huge = 1e35
 1     z=dsqrt(-2.d0*dlog(uniform()))*dcos(2.d0*pi*uniform())
       if (dabs(z).ge.huge) goto 1
       ran_nor = m + sqrt(v)*z
       return
       end

!================================================
! Uniform random number generator

        double precision function uniform()
        IMPLICIT NONE

        integer(4) z, k, is1, is2
        common /unif_seeds/ is1, is2
        save /unif_seeds/
!
        k = is1 / 53668.d0
        is1 = 40014.d0 * (is1 - k * 53668.d0) - k * 12211.d0
        if (is1 .lt. 0.d0) then
           is1 = is1 + 2147483563.d0
        endif
!
        k = is2 / 52774.d0
        is2 = 40692.d0 * (is2 - k * 52774.d0) - k * 3791.d0
        if (is2 .lt. 0.d0) then
           is2 = is2 + 2147483399.d0
        endif
!
        z = is1 - is2
        if (z .lt. 1.d0) then
           z = z + 2147483562.d0
        endif
!
        uniform = z / 2147483563.d0
        return
        end
        
!==================================================================
! Gamma random number generator
! Here (aa, bb) is in the (k, theta) form mentioned in wiki

          double precision function gam_ran(aa, bb)
          IMPLICIT NONE
          double precision aa, bb, d, c, w, u, v, x, z, yy, e, uniform
          logical acc
          integer(4) is1, is2
          common /unif_seeds/ is1, is2
          save /unif_seeds/
!         integer(4) l
          

          if (aa.lt.1.d0) then
!               { ********************* Johnk's algorithm***}
         x = 2.d0
         yy=0.d0
         do while((x+yy).gt.1.d0)
           u = uniform()
           v = uniform()
           x = exp(log(u)/aa)
           yy = exp(log(v)/(1.d0-aa))
           enddo
         e = -log(uniform())
         z = e*x/(x + yy)
         gam_ran = z * bb

        else if (aa.eq.1.d0) then
        gam_ran = -bb * log(uniform())
        else
!               { ************************* Best's algorithm***}
        d = aa-1.d0
        c = 3.d0*aa - (.75d0)
        acc = .false.
        do while (.not.acc)
          u = uniform()
          v = uniform()
          w = u*(1.d0-u)
          yy = sqrt(c/w)*(u-0.5d0)
          x = d+yy
          if (x.gt.0.d0) then
            z = 64.d0*w**3*v**2
            acc = (z.le.(1.d0-2.d0*yy**2/x))
            if (.not.acc) acc = (log(z).le.(2.d0 * (d*log(x/d)-yy) ))
            endif
          enddo
        gam_ran = x * bb
       endif
       return
       end
       
!======================================
! Multinomial sample
! check: 0<prob<1 and sum(prob)=1
! input: capn is the length of the probability vector
!        prob is the probability vector
! output: integer in the range [1, capn]


INTEGER FUNCTION multinom(capn, prob)
implicit none

integer(4):: capn, is1, is2, kk
real(8):: prob(capn), r_unif, sum, uniform

  common /unif_seeds/ is1, is2
  save /unif_seeds/

  r_unif= uniform()
  kk= 0
  sum= 0.d0
  
  do while(sum<r_unif)
    kk= kk+1
    sum= sum+prob(kk)
  end do
  
  multinom= kk

return
end


!======================================
! Dirichlet sample
! method mentioned in wiki article for Dirichlet (from gamma)
! input: vec of length cat
! output: out of the same length

subroutine ran_dir(cat, out, vec)
implicit none

integer(4):: cat, indx, is1, is2
real(8):: dir_den, gam_ran, gam_val(cat), out(cat), vec(cat)

  common /unif_seeds/ is1, is2
  save /unif_seeds/


  dir_den= 0.d0
  do indx= 1, cat
    gam_val(indx)= gam_ran(vec(indx), 1.d0)
    dir_den= dir_den+ gam_val(indx)
  end do
  
  out= gam_val/dir_den

return
end subroutine
!
!
!
! SIMULATION OF X AND MIXTURE PARAMETERS
!
!
subroutine mh_parkcmp(des_mat, isd, mm, n, ncolp, &
                      ncomp, ndim, capr, wbar, sumsqwdiff, &
                      a_sigmau, a_sigmax, b_sigmau, b_sigmax,& 
                      prmn, prvar, store_gamma, store_myx, store_pi,& 
                      store_sigma2u, store_sigma2x,  store_vec)
implicit none

integer(4):: capr, mm, n, ncolp, ncomp, ndim 
        
real(8)::  des_mat(n, ncolp), wbar(n), sumsqwdiff, a_sigmau, a_sigmax, & 
          b_sigmau, b_sigmax, prmn(ncomp, ncolp), prvar(ncomp, ncolp), & 
          store_gamma(capr, ndim), store_myx(1000, n), store_pi(capr, ncomp), & 
          store_sigma2u(capr), store_sigma2x(capr, ncomp), store_vec(capr, ncomp)  
          
real(8)::  indic(n, ncomp), mn, newgamma(ncolp), myx(n), nlk, olk,   & 
          scale_u, shape_u,  sigma2u, sigma2x(ncomp), ssdiag(n, ncolp), pi(ncomp), & 
          prob(ncomp), pout(ncomp), prop, sstemp1, tempo(n, ncomp), temp1(n), &
          temp100(ncomp), temp200(ncomp), vgamma(ndim), vec(ncomp), vr

integer(4) :: capr1, i,  is1, is2, isd, j, k, k0, k1,& 
             k2, myr, seed1, seed2  

integer(4):: multinom
real(8)::  gam_ran, ran_nor, uniform


  common /unif_seeds/ is1, is2
  save /unif_seeds/
  seed1=10005
  seed2=isd      
  CALL set_uniform(seed1, seed2)
  prob=1.d0/ncomp

   indic=0.d0
   do i=1, n
     indic(i, multinom(ncomp, prob))=1.d0
   end do
   do k=1, ncomp
     vec(k)=sum(indic(:, k))
   end do
 


   do j=1, ndim
     CALL RANDOM_NUMBER(mn)
     vgamma(j)=mn
   end do
   sigma2x=1.d0
   sigma2u=1.d0
!
   ssdiag(:, 1)=1.d0
   do j=2, ncolp
     ssdiag(:, j)= des_mat(:, j)*des_mat(:, j)    
   end do

    scale_u=1/(sumsqwdiff*0.5d0+1/b_sigmau)
    shape_u=0.5d0*n*mm + a_sigmau
!    
     if(capr>1000) then 
       capr1= capr-1000
     else 
       capr1=capr
     endif

! beginning of MCMC

     do myr=1, capr1

      do k=1, ncomp 
       k0=(k-1)*ncolp
       k2=k*ncolp        
       do j=1, ncolp
        newgamma=vgamma(k0+1:k2)
        newgamma(j)=0.d0
        vr=1/(sum(ssdiag(:, j)*indic(:, k))/(sigma2x(k)+sigma2u/mm)+1/prvar(k, j))
        mn=vr*(sum(indic(:, k)*(wbar-matmul(des_mat, newgamma))*des_mat(:,j))/&
        (sigma2x(k)+sigma2u/mm)+prmn(k, j)/prvar(k, j))
        vgamma(k0+j)=ran_nor(mn, vr)
       end do     
      end do
!
      store_gamma(myr, :)=vgamma
!
      do k=1, ncomp
       k0=(k-1)*ncolp
       k2=k*ncolp   
       tempo(:, k) =matmul(des_mat,vgamma(k0+1:k2))
!
       prop=1.d0/gam_ran(3.d0, 1.d0)
!      
!     
       temp1=indic(:, k)*(wbar- tempo(:, k))
       sstemp1=sum(temp1*temp1)
       olk=-0.5* sstemp1/(sigma2x(k)+sigma2u/mm)- & 
       0.5*vec(k)*log(sigma2x(k)+sigma2u/mm)
       nlk=-0.5* sstemp1/(prop+sigma2u/mm)-   &
       0.5*vec(k)*log(prop+sigma2u/mm)
       if(uniform()< exp(nlk-olk)) sigma2x(k)=prop
      end do
      store_sigma2x(myr, :)=sigma2x 
      indic=0.d0   
      do i=1, n 
        temp100=wbar(i)-tempo(i, :)
        temp200=0.5*temp100*temp100/(sigma2x+sigma2u/mm)
        prob= pi*exp(-temp200)/sqrt(sigma2x+sigma2u/mm)
        prob=prob/sum(prob)
        indic(i, multinom(ncomp, prob))=1.d0
      end do
       
      do k=1, ncomp
       vec(k)=sum(indic(:, k))
      end do 
      store_vec(myr, :)=vec
      call ran_dir(ncomp, pout, (vec+1.d0))
      pi=pout
      store_pi(myr, :)=pi
!
      sigma2u=1/gam_ran(shape_u, scale_u)
      store_sigma2u(myr)=sigma2u
!
     end do
! 
    if(capr>1000) then 

     do myr=(capr1+1), capr         
     
      do k=1, ncomp 
       k0=(k-1)*ncolp
       k2=k*ncolp        
       do j=1, ncolp
        newgamma=vgamma(k0+1:k2)
        newgamma(j)=0.d0
        vr=1/(sum(ssdiag(:, j)*indic(:, k))/(sigma2x(k)+sigma2u/mm)+1/prvar(k, j))
        mn=vr*(sum(indic(:, k)*(wbar-matmul(des_mat, newgamma))*des_mat(:,j))/&
        (sigma2x(k)+sigma2u/mm)+prmn(k, j)/prvar(k, j))
        vgamma(k0+j)=ran_nor(mn, vr)
       end do     
      end do
!
      store_gamma(myr, :)=vgamma
!
      do k=1, ncomp
       k0=(k-1)*ncolp
       k2=k*ncolp   
       tempo(:, k) =matmul(des_mat,vgamma(k0+1:k2))
!
       prop=1.d0/gam_ran(3.d0, 1.d0)
!      
!     
       temp1=indic(:, k)*(wbar- tempo(:, k))
       sstemp1=sum(temp1*temp1)
       olk=-0.5* sstemp1/(sigma2x(k)+sigma2u/mm)- & 
       0.5*vec(k)*log(sigma2x(k)+sigma2u/mm)
       nlk=-0.5* sstemp1/(prop+sigma2u/mm)-   &
       0.5*vec(k)*log(prop+sigma2u/mm)
       if(uniform()< exp(nlk-olk)) sigma2x(k)=prop
      end do
      store_sigma2x(myr, :)=sigma2x 
      indic=0.d0   
      do i=1, n 
        temp100=wbar(i)-tempo(i, :)
        temp200=0.5*temp100*temp100/(sigma2x+sigma2u/mm)
        prob= pi*exp(-temp200)/sqrt(sigma2x+sigma2u/mm)
        prob=prob/sum(prob)
        indic(i, multinom(ncomp, prob))=1.d0
      end do
       
      do k=1, ncomp
       vec(k)=sum(indic(:, k))
      end do 
       store_vec(myr, :)=vec
      call ran_dir(ncomp, pout, (vec+1.d0))
      pi=pout
!
      store_pi(myr, :)=pi
      sigma2u=1/gam_ran(shape_u, scale_u)
      store_sigma2u(myr)=sigma2u
 
     
      do i=1, n 
       vr= 1/ (sum(indic(i, :)/sigma2x) +mm/sigma2u)
       mn=vr*(wbar(i)*mm/sigma2u+ sum(indic(i, :)*tempo(i, :)/sigma2x))
       myx(i)=ran_nor(mn, vr)  
      enddo
      store_myx(myr-capr1, :)=myx
     end do
    endif
!
return
end subroutine
!
!
!
!
!
subroutine mh_parkcmpw1(des_mat, isd, mm, n, ncolp, &
                        ncomp, ndim, capr, wbar, &
                        a_sigmax, b_sigmax, prmn, prvar,&
                        store_gamma, store_myx, store_pi,& 
                        sigma2u, store_sigma2x, store_vec)
implicit none

integer(4):: capr, mm, n, ncolp, ncomp, ndim 
        
real(8)::  des_mat(n, ncolp), wbar(n), a_sigmax, & 
          b_sigmax, prmn(ncomp, ncolp), prvar(ncomp, ncolp), & 
          store_gamma(capr, ndim), store_myx(1000, n), store_pi(capr, ncomp), & 
          store_sigma2x(capr, ncomp), store_vec(capr, ncomp)  
          
real(8)::  indic(n, ncomp), mn, newgamma(ncolp), myx(n), nlk, olk,   & 
          sigma2u, sigma2x(ncomp), ssdiag(n, ncolp), pi(ncomp), & 
          prob(ncomp), pout(ncomp), prop, sstemp1, tempo(n, ncomp), temp1(n), &
          temp100(ncomp), temp200(ncomp), vgamma(ndim), vec(ncomp), vr

integer(4) :: capr1, i,  is1, is2, isd, j, k, k0, k1,& 
             k2, myr, seed1, seed2  

integer(4):: multinom
real(8)::  gam_ran, ran_nor, uniform


  common /unif_seeds/ is1, is2
  save /unif_seeds/
  seed1=10005
  seed2=isd      
  CALL set_uniform(seed1, seed2)
  prob=1.d0/ncomp

   indic=0.d0
   do i=1, n
     indic(i, multinom(ncomp, prob))=1.d0
   end do
   do k=1, ncomp
     vec(k)=sum(indic(:, k))
   end do
 


   do j=1, ndim
     CALL RANDOM_NUMBER(mn)
     vgamma(j)=mn
   end do
   sigma2x=1.d0
!   sigma2u=1.d0
!
   ssdiag(:, 1)=1.d0
   do j=2, ncolp
     ssdiag(:, j)= des_mat(:, j)*des_mat(:, j)    
   end do

!    scale_u=1/(sumsqwdiff*0.25d0+1/b_sigmau)
!    shape_u=0.5d0*n+a_sigmau
!    
     if(capr>1000) then 
       capr1= capr-1000
     else 
       capr1=capr
     endif

! begining of MCMC

     do myr=1, capr1

      do k=1, ncomp 
       k0=(k-1)*ncolp
       k2=k*ncolp        
       do j=1, ncolp
        newgamma=vgamma(k0+1:k2)
        newgamma(j)=0.d0
        vr=1/(sum(ssdiag(:, j)*indic(:, k))/(sigma2x(k)+sigma2u/mm)+1/prvar(k, j))
        mn=vr*(sum(indic(:, k)*(wbar-matmul(des_mat, newgamma))*des_mat(:,j))/&
        (sigma2x(k)+sigma2u/mm)+prmn(k, j)/prvar(k, j))
        vgamma(k0+j)=ran_nor(mn, vr)
       end do     
      end do
!
      store_gamma(myr, :)=vgamma
!
      do k=1, ncomp
       k0=(k-1)*ncolp
       k2=k*ncolp   
       tempo(:, k)=matmul(des_mat,vgamma(k0+1:k2))
!
       prop=1.d0/gam_ran(3.d0, 1.d0)
!      
!     
       temp1=indic(:, k)*(wbar- tempo(:, k))
       sstemp1=sum(temp1*temp1)
       olk=-0.5* sstemp1/(sigma2x(k)+sigma2u/mm)- & 
       0.5*vec(k)*log(sigma2x(k)+sigma2u/mm)
       nlk=-0.5* sstemp1/(prop+sigma2u/mm)-   &
       0.5*vec(k)*log(prop+sigma2u/mm)
       if(uniform()< exp(nlk-olk)) sigma2x(k)=prop
      end do
      store_sigma2x(myr, :)=sigma2x 
      indic=0.d0   
      do i=1, n 
        temp100=wbar(i)-tempo(i, :)
        temp200=0.5*temp100*temp100/(sigma2x+sigma2u/mm)
        prob= pi*exp(-temp200)/sqrt(sigma2x+sigma2u/mm)
        prob=prob/sum(prob)
        indic(i, multinom(ncomp, prob))=1.d0
      end do
       
      do k=1, ncomp
       vec(k)=sum(indic(:, k))
      end do 
      store_vec(myr, :)=vec
      call ran_dir(ncomp, pout, (vec+1.d0))
      pi=pout
      store_pi(myr, :)=pi
!
!      sigma2u=1/gam_ran(shape_u, scale_u)
!      store_sigma2u(myr)=sigma2u
!
     end do
! 
    if(capr>1000) then 

     do myr=(capr1+1), capr         
     
      do k=1, ncomp 
       k0=(k-1)*ncolp
       k2=k*ncolp        
       do j=1, ncolp
        newgamma=vgamma(k0+1:k2)
        newgamma(j)=0.d0
        vr=1/(sum(ssdiag(:, j)*indic(:, k))/(sigma2x(k)+sigma2u/mm)+1/prvar(k, j))
        mn=vr*(sum(indic(:, k)*(wbar-matmul(des_mat, newgamma))*des_mat(:,j))/&
        (sigma2x(k)+sigma2u/mm)+prmn(k, j)/prvar(k, j))
        vgamma(k0+j)=ran_nor(mn, vr)
       end do     
      end do
!
      store_gamma(myr, :)=vgamma
!
      do k=1, ncomp
       k0=(k-1)*ncolp
       k2=k*ncolp   
       tempo(:, k) =matmul(des_mat,vgamma(k0+1:k2))
!
       prop=1.d0/gam_ran(3.d0, 1.d0)
!      
!     
       temp1=indic(:, k)*(wbar- tempo(:, k))
       sstemp1=sum(temp1*temp1)
       olk=-0.5* sstemp1/(sigma2x(k)+sigma2u/mm)- & 
       0.5*vec(k)*log(sigma2x(k)+sigma2u/mm)
       nlk=-0.5* sstemp1/(prop+sigma2u/mm)-   &
       0.5*vec(k)*log(prop+sigma2u/mm)
       if(uniform()< exp(nlk-olk)) sigma2x(k)=prop
      end do
      store_sigma2x(myr, :)=sigma2x 
      indic=0.d0   
      do i=1, n 
        temp100=wbar(i)-tempo(i, :)
        temp200=0.5*temp100*temp100/(sigma2x+sigma2u/mm)
        prob= pi*exp(-temp200)/sqrt(sigma2x+sigma2u/mm)
        prob=prob/sum(prob)
        indic(i, multinom(ncomp, prob))=1.d0
      end do
       
      do k=1, ncomp
       vec(k)=sum(indic(:, k))
      end do 
       store_vec(myr, :)=vec
      call ran_dir(ncomp, pout, (vec+1.d0))
      pi=pout
!
      store_pi(myr, :)=pi
!      sigma2u=1/gam_ran(shape_u, scale_u)
!      store_sigma2u(myr)=sigma2u
 
     
      do i=1, n 
       vr= 1/ (sum(indic(i, :)/sigma2x) +mm/sigma2u)
       mn=vr*(wbar(i)*mm/sigma2u+ sum(indic(i, :)*tempo(i, :)/sigma2x))
       myx(i)=ran_nor(mn, vr)  
      enddo
      store_myx(myr-capr1, :)=myx
     end do
    endif
!
return
end subroutine
!

