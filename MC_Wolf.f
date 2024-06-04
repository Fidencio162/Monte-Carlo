! The diameter of the particles is used as unit length.
! Fidencio Pérez-Hernández 25.09.2022
 	implicit double precision(a-h,o-z)
	parameter(mp=1000000,mr=2**11,nmq=2**6,nvq=200)
	dimension x(mp),y(mp),z(mp)
	dimension qi(nmq),Sq(nmq),s(nmq)
	dimension qx(nmq,nvq),qy(nmq,nvq),qz(nmq,nvq)
	dimension r(mr),g11(mr),g12(mr),g22(mr)
	common/box/boxl,rc,np
	common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	common/parameters/diam,rho
	
	pi = 4.d0 * datan(1.d0)
	
	np = 8**3
	sigma1d = 1.53d0
	sigma2d = 1.d0
	sigma12d = 0.5d0 * (sigma1d + sigma2d)
	xi1 = 0.5d0
	xi2 = 1.d0 - xi1
	
	! packing fraction
	rho = 0.0285678d0
	d = (1.d0 / rho)**(1./3.)
	beta = 1.96d0
	temp = 1.d0 / beta
	
	! box length
	boxl = (dfloat(np) / rho)**(1.d0 / 3.d0)
	rc = boxl / 2.d0
	dr = rc / float(mr)
	nz1 = 1
	nz2 = -1
	
	print *, 'The length of the box is: ', boxl
	
	np1 = -nz2 * np / (nz1 - nz2)
	np2 = np - np1
	rho1 = np1 / boxl**3
	rho2 = np2 / boxl**3
	alpha = 4.d0 / boxl
	
	print *, 'rho', rho, 'rho1=', rho1, 'rho2=', rho2
	print *, 'n_c', np1, 'n_a', np2
	
	phi1 = pi * rho1 * sigma1d**3 / 6.d0
	phi2 = pi * rho2 * sigma2d**3 / 6.d0
	phi = np1 * pi * (sigma1d**3) / (6.d0 * boxl**3)
	phi = phi + pi * np2 * (sigma2d**3) / (6.d0 * boxl**3)
	
	print *, 'phi=', phi, 'phi1=', phi1, 'phi2=', phi2
	print *, 'rhoT=', rho, 'rho1=', rho1, 'rho2=', rho2
	
	call iniconf(x, y, z, d)
	
	open(10, file='iniconf.xyz', status='unknown')
	write(10, 49) np
	write(10, 50)
	do i = 1, np1
	    write(10, 51) x(i), y(i), z(i)
	enddo
	do i = np1 + 1, np
	    write(10, 52) x(i), y(i), z(i)
	enddo
	close(10)
	
	open(30, file='energyMC.dat', status='unknown')
	del = 0.1d0
	nattemp = 0
	nacc = 1
	iseed = 123456789
	
	do i = 1, 2000000
	    call mcmove(x, y, z, ener, nattemp, nacc, del, iseed)
	    call adjust(nattemp, nacc, del)
	    if (mod(i, 100) .eq. 0) write(30, *) i*1.d0, ener / np
	enddo
	
	print *, 'The system has thermalized :)'
	
	open(20, file='finalconf.xyz', status='unknown')
	do i = 1, np
	    write(20, *) x(i), y(i), z(i)
	enddo
	close(20)
	
	open(unit=13, file='gr_NaCl.dat')
	nacco = nacc
	ncp = 0
	do i = 1, mr
	    g11(i) = 0.d0
	    g12(i) = 0.d0
	    g22(i) = 0.d0
	enddo
	print *, 'Calculating g(r)...'
	do j = 1, 3000000
	    call mcmove(x, y, z, ener, nattemp, nacc, del, iseed)
	    call adjust(nattemp, nacc, del)
	    if (mod(j, np) .eq. 0) then
	        call average(x, y, z, g11, g12, g22, nattemp, nacc, del, dr, iseed)
	        ncp = ncp + 1
	    endif
	enddo
	
	print *, 'Average number:', ncp
	
	! Calculating g(r)
	do i = 2, mr
	    r(i) = (i-1) * dr
	    dv1 = (4.d0 * pi * r(i)**2 * dr) * rho1
	    dv12 = (4.d0 * pi * r(i)**2 * dr) * rho2
	    dv2 = (4.d0 * pi * r(i)**2 * dr) * rho2
	    g11(i) = g11(i) / (np1 * ncp * dv1)
	    g12(i) = g12(i) / (np1 * ncp * dv12)
	    g22(i) = g22(i) / (np2 * ncp * dv2)
	    write(13, *) r(i), g11(i), g12(i), g22(i)
	enddo
	close(13)
	
	print *, 'Calculating g(r) finished'
	print *, 'Program finished!'
	
	50 format(x, 'CONFIGURACION INICIAL')
	51 format(3X, 'C', 3f15.7)
	52 format(3X, 'N', 3f15.7)
	49 format(2x, i5)
	98 format(x, 'CONFIGURACION', 11x, i9)
	99 format(20x, i5)
	100 format(3f15.3)
	200 format(2f15.7)
	201 format(5x, '1', 2x, 'Sol', 4x, 'C', 2x, i5, 3f8.3, 2x, '0.0000', 2x, '0.0000', 2x, '0.0000')
	202 format(5x, '1', 2x, 'Sol', 4x, 'N', 2x, i5, 3f8.3, 2x, '0.0000', 2x, '0.0000', 2x, '0.0000')
	end	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c This subroutine calculates the initial configuration in 3D.		 
      	subroutine iniconf(xc,yc,zc,d)
      	implicit double precision(a-h,o-z)
	  parameter(mp=1000000,mr=2**11,nmq=2**6,nvq=200)
	  dimension xc(mp),yc(mp),zc(mp)
	  dimension xr(mp),yr(mp),zr(mp)
	  common/box/boxl,rc,np

	  xc(1)=-(boxl-d)/2.d0
	  yc(1)=-(boxl-d)/2.d0
	  zc(1)=-(boxl-d)/2.d0
     	do i=2,np
         xc(i)=xc(i-1)+d
         yc(i)=yc(i-1)
         zc(i)=zc(i-1)
         if (xc(i) .gt. boxl/2.) then
            xc(i)=xc(1)
            yc(i)=yc(i-1)+d
            if (yc(i) .gt. boxl/2.) then
               xc(i)=xc(1)
               yc(i)=yc(1)
               zc(i)=zc(i-1)+d
            endif
         endif
      	enddo

	  return
      	end

c This configuration calculates the energy of a given configuration
	subroutine energy(x,y,z,xj,yj,zj,ener,j)
	  implicit double precision(a-h,o-z)
	  parameter(mp=1000000)
	dimension x(mp),y(mp),z(mp)
      	  common/box/boxl,rc,np
	  common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	  common/parameters/diam,rho
	  pi=4.d0*datan(1.d0)
	  ener=0.d0
      	do i=1,j-1
            !uij=0.d0
            xij=x(i)-xj
            yij=y(i)-yj
	    zij=z(i)-zj
            xij=xij-boxl*dnint(xij/boxl) 
            yij=yij-boxl*dnint(yij/boxl) 
	    zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2) 
            if (rij .lt. rc) then
               call Wolf(rij,uij,i,j)
               ener=ener+uij
            endif 
      	enddo
      	
      	 do i=j+1,np
            xij=x(i)-xj
            yij=y(i)-yj
	    zij=z(i)-zj
            xij=xij-boxl*dnint(xij/boxl) 
            yij=yij-boxl*dnint(yij/boxl) 
	    zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2) 
            if (rij .lt. rc) then
               call Wolf(rij,uij,i,j)
               ener=ener+uij
            endif 
      	enddo
	  return
	  end

c This subroutine calculates the pair potential between particles i & j
subroutine Wolf(rij, uij, i, j)
  implicit double precision(a-h,o-z)
  parameter(mp=1000000, mr=2**11, nmq=2**6, nvq=200)
  dimension q(mp), sigma(mp)
  common/electric/np1, np2, nz1, nz2, alpha, sigma1d, sigma2d, sigma12d
  common/box/boxl, rc, np
  pi = dacos(-1.d0)
  
  ! Initialize charges and sigmas for particles
  do l = 1, np1
    q(l) = real(nz1)
    sigma(l) = sigma1d
  enddo
  
  do l = np1 + 1, np
    q(l) = real(nz2)
    sigma(l) = sigma2d
  enddo
  
  ! Calculate error functions and constants
  erfcij = derfc(alpha * rij)
  erfcc = derfc(alpha * rc)
  c2api = 2.d0 * alpha / dsqrt(pi)
  Bcut = erfcc / rc + c2api * dexp(-alpha**2 * rc**2) / rc
  
  ! Calculate Aij
  Aij = q(i) * q(j)
  sigmaij = 0.5d0 * (sigma(i) + sigma(j))
  
  ! Calculate Bij
  Bij = abs(Aij) * (sigmaij**8 * derfc(alpha * sigmaij) &
      + c2api * sigmaij**9 * dexp(-alpha**2 * sigmaij**2) &
      - sigmaij**10 * Bcut) / 9.d0
  
  ! Calculate the pair potential uij
  uij = Bij / rij**9 + Aij * (erfcij / rij - erfcc / rc)
  
  return
end

c This subroutine displace the system to a new configuration
	  subroutine mcmove(x,y,z,ener,nattemp,nacc,del,iseed)
	  implicit double precision(a-h,o-z)
	  parameter(mp=1000000)
	  dimension x(mp),y(mp),z(mp)
      	  common/box/boxl,rc,np
      	  beta=1.96d0
	  temp=1.d0/beta
	  nattemp=nattemp+1
	  no=int(ranf(iseed)*np)+1
	  xo=x(no)
	  yo=y(no)
	  zo=z(no)
	  call energy(x,y,z,xo,yo,zo,enero,no)
	  x(no)=x(no)+(ranf(iseed)-0.5d0)*del
	  y(no)=y(no)+(ranf(iseed)-0.5d0)*del
	  z(no)=z(no)+(ranf(iseed)-0.5d0)*del
c periodic boundary conditions
	  x(no)=x(no)-boxl*dnint(x(no)/boxl)
	  y(no)=y(no)-boxl*dnint(y(no)/boxl)
	  z(no)=z(no)-boxl*dnint(z(no)/boxl)
	  call energy(x,y,z,x(no),y(no),z(no),enern,no)
	  if (ranf(iseed) .lt. dexp(-(enern-enero)*beta)) then
	     ener=ener+0.5d0*(enern-enero)
	     nacc=nacc+1
	     else
	     x(no)=xo
	     y(no)=yo
	     z(no)=zo
	  endif
	  return
	  end

	 subroutine average(x,y,z,g11,g12,g22,nattemp,nacc,del,dr,iseed)
	 implicit double precision(a-h,o-z)
	 parameter(mp=1000000,mr=2**11)
        dimension x(mp),y(mp),z(mp),g11(mr),g12(mr),g22(mr)
	 common/box/boxl,rc,np
	 common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	 nattemp=nattemp+1
	 nbin=0
       ! g11
       do i=1,np1-1
          do j=i+1,np1
             xij=x(j)-x(i)
             yij=y(j)-y(i)
             zij=z(j)-z(i)
             xij=xij-boxl*dnint(xij/boxl)
             yij=yij-boxl*dnint(yij/boxl)
             zij=zij-boxl*dnint(zij/boxl)
             rij2=xij*xij+yij*yij+zij*zij
             rij=dsqrt(rij2)
             if (rij .lt. rc) then
                nbin=dint(rij/dr)
                if (nbin .lt. mr) then
                   g11(nbin)=g11(nbin)+2.
                endif
             endif
          enddo
       enddo
       
       ! g12
       do i=1,np1
          do j=np1+1,np
             xij=x(j)-x(i)
             yij=y(j)-y(i)
             zij=z(j)-z(i)
             xij=xij-boxl*dnint(xij/boxl)
             yij=yij-boxl*dnint(yij/boxl)
             zij=zij-boxl*dnint(zij/boxl)
             rij2=xij*xij+yij*yij+zij*zij
             rij=dsqrt(rij2)
             if (rij .lt. rc) then
                nbin=dint(rij/dr)
                if (nbin .lt. mr) then
                   g12(nbin)=g12(nbin)+1.
                endif
             endif
          enddo
       enddo
       
       !g22
       do i=np1+1,np-1
          do j=i+1,np
             xij=x(j)-x(i)
             yij=y(j)-y(i)
             zij=z(j)-z(i)
             xij=xij-boxl*dnint(xij/boxl)
             yij=yij-boxl*dnint(yij/boxl)
             zij=zij-boxl*dnint(zij/boxl)
             rij2=xij*xij+yij*yij+zij*zij
             rij=dsqrt(rij2)
             if (rij .lt. rc) then
                nbin=dint(rij/dr)
                if (nbin .lt. mr) then
                   g22(nbin)=g22(nbin)+2.
                endif
             endif
          enddo
       enddo
	     nacc=nacc+1
	  return
	  end   

c Random generator algorithm
    FUNCTION ranf(Idum) !Numerical Recipes
	implicit double precision(a-h,o-z)

      PARAMETER (IM1=2147483563, IM2=2147483399,
     & AM=1./IM1, IMM1=IM1-1,IA1=40014, IA2=40692,
     & IQ1=53668, IQ2=52774, IR1=12211,IR2=3791,
     & NTAb=32, NDIv=1+IMM1/NTAb, EPS=1.2E-7, RNMx=1.-EPS)
c
	dimension iv(ntab)
      SAVE iv, iy, idum2
      DATA idum2/123456789/, iv/NTAb*0/, iy/0/
      IF (Idum.LE.0) THEN
         Idum = MAX(-Idum, 1)
         idum2 = Idum
         DO j = NTAb + 8, 1, -1
            k = Idum/IQ1
            Idum = IA1*(Idum-k*IQ1) - k*IR1
            IF (Idum.LT.0) Idum = Idum + IM1
            IF (j.LE.NTAb) iv(j) = Idum
         END DO
         iy = iv(1)
      END IF
      k = Idum/IQ1
      Idum = IA1*(Idum-k*IQ1) - k*IR1
      IF (Idum.LT.0) Idum = Idum + IM1
      k = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2) - k*IR2
      IF (idum2.LT.0) idum2 = idum2 + IM2
      j = 1 + iy/NDIv
      iy = iv(j) - idum2
      iv(j) = Idum
      IF (iy.LT.1) iy = iy + IMM1
      ranf = MIN(AM*iy, RNMx)
      RETURN
      END

c This subroutine adjusts the displacement of particles
	  subroutine adjust(nattemp,nacc,dr)
	  implicit double precision(a-h,o-z)
	  if (mod(nattemp,nacc) .eq. 0) then
	     ratio=real(nacc)/real(nattemp)
	     if (ratio .gt. 0.5) then
	        dr=dr*1.05
	     else
	        dr=dr*0.95
	     endif
	  endif
      return
      end      
