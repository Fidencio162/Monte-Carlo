! The diameter of the particles is used as unit length.
c Ramon Castaneda Priego. 25.09.2018
      implicit double precision(a-h,o-z)
      logical :: has_overlap
      parameter(mp=1000,mr=2**11,nmq=2**6,nvq=200)
	  dimension x(mp),y(mp),z(mp)
	  dimension qi(nmq),Sq(nmq),s(nmq)
     	 dimension r(mr),g(mr)
     	 dimension sigma(mp)
	common/box/boxl,rc,np
      	common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	common/parameters/diam,rho
	common/poths/dlr,dla,a2,dT
	!common/vectorq/qx(nmq,nvq),qy(nmq,nvq),qz(nmq,nvq)
	!common/configxyz/ x(mp),y(mp),z(mp)
      	pi=4.d0*datan(1.d0)
c number of particles
	  np=10**3
	  
	  open(unit=1,file='sigmacv20.dat')
	   	   sd=0.2d0 !standard deviation
	do i=1,np
	read(1,*) sn,sigma(i)
	enddo
	sigma_max=maxval(sigma)
	   	   sigmain=0.d0
	   	   do i=1,np
	   	   sigmain=sigmain+sigma(i)
	   	   enddo
	   	   sigmap=sigmain/np
	   	   print*, sigmap
	   	   cvs=sd/sigmap !coefficient of variation
	   	   close(1)
c packing fraction
      !	phi=0.034265763d0
      !phi=0.014958d0 ! 1 [M]
      !phi=0.02991613d0 ! 2 [M]
      	rho=0.65d0 !(6.d0/pi)*phi
	  d=(1.d0/rho)**(1./3.)
	  dT=1.4737d0
	  beta=1.d0/dT
	  dlr=50.d0
      dla=49.d0
      dT=1.4737
      a2=dlr*(dlr/dla)**dla
c time step
	   deltat=0.001d0
c number of configurations to thermalize
	   nct=1
c box length
	   boxl=(dfloat(np)/rho)**(1.d0/3.d0)
	   rc=boxl/2.d0
	   dr=rc/mr
	   print*,'The box^3 : ',boxl**3
	   print*,'The mean interparticle distance is: ',d
	   print*,'Cut radius: ',rc
	   print*, 'particle number:', np
	   
	   sigma3i=0.d0
	do k=1,np
	  sigma3i=sigma3i+sigma(k)**3
	enddo
	phi=pi*rho*sigma3i/(6.d0*real(np))
	!phiV=sigma3i/boxl**3
	  phic=pi*rho/6.d0
	  print*, 'phi=',phi,'phi_HS:', phic !,'phiV:',phiV
	  print*,'Volumen total de partículas', sigma3i
	  !stop
!	square charge sum 
!	z2=npanion*(nzanion*nzanion)+npcation*(nzcation*nzcation)
c parameters of the potential (WCA: lambda, T*)
	!alpha=4.d0/boxl
!	square charge sum 
	!z2=npanion*(nzanion*nzanion)+npcation*(nzcation*nzcation)
c initial configuration
	!call init_fcc(x, y, z, sigma, np, np_out, boxl, d*dsqrt(2.38d0))
	  !call init_polydisperse(x, y, z, sigma, np, boxl, 0.85d0)
	  !call iniconfig(x,y,z,d)
	  call iniconfig_random(x,y,z,sigma,boxl,np)
c write the initial configuration
! Verificar exclusión de volumen
	ni=0
	do j=1,np
	xj=x(j)
	yj=y(j)
	zj=z(j)
       	  if (.not. has_overlap(j,xj,yj,zj,boxl,sigma,x,y,z,np)) then
    	  ni=ni+1
	  end if
	enddo
	print*, ni
	!stop
	  open(10,file='iniconf.dat',status='unknown')
	  write(10,49) np
	  write(10,50)
	  do i=1,np1
	  !print*, i
	     write(10,51)x(i),y(i),z(i)
	   enddo
        do i=np1+1,np
       ! print*, i
	     write(10,52)x(i),y(i),z(i)
      enddo
	  close(10)
	  !stop
	  enertot=0.d0
	  do i=1,np
           call energy(x,y,z,x(i),y(i),z(i),sigma,ener,i)
           !print*, i,x(i),y(i),z(i),ener
           enertot=enertot+ener
          enddo
          enertot=0.5d0*enertot
      print*,'Energy per particle of the initial configuration:'
     & ,enertot/real(np)
c MC cycle to thermalize the system
	!stop
      open(30,file='energy12.dat',status='unknown')
!      open(7,file='config.gro')
	  del=0.1d0
	  nattemp=0
	  nacc=1
	  iseed=-7767633
	  
	  do i=1,4000000
	     call mcmove(x,y,z,sigma,ener,nattemp,nacc,del,iseed)
	     call adjust(nattemp,nacc,del)
	     if (mod(i,100) .eq. 0) write(30,*)i*1.d0,ener/np
	     if (mod(i,1000) .eq. 0) then
	     accept=real(nacc)/real(nattemp)
            print*,i,del,ener/real(np),accept
!            write(7,98) i
!             write(7,99) np
!            do j=1,np1
!	     write(7,201)j,x(j),y(j),z(j)
!      	    enddo
!      	    do j=np1+1,np
!	     write(7,202)j,x(j),y(j),z(j)
!      	    enddo
 !     	    write(7,100)ener, ener, ener
         endif
      enddo
	  print*,'The system has thermalized :)'
	  !stop
c write the final configuration and the energy
	  open(20,file='confequi12.dat',status='unknown')
	  !write(20,49) np
	  !write(20,50)
	  do i=1,np
	  !print*, i
	     write(20,*)x(i),y(i),z(i)
	   enddo
        !do i=np1+1,np
       ! print*, i
	!     write(20,52)x(i),y(i),z(i)
      !enddo
	  close(20)
	!  close(30)
	  
	!stop
c MC cycle to calculate the g(r)
      nacco=nacc
      ncp=0
	  do i=1,mr
	  g(i)=0.d0
	  enddo
	  
!	  print*, 'calculating g(r)...'
!	  do i=1,10000000
!	  call mcmove(x,y,z,ener,nattemp,nacc,del,iseed)
!	  if (mod(i,100).eq.0) then
	     !ener=ener/np
	     !ener=ener/real(i)
!	     call average(x,y,z,g,nattemp,nacc,del,dr,iseed)
 !            call adjust(nattemp,nacc,del)
  !           ncp=ncp+1
 !          endif
!	  enddo
	  nav=nacc-nacco
!	  print*,'Average number:',ncp
!	  open(50,file='gr1.dat',status='unknown')
!	  do i=2,mr
!	     r(i)=(i-1)*dr
!	     dv=(4.d0*pi*r(i)**2.*dr)*rho
!	     g(i)=g(i)/(np*ncp*dv)
!	     write(50,200)r(i),g(i)
!	  enddo
!	  close(50)
!	  print*, 'calculating g(r) finished'
!	  stop
	  
	! calculate the structure factor 
	  ! k vector
	  
	  open(unit=22,file='Skran2.dat')
	  dq=pi/rc
	  do i=1,nmq
	    qi(i)=(i-1)*dq
	   ! print*, qi(i)
	  enddo
	  

	  !ncp1=0
	  
	  
	 do i=2,nmq
	    s(i)=0.d0
	 enddo
	 npmc=5000000		!numero de pasos Monte Carlo
	 nmcp=100			!numero de Monte Carlo para promediar
	 print*, 'calculating g(r)...'
	  do j=1,10000000
	    call mcmove(x,y,z,sigma,ener,nattemp,nacc,del,iseed)
	    call adjust(nattemp,nacc,del)
	    if(mod(j,np).eq. 0) then
	    call average(x,y,z,g,nattemp,nacc,del,dr,iseed)
	   ! call adjust(nattemp,nacc,del)
             ncp=ncp+1
	    !print*, s
	    !ncp1=ncp1+1
	    endif
           if (mod(i,500000) .eq. 0) print*,i,del,ener/real(np),accept
	  enddo
	  
	  do k=2,nmq
	    Sq(k)=s(k)/ncp
	    write(22,*) qi(k), Sq(k)
	  enddo
	  close(22)
	!print*, 'calculating g(r)...'
	  open(unit=13,file='gr12.dat')
	  print*,'Average number:',ncp
	  
c calculo de la g11(r),g12(r),g22(r)	 
	  do i=2,mr
	     r(i)=(i-1)*dr
	     dv=(4.d0*pi*r(i)**2.*dr)*rho
	     g(i)=g(i)/(np*ncp*dv)
	     write(13,*)r(i),g(i)
	  enddo
	  close(13)
	  
	  print*, 'calculating g(r) finished'
	  print*, 'finish!'
	  
	  print*, '---------------------------------------------------'
	  print*, '----------FINAL CONFIGURATION----------------------'
	  print*, '---------------------------------------------------'
	  
	 ! do i=1,np
	  !   print*, x(i),y(i),z(i)
	 ! enddo
	  print*, '---------------------------------------------------'
	  
	  print*, 'program finished!'
50	format(x,'CONFIGURACION INICIAL')
51	FORMAT(3X,'C',3f15.7)
52	FORMAT(3X,'N',3f15.7)
49	format(2x,i5)
98	format(x,'CONFIGURACION',11x,i9)
99	format(20x,i5)
	  
100   format(3f15.3)
200   format(2f15.7)
201 	format(5x,'1',2x,'Sol',4x,'C',2x,i5,3f8.3,2x,'0.0000',2x,
     &  '0.0000',2x,'0.0000')
202 	format(5x,'1',2x,'Sol',4x,'N',2x,i5,3f8.3,2x,'0.0000',2x,
     &  '0.0000',2x,'0.0000')
	  end
c
c This subroutine calculates the initial configuration in 3D.
      subroutine iniconfig(x,y,z,d)
      implicit double precision(a-h,o-z)
      parameter(mp=1000)
      dimension x(mp),y(mp),z(mp)
	!common/configxyz/ x(mp),y(mp),z(mp)
	  common/box/boxl,rc,np
	  common/parameters/diam,rho
	  d=rho**(-1.d0/3.d0)
      x(1)=-rc+d/2.d0
      y(1)=-rc+d/2.d0
	  z(1)=-rc+d/2.d0
      do i=1,np-1
         x(i+1)=x(i)+d
         y(i+1)=y(i)
	     z(i+1)=z(i)
         if (x(i+1) .gt. rc) then
            x(i+1)=-rc+d/2.
            y(i+1)=y(i+1)+d
	        z(i+1)=z(i)
	        if (y(i+1) .gt. rc) then
	           x(i+1)=-rc+d/2.
	           y(i+1)=-rc+d/2.
	           z(i+1)=z(i+1)+d
	        endif
	     endif
	   enddo
       return
       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine init_fcc(x, y, z, diam, np, np_out, boxl, a)
    	implicit none
    	integer, intent(in) :: np
    	real(8), intent(in) :: boxl, a
    	real(8), intent(out) :: x(np),y(np),z(np),diam(np)
    	integer, intent(out) :: np_out

    	integer :: nc, i, j, k, n, p
    	real(8), parameter :: rel_pos(3,4) = reshape([
     &   0.0d0, 0.0d0, 0.0d0,   ! posición 1
     &   0.5d0, 0.5d0, 0.0d0,   ! posición 2
     &   0.5d0, 0.0d0, 0.5d0,   ! posición 3
     &   0.0d0, 0.5d0, 0.5d0],  ! posición 4
     &   shape=[3,4])

    	nc = int(boxl / a)
    	np_out = 0
    	p = 0

    	do k = 0, nc-1
          do j = 0, nc-1
            do i = 0, nc-1
                do n = 1, 4
                    p = p + 1
                    if (p > np) then
                        print *, "superado límite de np_max =", np, p
                        stop
                    end if
                    x(p) = (i + rel_pos(1,n)) * a - 0.5d0 * boxl
                    y(p) = (j + rel_pos(2,n)) * a - 0.5d0 * boxl
                    z(p) = (k + rel_pos(3,n)) * a - 0.5d0 * boxl
                    diam(p) = a / sqrt(2.d0)  ! para contacto sin traslape
                end do
            end do
          end do
    	end do

    	np_out = p
	end subroutine init_fcc

	subroutine init_polydisperse(x, y, z, diam, np, boxl, a)
    	implicit none
    	integer, intent(in) :: np
    	real(8), intent(in) :: boxl, a
    	real(8), intent(inout) :: x(np), y(np), z(np), diam(np)

    	integer :: nx, ny, nz, ix, iy, iz, i, j, count, tries
    	real(8) :: dx, dy, dz, rij, rmin, max_tries
    	real(8) :: xi, yi, zi
    	logical :: ok

    	! Dimensiones de la malla
    	nx = int(boxl / a)
    	ny = int(boxl / a)
    	nz = int(boxl / a)
	print*, "nx,ny,nz:",nx,ny,nz
    	if (nx * ny * nz < np) then
        print *,"No hay suficientes sitios para ",np," partículas."
          stop 1
    	end if

    	count = 0
    	do iz = 0, nz-1
          do iy = 0, ny-1
            do ix = 0, nx-1
                xi = (ix + 0.5d0) * a - 0.5d0 * boxl
                yi = (iy + 0.5d0) * a - 0.5d0 * boxl
                zi = (iz + 0.5d0) * a - 0.5d0 * boxl

                ! Verifica que no haya traslape con partículas ya colocadas
                ok = .true.
                do j = 1, count
                    dx = xi - x(j)
                    dy = yi - y(j)
                    dz = zi - z(j)

                    ! Condiciones periódicas en x, y (opcional)
                    dx = dx - boxl * dnint(dx / boxl)
                    dy = dy - boxl * dnint(dy / boxl)
                    dz = dz - boxl * dnint(dz / boxl)

                    rij = sqrt(dx*dx + dy*dy + dz*dz)
                    rmin = 0.5d0 * (diam(j) + diam(count+1))  ! uso diam de la nueva partícula

                    if (rij < rmin) then
                        ok = .false.
                        exit
                    end if
                end do

                if (ok) then
                    count = count + 1
                    x(count) = xi
                    y(count) = yi
                    z(count) = zi
                end if

                if (count == np) exit
            end do
            if (count == np) exit
        end do
        if (count == np) exit
    	end do

    	if (count < np) then
	  print *,"Solo se colocaron",count, "de",np, "partículas"
          stop 1
    	end if
	end subroutine init_polydisperse

      	subroutine iniconfig_random(x,y,z,diam,boxl,np)
      	implicit double precision(a-h,o-z)
    	integer :: max_tries
    	integer, intent(in) :: np
    	real(8), intent(in) :: boxl, diam(np)
    	real(8), intent(inout) :: x(np), y(np), z(np)

    	integer :: i, j, tries
    	logical :: traslapa
    	real(8) :: xi, yi, zi
    	real(8) :: dx, dy, dz, rij, min_dist
    	real(8) :: u
    	max_tries = 50000000

    	do i = 1, np
        tries = 0

        do
            traslapa = .false.

            call random_number(u)
            xi = (u - 0.5d0) * boxl

            call random_number(u)
            yi = (u - 0.5d0) * boxl

            call random_number(u)
            zi = (u - 0.5d0) * boxl

            ! Verifica traslape con partículas ya colocadas
            do j = 1, i - 1
                dx = xi - x(j)
                dy = yi - y(j)
                dz = zi - z(j)

                ! PBC
                dx = dx - boxl * dnint(dx / boxl)
                dy = dy - boxl * dnint(dy / boxl)
                dz = dz - boxl * dnint(dz / boxl)

                rij = dsqrt(dx*dx + dy*dy + dz*dz)
                min_dist = 0.5d0 * (diam(i) + diam(j))

                if (rij < min_dist) then
                    traslapa = .true.
                    exit
                end if
            end do

            tries = tries + 1
            if (tries > max_tries) then
                print *, "ERROR: No se pudo insertar la partícula", i
                stop 1
            end if

            if (.not. traslapa) exit
	end do

        ! Posición aceptada
        x(i) = xi
        y(i) = yi
        z(i) = zi
    	end do
       return
       end
c
c This configuration calculates the energy of a given configuration
	subroutine energy(x,y,z,xj,yj,zj,sigma,ener,j)
	  implicit double precision(a-h,o-z)
	  parameter(mp=1000)
	dimension x(mp),y(mp),z(mp)
	dimension sigma(1000)
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
               call WolfMie(rij,uij,sigma,i,j)
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
               call WolfMie(rij,uij,sigma,i,j)
               ener=ener+uij
            endif 
      	enddo
      !print*, ener
	  return
	  end
c
c This configuration calculates the pair potential between particles i & j
     	subroutine WolfMie(rij,uij,sigma,i,j)
	  implicit double precision(a-h,o-z)
	  parameter(mp=1000,mr=2**11)
	  dimension q(mp),sigma(1000)
	  common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	  common/box/boxl,rc,np
	  common/poths/dlr,dla,a2,dT
	 pi=dacos(-1.d0)
	 dl=50.d0
      
      	sigmaij=0.5d0*(sigma(i)+sigma(j))
      	b2=dlr/dla
      	!print*,rij,b2*sigmaij
      	if (rij .lt. sigmaij*b2) then
         uij=(a2/dT)*((sigmaij/rij)**dlr-(sigmaij/rij)**dla)+1.d0/dT
      	else
         uij=0.d0
      	endif
	  return
	  end
c
c This subroutine displace the system to a new configuration
	  subroutine mcmove(x,y,z,sigma,ener,nattemp,nacc,del,iseed)
	  implicit double precision(a-h,o-z)
	  logical :: has_overlap
	  parameter(mp=1000)
	  dimension sigma(1000)
	  dimension x(mp),y(mp),z(mp)
      	  common/box/boxl,rc,np
      	  common/poths/dlr,dla,a2,dT
	  !beta=1.d0/dT
	  nattemp=nattemp+1
	  no=int(ranf(iseed)*np)+1
	  xo=x(no)
	  yo=y(no)
	  zo=z(no)
	  call energy(x,y,z,xo,yo,zo,sigma,enero,no)
	  xn=x(no)+(ranf(iseed)-0.5d0)*del
	  yn=y(no)+(ranf(iseed)-0.5d0)*del
	  zn=z(no)+(ranf(iseed)-0.5d0)*del
	  ! Verificar exclusión de volumen
       	if (.not. has_overlap(no,xn,yn,zn,boxl,sigma,x,y,z,np)) then
    	  ! Aplica movimiento y periodic boundary conditions
    	  x(no) = xn-boxl*dnint(x(no)/boxl)
    	  y(no) = yn-boxl*dnint(y(no)/boxl)
    	  z(no) = zn-boxl*dnint(z(no)/boxl)
    	  
	  call energy(x,y,z,x(no),y(no),z(no),sigma,enern,no)
	  if (ranf(iseed) .lt. dexp(-(enern-enero))) then
	     ener=ener+0.5d0*(enern-enero)
	     nacc=nacc+1
	     else
	     x(no)=xo
	     y(no)=yo
	     z(no)=zo
	  endif
	  endif
	  return
	  end
c
c This subroutine calculates the g(r)
	 subroutine average(x,y,z,g,nattemp,nacc,del,dr,iseed)
	 implicit double precision(a-h,o-z)
	 parameter(mp=1000,mr=2**11)
        dimension x(mp),y(mp),z(mp),g(mr)
	 common/box/boxl,rc,np
	 common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	 nattemp=nattemp+1
	 nbin=0
       do i=1,np-1
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
c jfix(x) obtiene el entero de x
                nbin=dint(rij/dr)+1
                if (nbin .lt. mr) then
                   g(nbin)=g(nbin)+2.
                endif
             endif
          enddo
       enddo
	     nacc=nacc+1
	  return
	  end

c Random generator algorithm
      FUNCTION ranf(Idum)
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
      
      	FUNCTION gasdev(idum)
	implicit double precision(a-h,o-z)
	SAVE iset,gset
	DATA iset/0/
!	  if (idum.lt.0) iset=0
	  if (iset.eq.0) then
   1    v1=2.d0*ranf(idum)-1.d0
	  v2=2.d0*ranf(idum)-1.d0
	  rsq=v1**2+v2**2
	  if((rsq.ge.(1.d0)).or.(rsq.eq.(0.d0)))goto 1 
	  fac=dsqrt(-2.d0*dlog(rsq)/rsq)
	  gset=v1*fac
	  gasdev=v2*fac
	  iset=1
	  else
	  gasdev=gset
	  iset=0
	endif
	return
	END  
      
c      	function ranf(idum)
c	implicit integer*8(i-n),real*8(a-h,o-z)
c
c      parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
c     *MASK=123459876)
c      idum=ieor(idum,MASK)
c      k=idum/IQ
c      idum=IA*(idum-k*IQ)-IR*k
c      if (idum.lt.0) idum=idum+IM
c      ranf=AM*idum
c      idum=ieor(idum,MASK)
c      return
c      end
c
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
      
      subroutine Sk(x,y,z,s,qx,qy,qz)
	implicit integer*8(i-n),real*8(a-h,o-z)
	parameter(mp=1000000,nmq=2**6,nvq=200)
	dimension x(mp),y(mp),z(mp),s(nmq)
	dimension auxc(nvq),auxs(nvq)
	dimension qx(nmq,nvq),qy(nmq,nvq),qz(nmq,nvq)
	common/box/boxl,rc,np

	do i=2,nmq
	  sum  = 0.d0
	  parti =0.d0

	  auxc = 0.d0
	  auxs = 0.d0

	  do k=1,np
	  !print*, x(k)
         xaux = x(k) - boxl*dnint(x(k)/boxl)
         yaux = y(k) - boxl*dnint(y(k)/boxl)
         zaux = z(k) - boxl*dnint(z(k)/boxl)
         rij2 = xaux*xaux + yaux*yaux + zaux*zaux
         rij  = dsqrt(rij2)

         if (rij .lt. rc) then
            parti=parti+1.d0
            do j = 1,nvq
               arg = qx(i,j)*xaux + qy(i,j)*yaux + qz(i,j)*zaux
               auxc(j) = auxc(j)+dcos(arg)
               auxs(j) = auxs(j)+dsin(arg)
            enddo
         endif
	  enddo


	  do l=1,nvq
            sum = sum + auxc(l)*auxc(l) + auxs(l)*auxs(l)
	  enddo

         auxsq = sum/(nvq*parti)
         s(i)=s(i)+auxsq
	enddo
	return
	end
	
	logical function has_overlap(i,xnew,ynew,znew,box,diam,x,y,z,N)
	implicit double precision(a-h,o-z)
    	real(8), intent(in) :: xnew, ynew, znew
    	real(8), intent(in) :: x(N), y(N), z(N), diam(N)
    	integer, intent(in) :: i, N
    	real(8) :: dx, dy, dz, rij, sig

    	has_overlap = .false.
    	do j = 1, N
          if (j /= i) then
            dx = xnew - x(j)
            dy = ynew - y(j)
            dz = znew - z(j)
            ! Aplica condiciones periódicas
            dx = dx - box * dnint(dx / box)
            dy = dy - box * dnint(dy / box)
            dz = dz - box * dnint(dz / box)
            rij = dsqrt(dx*dx + dy*dy + dz*dz)
            sig = 0.5d0*(diam(i) + diam(j))
            if (rij < sig- 1d-12) then
!          print *, "Traslape detectado entre partículas ",i, " y ",j,
!     &          " | r_ij = ", rij, " < sig = ", sig
                has_overlap = .true.
                return
            end if
          end if
        end do
	end function

