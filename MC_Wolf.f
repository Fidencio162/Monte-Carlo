! The diameter of the particles is used as unit length.
! Ramon Castaneda Priego. 25.09.2018
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

! number of particles
sigmaCl = 3.62d-10
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

! time step
deltat = 0.001d0
! number of configurations to thermalize
nct = 1

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

enertot = 0.d0
do i = 1, np
    call energy(x, y, z, x(i), y(i), z(i), ener, i)
    enertot = enertot + ener
enddo
enertot = 0.5d0 * enertot

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

open(unit=13, file='grnT1SCMCCaCl.dat')
nacco = nacc
ncp = 0
do i = 1, mr
    g11(i) = 0.d0
    g12(i) = 0.d0
    g22(i) = 0.d0
enddo

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
