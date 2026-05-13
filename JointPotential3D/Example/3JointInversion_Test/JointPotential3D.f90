! Joint Potential 1.0 Package
! Jonathan Carrillo, CICESE
!
! REVISION HISTORY:
! REFERENCES: CARRILLO & GALLARDO, 2018: GEOPHYSICAL JOURNAL INTERNATIONAL
!01-Sep-2018 Modificacion para alphas y Cg dinamicos 
!12-Agosto-2019 Modificacion para correr en computadora de escritorio 
!Se añadieron Modelos a priori, desv estándar por separado y gradientes cruzados          
Program inv3DconjGraviMagneto
implicit none
integer::i,j,k,nc,nc1,nc2,Nd1,Nd2,ncoef,conta
integer::h,l,Action,Method,Sensitivity,StdDev,Apriori,nX,nY,nZ
double precision :: rms,error,alpha1,error1,error2,error3,rms1,rms2,rms3
double precision :: error4,error5,error6,error7,error8,rms4,rms5,rms6,rms7,rms8,pi
double precision :: alpha1X,alpha1Y,alpha1Z,alpha2X,alpha2Y,alpha2Z
double precision ::delt_modelo,delt_m1rms,delt_m2rms,delt_arms,delt_mrms
double precision :: Cdd1,Cdd2,Cgg,x1,y1,z1,z2,x2,y2, mi,md,fi,fd,theta,gm1,gm2,factor1,factor2,conta_real
double precision , allocatable, dimension(:)::delt_m1,delt_m2,delt_a
double precision , allocatable, dimension(:)::x0
double precision ,allocatable, dimension(:)::x01,y01,z01,x02,y02,z02
double precision , allocatable, dimension(:,:)::Ad1,Ad2,Ah1,Ah2,Cd1,Cd2,Cmm,BB,PP
!A partir de aqui son vectores
double precision , allocatable, dimension(:,:)::a0,g_obs,d1,d2,d2aux,m01,m02,g0,d10,d20
double precision , allocatable, dimension(:,:)::vecCd1,vecCd2,vecCg
double precision , allocatable, dimension(:,:)::vecCmm,Ad1t,Ag1t,Ad2t,Ag2t,Aat

double precision , allocatable, dimension(:,:)::Dfm1,Dfm2,Cg,Ag1,Ag2,Aa
double precision , allocatable, dimension(:,:)::ter11a,ter11c
double precision , allocatable, dimension(:,:)::ter22a1,ter22a,ter22c
double precision , allocatable, dimension(:,:)::ter12a1,ter21a1
double precision , allocatable, dimension(:,:)::B11,B22,B33,B12,B13,B21,B23,B31,B32
double precision , allocatable, dimension(:,:)::Ad1m01,Ag1m01,Ag2m02,Aaa0,Ad2m02
double precision , allocatable, dimension(:,:)::dat1a,dat1c,dat2a,dat2c,dat3
double precision , allocatable, dimension(:,:)::P1a,P1c,P1
double precision , allocatable, dimension(:,:)::P2a,P2c,P2,P3
double precision , allocatable, dimension(:,:)::invB,param,sol_m1,sol_m2,sol_a
double precision , allocatable, dimension(:,:)::g1,d11,d21
double precision , allocatable, dimension(:,:)::m01pr,vecCm01pr,vecCm02,Cm01pr,Cm02,invCm01pr
double precision , allocatable, dimension(:,:)::m02pr,vecCm02pr,Cm02pr,invCm02pr
double precision ,allocatable, dimension(:,:)::Dz,Dzt,DztDz,Dy,Dyt,DytDy,Dx,Dxt,DxtDx, malla
character(len=20) :: str
real::tiempo

interface

subroutine const_mat_diag(n,vecDiag,Mat_diag)
!Construye una matriz diagonal del vector que se le da, vecDiag, de n elementos
implicit none
integer, intent(in)::n
double precision, intent(in)::vecDiag(n,1)
double precision ,intent(out)::Mat_diag(n,n)
integer::i,j
end subroutine const_mat_diag

Subroutine Agm1_sin_a00_pol2d(nc,m1,m2,coef,Dfm1)
implicit none
integer, intent(in)::nc
double precision, intent(in)::m1(nc,1),m2(nc,1),coef(8,1)
double precision ,intent(out)::Dfm1(nc,1)
double precision :: Agm1(nc,8)
integer::i,j
end subroutine Agm1_sin_a00_pol2d

Subroutine Agm2_sin_a00_pol2d(nc,m1,m2,coef,Dfm2)
implicit none
integer, intent(in)::nc
double precision, intent(in)::m1(nc,1),m2(nc,1),coef(8,1)
double precision ,intent(out)::Dfm2(nc,1)
double precision :: Agm2(nc,8)
integer::i,j
end subroutine Agm2_sin_a00_pol2d

Subroutine Aga_sin_a00_pol2d(nc,m1,m2,Aga)
implicit none
integer, intent(in)::nc
double precision, intent(in)::m1(nc,1),m2(nc,1)
double precision ,intent(out)::Aga(nc,8)
integer::i,j
end subroutine Aga_sin_a00_pol2d

Subroutine Const_B(nc1,nc2,ncoef,B11,B12,B13,B21,B22,B23,B31,B32,B33,B)
!Esta subrutina construye B, de la formulacion Bm=P
implicit none
integer, intent(in)::nc1,nc2,ncoef
double precision, intent(in)::B11(nc1,nc1),B22(nc2,nc2),B33(ncoef,ncoef),B12(nc1,nc2),B13(nc1,ncoef)
double precision, intent(in)::B21(nc2,nc1),B23(nc2,ncoef),B31(ncoef,nc2),B32(ncoef,nc2)
double precision ,intent(out)::B(nc1+nc2+ncoef,nc1+nc2+ncoef)
integer::i,j
end subroutine Const_B

subroutine Const_P(nc1,nc2,ncoef,P1,P2,P3,P)
!Esta subrutina construye P, de la formulacion Bm=P
implicit none
integer, intent(in)::nc1,nc2,ncoef
double precision, intent(in)::P1(nc1,1),P2(nc2,1),P3(ncoef,1)
double precision ,intent(out)::P(nc1+nc2+ncoef,1)
integer::i,j
end subroutine Const_P

subroutine inverse(a,c,n)
implicit none
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k
end subroutine inverse

subroutine traspuesta(n1,n2,A,B)
!Calcula la traspuesta B de una matriz A de n1xn2
implicit none
integer, intent(in)::n1,n2
double precision, intent(in)::A(n1,n2)
double precision ,intent(out)::B(n2,n1)
integer::i,j
end subroutine traspuesta


subroutine gbox(x0,y0,z0,Ndat,Nprism,malla,g)
!Modificado para modelo de capas
implicit none
! SUBRUTINA QUE CALCULA LA ATRACCI√ìN VERTICAL DE UN PRISMA RECT√ÅNGULAR
! Los lados del prisma son paralelos al eje x, y y z; donde z aumenta
! hacia abajo.
! Los par¬†metros a ingresar son:
! el punto de observaci¬¢n (x0,y0,z0).
! El prima se extiende de x1 a x2, de y1 a y2 y de z1 a z2.
! La densidad del prisma es rho en unidades kg/m^3
! Las distancias estan dadas en km.
! La atraccion vertical de gravedad g en mGal.
integer, intent(in)::Ndat,Nprism
double precision ,intent(in)::x0(Ndat),y0(Ndat),z0(Ndat),malla(Nprism,6)
double precision , dimension(:,:), intent(out)::g
integer::i,j,k,h,l
double precision ::sum,rijk,ijk,arg1,arg2,arg3
double precision ,dimension(2)::x,y,z,isign
double precision ,parameter::gamma=6.67e-11,twopi=6.2831853,si2mg=1.e5,km2m=1.e3
end subroutine gbox


subroutine oper_derX(n,Dx)
implicit none
integer, intent(in)::n
double precision ,intent(out)::Dx(n,n)
integer::i,j
end subroutine oper_DerX

subroutine oper_derY(n,nx,ny,nz,Dy)
implicit none
integer, intent(in)::n,nx,ny,nz
double precision ,intent(out)::Dy(n,n)
integer::i,j,k
end subroutine oper_derY

subroutine oper_derZ(n,nx,ny,nz,Dz)
implicit none
integer, intent(in)::n,nx,ny,nz
double precision ,intent(out)::Dz(n,n)
integer::i,j,k
end subroutine oper_derZ

 subroutine mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,theta,gm)
      double precision ,intent(in)::x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,theta
      double precision ,intent(out)::gm
      double precision :: alpha(2),beta(2),ma,mb,mc,t,fa,fb,fc
      double precision :: fm1,fm2,fm3,fm4,fm5,fm6,hsq,alphasq,sign
      double precision :: r0sq,r0h,alphabeta,arg3,arg4,tlog,tatan
      real :: arg1,arg2, h,r0 !alog solo permite reales
      !data cm/1.e-7/,t2nt/1.e9/
      double precision ,parameter::cm=1.e-7,t2nt=1.e9
end subroutine mbox


subroutine dircos(incl,decl,azim,a,b,c)
      double precision ,intent(in)::incl,decl,azim
      double precision ,intent(out)::a,b,c
      double precision :: xincl,xdecl,xazim
      double precision ,parameter::d2rad=.017453293
 end subroutine dircos

END INTERFACE

Write(*,*) '------------------------------------------------------------'
Write(*,*) 'JointPotential Inversion v1.0                    August 2019'
Write(*,*) '------------------------------------------------------------'
Write(*,*) '------------------------------------------------------------'
Write(*,*) '1 Forward Modelling Gravity and Magnetic data 3D'
Write(*,*) '------------------------------------------------------------'
Write(*,*) '2 Inverse Modelling Gravity and Magnetic data 3D'
Write(*,*) '------------------------------------------------------------'
Write(*,*) '3 Joint Modelling Gravity and Magnetic data 3D'
Write(*,*) '------------------------------------------------------------'
Write(*,*) '4 Cross Gradients Measure'
Write(*,*) '------------------------------------------------------------'
!STRUCTURE
! 1 LReading data from external files: 196:760
!                                      223 Startup FWD
!                                      275 Startup INV
!                                      353 Startup JOINT
! 2 Perform process selected: FWD lines 750:860
!                             INV lines 750:860
! Cobnstruyo B y P
 pi=4.0*atan(1.0)
call cpu_time(tiempo) 
write (*, *) tiempo 

!!!!!!!!!!!!!!!!!!!!!!FIRST EXTERNAL FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(*,*) 'Reading startup file from startup.dat'
 !STARTUP FILE
!!
!!!!Startup contains all information about the procedue to perform
!!!! Forward modelling, inverse modelling, joint inversion or cross gradients
open(unit=1,file='startup.dat')

read(1,*) Action !Type of modelling 
                 !1 Forward Modelling
                 !2 Inverse Modelling
                 !3 Joint inversion
                 !4 Cross Gradient measure           

select case (Action)                                         !!!!STARTS ACTION1
case(1)                                
!!!!1FORWARD MODELLING SECTION
write(*,*)'Reading Forward Modelling Sartup'
read(1,*) Method !Type of method 
                 !1 Gravity
                 !2 Magnetic
                 
!!!1.1 Gravity data
if (Method==1) then
write(*,*)'Gravity data'
read(1,*) Nd1
read(1,*) Nc1
read(1,*) Nx
read(1,*) Ny
read(1,*) Nz

write(*,*)'Allocate gravity forward matrixes'
allocate(m01(nc1,1))
allocate(d1(Nd1,1))
allocate(Ad1(Nd1,nc1))
allocate(Ad1t(nc1,nd1))
allocate(d11(nd1,1))
allocate(malla(Nc1,6))

allocate(x01(Nd1),y01(Nd1),z01(Nd1))

!!!1.2 Magnetic data
else if (Method==2) then
write(*,*)'Magnetic data'
read(1,*) mi !mi
read(1,*) md !md
read(1,*) fi !fi
read(1,*) fd !fd
read(1,*) theta !theta
read(1,*) Nd2
read(1,*) Nc2
read(1,*) Nx
read(1,*) Ny
read(1,*) Nz

write(*,*)'Nd2'
write(*,*) Nd2
write(*,*)'Nc2'
write(*,*) Nc2

write(*,*)'Allocate magnetic forward matrixes'
allocate(m02(nc2,1))
allocate(d2(Nd2,1))
allocate(Ad2(Nd2,nc2))
allocate(Ad2t(nc2,nd2))
allocate(d21(nd2,1))
allocate(malla(Nc2,6))

allocate(x02(Nd2),y02(Nd2),z02(Nd2))

else
write(*,*)'Invalid option for method'
end if

case(2)
!!! 2INVERSE MODELLING
write(*,*)' Reading Inverse Modelling Startup'
read(1,*) Method !Type of method 
                 !1 Gravity
                 !2 Magnetic
                 
!!!1.1 Gravity data
if (Method==1) then
write(*,*)'Gravity data'
read(1,*) Apriori !A priori
!Smothness factors for density model for the 3 directions
read(1,*) alpha1X !Alpha1X
read(1,*) alpha1Y !Alpha1Y
read(1,*) alpha1Z !Alpha1Z
read(1,*) Nd1
read(1,*) Nc1
read(1,*) Nx
read(1,*) Ny
read(1,*) Nz


!!!!Allocate gravity matrixes

write(*,*)'Allocate gravity inversion matrixes'
allocate(d1(Nd1,1))
allocate(Ad1(Nd1,nc1),Cd1(nd1,nd1))
allocate(ter11a(nc1,nc1),P1a(nc1,1),vecCd1(nd1,1))
allocate(BB(nc1,nc1),PP(nc1,1))
allocate(Ad1t(nc1,nd1))
allocate(invB(nc1,nc1),param(nc1,1),sol_m1(nc1,1))
allocate(d11(nd1,1))
allocate(m01pr(nc1,1),vecCm01pr(nc1,1),Cm01pr(nc1,nc1),invCm01pr(nc1,nc1))
allocate(malla(Nc1,6))

allocate(x01(Nd1),y01(Nd1),z01(Nd1))
allocate(Dx(Nc1,Nc1),Dy(Nc1,Nc1),Dz(Nc1,Nc1),Dxt(Nc1,Nc1),Dyt(Nc1,Nc1),Dzt(Nc1,Nc1))
allocate(DxtDx(nc1,nc1),DytDy(nc1,nc1),DztDz(nc1,nc1))

!!!1.2 Magnetic data
else if (Method==2) then
write(*,*)'Magnetic data'
read(1,*) Apriori !A priori
!Smothness factors for magnetization model for the 3 directions
read(1,*) alpha2X !Alpha2X
read(1,*) alpha2Y !Alpha2Y
read(1,*) alpha2Z !Alpha2Z
read(1,*) mi !mi
read(1,*) md !md
read(1,*) fi !fi
read(1,*) fd !fd
read(1,*) theta !theta
read(1,*) Nd2
read(1,*) Nc2
read(1,*) Nx
read(1,*) Ny
read(1,*) Nz

write(*,*)'FLAG1 Nz'
write(*,*) Nz
!!!!Allocate gravity matrixes

!!!allocate joint inversion parameters
write(*,*)'Allocate magnetic inversion matrixes'
allocate(d2(Nd2,1))
allocate(Ad2(Nd2,nc2),Cd2(nd2,nd2))
allocate(ter22a(nc2,nc2),P2a(nc2,1),vecCd2(nd2,1))
allocate(BB(nc2,nc2),PP(nc2,1))
allocate(Ad2t(nc2,nd2))
allocate(invB(nc2,nc2),param(nc2,1),sol_m2(nc2,1))
allocate(d21(nd2,1))
allocate(m02pr(nc2,1),vecCm02pr(nc2,1),Cm02pr(nc2,nc2),invCm02pr(nc2,nc2))
allocate(malla(Nc2,6))

allocate(x02(Nd2),y02(Nd2),z02(Nd2))
allocate(Dx(Nc2,Nc2),Dy(Nc2,Nc2),Dz(Nc2,Nc2),Dxt(Nc2,Nc2),Dyt(Nc2,Nc2),Dzt(Nc2,Nc2))
allocate(DxtDx(nc2,nc2),DytDy(nc2,nc2),DztDz(nc2,nc2))
else
write(*,*)'Invalid option for method'
end if     


!!! 3JOINT INVERSION
case(3)
write(*,*)'Reading Joint inversion Startup'
read(1,*) Apriori !A priori
!Smothness factors for density model for the 3 directions
read(1,*) alpha1X !Alpha1X
read(1,*) alpha1Y !Alpha1Y
read(1,*) alpha1Z !Alpha1Z
!Smothness factors for magnetization model for the 3 directions
read(1,*) alpha2X !Alpha2X
read(1,*) alpha2Y !Alpha2Y
read(1,*) alpha2Z !Alpha2Z

read(1,*) mi !mi
read(1,*) md !md
read(1,*) fi !fi
read(1,*) fd !fd
read(1,*) theta !theta
read(1,*) Nd1
read(1,*) Nd2 
read(1,*) nc1
read(1,*) nc2
read(1,*) nc
read(1,*) ncoef
read(1,*) Nx
read(1,*) Ny
read(1,*) Nz

!!!allocate joint inversion parameters
write(*,*)'Allocate Joint inversion matrixes'
allocate(d1(Nd1,1),d2(Nd2,1),g_obs(nc,1))
allocate(Ad1(Nd1,nc1),Ad2(Nd2,nc2),Cd1(nd1,nd1),Cd2(nd2,nd2))
allocate(a0(ncoef,1),m01(nc1,1),m02(nc2,1))
allocate(DxtDx(nc,nc),DytDy(nc,nc),DztDz(nc,nc))
allocate(vecCg(nc,1))
allocate(vecCd1(nd1,1),vecCd2(nd2,1))
allocate(Dfm1(nc,1),Dfm2(nc,1),Cg(nc,nc))
allocate(Ag1(nc,nc1),Ag2(nc,nc2),Aa(nc,ncoef))
allocate(g0(nc,1),d10(nd1,1),d20(nd2,1),BB(nc1+nc2+ncoef,nc1+nc2+ncoef),PP(nc1+nc2+ncoef,1))
allocate(Ad1t(nc1,nd1),Ag1t(nc1,nc))
allocate(Ad2t(nc2,nd2),Ag2t(nc2,nc),Aat(ncoef,nc))
allocate(ter11a(nc1,nc1),ter11c(nc1,nc1))
allocate(ter22a1(nc2,nd2),ter22a(nc2,nc2),ter22c(nc2,nc2))
allocate(B11(nc1,nc1),B22(nc2,nc2),B33(ncoef,ncoef),B12(nc1,nc2),B13(nc1,ncoef),B21(nc2,nc1),B23(nc2,ncoef))
allocate(B31(ncoef,nc1),B32(ncoef,nc2))
allocate(Ad1m01(nd1,1),Ag1m01(nc,1),Ag2m02(nc,1),Aaa0(nc,1),Ad2m02(nd2,1))
allocate(dat1a(nd1,1),dat1c(nc,1),dat2a(nd2,1),dat2c(nc,1),dat3(nc,1))
allocate(P1a(nc1,1),P1c(nc1,1),P1(nc1,1))
allocate(P2a(nc2,1),P2c(nc2,1),P2(nc2,1))
allocate(P3(ncoef,1))
allocate(invB(nc1+nc2+ncoef,nc1+nc2+ncoef),param(nc1+nc2+ncoef,1),sol_m1(nc1,1),sol_m2(nc2,1),sol_a(ncoef,1))
allocate(g1(nc,1),d11(nd1,1),d21(nd2,1))
allocate(m01pr(nc1,1),vecCm01pr(nc1,1),Cm01pr(nc1,nc1),invCm01pr(nc1,nc1))
allocate(m02pr(nc2,1),vecCm02pr(nc2,1),Cm02pr(nc2,nc2),invCm02pr(nc2,nc2))
allocate(malla(Nc,6))


!!DATA
allocate(x01(Nd1),y01(Nd1),z01(Nd1))
allocate(x02(Nd2),y02(Nd2),z02(Nd2))
allocate(Dx(Nc1,Nc1),Dy(Nc1,Nc1),Dz(Nc1,Nc1),Dxt(Nc1,Nc1),Dyt(Nc1,Nc1),Dzt(Nc1,Nc1))

!!! 4CROSS GRADIENTS MEASURE
case(4)
write(*,*)'Cross gradient measure'

CASE DEFAULT
      WRITE(*,*) 'Invalid option for module'
END SELECT                                          !!!!ENDS ACTION1

close(1)




!!!!!!!!!!!!!!!!!!!!!!SECOND EXTERNAL FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READING MESH!!!!!
 write(*,*) 'Reading mesh file from mesh.dat'
select case (Action)                                         !!!!STARTS ACTION2

!!!!1FORWARD MODELLING SECTION
case(1)
write(*,*)'Reading Forward Modelling Mesh'

!!!1.1 Gravity data
if (Method==1) then
write(*,*)'Gravity data'

open(unit=2,file='mesh.dat')
do i=1,Nc1
read(2,*) (malla(i,j), j=1,6)
end do
close(2)

!!!1.2 Magnetic data
else if (Method==2) then
write(*,*)'Magnetic data'

open(unit=3,file='mesh.dat')
do i=1,Nc2
read(3,*) (malla(i,j), j=1,6)
end do
close(3)

else
write(*,*)'Invalid option for method'
end if

!!! 2INVERSE MODELLING
case(2)
write(*,*)'Reading Inverse Modelling Mesh'
!!!1.1 Gravity data
if (Method==1) then
write(*,*)'Gravity data'

open(unit=4,file='mesh.dat')
do i=1,Nc1
read(4,*) (malla(i,j), j=1,6)
end do
close(4)

!!!1.2 Magnetic data
else if (Method==2) then
write(*,*)'Magnetic data'

open(unit=5,file='mesh.dat')
do i=1,Nc2
read(5,*) (malla(i,j), j=1,6)
end do
close(5)

else
write(*,*)'Invalid option for method'
end if
                 
!!! 3JOINT INVERSION
case(3)
write(*,*)'Reading Joint Inversion Mesh'

open(unit=66,file='mesh.dat')
do i=1,Nc
read(66,*) (malla(i,j), j=1,6)
end do
close(66)

!!! 4CROSS GRADIENTS MEASURE
case(4)
write(*,*)'Cross gradient measure mesh'

open(unit=7,file='mesh.dat')
do i=1,Nc
read(7,*) (malla(i,j), j=1,6)
end do
close(7)


CASE DEFAULT
      WRITE(*,*) 'Invalid option for module'
END SELECT                                         !!!!ENDS ACTION2

!!!!!READING MESH ENDS


!!!!!!!!!!!!!!!!!!!!!!THIRD EXTERNAL FILE: DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READING DATA!!!!!
 write(*,*) 'Reading data files (inversion) or models (forward modelling)'
select case (Action)                                         !!!!STARTS ACTION3
!!!!1FORWARD MODELLING SECTION
case(1)
write(*,*)'Reading Forward Modelling Models'

!!!1.1 Gravity data
if (Method==1) then
write(*,*)'m01.dat'

open(unit=8,file='m01.dat')
do i=1,Nc1
read(8,*) m01(i,1)
end do
close(8)

write(*,*)'Gravity UTM locations (km)'
write(*,*)'Important: the program takes the original UTM axes' 
write(*,*) 'but internally X is the North, Y is the Easth and Z positive at depth'
!!!en km
open(unit=99,file='X_UTM_grav.dat')
do i=1,Nd1
read(99,*) y01(i)
end do
close(99)
  
  !POR QUE NO FUNCIONA LA UNIDAD 6?
open(unit=10,file='Y_UTM_grav.dat')
do i=1,Nd1
read(10,*) x01(i)
end do
close(10)

open(unit=11,file='Z_UTM_grav.dat')
do i=1,Nd1
read(11,*) z01(i)
end do
close(11)

z01=-z01 !negative to depth


!!!1.2 Magnetic data
else if (Method==2) then
write(*,*)'m02.dat'
open(unit=12,file='m02.dat')
do i=1,Nc2
read(12,*) m02(i,1)
end do
close(12)


write(*,*)'Important: the program takes the original UTM axes' 
write(*,*) 'but internally X is the North, Y is the Easth and Z positive at depth'
!!!en km
open(unit=13,file='X_UTM_mag.dat')
do i=1,Nd2
read(13,*) y02(i)
end do
close(13)

open(unit=14,file='Y_UTM_mag.dat')
do i=1,Nd2
read(14,*) x02(i)
end do
close(14)

open(unit=15,file='Z_UTM_mag.dat')
do i=1,Nd2
read(15,*) z02(i)
end do
close(15)

z02=-z02 !negative to depth

else
write(*,*)'Invalid option for method'
end if

!!! 2INVERSE MODELLING
case(2)
write(*,*)'Reading Inverse Modelling Data'
                 
!!!1.1 Gravity data
if (Method==1) then
write(*,*)'Gravity.dat'

!!!!!!!!!!!!!! 2 READING GRAVITY DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=16,file='Gravity.dat')
do i=1,Nd1
read(16,*) d1(i,1)
end do
close(16)

open(unit=17,file='Cd1.dat')
do i=1,Nd1
read(17,*) vecCd1(i,1)
end do
close(17)
write(*,*)'Gravity UTM locations (km)'
write(*,*)'Important: the program takes the original UTM axes' 
write(*,*) 'but internally X is the North, Y is the Easth and Z positive at depth'
!!!en km
open(unit=18,file='X_UTM_grav.dat')
do i=1,Nd1
read(18,*) y01(i)
end do
close(18)
  
  !POR QUE NO FUNCIONA LA UNIDAD 6?
open(unit=19,file='Y_UTM_grav.dat')
do i=1,Nd1
read(19,*) x01(i)
end do
close(19)

open(unit=20,file='Z_UTM_grav.dat')
do i=1,Nd1
read(20,*) z01(i)
end do
close(20)

z01=-z01 !negative to depth


!!!1.2 Magnetic data
else if (Method==2) then
write(*,*)'Magnetic.dat'
!!!!!!!!!!!!!! 3 READING MAGNETIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=21,file='Magnetic.dat')
do i=1,Nd2
read(21,*) d2(i,1)
end do
close(21)


open(unit=22,file='Cd2.dat')
do i=1,Nd2
read(22,*) vecCd2(i,1)
end do
close(22)

write(*,*)'FLAG2 Cd2'
write(*,*) vecCd2(3,1)

write(*,*)'Magnetic UTM locations (km)'
write(*,*)'Important: the program takes the original UTM axes' 
write(*,*) 'but internally X is the North, Y is the Easth and Z positive at depth'
!!!en km
open(unit=23,file='X_UTM_mag.dat')
do i=1,Nd2
read(23,*) y02(i)
end do
close(23)

open(unit=24,file='Y_UTM_mag.dat')
do i=1,Nd2
read(24,*) x02(i)
end do
close(24)

open(unit=25,file='Z_UTM_mag.dat')
do i=1,Nd2
read(25,*) z02(i)
end do
close(25)

z02=-z02 !negative to depth

else
write(*,*)'Invalid option for method'
end if
                 


!!! 3JOINT INVERSION
case(3)
write(*,*)'Reading Joint Inversion Data'
!!!!!!!!!!!!!! 2 READING GRAVITY DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)'Gravity.dat'
open(unit=26,file='Gravity.dat')
do i=1,Nd1
read(26,*) d1(i,1)
end do
close(26)

open(unit=27,file='Cd1.dat')
do i=1,Nd1
read(27,*) vecCd1(i,1)
end do
close(27)

write(*,*)'Gravity UTM locations (km)'
write(*,*)'Important: the program takes the original UTM axes' 
write(*,*) 'but internally X is the North, Y is the Easth and Z positive at depth'
!!!en km
open(unit=28,file='X_UTM_grav.dat')
do i=1,Nd1
read(28,*) y01(i)
end do
close(28)
  
  !POR QUE NO FUNCIONA LA UNIDAD 6?
open(unit=29,file='Y_UTM_grav.dat')
do i=1,Nd1
read(29,*) x01(i)
end do
close(29)

 write(*,*) 'x01(2)'
  write(*,*) x01(2)

open(unit=30,file='Z_UTM_grav.dat')
do i=1,Nd1
read(30,*) z01(i)
end do
close(30)

z01=-z01 !negative to depth


!!!!!!!!!!!!!! 3 READING MAGNETIC DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)'Magnetic.dat'
open(unit=31,file='Magnetic.dat')
do i=1,Nd2
read(31,*) d2(i,1)
end do
close(31)

open(unit=32,file='Cd2.dat')
do i=1,Nd2
read(32,*) vecCd2(i,1)
end do
close(32)

write(*,*)'Magnetic UTM locations (km)'
write(*,*)'Important: the program takes the original UTM axes' 
write(*,*) 'but internally X is the North, Y is the Easth and Z positive at depth'
!!!en km
open(unit=33,file='X_UTM_mag.dat')
do i=1,Nd2
read(33,*) y02(i)
end do
close(33)

open(unit=34,file='Y_UTM_mag.dat')
do i=1,Nd2
read(34,*) x02(i)
end do
close(34)

open(unit=35,file='Z_UTM_mag.dat')
do i=1,Nd2
read(35,*) z02(i)
end do
close(35)

write(*,*) 'Nc'
write(*,*) Nc

open(unit=36,file='Cgg.dat')
do i=1,Nc
read(36,*) vecCg(i,1)
end do
close(36)

z02=-z02 !negative to depth

!!! 4CROSS GRADIENTS MEASURE
case(4)
write(*,*)'Cross gradient measure mesh'


CASE DEFAULT
      WRITE(*,*) 'Invalid option for module'
END SELECT                                           !!!!ENDS ACTION3
!!!!!READING DATA ENDS READING DATA ENDS
!!!!!READING DATA ENDS READING DATA ENDS
!!!!!READING DATA ENDS READING DATA ENDS
!!!!!READING DATA ENDS READING DATA ENDS
!!!!!READING DATA ENDS READING DATA ENDS


Write(*,*) '------------------------------------------------------------'
Write(*,*) 'JointPotential Inversion v1.0                    August 2019'
Write(*,*) '------------------------------------------------------------'
Write(*,*) '------------------------------------------------------------'
Write(*,*) 'READING STEP FINISHED'
Write(*,*) '------------------------------------------------------------'
Write(*,*) 'STARTING THE PROCESS SELECTED'
Write(*,*) '------------------------------------------------------------'

select case (Action)                                         !!!!STARTS ACTION4
case(1)
Write(*,*) 'FWD'
if (Method==1) then

WRITE(*,*) 'Constructing gravity sensitivity matrix'
call gbox(x01,y01,z01,Nd1,Nc1,malla,Ad1)
!end do
 Ad1=Ad1*1000 !to gm/cm3

!WRITE(*,*) 'Writting GG.dat'

!open(unit=9,file='GG.dat',status='replace',action='write')
!do i=1,Nd1
!write(9,*) (Ad1(i,j), j=1,Nc1)
!end do
!close(9)

WRITE(*,*) 'Writing gravity data'

d1=matmul(Ad1,m01)

open(unit=38,file='Synthetic_gravity.dat',status='replace',action='write')
do i=1,nd1
write(38,*) d1(i,1)
end do
close(38)


write(*,*)'Deallocate gravity forward matrixes'
deallocate(m01)
deallocate(d1)
deallocate(Ad1)
deallocate(Ad1t)
deallocate(d11)
deallocate(malla)
deallocate(x01,y01,z01)

else if (Method==2) then

do h=1,Nd2

do l=1,Nc2

x1=malla(l,1)
y1=malla(l,3)
z1=malla(l,5)
x2=malla(l,2)
y2=malla(l,4)
z2=malla(l,6)

call mbox(x02(h),y02(h),z02(h),x1,y1,z1,x2,y2,mi,md,fi,fd,theta,gm1)
call mbox(x02(h),y02(h),z02(h),x1,y1,z2,x2,y2,mi,md,fi,fd,theta,gm2)
Ad2(h,l)=gm1-gm2


end do !Nprism

!write(*,*) 'Ndat mag'
!write(*,*) h

end do !h datos

!WRITE(*,*)'Writting GM.dat'

!open(unit=16,file='GM.dat',status='replace',action='write')
!do i=1,Nd2
!write(16,*) (Ad2(i,j), j=1,Nc2)
!end do
!close(16)


WRITE(*,*) 'Writing magnetic data'

d2=matmul(Ad2,m02)

open(unit=39,file='Synthetic_magnetic.dat',status='replace',action='write')
do i=1,nd2
write(39,*) d2(i,1)
end do
close(39)



write(*,*)'Deallocate magnetic forward matrixes'
deallocate(m02)
deallocate(d2)
deallocate(Ad2)
deallocate(Ad2t)
deallocate(d21)
deallocate(malla)

deallocate(x02,y02,z02)

else
write(*,*)'Invalid option for method'
!NOTHING
end if

!!!!!!!!!!!FORWARD MODELLING FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!INVERSION MODELLING SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case(2)
write(*,*)'INV'

if (Method==1) then

WRITE(*,*) 'Constructing gravity sensitivity matrix'
call gbox(x01,y01,z01,Nd1,Nc1,malla,Ad1)
!end do
 Ad1=Ad1*1000 

write(*,*)'Calculating Spatial derivates using the grid information in density model'

call oper_derX(Nc1,Dx)
call oper_derY(Nc1,nx,ny,nz,Dy)
call oper_derZ(Nc1,nx,ny,nz,Dz)

!Weighthin operators
DO i= 2,NC1
Dx(i,:) = Dx(i,:) / ((malla(i,2)-malla(i-1,1))/2)
enddo

do i=1,Nz-1
do j=2,Ny-1
do k=1,Nx-1
Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)=Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)&
&/((malla(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),4)-malla(k+(j-2)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),3))/2)
end do
end do
end do

do i=2,Nz-1
do j=1,Ny-1
do k=1,Nx-1
Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)=Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)&
&/((malla(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),6)-malla(k+(j-1)*(Nx-1)+(i-2)*(Nx-1)*(Ny-1),5))/2)
end do
end do
end do

do i=1,Nc1
do j=1,Nc1
Dxt(i,j)=Dx(j,i)
Dyt(i,j)=Dy(j,i)
Dzt(i,j)=Dz(j,i)
end do
end do

Write(*,*)'Multiplicamos Dzt*Dz'
DxtDx=matmul(Dxt,Dx)
DytDy=matmul(Dyt,Dy)
DztDz=matmul(Dzt,Dz)


!!!!!!!!!!!!GRAVITY INVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)'Performing gravity inversion'

write(*,*)'Weighting matrixes'
do i= 1,Nd1
 Ad1(i,:) = Ad1(i,:) / vecCd1(i,1) ! Array math - possible alloc problem
 enddo
 
  do i = 1,Nd1
 d1(i,1) = d1(i,1) / vecCd1(i,1) ! Array math - possible alloc problem
 enddo

select case (Apriori)
case(1)      
write(*,*) 'Using density a priori model'
open(unit=40,file='m01pr_density.dat')
do i=1,Nc1
    read(40,*) m01pr(i,1)
end do
close(40)

open(unit=41,file='cov_m01pr.dat')
do i=1,Nc1
    read(41,*) vecCm01pr(i,1)
end do
close(41)

!Debo calcular la inversa de la diagonal dividiendo, lo hago de desde ahora para evitar la división entre cero después
 !Solo debemos tomar en cuenta que ahora vecCm01pr es en realidad el inverso cuadrado
 vecCm01pr=1./vecCm01pr**2
   call const_mat_diag(nc1,vecCm01pr,invCm01pr)

 write(*,*)'traspuesta Ad1:'
call traspuesta(nd1,nc1,Ad1,Ad1t)

write(*,*)'mult Ad1t*nvCd1*Ad1:'
ter11a=matmul(Ad1t,Ad1)
P1a=matmul(Ad1t,d1)

BB=ter11a+alpha1X*DxtDx+alpha1Y*DytDy+alpha1Z*DztDz+invCm01pr
PP=P1a+matmul(invCm01pr,m01pr)

case(2)
write(*,*)'You are NOT using a priori density model'

 write(*,*)'traspuesta Ad1:'
call traspuesta(nd1,nc1,Ad1,Ad1t)

write(*,*)'mult Ad1t*nvCd1*Ad1:'
ter11a=matmul(Ad1t,Ad1)
P1a=matmul(Ad1t,d1)

BB=ter11a+alpha1X*DxtDx+alpha1Y*DytDy+alpha1Z*DztDz
PP=P1a
CASE DEFAULT
WRITE(*,*) 'Invalid option'
END SELECT


call inverse(BB,invB,nc1)  !!Inversion with NO priori models
param=matmul(invB,PP)

do i=1,nc1
sol_m1(i,1)=param(i,1)
end do

do i= 1,Nd1
 Ad1(i,:) = Ad1(i,:) * vecCd1(i,1) ! Array math - possible alloc problem
enddo

  do i = 1,Nd1
 d1(i,1) = d1(i,1) * vecCd1(i,1) ! Array math - possible alloc problem
 enddo
 
d11=matmul(Ad1,sol_m1)
!ERROR RMS
error1=0


error1=0
do i=1,nd1
   error1=((d1(i,1)-d11(i,1))/(vecCd1(i,1)))**2+error1
end do
rms1=sqrt(error1/nd1)


open(unit=44,file='m11_separate.dat',status='replace',action='write')
do i=1,nc1
write(44,*) sol_m1(i,1)
end do
close(44)

open(unit=45,file='d11_separate.dat',status='replace',action='write')
do i=1,Nd1
write(45,*) d11(i,1)
end do
close(45)

open(unit=46,file='RMS_grav.dat',status='replace',action='write')
write(46,*) rms1
close(46)


write(*,*)'Deallocate gravity inversion matrixes'
deallocate(d1)
deallocate(Ad1,Cd1)
deallocate(ter11a,P1a,vecCd1)
deallocate(BB,PP)
deallocate(Ad1t)
deallocate(invB,param,sol_m1)
deallocate(d11)
deallocate(m01pr,vecCm01pr,Cm01pr,invCm01pr)
deallocate(malla)
deallocate(x01,y01,z01)
deallocate(Dx,Dy,Dz,Dxt,Dyt,Dzt)
deallocate(DxtDx,DytDy,DztDz)

!!!!!!!MAGNETIC INVERSION!!!!!!!
else if (Method==2) then

do h=1,Nd2

do l=1,Nc2

x1=malla(l,1)
y1=malla(l,3)
z1=malla(l,5)
x2=malla(l,2)
y2=malla(l,4)
z2=malla(l,6)

call mbox(x02(h),y02(h),z02(h),x1,y1,z1,x2,y2,mi,md,fi,fd,theta,gm1)
call mbox(x02(h),y02(h),z02(h),x1,y1,z2,x2,y2,mi,md,fi,fd,theta,gm2)
Ad2(h,l)=gm1-gm2


end do !Nprism

!write(*,*) 'Ndat mag'
!write(*,*) h

end do !h datos



write(*,*)'Calculating Spatial derivates using the grid information in magnetization model'

call oper_derX(Nc2,Dx)
call oper_derY(Nc2,nx,ny,nz,Dy)
call oper_derZ(Nc2,nx,ny,nz,Dz)

!Weighthin operators
DO i= 2,NC2
Dx(i,:) = Dx(i,:) / ((malla(i,2)-malla(i-1,1))/2)
enddo

write(*,*)'FLAG3.5 Dx(3,3)'
write(*,*) Dx(3,3)

do i=1,Nz-1
do j=2,Ny-1
do k=1,Nx-1
Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)=Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)&
&/((malla(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),4)-malla(k+(j-2)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),3))/2)
end do
end do
end do

do i=2,Nz-1
do j=1,Ny-1
do k=1,Nx-1
Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)=Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)&
&/((malla(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),6)-malla(k+(j-1)*(Nx-1)+(i-2)*(Nx-1)*(Ny-1),5))/2)
end do
end do
end do

do i=1,Nc2
do j=1,Nc2
Dxt(i,j)=Dx(j,i)
Dyt(i,j)=Dy(j,i)
Dzt(i,j)=Dz(j,i)
end do
end do

Write(*,*)'Multiplicamos Dzt*Dz'
DxtDx=matmul(Dxt,Dx)
DytDy=matmul(Dyt,Dy)
DztDz=matmul(Dzt,Dz)


write(*,*)'FLAG3.55 DxtDx(3,3)'
write(*,*) DxtDx(3,3)
!!!!!!!!!!!!MAGNETIC INVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)'Performing magnetic inversion'
write(*,*)'Weighting matrixes'

do i= 1,Nd2
 Ad2(i,:) = Ad2(i,:) / vecCd2(i,1) ! Array math - possible alloc problem
enddo

write(*,*)'FLAG3.6 Ad2(3,3)'
write(*,*) Ad2(3,3)

  do i = 1,Nd2
 d2(i,1) = d2(i,1) / vecCd2(i,1) ! Array math - possible alloc problem
 enddo


select case(Apriori)
case(1)      
write(*,*) 'Using magnetization a priori model'
open(unit=42,file='m02pr_magnetization.dat')
do i=1,Nc2
    read(42,*) m02pr(i,1)
end do
close(42)


write(*,*)'FLAG3.8 m02pr(3,1)'
write(*,*) m02pr(3,1)

open(unit=43,file='cov_m02pr.dat')
do i=1,Nc2
    read(43,*) vecCm02pr(i,1)
end do
close(43)


!Debo calcular la inversa de la diagonal dividiendo, lo hago de desde ahora para evitar la división entre cero después
 !Solo debemos tomar en cuenta que ahora vecCm01pr es en realidad el inverso cuadrado
 vecCm02pr=1./vecCm02pr**2

  call const_mat_diag(nc2,vecCm02pr,invCm02pr)
  
   write(*,*)'traspuesta Ad2:'
call traspuesta(nd2,nc2,Ad2,Ad2t)

write(*,*)'mult Ad2t*nvCd2*Ad2:'
ter22a=matmul(Ad2t,Ad2)
P2a=matmul(Ad2t,d2)

BB=ter22a+alpha2X*DxtDx+alpha2Y*DytDy+alpha2Z*DztDz+invCm02pr
PP=P2a+matmul(invCm02pr,m02pr)
write(*,*)'FLAG3 BB(3,1)'
write(*,*) BB(3,3)
write(*,*)'FLAG4 PP(3,1)'
write(*,*) PP(3,1)

case(2)
write(*,*)'You are NOT using a priori magnetization model'
   write(*,*)'traspuesta Ad2:'
call traspuesta(nd2,nc2,Ad2,Ad2t)

write(*,*)'mult Ad2t*nvCd2*Ad2:'
ter22a=matmul(Ad2t,Ad2)
P2a=matmul(Ad2t,d2)

BB=ter22a+alpha2X*DxtDx+alpha2Y*DytDy+alpha2Z*DztDz
PP=P2a

CASE DEFAULT
WRITE(*,*) 'Invalid option'
END SELECT

call inverse(BB,invB,nc2)  !!Inversion with NO priori models
param=matmul(invB,PP)


do i=1,nc2
sol_m2(i,1)=param(i,1)
end do

do i= 1,Nd2
 Ad2(i,:) = Ad2(i,:) * vecCd2(i,1) ! Array math - possible alloc problem
enddo

  do i = 1,Nd2
 d2(i,1) = d2(i,1) * vecCd2(i,1) ! Array math - possible alloc problem
 enddo
 
d21=matmul(Ad2,sol_m2)
!ERROR RMS
error2=0


error2=0
do i=1,nd2
   error2=((d2(i,1)-d21(i,1))/(vecCd2(i,1)))**2+error2
end do
rms2=sqrt(error2/nd2)


open(unit=47,file='m21_separate.dat',status='replace',action='write')
do i=1,nc2
write(47,*) sol_m2(i,1)
end do
close(47)

open(unit=48,file='d21_separate.dat',status='replace',action='write')
do i=1,Nd2
write(48,*) d21(i,1)
end do
close(48)

open(unit=49,file='RMS_mag.dat',status='replace',action='write')
write(49,*) rms2
close(49)


write(*,*)'Deallocate magnetic inversion matrixes'
deallocate(d2)
deallocate(Ad2,Cd2)
deallocate(DxtDx,DytDy,DztDz)
deallocate(ter22a,P2a,vecCd2)
deallocate(BB,PP)
deallocate(Ad2t)
deallocate(invB,param,sol_m2)
deallocate(d21)
deallocate(m02pr,vecCm02pr,Cm02pr,invCm02pr)
deallocate(malla)
deallocate(x02,y02,z02)
deallocate(Dx,Dy,Dz,Dxt,Dyt,Dzt)

else
write(*,*)'Invalid option for method'
!NOTHING
end if






!!!!!!!!!!!!!!!!!!!!!JOINT INVERSION MODULE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ENTRAN COORDENADAS UTM, PERO INTERNAMENTE SE LEEN AL REVÉS (X AL NORTE, Y AL ESTE)


write(*,*)'JOINT'
case(3) !!! 	Case 3 to action

select case(Apriori)

case(1)        !!!!!OPTION 1: BOTH A PRIORI
write(*,*) 'Using both a priori models'
open(unit=50,file='m01pr_density.dat')
do i=1,Nc1
    read(50,*) m01pr(i,1)
end do
close(50)

open(unit=51,file='cov_m01pr.dat')
do i=1,Nc1
    read(51,*) vecCm01pr(i,1)
end do
close(51)

open(unit=52,file='m02pr_magnetization.dat')
do i=1,Nc2
    read(52,*) m02pr(i,1)
end do
close(52)


open(unit=53,file='cov_m02pr.dat')
do i=1,Nc2
    read(53,*) vecCm02pr(i,1)
end do
close(53)


!Debo calcular la inversa de la diagonal dividiendo, lo hago de desde ahora para evitar la división entre cero después
 !Solo debemos tomar en cuenta que ahora vecCm01pr es en realidad el inverso cuadrado
 vecCm01pr=1./vecCm01pr**2
  vecCm02pr=1./vecCm02pr**2


  call const_mat_diag(nc1,vecCm01pr,invCm01pr)
    call const_mat_diag(nc2,vecCm02pr,invCm02pr)
   
   write(*,*)'invCm01pr (3,3)'
write(*,*) invCm01pr(3,3)

 write(*,*)'invCm02pr (3,3)'
write(*,*) invCm02pr(3,3)
    
    case(2)        !!!!!!!!OPTION 2: NO
WRITE(*,*) 'Not using a priori models'
!else if(Apriori==2)       then    !!!!OPTION 2: NO
!NOTHING
  case (3)      !!!!!!!!!!!!!!!!!!! !OPION 3: ONLY GRAV A PRIORI
  WRITE(*,*) 'Using density a priori model'
open(unit=54,file='m01pr_density.dat')
do i=1,Nc1
    read(54,*) m01pr(i,1)
end do
close(54)

open(unit=55,file='cov_m01pr.dat')
do i=1,Nc1
    read(55,*) vecCm01pr(i,1)
end do
close(55)


 !Debo calcular la inversa de la diagonal dividiendo, lo hago de desde ahora para evitar la división entre cero después
 !Solo debemos tomar en cuenta que ahora vecCm01pr es en realidad el inverso cuadrado
 vecCm01pr=1./vecCm01pr**2

 call const_mat_diag(nc1,vecCm01pr,invCm01pr)
 
 case(4)      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OPION 4: ONLY MAG A PRIORI
WRITE(*,*) 'Using magnetization a priori model'

open(unit=56,file='m02pr_magnetization.dat')
do i=1,Nc2
    read(56,*) m02pr(i,1)
end do
close(56)


open(unit=57,file='cov_m02pr.dat')
do i=1,Nc2
    read(57,*) vecCm02pr(i,1)
end do
close(57)

  !Debo calcular la inversa de la diagonal dividiendo, lo hago de desde ahora para evitar la división entre cero después
 !Solo debemos tomar en cuenta que ahora vecCm01pr es en realidad el inverso cuadrado
  vecCm02pr=1./vecCm02pr**2

 call const_mat_diag(nc2,vecCm02pr,invCm02pr)

CASE DEFAULT
      WRITE(*,*) 'Invalid option'
END SELECT 

!!!!!!!!SENSITIVITY
WRITE(*,*) 'Constructing gravity sensitivity matrix'
call gbox(x01,y01,z01,Nd1,Nc1,malla,Ad1)
!end do
 Ad1=Ad1*1000 


WRITE(*,*) 'Constructing magnetic sensitivity matrix'

do h=1,Nd2

do l=1,Nc2

x1=malla(l,1)
y1=malla(l,3)
z1=malla(l,5)
x2=malla(l,2)
y2=malla(l,4)
z2=malla(l,6)

call mbox(x02(h),y02(h),z02(h),x1,y1,z1,x2,y2,mi,md,fi,fd,theta,gm1)
call mbox(x02(h),y02(h),z02(h),x1,y1,z2,x2,y2,mi,md,fi,fd,theta,gm2)
Ad2(h,l)=gm1-gm2


end do !Nprism

!write(*,*) 'Ndat mag'
!write(*,*) h

end do !h datos




!!!!SPATIAL DERIVATES
write(*,*)'Calculating Spatial derivates using the grid information in models'

call oper_derX(Nc1,Dx)
call oper_derY(Nc1,nx,ny,nz,Dy)
call oper_derZ(Nc1,nx,ny,nz,Dz)

!Weighthin operators
DO i= 2,NC1
Dx(i,:) = Dx(i,:) / ((malla(i,2)-malla(i-1,1))/2)
enddo

do i=1,Nz-1
do j=2,Ny-1
do k=1,Nx-1
Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)=Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)&
&/((malla(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),4)-malla(k+(j-2)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),3))/2)
end do
end do
end do

do i=2,Nz-1
do j=1,Ny-1
do k=1,Nx-1
Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)=Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),:)&
&/((malla(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),6)-malla(k+(j-1)*(Nx-1)+(i-2)*(Nx-1)*(Ny-1),5))/2)
end do
end do
end do

do i=1,Nc1
do j=1,Nc1
Dxt(i,j)=Dx(j,i)
Dyt(i,j)=Dy(j,i)
Dzt(i,j)=Dz(j,i)
end do
end do

Write(*,*)'Multiplicamos Dzt*Dz'
DxtDx=matmul(Dxt,Dx)
DytDy=matmul(Dyt,Dy)
DztDz=matmul(Dzt,Dz)

!!!!!!!!!!!!!! 4 RELACION DE CORRESPONDENCIA !!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nc
    g_obs(i,1)=-1
end do

 !!PESO TODAS LAS MATRICES DE SENSIBILIDAD   Y DATOS


!Modelo inicial de coeficientes
!do i=1,ncoef
 a0(1,1)=1.5
  a0(2,1)=-1.5
!a0(3,1)=1.5
!end do


!!! SECCION DEL CICLO
 rms=10
 open(unit=101,file='Xout_sal.dat',status='replace',action='write')
 open(unit=500,file='rms1_2_5.dat',status='replace',action='write')
open(unit=501,file='Lagrange_mult.dat',status='replace',action='write')
open(unit=502,file='delta_modelo.dat',status='replace',action='write')
open(unit=503,file='delta_m1.dat',status='replace',action='write')
open(unit=504,file='delta_m2.dat',status='replace',action='write')
open(unit=505,file='delta_a.dat',status='replace',action='write')

conta=0

do i=1,nc1
m01(i,1)=cos((0.01+(2*pi-0.01)*rand()))
end do

do i=1,nc2
m02(i,1)=sin((0.01+(2*pi-0.01)*rand()))
end do

open(unit=888,file='m01.dat',status='replace',action='write')
open(unit=988,file='m02.dat',status='replace',action='write')
do i=1,nc1
write(888,*) m01(i,1)
   end do

do i=1,nc2
write(988,*) m02(i,1)
   end do

 close(888)
 close(988)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ITERA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do while(abs(rms)>=1.5)

    conta=conta+1

 do i= 1,Nd1
 Ad1(i,:) = Ad1(i,:) / vecCd1(i,1) ! Array math - possible alloc problem
 enddo

 do i = 1,Nd1
 d1(i,1) = d1(i,1) / vecCd1(i,1) ! Array math - possible alloc problem
 enddo


 do i= 1,Nd2
 Ad2(i,:) = Ad2(i,:) / vecCd2(i,1) ! Array math - possible alloc problem
 enddo

 do i = 1,Nd2
 d2(i,1) = d2(i,1) / vecCd2(i,1) ! Array math - possible alloc problem
 enddo

  do i = 1,nc
 g_obs(i,1) = g_obs(i,1) / vecCg(i,1) ! Array math - possible alloc problem
 enddo
!  do i = 1,nc1
! m01pr(i,1) = m01pr(i,1) / vecCm01pr(i,1) ! Array math - possible alloc problem
 !enddo
 

!SECCION PARA CONSTRUIR LAS MATRICES DE SENSIBILIDAD DE LA RELACION
write(*,*)'Voy a calcular Agm1 y Agm2'
!call Agm1_sin_a00_pol2d(nc,m01,m02,a0,Dfm1)
!call Agm2_sin_a00_pol2d(nc,m01,m02,a0,Dfm2)


do i=1,nc
Dfm1(i,1)=a0(1,1)
    end do

do i=1,nc
Dfm2(i,1)=a0(2,1)
    end do


call const_mat_diag(nc,Dfm1,Ag1)
call const_mat_diag(nc,Dfm2,Ag2)


! open(unit=298,file='Ag1.dat',status='replace',action='write')
!do i=1,nc
!write(298,*) (Ag1(i,j), j=1,nc)
!   end do
! close(298)

!open(unit=299,file='Ag2.dat',status='replace',action='write')
!do i=1,nc
!write(299,*) (Ag2(i,j), j=1,nc)
!   end do
! close(299)


write(*,*)'Voy a calcular Aga'
!call Aga_sin_a00_pol2d(nc,m01,m02,Aa)
do i=1,nc
Aa(i,1)=m01(i,1)
Aa(i,2)=m02(i,1)
!Aa(i,3)=m02(i,1)**2
end do

open(unit=300,file='Aga.dat',status='replace',action='write')
do i=1,nc
write(300,*) (Aa(i,j), j=1,ncoef)
   end do
 close(300)
 
 !PESO LAS MATRICES DE SENSIBILIDAD RESPECTO A g
 
  do i= 1,nc
 Ag1(i,:) = Ag1(i,:) / vecCg(i,1) ! Array math - possible alloc problem
 enddo

  do i= 1,nc
 Ag2(i,:) = Ag2(i,:) / vecCg(i,1) ! Array math - possible alloc problem
 enddo

 do i= 1,nc
 Aa(i,:) = Aa(i,:) / vecCg(i,1) ! Array math - possible alloc problem
 enddo

 !Respuesta del modelo inicial
g0=matmul(Aa,a0)
d10=matmul(Ad1,m01)
d20=matmul(Ad2,m02)
 !Traspongo las matrices
 
 write(*,*)'traspuesta Ad1:'
call traspuesta(nd1,nc1,Ad1,Ad1t)

write(*,*)'mult Ad1t*nvCd1*Ad1:'
ter11a=matmul(Ad1t,Ad1)

write(*,*)'ter11a 2,2'
WRITE(*,*) ter11a(2,2)

!Construccion B
!Momentaneamente se definira como identidad al termino Ad2
write(*,*)'Voy a construir B'
write(*,*)'Voy a constuir B11 y B22'
call traspuesta(nc,nc1,Ag1,Ag1t)
call traspuesta(nc,nc2,Ag2,Ag2t)
call traspuesta(nc,ncoef,Aa,Aat)

!ter11c1=matmul(Ag1t,invCg) !!SE REPITE EN B11, B12 Y B13
ter11c=matmul(Ag1t,Ag1)

write(*,*)'ter11c(4,1)'
write(*,*) ter11c(4,1)


 write(*,*)'traspuesta Ad2:'
call traspuesta(nd2,nc2,Ad2,Ad2t)

 write(*,*)'Ad2t(1,1)'
write(*,*) Ad2t(1,1)

write(*,*)'mult Ad2t*nvCd2*Ad2:'
ter22a=matmul(Ad2t,Ad2)

write(*,*)'ter22a 1,1'
WRITE(*,*) ter22a(1,1)

!ter22c1=matmul(Ag2t,invCg)  !!SE REPITE EN B21, B22 Y B23
ter22c=matmul(Ag2t,Ag2)

!!B22=ter22a+alpha2*PTPf+ter22c!+invCm02
select case (Apriori)
case (1)
B11=ter11a+alpha1X*DxtDx+alpha1Y*DytDy+alpha1Z*DztDz+ter11c+invCm01pr
B22=ter22a+alpha2X*DxtDx+alpha2Y*DytDy+alpha2Z*DztDz+ter22c+invCm02pr
case (2)
B11=ter11a+alpha1X*DxtDx+alpha1Y*DytDy+alpha1Z*DztDz+ter11c
B22=ter22a+alpha2X*DxtDx+alpha2Y*DytDy+alpha2Z*DztDz+ter22c
case (3)
B11=ter11a+alpha1X*DxtDx+alpha1Y*DytDy+alpha1Z*DztDz+ter11c+invCm01pr
B22=ter22a+alpha2X*DxtDx+alpha2Y*DytDy+alpha2Z*DztDz+ter22c
case (4)
B11=ter11a+alpha1X*DxtDx+alpha1Y*DytDy+alpha1Z*DztDz+ter11c
B22=ter22a+alpha2X*DxtDx+alpha2Y*DytDy+alpha2Z*DztDz+ter22c+invCm02pr
!B11=gtg+alpha1*DxtDx+alpha2*DytDy+alpha3*DztDz
CASE DEFAULT
WRITE(*,*) 'Invalid option'
END SELECT

!WJTWJ+LagrangeMu*PTP+ter22c
write(*,*)'B11 2,2'
WRITE(*,*) B11(2,2)

write(*,*)'B22 2,2'
WRITE(*,*) B22(2,2)
write(*,*)'Voy a constuir B33'

!ter33a1=matmul(Aat,invCg)  !! SE REPITE EN B31, B32, B33
B33=matmul(Aat,Aa)

write(*,*)'B33 2,2'
WRITE(*,*) B33(2,2)

write(*,*)'Voy a constuir B12, B21, B13, B31, B23, B32'


B12=matmul(Ag1t,Ag2)

write(*,*)'B12 2,2'
WRITE(*,*) B12(2,2)

B21=matmul(Ag2t,Ag1)

write(*,*)'B21 2,2'
WRITE(*,*) B21(2,2)

B13=matmul(Ag1t,Aa)

write(*,*)'B13 2,2'
WRITE(*,*) B13(2,2)

B31=matmul(Aat,Ag1)

write(*,*)'B31 2,2'
WRITE(*,*) B31(2,2)

B23=matmul(Ag2t,Aa)

write(*,*)'B23 2,2'
WRITE(*,*) B23(2,2)

B32=matmul(Aat,Ag2)

write(*,*)'B32 2,2'
WRITE(*,*) B32(2,2)

call Const_B(nc1,nc2,ncoef,B11,B12,B13,B21,B22,B23,B31,B32,B33,BB)

write(*,*)'B(3,4)'
write(*,*) BB(3,4)

!open(unit=199,file='BB.dat',status='replace',action='write')
!do i=1,nc1+nc2+ncoef
!write(199,*) (BB(i,j), j=1,nc1+nc2+ncoef)
!   end do
! close(199)

write(*,*)'Voy a constuir P'

Ad1m01=matmul(Ad1,m01)
Ad2m02=matmul(Ad2,m02)
Ag1m01=matmul(Ag1,m01)
Ag2m02=matmul(Ag2,m02)
Aaa0=matmul(Aa,a0)

dat2c=g_obs-g0+Ag1m01+Ag2m02+Aaa0    !CUIDADO

!!!!!!dat3=g_obs-g0+Ag1m01+Ag2m02+Aaa0

P1a=matmul(Ad1t,d1)

write(*,*)'P1a(4,1)'
write(*,*) P1a(4,1)

P1c=matmul(Ag1t,dat2c)


P2c=matmul(Ag2t,dat2c)
P2=matmul(Ad2t,d2)


select case(Apriori)
case(1)
P1=P1a+P1c+matmul(invCm01pr,m01pr)
P2=P2+P2c+matmul(invCm02pr,m02pr)
case (2)
P1=P1a+P1c!+matmul(invCm01pr,m01pr)
P2=P2+P2c!+matmul(invCm02pr,m02pr)
case (3)
P1=P1a+P1c+matmul(invCm01pr,m01pr)
P2=P2+P2c!+matmul(invCm02pr,m02pr)
case (4)
P1=P1a+P1c!+matmul(invCm01pr,m01pr)
P2=P2+P2c+matmul(invCm02pr,m02pr)
CASE DEFAULT
WRITE(*,*) 'Invalid option'
END SELECT

write(*,*)'     P1 2,2'
WRITE(*,*) P1(2,2)
write(*,*)'P2 2,2'
WRITE(*,*) P2(2,2)

P3=matmul(Aat,dat2c)
write(*,*)'P3 2,2'
WRITE(*,*) P3(2,2)
call Const_P(nc1,nc2,ncoef,P1,P2,P3,PP)

write(*,*)'PP(4,1)'
write(*,*) PP(4,1)

write(*,*)'Voy a invertir B'
write(*,*)'BB(2664,2664)'
write(*,*) BB(2664,2664)

write(*,*) 'Resuelve el sistema de ecuaciones'

!!!!AQUI OEGAR PARDISO


!!AQUI TERMINA PARDISO
call inverse(BB,invB,nc1+nc2+ncoef)
param=matmul(invB,PP)

!CALCULO EL PASO DEL MODELO GENERAL, Y UNO POR UNO     !!!!CUIDADO, NO HAY PASOS PORQUE NO ES ITERATIVO
delt_m1rms=0
delt_m2rms=0
delt_arms=0

do i=1,nc1
 delt_m1rms=(param(i,1)-m01(i,1))**2+delt_m1rms
end do

delt_m1rms=sqrt(delt_m1rms/nc1)

do i=1,nc2
 delt_m2rms=(param(i+nc1,1)-m02(i,1))**2+delt_m2rms
end do

delt_m2rms=sqrt(delt_m2rms/nc2)

do i=1,ncoef
 delt_arms=(param(i+nc1+nc2,1)-a0(i,1))**2+delt_arms
end do

delt_arms=sqrt(delt_arms/ncoef)



delt_mrms=delt_m1rms+delt_m2rms+delt_arms

do i=1,nc1
sol_m1(i,1)=param(i,1)
end do

do i=1,nc2
sol_m2(i,1)=param(i+nc1,1)
end do

do i=1,ncoef
sol_a(i,1)=param(i+nc1+nc2,1)
end do

!RESPUESTA DEL NUEVO MODELO,(sol)

 do i= 1,Nd1
 Ad1(i,:) = Ad1(i,:) * vecCd1(i,1) ! Array math - possible alloc problem
 enddo

 do i = 1,Nd1
 d1(i,1) = d1(i,1) * vecCd1(i,1) ! Array math - possible alloc problem
 enddo


 do i= 1,Nd2
 Ad2(i,:) = Ad2(i,:) * vecCd2(i,1) ! Array math - possible alloc problem
 enddo

 do i = 1,Nd2
 d2(i,1) = d2(i,1) * vecCd2(i,1) ! Array math - possible alloc problem
 enddo
 

 do i = 1,nc
 g_obs(i,1) = g_obs(i,1) * vecCg(i,1) ! Array math - possible alloc problem
 enddo
 
   do i= 1,nc
 Ag1(i,:) = Ag1(i,:) * vecCg(i,1) ! Array math - possible alloc problem
 enddo

  do i= 1,nc
 Ag2(i,:) = Ag2(i,:) * vecCg(i,1) ! Array math - possible alloc problem
 enddo

 do i= 1,nc
 Aa(i,:) = Aa(i,:) * vecCg(i,1) ! Array math - possible alloc problem
 enddo

!RESPUESTA DEL NUEVO MODELO

g1=matmul(Aa,sol_a)
d11=matmul(Ad1,sol_m1)
d21=matmul(Ad2,sol_m2)


!!!ERRORES
!ERROR RMS
error1=0
error4=0
error5=0
error6=0
error7=0
error8=0

do i=1,nc
 error1=((g_obs(i,1)-g1(i,1))/(vecCg(i,1)))**2+error1
end do

rms1=sqrt(error1/nc)


error2=0
do i=1,nd1
   error2=((d1(i,1)-d11(i,1))/(vecCd1(i,1)))**2+error2
end do
rms2=sqrt(error2/nd1)

error3=0
do i=1,nd2
   error3=((d2(i,1)-d21(i,1))/(vecCd2(i,1)))**2+error3
end do
rms3=sqrt(error3/nd2)

rms=rms1+rms2+rms3!+rms7+rms8


write(*,*)'RMS 1 g, RMS2 d1, RMS3 d2'
WRITE(*,*) rms1,rms2,rms3!,rms7,rms8


write(500,*) rms1,rms2,rms3
write(502,*) delt_modelo
write(503,*) delt_m1rms
write(504,*) delt_m2rms
write(505,*) delt_arms


m01=sol_m1
m02=sol_m2
a0=sol_a

!DATOS PREDICHOS
open(unit=14,file='d11_'//trim(str(conta))//'.dat',status='replace',action='write')
do i=1,nd1
write(14,*) d11(i,1)
    end do
close(14)

open(unit=15,file='d21_'//trim(str(conta))//'.dat',status='replace',action='write')
do i=1,nd2
write(15,*) d21(i,1)
    end do
close(15)

open(unit=16,file='g_pred_'//trim(str(conta))//'.dat',status='replace',action='write')
do i=1,nc
write(16,*) g1(i,1)
    end do
close(16)

open(unit=17,file='Sol_a__VECTOR_'//trim(str(conta))//'.dat',status='replace',action='write')
do i=1,ncoef
write(17,*) sol_a(i,1)
    end do
close(17)

open(unit=18,file='m11_'//trim(str(conta))//'.dat',status='replace',action='write')
open(unit=19,file='m21_'//trim(str(conta))//'.dat',status='replace',action='write')
do i=1,nc1
write(18,*) sol_m1(i,1)
   end do

do i=1,nc2
write(19,*) sol_m2(i,1)
   end do

 close(18)
 close(19)

   IF((rms1)<=2) THEN
         EXIT
      END IF

   !IF((delt_mrms)<=0.0005) THEN
    !     EXIT
     ! END IF
      
  end do


close(500)
close(501)
close(502)
close(503)
close(504)
close(505)

!DATOS PREDICHOS

write(*,*)'Ag1 3,3'
WRITE(*,*) Ag1(3,3)

write(*,*)'Ag2 3,3'
WRITE(*,*) Ag2(3,3)

write(*,*)'Aa 4,3'
WRITE(*,*) Aa(4,3)

write(*,*)'B 3,3'
WRITE(*,*) BB(3,3)

write(*,*)'P 5,1'
WRITE(*,*) PP(5,1)

write(*,*)'sol_m1 5,1'
WRITE(*,*) sol_m1(5,1)

write(*,*)'sol_m2 5,1'
WRITE(*,*) sol_m2(5,1)

write(*,*)'sol_a 5,1'
WRITE(*,*) sol_a(5,1)

write(*,*)'RMS'
WRITE(*,*) rms

call cpu_time(tiempo) 
write (*, *) tiempo


! libera espacio de arreglos dinamicos 
deallocate(x01,y01,z01,x02,y02,z02)
deallocate(Dz,Dzt,Dy,Dyt,Dx,Dxt)
deallocate(d1,d2,g_obs)
deallocate(Ad1,Ad2,Cd1,Cd2)
deallocate(a0,m01,m02)
deallocate(DxtDx,DytDy,DztDz)
deallocate(vecCg)
deallocate(vecCd1,vecCd2)
deallocate(Dfm1,Dfm2,Cg)
deallocate(Ag1,Ag2,Aa)
deallocate(g0,d10,d20,BB,PP)
deallocate(Ad1t,Ag1t)
deallocate(Ad2t,Ag2t,Aat)
deallocate(ter11a,ter11c)
deallocate(ter22a1,ter22a,ter22c)
deallocate(B11,B22,B33,B12,B13,B21,B23)
deallocate(B31,B32)
deallocate(Ad1m01,Ag1m01,Ag2m02,Aaa0,Ad2m02)
deallocate(dat1a,dat1c,dat2a,dat2c,dat3)
deallocate(P1a,P1c,P1)
deallocate(P2a,P2c,P2)
deallocate(P3)
deallocate(invB,param,sol_m1,sol_m2,sol_a)
deallocate(g1,d11,d21)

CASE DEFAULT
WRITE(*,*) 'Invalid option Action' !!!!!!ENDS ACTION 4
END SELECT


END PROGRAM



















!!!SECCION DE SUBRUTINAS
Subroutine Agm1_sin_a00_pol2d(nc,m1,m2,coef,Dfm1)
!Esta funcion calcula el jacobiano respecto a m1, tomando en cuenta que
!el coeficiente a00 no se toma en cuenta para evitar soluci?n trivial
!La entrada son los dos modelos, m1 y m2, y los coeficientes propuestos coef del polinomio
implicit none
integer, intent(in)::nc
double precision , intent(in)::m1(nc,1),m2(nc,1),coef(8,1)
double precision ,intent(out)::Dfm1(nc,1)
double precision :: Agm1(nc,8)
integer::i,j

do i=1,nc
Agm1(i,1)=0
Agm1(i,2)=0
Agm1(i,3)=1*coef(3,1)
Agm1(i,4)=1*coef(4,1)*(m2(i,1)**(1))
Agm1(i,5)=1*coef(5,1)*(m2(i,1)**(2))
Agm1(i,6)=2*coef(6,1)*(m1(i,1)**(2-1))
Agm1(i,7)=2*coef(7,1)*(m1(i,1)**(2-1))*(m2(i,1)**(1))
Agm1(i,8)=2*coef(8,1)*(m1(i,1)**(2-1))*(m2(i,1)**(2))
end do

do i=1,nc
Dfm1(i,1)=Agm1(i,1)+Agm1(i,2)+Agm1(i,3)+Agm1(i,4)+Agm1(i,5)+Agm1(i,6)+Agm1(i,7)+Agm1(i,8)
end do

return

end subroutine Agm1_sin_a00_pol2d


Subroutine Agm2_sin_a00_pol2d(nc,m1,m2,coef,Dfm2)
!Esta funcion calcula el jacobiano respecto a m1, tomando en cuenta que
!el coeficiente a00 no se toma en cuenta para evitar soluci?n trivial
!La entrada son los dos modelos, m1 y m2, y los coeficientes propuestos coef del polinomio
implicit none
integer, intent(in)::nc
double precision, intent(in)::m1(nc,1),m2(nc,1),coef(8,1)
double precision ,intent(out)::Dfm2(nc,1)
double precision :: Agm2(nc,8)
integer::i,j

do i=1,nc


Agm2(i,1)=1*coef(1,1)
Agm2(i,2)=2*coef(2,1)*(m2(i,1)**(2-1))
Agm2(i,3)=0
Agm2(i,4)=1*coef(4,1)*(m1(i,1)**(1))
Agm2(i,5)=2*coef(5,1)*(m1(i,1)**(1))*(m2(i,1)**(2-1))
Agm2(i,6)=0
Agm2(i,7)=1*coef(7,1)*(m1(i,1)**(2))
Agm2(i,8)=2*coef(8,1)*(m1(i,1)**(2))*(m2(i,1)**(2-1))
end do

do i=1,nc
Dfm2(i,1)=Agm2(i,1)+Agm2(i,2)+Agm2(i,3)+Agm2(i,4)+Agm2(i,5)+Agm2(i,6)+Agm2(i,7)+Agm2(i,8)
end do

return

end subroutine Agm2_sin_a00_pol2d

Subroutine Aga_sin_a00_pol2d(nc,m1,m2,Aga)
!Esta funcion calcula el jacobiano respecto a coef, tomando en cuenta que
!el coeficiente a00 no se toma en cuenta para evitar soluci?n trivial
implicit none
integer, intent(in)::nc
double precision, intent(in)::m1(nc,1),m2(nc,1)
double precision ,intent(out)::Aga(nc,8)
integer::i,j

do i=1,nc
Aga(i,1)=(m2(i,1))
Aga(i,2)=(m2(i,1)**(2))
Aga(i,3)=(m1(i,1))
Aga(i,4)=(m1(i,1))*(m2(i,1))
Aga(i,5)=(m1(i,1))*(m2(i,1)**(2))
Aga(i,6)=(m1(i,1)**(2))
Aga(i,7)=(m1(i,1)**(2))*(m2(i,1))
Aga(i,8)=(m1(i,1)**(2))*(m2(i,1)**(2))
end do

return

end subroutine Aga_sin_a00_pol2d

!Esta funcion construye B, de la ecuacion Bm=P
Subroutine Const_B(nc1,nc2,ncoef,B11,B12,B13,B21,B22,B23,B31,B32,B33,B)
!Esta subrutina construye B, de la formulacion Bm=P
implicit none
integer, intent(in)::nc1,nc2,ncoef
double precision, intent(in)::B11(nc1,nc1),B22(nc2,nc2),B33(ncoef,ncoef),B12(nc1,nc2),B13(nc1,ncoef)
double precision, intent(in)::B21(nc2,nc1),B23(nc2,ncoef),B31(ncoef,nc2),B32(ncoef,nc2)
double precision ,intent(out)::B(nc1+nc2+ncoef,nc1+nc2+ncoef)
integer::i,j

do i=1,nc1
    do j=1,nc1
B(i,j)=B11(i,j)
!B(i,j+nc1)=B12(i,j)
!B(i+nc1,j)=B21(i,j)
!B(i+nc1,j+nc1)=B22(i,j)
    end do
    end do
do i=1,nc2
    do j=1,nc2
!B(i,j)=B11(i,j)
!B(i,j+nc1)=B12(i,j)
!B(i+nc1,j)=B21(i,j)
B(i+nc1,j+nc1)=B22(i,j)
    end do
    end do
do i=1,nc1
    do j=1,nc2
!B(i,j)=B11(i,j)
B(i,j+nc1)=B12(i,j)
!B(i+nc1,j)=B21(i,j)
!B(i+nc1,j+nc1)=B22(i,j)
    end do
    end do

do i=1,nc2
    do j=1,nc1
!B(i,j)=B11(i,j)
!B(i,j+nc1)=B12(i,j)
B(i+nc1,j)=B21(i,j)
!B(i+nc1,j+nc1)=B22(i,j)
    end do
    end do

do i=1,nc1
    do j=1,ncoef
B(i,j+nc1+nc2)=B13(i,j)
!B(i+nc1,j+nc1+nc2)=B23(i,j)
    end do
end do

do i=1,nc2
    do j=1,ncoef
!B(i,j+nc1+nc2)=B13(i,j)
B(i+nc1,j+nc1+nc2)=B23(i,j)
    end do
end do

do i=1,ncoef
    do j=1,nc1
B(i+nc1+nc2,j)=B31(i,j)
!B(i+2*nc,j+nc)=B32(i,j)
    end do
end do

do i=1,ncoef
    do j=1,nc2
!B(i+2*nc,j)=B31(i,j)
B(i+nc1+nc2,j+nc1)=B32(i,j)
    end do
end do

do i=1,ncoef
    do j=1,ncoef
B(i+nc1+nc2,j+nc1+nc2)=B33(i,j)
end do
end do

return
end subroutine Const_B

subroutine Const_P(nc1,nc2,ncoef,P1,P2,P3,P)
!Esta subrutina construye P, de la formulacion Bm=P
implicit none
integer, intent(in)::nc1,nc2,ncoef
double precision, intent(in)::P1(nc1,1),P2(nc2,1),P3(ncoef,1)
double precision ,intent(out)::P(nc1+nc2+ncoef,1)
integer::i,j


do i=1,nc1
P(i,1)=P1(i,1)
!P(i+nc,1)=P2(i,1)
    end do

do i=1,nc2
!P(i,1)=P1(i,1)
P(i+nc1,1)=P2(i,1)
    end do

do i=1,ncoef
    P(i+nc1+nc2,1)=P3(i,1)
    end do

return

end subroutine Const_P


 subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do

end subroutine inverse

subroutine traspuesta(n1,n2,A,B)
!Calcula la traspuesta B de una matriz A de n1xn2
implicit none
integer, intent(in)::n1,n2
double precision, intent(in)::A(n1,n2)
double precision ,intent(out)::B(n2,n1)
integer::i,j

do i=1,n2
do j=1,n1
B(i,j)=A(j,i)
end do
end do

return

end subroutine traspuesta


subroutine const_mat_diag(n,vecDiag,Mat_diag)
!Construye una matriz diagonal del vector que se le da, vecDiag, de n elementos
implicit none
integer, intent(in)::n
double precision, intent(in)::vecDiag(n,1)
double precision ,intent(out)::Mat_diag(n,n)
integer::i,j


do i=1,n
    do j=1,n
        Mat_diag(i,j)=0;
        end do
        end do

        do i=1,n
            do j=1,n
                if(i.eq.j) then
                    Mat_diag(i,j)=vecDiag(i,1)
                    end if

                end do
                end do
return
end subroutine const_mat_diag

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str




subroutine gbox(x0,y0,z0,Ndat,Nprism,malla,g)
!Modificado para modelo de capas
implicit none
! SUBRUTINA QUE CALCULA LA ATRACCI√ìN VERTICAL DE UN PRISMA RECT√ÅNGULAR
! Los lados del prisma son paralelos al eje x, y y z; donde z aumenta
! hacia abajo.
! Los par¬†metros a ingresar son:
! el punto de observaci¬¢n (x0,y0,z0).
! El prima se extiende de x1 a x2, de y1 a y2 y de z1 a z2.
! La densidad del prisma es rho en unidades kg/m^3
! Las distancias estan dadas en km.
! La atraccion vertical de gravedad g en mGal.
integer, intent(in)::Ndat,Nprism
double precision ,intent(in)::x0(Ndat),y0(Ndat),z0(Ndat),malla(Nprism,6)
double precision , dimension(:,:), intent(out)::g
integer::i,j,k,h,l
double precision ::sum,rijk,ijk,arg1,arg2,arg3
double precision ,dimension(2)::x,y,z,isign
double precision ,parameter::gamma=6.67e-11,twopi=6.2831853,si2mg=1.e5,km2m=1.e3
isign=(/-1,1/)

do h=1,Ndat

!do l=1,Nz-1
 !do m=1,Ny-1
  !do n=1,Nx-1
!conta=(l-1)*(Nx-1)*(Ny-1)+(m-1)*(Nx-1)+n
 do l=1,Nprism
 
y(1)=x0(h)-malla(l,1)
x(1)=y0(h)-malla(l,3)
z(1)=z0(h)-malla(l,5)
y(2)=x0(h)-malla(l,2)
x(2)=y0(h)-malla(l,4)
z(2)=z0(h)-malla(l,6)

sum=0.0
do i=1,2
   do j=1,2
      do k=1,2
         rijk=sqrt(x(i)**2+y(j)**2+z(k)**2)
         ijk=isign(i)*isign(j)*isign(k)
         arg1=atan2((x(i)*y(j)),(z(k)*rijk))
         if(arg1<0.0) then
            arg1=arg1+twopi
         end if
         arg2=rijk+y(j)
         arg3=rijk+x(i)
         if(arg2<=0.0) then
            write(6,*) 'GBOX: punto incorrecto del campo (arg2)'
         end if
         if(arg3<=0.0) then
            write(6,*) 'GBOX: punto incorrecto del campo (arg3)'
         endif
         arg2=log(arg2)
         arg3=log(arg3)
         sum=sum+ijk*(z(k)*arg1-x(i)*arg2-y(j)*arg3)
      end do
   end do
end do
g(h,l)=gamma*sum*si2mg*km2m

!end do !xM
!end do !yM
!end do !zM
end do !Nprism
!write(*,*) 'Ndat'
!write(*,*) h
end do !Datos h

return

end subroutine gbox

subroutine oper_derX(n,Dx)
implicit none
integer, intent(in)::n
double precision ,intent(out)::Dx(n,n)
integer::i,j
do i=1,n
do j=1,n
Dx(i,j)=0;
end do
end do
do i=2,n
do j=1,n
if(i.eq.j) then
Dx(i,j)=1
Dx(i,j-1)=-1
end if
end do
end do
return
end subroutine oper_DerX

subroutine oper_derY(n,nx,ny,nz,Dy)
implicit none
integer, intent(in)::n,nx,ny,nz
double precision ,intent(out)::Dy(n,n)
integer::i,j,k
do i=1,n
do j=1,n
Dy(i,j)=0
end do
end do


do i=1,Nz-1
do j=2,Ny-1
do k=1,Nx-1
Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),k+(j-2)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1))=-1
Dy(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1))=1
end do
end do
end do


return
end subroutine oper_derY

subroutine oper_derZ(n,nx,ny,nz,Dz)
implicit none
integer, intent(in)::n,nx,ny,nz
double precision ,intent(out)::Dz(n,n)
integer::i,j,k

do i=1,n
do j=1,n
Dz(i,j)=0
end do
end do


do i=2,Nz-1
do j=1,Ny-1
do k=1,Nx-1
Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),k+(j-1)*(Nx-1)+(i-2)*(Nx-1)*(Ny-1))=-1
Dz(k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1),k+(j-1)*(Nx-1)+(i-1)*(Nx-1)*(Ny-1))=1
end do
end do
end do

return
end subroutine oper_derZ





 subroutine mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,theta,gm)
  !   MODIFICADO PARA QUE SALGA SIN MULTIPLICAR POR LA MAGNETIZACION
!SOLO SALE LA SENSIBILIDAD
!  Subroutine MBOX computes the total field anomaly of an infinitely
!  extended rectangular prism.  Sides of prism are parallel to x,y,z
!  axes, and z is vertical down.  Bottom of prism extends to infinity.
!  Two calls to mbox can provide the anomaly of a prism with finite
!  thickness; e.g.,
!
!     call mbox(x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,m,theta,t1)
!     call mbox(x0,y0,z0,x1,y1,z2,x2,y2,mi,md,fi,fd,m,theta,t2)
!     t=t1-t2
!
!  Requires subroutine DIRCOS.  Method from Bhattacharyya (1964).
!
!  Input parameters:
!    Observation point is (x0,y0,z0).  Prism extends from x1 to
!    x2, y1 to y2, and z1 to infinity in x, y, and z directions,
!    respectively.  Magnetization defined by inclination mi,
!    declination md, intensity m.  Ambient field defined by
!    inclination fi and declination fd.  X axis has declination
!    theta. Distance units are irrelevant but must be consistent.
!    Angles are in degrees, with inclinations positive below
!    horizontal and declinations positive east of true north.
!    Magnetization in A/m.
!
!  Output paramters:
!    Total field anomaly t, in nT.
      double precision ,intent(in)::x0,y0,z0,x1,y1,z1,x2,y2,mi,md,fi,fd,theta
      double precision ,intent(out)::gm
      double precision :: alpha(2),beta(2),ma,mb,mc,t,fa,fb,fc
      double precision :: fm1,fm2,fm3,fm4,fm5,fm6,hsq,alphasq,sign
      double precision :: r0sq,r0h,alphabeta,arg3,arg4,tlog,tatan
      real :: arg1,arg2, h,r0 !alog solo permite reales
      !data cm/1.e-7/,t2nt/1.e9/
      double precision ,parameter::cm=1.e-7,t2nt=1.e9
      call dircos(mi,md,theta,ma,mb,mc)
      call dircos(fi,fd,theta,fa,fb,fc)
      fm1=ma*fb+mb*fa
      fm2=ma*fc+mc*fa
      fm3=mb*fc+mc*fb
      fm4=ma*fa
      fm5=mb*fb
      fm6=mc*fc
      alpha(1)=x1-x0
      alpha(2)=x2-x0
      beta(1)=y1-y0
      beta(2)=y2-y0
      h=z1-z0
      t=0.0
      hsq=h**2
      do 1 i=1,2
         alphasq=alpha(i)**2
         do 1 j=1,2
            sign=1.0
            if(i.ne.j)sign=-1.
            r0sq=alphasq+beta(j)**2+hsq
            r0=sqrt(r0sq)
            r0h=r0*h
            alphabeta=alpha(i)*beta(j)
            arg1=(r0-alpha(i))/(r0+alpha(i))
            arg2=(r0-beta(j))/(r0+beta(j))
            arg3=alphasq+r0h+hsq
            arg4=r0sq+r0h-alphasq
            
             if (arg1<0) then
            arg1=arg1*(-1)
            end if

             if (arg2<0) then
            arg2=arg2*(-1)
            end if
            
            tlog=fm3*alog(arg1)/2.+fm2*alog(arg2)/2.-fm1*alog(r0+h)
            tatan=-fm4*atan2(alphabeta,arg3)-fm5*atan2(alphabeta,arg4)+fm6*atan2(alphabeta,r0h)
    1 t=t+sign*(tlog+tatan)
      !t=t*m*cm*t2nt
      gm=t*cm*t2nt
      
      return
      end subroutine mbox



subroutine dircos(incl,decl,azim,a,b,c)
!
!  Subroutine DIRCOS computes direction cosines from inclination
!  and declination.
!
!  Input parameters:
!    incl:  inclination in degrees positive below horizontal.
!    decl:  declination in degrees positive east of true north.
!    azim:  azimuth of x axis in degrees positive east of north.
!
!  Output parameters:
!    a,b,c:  the three direction cosines.
      double precision ,intent(in)::incl,decl,azim
      double precision ,intent(out)::a,b,c
      double precision :: xincl,xdecl,xazim
      double precision ,parameter::d2rad=.017453293
      !data d2rad/.017453293/
      xincl=incl*d2rad
      xdecl=decl*d2rad
      xazim=azim*d2rad
      a=cos(xincl)*cos(xdecl-xazim)
      b=cos(xincl)*sin(xdecl-xazim)
      c=sin(xincl)
      return
      end subroutine dircos