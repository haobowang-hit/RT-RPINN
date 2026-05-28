       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
!  --------------------------------------------------------------------------------

      INTEGER nProny, counter, i, j
      REAL*8 C11inf, C12inf, C22inf, C23inf, C66inf
      REAL*8 c111,c122,c123,c121,c222,c233,c232,c223,c444,c665,c666
      REAL*8 c1, c2, Tr, aT, dtr
      REAL*8 rhoi(6),strainE(6),Tp(6)
      PARAMETER (EPS=2.22D-16) !SMALLEST NUMBER REAL*8 CAN STORE
      DIMENSION C11i(6),C12i(6),C22i(6),C23i(6),C66i(6)
      DIMENSION q(6,6),qold(6,6)
      

!      A = 1.
!      EPS   = EPSILON(A)        !NEGLIGIBLE COMPARED TO 1.0, SAME TYPE AS A
!      TINY2 = TINY(A)           !SMALLEST NUMBER OF SAME TYPE AS A
!      HUGE2 = HUGE(A)           !LARGEST NUMBER OF SAME TYPE AS A
! -----------------------------------------------------------
!     UMAT FOR 3D SOLID ELEMENTS
!     F77 IMPLICIT NAME CONVENTION
! -----------------------------------------------------------
!       NDI: # of direct components (11,...) of DDSDDE, DDSDDT, and DRPLDE
!       NSHR: # of engineering shear components (12,...) of DDSDDE, DDSDDT, and DRPLDE
!       NTENS = NDI + NSHR: Size of the stress or strain component array
!       TIME(1):    Value of step time at the beginning of the current increment.
!       TIME(2):    Value of total time at the beginning of the current increment.
!       DTIME:      Time increment.
!       STRESS(NTENS):  passed in as the stress tensor at the beginning of the increment
!                       must be updated to be the stress tensor at the end of the increment
! -----------------------------------------------------------
C-----material parameters---------C
C Number of Prony terms
      nProny=6
      
  
C  Obtain Prony coefficients
      do i=1,nProny
         rhoi(i)=PROPS(i)
         C11i(i)=PROPS(7+i)
         C12i(i)=PROPS(14+i)
         C22i(i)=PROPS(21+i)
         C23i(i)=PROPS(28+i)
         C66i(i)=PROPS(35+i)
      end do
      
C Obtain long term modulus
      C11inf=PROPS(7)
      C12inf=PROPS(14)
      C22inf=PROPS(21)
      C23inf=PROPS(28)
      C66inf=PROPS(35) 
      Tp(1)=PROPS(42) 
      Tp(2)=PROPS(43)
      Tp(3)=PROPS(43)
      Tp(4)=0
      Tp(5)=0
      Tp(6)=0
C WLF parameters
      Tr=323
      c1=14.8
      c2=45.6
C----------solution-dependent variables----------C


      
C     old heredity integrals
         counter = 1
      do j = 1,NTENS
         do i = 1,nProny
        qold(j,i) = STATEV(counter)
        counter = counter+1
      end do
         end do
      
      
C----------update stress--------------C
      
      
C     generalized section strains at the end of increment
      strainE(1)=STRAN(1) + DSTRAN(1)
      strainE(2)=STRAN(2) + DSTRAN(2)
      strainE(3)=STRAN(3) + DSTRAN(3)
      strainE(4)=STRAN(4) + DSTRAN(4)
      strainE(5)=STRAN(5) + DSTRAN(5)
      strainE(6)=STRAN(6) + DSTRAN(6)
      
      
C     shift factor
      if (Temp+Dtemp.GT.317.4) then
      aT=10**(-c1*(TEMP+DTEMP-Tr)/(c2+TEMP+DTEMP-Tr))
      else
      aT=2.71828**(27403.3*(1/(TEMP+DTEMP)-1/336.0))  
      end if
      
      
C     reduced time increment
      dtr=DTIME/aT

      
C     compute q
      do j = 1, NTENS
        do i = 1, nProny
          q(j,i)=exp(-dtr/rhoi(i))*qold(j,i)+
     1 (DSTRAN(j)-Tp(j)*DTEMP)*(1-exp(-dtr/rhoi(i)))/(dtr/rhoi(i))
        end do
      end do

C     compute a
      c111=0
      c122=0
      c123=0
      c121=0
      c222=0
      c233=0
      c232=0
      c223=0
      c444=0
      c665=0
      c666=0

        c111=c111+C11i(1)*q(1,1)+C11i(2)*q(1,2)
     1 +C11i(3)*q(1,3)+C11i(4)*q(1,4)+C11i(5)*q(1,5)
     1 +C11i(6)*q(1,6)
        c122=c122+C12i(1)*q(2,1)+C12i(2)*q(2,2)
     1 +C12i(3)*q(2,3)+C12i(4)*q(2,4)+C12i(5)*q(2,5)
     1 +C12i(6)*q(2,6)
        c123=c123+C12i(1)*q(3,1)+C12i(2)*q(3,2)
     1 +C12i(3)*q(3,3)+C12i(4)*q(3,4)+C12i(5)*q(3,5)
     1 +C12i(6)*q(3,6)
        c121=c121+C12i(1)*q(1,1)+C12i(2)*q(1,2)
     1 +C12i(3)*q(1,3)+C12i(4)*q(1,4)+C12i(5)*q(1,5)
     1 +C12i(6)*q(1,6)
        c222=c222+C22i(1)*q(2,1)+C22i(2)*q(2,2)
     1 +C22i(3)*q(2,3)+C22i(4)*q(2,4)+C22i(5)*q(2,5)
     1 +C22i(6)*q(2,6)
        c233=c233+C23i(1)*q(3,1)+C23i(2)*q(3,2)
     1 +C23i(3)*q(3,3)+C23i(4)*q(3,4)+C23i(5)*q(3,5)
     1 +C23i(6)*q(3,6)
        c232=c232+C23i(1)*q(2,1)+C23i(2)*q(2,2)
     1 +C23i(3)*q(2,3)+C23i(4)*q(2,4)+C23i(5)*q(2,5)
     1 +C23i(6)*q(2,6)
        c223=c223+C22i(1)*q(3,1)+C22i(2)*q(3,2)
     1 +C22i(3)*q(3,3)+C22i(4)*q(3,4)+C22i(5)*q(3,5)
     1 +C22i(6)*q(3,6)
        c444=c444+(C22i(1)-C23i(1))/2*q(4,1)+(C22i(2)-C23i(2))/2*q(4,2)
     1 +(C22i(3)-C23i(3))/2*q(4,3)+(C22i(4)-C23i(4))/2*q(4,4)
     1 +(C22i(5)-C23i(5))/2*q(4,5)+(C22i(6)-C23i(6))/2*q(4,6)
        c665=c665+C66i(1)*q(5,1)+C66i(2)*q(5,2)
     1 +C66i(3)*q(5,3)+C66i(4)*q(5,4)+C66i(5)*q(5,5)+C66i(6)*q(5,6)
        c666=c666+C66i(1)*q(6,1)+C66i(2)*q(6,2)
     1 +C66i(3)*q(6,3)+C66i(4)*q(6,4)+C66i(5)*q(6,5)+C66i(6)*q(6,6)

      
      STRESS(1) = C11inf*strainE(1)+c111+C12inf*strainE(2)+c122+
     1 C12inf*strainE(3)+c123
      STRESS(2) = C12inf*strainE(1)+c121+C22inf*strainE(2)+c222+
     1 C23inf*strainE(3)+c233
      STRESS(3) = C12inf*strainE(1)+c121+C23inf*strainE(2)+c232+
     1 C22inf*strainE(3)+c223
      STRESS(4) = (C22inf-C23inf)/2*strainE(4)+c444
      STRESS(5) = C66inf*strainE(5)+c665
      STRESS(6) = C66inf*strainE(6)+c666
      
C-----update solution-dependent variables-----C
       
C     hereditary integrals
      counter = 1
      do j = 1,NTENS
        do i = 1, nProny
          STATEV(counter)=q(j,i)
          counter=counter+1
        end do
      end do     
            

C-----update Jacobian (tangent stiffess)-----C
      do i = 1, NTENS
        do j = 1, NTENS
          DDSDDE(i,j)=0.0
        end do
      end do
       
       

        DDSDDE(1,1)=DDSDDE(1,1)+C11inf+
     1 C11i(1)*(1-exp(-dtr/rhoi(1)))/(dtr/rhoi(1))
     1 +C11i(2)*(1-exp(-dtr/rhoi(2)))/(dtr/rhoi(2))
     1 +C11i(3)*(1-exp(-dtr/rhoi(3)))/(dtr/rhoi(3))
     1 +C11i(4)*(1-exp(-dtr/rhoi(4)))/(dtr/rhoi(4))
     1 +C11i(5)*(1-exp(-dtr/rhoi(5)))/(dtr/rhoi(5))
     1 +C11i(6)*(1-exp(-dtr/rhoi(6)))/(dtr/rhoi(6))  
        
        DDSDDE(1,2)=DDSDDE(1,2)+C12inf+
     1 C12i(1)*(1-exp(-dtr/rhoi(1)))/(dtr/rhoi(1))
     1 +C12i(2)*(1-exp(-dtr/rhoi(2)))/(dtr/rhoi(2))
     1 +C12i(3)*(1-exp(-dtr/rhoi(3)))/(dtr/rhoi(3))
     1 +C12i(4)*(1-exp(-dtr/rhoi(4)))/(dtr/rhoi(4))
     1 +C12i(5)*(1-exp(-dtr/rhoi(5)))/(dtr/rhoi(5))
     1 +C12i(6)*(1-exp(-dtr/rhoi(6)))/(dtr/rhoi(6))
           
        DDSDDE(2,2)=DDSDDE(2,2)+C22inf+
     1 C22i(1)*(1-exp(-dtr/rhoi(1)))/(dtr/rhoi(1))
     1 +C22i(2)*(1-exp(-dtr/rhoi(2)))/(dtr/rhoi(2))
     1 +C22i(3)*(1-exp(-dtr/rhoi(3)))/(dtr/rhoi(3))
     1 +C22i(4)*(1-exp(-dtr/rhoi(4)))/(dtr/rhoi(4))
     1 +C22i(5)*(1-exp(-dtr/rhoi(5)))/(dtr/rhoi(5))
     1 +C22i(6)*(1-exp(-dtr/rhoi(6)))/(dtr/rhoi(6))
        
        DDSDDE(2,3)=DDSDDE(2,3)+C23inf+
     1 C23i(1)*(1-exp(-dtr/rhoi(1)))/(dtr/rhoi(1))
     1 +C23i(2)*(1-exp(-dtr/rhoi(2)))/(dtr/rhoi(2))
     1 +C23i(3)*(1-exp(-dtr/rhoi(3)))/(dtr/rhoi(3))
     1 +C23i(4)*(1-exp(-dtr/rhoi(4)))/(dtr/rhoi(4))
     1 +C23i(5)*(1-exp(-dtr/rhoi(5)))/(dtr/rhoi(5))
     1 +C23i(6)*(1-exp(-dtr/rhoi(6)))/(dtr/rhoi(6)) 
        
        DDSDDE(6,6)=DDSDDE(6,6)+C66inf+
     1 C66i(1)*(1-exp(-dtr/rhoi(1)))/(dtr/rhoi(1))
     1 +C66i(2)*(1-exp(-dtr/rhoi(2)))/(dtr/rhoi(2))
     1 +C66i(3)*(1-exp(-dtr/rhoi(3)))/(dtr/rhoi(3))
     1 +C66i(4)*(1-exp(-dtr/rhoi(4)))/(dtr/rhoi(4))
     1 +C66i(5)*(1-exp(-dtr/rhoi(5)))/(dtr/rhoi(5))
     1 +C66i(6)*(1-exp(-dtr/rhoi(6)))/(dtr/rhoi(6))
      
      
      DDSDDE(2,1)=DDSDDE(1,2)
      DDSDDE(3,1)=DDSDDE(1,2)
      DDSDDE(3,2)=DDSDDE(2,3)
      DDSDDE(1,3)=DDSDDE(1,2)
      DDSDDE(3,3)=DDSDDE(2,2)
      DDSDDE(4,4)=(DDSDDE(2,2)-DDSDDE(2,3))/2
      DDSDDE(5,5)=DDSDDE(6,6)

       
C----------update energy----------C
      
      
C  update elastic strain energy in SSE
      SSE = 0.0
      
      
C  update plastic dissipation in SPD
      SPD = 0.0
      

      
      RETURN
      END   
C***********************************************************************
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
      
      
      INCLUDE 'ABA_PARAM.INC'
      
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
      
      
C     Initialize state variables
      do i = 1,NSTATV
        STATEV(i) = 0.0
      end do
   
   
      RETURN
      END