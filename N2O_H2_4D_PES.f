      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION R,THN2O,THH2,PHI,POT

C     Initialized fitting parameters before calculation
      CALL PARAREAD(0)    
C	SUBROUTINE PARAREAD SHOULD BE CALLED BEFORE CALCULATING POTENTIAL ENERGY
C	PARAREAD(IV),IV=0 FOR v3=0 STATE OF N2O AND IV=1 FOR v3=1 STATE OF N2O
      PI=DACOS(-1.D0)
 
      DO R=3.0D0,3.6D0,0.2D0
       DO THN2O=0.D0,180.D0,30.D0
        DO THH2=0.D0,180.D0,30.D0
         DO PHI=0.D0,90.D0,30.D0
          CALL H2N2O4DPES(R,THN2O*PI/180.D0,THH2*PI/180.D0,
     & PHI*PI/180.D0,POT)
          WRITE(11,663) R,THN2O,THH2,PHI,POT
         ENDDO
        ENDDO
       ENDDO
      ENDDO
 663  FORMAT(1X,f7.1,3f9.1,f20.8)
      STOP
      END 

c***********************************************************************
      SUBROUTINE H2N2O4DPES(RTP,THN2O,THH2,ANPHI,YC)
c******************************************************************************
c	Subroutine to generate values of the 4D-MLR analyic
c 	potential energy surfaces for complexes formed between 
c  	H2 and N2O, as determined by:
c    	L. Wang, D. Xie, Robert J. Le Roy and P.-N. Roy,
c     in the paper 
c     "A new six-dimensional potential energy surface for H2--N2O
c     and its adiabatic-hindered-rotor treatment"
 
c  	before first call, call subroutine pararead(iv) for 
c  	the initialization of parameters
c-------------------
c	 Input variables:
c-------------------
c 	RTP - distance between N2O and He centre of mass in [Angst], 
c       	pointing from the center of mass of N2O to He.
c 	THN2O - Jacobi angular coordinate in arc, 
c         	which is the angle between the vector RTP 
c         	pointing from the center of mass of N2O to H2 
c         	and the vector pointing from O atom to N. 
c 	THH2 - Jacobi angular coordinate in arc,
c         	which is the angle between RTP and the axis of H2
c 	ANPHI - Jacobi angular coorinates of the dihedral angle
c---------------------
c Output: V [cm-1]  is the calculated interaction energy '\Delta{V}'.
c--------------------------------------------------------------------------
      INTEGER MXDATA, MXPARM, NPMAX, MXN, MXL
      PARAMETER (MXDATA=15000, MXPARM=300, NPMAX=30,MXN=30,MXL=30)

      INTEGER IFXP(MXPARM),NPHI(NPMAX),LCM(0:MXN),I,J,K,L,M,IP,
     1 IDAT,NDATA,NPARM,p,q,NDE,NRE,NCN,MCM,MMN,NS,NL,NPOW,NPS,
     2  L1DE(MXPARM),L2DE(MXPARM),LDE(MXPARM),L1RE(MXPARM),
     3  L2RE(MXPARM),LRE(MXPARM),L1PHI(MXPARM,MXPARM),
     4  L2PHI(MXPARM,MXPARM),LPHI(MXPARM,MXPARM)

      REAL*8  YC,PV(MXPARM),Re,De,Vasy,RREF,AREF,AREFp,AREFq,Rep,CN,RC6,
     1 RCN,VLRe,phiINF,RTPp,RTPq,yp,yq,ype,yPOW,XP,SUM,VLR,XPW,CL(MXN)

      REAL*8 DIN2O,DISN2O,POAVGN2O,PODELN2O,QUAN2O,
     1  POAVGH2,POPAH2,POPEH2,QUAH2

      REAL*8 AL1L2L0,PI
      REAL*8  RTP,THN2O,RCM(MXN,0:MXL),PHI(NPMAX),
     1  THH2,ANPHI
      COMMON /DATABLK/p,q,NS,NL,RREF,
     1 NPHI,L1PHI,L2PHI,LPHI,NDE,L1DE,L2DE,LDE,NRE,L1RE,L2RE,LRE,PV
c=======================================================================
c  For the case of an  MLR_{p}  potential ...
c-----------------------------------------------------------------------
       CALL GETLRPARA(THN2O,THH2,ANPHI,CL(4),CL(5),CL(6),CL(7),CL(8))
       pi=dacos(-1.d0)
       De=0.0d0
       DO I=0,NDE-1
         De=De+AL1L2L0(L1DE(I+1),L2DE(I+1),LDE(I+1),THN2O,
     1 THH2,ANPHI)*PV(I+1)
       ENDDO
c       write(*,*) 'detest',de 
       
c caculate the derivative of the parameters of Re not including
c the coefficient before, only the Legendre expansion.

       Re=0.0d0
       IP=NDE
       DO I=0,NRE-1
          IP=IP+1
          Re=Re+AL1L2L0(L1RE(I+1),L2RE(I+1),LRE(I+1),THN2O,
     1 THH2,ANPHI)*PV(IP)
       ENDDO
c      write(*,*) 'retest:',re  

       AREF= RREF*Re
       IF(RREF.LE.0.d0) AREF= Re
       AREFp= AREF**p
       AREFq= AREF**q
       Rep= Re**p

c included the higher order coefficients such as C7, C8, C9,C10 etc.
       RCN=0.0D0
       DO M=4,8
        RCN=RCN+CL(M)/Re**M
       END DO
       VLRe= RCN
C       write(*,*) de,vlre
       phiINF= DLOG(2.d0*De/VLRe)

       RTPp= RTP**p
       RTPq= RTP**q
       yp= (RTPp - AREFp)/(RTPp + AREFp)
       yq= (RTPq - AREFq)/(RTPq + AREFq)
       ype= (RTPp - Rep)/(RTPp + Rep)
 
c caculate the derivative of the parameters of PHI(N) not including
c the coefficient before, the Legendre expansion and exponent expansion.

       NPOW= NS+1
       IF(RTP.GE.Re) NPOW= NL+1
        yPOW= 1.d0 - yp
        SUM=0.0 
        NPS=0
        IP=NDE+NRE
        DO J=1,NPOW
          IP=IP+1
          PHI(J)= PV(IP)*AL1L2L0(L1PHI(J,1),L2PHI(J,1),LPHI(J,1),
     1 THN2O,THH2,ANPHI)
c          write(*,*) ip,L1PHI(J,1),L2PHI(J,1),LPHI(J,1)
           DO  K=2,NPHI(J)
             IP=IP+1
            PHI(J)=PHI(J)+ PV(IP)*AL1L2L0(L1PHI(J,K),L2PHI(J,K),
     1 LPHI(J,K),THN2O,THH2,ANPHI)
c           write(*,*) ip,L1PHI(J,K),L2PHI(J,K),LPHI(J,K)
           ENDDO
c        write(*,*) 'phitest:',j,phi(j)
        NPS=NPS+NPHI(J)
        SUM=SUM+PHI(J)*yq**(J-1)
        ENDDO

caculate the derivative of the parameters of Vasy 

        IP=NDE+NRE+NPS
        Vasy=PV(IP+1) 

        XP= SUM*yPOW+ phiINF*yp


c included the higher order coefficients such as C4-C7 etc.
       RCN=0.0D0
       DO M=4,8
             RCN=RCN+CL(M)/RTP**M
       ENDDO
         VLR= RCN
c         write(*,99) rtp(idat),thn2o(idat)*180.d0/pi,
c     &thh2(idat)*180.d0/pi,anphi(idat)*180.d0/pi,rcn
         XPW= DEXP(-XP*ype) * VLR/VLRe
         YC= De*(1.d0 - XPW)**2-De+Vasy 

      RETURN
      END


c***********************************************************************
c    
c This subroutine fills up the matrix w3j for C6 with values of 3-j coefficient
c     
      subroutine fill3j(l1max,l2max,lmax)
      implicit real*8 (a-h,o-z)
      common /w3jcg/w3j(0:50,0:30,0:30,0:60)
      dimension x(60)

      do l1=0,l1max
       do l2=0,l2max
        lmin=iabs(l1-l2)
        mm=min0(l1,l2)
        llmax=min0(lmax,l1+l2)

        do l=lmin,llmax
         do m=0,mm
          m1=m
          m2=-m
          mmm=0
          call cgc(l1,m1,l2,m2,l,mmm,c,1)
          w3j(m,l1,l2,l)=c
         end do
        end do
       end do
      end do
c      write(*,*) w3j(2,2,2,2)
      return
      end
C ----------------------------------------------------------------------------
c
c Calculate the function C6Al1l2L for a given set of angles...
c It is assumed that the th1 and th2 angles are between the monomer bond
c and the "inner" part of the intermolecular axis.
c
       function al1l2l0(l1,l2,l,th1,th2,phi)
       implicit real*8 (a-h,o-z)
       dimension p1(0:50,0:50), p2(0:50,0:50)
       common /w3jcg/w3j(0:50,0:30,0:30,0:60)
       data izer/0/, ione/1/, pifact/12.566370614359d0/
c
       c1 = dcos(th1)
       c2 = dcos(th2)
c       write(*,*) w3j(0,0,0,0)
       call plmrb(p1,c1,l1)
       call plmrb(p2,c2,l2)
       mmax = min(l1,l2)
       sum = 0.d0
       do m=1,mmax
        value=w3j(m,l1,l2,l)
        sum = sum + (-1)**m*value*p1(l1,m)*p2(l2,m)*dcos(m*phi)
       end do
       value=w3j(0,l1,l2,l)
       sum = 2*sum + value*p1(l1,0)*p2(l2,0)
c
       al1l2l0 = sum*pifact*(-1.d0)**(l1+l2+l)/dsqrt((2.d0*l1+1.d0)*
     1           (2.d0*l2+1.d0))
c
c       write(*,1001) l1,l2,l,value
c       write(*,*) 'ytest'
c       write(*,1002) l1,c2,p2(l1,0)
c       write(*,*) 'atest'
c       write(*,1003) l1,l2,l,th1,th2,phi,al1l2l0
1001   format(1x,3i4,f14.6)
1002   format(1x,i4,f10.4,f14.6)
1003   format(1X,3I4,3F8.4,F12.6)

       return
       end

C --------------------------------------------------------------------------
c
c Compute the set of associated Legendre polynomials P_lm
c for l=0,1,...,lmax, and m=0,1,...,l. First the standard
c polynomials
c
c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
c
c are computed, and then multiplied by
c
c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
c
c to get the P_lm polynomials....
c
        subroutine plmrb(p,x,lmax)
        implicit real*8 (a-h,o-z)
        dimension p(0:50,0:50)
        common/factorial/ fact(0:40)
c inverse of dsqrt(2Pi)
        data twopinv /0.3989422804014d0/
c
c starting value
c
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)
c
c compute the diagonal elements
c
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do
c
c compute P_lm along the columns with fixed m

c
        do m = 0,lmax-1
        do l = m,lmax-1
         if((l-1).lt.m) then
           pp = 0
         else
           pp = p(l-1,m)
         endif
         p(l+1,m) = ((2*l+1)*x*p(l,m)-(l+m)*pp)/(l-m+1)
        end do
        end do
c
c Renormalize values...
c
        do l=0,lmax
        mm = 1
        do m=0,l
         dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
         p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
         mm = -mm
        end do
        end do
c
        return
        end

C -------------------------------------------------------------------------
c
c compute the matrix of N!
c
        subroutine fct(nmax)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c
        f(0) = 1.d0
        do i=1,nmax
         f(i) = f(i-1)*i
        end do
        return
        end
C -------------------------------------------------------------------------
c
Calculate the Clebsh-Gordan coefficient (or the 3-j symbol)
c The parameter ind3j.eq.1 indicates that the 3-J symbol is returned
c
        subroutine cgc(j1,m1,j2,m2,j,m,value,ind3j)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c
        d3jfact = 1.d0
        if(ind3j.eq.1) then
         d3jfact = ((-1.d0)**(j1-j2-m))/dsqrt(dfloat(2*j+1))
         m = -m
        endif

c
c Check the triangle conditions
c
        if(j.gt.(j1+j2)) write(6,*)'triangle violated'
        if(j.lt.abs(j1-j2)) write(6,*)'triangle violated'
        if((m1+m2).ne.m) then
          value = 0.d0
          return
        endif


c Calculation proper... the pre-sum factor....
c
        facn = (2*j+1)*f(j1+j2-j)*f(j1-m1)*f(j2-m2)*f(j+m)*f(j-m)
        facd = f(j1+j2+j+1)*f(j+j1-j2)*f(j+j2-j1)*f(j1+m1)*f(j2+m2)
        fac = dsqrt(facn/facd)

c
c determine the limit of k summation...
c
        kmax = min(j2+j-m1,j-m,j1-m1)
        if(kmax.lt.0) kmax = 0
        kmin = max(-j1-m1,-j2+j-m1,0)

c
c perform the summation (at least one cycle must be completed...
c
        sum = 0.d0
        do k = kmin,kmax
         facn = f(j1+m1+k)*f(j2+j-m1-k)
         facd = f(k)*f(j-m-k)*f(j1-m1-k)*f(j2-j+m1+k)
         sum = sum + (facn/facd)*(-1)**k
        end do
        value = d3jfact*fac*sum*(-1)**(j1-m1)
       return
       end
C -------------------------------------------------------------------------
       SUBROUTINE flush(nunit)
       endfile nunit
       backspace nunit
       end
C -------------------------------------------------------------------------


      SUBROUTINE PARAREAD(IV)
C	INITIAL ALL THE PARAS REQUIRED FOR PES
C	INPUT: IV=0, v3=0 STATE OF N2O; IV=1, v3=1 STATE OF N2O
	INTEGER IV
      INTEGER MXDATA, MXPARM, NPMAX, MXN, MXL
      PARAMETER (MXDATA=15000, MXPARM=300, NPMAX=30,MXN=30,MXL=30)
      INTEGER NPARM,ND,NR,NP1,NP2,NP3,NP4,NP5
	PARAMETER (NPARM=286,ND=82,NR=72,NP1=48,NP2=40,NP3=17,NP4=17,
     1 NP5=9)
      REAL*8 DIN2O,DISN2O,POAVGN2O,PODELN2O,QUAN2O,
     1  POAVGH2,POPAH2,POPEH2,QUAH2
c-----------------------------------------------------------------
C***FOR TEST
      INTEGER ITEST

      INTEGER NPHI(NPMAX),I,J,K,L,M,IP,
     1  p,q,NDE,NRE,NCN,MCM,NS,NL,NPOW,NPS,
     1  L1DE(MXPARM),L2DE(MXPARM),LDE(MXPARM),L1RE(MXPARM),
     1  L2RE(MXPARM),LRE(MXPARM),L1PHI(MXPARM,MXPARM),
     1  L2PHI(MXPARM,MXPARM),LPHI(MXPARM,MXPARM)
      REAL*8  PV(MXPARM)
      REAL*8  RREF,CN,PI,AL1L2L0
      INTEGER ITMP
      REAL*8  RCM(30,0:30),
     1  THH2(MXDATA),ANPHI(MXDATA)
C	FOLLOWING VALUES ARE DIFFERENT FOR v3=0 AND 1 STATE
	REAL*8  TCN(2),TRCM(2,30,0:30),TPV(2,MXPARM),TDIN2O(2),
     1 TDISN2O(2),TPOAVGN2O(2),TPODELN2O(2),TQUAN2O(2)
	
      COMMON /DATABLK/p,q,NS,NL,RREF,
     1 NPHI,L1PHI,L2PHI,LPHI,NDE,L1DE,L2DE,LDE,
     1 NRE,L1RE,L2RE,LRE,PV
      COMMON /LONGR/CN,RCM,NCN,MCM,DIN2O,DISN2O,POAVGN2O,
     1  PODELN2O,QUAN2O,POAVGH2,POPAH2,POPEH2,QUAH2
 
c-----------------------------------------------------------------------
      CALL FCT(40)
      CALL FILL3J(14,14,14)

C	PES BEGIN HERE
C	EXPANSION LENGTH BEGIN HERE
      DATA NDE/82/
	DATA NRE/72/
	DATA p/5/
	DATA q/3/
	DATA NS/4/
	DATA NL/4/
	DATA RREF/0.D0/
	DATA (NPHI(I),I=1,5) /48,40,17,17,9/

C	LONG RANGE PARAMETERS BEGIN HERE	
	DATA NCN/6/
	DATA (TCN(I),I=1,2) /2.3132926d5,2.3298419d5/
	DATA MCM/8/
      DATA (TRCM(1,6,J),J=0,2) /1.0D0,0.326423D0,0.373827D0/
	DATA (TRCM(1,7,J),J=1,3,2) /0.56341392D0,0.03768641D0/
	DATA (TRCM(1,8,J),J=0,3) /10.24082D0,19.32325D0,1.33809D0,
     & 3.00158D0/
      DATA (TRCM(2,6,J),J=0,2) /1.0D0,0.3266091D0,0.373827D0/
      DATA (TRCM(2,7,J),J=1,3,2) /0.56341392D0,0.03768641D0/
      DATA (TRCM(2,8,J),J=0,3) /10.24082D0,19.32325D0,1.33809D0,
     & 3.00158D0/
	DATA (TDIN2O(I),I=1,2) /0.064252D0,0.065108d0/
 	DATA (TDISN2O(I),I=1,2) /0.017137D0,0.03747d0/
	DATA (TPOAVGN2O(I),I=1,2) /19.663929D0,19.804605d0/
	DATA (TPODELN2O(I),I=1,2) /19.256302D0,19.405092d0/
	DATA (TQUAN2O(I),I=1,2) /-2.6381d0,-2.58948d0/
	DATA POAVGH2 /5.414D0/
	DATA POPAH2 /6.7632D0/
	DATA POPEH2 /4.7393D0/
	DATA QUAH2 /0.481D0/

C	PES L1 L2 L START HERE
	DATA (L1DE(I),I=1,ND)/
     &    0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2,
     &    2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  4,
     &    4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,
     &    5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,
     &    6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,
     &    8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,
     &    10, 10, 10, 10, 11, 11, 11, 12, 12, 12/
	DATA (L2DE(I),I=1,ND)/
     &    0,  2,  4,  0,  2,  2,  4,  4,  0,  2,  2,  2,
     &    4,  4,  4,  0,  2,  2,  2,  4,  4,  4,  4,  0,
     &    2,  2,  2,  4,  4,  4,  4,  4,  0,  2,  2,  2,
     &    4,  4,  4,  4,  4,  0,  2,  2,  2,  4,  4,  4,
     &    4,  4,  0,  2,  2,  2,  4,  4,  4,  4,  4,  0,
     &    2,  2,  2,  4,  4,  4,  4,  4,  0,  2,  2,  2,
     &    0,  2,  2,  2,  0,  2,  2,  0,  2,  2/
	DATA (LDE(I),I=1,ND)/
     &    0,  2,  4,  1,  1,  3,  3,  5,  2,  0,  2,  4,
     &    2,  4,  6,  3,  1,  3,  5,  1,  3,  5,  7,  4,
     &    2,  4,  6,  0,  2,  4,  6,  8,  5,  3,  5,  7,
     &    1,  3,  5,  7,  9,  6,  4,  6,  8,  2,  4,  6,
     &    8, 10,  7,  5,  7,  9,  3,  5,  7,  9, 11,  8,
     &    6,  8, 10,  4,  6,  8, 10, 12,  9,  7,  9, 11,
     &    10, 8, 10, 12, 11,  9, 11, 12, 10, 12/
	DATA (L1RE(I),I=1,NR)/
     &    0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2,
     &    2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  4,
     &    4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,
     &    5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,
     &    6,  6,  7,  7,  7,  7,  8,  8,  8,  8,  9,  9,
     &    9,  9, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12/
	DATA (L2RE(I),I=1,NR)/
     &    0,  2,  4,  0,  2,  2,  4,  4,  0,  2,  2,  2,
     &    4,  4,  4,  0,  2,  2,  2,  4,  4,  4,  4,  0,
     &    2,  2,  2,  4,  4,  4,  4,  4,  0,  2,  2,  2,
     &    4,  4,  4,  4,  4,  0,  2,  2,  2,  4,  4,  4,
     &    4,  4,  0,  2,  2,  2,  0,  2,  2,  2,  0,  2,
     &    2,  2,  0,  2,  2,  2,  0,  2,  2,  0,  2,  2/
	DATA (LRE(I),I=1,NR)/
     &    0,  2,  4,  1,  1,  3,  3,  5,  2,  0,  2,  4,
     &    2,  4,  6,  3,  1,  3,  5,  1,  3,  5,  7,  4,
     &    2,  4,  6,  0,  2,  4,  6,  8,  5,  3,  5,  7,
     &    1,  3,  5,  7,  9,  6,  4,  6,  8,  2,  4,  6,
     &    8, 10,  7,  5,  7,  9,  8,  6,  8, 10,  9,  7,
     &    9, 11, 10,  8, 10, 12, 11,  9, 11, 12, 10, 12/
	DATA (L1PHI(1,I),I=1,NP1)/
     &    0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2,
     &    2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  4,
     &    4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,
     &    6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8/
	DATA (L2PHI(1,I),I=1,NP1)/
     &    0,  2,  4,  0,  2,  2,  4,  4,  0,  2,  2,  2,
     &    4,  4,  4,  0,  2,  2,  2,  4,  4,  4,  4,  0,
     &    2,  2,  2,  4,  4,  4,  4,  4,  0,  2,  2,  2,
     &    0,  2,  2,  2,  0,  2,  2,  2,  0,  2,  2,  2/
	DATA (LPHI(1,I),I=1,NP1)/
     &    0,  2,  4,  1,  1,  3,  3,  5,  2,  0,  2,  4,
     &    2,  4,  6,  3,  1,  3,  5,  1,  3,  5,  7,  4,
     &    2,  4,  6,  0,  2,  4,  6,  8,  5,  3,  5,  7,
     &    6,  4,  6,  8,  7,  5,  7,  9,  8,  6,  8, 10/
	DATA (L1PHI(2,I),I=1,NP2)/
     &    0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2,
     &    2,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  4,
     &    4,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,
     &    6,  6,  6,  6/
	DATA (L2PHI(2,I),I=1,NP2)/
     &    0,  2,  4,  0,  2,  2,  4,  4,  0,  2,  2,  2,
     &    4,  4,  4,  0,  2,  2,  2,  4,  4,  4,  4,  0,
     &    2,  2,  2,  4,  4,  4,  4,  4,  0,  2,  2,  2,
     &    0,  2,  2,  2/
	DATA (LPHI(2,I),I=1,NP2)/
     &    0,  2,  4,  1,  1,  3,  3,  5,  2,  0,  2,  4,
     &    2,  4,  6,  3,  1,  3,  5,  1,  3,  5,  7,  4,
     &    2,  4,  6,  0,  2,  4,  6,  8,  5,  3,  5,  7,
     &    6,  4,  6,  8/
	DATA (L1PHI(3,I),I=1,NP3)/0,0,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4/
	DATA (L2PHI(3,I),I=1,NP3)/0,2,0,2,2,0,2,2,2,0,2,2,2,0,2,2,2/
	DATA (LPHI(3,I),I=1,NP3) /0,2,1,1,3,2,0,2,4,3,1,3,5,4,2,4,6/
	DATA (L1PHI(4,I),I=1,NP4)/0,0,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4/
	DATA (L2PHI(4,I),I=1,NP4)/0,2,0,2,2,0,2,2,2,0,2,2,2,0,2,2,2/
	DATA (LPHI(4,I),I=1,NP4) /0,2,1,1,3,2,0,2,4,3,1,3,5,4,2,4,6/
	DATA (L1PHI(5,I),I=1,NP5)/0,0,1,1,1,2,2,2,2/
	DATA (L2PHI(5,I),I=1,NP5)/0,2,0,2,2,0,2,2,2/
	DATA (LPHI(5,I),I=1,NP5) /0,2,1,1,3,2,0,2,4/

C	PES PARAS START HERE
	DATA (TPV(1,I),I=1,NPARM)/
     &    8.45100D+01, -3.71400D+01,  5.80000D+00,  1.26200D+01,
     &    1.65000D+01,  6.04000D+01, -2.70000D+00, -4.00000D+00,
     &   -1.18040D+02,  2.42500D+01,  3.72000D+01,  4.05800D+02,
     &   -4.30000D+00,  5.20000D+00, -8.00000D+00, -9.80000D+00,
     &    1.32000D+01, -2.20000D+00, -6.83000D+01,  2.10000D+00,
     &    3.00000D-01,  4.20000D+00,  1.39000D+01,  1.26100D+02,
     &   -2.89000D+01,  4.30000D+00, -2.44900D+02,  1.60000D+00,
     &    0.00000D+00,  6.80000D+00,  5.80000D+00,  4.00000D+01,
     &    3.02000D+01, -1.57000D+01, -8.20000D+00, -3.73000D+01,
     &    1.80000D+00, -5.00000D-01,  9.00000D-01, -2.10000D+00,
     &   -1.80000D+01, -8.37000D+01,  2.35000D+01, -8.30000D+00,
     &    1.16400D+02, -2.50000D+00,  1.40000D+00, -3.80000D+00,
     &   -7.00000D-01, -3.49000D+01, -2.54000D+01,  1.25000D+01,
     &    3.90000D+00,  2.20000D+01, -1.40000D+00,  0.00000D+00,
     &   -7.00000D-01, -1.30000D+00, -3.30000D+00,  4.94000D+01,
     &   -1.50000D+01,  6.00000D+00, -6.79000D+01,  1.60000D+00,
     &   -6.00000D-01,  1.80000D+00,  0.00000D+00,  1.49000D+01,
     &    1.68000D+01, -5.60000D+00, -2.60000D+00, -1.77000D+01,
     &   -2.52000D+01,  1.12000D+01, -3.40000D+00,  2.48000D+01,
     &   -4.20000D+00,  6.00000D-01, -1.30000D+00,  8.20000D+00,
     &   -1.01000D+01,  1.80000D+00,  3.73612D+00,  2.00300D-01,
     &    6.00000D-03, -2.57700D-01, -2.76000D-02, -1.30800D-01,
     &   -5.00000D-03, -8.00000D-03,  2.14210D+00, -4.98000D-02,
     &   -3.00000D-02, -1.17700D+00,  7.00000D-03,  3.00000D-03,
     &   -1.20000D-02, -2.11000D-02,  0.00000D+00,  4.00000D-03,
     &    2.62000D-01,  0.00000D+00,  0.00000D+00,  8.00000D-03,
     &    4.20000D-02, -8.18900D-01,  7.90000D-02, -5.40000D-02,
     &    4.87000D-01, -2.00000D-03,  7.00000D-03,  8.00000D-03,
     &    3.40000D-02,  2.10000D-01, -1.87400D-01,  1.60000D-02,
     &    2.10000D-02,  2.56000D-01, -2.00000D-03,  0.00000D+00,
     &    0.00000D+00, -4.00000D-03, -8.20000D-02,  2.91000D-01,
     &   -5.70000D-02,  3.00000D-02,  1.05000D-01,  8.00000D-03,
     &   -3.00000D-03,  1.00000D-02,  0.00000D+00, -1.20000D-01,
     &    1.07000D-01, -2.10000D-02,  7.00000D-03, -3.90000D-02,
     &   -8.00000D-02,  2.50000D-02,  0.00000D+00, -9.80000D-02,
     &   -5.10000D-02,  2.10000D-02, -5.00000D-03,  9.00000D-03,
     &   -8.00000D-03,  7.00000D-03, -1.50000D-02,  6.20000D-02,
     &    1.00000D-02, -5.00000D-03,  0.00000D+00,  1.20000D-02,
     &   -2.40000D-02,  7.00000D-03, -3.05000D-01,  2.69000D-01,
     &   -8.10000D-02, -9.80000D-02, -3.10000D-02,  2.00000D-02,
     &    4.50000D-02, -4.00000D-02,  1.25100D+00,  1.78000D-01,
     &   -5.92000D-01, -7.60000D-01, -5.00000D-02,  2.22000D-01,
     &   -2.00000D-01,  0.00000D+00, -2.90000D-02,  3.00000D-02,
     &    1.50000D-01, -3.10000D-02, -2.00000D-02,  8.00000D-02,
     &    2.20000D-01, -1.45000D-01,  5.00000D-02,  0.00000D+00,
     &    2.60000D-01, -2.00000D-02,  6.00000D-02, -7.00000D-02,
     &    2.50000D-01,  7.20000D-01, -6.00000D-02, -2.00000D-02,
     &    0.00000D+00,  3.00000D-02, -8.00000D-02, -3.00000D-02,
     &    3.00000D-02, -1.20000D-01,  7.00000D-02,  3.00000D-02,
     &    0.00000D+00, -1.80000D-01,  5.00000D-02,  0.00000D+00,
     &    0.00000D+00,  6.00000D-02,  1.97000D-01,  8.70000D-01,
     &   -8.00000D-02,  3.20000D-01,  0.00000D+00, -2.10000D-01,
     &   -7.00000D-02, -5.70000D-01, -4.00000D-01,  7.60000D-01,
     &   -2.21000D+00, -4.30000D+00, -2.40000D-01,  6.60000D-01,
     &   -1.12000D+00, -9.00000D-02,  1.00000D-01,  0.00000D+00,
     &    4.00000D-01,  0.00000D+00,  0.00000D+00,  4.00000D-01,
     &    1.10000D+00,  1.70000D-01,  3.00000D-01,  0.00000D+00,
     &    1.00000D+00, -1.10000D-01,  2.10000D-01,  0.00000D+00,
     &    1.10000D+00,  3.60000D+00,  0.00000D+00,  0.00000D+00,
     &    0.00000D+00, -3.00000D-01, -3.30000D-01,  0.00000D+00,
     &    4.00000D-01, -3.00000D-01, -9.00000D-02,  1.20000D+00,
     &    4.60000D-01,  2.00000D-01, -9.00000D-01, -1.00000D+00,
     &    1.50000D+00, -1.10000D+00, -6.50000D+00, -6.00000D-01,
     &   -3.00000D-01,  3.00000D-01,  4.50000D+00,  4.00000D-01,
     &   -6.00000D-01,  7.00000D-01,  4.70000D+00,  6.90000D-01,
     &    2.60000D+00,  1.20000D+00,  3.00000D-01, -8.00000D-01,
     &   -2.60000D+00,  4.70000D+00, -9.30000D+00, -2.24000D+01,
     &   -1.10000D+00, -9.00000D-01,  1.30000D+00,  8.50000D+00,
     &    6.00000D-01, -3.60000D+00,  0.00000D+00,  4.00000D+00,
     &   -1.04000D+00,  4.50000D+00,  1.20000D+00,  0.00000D+00,
     &    0.00000D+00, -1.00000D+00,  3.40000D+00, -2.17000D+01,
     &   -3.78000D+01,  0.00000D+00/
      DATA (TPV(2,I),I=1,NPARM)/
     &    8.45100D+01, -3.61400D+01,  5.60000D+00,  1.28000D+01,
     &    1.63000D+01,  5.99000D+01, -2.70000D+00, -3.90000D+00,
     &   -1.16900D+02,  2.38000D+01,  3.70000D+01,  3.99300D+02,
     &   -4.20000D+00,  5.10000D+00, -7.60000D+00, -1.00000D+01,
     &    1.31000D+01, -1.90000D+00, -6.65000D+01,  2.00000D+00,
     &    3.00000D-01,  4.20000D+00,  1.35000D+01,  1.25000D+02,
     &   -2.82000D+01,  3.90000D+00, -2.44100D+02,  1.50000D+00,
     &    0.00000D+00,  6.70000D+00,  5.60000D+00,  3.89000D+01,
     &    3.06000D+01, -1.56000D+01, -8.20000D+00, -3.85000D+01,
     &    1.80000D+00, -5.00000D-01,  9.00000D-01, -2.10000D+00,
     &   -1.75000D+01, -8.24000D+01,  2.29000D+01, -8.10000D+00,
     &    1.14200D+02, -2.40000D+00,  1.40000D+00, -3.70000D+00,
     &   -8.00000D-01, -3.48000D+01, -2.56000D+01,  1.24000D+01,
     &    3.90000D+00,  2.24000D+01, -1.40000D+00,  0.00000D+00,
     &   -7.00000D-01, -1.30000D+00, -3.40000D+00,  4.85000D+01,
     &   -1.45000D+01,  5.90000D+00, -6.67000D+01,  1.50000D+00,
     &   -6.00000D-01,  1.70000D+00,  0.00000D+00,  1.48000D+01,
     &    1.69000D+01, -5.60000D+00, -2.60000D+00, -1.79000D+01,
     &   -2.45000D+01,  1.08000D+01, -3.40000D+00,  2.42000D+01,
     &   -4.20000D+00,  6.00000D-01, -1.30000D+00,  7.90000D+00,
     &   -9.70000D+00,  1.80000D+00,  3.73772D+00,  1.98300D-01,
     &    6.00000D-03, -2.58900D-01, -2.77000D-02, -1.30100D-01,
     &   -5.00000D-03, -7.00000D-03,  2.14010D+00, -4.93000D-02,
     &   -3.00000D-02, -1.15800D+00,  7.00000D-03,  3.00000D-03,
     &   -1.20000D-02, -1.89000D-02,  0.00000D+00,  4.00000D-03,
     &    2.54000D-01,  0.00000D+00,  0.00000D+00,  9.00000D-03,
     &    4.30000D-02, -8.16800D-01,  7.90000D-02, -5.20000D-02,
     &    4.98000D-01, -2.00000D-03,  7.00000D-03,  8.00000D-03,
     &    3.30000D-02,  2.01000D-01, -1.89500D-01,  1.60000D-02,
     &    2.20000D-02,  2.57000D-01, -2.00000D-03,  0.00000D+00,
     &    0.00000D+00, -4.00000D-03, -8.10000D-02,  2.87000D-01,
     &   -5.70000D-02,  2.80000D-02,  1.05000D-01,  8.00000D-03,
     &   -3.00000D-03,  1.10000D-02,  0.00000D+00, -1.19000D-01,
     &    1.07000D-01, -2.10000D-02,  7.00000D-03, -3.70000D-02,
     &   -7.90000D-02,  2.70000D-02,  0.00000D+00, -9.90000D-02,
     &   -5.10000D-02,  2.20000D-02, -5.00000D-03,  8.00000D-03,
     &   -9.00000D-03,  6.00000D-03, -1.40000D-02,  6.40000D-02,
     &    1.10000D-02, -5.00000D-03,  0.00000D+00,  1.20000D-02,
     &   -2.40000D-02,  7.00000D-03, -3.04500D-01,  2.63000D-01,
     &   -7.90000D-02, -9.80000D-02, -3.00000D-02,  2.00000D-02,
     &    4.30000D-02, -4.00000D-02,  1.25700D+00,  1.75000D-01,
     &   -5.74000D-01, -7.50000D-01, -5.40000D-02,  2.17000D-01,
     &   -1.90000D-01,  0.00000D+00, -3.00000D-02,  3.00000D-02,
     &    1.50000D-01, -2.90000D-02, -1.00000D-02,  8.00000D-02,
     &    2.10000D-01, -1.41000D-01,  6.00000D-02,  0.00000D+00,
     &    2.60000D-01, -2.10000D-02,  5.80000D-02, -7.00000D-02,
     &    2.40000D-01,  7.10000D-01, -7.00000D-02, -2.00000D-02,
     &    0.00000D+00,  4.00000D-02, -8.00000D-02, -3.00000D-02,
     &    3.00000D-02, -1.20000D-01,  7.00000D-02,  3.00000D-02,
     &    0.00000D+00, -1.80000D-01,  5.00000D-02,  0.00000D+00,
     &    0.00000D+00,  6.00000D-02,  2.01000D-01,  8.40000D-01,
     &   -7.00000D-02,  3.10000D-01,  0.00000D+00, -2.20000D-01,
     &   -7.00000D-02, -5.30000D-01, -3.90000D-01,  7.70000D-01,
     &   -2.14000D+00, -4.29000D+00, -2.50000D-01,  6.40000D-01,
     &   -1.08000D+00, -1.00000D-01,  1.00000D-01,  0.00000D+00,
     &    4.00000D-01,  0.00000D+00,  0.00000D+00,  4.00000D-01,
     &    1.10000D+00,  1.80000D-01,  4.00000D-01,  0.00000D+00,
     &    1.00000D+00, -1.20000D-01,  2.10000D-01,  0.00000D+00,
     &    1.00000D+00,  3.40000D+00,  0.00000D+00,  0.00000D+00,
     &    0.00000D+00, -3.00000D-01, -3.60000D-01,  0.00000D+00,
     &    3.00000D-01, -3.00000D-01, -9.00000D-02,  1.20000D+00,
     &    4.80000D-01,  2.00000D-01, -8.00000D-01, -1.00000D+00,
     &    1.40000D+00, -1.10000D+00, -6.40000D+00, -6.20000D-01,
     &   -3.00000D-01,  3.00000D-01,  4.40000D+00,  4.00000D-01,
     &   -6.00000D-01,  6.00000D-01,  4.70000D+00,  7.20000D-01,
     &    2.60000D+00,  1.20000D+00,  3.00000D-01, -7.00000D-01,
     &   -2.60000D+00,  4.50000D+00, -9.10000D+00, -2.18000D+01,
     &   -1.10000D+00, -8.00000D-01,  1.30000D+00,  8.30000D+00,
     &    7.00000D-01, -3.90000D+00,  0.00000D+00,  4.00000D+00,
     &   -1.00000D+00,  4.40000D+00,  1.10000D+00,  0.00000D+00,
     &    0.00000D+00, -1.00000D+00,  3.30000D+00, -2.12000D+01,
     &   -3.71000D+01,  0.00000D+00/

C	PES DATAS ENDS HERE

C	FOR GROUND AND EXCITED ROTATIONAL STATE OF N2O
	IF(IV.EQ.0) THEN
	 CN=TCN(1)
	 DO I=0,2
	  RCM(6,I)=TRCM(1,6,I)
	 ENDDO
	 DO I=1,3,2
	  RCM(7,I)=TRCM(1,7,I)
	 ENDDO
	 DO I=0,3
	  RCM(8,I)=TRCM(1,8,I)
	 ENDDO
	 DIN2O=TDIN2O(1)
	 DISN2O=TDISN2O(1)
	 POAVGN2O=TPOAVGN2O(1)
	 PODELN2O=TPODELN2O(1)
	 QUAN2O=TQUAN2O(1)
	DO I=1,NPARM
	  PV(I)=TPV(1,I)
	 ENDDO
	ELSE IF (IV.EQ.1) THEN
       CN=TCN(2)
       DO I=0,2
        RCM(6,I)=TRCM(2,6,I)
       ENDDO
       DO I=1,3,2
        RCM(7,I)=TRCM(2,7,I)
       ENDDO
       DO I=0,3
        RCM(8,I)=TRCM(2,8,I)
       ENDDO
       DIN2O=TDIN2O(2)
       DISN2O=TDISN2O(2)
       POAVGN2O=TPOAVGN2O(2)
       PODELN2O=TPODELN2O(2)
       QUAN2O=TQUAN2O(2)
       DO I=1,NPARM
        PV(I)=TPV(2,I)
       ENDDO
	ELSE
	 WRITE(*,*) "ERROR IN THE VIBRATIONAL STATE OF N2O!"
	 STOP
      ENDIF
      RETURN
      END

c*******************************
      SUBROUTINE GETLRPARA(THN,THH,PHII,C4,C5,C6,C7,C8)
C***** TO GET LONG RANGE PARAMETERS WITH ANGLES **************
C***** INPUT: ALL IN ARC(0 TO PI) ****************
C***** OUTPUT: C4-C8 IN ORDER OF ANGSTROM*CM-1 ****************
      PARAMETER (MXDATA=15000, MXPARM=300, NPMAX=30,MXM=30,MXL=30)
      REAL*8 PI,AL1L2L0
      INTEGER NCN,MCM
      REAL*8 THN,THH,THH1,PHII,C4,C5,C6,C7,C8,C6IND,C6DISP,C7IND,C7DISP
      REAL*8 CN,RCM(30,0:30),DIN2O,DISN2O,POAVGN2O,tmp,
     1  PODELN2O,QUAN2O,POAVGH2,POPAH2,POPEH2,QUAH2
      COMMON /LONGR/CN,RCM,NCN,MCM,DIN2O,DISN2O,POAVGN2O,
     1  PODELN2O,QUAN2O,POAVGH2,POPAH2,POPEH2,QUAH2
      real*8 c6ind1

      PI=DACOS(-1.0D0)
      COSN2O=DCOS(THN)
      SINN2O=DSIN(THN)
      THH1=PI-THH
      COSH2=DCOS(THH1)
      SINH2=DSIN(THH1)
      COSPHI=DCOS(PHII)
      C4=(3.0D0/2.0D0)*DIN2O*QUAH2*(COSN2O*(3.0D0*COSH2**2-1.0D0)+
     1 2.0D0*SINH2*COSH2*SINN2O*COSPHI)*219474.D0*(0.5291771D0**4)
      C5=0.75D0*QUAN2O*QUAH2*9107.307339*
     & (1.0-5.0*COSN2O**2-5.0*COSH2**2
     & -15.0*COSN2O**2*COSH2**2+
     & 2.0*(4.0*COSN2O*COSH2-SINN2O*SINH2*CPHI)**2)

C      C5=(3.0D0/4.0D0)*QUAN2O*QUAH2*(1.0D0-5.0D0*COSN2O**2.D0-5.0D0
C     1 *COSH2**2.D0+17.0D0*COSN2O**2.D0*COSH2**2.D0+2.D0*SINN2O**2
C     1 *SINH2**2*COSPHI**2+16.D0*SINN2O*COSN2O*SINH2*COSH2*COSPHI)
C     1 *219474*(0.5291771D0**5)
      C6DISP=CN+CN*RCM(6,1)*AL1L2L0(2,0,2,THN,THH1,PHII)
     1 +CN*RCM(6,2)*AL1L2L0(0,2,2,THN,THH1,PHII)
      C7DISP=CN*RCM(7,1)*AL1L2L0(1,0,1,THN,THH1,PHII)+
     1 CN*RCM(7,3)*AL1L2L0(3,0,3,THN,THH1,PHII)
      C6IND=(POAVGH2*DISN2O*(3.D0*COSN2O**2+1.D0)/2.D0+(1.0d0/6.0d0)*
     1 (POPAH2-POPEH2)*DISN2O*(12.D0*COSN2O**2*COSH2**2+3.D0*
     1 SINN2O**2*SINH2**2*COSPHI**2-3.D0*COSN2O**2-1.D0+
     1 12.D0*SINN2O*COSN2O*SINH2*COSH2*COSPHI))
     1 *219474*(0.5291771**6)
      C6ind1=4819.379496d0*DISN2O*(0.5d0*POAVGH2*(3.0d0*COSN2O**2+1.0d0)
     & +(1.0d0/6.0d0)*(popah2-popeh2)*(12.0d0*Cosn2o**2*cosh2**2
     & +3.0D0*Sinn2o**2*Sinh2**2*Cosphi**2-3.0D0*Cosn2o**2-1.0D0
     & +12.0D0*Sinh2*Cosh2*Sinn2o*Cosn2o*cosphi))
c      write(*,*) c6ind,c6ind1
      C7IND=6.D0*POAVGH2*DIN2O*QUAH2*COSN2O**3*219474*(0.5291771**7)
      C4=-1.0D0*C4
      C5=-1.0D0*C5
      C6=C6DISP+C6IND
      C7=C7DISP+C7IND
      C8=CN*RCM(8,0)+CN*RCM(8,1)*AL1L2L0(2,0,2,THN,THH1,PHII)+
     & CN*RCM(8,2)*AL1L2L0(4,0,4,THN,THH1,PHII)+
     & CN*RCM(8,3)*AL1L2L0(0,2,2,THN,THH1,PHII)
c      WRITE(*,*) 'lrtest:'
c      WRITE(*,410) c4,c5,c6,c7,c8
c  410 FORMAT(1X,4f18.6)
      RETURN
      END
