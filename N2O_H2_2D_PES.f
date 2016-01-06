      IMPLICIT REAL*8 (A-H,O-Z)
      iv3=1;isoh2=1  
C	iv3=0 and 1 for the ground and excited vibrational state of N2O
C	isoh2=0 and 1 for para-H2 and ortho-D2 

C     Initialized fitting parameters before calculation
      CALL PARAREAD(iv3,isoh2)

      OPEN(11,FILE='h2-n2o-pes.chk')
      DO R=2.0D0, 10.0D0, 0.05D0
      DO TH=0.0D0,180.D0,1.D0
        CALL  H2N2OPES(R,TH,VCAL)
        IF(VCAL.GT.1000.D0) VCAL=1000.D0
        WRITE(11,663)R,TH,VCAL
      ENDDO
      ENDDO
 663  FORMAT(1x,f8.2,f8.2,f16.8)
      END
c=============================================
      
   
c***************************************************
      SUBROUTINE H2N2OPES(R,TH,YC)
c***************************************************
c 	Subroutine to generate values of the vibrationally averaged 
c  	2D-MLR analyic potential energy surfaces for complexes formed 
c  	between para-H2 or ortho-H2 and N2O isotopologue {14}N2{16}O 
c	in vibrational level v3= 0 or 1, as determined by:
c    	L. Wang, D. Xie, Robert J. Le Roy and P.-N. Roy,
c     in the paper 
c     "A new six-dimensional potential energy surface for H2--N2O
c     and its adiabatic-hindered-rotor treatment"

c  	Before first call, call input subroutine pararead(iv3,isoh2) for 
c  	parameters(depend on the vibrational state (iv3) of N2O
c	and the isotop of H2: isoh2=0 for para-H2 and 1 for ortho-H2)
c  	of MLR functions; 
c-------------------
c** Input variables:
c-------------------
c 	R - distance between the center of mass of N2O and H2 in [Angst], 
c 	pointing from the center of mass of N2O to H2.
c 	TH - Jacobi angular coordinate 'theta' in degrees, which is the 
c 	angle between the vector R pointing from the center of mass 
c 	of N2O to H2 and the vector pointing from O atom to N. 
c---------------------
c** 	Output:   YC [cm-1]  is the calculated interaction energy '\Delta{V}'.
c----------------------------------------------------

      INTEGER MXDATA, MXPARM, NPMAX, MXN, MXL
      PARAMETER (MXDATA=15000, MXPARM=80, NPMAX=20,MXN=20,MXL=20)

      INTEGER NPHI(NPMAX),LCM(0:MXN),I,J,K,L,M,IP,
     1 IDAT,NDATA,p,q,NDE,NRE,NCN,MCM,MMN,NS,NL,NPOW,NPS

      REAL*8  R,TH,YC,PV(MXPARM),
     1  Re,De,Vasy,RREF,AREF,AREFp,AREFq,Rep,CN,RC6,RCN,VLRe,
     2  phiINF,RTPp,RTPq,yp,yq,ype,yPOW,XP,SUM,VLR,XPW

      REAL*8 Pn(0:NPMAX+1),PHI(NPMAX),RCM(MXN,0:MXL) 
      COMMON /DATABLK/RCM,RREF,CN,
     1 NPHI,LCM,p,q,NDE,NRE,NCN,MCM,NS,NL,PV
c=======================================================================
c  For the case of an  MLR_{p}  potential ...
c-----------------------------------------------------------------------

       PI=DACOS(-1.0D0)
       CTH=DCOS(TH*PI/180.0D0)

       Pn(0)=1.
       Pn(1)=CTH
       DO I=1,NPMAX
        Pn(I+1)=((2.0d0*I+1.0d0)*CTH*Pn(I)
     1                         -dfloat(I)*Pn(I-1))/dfloat(I+1)
       ENDDO 
c caculate the derivative of the parameters of De not including
c the coefficient before, only the Legendre expansion.

       De=0.0d0
       DO I=0,NDE-1
         De=De+Pn(I)*PV(I+1)
       ENDDO 
       
c caculate the derivative of the parameters of Re not including
c the coefficient before, only the Legendre expansion.

       Re=0.0d0
       IP=NDE
       DO I=0,NRE-1
          IP=IP+1
          Re=Re+Pn(I)*PV(IP)
       ENDDO  

       AREF= RREF*Re
       IF(RREF.LE.0.d0) AREF= Re
       AREFp= AREF**p
       AREFq= AREF**q
       Rep= Re**p
c       write(7,*) Re,AREF,AREFp,AREFq
c only included the C6 coefficient 
       RC6=0.0D0
       DO L=0,LCM(NCN),2
        RC6=RC6+RCM(NCN,L)*Pn(L)
       ENDDO  
       VLRe= CN*RC6/Re**NCN
c included the higher order coefficients such as C7, C8, C9,C10 etc.
      IF(MCM.GT.NCN) THEN
       MMN=MCM-NCN
c      IF(p.LE.MMN)THEN MMN=0
      IF(MMN.GT.0) THEN
       RCN=0.0D0
       DO M=NCN,MCM
         IF (MOD(M,2).eq.0) then
           DO L=0,LCM(M),2
             RCN=RCN+(RCM(M,L)*Pn(L))/Re**(M-NCN)
           ENDDO 
         ELSE
           DO L=1,LCM(M),2
             RCN=RCN+(RCM(M,L)*Pn(L))/Re**(M-NCN)
           ENDDO
         ENDIF
         ENDDO 
         VLRe= CN*RCN/Re**NCN
      ENDIF
      ENDIF  

       phiINF= DLOG(2.d0*De/VLRe)

       RTPp= R**p
       RTPq= R**q
       yp= (RTPp - AREFp)/(RTPp + AREFp)
       yq= (RTPq - AREFq)/(RTPq + AREFq)
       ype= (RTPp - Rep)/(RTPp + Rep)
C       write(7,*) yp,yq,ype
 
c caculate the derivative of the parameters of PHI(N) not including
c the coefficient before, the Legendre expansion and exponent expansion.

       NPOW= NS+1
       IF(R.GE.Re) NPOW= NL+1
        yPOW= 1.d0 - yp
        SUM=0.0 
        NPS=0
        IP=NDE+NRE
        DO J=1,NPOW
          IP=IP+1
          PHI(J)= PV(IP)*Pn(0)
           DO  K=2,NPHI(J)
             IP=IP+1
            PHI(J)=PHI(J)+ PV(IP)*Pn(K-1)
           ENDDO
        NPS=NPS+NPHI(J)
        SUM=SUM+PHI(J)*yq**(J-1)
        ENDDO

caculate the derivative of the parameters of Vasy 

        IP=NDE+NRE+NPS
        Vasy=PV(IP+1) 

        XP= SUM*yPOW+ phiINF*yp


c only included the C6 coefficient 
       VLR= CN*RC6/R**NCN
c included the higher order coefficients such as C7, C8, C9,C10 etc.
      IF(MCM.GT.NCN) THEN
       RCN=0.0D0
       DO M=NCN,MCM
         IF (MOD(M,2).eq.0) then
           DO L=0,LCM(M),2
             RCN=RCN+(RCM(M,L)*Pn(L))/R**(M-NCN)
             ENDDO
         ELSE
           DO L=1,LCM(M),2
             RCN=RCN+(RCM(M,L)*Pn(L))/R**(M-NCN)
             ENDDO
         ENDIF
         ENDDO
         VLR= CN*RCN/R**NCN
       ENDIF 

         XPW= DEXP(-XP*ype) * VLR/VLRe
         YC= De*(1.d0 - XPW)**2-De+Vasy 

      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE PARAREAD(iv3,isoh2)
C	Subroutine to intial parameters defining the angle-averaged
C	two dimensional PES for para-H2--N2O and ortho-D2--N2O
C	Input parameters:
C	iv3: v3 vibrational level of N2O: iv3=0 for v3=0 and iv3=1 for v3=1
C	isoh2: isoh2=0 for para-H2 and isoh2=1 for ortho-D2 
      INTEGER MXDATA, MXPARM, NPMAX, MXN, MXL
      PARAMETER (MXDATA=15000, MXPARM=80, NPMAX=20,MXN=20,MXL=20)
      INTEGER iv3,isoh2
	INTEGER NPARM
	PARAMETER(NPARM=71)
      INTEGER NPHI(NPMAX),LCM(0:MXN),I,J,K,L,M,IP,
     1  p,q,NDE,NRE,NCN,MCM,NS,NL,NPOW,NPS
      REAL*8  PV(MXPARM),RCM(MXN,0:MXL)
      REAL*8  RREF,CN,PI
C     FOR GROUND AND EXCITED STATE OF N2O, AND pH2 and oD2
      REAL*8  TC62(4)
      REAL*8  TCN(4)
      REAL*8  TPV(4,MXPARM)
   
      COMMON /DATABLK/RCM,RREF,CN,
     1 NPHI,LCM,p,q,NDE,NRE,NCN,MCM,NS,NL,PV

C	PES DATA BEGIN HERE
C	EXPANSION LENGTH OF PARAMETERS
	DATA NDE/19/
	DATA NRE/15/
	DATA p/5/
	DATA q/3/
	DATA NS/5/
	DATA NL/5/
	DATA RREF/0.d0/
	DATA (NPHI(I),I=1,6)/11,9,7,5,3,1/

C	LONG RANGE PARAMETERS BEGIN HERE
	DATA NCN/6/
	DATA MCM/8/
	DATA (LCM(I),I=6,8)/2,3,4/
	DATA (TCN(I),I=1,4)/2.3132926D5,2.3298419D5,
     & 2.3132926D5,2.3298419D5/
	DATA (TC62(I),I=1,4)/0.326423D0,0.3266091D0,
     & 0.326423D0,0.3266091D0/
	DATA RCM(6,0)/1.D0/
	DATA (RCM(7,I),I=1,3,2)/0.56341392D0,0.03768641/
	DATA (RCM(8,I),I=0,4,2)/10.24082D0,19.32325D0,1.33809D0/

C	FITTING PARAMETERS BEGIN HERE
	DATA (TPV(1,I),I=1,NPARM)/
     &    8.51780D+01, -7.73000D+00, -5.26900D+01,  4.22000D+00,
     &    4.23600D+01, -9.50000D+00, -2.35900D+01,  6.80000D+00,
     &    1.22100D+01, -3.83000D+00, -5.70000D+00,  8.20000D-01,
     &    1.85000D+00,  7.00000D-02, -2.30000D-01, -2.00000D-02,
     &    8.00000D-02,  0.00000D+00,  0.00000D+00,  3.71400D+00,
     &    1.51480D-01,  9.69180D-01, -2.60000D-03, -2.79900D-01,
     &    5.49000D-02,  9.03000D-02, -2.78000D-02, -2.57000D-02,
     &    1.15000D-02,  1.00000D-03, -2.40000D-03,  2.40000D-03,
     &    2.00000D-04, -6.00000D-04, -2.33300D-01,  5.37000D-02,
     &    4.09300D-01,  1.10000D-02,  5.90000D-02,  1.10000D-02,
     &   -6.20000D-02, -1.60000D-02,  2.50000D-02,  3.00000D-03,
     &   -1.00000D-02,  2.03000D-01, -2.00000D-02,  4.36000D-01,
     &    0.00000D+00, -2.00000D-02,  0.00000D+00,  2.00000D-02,
     &   -1.00000D-02, -7.00000D-02, -2.70000D-01, -2.10000D-01,
     &    0.00000D+00,  0.00000D+00,  9.00000D-02,  1.30000D-01,
     &    1.60000D-01, -1.00000D-01, -4.30000D-01,  0.00000D+00,
     &   -4.40000D-01, -4.50000D-01, -5.00000D-01,  1.20000D+00,
     &    3.00000D+00,  2.91000D+00,  0.00000D+00/

	DATA (TPV(2,I),I=1,NPARM)/
     &    8.51430D+01, -7.82000D+00, -5.21400D+01,  4.29000D+00,
     &    4.19400D+01, -9.61000D+00, -2.32000D+01,  6.84000D+00,
     &    1.19800D+01, -3.85000D+00, -5.53000D+00,  8.20000D-01,
     &    1.78000D+00,  7.00000D-02, -2.20000D-01, -2.00000D-02,
     &    7.00000D-02,  0.00000D+00,  0.00000D+00,  3.71624D+00,
     &    1.52420D-01,  9.68300D-01, -3.40000D-03, -2.78600D-01,
     &    5.56000D-02,  8.87000D-02, -2.78000D-02, -2.55000D-02,
     &    1.15000D-02,  8.00000D-04, -2.60000D-03,  2.40000D-03,
     &    2.00000D-04, -6.00000D-04, -2.33800D-01,  5.38000D-02,
     &    4.09600D-01,  1.10000D-02,  6.00000D-02,  1.30000D-02,
     &   -6.10000D-02, -1.60000D-02,  2.60000D-02,  3.00000D-03,
     &   -1.00000D-02,  1.97000D-01, -1.00000D-02,  4.35000D-01,
     &    0.00000D+00, -1.00000D-02,  0.00000D+00,  1.00000D-02,
     &   -1.00000D-02, -7.00000D-02, -2.90000D-01, -2.20000D-01,
     &    0.00000D+00,  0.00000D+00,  9.00000D-02,  1.30000D-01,
     &    1.50000D-01,  0.00000D+00, -4.70000D-01,  0.00000D+00,
     &   -4.00000D-01, -4.90000D-01, -4.00000D-01,  1.30000D+00,
     &    3.00000D+00,  2.63000D+00,  0.00000D+00/

	DATA (TPV(3,I),I=1,NPARM)/
     &    8.81300D+01, -8.44000D+00, -5.60200D+01,  5.52000D+00,
     &    4.52900D+01, -1.03000D+01, -2.57900D+01,  7.52000D+00,
     &    1.35300D+01, -4.28000D+00, -6.42000D+00,  1.08000D+00,
     &    2.22000D+00, -5.00000D-02, -4.20000D-01,  3.00000D-02,
     &    1.10000D-01,  0.00000D+00,  0.00000D+00,  3.69893D+00,
     &    1.53500D-01,  9.83500D-01, -8.90000D-03, -2.89600D-01,
     &    5.58000D-02,  9.74000D-02, -2.89000D-02, -2.99000D-02,
     &    1.23000D-02,  3.00000D-03, -3.30000D-03,  2.90000D-03,
     &    6.00000D-04, -1.60000D-03, -2.39400D-01,  5.20000D-02,
     &    4.33000D-01,  1.10000D-02,  4.80000D-02,  1.70000D-02,
     &   -5.70000D-02, -1.80000D-02,  2.60000D-02,  6.00000D-03,
     &   -1.40000D-02,  1.80000D-01, -1.00000D-02,  4.35000D-01,
     &    0.00000D+00, -5.00000D-02,  0.00000D+00,  4.00000D-02,
     &   -1.00000D-02, -9.00000D-02, -2.90000D-01, -1.90000D-01,
     &    9.00000D-02,  0.00000D+00,  0.00000D+00,  1.20000D-01,
     &    1.90000D-01,  1.00000D-01, -5.70000D-01,  3.00000D-01,
     &   -4.70000D-01, -3.00000D-01, -7.00000D-01,  1.30000D+00,
     &    2.30000D+00,  2.72000D+00,  0.00000D+00/

	DATA (TPV(4,I),I=1,NPARM)/
     &    8.80040D+01, -8.54000D+00, -5.53700D+01,  5.58000D+00,
     &    4.47500D+01, -1.04100D+01, -2.53300D+01,  7.56000D+00,
     &    1.32400D+01, -4.30000D+00, -6.21000D+00,  1.08000D+00,
     &    2.13000D+00, -5.00000D-02, -4.10000D-01,  3.00000D-02,
     &    1.10000D-01,  0.00000D+00,  0.00000D+00,  3.70165D+00,
     &    1.54700D-01,  9.82700D-01, -9.70000D-03, -2.88300D-01,
     &    5.64000D-02,  9.59000D-02, -2.86000D-02, -2.95000D-02,
     &    1.22000D-02,  2.70000D-03, -3.60000D-03,  2.80000D-03,
     &    7.00000D-04, -1.60000D-03, -2.40800D-01,  5.10000D-02,
     &    4.34000D-01,  1.10000D-02,  4.80000D-02,  1.80000D-02,
     &   -5.80000D-02, -1.80000D-02,  2.60000D-02,  6.00000D-03,
     &   -1.30000D-02,  1.75000D-01,  0.00000D+00,  4.37000D-01,
     &    0.00000D+00, -5.00000D-02, -1.00000D-02,  3.00000D-02,
     &    0.00000D+00, -9.00000D-02, -2.80000D-01, -1.90000D-01,
     &    8.00000D-02,  0.00000D+00,  0.00000D+00,  1.50000D-01,
     &    2.00000D-01,  2.00000D-01, -6.10000D-01,  3.00000D-01,
     &   -4.60000D-01, -3.00000D-01, -7.00000D-01,  1.40000D+00,
     &    2.30000D+00,  2.46000D+00,  0.00000D+00/

C	FOR DIFFERENT ISOTOPES OF H2 AND VIBRATIONAL STATES OF N2O
	IF((IV3.EQ.0).AND.(ISOH2.EQ.0)) THEN
	 CN=TCN(1)
	 RCM(6,2)=TC62(1)
	 DO I=1,NPARM
	  PV(I)=TPV(1,I)
	 ENDDO
	ELSE IF((IV3.EQ.1).AND.(ISOH2.EQ.0)) THEN
       CN=TCN(2)
       RCM(6,2)=TC62(2)
       DO I=1,NPARM
        PV(I)=TPV(2,I)
       ENDDO
      ELSE IF((IV3.EQ.0).AND.(ISOH2.EQ.1)) THEN
       CN=TCN(3)
       RCM(6,2)=TC62(3)
       DO I=1,NPARM
        PV(I)=TPV(3,I)
       ENDDO
      ELSE IF((IV3.EQ.1).AND.(ISOH2.EQ.1)) THEN
       CN=TCN(4)
       RCM(6,2)=TC62(4)
       DO I=1,NPARM
        PV(I)=TPV(4,I)
       ENDDO
	ELSE
	 WRITE(*,*) 'ERROR IN VIBRATIONAL STATE OF N2O OR ISOTOPES OF H2'
	 STOP
	ENDIF
      RETURN
      END
