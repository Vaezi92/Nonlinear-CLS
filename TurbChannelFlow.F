      PROGRAM SBLRE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C  1.    A 3-DIMENSIONAL FULLY ELLITPIC K-E TURBULENT FLOW SOLVER     C
C  2.    GENERAL NONORTHOGONAL GRID                                   C
C  3.    NONSTAGGERED FININTE VOLUME FORMULATION                      C
C  4.    SIMPLE PRESSURE CORRECTION ALGORITHM                         C
C  5.    LAUNDER AND SHARMA K-EPSILON  MODEL                          C 
C  6.    LOW-RE DSM MODELS                                            C
C  7.    TEMPERATURE EQUATION                                         C
C  8     ROTATION TERMS                                               C
C  9.    BOUNDED QUICK SCHEME                                         C
C                                                                     C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C*	
C*
	PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C*
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
C	CHARACTER KEY
C	LOGICAL(4) CHK_NAN /.FALSE./
C	LOGICAL(4) PRESSED /.FALSE./
C	INTEGER FINAL_IT,ISTAT1,ISTAT2
C      INTEGER(2) tmphour1, tmpminute1, tmpsecond1, tmphund1
C	INTEGER(2) tmphour2, tmpminute2, tmpsecond2, tmphund2
C	INTEGER(2) HOUR,MINUTE,SECOND
C*
C
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
C
      DIMENSION XLM(NZ),IG(4),KG(4),UYZ(4,NY,NZ),VYZ(4,NY,NZ),
     1   WYZ(4,NY,NZ),UXY(NX,NY,4),VXY(NX,NY,4),WXY(NX,NY,4)
     1  ,UUXY(NX,NY,4),VVXY(NX,NY,4),UVXY(NX,NY,4)
      DIMENSION JG(4),UXZ(NX,4,NZ),VXZ(NX,4,NZ),
     1   WXZ(NX,4,NZ),UUXZ(NX,4,NZ),VVXZ(NX,4,NZ),UVXZ(NX,4,NZ)
C
      dimension uex(ny,NZ),vex(ny,NZ),wex(ny,NZ),tkex(ny,NZ),
     1          tedx(ny,NZ),uuex(ny,NZ),vvex(NY,NZ),wwex(ny,NZ),
     1          visex(ny,NZ),pex(ny,NZ),TEX(NY,NZ),uvex(ny,nz),
     1          uwex(ny,nz),vwex(ny,nz),cxex(ny,nz)
c
      LOGICAL INCALU,INCALV,INCALW,INCALP,INCALK,INCALE,INPROP
      OPEN(16,FILE='UIND.DAT')
      OPEN(17,FILE='VIND.DAT')
      OPEN(18,FILE='TEIND.DAT')
      OPEN(7,FILE='UVLFL.DAT')
      OPEN(8,FILE='VVLFL.DAT')
      OPEN(9,FILE='WVLFL.DAT')
      OPEN(10,FILE='PRSFL.DAT')
      OPEN(11,FILE='TKEFL.DAT')
      OPEN(12,FILE='TEDFL.DAT')
      OPEN(13,FILE='XGRID.DAT')
      OPEN(14,FILE='YGRID.DAT')
      OPEN(15,FILE='ZGRID.DAT')
      OPEN(21,FILE='CVXFL.DAT')
      OPEN(22,FILE='CVYFL.DAT')
      OPEN(23,FILE='CVZFL.DAT')
      OPEN(24,FILE='VISFL.DAT')
      OPEN(39,FILE='TMPFL.DATA')
      OPEN(25,FILE='YCGRD.DATA')
      OPEN(26,FILE='ZCGRD.DATA')
      OPEN(27,FILE='XCGRD.DATA')
C*
      OPEN(28,FILE='DATNS.DATA')
      OPEN(29,FILE='DATSS.DATA')
      OPEN(30,FILE='UUDAT.DATA')
      OPEN(31,FILE='VVDAT.DATA')
      OPEN(32,FILE='WWDAT.DATA')
      OPEN(33,FILE='UVDAT.DATA')
      OPEN(34,FILE='UWDAT.DATA')
      OPEN(35,FILE='VWDAT.DATA')
C*
      OPEN(36,FILE='UIN1.DATA')
      OPEN(37,FILE='TEIN1.DATA')
      OPEN(38,FILE='VIN1.DATA')
      OPEN(40,FILE='UGRP.DATA')
      OPEN(41,FILE='VGRP.DATA')
      OPEN(401,FILE='WGRP.DATA')
      OPEN(42,FILE='ELXY.DAT')
      OPEN(45,FILE='ELXYC.DAT')
      OPEN(51,FILE='EXPLD.DAT')
      OPEN(52,FILE='RESIDUAL.DAT')
      OPEN(53,FILE='TEMP.DAT')
      OPEN(54,FILE='YPLUS.TXT')
      OPEN(55,FILE='ZPLUS.TXT')
      OPEN(56,FILE='YSTAR.TXT')
      OPEN(57,FILE='CMU.TXT')      
      OPEN(58,FILE='S.TXT')      
      OPEN(59,FILE='OMG.TXT')                  
      OPEN(60,FILE='rt.TXT')                  
      OPEN(61,FILE='ETA.TXT') 
      OPEN(62,FILE='xl3.TXT')                        
      AU=0.016
      AD=0.263
      ETA=-1.
      RATIO=0.5
      NI=NX
      NJ=NY
      NK=NZ
      READ(13,*)NI,UBL
      READ(13,*)LSCN,LSCS,LSCE,LSCW
      READ(13,*)LSTPN,LSTPS,LSTPE,LSTPW
*      UBL=UBL*0.9
      lstpn(1)=lstpn(1)+1
      lstps(1)=lstps(1)+1
      lstpn(3)=lstpn(3)+1
      lstps(3)=lstps(3)+1
C      NIRMS=NI-3
C      NIRMF=NI-1
      READ(14,*)NJ
C
      NJR=1
      NJRN=NJ-1
      I180=1
      I90 =0
      NIM1=NI-1
      NJM1=NJ-1
      NKM1=NK-1
C      NIR=10
      NKR=62
      NKRH=32
C------------- TURBULENCE CLOSURE SELECTION ---------------------------
C--------  IEVM=1  EFFECTIVE VISCOSITY K-EPSILON CLOSURE              *
C--------  IEVM=0  SECOND-MOMENT DSM CLOSURE                          *
C----------------------------------------------------------------------
      IEVM=1
C-------------NEAR WALL MODEL SELECTION -------------------------------
C--------  INWM=1  ONE-EQUATION MODEL                                 *
C--------  INWM=2  TWO-EQUATION MODEL                                 *
C----------------------------------------------------------------------
C*
      INWM=2
C*
C*-------------K-EPSILON MODEL SELECTION ------------------------------
C*-------		INLEVM=0	LINEAR K-EPSILON MODEL
C*-------		INLEVM=1	NON-LINEAR K-EPSILON MODEL
C*---------------------------------------------------------------------
	INLEVM=1
C*---------------------------------------------------------------------
C*-------		ICRAFT=0	SUGA FORMULA
C*-------		ICRAFT=1	CRAFT FORMULA
	ICRAFT=1
C*---------------------------------------------------------------------
C*-------		IYAP=0	YAP CORRECTION TERM
C*-------		IYAP=1	NEW YAP CORRECTION TERM
	IYAP=1
C-------- LW-RE NEAR-WALL ZONES   --------------------------------------
C------  INTEGER VARIABLES :        LWN, LWS, LWT AND LWB               *
C-------- REFER TO : NORTH, SOUTH, TOP AND BOTTOM BOUNDARIES RESPECTIVELY
C-------  TO CREATE A LOW-RE NEAR WALL ZONE SET RELEVANT VARIABLE TO 1  *
C------- OTHERWISE SET ABOVE INTEGER VARIABLES TO ZERO                  *
C------------------------------------------------------------------------
      LWN=1
      LWS=1
      LWT=1
C*      LWB=1
C*
	LWB=0
C*
C--------- NUMBER OF NODES IN THE NEAR-WALL ZONES -----------------------
      INWP=15
C-------------------------------------------------------------------------
C---- CHOOSE WHETHER TO RESTART EXISTING SOLUTION (NREAD=1)   ------------
C------- OR  WHETHER TO START FROM AN INITIAL FLOWFIELD (NREAD=0) --------
      NREAD=0
C-------------------------------------------------------------------------
      NJSL=INWP+1
      NJSL1=NJSL+1
      NJNH=NJ-INWP
      NJNH1=NJNH-1
      NKBL=INWP+1
      NKBL1=NKBL+1
      NKTH=NK-INWP
      NKTH1=NKTH-1
C     WRITE(*,*)NJSL,NJNH,NKBL,NKTH
C------------------------------------------------
      RESORM=0.
      RESORU=0.
      RESORV=0.
      RESORW=0.
      RESORK=0.
      RESORE=0.
C
      INCALU=.TRUE.
      INCALV=.TRUE.
      INCALW=.TRUE.
      INCALP=.TRUE.
      INCALK=.TRUE.
      INCALE=.TRUE.
      INPROP=.TRUE.
C
      CALL GRID1
      CALL GRID2
C
      TIME=0.
      CALL INIT
C
CMV     GO TO 112
      IF(NREAD.EQ.1)GO TO 112
      DO K=1,NK
      DO J=1,NJ
      DO I=1,NI
      VIS(I,J,K)=VISCOS
      ENDDO
      ENDDO
      ENDDO
      DO K=1,NKM1
      DO J=2,NJM1
      DO I=1,NI
      U(I,J,K)=UBL
      V(I,J,K)=0
      W(I,J,K)=0
      TKE(I,J,K)=0.02*UBL*UBL
      TED(I,J,K)=DEN(I,J,K)*0.09*(TKE(I,J,K)**2)/(10.*VISCOS)
      VIS(I,J,K)=11.*VISCOS
      UU(I,J,K)=2.*TKE(I,J,K)/3.
      VV(I,J,K)=2.*TKE(I,J,K)/3.
      WW(I,J,K)=2.*TKE(I,J,K)/3.
      UV(I,J,K)=0.
      UW(I,J,K)=0.
      VW(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO

  112 CONTINUE

       DO K=1,NKM1
       DO J=2,NJM1
       U(1,J,K)=UBL
       ENDDO
       ENDDO
       UIN=UBL
      CAREA=0.
      FLOWIN=0.
      XMONIN=0.
      DO 5 K=2,NKM1
      DO 5 J=2,NJ-1
C      DO 5 K=NKR+1,NKM1
      X1=X(4,J-1,K-1)
      X2=X(4,J,K-1)
      X5=X(4,J-1,K)
      X6=X(4,J,K)
C
      Y1=Y(4,J-1,K-1)
      Y2=Y(4,J,K-1)
      Y5=Y(4,J-1,K)
      Y6=Y(4,J,K)
C
      Z1=Z(4,J-1,K-1)
      Z2=Z(4,J,K-1)
      Z5=Z(4,J-1,K)
      Z6=Z(4,J,K)
C
      CALL FLORAT(U(1,J,K) ,0. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,FLOW)
      CALL FLORAT(1.,0. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,AREA)
       CAREA=CAREA+AREA
      FLOWIN=FLOWIN+DEN(1,J,K)*FLOW
      CALL FLORAT(UBL**2,0. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,XMON)
      XMONIN=XMONIN+DEN(1,J,K)*XMON
 5    CONTINUE
C

      VELM=FLOWIN/CAREA
      REY = VELM *(Y(2,NJ-1,2)-Y(2,1,2))/ VISCOS
      DELT=10000.*(Y(2,NJ-1,2)-Y(2,1,2))/VELM
      
      PRINT *,' FLOWIN, XMONIN, UBLK, REYNOD NO,   DELT, CROSS-AREA'
      WRITE(*,*) FLOWIN, XMONIN, VELM, REY, DELT, CAREA
      IF(NREAD.EQ.0)THEN
      DO 6002 K=1,NK
      DO 6002 J=1,NJ
      DO 6002 I=1,NI
 6002 T(I,J,K)=300.
      ENDIF
C
c      call props
      CALL ADVECT(1.,1.)
      if(nread.eq.0)then
      NSWPP1=NSWPP
      NSWPP=0
      CALL CALCP
      NSWPP=NSWPP1
      endif
C
      DO K=1,NK
      DO J=1,NJ
      DO I=1,NI
      U0(I,J,K)=U(I,J,K)
      V0(I,J,K)=V(I,J,K)
      W0(I,J,K)=W(I,J,K)
      T0(I,J,K)=T(I,J,K)
      TKE0(I,J,K)=TKE(I,J,K)
      TED0(I,J,K)=TED(I,J,K)
      UU0(I,J,K)=UU(I,J,K)
      VV0(I,J,K)=VV(I,J,K)
      WW0(I,J,K)=WW(I,J,K)
      UV0(I,J,K)=UV(I,J,K)
      UW0(I,J,K)=UW(I,J,K)
      VW0(I,J,K)=VW(I,J,K)
      ENDDO
      ENDDO
      ENDDO
C**********
      DO 20 IT=1,MAXIT
      CALL STRAIN
C
C      CALL STRAIN
       KK=0
      IF(INCALU) CALL CALCU
      IF(INCALV) CALL CALCV
      IF(INCALW) CALL CALCW
      IF(INCALP) CALL CALCP
      CALL STRAIN
c      call modtke
c      call calcuu
c      call calcvv
c      call calcww
c      call calcuv
c      call calcuw
c      call calcvw
c      DO  K=1,NK
c      DO  J=1,NJ
c      DO  I=1,NI
c      TKE(I,J,K)=.5*(UU(I,J,K)+VV(I,J,K)+WW(I,J,K))
c      enddo
c      enddo
c      enddo
      IF(IEVM.EQ.1.AND.INCALK) CALL CALCKE
C------------  DSM PART ------------------------------------------------
      IF(IEVM.EQ.0)THEN
      CALL MODTKE
      CALL CALCUU
      CALL CALCWW
      CALL CALCVV
      CALL CALCUV
      CALL CALCUW
      CALL CALCVW

      DO 30 K=1,NK
      DO 30 J=1,NJ
      DO 30 I=1,NI
   30 TKE(I,J,K)=.5*(UU(I,J,K)+VV(I,J,K)+WW(I,J,K))
      WRITE(11,*)(((TKE(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      ENDIF
C-----------------------------------------------------------------------
      IF(INCALE) CALL CALCED
      IF(INPROP) CALL PROPS
C
*      N1=0
*      DO 821 N=1,1
*      RESORT=RESORT/(FLOWIN*10.)
*      N1=N1+1
*      WRITE(6,* ) N1,RESORT
  821 CALL CALCT
*  821  CONTINUE 
C
      RESORM=RESORM/FLOWIN
      RESORU=RESORU/XMONIN
      RESORV=RESORV/XMONIN
      RESORW=RESORW/XMONIN
      RESORT=RESORT/(FLOWIN*10.)
C      RESORU=RESORU/(FLOWIN*TKEIN)
C      RESORV=RESORV/(FLOWIN*TKEIN)
C      RESORW=RESORW/(FLOWIN*TKEIN)
      RESORK=RESORK/(FLOWIN*TKEIN)
      RESORE=RESORE/(FLOWIN*TEDIN)
C
      SORCE=AMAX1(RESORM,RESORU,RESORV,RESORW)
C*
C      NAN_CHK=ISNAN(RESORE)
C      IF (NAN_CHK) THEN
C      WRITE(*,*) "EPSILON HAS BEEN NAN."
C      STOP
C      ENDIF
      IF(MOD(IT,1).EQ.0)
     1WRITE(52,* ) IT+FINAL_IT,RESORM,RESORU,RESORV,RESORW,RESORT,
     1               RESORK,RESORE

      IF(MOD(IT,1).EQ.0)
     1WRITE(*,* ) IT+FINAL_IT,RESORM,RESORU,RESORV,RESORW,RESORT,
     1               RESORK,RESORE

 89   FORMAT(I5,6F10.4)
C*      
	IF(SORCE.LE.SORMAX) GO TO 1000
C*
C**********	EMERGENCY EXIT 
C      PRESSED=PEEKCHARQQ()
C      IF (PRESSED) THEN 
C      KEY=GETCHARQQ()
C      IF(KEY.EQ.'q' .OR. KEY.EQ.'Q') THEN 
C      IT_NO=IT+1
C      GOTO 1000	  
C      ENDIF
C      PRESSED=.FALSE.
C      GOTO 20
C      ENDIF
C**********
   20 CONTINUE
C      IT_NO=IT   
 1000 CONTINUE
      NJH=NJ/2
      UIN2=UIN*UIN
      CALL MODTBL

C      DO 1001 K=2,NKM1
C      WRITE(*,*)K,U(3,NJH,K),U(3,K,2),P(3,NI+1-K,2)
C 1001 CONTINUE
C
CCC
C
      REWIND 7
      REWIND 8
      REWIND 9
      REWIND 10
      REWIND 11
      REWIND 12
      REWIND 21
      REWIND 22
      REWIND 23
      REWIND 24
      REWIND 30
      REWIND 31
      REWIND 32
      REWIND 33
      REWIND 34
      REWIND 35
      REWIND 39
C
c      GO TO 5002
      WRITE( 7,1999)(((U(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE( 8,*)(((V(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE( 9,*)(((W(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(10,*)DPEX,(((P(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(11,*)(((TKE(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(12,*)(((TED(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(21,*)(((CX(I,J,K),I=1,NIM1),J=1,NJM1),K=1,NKM1)
      WRITE(22,*)(((CY(I,J,K),I=1,NIM1),J=1,NJM1),K=1,NKM1)
      WRITE(23,*)(((CZ(I,J,K),I=1,NIM1),J=1,NJM1),K=1,NKM1)
      WRITE(24,*)(((VIS(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(25,999)(((YC(I,J,K),I=1,NI),J=1,NJ),K=1,NK)
      WRITE(26,999)(((ZC(I,J,K),I=1,NI),J=1,NJ),K=1,NK)
      WRITE(27,999)(((XC(I,J,K),I=1,NI),J=1,NJ),K=1,NK)

C      IF (INLEVM.EQ.0) CMU1=0.09
C      write(57,*) CMU1      
      WRITE(58,*)STLD
      WRITE(59,*)OMGTLD
   
      write(60,*) RTL    
      IF (INLEVM.EQ.0) THEN
 	DO 3021 I=2,NI-1
	DO 3021 J=2,NJ-1
	DO 3021 K=2,NK-1
      S11=2.0*DUDXC(I,J,K)
      S12=DUDYC(I,J,K)+DVDXC(I,J,K)
      S13=DUDZC(I,J,K)+DWDXC(I,J,K)
      S22=2.0*DVDYC(I,J,K)
      S23=DVDZC(I,J,K)+DWDYC(I,J,K)
      S33=2.0*DWDZC(I,J,K)
      OMG12=DUDYC(I,J,K)-DVDXC(I,J,K)
      OMG13=DUDZC(I,J,K)-DWDXC(I,J,K)
      OMG23=DVDZC(I,J,K)-DWDYC(I,J,K)
      ST2=S11**2+S22**2+S33**2+2.0*S12**2+2.0*S13**2+2.0*S23**2
      OMGT2=2.0*OMG12**2+2.0*OMG13**2+2.0*OMG23**2
	STLDP=TKE(I,J,K)*SQRT(0.5*ST2  )/(TED(I,J,K)+tiny)
	OMGTLDP=TKE(I,J,K)*SQRT(0.5*OMGT2)/(TED(I,J,K)+tiny)
	ETAL(I,J,K)=AMAX1(STLDP,OMGTLDP)
 3021 CONTINUE	      
      ENDIF
      WRITE(61,*) ETAL 
          
      
C      WRITE(28,1999)((UU(3,J,K),J=1,NJ),K=1,NK)
C      WRITE(28,1999)((VV(3,J,K),J=1,NJ),K=1,NK)
C      WRITE(28,1999)((WW(3,J,K),J=1,NJ),K=1,NK)
C      WRITE(29,1999)((UV(3,J,K),J=1,NJ),K=1,NK)
C      WRITE(29,1999)((UW(3,J,K),J=1,NJ),K=1,NK)
C      WRITE(29,1999)((VW(3,J,K),J=1,NJ),K=1,NK)
      WRITE(30,*)(((UU(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(31,*)(((VV(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(32,*)(((WW(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(33,*)(((UV(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(34,*)(((UW(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(35,*)(((VW(I,J,K),I=1,NI  ),J=1,NJ  ),K=1,NK  )
      WRITE(36,1999)((U(3,J,K),J=1,NJ),K=1,NK)
      WRITE(37,1999)((W(3,J,K),J=1,NJ),K=1,NK)
      WRITE(38,1999)((V(3,J,K),J=1,NJ),K=1,NK)
C
      WRITE(53,*)U,V,W 
 5002 WRITE(39,*)QWOCP,PRLTM,TTBLK,FLOWEN,T,QS,QN
 999  FORMAT(5F14.5)
 1999 FORMAT(5F18.5)
CCC
C  ----  GRAPHICAL OUTPUT --------------------------
C
      REWIND 40 
      REWIND 41
      NIG=4
      IG(1)=1
      IG(2)=3
      IG(3)=10
      IG(4)=21
      NJG=4
      NKG=4
      KG(1)=10
      KG(2)=27
      KG(3)=28
      KG(4)=32
      JG(1)=5
      JG(2)=25
      JG(3)=35
      JG(4)=57
      DO 2999 N=1,NJG
      DO 2999 I=1,NI
      DO 2999 K=1,NK
      VXZ(I,N,K)=V(I,JG(N),K)
      WXZ(I,N,K)=W(I,JG(N),K)
      VVXZ(I,N,K)=VV(I,JG(N),K)
      UVXZ(I,N,K)=UV(I,JG(N),K)
      UUXZ(I,N,K)=UU(I,JG(N),K)
 2999 UXZ(I,N,K)=U(I,JG(N),K)
      DO 3000 N=1,NIG
      DO 3000 J=1,NJ
      DO 3000 K=1,NK
      VYZ(N,J,K)=V(IG(N),J,K)
      WYZ(N,J,K)=W(IG(N),J,K)
 3000 UYZ(N,J,K)=U(IG(N),J,K)
      DO 3001 N=1,NKG
      DO 3001 I=1,NI
      DO 3001 J=1,NJ
      VXY(I,J,N)=V(I,J,KG(N))
      WXY(I,J,N)=W(I,J,KG(N))
      VVXY(I,J,N)=VV(I,J,KG(N))
      UVXY(I,J,N)=UW(I,J,KG(N))
      UUXY(I,J,N)=UU(I,J,KG(N))
 3001 UXY(I,J,N)=U(I,J,KG(N))
      WRITE(40,*)IG,UYZ,VYZ,WYZ
      WRITE(41,*)KG,UXY,VXY,WXY,UUXY,VVXY,UVXY
      WRITE(401,*)JG,UXZ,VXZ,WXZ,UUXZ,VVXZ,UVXZ
      do 3002 i=2,ni-1
      YSW=YC(I,2,2)-Y(I,1,2)
      YSTRS=(SQRT(TKE(I,2,2)))*xl3(i,2,2)/VISCOS
      YSTRN=(SQRT(TKE(I,NJ-1,2)))*xl3(i,NJ-1,2)/VISCOS
      WRITE(56,*) I,YSTRS,YSTRN
 3002 CONTINUE

C*
      OPEN(68,FILE='UU')
       DO J=1,NZ
       DO I=1,NY
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(68,*) UU(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(68)
      OPEN(69,FILE='UV')     
       DO J=1,NZ
       DO I=1,NY
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(69,*) UV(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(69)
       OPEN(70,FILE='VV')
       DO J=1,NZ
       DO I=1,NY
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(70,*) VV(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(70)
       OPEN(70,FILE='WW')
       DO J=1,NZ
       DO I=1,NY
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(70,*) WW(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(70)
      OPEN(71,FILE='U')   
       DO J=1,NZ
       DO I=1,NY
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(71,*) U(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(71)
      OPEN(72,FILE='V')     
       DO J=1,NZ
       DO I=1,NY
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(72,*) V(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(72)
       OPEN(73,FILE='W')
       DO J=1,NZ
       DO I=1,NY
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(73,*) W(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(73)

       WRITE(77,*) "VARIABLES=X,Y,Z,U"
       WRITE(77,*) 'ZONE T="1"'
       WRITE(77,*) "I=",NX,"J=",NZ,"K=",NY 
       DO I=1,NY
       DO J=1,NZ
       DO K=1,NX
       KK=J
       JJ=I
       II=K
       WRITE(77,*) X(II,JJ,KK),Y(II,JJ,KK),Z(II,JJ,KK),U(II,JJ,KK)
       ENDDO
       ENDDO
       ENDDO
       CLOSE(77)

	   
       DO I=2,NY-1
       DO J=2,NZ-1
       DO K=2,NX-1
       KK=J
       JJ=I
       II=K
       YPLUS(II,JJ,KK)=0.
       ENDDO
       ENDDO
       ENDDO
	   
      DO K=2,NK-1
      DO I=2,NI-1
      TAUWX=VISCOS*(DUDYC(I,2,K)+DVDXC(I,2,K))
      TAUWZ=VISCOS*(DVDZC(I,2,K)+DWDYC(I,2,K))
      TAUW=SQRT(TAUWX*TAUWX+TAUWZ*TAUWZ)
      UTAU=SQRT(TAUW/DEN(I,2,K))
      YPLUSS=XL3(I,2,K)*UTAU*DENSIT/VISCOS

      TAUWXN=VISCOS*(DUDYC(I,NJ-1,K)+DVDXC(I,NJ-1,K))
      TAUWZN=VISCOS*(DVDZC(I,NJ-1,K)+DWDYC(I,NJ-1,K))
      TAUWN=SQRT(TAUWXN*TAUWXN+TAUWZN*TAUWZN)
      UTAUN=SQRT(TAUWN/DEN(I,NJ-1,K))
      YPLUSN=XL3(I,NJ-1,K)*UTAUN*DENSIT/VISCOS
      YPLUS(I,2,K)=YPLUSS
      YPLUS(I,NJ-1,K)=YPLUSN
      ENDDO
      ENDDO

      DO J=2,NJ-1
      DO I=2,NI-1
      TAUXZ=VISCOS*(DWDXC(I,J,NK-1)+DUDZC(I,J,NK-1))
      TAUYZ=VISCOS*(DVDZC(I,J,NK-1)+DWDYC(I,J,NK-1))
      TAUW=SQRT(TAUXZ**2+TAUYZ**2)
      UTAU=SQRT(TAUW/DENSIT)
      ZPLUS=XL3(I,J,NK-1)*UTAU*DENSIT/VISCOS
      YPLUS(I,J,NK-1)=ZPLUS
      ENDDO
      ENDDO
      OPEN(55,FILE='YPLUS.PLT')
      WRITE(55,*) "VARIABLES=X,Y,Z,YPLUS"
      WRITE(55,*) 'ZONE'
      WRITE(55,*) "I=",NX-2,"J=",NY-2,"K=",NZ-2
      DO J=2,NZ-1
      DO I=2,NY-1
      DO K=2,NX-1
      KK=J
      JJ=I
      II=K
      WRITE(55,*) X(II,JJ,KK),Y(II,JJ,KK),Z(II,JJ,KK),YPLUS(II,JJ,KK)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(55)
	  
      TBL(1)=T(1,1,1)
      DO I=2,NX
      TBL(I)=TBL(I-1)+(QWOCP*(Y(1,NY,1)-Y(1,1,1))*(sqrt((X(I,NY,1)-
	1    X(I-1,NY,1))**2+(Y(I,NY,1)-Y(I-1,NY,1))**2)+sqrt((X(I,1,1)-
	1    X(I-1,1,1))**2+(Y(I,1,1)-Y(I-1,1,1))**2))/(FLOWIN))	
      ENDDO
C      TBL(1)=T(1,1,1)
C      DO I=2,NX
C      SUMM=0
C      DO J=10,NY-10
C      DO K=1,NZ-5
C      SUMM=SUMM+T(I,J,K)
C      ENDDO
C      ENDDO
C      TBL(I)=SUMM/((NZ-5)*(NY-19))
C      ENDDO
      OPEN(85,FILE="TBL.PLT")
      DO I=1,NX
      WRITE(85,*)X(I,1,1),TBL(I)
      ENDDO
      CLOSE(85)
	  
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
      NU(I,J,K)=0.
      IF(J.EQ.1.)THEN
      NU(I,J,K)=(QWOCP*PRLTM*(Y(1,NY,1)-Y(1,1,1)))/(VISCOS*(T(I,J,K)
	1    -TBL(I)))
      NU(1,J,K)=0.
      ENDIF
      IF(J.EQ.NY)THEN
      NU(I,J,K)=(QWOCP*PRLTM*(Y(1,NY,1)-Y(1,1,1)))/(VISCOS*(T(I,J,K)
	1    -TBL(I)))
      NU(1,J,K)=0.
      ENDIF
      IF(K.EQ.NZ)THEN
      NU(I,J,K)=(QWOCP*PRLTM*(Y(1,NY,1)-Y(1,1,1)))/(VISCOS*(T(I,J,K)
	1    -TBL(I)))
C      NU(1,J,K)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO
	  
      OPEN(80,FILE='NU')
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
      WRITE(80,*)NU(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      CLOSE(80)
c-------- pressure drop --------------------------------------------
      UBL=FLOWIN/(DENSIT*0.5)
      NJH=NJ/2
      DPPR=P(NI-1,2,2)-P(3,2,2)
      DPSC=P(NI-1,NJM1,2)-P(3,NJM1,2)
      DPTP=P(NI-1,NJH,NKM1)-P(3,NJH,NKM1)
      HXX1=HX+XC(NI-1,2,2)
      HXX2=HX+XC(3,2,2)
      DPPR=DPPR-.5*OMZ*OMZ*DENSIT*(HXX1*HXX1-HXX2*HXX2)
      HXX1=HX+XC(NI-1,NJM1,2)
      HXX2=HX+XC(3,NJM1,2)
      DPSC=DPSC-.5*OMZ*OMZ*DENSIT*(HXX1*HXX1-HXX2*HXX2)
      HXX1=HX+XC(NI-1,NJH,NKM1)
      HXX2=HX+XC(3,NJH,NKM1)
      DPTP=DPTP-.5*OMZ*OMZ*DENSIT*(HXX1*HXX1-HXX2*HXX2)
      DNM=0.5*DENSIT*UBL*UBL*(XC(NIM1,NJH,NKM1)-XC(3,NJH,NKM1))
      DPPR=DPPR/DNM
      DPSC=DPSC/DNM
      DPTP=DPTP/DNM
      WRITE(*,*)OMZ,DPPR,DPSC,DPTP
C--------------------------------------------------------------------
      DO K=1,NK
      DO J=1,NJ
      UEX(J,K)=U(NI,J,K)
      VEX(J,K)=V(NI,J,K)
      WEX(J,K)=W(NI,J,K)
      PEX(J,K)=P(NI,J,K)
      TKEX(J,K)=TKE(NI,J,K)
      TEDX(J,K)=TED(NI,J,K)
      TEX(J,K)=T(NI,J,K)
      VISEX(J,K)=VIS(NI,J,K)
      UUEX(J,K)=UU(NI,J,K)
      VVEX(J,K)=VV(NI,J,K)
      WWEX(J,K)=WW(NI,J,K)
      UVEX(J,K)=UV(NI,J,K)
      UWEX(J,K)=UW(NI,J,K)
      VWEX(J,K)=VW(NI,J,K)
      cxex(j,k)=cx(4,j,k)
      ENDDO
      ENDDO
C
      WRITE(51,*)uex,vex,wex,pex,tex,tkex,tedx,visex,uuex,vvex,wwex,
     1           uvex,uwex,vwex,cxex
c----------------------------------------------------------------------
C*
C      CALL GETTIM(tmphour2, tmpminute2, tmpsecond2, tmphund2)
C      IF(tmphour2.LT.tmphour1) tmphour2=tmphour2+24
C      TIME2=TMPHOUR2*3600+TMPMINUTE2*60+TMPSECOND2
C      HOUR=FLOOR(ABS(TIME2-TIME1))/3600
C      MINUTE=FLOOR(ABS(TIME2-TIME1-HOUR*3600))/60
C      SECOND=FLOOR(ABS(TIME2-TIME1-HOUR*3600-MINUTE*60))
C      WRITE(*,1001) HOUR," HOURS, ",MINUTE," MINUTES AND ",SECOND,
C    1  " SECONDS FOR ",IT_NO-1," INTERATIONS"
C1001   FORMAT(I2,A7,I2,A13,I2,A13,I5,A12)
C*
      STOP
      END
C
      SUBROUTINE BCMOD
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
	
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C------------------------------------------------------------------
      ENTRY MODP
C---- EAST-WEST
      PCN=P(2,3,3)
      DO 10 J=1,NJ
      DO 10 K=1,NK
      HXX2=HX+(1.-FX(2,J,K))*XC(2,J,K)+FX(2,J,K)*XC(3,J,K)
      DEN2=(1.-FX(2,J,K))*DEN(2,J,K)+FX(2,J,K)*DEN(3,J,K)
      HXX1=HX+(1.-FX(3,J,K))*XC(3,J,K)+FX(3,J,K)*XC(4,J,K)
      DEN1=(1.-FX(3,J,K))*DEN(3,J,K)+FX(3,J,K)*DEN(4,J,K)
      DPXS=0.5*OMZ*OMZ*(DEN1*HXX1*HXX1-DEN2*HXX2*HXX2)
      AW(2   ,J,K)=0.
  10   CONTINUE
C---- NORTH-SOUTH
      DO 20 K=1,NK
      DO 20 I=1,NI
      AN(I,NJM1,K)=0.
      AS(I,2   ,K)=0.
      P(I,1,K)=P(I,2,K)
      P(I,NJ,K)=P(I,NJM1,K)
 20   CONTINUE
      DO 30 I=1,NI
      DO 30 J=1,NJ
      AT(I,J,NK-1)=0.
      AB(I,J,2   )=0.
      P(I,J,1)=P(I,J,2)-.0*(P(I,J,3)-P(I,J,2))
       P(I,J,NK)=P(I,J,NK-1)+0.*(P(I,J,NK-1)-P(I,J,NK-2))
 30   CONTINUE
C
      RETURN
C
C-----------------------------------------------------------------------
      ENTRY MODTT
      EXCHAT=(VISCOS/PRLTM)
C----- EAST-WEST
      NJH=NJ/2
      DELTIN=T(NI-2,NJH,2)-300.
      DTDX=2.*(Y(2,NJ,2)-Y(2,1,2))*QWOCP/FLOWIN
      SAREA= (0.5+0.5+1.)*(X(NI,5,5)-X(4,5,5))
     1     - (2.*(Y(1,1,5))**2)
     1     +  2.*0.5*(Y(1,1,5))
      DO 210 K=1,NK
      DO 210 J=1,NJ
      AE(NI-1,J,K)=0.
      T(NI,J,K)=T(NI-1,J,K)
  210 CONTINUE
C------ NORTH-SOUTH
      DO 220 K=2,NKM1
      DO 220 I=2,NIM1
      J=2
C
      RAREAS=0.5*SQRT(
     1      ( ((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1       -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))**2)
     1     +( ((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1       -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))**2)
     1     +( ((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1       -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))**2)   )
      SU(I,2,K)=SU(I,2,K)+QWOCP*RAREAS
C
      J=NJ-1
      RAREAN=0.5*SQRT(
     1      ( ((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1       -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))**2)
     1     +( ((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1       -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))**2)
     1     +( ((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1       -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))**2)   )
      SU(I,NJM1,K)=SU(I,NJM1,K)+QWOCP*RAREAN     
      EXCHTT=EXCHAT+.5*(VIS(I,2,K)-VISCOS)/PRTTM
      T(I,1,K)=T(I,2,K)+.5*QWOCP*VOL(I,2,K)/(EXCHTT*RAREAS)      
      EXCHTT=EXCHAT+.5*(VIS(I,NJ-1,K)-VISCOS)/PRTTM
      T(I,NJ,K)=T(I,NJ-1,K)+.5*QWOCP*VOL(I,NJ-1,K)/(EXCHTT*RAREAN)   
  220 CONTINUE
C
C
      DO 223 K=2,NKM1
      DO 223 I=2,NIM1
      AS(I,2,K)=0.
      AN(I,NJM1,K)=0.
  223 CONTINUE
C
C------ TOP-BOTTOM
      DO 230 J=1,NJ
      DO 230 I=1,NI
      K=2
       T(I,J,K-1)=T(I,J,K)
C
      K=1
c      RAREAB=0.5*SQRT(
c     1      ( ((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
c     1       -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))**2)
c     1     +( ((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
c     1       -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))**2)
c     1     +( ((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
c     1       -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))**2)  )
c      SU(I,J,2)=SU(I,J,2)+RAREAB*0.0
C
      K=NK-1
      RAREAT=0.5*SQRT(
     1      ( ((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1       -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))**2)
     1     +( ((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1       -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))**2)
     1     +( ((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1       -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))**2)  )
      SU(I,J,NK-1)=SU(I,J,NK-1)+RAREAT*QWOCP
      EXCHTT=EXCHAT+.5*(VIS(I,J,NK-1)-VISCOS)/PRTTM
      T(I,J,NK)=T(I,J,NK-1)+0.5*QWOCP*VOL(I,J,NK-1)/(EXCHTT*RAREAT)
C      T(I,J,1)=T(I,J,2)+.0*QWOCP*VOL(I,J,2)/(EXCHAT*RAREAB)
  230 CONTINUE
C
C
      DO 221 I=1,NI
      T(I,1,1)=T(I,1,2)
      T(I,NJ,1)=T(I,NJ,2)
      T(I,1,NK)=.5*(T(I,1,NKM1)+T(I,2,NK))
      T(I,NJ,NK)=0.5*(T(I,NJ,NKM1)+T(I,NJM1,NK))
  221 CONTINUE
C
      DO 233 J=2,NJM1
      DO 233 I=2,NIM1
      AB(I,J,2)=0.
      AT(I,J,NK-1)=0.
  233 CONTINUE
c
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=1,NKM1
      DO I=1,NI
      T(I,J+1,K)=TWFNN(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN  
      J=2
      DO K=1,NKM1
      DO I=1,NI
      T(I,J-1,K)=TWFNS(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN  
      K=NKM1
      DO J=2,NJM1
      DO I=1,NI
      T(I,J,K+1)=TWFNT(I,J)
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=1,NI
      T(I,J,K-1)=TWFNB(I,J)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      RETURN
C--------------------------------------------------------------------
      ENTRY MODTBL
      I=3
      JRN=1
      NJRN=NJM1
      DTDX=2.*(Y(2,NJ,2)-Y(2,1,2))*QWOCP/FLOWIN
      TTBLK=0.
      FLOWEN=0.
      DO 725 J=2,NJM1
c      DO 725 J=njr+1,njrn
      DO 725 K=2,NKM1
C
      X1=(X(I,J-1,K-1)+X(I-1,J-1,K-1))/2.
      X2=(X(I,J,K-1)+X(I-1,J,K-1))/2.
      X5=(X(I,J-1,K)+X(I-1,J-1,K))/2.
      X6=(X(I,J,K)+X(I-1,J,K))/2.
C
      Y1=(Y(I,J-1,K-1)+Y(I-1,J-1,K-1))/2.
      Y2=(Y(I,J,K-1)+Y(I-1,J,K-1))/2.
      Y5=(Y(I,J-1,K)+Y(I-1,J-1,K))/2.
      Y6=(Y(I,J,K)+Y(I-1,J,K))/2.
C
      Z1=(Z(I,J-1,K-1)+Z(I-1,J-1,K-1))/2.
      Z2=(Z(I,J,K-1)+Z(I-1,J,K-1))/2.
      Z5=(Z(I,J-1,K)+Z(I-1,J-1,K))/2.
      Z6=(Z(I,J,K)+Z(I-1,J,K))/2.
C
      TCR=T(3,J,K)+(XC(2,J,K)-XC(2,2,2))*DTDX
      VELTM=U(I,J,K)*TCR
      CALL FLORAT(VELTM,0. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,TBLK)
      CALL FLORAT(U(I,J,K),0. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,FLOW)
      TTBLK=TTBLK+TBLK
      FLOWEN=FLOWEN+FLOW
  725 CONTINUE
      TTBLK=TTBLK/FLOWEN
      RETURN
C-------------------------------------------------------------------------
      ENTRY MODU
C
C---- TOP-BOTTOM
C      DO 1002 K=1,NK
      DO 1002 I=1,NI
      DO 1002 J=1,NJ
      AB(I,J,2)=0.
      U(I,J,1)=U(I,J,2)
      U(I,J,NK)=0.
 1002 CONTINUE
C---- NORTH-SOUTH
      DO 36 K=2,NKM1
      DO 36 I=2,NIM1
      U(I,1,K)=0.
      U(I,NJ,K)=0.
   36 CONTINUE
c----  EAST-WEST 
      DO K=1,NKM1
      DO J=2,NJM1
c      SP(NIM1,J,K)=SP(NIM1,J,K)-MAX(SU(NIM1,J,K),0.)/
c      1                          (ABS(U(NIM1,J,K))+TINY)
c      SU(NIM1,J,K)=MIN(SU(NIM1,J,K),0.)
      AE(NIM1,J,K)=0.
C      U(1,J,K)=UIN
      ENDDO
      ENDDO
      OPEN(90,FILE="UVW1")
      DO I=1,NY
      DO J=1,NZ
      KK=J
      JJ=I
      READ(90,*) U(1,JJ,KK),V(1,JJ,KK),W(1,JJ,KK),TKE(1,JJ,KK),
	1    TED(1,JJ,KK),UU(1,JJ,KK),VV(1,JJ,KK),WW(1,JJ,KK),UV(1,JJ,KK)
	1    ,UW(1,JJ,KK),VW(1,JJ,KK),CMU1(1,JJ,KK),VIS(1,JJ,KK)
      ENDDO
      ENDDO
      CLOSE(90)	  
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPNX(I,K)
      SP(I,J,K)=SP(I,J,K)-SPUWN(I,K)
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPSX(I,K)
      SP(I,J,K)=SP(I,J,K)-SPUWS(I,K)
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPTX(I,J)
      SP(I,J,K)=SP(I,J,K)-SPUWT(I,J)
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPBX(I,J)
      SP(I,J,K)=SP(I,J,K)-SPUWB(I,J)
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      IF(IWFN.GE.1.AND.IEVM.EQ.0)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPNX(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPSX(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPTX(I,J)
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPBX(I,J)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      RETURN
C---------------------------------------------------------------------------
      ENTRY MODXU
C
C     EAST
C
C      IF(I180.EQ.1) THEN
      I=NI-1
      NJR=1
      NJRN=NJM1
      FLOWX =0.
      DENAX =0.
      NJH=NJ/2
      NKH=NK/2
      DO 25 K=2,NKM1
      DO 25 J=2,NJ-1
C
      X1=(X(I,J-1,K-1)+X(I-1,J-1,K-1))/2.
      X2=(X(I,J,K-1)+X(I-1,J,K-1))/2.
      X5=(X(I,J-1,K)+X(I-1,J-1,K))/2.
      X6=(X(I,J,K)+X(I-1,J,K))/2.
C
      Y1=(Y(I,J-1,K-1)+Y(I-1,J-1,K-1))/2.
      Y2=(Y(I,J,K-1)+Y(I-1,J,K-1))/2.
      Y5=(Y(I,J-1,K)+Y(I-1,J-1,K))/2.
      Y6=(Y(I,J,K)+Y(I-1,J,K))/2.
C
      Z1=(Z(I,J-1,K-1)+Z(I-1,J-1,K-1))/2.
      Z2=(Z(I,J,K-1)+Z(I-1,J,K-1))/2.
      Z5=(Z(I,J-1,K)+Z(I-1,J-1,K))/2.
      Z6=(Z(I,J,K)+Z(I-1,J,K))/2.
C
      CALL FLORAT(U(I,J,K),0. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,FLOW)
      FLOWX =FLOWX +DEN(I,J,K)*FLOW
      CALL FLORAT(1., 0. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,AREA)
      PRESE=P(I-1,J,K)*(1.-FX(I-1,J,K))+P(I,J,K)*FX(I-1,J,K)
      PRESW=P(I-2,J,K)*(1.-FX(I-2,J,K))+P(I-1,J,K)*FX(I-2,J,K)
c      DENAX =DENAX +DEN(I,J,K)*AREA*vol(i,j,k)*DU1(I,J,K)
c     1       /xx(i,j,k)
      DENAX =DENAX +DEN(I,J,K)*AREA*DU1(I,J,K)
c      DENAX =DENAX +DEN(I,J,K)*AREA*DU1(I,J,K)*(PRESW-PRESE)
c      DENAX =DENAX +DEN(I-2,J,K)*AREA
 25   CONTINUE
C
      PPD=(FLOWIN-FLOWX)/DENAX
      DPEX=DPEX-URFP*PPD
      RATM=(FLOWIN-FLOWX)/FLOWIN
      WRITE(*,*)FLOWIN,FLOWX,RATM
C
      I=NIM1
      DO 351 J=2,NJm1
      DO 351 K=2,NKm1
      PRESE=P(I-1,J,K)*(1.-FX(I-1,J,K))+P(I,J,K)*FX(I-1,J,K)
      PRESW=P(I-2,J,K)*(1.-FX(I-2,J,K))+P(I-1,J,K)*FX(I-2,J,K)
      dpexl=dpex*vol(ni-1,j,k)/xx(nim1,j,k)
c      dpexl=dpex*(PRESE-PRESW)
      HXX1=HX+(1.-FX(NI-1,J,K))*XC(NI-1,J,K)+FX(NI-1,J,K)*
     1     XC(NI,J,K)
      DEN1=(1.-FX(NI-1,J,K))*DEN(NI-1,J,K)+FX(NI-1,J,K)*
     1     DEN(NI,J,K)
      HXX2=HX+(1.-FX(NI-2,J,K))*XC(NI-2,J,K)+FX(NI-2,J,K)*
     1     XC(NI-1,J,K)
      DEN2=(1.-FX(NI-2,J,K))*DEN(NI-2,J,K)+FX(NI-2,J,K)*
     1     DEN(NI-1,J,K)
c      HXX2=(XC(4,J,K))
c      DEN2=(DEN(4,J,K))
C      DPXS=0.5*OMZ*OMZ*(DEN1*HXX1*HXX1-DEN2*HXX2*HXX2)
C      P(NI,J,K)=P(NI-1,J,K)+(P(4,J,K)-P(3,J,K))
C     1 +0.5*OMZ*OMZ*(DEN(NI,J,K)*(XC(NI,J,K)**2)-(XC(NI-1,J,K)**2)
C     1 -.5*((XC(4,J,K)**2)-(XC(3,J,K)**2)) )+0.*DPEX
C      P(NI,J,K)=p(NI-1,J,K)+DPEXL
C     1 +0.5*OMZ*OMZ*(DEN(NI,J,K)*(XC(NI,J,K)**2)-DEN2*
C     2  (xc(NI-1,J,K)**2) ) +DPEX
C      p(ni-1,j,k)=p(ni,j,k)
C      P(NI-1,J,K)=((1.-FX(NI-3,J,K))*P(NI-3,J,K)+FX(NI-3,J,K)*
C     1   P(NI-2,J,K)+DPXS+DPEX-(1.-FX(NI-2,J,K))*P(NI-2,J,K))
C     2       /FX(NI-2,J,K)
      P(NI,J,K)=( (1.-FX(NI-2,J,K))*P(NI-2,J,K)+FX(NI-2,J,K)*P(NI-1,J,K)
     1           -(1.0-FX(NI-1,J,K))*P(NI-1,J,K) + DPEX )/FX(NI-1,J,K)
  351 CONTINUE
C
      DO 35 J=2,NJ-1
      DO 35 K=2,NKM1
      PRESE=P(I-1,J,K)*(1.-FX(I-1,J,K))+(
     1      P(I,J,K)-PPD)*FX(I-1,J,K)
      PRESW=P(I-2,J,K)*(1.-FX(I-2,J,K))+P(I-1,J,K)*FX(I-2,J,K)
c      U(NIM1,J,K)=U(NIM1,J,K)+DU1(NI-1,J,K)*PPD*vol(nim1,j,k)/
c     1          xx(nim1,j,k)
      U(NIM1,J,K)=U(NIM1,J,K)+DU1(NI-1,J,K)*PPD     
C     1           * vol(nim1,j,k)/xx(nim1,j,k)
      U(ni,j,k)=u(ni-1,j,k)
c      U(NI,J,K)=U(NIM1,J,K)+DU1(NI-1,J,K)*PPD*(PRESW-PRESE)
c      U(NI-1,J,K)=U(NIM1,J,K)+DU1(NI-1,J,K)*PPD
C      U(2,J,K)=U(NI-2,J,K)-0.*PPD
C      U(2,J,K)=U(NI-2,J,K)
C      U(3,J,K)=U(NI-1,J,K)
c      u(ni,j,k)=u(ni-1,j,k)
C      U(NI,J,K)=U(4,J,K)
C      U(1,J,K)=0.5*(U(NI-3,J,K)+U(NI-2,J,K))
C      DU1(3,J,K)=DU1(NI-1,J,K)
C      DU1(2,J,K)=DU1(NI-2,J,K)
C      DU1(NI,J,K)=DU1(4,J,K)
C      P(NI,J,K)=P(NIM1,J,K)+0.5*(P(NIM1,J,K)-P(NI-2,J,K))
 35   CONTINUE
      DO  J=1,NJ
      DO  K=1,NK
c      DO  I=1,NIR
      U(nI,J,K)=min(U(nI,J,K),0.)
      U(nim1,j,k)=min(u(nim1,j,k),0.)
       enddo
       enddo
C      END IF
      RETURN
C---------------------------------------------------------------------------
      ENTRY MODV
C---- TOP-BOTTOM
C      DO 1004 K=1,NK
      DO 1004 I=1,NI
      DO 1004 J=1,NJ
      AB(I,J,2)=0.
      V(I,J,1)=V(I,J,2)
      V(I,J,NK)=0.
 1004 CONTINUE
C
C     NORTH
C
      IF(I90.EQ.1) THEN
      I=NIM1
      FLOWY =0.
      DENAY =0.
      DO 26 J=2,NJM1
      DO 26 K=2,NKM1
C
      X1=(X(I,J-1,K-1)+X(I-1,J-1,K-1))/2.
      X2=(X(I,J,K-1)+X(I-1,J,K-1))/2.
      X5=(X(I,J-1,K)+X(I-1,J-1,K))/2.
      X6=(X(I,J,K)+X(I-1,J,K))/2.
C
      Y1=(Y(I,J-1,K-1)+Y(I-1,J-1,K-1))/2.
      Y2=(Y(I,J,K-1)+Y(I-1,J,K-1))/2.
      Y5=(Y(I,J-1,K)+Y(I-1,J-1,K))/2.
      Y6=(Y(I,J,K)+Y(I-1,J,K))/2.
C
      Z1=(Z(I,J-1,K-1)+Z(I-1,J-1,K-1))/2.
      Z2=(Z(I,J,K-1)+Z(I-1,J,K-1))/2.
      Z5=(Z(I,J-1,K)+Z(I-1,J-1,K))/2.
      Z6=(Z(I,J,K)+Z(I-1,J,K))/2.
C
      CALL FLORAT(0.,V(NIM1,J,K),0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,FLOW)
      FLOWY =FLOWY +DEN(NIM1,J,K)*FLOW
      CALL FLORAT(0., 1. ,0. ,X1,Y1,Z1,
     1 X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,AREA)
      DENAY =DENAY +DEN(NIM1,J,K)*AREA
 26   CONTINUE
C
      VINC=(FLOWIN-FLOWY)/DENAY
C
      WRITE(6,*) FLOWIN,FLOWY,DENAY
      END IF
C---- NORTH-SOUTH
      DO 1115 K=1,NK
      DO 1115 I=1,NI
      V(I,1,K)=0.0
      V(I,NJ,K)=0.
 1115 CONTINUE
C---- EAST-WEST
      DO 2004 J=2,NJ-1
      DO 2004 K=1,NK
      AE(NIM1,J,K)=0.
      V(NI,J,K)=V(NI-1,J,K)
 2004 CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPNY(I,K)
      SP(I,J,K)=SP(I,J,K)-SPUWN(I,K)
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPSY(I,K)
      SP(I,J,K)=SP(I,J,K)-1.*SPUWS(I,K)
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPTY(I,J)
      SP(I,J,K)=SP(I,J,K)-1.*SPUWT(I,J)
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPBY(I,J)
      SP(I,J,K)=SP(I,J,K)-SPUWT(I,J)
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      IF(IWFN.GE.1.AND.IEVM.EQ.0)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPNY(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPSY(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPTY(I,J)
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPBY(I,J)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      RETURN
C----------------------------------------------------------------------------
      ENTRY MODW
C---- TOP-BOTTOM
C      DO 2014 K=1,NK
      DO 2014 J=1,NJ
      DO 2014 I=1,NI
      W(I,J,NK)=0.
      W(I,J,1)=0.
 2014 CONTINUE
C---- NORTH-SOUTH
      DO 1116 K=2,NKM1
      DO 1116 I=2,NIM1
C      AN(I,NJM1,K)=0
C      AS(I,2,K)=0.
      W(I,NJ,K)=0.
      W(I,1,K)=0.
 1116 CONTINUE
C---- EAST-WEST
      DO 1117  K=2,NKM1
      DO 1117 J=2,NJM1
      AE(NIM1,J,K)=0.
      W(NI,J,K)=W(NI-1,J,K)
 1117 CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPNZ(I,K)
      SP(I,J,K)=SP(I,J,K)-SPUWN(I,K)
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPSZ(I,K)
      SP(I,J,K)=SP(I,J,K)-1.*SPUWS(I,K)
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPTZ(I,J)
      SP(I,J,K)=SP(I,J,K)-1.*SPUWT(I,J)
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)-TAWPBZ(I,J)
      SP(I,J,K)=SP(I,J,K)-SPUWB(I,J)
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      IF(IWFN.GE.1.AND.IEVM.EQ.0)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPNZ(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPSZ(I,K)
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPTZ(I,J)
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+PAWPBZ(I,J)
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      RETURN
C----------------------------------------------------------------------------
      RETURN
C----------------------------------------------------------------------------
      ENTRY MODTKE
C---- TOP
C      DO 70 K=1,NK
      DO 70 I=1,NI
      DO 70 J=1,NJ
      TKE(I,J,NK)=0.
 70   CONTINUE
C---- BOTTOM
      K=2
      DO 71 I=1,NI
      DO 71 J=1,NJ
      AB(I,J,K)=0.
      TKE(I,J,K-1)=TKE(I,J,K)
 71   CONTINUE
C---- NORTH
      J=NJM1
      DO 80 K=1,NK
      DO 80 I=1,NI
      TKE(I,J+1,K)=0.
C      AN(I,J,K)=0.
 80   CONTINUE
C---- SOUTH
      J=2
      DO 81 K=1,NK
      DO 81 I=1,NI
      TKE(I,J-1,K)=0.
C      AS(I,J,K)=0.
 81   CONTINUE
C---- EAST-WEST
C      I=NIM1
      DO 90 J=2,NJM1
      DO 90 K=2,NKM1
      AE(nI-1,J,K)=0.
      TKE(NI,J,K)=TKE(NI-1,J,K)
 90   CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NK-2
      RCTN(K)=1.
      ENDDO
C      RCTN(2)=0.
      RCTN(NKM1)=0.
      DO K=2,NKM1
      DO I=2,NIM1
      SU1=DEN(I,J,K)*(PK1(I,J,K)-TED(I,J,K))*VOL(I,J,K)
      SU(I,J,K)=SU(I,J,K)-MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)+MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      SU1=( (PKWN(I,K)-DEN(I,J,K)*TEDAVWN(I,K))*
     1         VOL(I,J,K) )
      SU(I,J,K)=SU(I,J,K)+MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NK-2
      RCTS(K)=1.
      ENDDO
C      RCTS(2)=0.
      RCTS(NKM1)=0.
      DO K=2,NKM1
      DO I=2,NIM1
      SU1=DEN(I,J,K)*(PK1(I,J,K)-TED(I,J,K))*VOL(I,J,K)
      SU(I,J,K)=SU(I,J,K)-MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)+MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      SU1=( (PKWS(I,K)-DEN(I,J,K)*TEDAVWS(I,K))*
     1       VOL(I,J,K) )
      SU(I,J,K)=SU(I,J,K)+MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      RCNT(J)=1.
      ENDDO
      IF(LWS.EQ.1)RCNT(2)=RZER
      IF(LWN.EQ.1)RCNT(NJM1)=RZER
      DO I=2,NIM1
      DO J=2,NJM1
      SU1=DEN(I,J,K)*(PK1(I,J,K)-TED(I,J,K))*VOL(I,J,K)
      SU(I,J,K)=SU(I,J,K)-RCNT(J)*MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)+RCNT(J)*MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      SU1=( (PKWT(I,J)-DEN(I,J,K)*TEDAVWT(I,J))*
     1     VOL(I,J,K) )
      SU(I,J,K)=SU(I,J,K)+MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      AT(I,J,K)=0.
      ENDDO
C      IF(IT.EQ.1)PRINT*,I,PKWT(I,2),PKWT(I,NJ/2),PKWT(I,NJM1)
C      IF(IT.EQ.1)PRINT*,I,TEDAVWT(I,2),TEDAVWT(I,NJ/2),TEDAVWT(I,NJM1)
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      RCNB(J)=1.
      ENDDO
      IF(LWS.EQ.1)RCNB(2)=RZER
      IF(LWN.EQ.1)RCNB(NJM1)=RZER
      DO J=2,NJM1
      DO I=2,NIM1
      SU1=DEN(I,J,K)*(PK1(I,J,K)-TED(I,J,K))*VOL(I,J,K)
      SU(I,J,K)=SU(I,J,K)-RCNB(J)*MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)+RCNB(J)*MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      SU1=(PKWB(I,J)-DEN(I,J,K)*TEDAVWB(I,J))*VOL(I,J,K)
      SU(I,J,K)=SU(I,J,K)+MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      IF(IEVM.EQ.0.AND.IWFN.GE.1.AND.IWFD.EQ.1)THEN
      JST=2
      IF(LWS.EQ.1)JST=3
      JFN=NJM1
      IF(LWN.EQ.1)JFN=NJ-2
      KST=2
      IF(LWB.EQ.1)KST=3
      KFN=NKM1
      IF(LWT.EQ.1)KFN=NK-2
      DO K=KST,KFN
      DO J=JST,JFN
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*TKE(I,J,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      RETURN
C----------------------------------------------------------------------------
      ENTRY MODTED
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
       DE=FLOAT(IEVM)*0.236+FLOAT(1-IEVM)*0.263
C---- TOP 1-EQN REGION
C      DO 72 K=NKED+1,NKM1
C      DO 72 J=2,NJM1
C      DO 72 I=2,NIM1
C      RT=DEN(I,J,K)*TKE(I,J,K)**2/(TED(I,J,K)+TINY)/VISCOS
C      F1=1.
C      F2=1.-0.3*EXP(-RT**2)
C      SU(I,J,K)= GREAT*F1*C1*PK2(I,J,K)/F2/C2
C      SP(I,J,K)=-GREAT
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      FMU=(1.-EXP(-0.0198*XLST))*(1.+5.29/XLST)
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-DE*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      TED(I,J,K)= F1*C1*PK2(I,J,K)/F2/C2
C   72 CONTINUE
C---- TOP-BOTTOM
      DO 73 K=1,NK
      DO 73 I=1,NI
      DO 73 J=1,NJ
      AB(I,J,2)=0.
      TED(I,J,1)=TED(I,J,2)
      TED(I,J,NK)=0.
 73   CONTINUE
C---- BOTTOM 1-EQN REGION
C      DO 84 K=2,KED-1
C      DO 84 J=2,NJM1
C      DO 84 I=2,NIM1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      FMU=(1.-EXP(-0.0198*XLST))*(1.+5.29/XLST)
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-DE*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      RT=DEN(I,J,K)*TKE(I,J,K)**2/(TED(I,J,K)+TINY)/VISCOS
C      F1=1.
C   84 CONTINUE
C---- NORTH 1-EQN REGION
C      DO 382 K=2,NKM1
C      DO 382 J=NJED+1,NJM1
C      DO 382 I=2,NIM1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      FMU=(1.-EXP(-0.0198*XLST))*(1.+5.29/XLST)
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-DE*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      RT=DEN(I,J,K)*TKE(I,J,K)**2/(TED(I,J,K)+TINY)/VISCOS
C      F1=1.
C      F2=1.-0.3*EXP(-RT**2)
C      SU(I,J,K)= GREAT*F1*C1*PK2(I,J,K)/F2/C2
C      SP(I,J,K)=-GREAT
C      TED(I,J,K)= F1*C1*PK2(I,J,K)/F2/C2
C      AN(I,J,K)=0.
C 382  CONTINUE
C---- SOUTH 1-EQN REGION
C      DO 83 K=2,NKM1
C      DO 83 J=2,JED-1
C      DO 83 I=2,NIM1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      FMU=(1.-EXP(-0.0198*XLST))*(1.+5.29/XLST)
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-DE*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      RT=DEN(I,J,K)*TKE(I,J,K)**2/(TED(I,J,K)+TINY)/VISCOS
C      F1=1.
C      F2=1.-0.3*EXP(-RT**2)
C      SU(I,J,K)= GREAT*F1*C1*PK2(I,J,K)/F2/C2
C      SP(I,J,K)=-GREAT
C      TED(I,J,K)= F1*C1*PK2(I,J,K)/F2/C2
C      AS(I,J,K)=0.
C 83   CONTINUE
C---- EAST-WEST
C      I=NIM1
      DO 92 J=2,NJM1
      DO 92 K=2,NKM1
      AE(nI-1,J,K)=0.
      TED(NI,J,K)=TED(ni-1,J,K)
c      TED(1,J,K)=TED(NI-3,J,K)
c      TED(2,J,K)=TED(NI-2,J,K)
c      TED(NI-1,J,K)=TED(3,J,K)
 92   CONTINUE
C---- NORTH-SOUTH 
      DO 1119 K=2,NKM1
      DO 1119 I=2,NIM1
C      AN(I,NJM1,K)=0.
C      AS(I,2,K)=0.
      TED(I,1,K)=0.
      TED(I,NJ,K)=0.
 1119 CONTINUE
C
C---- WALL-FUNCTIONS
      IF(IWFN.EQ.1.OR.IWFN.EQ.2.OR.IWFN.EQ.4)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* TEDWN(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* TEDWS(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* TEDWT(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* TEDWB(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      IF(IWFN.EQ.3)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*MAX(SUTEDAN(I,K),RZER)
      SP(I,J,K)=
     1  VOL(I,J,K)*MIN(SUTEDAN(I,K),RZER)/(TED(I,J,K)+TINY)
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*MAX(SUTEDAS(I,K),RZER)
      SP(I,J,K)=
     1      VOL(I,J,K)*MIN(SUTEDAS(I,K),RZER)/(TED(I,J,K)+TINY)
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      RCNT(J)=0.
      ENDDO
      IF(LWS.EQ.1)RCNT(2)=1.
      IF(LWN.EQ.1)RCNT(NJM1)=1.
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=RCNT(J)*SU(I,J,K)+
     1            VOL(I,J,K)*MAX(SUTEDAT(I,J),RZER)
      SP(I,J,K)=RCNT(J)*SP(I,J,K)+
     1       VOL(I,J,K)*MIN(SUTEDAT(I,J),RZER)/(TED(I,J,K)+TINY)
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+VOL(I,J,K)*MAX(SUTEDAB(I,J),RZER)
      SP(I,J,K)=SP(I,J,K)+
     1      VOL(I,J,K)*MIN(SUTEDAB(I,J),RZER)/(TED(I,J,K)+TINY)
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      RETURN
C--------------------------------------------------------------------------
      ENTRY MODTUU
C---- TOP
      K=NKM1
      DO 101 I=1,NI
      DO 101 J=1,NJ
      UU(I,J,NK)=0.
 101  CONTINUE
C---- BOTTOM
      K=2
      DO 102 I=1,NI
      DO 102 J=1,NJ
      AB(I,J,K)=0.
      UU(I,J,K-1)=UU(I,J,K)
 102  CONTINUE
C---- NORTH
      J=NJM1
      DO 103 K=1,NK
      DO 103 I=1,NI
      UU(I,J+1,K)=0.
C      AN(I,J,K)=0.
 103  CONTINUE
C---- SOUTH
      J=2
      DO 104 K=1,NK
      DO 104 I=1,NI
      UU(I,J-1,K)=0.
C      AS(I,J,K)=0.
 104  CONTINUE
C---- EAST-WEST
      I=NIM1
      DO 105 J=1,NJ
      DO 105 K=1,NK
      AE(I,J,K)=0.
      UU(NI,J,K)=UU(ni-1,J,K)
c      UU(1,J,K)=UU(NI-3,J,K)
c      UU(NI-1,J,K)=UU(3,J,K)
c      UU(2,J,K)=UU(NI-2,J,K)
 105  CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UUWFN(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* UUWFS(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UUWFT(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UUWFB(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
      RETURN
C-----------------------------------------------------------------------------
      ENTRY MODTVV
C---- TOP
      K=NKM1
      DO 106 I=1,NI
      DO 106 J=1,NJ
      VV(I,J,K+1)=0.
 106  CONTINUE
C---- BOTTOM
      K=2
      DO 107 I=1,NI
      DO 107 J=1,NJ
      AB(I,J,K)=0.
      VV(I,J,K-1)=VV(I,J,K)
 107  CONTINUE
C---- NORTH
      J=NJM1
      DO 108 K=1,NK
      DO 108 I=1,NI
      VV(I,J+1,K)=0.
C      AN(I,J,K)=0.
 108  CONTINUE
C---- SOUTH
      J=2
      DO 109 K=1,NK
      DO 109 I=1,NI
      VV(I,J-1,K)=0.
C      AS(I,J,K)=0.
 109  CONTINUE
C---- EAST-WEST
C     I=NIM1
      DO 110 J=1,NJ
      DO 110 K=1,NK
      VV(NI,J,K)=VV(ni-1,J,K)
c      VV(1,J,K)=VV(NI-3,J,K)
      AE(nI-1,J,K)=0.
c      VV(NI-1,J,K)=VV(3,J,K)
c      VV(2,J,K)=VV(NI-2,J,K)
 110  CONTINUE
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*VVWFN(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* VVWFS(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*VVWFT(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*VVWFB(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      RETURN
C------------------------------------------------------------------------------
      ENTRY MODTWW
C---- TOP
      K=NKM1
      DO 111 I=1,NI
      DO 111 J=1,NJ
      WW(I,J,K+1)=0.
 111  CONTINUE
C---- BOTTOM
      K=2
      DO 112 I=1,NI
      DO 112 J=1,NJ
      AB(I,J,K)=0.
      WW(I,J,K-1)=WW(I,J,K)
 112  CONTINUE
C---- NORTH
      J=NJM1
      DO 113 K=1,NK
      DO 113 I=1,NI
      WW(I,J+1,K)=0.
C      AN(I,J,K)=0.
 113  CONTINUE
C---- SOUTH
      J=2
      DO 114 K=1,NK
      DO 114 I=1,NI
      WW(I,J-1,K)=0.
C      AS(I,J,K)=0.
 114  CONTINUE
C---- EAST-WEST
      I=NIM1
      DO 115 J=1,NJ
      DO 115 K=1,NK
      WW(NI,J,K)=WW(ni-1,J,K)
c      WW(1,J,K)=WW(NI-3,J,K)
      AE(I,J,K)=0.
c      WW(NI-1,J,K)=WW(3,J,K)
c      WW(2,J,K)=WW(NI-2,J,K)
 115  CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*WWWFN(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* WWWFS(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*WWWFT(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*WWWFB(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      RETURN
C----------------------------------------------------------------------------
      ENTRY MODTUV
C---- TOP
      K=NKM1
      DO 116 I=1,NI
      DO 116 J=1,NJ
      UV(I,J,K+1)=0.
 116  CONTINUE
C---- BOTTOM
      K=2
      DO 117 I=1,NI
      DO 117 J=1,NJ
      AB(I,J,K)=0.
      UV(I,J,K-1)=UV(I,J,K)
 117  CONTINUE
C---- NORTH
      J=NJM1
      DO 118 K=1,NK
      DO 118 I=1,NI
      UV(I,J+1,K)=0.
 118  CONTINUE
C---- SOUTH
      J=2
      DO 119 K=1,NK
      DO 119 I=1,NI
      UV(I,J-1,K)=0.
 119  CONTINUE
C---- EAST
      I=NIM1
      DO 120 J=1,NJ
      DO 120 K=1,NK
      UV(NI,J,K)=UV(ni-1,J,K)
c      UV(1,J,K)=UV(NI-3,J,K)
      AE(I,J,K)=0.
c      UV(NI-1,J,K)=UV(3,J,K)
c      UV(2,J,K)=UV(NI-2,J,K)
 120  CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UVWFN(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* UVWFS(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UVWFT(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UVWFB(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      RETURN
C----------------------------------------------------------------------------
      ENTRY MODTUW
C---- TOP
      K=NKM1
      DO 121 I=1,NI
      DO 121 J=1,NJ
      UW(I,J,K+1)=0.
 121  CONTINUE
C---- BOTTOM
      K=2
      DO 122 I=1,NI
      DO 122 J=1,NJ
      UW(I,J,K-1)=0.
 122  CONTINUE
C---- NORTH
      J=NJM1
      DO 123 K=1,NK
      DO 123 I=1,NI
      UW(I,J+1,K)=0.
C      AN(I,J,K)=0.
 123  CONTINUE
C---- SOUTH
      J=2
      DO 124 K=1,NK
      DO 124 I=1,NI
      UW(I,J-1,K)=0.
C      AS(I,J,K)=0.
 124  CONTINUE
C---- EAST-WEST
      I=NIM1
      DO 125 J=1,NJ
      DO 125 K=1,NK
      UW(NI,J,K)=UW(ni-1,J,K)
c      UW(1,J,K)=UW(NI-3,J,K)
      AE(I,J,K)=0.
c      UW(NI-1,J,K)=UW(3,J,K)
c      UW(2,J,K)=UW(NI-2,J,K)
 125  CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UWWFN(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* UWWFS(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UWWFT(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*UWWFB(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      RETURN
C----------------------------------------------------------------------------
      ENTRY MODTVW
C---- TOP
      K=NKM1
      DO 126 I=1,NI
      DO 126 J=1,NJ
      VW(I,J,K+1)=0.
 126  CONTINUE
C---- BOTTOM
      K=2
      DO 127 I=1,NI
      DO 127 J=1,NJ
      VW(I,J,K-1)=0.
 127  CONTINUE
C---- NORTH
      J=NJM1
      DO 128 K=1,NK
      DO 128 I=1,NI
      VW(I,J+1,K)=0.
 128  CONTINUE
C---- SOUTH
      J=2
      DO 129 K=1,NK
      DO 129 I=1,NI
      VW(I,J-1,K)=0.
 129  CONTINUE
C---- EAST-WEST
      I=NIM1
      DO 130 K=1,NK
      DO 130 J=1,NJ
      VW(NI,J,K)=VW(ni-1,J,K)
c      VW(1,J,K)=VW(NI-3,J,K)
      AE(I,J,K)=0.
c      VW(NI-1,J,K)=VW(3,J,K)
c      VW(2,J,K)=VW(NI-2,J,K)
 130  CONTINUE
C
C---- WALL-FUNCTION
      IF(IWFN.GE.1)THEN
      IF(LWN.EQ.1)THEN
      J=NJM1
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*VWWFN(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AN(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWS.EQ.1)THEN
      J=2
      DO K=2,NKM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)* VWWFS(I,K)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWT.EQ.1)THEN
      K=NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*VWWFT(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AT(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      IF(LWB.EQ.1)THEN
      K=2
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=VOL(I,J,K)*VWWFB(I,J)*GREAT
      SP(I,J,K)=-VOL(I,J,K)*GREAT
      AB(I,J,K)=0.
      ENDDO
      ENDDO
      ENDIF
      ENDIF
C
      RETURN
      END
C----------------------------------------------------------------------------
      SUBROUTINE INIT
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10,NSG=67 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),DELD(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
c
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
C
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
C
      COMMON /WFNBL/ TAWPNX(NX,NZ),TAWPNY(NX,NZ),TAWPNZ(NX,NZ),
     1        PKWN(NX,NZ),TEDAVWN(NX,NZ),TEDWN(NX,NZ)
     2        ,TAWPSX(NX,NZ),TAWPSY(NX,NZ),TAWPSZ(NX,NZ),PKWS(NX,NZ),
     3        TEDAVWS(NX,NZ),TEDWS(NX,NZ)
     4        ,TAWPTX(NX,NY),TAWPTY(NX,NY),TAWPTZ(NX,NY),PKWT(NX,NY),
     5        TEDAVWT(NX,NY),TEDWT(NX,NY)
     4        ,TAWPBX(NX,NY),TAWPBY(NX,NY),TAWPBZ(NX,NY),PKWB(NX,NY),
     5        TEDAVWB(NX,NY),TEDWB(NX,NY)
     6        ,UUWFN(NX,NZ),VVWFN(NX,NZ),WWWFN(NX,NZ),UVWFN(NX,NZ),
     7        UWWFN(NX,NZ),VWWFN(NX,NZ)
     6        ,UUWFS(NX,NZ),VVWFS(NX,NZ),WWWFS(NX,NZ),UVWFS(NX,NZ),
     7        UWWFS(NX,NZ),VWWFS(NX,NZ)
     8        ,UUWFT(NX,NY),VVWFT(NX,NY),WWWFT(NX,NY),UVWFT(NX,NY),
     9        UWWFT(NX,NY),VWWFT(NX,NY)
     8        ,UUWFB(NX,NY),VVWFB(NX,NY),WWWFB(NX,NY),UVWFB(NX,NY),
     9        UWWFB(NX,NY),VWWFB(NX,NY)
     1    ,SUTEDAN(NX,NZ),SUTEDAS(NX,NZ),SUTEDAT(NX,NY),SUTEDAB(NX,NY)
C
      COMMON /WFWNSTR/ PAWPNX(NX,NZ),PAWPNY(NX,NZ),PAWPNZ(NX,NZ),
     2        PAWPSX(NX,NZ),PAWPSY(NX,NZ),PAWPSZ(NX,NZ),
     4        PAWPTX(NX,NY),PAWPTY(NX,NY),PAWPTZ(NX,NY),
     4        PAWPBX(NX,NY),PAWPBY(NX,NY),PAWPBZ(NX,NY)
     5       , SPUWN(NX,NZ),SPUWS(NX,NZ),SPUWT(NX,NY),SPUWB(NX,NY)
     6       ,TWFNN(NX,NZ),TWFNS(NX,NZ),TWFNT(NX,NY),TWFNB(NX,NY)
      COMMON /SGBLK/
     1   URSGN(NX,NZ,NSG),UNSGN(NX,NZ,NSG),UTSGN(NX,NZ,NSG),
     2   TKESGN(NX,NZ,NSG),TEDSGN(NX,NZ,NSG), VISSGN(NX,NZ,NSG),
     3   TSGN(NX,NZ,NSG),UR2SGN(NX,NZ,NSG),UN2SGN(NX,NZ,NSG),
     4   UT2SGN(NX,NZ,NSG),URUNSGN(NX,NZ,NSG),URUTSGN(NX,NZ,NSG),
     5   UNUTSGN(NX,NZ,NSG),
     1   URSGS(NX,NZ,NSG),UNSGS(NX,NZ,NSG),UTSGS(NX,NZ,NSG),
     2   TKESGS(NX,NZ,NSG),TEDSGS(NX,NZ,NSG), VISSGS(NX,NZ,NSG),
     3   TSGS(NX,NZ,NSG),UR2SGS(NX,NZ,NSG),UN2SGS(NX,NZ,NSG),
     4   UT2SGS(NX,NZ,NSG),URUNSGS(NX,NZ,NSG),URUTSGS(NX,NZ,NSG),
     5   UNUTSGS(NX,NZ,NSG),
     1   URSGT(NX,NY,NSG),UNSGT(NX,NY,NSG),UTSGT(NX,NY,NSG),
     2   TKESGT(NX,NY,NSG),TEDSGT(NX,NY,NSG), VISSGT(NX,NY,NSG),
     3   TSGT(NX,NY,NSG),UR2SGT(NX,NY,NSG),UN2SGT(NX,NY,NSG),
     4   UT2SGT(NX,NY,NSG),URUNSGT(NX,NY,NSG),URUTSGT(NX,NY,NSG),
     5   UNUTSGT(NX,NY,NSG),FRCFN(NX,NZ),FRCFS(NX,NZ)
       COMMON /SGBLK1/
     1   URSGB(NX,NY,NSG),UNSGB(NX,NY,NSG),UTSGB(NX,NY,NSG),
     2   TKESGB(NX,NY,NSG),TEDSGB(NX,NY,NSG), VISSGB(NX,NY,NSG),
     3   TSGB(NX,NY,NSG),UR2SGB(NX,NY,NSG),UN2SGB(NX,NY,NSG),
     4   UT2SGB(NX,NY,NSG),URUNSGB(NX,NY,NSG),URUTSGB(NX,NY,NSG),
     5   UNUTSGB(NX,NY,NSG),
     6   EDNWSG(NSG),PKSG(NSG),PURSG(NSG),PUNSG(NSG),PURNSG(NSG),
     7   PURTSG(NSG),PUNTSG(NSG),TRESG(NSG),TEL(NSG),ETSG(NSG),
     8   YAPSG(NSG),PUTSG(NSG)
       COMMON /SGBLK2/
     1   ANSG(NSG),ASSG(NSG),APSG(NSG),SUSG(NSG),SPSG(NSG),SURRSG(NSG),
     2   SUNNSG(NSG),SUTTSG(NSG),SURNSG(NSG),SURTSG(NSG),SUNTSG(NSG),
     3   APRRSG(NSG),APNNSG(NSG),APTTSG(NSG),APRNSG(NSG),APRTSG(NSG),
     4   APNTSG(NSG),SURRAV,SUNNAV,SUTTAV,SURNAV,SUNTAV,SURTAV,TAWPN,
     5   TAWPTN,TAWPS,TAWPTS,TAWPT,TAWPTT,TAWPB,TWAPTB
C
C
      DIMENSION QUIN(NY,NZ),QVIN(NY,NZ),QTKIN(NY,NZ),QTEIN(NY,NZ)
      DIMENSION ELN2(NX,NY),ELS2(NX,NY),COSNX2(NX,NY),COSNY2(NX,NY)
      DIMENSION COSSX2(NX,NY),COSSY2(NX,NY)
      DIMENSION ELNI(NX,NY), ELSI(NX,NY), ELEI(NX,NY), ELWI(NX,NY),
     1          ELTI(NX,NY),X2D(NX,NY),Y2D(NX,NY)
C
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
	 
      dimension uex(ny,NZ),vex(ny,NZ),wex(ny,NZ),tkex(ny,NZ),
     1          tedx(ny,NZ),uuex(ny,NZ),vvex(NY,NZ),wwex(ny,NZ),
     1          visex(ny,NZ),pex(ny,NZ),TEX(NY,NZ),uvex(ny,nz),
     1          uwex(ny,nz),vwex(ny,nz),cxex(ny,nz),YCO(NY,NZ)
     1          ,ZCO(NY,NZ)
C------------------------------------------------------------------------
C          SPECIFY CODE CONSTANTS  GREAT AND TINY                       *
C------------------------------------------------------------------------
      GREAT=1.E+10
      TINY=1.E-10
C------- SPECIFY NUMBER OF TDMA SWEEPS  ---------------------------------
      NSWPU=1
      NSWPV=1
      NSWPW=1
      NSWPP=5
      NSWPK=1
      NSWPE=1
      NSWPT=1
C
      ISCMU=2
      ISCMK=2
      ISCME=2
C------------------------------------------------------------------------
C      SPECIFY UNDER-RELAXATION FACTORS                                 *
C------------------------------------------------------------------------
      URFU=0.38
      URFV=0.38
      URFW=0.38
      URFP=0.38
      URFK=0.5
      URFE=0.5
      URFT=0.75
      URFVIS=0.5
C------ SPECIFY REFERENCE LOCATION  (IPREF,JPREF,KPREF)  ----------------
      IPREF=1
      JPREF=NJ/2
      KPREF=2
C-------------------------------------------------------------------------
C    FLUID PROPERTIES AND FLOW QUANTIDIES                                *
C-------------------------------------------------------------------------
      UIN=39.28
      SORMAX=0.0001
      RENO=100000.
      VISCOS=2.044/100000.
      DENSIT=1.02461
      PRLTM=0.7
      QWOC=100.
      RZER=0.
C-------------------------------------------------------------------------
C                  ROTATION  PARAMETERS                                  *
C    1. DISTANCES BETWEEN THE ORIGIN OF THE ROTATING CO-ORDINATE AND     *
C  THE CENTER OF ROTATION IN THE X,Y,AND Z DIRECTIONS:    HX,HY AND HZ   *
C    2. ROTATION VECTOR OF THE ROTATING FRAME OF REFERENCE: (OMX,OMY,OMZ)*
C-------------------------------------------------------------------------
      HX=-0.2225
      HY=0.
      HZ=0.
      OMX=0.
      OMY=0.
C     OMZ=(-0.2*UBL)/(Y(1,NY,1)-Y(1,1,1))
      OMZ=0.
C-------------------------------------------------------------------------
C     TURBULENCE   MODEL  CONSTANTS                                      *
C-------------------------------------------------------------------------
      C1=1.44
      C2=1.92
      CMU=0.09
      CAPPA=0.4187
      ELOG=9.793
      ELC=9.793
      ELCMU=ELC*CMU**0.25
      RKCMU=CAPPA*CMU**0.25
      PRTKE=1.
      PRTED=1.3
      PRTTM=0.9
      TKEIN=0.005*UIN**2
c      TEDIN=TKEIN**1.5/(0.11*1.34)
      TEDIN=CMU*DENSIT*TKEIN*TKEIN/(5.*VISCOS)
      CA1=1.8
      CA2=0.6
      CA1W=0.5
      CA2W=0.3
C-----------------------------------------------------------------------
C      INITIALISATION OF FLOW VARIABLES                                *
C-----------------------------------------------------------------------
      DO 10 I=1,NI
      DO 10 J=1,NJ
      DO 10 K=1,NK
      U(I,J,K)=0.
      V(I,J,K)=0.
      W(I,J,K)=0.
      TKE(I,J,K)=TKEIN
      TED(I,J,K)=TEDIN
      EDNW(I,J,K)=0.
      P(I,J,K)=0.
      PP(I,J,K)=0.
      DEN(I,J,K)=DENSIT
C
C     VIS(I,J,K)=VISCOS
C
      VIS(I,J,K)=VISCOS+CMU*DENSIT*TKEIN**2/(TEDIN+TINY)
C
      PK1(I,J,K)=0.
      PK2(I,J,K)=0.
      AE(I,J,K)=0.
      AW(I,J,K)=0.
      AN(I,J,K)=0.
      AS(I,J,K)=0.
      AT(I,J,K)=0.
      AB(I,J,K)=0.
      AP(I,J,K)=0.
      SU(I,J,K)=0.
      SP(I,J,K)=0.
      DU1(I,J,K)=0.
      DU2(I,J,K)=0.
      DU3(I,J,K)=0.
      DV1(I,J,K)=0.
      DV2(I,J,K)=0.
      DV3(I,J,K)=0.
      DW1(I,J,K)=0.
      DW2(I,J,K)=0.
      DW3(I,J,K)=0.
      CX(I,J,K)=0.
      CY(I,J,K)=0.
      CZ(I,J,K)=0.
      CXM(I,J,K)=0.
      CYM(I,J,K)=0.
      CZM(I,J,K)=0.
10    CONTINUE
C
      IF(NREAD.EQ.1)THEN
C
      PRINT*,'7'
      READ(7,*)TIME,U,V,W,DPEX,P,TKE,TED,CX,CY,CZ,VIS,UU,VV,WW,UV,UW,
     1              VW
C      GO TO 77
      IF(IWFN.EQ.3.AND.LWN.EQ.1.AND.NRDSG.EQ.1)
     1 READ(7,*)
     2   URSGN,UNSGN,UTSGN,TKESGN,TEDSGN,VISSGN,TSGN,UR2SGN,UN2SGN,
     3   UT2SGN,URUNSGN,URUTSGN,UNUTSGN
      IF(IWFN.EQ.3.AND.LWS.EQ.1.AND.NRDSG.EQ.1)
     1 READ(7,*)
     2   URSGS,UNSGS,UTSGS,TKESGS,TEDSGS,VISSGS,TSGS,UR2SGS,UN2SGS,
     3   UT2SGS,URUNSGS,URUTSGS,UNUTSGS
      IF(IWFN.EQ.3.AND.LWT.EQ.1.AND.NRDSG.EQ.1)
     1 READ(7,*)
     2   URSGT,UNSGT,UTSGT,TKESGT,TEDSGT,VISSGT,TSGT,UR2SGT,UN2SGT,
     3   UT2SGT,URUNSGT,URUTSGT,UNUTSGT
      IF(IWFN.EQ.3.AND.LWB.EQ.1.AND.NRDSG.EQ.1)
     1 READ(7,*)
     2   URSGB,UNSGB,UTSGB,TKESGB,TEDSGB,VISSGB,TSGB,UR2SGB,UN2SGB,
     3   UT2SGB,URUNSGB,URUTSGB,UNUTSGB
   77 CONTINUE
C
      PRINT*,'42'
      READ(42,*)ELN2,ELS2,COSNX2,COSNY2,COSSX2,COSSY2
C
      READ(39,*)QWOCP,PRLTM,TTBLK,FLOWEN,VISCOS,T
      QWOCP=1000.
c      READ(47,*)YCO,ZCO
C
      ENDIF
c--------------------------------------------------------------------
      QWOCP=1000.
c---------------------------------------------------------------------
      NJ2=NJ/2
      NK2=NK/2
      HXX1=HX+(1.-FX(NI-2,NJ2,NK2))*XC(NI,NJ2,NK2)+
     1         FX(NI-2,NJ2,NK2)*XC(NI-1,NJ2,NK2)
      DEN1=(1.-FX(NI-2,NJ2,NK2))*DEN(NI-2,NJ2,NK2)+
     1     FX(NI-2,NJ2,NK2)*DEN(NI-1,NJ2,NK2)
      DELX=(1.-FX(NI-2,NJ2,NK2))*XC(NI-2,NJ2,NK2)+
     1     FX(NI-2,NJ2,NK2)*XC(NI-1,NJ2,NK2)
     1     -(1.-FX(NI-3,NJ2,NK2))*XC(NI-3,NJ2,NK2)-
     1     FX(NI-3,NJ2,NK2)*XC(NI-2,NJ2,NK2)
      HXX2=HX+(1.-FX(NI-3,NJ2,NK2))*XC(NI-3,NJ2,NK2)+
     1     FX(NI-3,NJ2,NK2)*XC(NI-2,NJ2,NK2)
      DEN2=(1.-FX(NI-3,NJ2,NK2))*DEN(NI-3,NJ2,NK2)+
     1     FX(NI-3,NJ2,NK2)*DEN(NI-2,NJ2,NK2)
      DPXS=DEN1*OMZ*OMZ*HXX*DELX
      DPXS=.5*OMZ*OMZ*(DEN1*HXX1*HXX1-DEN2*HXX2*HXX2)
 1999 FORMAT(5F18.5)
 999  FORMAT(5F14.5)          
      print*,"1"
C     GO  TO  118

      READ(51,*)uex,vex,wex,tex,tkex,tedx,visex
      1          ,UUEX,VVEX,WWEX,UVEX,UWEX,VWEX

      DO J=1,NJ
      VISEX(J,1)=VISEX(J,2)
      UEX(J,1)=UEX(J,2)
      VEX(J,1)=VEX(J,2)
      TEX(J,1)=TEX(J,2)
      TKEX(J,1)=TKEX(J,2)
      TEDX(J,1)=TEDX(J,2)
      ENDDO
      do k=1,nk
      do j=1,nj
c      cx(1,j,k)=cxex(j,k)
      tke(1,j,k)=tkex(j,k)
      ted(1,j,k)=tedx(j,k)
      vis(1,j,k)=visex(j,k)
      uu(1,j,k)=UUEX(j,k)
      vv(1,j,k)=VVex(j,k)
      ww(1,j,k)=WWex(j,k)
      uv(1,j,k)=UVex(j,k)
      uw(1,j,k)=UWex(j,k)
      vw(1,j,k)=VWex(j,k)
      t(1,j,k)=tex(j,k)
      u(1,j,k)=uex(j,k)
      v(1,j,k)=vex(j,k)
      w(1,j,k)=wex(j,k)
C      p(1,j,k)=pex(j,k)
      enddo
c
      enddo

  324 continue
C
  118 CONTINUE

      IF(NREAD.EQ.0)THEN
C
      print*, '2'
C
      print*, '3'
      WRITE(*,*)NI,NJ,NK
       dth=3.142/31.
       DO 20 K=1,NK
       DO 20 J=1,NJ
       DO 20 I=1,NI
C       I=1
C
      u(i,j,nk)=0.
      u(i,1,k)=0.
      u(i,nj,k)=0.
C      v(i,j,k)=vex(j,k)
      w(i,j,k)=wex(j,k)
C      p(i,j,k)=pex(j,k)
      tke(i,j,k)=tkex(j,k)
      ted(i,j,k)=tedx(j,k)
      vis(i,j,k)=visex(j,k)
      t(i,j,k)=tex(j,k)
      uu(i,j,k)=tkex(j,k)
      vv(i,j,k)=tkex(j,k)
      ww(i,j,k)=tkex(j,k)
      uv(i,j,k)=0.
      vw(i,j,k)=0.
      uw(i,j,k)=0.
      HPTNS=SQRT ( (YC(I,NJ-1,K)-YC(I,2,K))**2 +
     1             (XC(I,NJ-1,K)-XC(I,2,K))**2 )
      COSTHT=(YC(I,NJ-1,K)-YC(I,2,K))/HPTNS
      SINTHT= (XC(I,NJ-1,K)-XC(I,2,K))/HPTNS
      u(i,j,k)=uex(j,k)*costht+vex(j,k)*sintht
      v(i,j,k)=vex(j,k)*costht-uex(j,k)*sintht
C
 20   CONTINUE
C
C      
      ENDIF
C
      DO 21 J=1,NJ
      DO 21 I=1,NI
C      TED(I,J,1)=0.
      TED(I,J,NK)=0.
      TKE(I,J,NK)=0.
C      TKE(I,J,1)=0.
   21 CONTINUE
C
C--- GEOMETRICAL INFORMATION RELATED TO WALL REFLECTION TERMS ---------------
C******** THIS PART  DEPENDS ON THE GEOMETRY OF WALL BOUNDARIES *********
C----------------------------------------------------------------------------
      DO 41 K=2,NKM1
      DO 41 J=2,NJM1
      DO 41 I=2,NIM1
C-----------------------------------------------------------------------------
C        MINIMUM WALL DISTANCE                                               *
C-----------------------------------------------------------------------------
C
      ZT=0.25*(Z(I,J,NKM1)+Z(I-1,J,NKM1)+Z(I,J-1,NKM1)+Z(I-1,J-1,NKM1))
      ZB=0.25*(Z(I,J,1)+Z(I-1,J,1)+Z(I,J-1,1)+Z(I-1,J-1,1))
      XLT=ZT-ZC(I,J,K)
      XLB=ZC(I,J,K)-ZB
      XE=SQRT(GREAT)
      YE=SQRT(GREAT)
      XW=SQRT(GREAT)
      YW=SQRT(GREAT)
      XLW=GREAT
      XLE=GREAT
C
      XS=0.25*(X(I,1,K)+X(I-1,1,K)+X(I-1,1,K-1)+X(I,1,K-1))
      XN=0.25*(X(I,NJM1,K)+X(I-1,NJM1,K)+X(I-1,NJM1,K-1)+X(I,NJM1,K-1))
      YS=0.25*(Y(I,1,K)+Y(I-1,1,K)+Y(I-1,1,K-1)+Y(I,1,K-1))
      YN=0.25*(Y(I,NJM1,K)+Y(I-1,NJM1,K)+Y(I-1,NJM1,K-1)+Y(I,NJM1,K-1))
      xs=xc(i,1,k)
      ys=yc(i,1,k)
      xn=xc(i,nj,k)
      yn=yc(i,nj,k)
      XLS=SQRT( (XC(I,J,K)-XS)**2 + (YC(I,J,K)-YS)**2 )
      XLN=SQRT( (XN-XC(I,J,K))**2 + (YN-YC(I,J,K))**2 )
      COSXN(I,J,K)=-(XN-XC(I,J,K))/XLN
      COSYN(I,J,K)= -(YN-YC(I,J,K))/XLN
      COSXS(I,J,K)= (XC(I,J,K)-XS)/XLS
      COSYS(I,J,K)= (YC(I,J,K)-YS)/XLS
C
      XLS=XLS-FLOAT(LWS-1)*GREAT
      XLN=XLN-FLOAT(LWN-1)*GREAT
      XLT=XLT-FLOAT(LWT-1)*GREAT
      XLB=XLB-FLOAT(LWB-1)*GREAT
C
      XL3(I,J,K)=AMIN1(XLS,XLN,XLT,XLB,XLE,XLW)
C 
C
C----------------------------------------------------------------------------
C         EAST WALL
C----------------------------------------------------------------------------
C        INVERSE OF MINIMUM WALL DISTANCE
C      ELE(I,J,K)=MAX((-COSXN(I,J,K)/XLN),(-COSXS(I,J,K)/XLS),0.)
      ELE(I,J,K)=0.
C        DIRECTION  COSINES
      COSXE(I,J,K)=0.
      COSYE(I,J,K)=0.
      COSZE(I,J,K)=0. 
C----------------------------------------------------------------------------
C         WEST WALL
C----------------------------------------------------------------------------
C        INVERSE OF MINIMUM WALL DISTANCE
      ELW(I,J,K)=0.
C      ELW(I,J,K)=MAX((COSXN(I,J,K)/XLN),(COSXS(I,J,K)/XLS),0.)
C        DIRECTION  COSINES
      COSXW(I,J,K)=0.
      COSYW(I,J,K)=0.
      COSZW(I,J,K)=0. 
C----------------------------------------------------------------------------
C         NORTH WALL
C----------------------------------------------------------------------------
C        INVERSE OF MINIMUM WALL DISTANCE
      ELN(I,J,K)=1./XLN
C        DIRECTION  COSINES
c      COSXN(I,J,K)=0.
c      COSYN(I,J,K)=1.
      COSZN(I,J,K)=0. 
C
C----------------------------------------------------------------------------
C         SOUTH WALL
C----------------------------------------------------------------------------
C        INVERSE OF MINIMUM WALL DISTANCE
      ELS(I,J,K)=1./XLS
c      COSXS(I,J,K)=0.
c      COSYS(I,J,K)=1.
      COSZS(I,J,K)=0. 
C
c      
C----------------------------------------------------------------------------
C         TOP WALL
C----------------------------------------------------------------------------
C        INVERSE OF MINIMUM WALL DISTANCE
      ELT(I,J,K)=1./XLT
C        DIRECTION  COSINES
      COSXT(I,J,K)=0.
      COSYT(I,J,K)=0.
      COSZT(I,J,K)=-1. 
C----------------------------------------------------------------------------
C         BOTTOM WALL
C----------------------------------------------------------------------------
C        INVERSE OF MINIMUM WALL DISTANCE
      ELB(I,J,K)=1./XLB
C        DIRECTION  COSINES
      COSXB(I,J,K)=0.
      COSYB(I,J,K)=0.
      COSZB(I,J,K)=0. 
C
C
   41 CONTINUE
C
      nkh=nk/2
      DO J=1,NJ
      DO I=1,NI
      ELNI(I,J)=ELN(I,J,2)
      ELSI(I,J)=ELS(I,J,2)   
      ELEI(I,J)=ELE(I,J,2)
      ELWI(I,J)=ELW(I,J,2)
      ELTI(I,J)=ELT(I,J,nkh)
      X2D(I,J)=XC(I,J,2)
      Y2D(I,J)=YC(I,J,2)
      ENDDO
      enddo
      WRITE(45,*)X2D,Y2D,ELNI,ELSI,ELEI,ELWI,ELTI
      RETURN
      END
C
      SUBROUTINE CALCU
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
      PR=1.
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*U(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFICIENTS  ----------------------------
C
      CALL ADVECT(1.,PR)
C
      CALL BQUICK(0,1.,U)
C
      DO 10 K=2,NKM1
      DO 10 J=2,NJM1
      DO 10 I=2,NIM1
C
C---------- CELL FACE DENSITIES ------------------------------------
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
C     NON-MUSCL APPROACH
C
      CE=DENE*CX(I  ,J,K)
      CW=DENW*CX(I-1,J,K)
      CN=DENN*CY(I,J  ,K)
      CS=DENS*CY(I,J-1,K)
      CT=DENT*CZ(I,J,K  )
      CB=DENB*CZ(I,J,K-1)
C
C      CEPLUS=0.5*(1.+SIGN(1.,CE))
C      CEMINS=0.5*(1.-SIGN(1.,CE))
C      CWPLUS=0.5*(1.+SIGN(1.,CW))
C      CWMINS=0.5*(1.-SIGN(1.,CW))
C      CNPLUS=0.5*(1.+SIGN(1.,CN))
C      CNMINS=0.5*(1.-SIGN(1.,CN))
C      CSPLUS=0.5*(1.+SIGN(1.,CS))
C      CSMINS=0.5*(1.-SIGN(1.,CS))
C      CTPLUS=0.5*(1.+SIGN(1.,CT))
C      CTMINS=0.5*(1.-SIGN(1.,CT))
C      CBPLUS=0.5*(1.+SIGN(1.,CB))
C      CBMINS=0.5*(1.-SIGN(1.,CB))
C
C      DIP3O2= U(  I+2,J,K  )-U(  I+1,J,K  )
C      DIP1O2= U(  I+1,J,K  )-U(  I  ,J,K  )
C      DIM1O2= U(  I  ,J,K  )-U(  I-1,J,K  )
C      DIM3O2= U(  I-1,J,K  )-U(  I-2,J,K  )
C
C      DJP3O2= U(  I,J+2,K  )-U(  I,J+1,K  )
C      DJP1O2= U(  I,J+1,K  )-U(  I,J  ,K  )
C      DJM1O2= U(  I,J  ,K  )-U(  I,J-1,K  )
C      DJM3O2= U(  I,J-1,K  )-U(  I,J-2,K  )
C
C      DKP3O2= U(  I,  J,K+2)-U(  I,  J,K+1)
C      DKP1O2= U(  I,  J,K+1)-U(  I,  J,K  )
C      DKM1O2= U(  I,  J,K  )-U(  I,  J,K-1)
C      DKM3O2= U(  I,  J,K-1)-U(  I,  J,K-2)
C
C     B.C.
C
C      IF(I.EQ.2) THEN
C      DIM1O2=2.*DIM1O2
C      DIM3O2=0.
C      END IF
C
C      IF(J.EQ.2) THEN
C      DJM1O2=2.*DJM1O2
C      DJM3O2=0.
C      END IF
C
C      IF(K.EQ.2) THEN
C      DKM1O2=2.*DKM1O2
C      DKM3O2=0.
C      END IF
C
C      IF(I.EQ.NIM1) THEN
C      DIP1O2=2.*DIP1O2
C      DIP3O2=0.
C      END IF
C
C      IF(J.EQ.NJM1) THEN
C      DJP1O2=2.*DJP1O2
C      DJP3O2=0.
C      END IF
C
C      IF(K.EQ.NKM1) THEN
C      DKP1O2=2.*DKP1O2
C      DKP3O2=0.
C      END IF
C
C       SUE=
C     1 CE*CEPLUS*
C     1 ( (1.-ETA)*DIM1O2+(1.+ETA)*DIP1O2 )
C     1-CE*CEMINS*
C     1 ( (1.-ETA)*DIP3O2+(1.+ETA)*DIP1O2 )
C
C       SUW=
C     1 CW*CWPLUS*
C     1 ( (1.-ETA)*DIM3O2+(1.+ETA)*DIM1O2 )
C     1-CW*CWMINS*
C     1 ( (1.-ETA)*DIP1O2+(1.+ETA)*DIM1O2 )
C
C       SUN=
C     1 CN*CNPLUS*
C     1 ( (1.-ETA)*DJM1O2+(1.+ETA)*DJP1O2 )
C     1-CN*CNMINS*
C     1 ( (1.-ETA)*DJP3O2+(1.+ETA)*DJP1O2 )
C
C       SUS=
C     1 CS*CSPLUS*
C     1 ( (1.-ETA)*DJM3O2+(1.+ETA)*DJM1O2 )
C     1-CS*CSMINS*
C     1 ( (1.-ETA)*DJP1O2+(1.+ETA)*DJM1O2 )
C
C       SUT=
C     1 CT*CTPLUS*
C     1 ( (1.-ETA)*DKM1O2+(1.+ETA)*DKP1O2 )
C     1-CT*CTMINS*
C     1 ( (1.-ETA)*DKP3O2+(1.+ETA)*DKP1O2 )
C
C       SUB=
C     1 CB*CBPLUS*
C     1 ( (1.-ETA)*DKM3O2+(1.+ETA)*DKM1O2 )
C     1-CB*CBMINS*
C     1 ( (1.-ETA)*DKP1O2+(1.+ETA)*DKM1O2 )
C
C      IF(I.EQ.NIM1) SUE=0.
C      IF(J.EQ.NJM1) SUN=0.
C      IF(K.EQ.NKM1) SUT=0.
C      IF(I.EQ.2   ) SUW=0.
C      IF(J.EQ.2   ) SUS=0.
C      IF(K.EQ.2   ) SUB=0.
C
C      SU(I,J,K)=SU(I,J,K)-0.25*( SUE-SUW+SUN-SUS+SUT-SUB )*RATIO
C
C------ REYNOLDS STRESS GRADIENT TERMS ------------------------------
C
C      EAST-WEST FACE AREAS
C
      AREXE=(1.-FX(I,J,K))*  XX(I,J,K)+  FX(I,J,K)*  XX(I+1,J,K)
      AREXW=(1.-FX(I-1,J,K))*XX(I-1,J,K)+FX(I-1,J,K)*XX(I,J,K)
      AREYE=(1.-FX(I,J,K))*  XY(I,J,K)+  FX(I,J,K)*  XY(I+1,J,K)
      AREYW=(1.-FX(I-1,J,K))*XY(I-1,J,K)+FX(I-1,J,K)*XY(I,J,K)
      AREZE=(1.-FX(I,J,K))*  XZ(I,J,K)+  FX(I,J,K)*  XZ(I+1,J,K)
      AREZW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)+FX(I-1,J,K)*XZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      UEST=(1-FX(I,J,K))*  U(I,J,K)+  FX(I,J,K)*  U(I+1,J,K)
      UWST=(1-FX(I-1,J,K))*U(I-1,J,K)+FX(I-1,J,K)*U(I,J,K)
      UEEST=UEST+(U(I+1,J,K)-U(I,J,K))
      IF(I.NE.NIM1)
     1 UEEST=(1-FX(I+1,J,K))*U(I+1,J,K)+FX(I+1,J,K)*U(I+2,J,K)
      UWWST=UWST-(U(I,J,K)-U(I-1,J,K))
      IF(I.EQ.NIM1)
     1 UEEST=(1-FX(4,J,K))*U(4,J,K)+FX(4,J,K)*U(5,J,K)
      UWWST=UWST-(U(I,J,K)-U(I-1,J,K))
      IF(I.NE.2)
     1 UWWST=(1-FX(I-2,J,K))*U(I-2,J,K)+FX(I-2,J,K)*U(I-1,J,K)
      IF(I.EQ.2)
     1 UWWST=(1-FX(NI-4,J,K))*U(NI-4,J,K)+FX(NI-4,J,K)*U(NI-3,J,K)
C
C      EAST-WEST FACE STRESSES
C
C      FX(1,J,K)=0.5
C      FX(NIM1,J,K)=0.5
C
      UUE=(1.-FX(I,J,K))*  UU(I,J,K)+  FX(I,J,K)*  UU(I+1,J,K)
      UUW=(1.-FX(I-1,J,K))*UU(I-1,J,K)+FX(I-1,J,K)*UU(I,J,K)
      UVE=(1.-FX(I,J,K))*  UV(I,J,K)+  FX(I,J,K)*  UV(I+1,J,K)
      UVW=(1.-FX(I-1,J,K))*UV(I-1,J,K)+FX(I-1,J,K)*UV(I,J,K)
      UWE=(1.-FX(I,J,K))*  UW(I,J,K)+  FX(I,J,K)*  UW(I+1,J,K)
      UWW=(1.-FX(I-1,J,K))*UW(I-1,J,K)+FX(I-1,J,K)*UW(I,J,K)
C
C        EAST-WEST CONTRIBUTION
C
      SUSTW=DENW*(AREXW*UUW+AREYW*UVW+AREZW*UWW)
      SUSTE=DENE*(AREXE*UUE+AREYE*UVE+AREZE*UWE)
C --- RHIE AND CHOW ------------------
C      if(j.eq.2.or.j.eq.njm1.or.k.eq.nkm1)go to 128
C      GO TO 128
      SUSTW=SUSTW+
     1            (1.-FX(I-1,J,K))*(UWST-UWWST)*(
     1     ((XX(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS)+
     1      (XY(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS)+
     1      (XZ(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS))/VOL(I-1,J,K))
     1          + FX(I-1,J,K)*(UEST-UWST)*(
     1     ((XX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTE=SUSTE+
     1            (1.-FX(I,J,K))*(UEST-UWST)*(
     1     ((XX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FX(I,J,K)*(UEEST-UEST)*(
     1     ((XX(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS)+
     1      (XY(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS)+
     1      (XZ(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS))/VOL(I+1,J,K))
C
  128 continue
C      FX(1,J,K)=0.
C      FX(NIM1,J,K)=1.
C
C      NORTH-SOUTH FACE AREAS
C
      AREXN=(1.-FY(I,J,K))*  YX(I,J,K)+  FY(I,J,K)*  YX(I,J+1,K)
      AREXS=(1.-FY(I,J-1,K))*YX(I,J-1,K)+FX(I,J-1,K)*YX(I,J,K)
      AREYN=(1.-FY(I,J,K))*  YY(I,J,K)+  FY(I,J,K)*  YY(I,J+1,K)
      AREYS=(1.-FY(I,J-1,K))*YY(I,J-1,K)+FY(I,J-1,K)*YY(I,J,K)
      AREZN=(1.-FY(I,J,K))*  YZ(I,J,K)+  FY(I,J,K)*  YZ(I,J+1,K)
      AREZS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)+FY(I,J-1,K)*YZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      UNRT=(1-FY(I,J,K))*  U(I,J,K)+  FY(I,J,K)*  U(I,J+1,K)
      USTH=(1-FY(I,J-1,K))*U(I,J-1,K)+FY(I,J-1,K)*U(I,J,K)
      UNNRT=UNRT+(U(I,J+1,K)-U(I,J,K))
      IF(J.NE.NJM1)
     1 UNNRT=(1-FY(I,J+1,K))*U(I,J+1,K)+FY(I,J+1,K)*U(I,J+2,K)
      USSTH=USTH-(U(I,J,K)-U(I,J-1,K))
      IF(J.NE.2)
     1 USSTH=(1-FY(I,J-2,K))*U(I,J-2,K)+FY(I,J-2,K)*U(I,J-1,K)
C
C      NORTH-SOUTH FACE STRESSES
C
C      FY(I,1,K)=0.5
C     FY(I,NJM1,K)=.5
C
      UUN=(1.-FY(I,J,K))*  UU(I,J,K)+  FY(I,J,K)*  UU(I,J+1,K)
      UUS=(1.-FY(I,J-1,K))*UU(I,J-1,K)+FY(I,J-1,K)*UU(I,J,K)
      UVN=(1.-FY(I,J,K))*  UV(I,J,K)+  FY(I,J,K)*  UV(I,J+1,K)
      UVS=(1.-FY(I,J-1,K))*UV(I,J-1,K)+FY(I,J-1,K)*UV(I,J,K)
      UWN=(1.-FY(I,J,K))*  UW(I,J,K)+  FY(I,J,K)*  UW(I,J+1,K)
      UWS=(1.-FY(I,J-1,K))*UW(I,J-1,K)+FY(I,J-1,K)*UW(I,J,K)
C
C        NORTH-SOUTH CONTRIBUTION
C
      SUSTS=DENS*(AREXS*UUS+AREYS*UVS+AREZS*UWS)
      SUSTN=DENN*(AREXN*UUN+AREYN*UVN+AREZN*UWN)
C --- RHIE AND CHOW ------------------
C      if(j.eq.2.or.j.eq.njm1.or.k.eq.nkm1)go to 129
C      GO TO 129
      SUSTS=SUSTS+
     1            (1.-FY(I,J-1,K))*(USTH-USSTH)*(
     1     ((YX(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS)+
     1      (YY(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS)+
     1      (YZ(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS))/VOL(I,J-1,K))
     1          + FY(I,J-1,K)*(UNRT-USTH)*(
     1     ((YX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTN=SUSTN+
     1            (1.-FY(I,J,K))*(UNRT-USTH)*(
     1     ((YX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FY(I,J,K)*(UNNRT-UNRT)*(
     1     ((YX(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS)+
     1      (YY(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS)+
     1      (YZ(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS))/VOL(I,J+1,K))
C
c      SUSTS=DENS*(AREXS*UUS+AREYS*UVS+AREZS*UWS)
c      SUSTN=DENN*(AREXN*UUN+AREYN*UVN+AREZN*UWN)
  129 continue
C      FY(I,1,K)=0.
C      FY(I,NJM1,K)=1.
C
C      TOP-BOTTOM FACE AREAS
C
      AREXT=(1.-FZ(I,J,K))*  ZX(I,J,K)+  FZ(I,J,K)*  ZX(I,J,K+1)
      AREXB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)+FZ(I,J,K-1)*ZX(I,J,K)
      AREYT=(1.-FZ(I,J,K))*  ZY(I,J,K)+  FZ(I,J,K)*  ZY(I,J,K+1)
      AREYB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)+FZ(I,J,K-1)*ZY(I,J,K)
      AREZT=(1.-FZ(I,J,K))*  ZZ(I,J,K)+  FZ(I,J,K)*  ZZ(I,J,K+1)
      AREZB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)+FZ(I,J,K-1)*ZZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      UTOP=(1-FZ(I,J,K))*  U(I,J,K)+  FZ(I,J,K)*  U(I,J,K+1)
      UBOT=(1-FZ(I,J,K-1))*U(I,J,K-1)+FZ(I,J,K-1)*U(I,J,K)
      UTTOP=UTOP+(U(I,J,K+1)-U(I,J,K))
      IF(K.NE.NKM1)
     1 UTTOP=(1-FZ(I,J,K+1))*U(I,J,K+1)+FZ(I,J,K+1)*U(I,J,K+2)
      UBBOT=UBOT-(U(I,J,K)-U(I,J,K-1))
      IF(K.NE.2)
     1 UBBOT=(1-FZ(I,J,K-2))*U(I,J,K-2)+FZ(I,J,K-2)*U(I,J,K-1)
C
C      TOP-BOTTOM FACE STRESSES
C
C      FZ(I,J,1)=0.5
c      FZ(I,J,NKM1)=.5
C
      UUT=(1.-FZ(I,J,K))*  UU(I,J,K)+  FZ(I,J,K)*  UU(I,J,K+1)
      UUB=(1.-FZ(I,J,K-1))*UU(I,J,K-1)+FZ(I,J,K-1)*UU(I,J,K)
      UVT=(1.-FZ(I,J,K))  *UV(I,J,K)+  FZ(I,J,K)*  UV(I,J,K+1)
      UVB=(1.-FZ(I,J,K-1))*UV(I,J,K-1)+FZ(I,J,K-1)*UV(I,J,K)
      UWT=(1.-FZ(I,J,K))*  UW(I,J,K)+  FZ(I,J,K)*  UW(I,J,K+1)
      UWB=(1.-FZ(I,J,K-1))*UW(I,J,K-1)+FZ(I,J,K-1)*UW(I,J,K)
C
C        NORTH-SOUTH CONTRIBUTION
C
      SUSTB=DENB*(AREXB*UUB+AREYB*UVB+AREZB*UWB)
      SUSTT=DENT*(AREXT*UUT+AREYT*UVT+AREZT*UWT)
C      if(k.eq.2.or.k.eq.nkm1.or.j.eq.2.or.j.eq.njm1)go to 130
C      GO TO 130
C --- RHIE AND CHOW ------------------
      SUSTB=SUSTB+
     1            (1.-FZ(I,J,K-1))*(UBOT-UBBOT)*(
     1     ((ZX(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS)+
     1      (ZY(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS)+
     1      (ZZ(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS))/VOL(I,J,K-1))
     1          + FZ(I,J,K-1)*(UTOP-UBOT)*(
     1     ((ZX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTT=SUSTT+
     1            (1.-FZ(I,J,K))*(UTOP-UBOT)*(
     1     ((ZX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FZ(I,J,K)*(UTTOP-UTOP)*(
     1     ((ZX(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS)+
     1      (ZY(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS)+
     1      (ZZ(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS))/VOL(I,J,K+1))
C   
  130 continue
C      FZ(I,J,1)=0.
C      FZ(I,J,NKM1)=1.
C
C   OVERALL STRESS GRADIENT TERM
C
      SU(I,J,K)=SU(I,J,K)+FLOAT(1-IEVM)*
     1          (SUSTW-SUSTE+SUSTS-SUSTN+SUSTB-SUSTT)
   10 CONTINUE
C
C -------- PRESSURE GRADIENT SOURCE TERM  ----------------------------
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C      CELL FACE PRESSURES
C
      PRESE=P(I,J,K)*(1.-FX(I,J,K))+P(I+1,J,K)*FX(I,J,K)
      PRESW=P(I-1,J,K)*(1.-FX(I-1,J,K))+P(I,J,K)*FX(I-1,J,K)
      PRESN=P(I,J,K)*(1.-FY(I,J,K))+P(I,J+1,K)*FY(I,J,K)
      PRESS=P(I,J-1,K)*(1.-FY(I,J-1,K))+P(I,J,K)*FY(I,J-1,K)
      PREST=P(I,J,K)*(1.-FZ(I,J,K))+P(I,J,K+1)*FZ(I,J,K)
      PRESB=P(I,J,K-1)*(1.-FZ(I,J,K-1))+P(I,J,K)*FZ(I,J,K-1)
C
C       AREA VECTOR AT CELL CENTER 
C
      DU1(I,J,K)=XX(I,J,K)
      DU2(I,J,K)=YX(I,J,K)
      DU3(I,J,K)=ZX(I,J,K)
C
      SU(I,J,K)=SU(I,J,K)+DU1(I,J,K)*(PRESW-PRESE)
     1 +DU2(I,J,K)*(PRESS-PRESN)+DU3(I,J,K)*(PRESB-PREST)
   50 CONTINUE
C
C------ DIFFUSION  SOURCE TERMS ------------------------------------
C
      DO 60 K=2,NKM1
      DO 60 J=2,NJM1
      DO 60 I=2,NIM1
C
C ------ CONTROL VOLUMES FOR DUMMY (STAGGERED) VELOCITY NODES ----
C
       VOLE=VOL(I,J,K)*(1.-FX(I,J,K))+VOL(I+1,J,K)*FX(I,J,K)
       VOLW=VOL(I-1,J,K)*(1.-FX(I-1,J,K))+VOL(I,J,K)*FX(I-1,J,K)
       VOLN=VOL(I,J,K)*(1.-FY(I,J,K))+VOL(I,J+1,K)*FY(I,J,K)
       VOLS=VOL(I,J-1,K)*(1.-FY(I,J-1,K))+VOL(I,J,K)*FY(I,J-1,K)
       VOLT=VOL(I,J,K)*(1.-FZ(I,J,K))+VOL(I,J,K+1)*FZ(I,J,K)
       VOLB=VOL(I,J,K-1)*(1.-FZ(I,J,K-1))+VOL(I,J,K)*FZ(I,J,K-1)
C
C ------ CELL FACE VISCOCITIES ------------------------------------
C
       VISE=VIS(I,J,K)*(1.-FX(I,J,K))+VIS(I+1,J,K)*FX(I,J,K)
       VISW=VIS(I-1,J,K)*(1.-FX(I-1,J,K))+VIS(I,J,K)*FX(I-1,J,K)
       VISN=VIS(I,J,K)*(1.-FY(I,J,K))+VIS(I,J+1,K)*FY(I,J,K)
       VISS=VIS(I,J-1,K)*(1.-FY(I,J-1,K))+VIS(I,J,K)*FY(I,J-1,K)
       VIST=VIS(I,J,K)*(1.-FZ(I,J,K))+VIS(I,J,K+1)*FZ(I,J,K)
       VISB=VIS(I,J,K-1)*(1.-FZ(I,J,K-1))+VIS(I,J,K)*FZ(I,J,K-1)
C
       VISE=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISE
       VISW=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISW
       VISN=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISN
       VISS=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISS
       VIST=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VIST
       VISB=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISB
C
C     CALCULATE THE CORNER VALUES FOR U,V,W
C
      UTNE=
     1 U(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+U(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      UTNW=
     1 U(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+U(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      UBNE=
     1 U(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+U(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      UBNW=
     1 U(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+U(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+U(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      UTSE=
     1 U(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+U(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      UTSW=
     1 U(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+U(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      UBSE=
     1 U(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+U(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+U(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      UBSW=
     1 U(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+U(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      VTNE=
     1 V(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+V(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      VTNW=
     1 V(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+V(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      VBNE=
     1 V(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+V(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      VBNW=
     1 V(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+V(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+V(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      VTSE=
     1 V(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+V(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      VTSW=
     1 V(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+V(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      VBSE=
     1 V(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+V(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+V(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      VBSW=
     1 V(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+V(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      WTNE=
     1 W(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+W(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      WTNW=
     1 W(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+W(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      WBNE=
     1 W(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+W(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      WBNW=
     1 W(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+W(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+W(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      WTSE=
     1 W(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+W(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      WTSW=
     1 W(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+W(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      WBSE=
     1 W(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+W(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+W(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      WBSW=
     1 W(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+W(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
C ---- CROSS-DIFFUSION SOURCE TERMS  ------------------------------
C
C         EAST-WEST FACE CROSS-DIFFUSION
C
      Q12E=(XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*YX(I+1,J,K)+XY(I+1,J,K)*YY(I+1,J,K)+
     1 XZ(I+1,J,K)*YZ(I+1,J,K))*FX(I,J,K)
C
      Q12W=(XX(I-1,J,K)*YX(I-1,J,K)+XY(I-1,J,K)*YY(I-1,J,K)+
     1 XZ(I-1,J,K)*YZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*FX(I-1,J,K)
C
      Q13E=(XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*ZX(I+1,J,K)+XY(I+1,J,K)*ZY(I+1,J,K)+
     1 XZ(I+1,J,K)*ZZ(I+1,J,K))*FX(I,J,K)
C
      Q13W=(XX(I-1,J,K)*ZX(I-1,J,K)+XY(I-1,J,K)*ZY(I-1,J,K)+
     1 XZ(I-1,J,K)*ZZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*FX(I-1,J,K)
C        
      SCX=0.5*VISE*Q12E/VOLE*(UTNE+UBNE-UTSE-UBSE)
     1   +0.5*VISE*Q13E/VOLE*(UTNE+UTSE-UBNE-UBSE)
     1   -0.5*VISW*Q12W/VOLW*(UTNW+UBNW-UTSW-UBSW)
     1   -0.5*VISW*Q13W/VOLW*(UTNW+UTSW-UBNW-UBSW)
C
C        NORTH-SOUTH FACE CROSS-DIFFUSION
C
      Q21N=(YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*XX(I,J+1,K)+YY(I,J+1,K)*XY(I,J+1,K)+
     1 YZ(I,J+1,K)*XZ(I,J+1,K))*FY(I,J,K)
C
      Q21S=(YX(I,J-1,K)*XX(I,J-1,K)+YY(I,J-1,K)*XY(I,J-1,K)+
     1 YZ(I,J-1,K)*XZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*FY(I,J-1,K)
C
      Q23N=(YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*ZX(I,J+1,K)+YY(I,J+1,K)*ZY(I,J+1,K)+
     1 YZ(I,J+1,K)*ZZ(I,J+1,K))*FY(I,J,K)
C
      Q23S=(YX(I,J-1,K)*ZX(I,J-1,K)+YY(I,J-1,K)*ZY(I,J-1,K)+
     1 YZ(I,J-1,K)*ZZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*FY(I,J-1,K)
C
      SCY=0.5*VISN*Q21N/VOLN*(UTNE+UBNE-UTNW-UBNW)
     1   +0.5*VISN*Q23N/VOLN*(UTNE+UTNW-UBNE-UBNW)
     1   -0.5*VISS*Q21S/VOLS*(UTSE+UBSE-UTSW-UBSW)
     1   -0.5*VISS*Q23S/VOLS*(UTSE+UTSW-UBSE-UBSW)
C
C     TOP-BOTTOM FACE CROSS-DIFUSION
C
      Q31T=(ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*XX(I,J,K+1)+ZY(I,J,K+1)*XY(I,J,K+1)+
     1 ZZ(I,J,K+1)*XZ(I,J,K+1))*FZ(I,J,K)
C
      Q31B=(ZX(I,J,K-1)*XX(I,J,K-1)+ZY(I,J,K-1)*XY(I,J,K-1)+
     1 ZZ(I,J,K-1)*XZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*FZ(I,J,K-1)
C
      Q32T=(ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*YX(I,J,K+1)+ZY(I,J,K+1)*YY(I,J,K+1)+
     1 ZZ(I,J,K+1)*YZ(I,J,K+1))*FZ(I,J,K)
C
      Q32B=(ZX(I,J,K-1)*YX(I,J,K-1)+ZY(I,J,K-1)*YY(I,J,K-1)+
     1 ZZ(I,J,K-1)*YZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*FZ(I,J,K-1)
C
      SCZ=0.5*VIST*Q31T/VOLT*(UTNE+UTSE-UTNW-UTSW)
     1   +0.5*VIST*Q32T/VOLT*(UTNE+UTNW-UTSE-UTSW)
     1   -0.5*VISB*Q31B/VOLB*(UBNE+UBSE-UBNW-UBSW)
     1   -0.5*VISB*Q32B/VOLB*(UBNE+UBNW-UBSE-UBSW)
C
      SU(I,J,K)=SU(I,J,K)+SCX+SCY+SCZ
C
C------ VARIABLE VISCOSITY SOURCE TERMS  --------------------------
C
C        EAST WEST U,V,W GRADIENT TERMS
C
      XXXXE=(1.-FX(I,J,K))*XX(I,J,K)*XX(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*XX(I+1,J,K)
      XYXXE=(1.-FX(I,J,K))*XY(I,J,K)*XX(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*XX(I+1,J,K)
      XZXXE=(1.-FX(I,J,K))*XZ(I,J,K)*XX(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*XX(I+1,J,K)
      XXYXE=(1.-FX(I,J,K))*XX(I,J,K)*YX(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*YX(I+1,J,K)
      XYYXE=(1.-FX(I,J,K))*XY(I,J,K)*YX(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*YX(I+1,J,K)
      XZYXE=(1.-FX(I,J,K))*XZ(I,J,K)*YX(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*YX(I+1,J,K)
      XXZXE=(1.-FX(I,J,K))*XX(I,J,K)*ZX(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*ZX(I+1,J,K)
      XYZXE=(1.-FX(I,J,K))*XY(I,J,K)*ZX(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*ZX(I+1,J,K)
      XZZXE=(1.-FX(I,J,K))*XZ(I,J,K)*ZX(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*ZX(I+1,J,K)
C
      XXXXW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*XX(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*XX(I,J,K)
      XYXXW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*XX(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*XX(I,J,K)
      XZXXW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*XX(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*XX(I,J,K)
      XXYXW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*YX(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*YX(I,J,K)
      XYYXW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*YX(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*YX(I,J,K)
      XZYXW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*YX(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*YX(I,J,K)
      XXZXW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*ZX(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*ZX(I,J,K)
      XYZXW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*ZX(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*ZX(I,J,K)
      XZZXW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*ZX(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*ZX(I,J,K)
C
      SC1=VISE/VOLE*( (U(I+1,J,K)-U(I,J,K))*XXXXE+(V(I+1,J,K)-V(I,J,K))
     1 *XYXXE+(W(I+1,J,K)-W(I,J,K))*XZXXE )
     1   -VISW/VOLW*( (U(I,J,K)-U(I-1,J,K))*XXXXW+(V(I,J,K)-V(I-1,J,K))
     1 *XYXXW+(W(I,J,K)-W(I-1,J,K))*XZXXW )
C
      SC2=VISE/VOLE*( (UTNE+UBNE-UTSE-UBSE)*XXYXE+(VTNE+VBNE-VTSE-VBSE)
     1 *XYYXE+(WTNE+WBNE-WTSE-WBSE)*XZYXE )
     1   -VISW/VOLW*( (UTNW+UBNW-UTSW-UBSW)*XXYXW+(VTNW+VBNW-VTSW-VBSW)
     1 *XYYXW+(WTNW+WBNW-WTSW-WBSW)*XZYXW )
      SC2=0.5*SC2
C
      SC3=VISE/VOLE*( (UTNE+UTSE-UBNE-UBSE)*XXZXE+(VTNE+VTSE-VBNE-VBSE)
     1 *XYZXE+(WTNE+WTSE-WBNE-WBSE)*XZZXE )
     1   -VISW/VOLW*( (UTNW+UTSW-UBNW-UBSW)*XXZXW+(VTNW+VTSW-VBNW-VBSW)
     1 *XYZXW+(WTNW+WTSW-WBNW-WBSW)*XZZXW )
      SC3=0.5*SC3
C
C      NORTH-SOUTH U,V,W GRADIENT TERMS
C
      YXXXN=(1.-FY(I,J,K))*YX(I,J,K)*XX(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*XX(I,J+1,K)
      YYXXN=(1.-FY(I,J,K))*YY(I,J,K)*XX(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*XX(I,J+1,K)
      YZXXN=(1.-FY(I,J,K))*YZ(I,J,K)*XX(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*XX(I,J+1,K)
      YXYXN=(1.-FY(I,J,K))*YX(I,J,K)*YX(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*YX(I,J+1,K)
      YYYXN=(1.-FY(I,J,K))*YY(I,J,K)*YX(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*YX(I,J+1,K)
      YZYXN=(1.-FY(I,J,K))*YZ(I,J,K)*YX(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*YX(I,J+1,K)
      YXZXN=(1.-FY(I,J,K))*YX(I,J,K)*ZX(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*ZX(I,J+1,K)
      YYZXN=(1.-FY(I,J,K))*YY(I,J,K)*ZX(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*ZX(I,J+1,K)
      YZZXN=(1.-FY(I,J,K))*YZ(I,J,K)*ZX(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*ZX(I,J+1,K)
C
      YXXXS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*XX(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*XX(I,J,K)
      YYXXS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*XX(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*XX(I,J,K)
      YZXXS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*XX(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*XX(I,J,K)
      YXYXS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*YX(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*YX(I,J,K)
      YYYXS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*YX(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*YX(I,J,K)
      YZYXS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*YX(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*YX(I,J,K)
      YXZXS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*ZX(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*ZX(I,J,K)
      YYZXS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*ZX(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*ZX(I,J,K)
      YZZXS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*ZX(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*ZX(I,J,K)
C
      SC4=VISN/VOLN*( (UTNE+UBNE-UTNW-UBNW)*YXXXN+(VTNE+VBNE-VTNW-VBNW)
     1 *YYXXN+(WTNE+WBNE-WTNW-WBNW)*YZXXN )
     1   -VISS/VOLS*( (UTSE+UBSE-UTSW-UBSW)*YXXXS+(VTSE+VBSE-VTSW-VBSW)
     1 *YYXXS+(WTSE+WBSE-WTSW-WBSW)*YZXXS )
      SC4=0.5*SC4
C
      SC5=VISN/VOLN*( (U(I,J+1,K)-U(I,J,K))*YXYXN+(V(I,J+1,K)-V(I,J,K))
     1 *YYYXN+(W(I,J+1,K)-W(I,J,K))*YZYXN )
     1   -VISS/VOLS*( (U(I,J,K)-U(I,J-1,K))*YXYXS+(V(I,J,K)-V(I,J-1,K))
     1 *YYYXS+(W(I,J,K)-W(I,J-1,K))*YZYXS )
C
      SC6=VISN/VOLN*( (UTNE+UTNW-UBNE-UBNW)*YXZXN+(VTNE+VTNW-VBNE-VBNW)
     1 *YYZXN+(WTNE+WTNW-WBNE-WBNW)*YZZXN )
     1   -VISS/VOLS*( (UTSE+UTSW-UBSE-UBSW)*YXZXS+(VTSE+VTSW-VBSE-VBSW)
     1 *YYZXS+(WTSE+WTSW-WBSE-WBSW)*YZZXS )
      SC6=0.5*SC6
C
C    TOP-BOTTOM U,V,W GRADIENT TERMS
C
      ZXXXT=(1.-FZ(I,J,K))*ZX(I,J,K)*XX(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*XX(I,J,K+1)
      ZYXXT=(1.-FZ(I,J,K))*ZY(I,J,K)*XX(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*XX(I,J,K+1)
      ZZXXT=(1.-FZ(I,J,K))*ZZ(I,J,K)*XX(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*XX(I,J,K+1)
      ZXYXT=(1.-FZ(I,J,K))*ZX(I,J,K)*YX(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*YX(I,J,K+1)
      ZYYXT=(1.-FZ(I,J,K))*ZY(I,J,K)*YX(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*YX(I,J,K+1)
      ZZYXT=(1.-FZ(I,J,K))*ZZ(I,J,K)*YX(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*YX(I,J,K+1)
      ZXZXT=(1.-FZ(I,J,K))*ZX(I,J,K)*ZX(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*ZX(I,J,K+1)
      ZYZXT=(1.-FZ(I,J,K))*ZY(I,J,K)*ZX(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*ZX(I,J,K+1)
      ZZZXT=(1.-FZ(I,J,K))*ZZ(I,J,K)*ZX(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*ZX(I,J,K+1)
C
      ZXXXB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*XX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*XX(I,J,K)
      ZYXXB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*XX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*XX(I,J,K)
      ZZXXB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*XX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*XX(I,J,K)
      ZXYXB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*YX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*YX(I,J,K)
      ZYYXB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*YX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*YX(I,J,K)
      ZZYXB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*YX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*YX(I,J,K)
      ZXZXB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*ZX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*ZX(I,J,K)
      ZYZXB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*ZX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*ZX(I,J,K)
      ZZZXB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*ZX(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*ZX(I,J,K)
C
      SC7=VIST/VOLT*( (UTNE+UTSE-UTNW-UTSW)*ZXXXT+(VTNE+VTSE-VTNW-VTSW)
     1 *ZYXXT+(WTNE+WTSE-WTNW-WTSW)*ZZXXT )
     1   -VISB/VOLB*( (UBNE+UBSE-UBNW-UBSW)*ZXXXB+(VBNE+VBSE-VBNW-VBSW)
     1 *ZYXXB+(WBNE+WBSE-WBNW-WBSW)*ZZXXB )
      SC7=0.5*SC7
C
      SC8=VIST/VOLT*( (UTNE+UTNW-UTSE-UTSW)*ZXYXT+(VTNE+VTNW-VTSE-VTSW)
     1 *ZYYXT+(WTNE+WTNW-WTSE-WTSW)*ZZYXT )
     1   -VISB/VOLB*( (UBNE+UBNW-UBSE-UBSW)*ZXYXB+(VBNE+VBNW-VBSE-VBSW)
     1 *ZYYXB+(WBNE+WBNW-WBSE-WBSW)*ZZYXB )
      SC8=0.5*SC8
C
      SC9=VIST/VOLT*( (U(I,J,K+1)-U(I,J,K))*ZXZXT+(V(I,J,K+1)-V(I,J,K))
     1 *ZYZXT+(W(I,J,K+1)-W(I,J,K))*ZZZXT )
     1   -VISB/VOLB*( (U(I,J,K)-U(I,J,K-1))*ZXZXB+(V(I,J,K)-V(I,J,K-1))
     1 *ZYZXB+(W(I,J,K)-W(I,J,K-1))*ZZZXB )
C
      SU(I,J,K)=SU(I,J,K)+SC1+SC2+SC3+SC4+SC5+SC6+SC7+SC8+SC9
C
 60   CONTINUE
C---------------------------------------------------------------------
C          ROTATION TERMS                                            *
C---------------------------------------------------------------------
      DO 70 K=2,NKM1
      DO 70 J=2,NJM1
      DO 70 I=2,NIM1
C----  CORIOLIS TERMS -----------------------------------------------
      SU(I,J,K)=SU(I,J,K)+2.*DEN(I,J,K)*VOL(I,J,K)*
     +                    (OMZ*V(I,J,K)-OMY*W(I,J,K))
C----  CENTRIFUGAL TERMS --------------------------------------------
      HXX=HX+XC(I,J,K)
      HYY=HY+YC(I,J,K)
      HZZ=HZ+ZC(I,J,K)
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*
     +                    (HXX*(OMY*OMY+OMZ*OMZ)-OMX*(HYY*OMY+HZZ*OMZ))
   70 CONTINUE
C------ BOUNDARY CONDITIONS   ---------------------------------------
      CALL MODU
C
      RESORU=0.
      DO 30 K=2,NKM1
      DO 30 J=2,NJM1
      DO 30 I=2,NIM1
      AN(I,J,K)=ABS(AN(I,J,K))
      AS(I,J,K)=ABS(AS(I,J,K))
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*U(I+1,J,K)+AW(I,J,K)*U(I-1,J,K)
     1     +AN(I,J,K)*U(I,J+1,K)+AS(I,J,K)*U(I,J-1,K)
     1     +AT(I,J,K)*U(I,J,K+1)+AB(I,J,K)*U(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*U(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
      RESORU=RESORU+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFU
      SU(I,J,K)=SU(I,J,K)+(1.-URFU)*AP(I,J,K)*U(I,J,K)
      DU1(I,J,K)=DU1(I,J,K)/AP(I,J,K)
      DU2(I,J,K)=DU2(I,J,K)/AP(I,J,K)
      DU3(I,J,K)=DU3(I,J,K)/AP(I,J,K)
 30   CONTINUE
   
   36 CONTINUE
C------ EXIT CONDITIONS ----------------------------------------------
      CALL MODXU
C
      DO 40 N=1,NSWPU
      CALL LISOLV(4,2,2,NI,NJ,NK,U)
C      CALL LISOLV(4,NJR+1,2,NI,NJRN+1,NK,U)
C      CALL LISOLV(NIR+1,2,2,NI-2,NJ,NK,U)
 40   CONTINUE
C      DO 77 J=1,NJ
C      DO 77 K=1,NK
C      DO 78 I=1,NIR
C      U(I,J,K)=MAX(U(I,J,K),0.)
C   78 CONTINUE
C      DO 79 I=NI-2,NI
C      U(I,J,K)=MAX(U(I,J,K),0.)
C   79 CONTINUE
C   77 CONTINUE
      CALL MODU
      CALL MODXU
      RETURN
      END
C
      SUBROUTINE CALCV
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
      PR=1
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*V(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFICIENTS  ----------------------------
C
      CALL ADVECT(1.,PR)
      CALL BQUICK(0,1.,V)
C
      DO 10 K=2,NKM1
      DO 10 J=2,NJM1
      DO 10 I=2,NIM1

C
C---------- CELL FACE DENSITIES ------------------------------------
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
C     NON-MUSCL APPROACH
C
      CE=DENE*CX(I  ,J,K)
      CW=DENW*CX(I-1,J,K)
      CN=DENN*CY(I,J  ,K)
      CS=DENS*CY(I,J-1,K)
      CT=DENT*CZ(I,J,K  )
      CB=DENB*CZ(I,J,K-1)
C
      CEPLUS=0.5*(1.+SIGN(1.,CE))
      CEMINS=0.5*(1.-SIGN(1.,CE))
      CWPLUS=0.5*(1.+SIGN(1.,CW))
      CWMINS=0.5*(1.-SIGN(1.,CW))
      CNPLUS=0.5*(1.+SIGN(1.,CN))
      CNMINS=0.5*(1.-SIGN(1.,CN))
      CSPLUS=0.5*(1.+SIGN(1.,CS))
      CSMINS=0.5*(1.-SIGN(1.,CS))
      CTPLUS=0.5*(1.+SIGN(1.,CT))
      CTMINS=0.5*(1.-SIGN(1.,CT))
      CBPLUS=0.5*(1.+SIGN(1.,CB))
      CBMINS=0.5*(1.-SIGN(1.,CB))
C
C      DIP3O2= V(  I+2,J,K  )-V(  I+1,J,K  )
C      DIP1O2= V(  I+1,J,K  )-V(  I  ,J,K  )
C      DIM1O2= V(  I  ,J,K  )-V(  I-1,J,K  )
C      DIM3O2= V(  I-1,J,K  )-V(  I-2,J,K  )
C
C      DJP3O2= V(  I,J+2,K  )-V(  I,J+1,K  )
C      DJP1O2= V(  I,J+1,K  )-V(  I,J  ,K  )
C      DJM1O2= V(  I,J  ,K  )-V(  I,J-1,K  )
C      DJM3O2= V(  I,J-1,K  )-V(  I,J-2,K  )
C
C      DKP3O2= V(  I,  J,K+2)-V(  I,  J,K+1)
C      DKP1O2= V(  I,  J,K+1)-V(  I,  J,K  )
C      DKM1O2= V(  I,  J,K  )-V(  I,  J,K-1)
C      DKM3O2= V(  I,  J,K-1)-V(  I,  J,K-2)
C
C     B.C.
C
C      IF(I.EQ.2) THEN
C      DIM1O2=2.*DIM1O2
C      DIM3O2=0.
C      END IF
C
C      IF(J.EQ.2) THEN
C      DJM1O2=2.*DJM1O2
C      DJM3O2=0.
C      END IF
C
C      IF(K.EQ.2) THEN
C      DKM1O2=2.*DKM1O2
C      DKM3O2=0.
C      END IF
C
C      IF(I.EQ.NIM1) THEN
C      DIP1O2=2.*DIP1O2
C      DIP3O2=0.
C      END IF
C
C      IF(J.EQ.NJM1) THEN
C      DJP1O2=2.*DJP1O2
C      DJP3O2=0.
C      END IF
C
C      IF(K.EQ.NKM1) THEN
C      DKP1O2=2.*DKP1O2
C      DKP3O2=0.
C      END IF
C
C       SUE=
C     1 CE*CEPLUS*
C     1 ( (1.-ETA)*DIM1O2+(1.+ETA)*DIP1O2 )
C     1-CE*CEMINS*
C     1 ( (1.-ETA)*DIP3O2+(1.+ETA)*DIP1O2 )
C
C       SUW=
C     1 CW*CWPLUS*
C     1 ( (1.-ETA)*DIM3O2+(1.+ETA)*DIM1O2 )
C     1-CW*CWMINS*
C     1 ( (1.-ETA)*DIP1O2+(1.+ETA)*DIM1O2 )
C
C       SUN=
C     1 CN*CNPLUS*
C     1 ( (1.-ETA)*DJM1O2+(1.+ETA)*DJP1O2 )
C     1-CN*CNMINS*
C     1 ( (1.-ETA)*DJP3O2+(1.+ETA)*DJP1O2 )
C
C       SUS=
C     1 CS*CSPLUS*
C     1 ( (1.-ETA)*DJM3O2+(1.+ETA)*DJM1O2 )
C     1-CS*CSMINS*
C     1 ( (1.-ETA)*DJP1O2+(1.+ETA)*DJM1O2 )
C
C       SUT=
C     1 CT*CTPLUS*
C     1 ( (1.-ETA)*DKM1O2+(1.+ETA)*DKP1O2 )
C     1-CT*CTMINS*
C     1 ( (1.-ETA)*DKP3O2+(1.+ETA)*DKP1O2 )
C
C       SUB=
C     1 CB*CBPLUS*
C     1 ( (1.-ETA)*DKM3O2+(1.+ETA)*DKM1O2 )
C     1-CB*CBMINS*
C     1 ( (1.-ETA)*DKP1O2+(1.+ETA)*DKM1O2 )
C
C      IF(I.EQ.NIM1) SUE=0.
C      IF(J.EQ.NJM1) SUN=0.
C      IF(K.EQ.NKM1) SUT=0.
C      IF(I.EQ.2   ) SUW=0.
C      IF(J.EQ.2   ) SUS=0.
C      IF(K.EQ.2   ) SUB=0.
C
C      SU(I,J,K)=SU(I,J,K)-0.25*( SUE-SUW+SUN-SUS+SUT-SUB )*RATIO
C
C------ REYNOLDS STRESS GRADIENT TERMS ------------------------------
C
C      EAST-WEST FACE AREAS
C
      AREXE=(1.-FX(I,J,K))*  XX(I,J,K)+  FX(I,J,K)*  XX(I+1,J,K)
      AREXW=(1.-FX(I-1,J,K))*XX(I-1,J,K)+FX(I-1,J,K)*XX(I,J,K)
      AREYE=(1.-FX(I,J,K))*  XY(I,J,K)+  FX(I,J,K)*  XY(I+1,J,K)
      AREYW=(1.-FX(I-1,J,K))*XY(I-1,J,K)+FX(I-1,J,K)*XY(I,J,K)
      AREZE=(1.-FX(I,J,K))*  XZ(I,J,K)+  FX(I,J,K)*  XZ(I+1,J,K)
      AREZW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)+FX(I-1,J,K)*XZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      VEST=(1-FX(I,J,K))*  V(I,J,K)+  FX(I,J,K)*  V(I+1,J,K)
      VWST=(1-FX(I-1,J,K))*V(I-1,J,K)+FX(I-1,J,K)*V(I,J,K)
      VEEST=VEST+(V(I+1,J,K)-V(I,J,K))
      IF(I.NE.NIM1)
     1 VEEST=(1-FX(I+1,J,K))*V(I+1,J,K)+FX(I+1,J,K)*V(I+2,J,K)
      IF(I.EQ.NIM1)
     1 VEEST=(1-FX(4,J,K))*V(4,J,K)+FX(4,J,K)*V(5,J,K)
      VWWST=VWST-(V(I,J,K)-V(I-1,J,K))
      IF(I.NE.2)
     1 VWWST=(1-FX(I-2,J,K))*V(I-2,J,K)+FX(I-2,J,K)*V(I-1,J,K)
      IF(I.EQ.2)
     1 VWWST=(1-FX(NI-4,J,K))*V(NI-4,J,K)+FX(NI-4,J,K)*V(NI-3,J,K)
C
C      EAST-WEST FACE STRESSES
C
C      FX(1,J,K)=0.5
C      FX(NIM1,J,K)=0.5
C
      UVE=(1.-FX(I,J,K))*  UV(I,J,K)+  FX(I,J,K)*  UV(I+1,J,K)
      UVW=(1.-FX(I-1,J,K))*UV(I-1,J,K)+FX(I-1,J,K)*UV(I,J,K)
      VVE=(1.-FX(I,J,K))*  VV(I,J,K)+  FX(I,J,K)*  VV(I+1,J,K)
      VVW=(1.-FX(I-1,J,K))*VV(I-1,J,K)+FX(I-1,J,K)*VV(I,J,K)
      VWE=(1.-FX(I,J,K))*  VW(I,J,K)+  FX(I,J,K)*  VW(I+1,J,K)
      VWW=(1.-FX(I-1,J,K))*VW(I-1,J,K)+FX(I-1,J,K)*VW(I,J,K)
C
C        EAST-WEST CONTRIBUTION
C
      SUSTW=DENW*(AREXW*UVW+AREYW*VVW+AREZW*VWW)
      SUSTE=DENE*(AREXE*UVE+AREYE*VVE+AREZE*VWE)
C --- RHIE AND CHOW ------------------
C      if(j.eq.2.or.j.eq.njm1.or.k.eq.nkm1)go to 128
C      GO TO 128
      SUSTW=SUSTW+
     1            (1.-FX(I-1,J,K))*(VWST-VWWST)*(
     1     ((XX(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS)+
     1      (XY(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS)+
     1      (XZ(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS))/VOL(I-1,J,K))
     1          + FX(I-1,J,K)*(VEST-VWST)*(
     1     ((XX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTE=SUSTE+
     1            (1.-FX(I,J,K))*(VEST-VWST)*(
     1     ((XX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FX(I,J,K)*(VEEST-VEST)*(
     1     ((XX(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS)+
     1      (XY(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS)+
     1      (XZ(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS))/VOL(I+1,J,K))
C
  128 continue
C      FX(1,J,K)=0.
C      FX(NIM1,J,K)=1.
C
C      NORTH-SOUTH FACE AREAS
C
      AREXN=(1.-FY(I,J,K))*  YX(I,J,K)+  FY(I,J,K)*  YX(I,J+1,K)
      AREXS=(1.-FY(I,J-1,K))*YX(I,J-1,K)+FX(I,J-1,K)*YX(I,J,K)
      AREYN=(1.-FY(I,J,K))*  YY(I,J,K)+  FY(I,J,K)*  YY(I,J+1,K)
      AREYS=(1.-FY(I,J-1,K))*YY(I,J-1,K)+FY(I,J-1,K)*YY(I,J,K)
      AREZN=(1.-FY(I,J,K))*  YZ(I,J,K)+  FY(I,J,K)*  YZ(I,J+1,K)
      AREZS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)+FY(I,J-1,K)*YZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      VNRT=(1-FY(I,J,K))*  V(I,J,K)+  FY(I,J,K)*  V(I,J+1,K)
      VSTH=(1-FY(I,J-1,K))*V(I,J-1,K)+FY(I,J-1,K)*V(I,J,K)
      VNNRT=VNRT+(V(I,J+1,K)-V(I,J,K))
      IF(J.NE.NJM1)
     1 VNNRT=(1-FY(I,J+1,K))*V(I,J+1,K)+FY(I,J+1,K)*V(I,J+2,K)
      VSSTH=VSTH-(V(I,J,K)-V(I,J-1,K))
      IF(J.NE.2)
     1 VSSTH=(1-FY(I,J-2,K))*V(I,J-2,K)+FY(I,J-2,K)*V(I,J-1,K)
C
C      NORTH-SOUTH FACE STRESSES
C
C      FY(I,1,K)=.5
C      FY(I,NJM1,K)=0.5
C
      UVN=(1.-FY(I,J,K))*  UV(I,J,K)+  FY(I,J,K)*  UV(I,J+1,K)
      UVS=(1.-FY(I,J-1,K))*UV(I,J-1,K)+FY(I,J-1,K)*UV(I,J,K)
      VVN=(1.-FY(I,J,K))*  VV(I,J,K)+  FY(I,J,K)*  VV(I,J+1,K)
      VVS=(1.-FY(I,J-1,K))*VV(I,J-1,K)+FY(I,J-1,K)*VV(I,J,K)
      VWN=(1.-FY(I,J,K))*  VW(I,J,K)+  FY(I,J,K)*  VW(I,J+1,K)
      VWS=(1.-FY(I,J-1,K))*VW(I,J-1,K)+FY(I,J-1,K)*VW(I,J,K)
C
C        NORTH-SOUTH CONTRIBUTION
C
      SUSTS=DENS*(AREXS*UVS+AREYS*VVS+AREZS*VWS)
      SUSTN=DENN*(AREXN*UVN+AREYN*VVN+AREZN*VWN)
C --- RHIE AND CHOW ------------------
C      if(j.eq.2.or.j.eq.njm1.or.k.eq.nkm1)go to 129
C      GO TO 129
      SUSTS=SUSTS+
     1            (1.-FY(I,J-1,K))*(VSTH-VSSTH)*(
     1     ((YX(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS)+
     1      (YY(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS)+
     1      (YZ(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS))/VOL(I,J-1,K))
     1          + FY(I,J-1,K)*(VNRT-VSTH)*(
     1     ((YX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTN=SUSTN+
     1            (1.-FY(I,J,K))*(VNRT-VSTH)*(
     1     ((YX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FY(I,J,K)*(VNNRT-VNRT)*(
     1     ((YX(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS)+
     1      (YY(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS)+
     1      (YZ(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS))/VOL(I,J+1,K))
C      IF(J.EQ.2)SUSTS=SUSTN
C      IF(J.EQ.NJM1)SUSTN=SUSTS
C
c      SUSTS=DENS*(AREXS*UVS+AREYS*VVS+AREZS*VWS)
c      SUSTN=DENN*(AREXN*UVN+AREYN*VVN+AREZN*VWN)
  129 continue
C      FY(I,1,K)=0.
C      FY(I,NJM1,K)=1.
C
C      TOP-BOTTOM FACE AREAS
C
      AREXT=(1.-FZ(I,J,K))*  ZX(I,J,K)+  FZ(I,J,K)*  ZX(I,J,K+1)
      AREXB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)+FZ(I,J,K-1)*ZX(I,J,K)
      AREYT=(1.-FZ(I,J,K))*  ZY(I,J,K)+  FZ(I,J,K)*  ZY(I,J,K+1)
      AREYB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)+FZ(I,J,K-1)*ZY(I,J,K)
      AREZT=(1.-FZ(I,J,K))*  ZZ(I,J,K)+  FZ(I,J,K)*  ZZ(I,J,K+1)
      AREZB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)+FZ(I,J,K-1)*ZZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      VTOP=(1-FZ(I,J,K))*  V(I,J,K)+  FZ(I,J,K)*  V(I,J,K+1)
      VBOT=(1-FZ(I,J,K-1))*V(I,J,K-1)+FZ(I,J,K-1)*V(I,J,K)
      VTTOP=VTOP+(V(I,J,K+1)-V(I,J,K))
      IF(K.NE.NKM1)
     1 VTTOP=(1-FZ(I,J,K+1))*V(I,J,K+1)+FZ(I,J,K+1)*V(I,J,K+2)
      VBBOT=VBOT-(V(I,J,K)-V(I,J,K-1))
      IF(K.NE.2)
     1 VBBOT=(1-FZ(I,J,K-2))*V(I,J,K-2)+FZ(I,J,K-2)*V(I,J,K-1)
C
C      TOP-BOTTOM FACE STRESSES
C
C      FZ(I,J,1)=.5
c      FZ(I,J,NKM1)=.5
C
      UVT=(1.-FZ(I,J,K))*  UV(I,J,K)+  FZ(I,J,K)*  UV(I,J,K+1)
      UVB=(1.-FZ(I,J,K-1))*UV(I,J,K-1)+FZ(I,J,K-1)*UV(I,J,K)
      VVT=(1.-FZ(I,J,K))  *VV(I,J,K)+  FZ(I,J,K)*  VV(I,J,K+1)
      VVB=(1.-FZ(I,J,K-1))*VV(I,J,K-1)+FZ(I,J,K-1)*VV(I,J,K)
      VWT=(1.-FZ(I,J,K))*  VW(I,J,K)+  FZ(I,J,K)*  VW(I,J,K+1)
      VWB=(1.-FZ(I,J,K-1))*VW(I,J,K-1)+FZ(I,J,K-1)*VW(I,J,K)
C
C        TOP-BOTTOM CONTRIBUTION
C
      SUSTB=DENB*(AREXB*UVB+AREYB*VVB+AREZB*VWB)
      SUSTT=DENT*(AREXT*UVT+AREYT*VVT+AREZT*VWT)
C --- RHIE AND CHOW ------------------
C      if(k.eq.2.or.k.eq.nkm1.or.j.eq.2.or.j.eq.njm1)go to 130
C      GO TO 130
      SUSTB=SUSTB+
     1            (1.-FZ(I,J,K-1))*(VBOT-VBBOT)*(
     1     ((ZX(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS)+
     1      (ZY(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS)+
     1      (ZZ(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS))/VOL(I,J,K-1))
     1          + FZ(I,J,K-1)*(VTOP-VBOT)*(
     1     ((ZX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTT=SUSTT+
     1            (1.-FZ(I,J,K))*(VTOP-VBOT)*(
     1     ((ZX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FZ(I,J,K)*(VTTOP-VTOP)*(
     1     ((ZX(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS)+
     1      (ZY(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS)+
     1      (ZZ(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS))/VOL(I,J,K+1))
C
  130 continue
C      FZ(I,J,1)=0.
c      FZ(I,J,NKM1)=1.
C
C   OVERALL STRESS GRADIENT TERM
C
      SU(I,J,K)=SU(I,J,K)+FLOAT(1-IEVM)*
     1          (SUSTW-SUSTE+SUSTS-SUSTN+SUSTB-SUSTT)
   10 CONTINUE
C
C -------- PRESSURE GRADIENT SOURCE TERM  ----------------------------
C
C      DO  K=2,NKM1
C
C      DO  I=2,NIR
C      P(I,NJR,K)=1.*P(I,NJR+1,K)+0.*P(I,NJR+2,K)
C      P(I,NJRN+1,K)=1.*P(I,NJRN,K)+0.*P(I,NJRN-1,K)
C      END DO
C      DO 85 I=NIRMS+1,NIRMF
C      P(I,NJR,K)=1.*P(I,NJR+1,K)-0.*P(I,NJR+2,K)
C      P(I,NJRN+1,K)=1.*P(I,NJRN,K)-0.*P(I,NJRN-1,K)
C   85 CONTINUE
C      DO 85 I=NI-2,NI-1
C      P(I,NJR,K)=1.*P(I,NJR+1,K)+0.*P(I,NJR+2,K)
C      P(I,NJRN+1,K)=1.*P(I,NJRN,K)+0.*P(I,NJRN-1,K)
C   85 CONTINUE
C      END DO
C
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
      PRESE=P(I,J,K)*(1.-FX(I,J,K))+P(I+1,J,K)*FX(I,J,K)
      PRESW=P(I-1,J,K)*(1.-FX(I-1,J,K))+P(I,J,K)*FX(I-1,J,K)
      PRESN=P(I,J,K)*(1.-FY(I,J,K))+P(I,J+1,K)*FY(I,J,K)
      PRESS=P(I,J-1,K)*(1.-FY(I,J-1,K))+P(I,J,K)*FY(I,J-1,K)
      PREST=P(I,J,K)*(1.-FZ(I,J,K))+P(I,J,K+1)*FZ(I,J,K)
      PRESB=P(I,J,K-1)*(1.-FZ(I,J,K-1))+P(I,J,K)*FZ(I,J,K-1)
C
      DV1(I,J,K)=XY(I,J,K)
      DV2(I,J,K)=YY(I,J,K)
      DV3(I,J,K)=ZY(I,J,K)
C
      SU(I,J,K)=SU(I,J,K)+DV1(I,J,K)*(PRESW-PRESE)
     1 +DV2(I,J,K)*(PRESS-PRESN)+DV3(I,J,K)*(PRESB-PREST)
C
 50   CONTINUE
C
C------ DIFFUSION  SOURCE TERMS ------------------------------------
C
      DO 60 K=2,NKM1
      DO 60 J=2,NJM1
      DO 60 I=2,NIM1
C
C ------ CELL FACE VISCOCITIES ------------------------------------
C
       VISE=VIS(I,J,K)*(1.-FX(I,J,K))+VIS(I+1,J,K)*FX(I,J,K)
       VISW=VIS(I-1,J,K)*(1.-FX(I-1,J,K))+VIS(I,J,K)*FX(I-1,J,K)
       VISN=VIS(I,J,K)*(1.-FY(I,J,K))+VIS(I,J+1,K)*FY(I,J,K)
       VISS=VIS(I,J-1,K)*(1.-FY(I,J-1,K))+VIS(I,J,K)*FY(I,J-1,K)
       VIST=VIS(I,J,K)*(1.-FZ(I,J,K))+VIS(I,J,K+1)*FZ(I,J,K)
       VISB=VIS(I,J,K-1)*(1.-FZ(I,J,K-1))+VIS(I,J,K)*FZ(I,J,K-1)
C
       VISE=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISE
       VISW=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISW
       VISN=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISN
       VISS=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISS
       VIST=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VIST
       VISB=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISB
C
C ------ CONTROL VOLUMES FOR DUMMY (STAGGERED) VELOCITY NODES ----
C
       VOLE=VOL(I,J,K)*(1.-FX(I,J,K))+VOL(I+1,J,K)*FX(I,J,K)
       VOLW=VOL(I-1,J,K)*(1.-FX(I-1,J,K))+VOL(I,J,K)*FX(I-1,J,K)
       VOLN=VOL(I,J,K)*(1.-FY(I,J,K))+VOL(I,J+1,K)*FY(I,J,K)
       VOLS=VOL(I,J-1,K)*(1.-FY(I,J-1,K))+VOL(I,J,K)*FY(I,J-1,K)
       VOLT=VOL(I,J,K)*(1.-FZ(I,J,K))+VOL(I,J,K+1)*FZ(I,J,K)
       VOLB=VOL(I,J,K-1)*(1.-FZ(I,J,K-1))+VOL(I,J,K)*FZ(I,J,K-1)
C
C
C     CALCULATE THE CORNER VALUES FOR U,V,W
C
      UTNE=
     1 U(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+U(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      UTNW=
     1 U(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+U(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      UBNE=
     1 U(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+U(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      UBNW=
     1 U(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+U(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+U(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      UTSE=
     1 U(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+U(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      UTSW=
     1 U(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+U(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      UBSE=
     1 U(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+U(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+U(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      UBSW=
     1 U(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+U(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      VTNE=
     1 V(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+V(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      VTNW=
     1 V(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+V(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      VBNE=
     1 V(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+V(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      VBNW=
     1 V(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+V(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+V(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      VTSE=
     1 V(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+V(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      VTSW=
     1 V(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+V(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      VBSE=
     1 V(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+V(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+V(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      VBSW=
     1 V(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+V(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      WTNE=
     1 W(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+W(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      WTNW=
     1 W(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+W(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      WBNE=
     1 W(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+W(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      WBNW=
     1 W(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+W(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+W(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      WTSE=
     1 W(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+W(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      WTSW=
     1 W(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+W(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      WBSE=
     1 W(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+W(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+W(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      WBSW=
     1 W(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+W(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
C ---- CROSS-DIFFUSION SOURCE TERMS  ------------------------------
C
C         EAST-WEST FACE CROSS-DIFFUSION
C
      Q12E=(XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*YX(I+1,J,K)+XY(I+1,J,K)*YY(I+1,J,K)+
     1 XZ(I+1,J,K)*YZ(I+1,J,K))*FX(I,J,K)
C
      Q12W=(XX(I-1,J,K)*YX(I-1,J,K)+XY(I-1,J,K)*YY(I-1,J,K)+
     1 XZ(I-1,J,K)*YZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*FX(I-1,J,K)
C
      Q13E=(XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*ZX(I+1,J,K)+XY(I+1,J,K)*ZY(I+1,J,K)+
     1 XZ(I+1,J,K)*ZZ(I+1,J,K))*FX(I,J,K)
C
      Q13W=(XX(I-1,J,K)*ZX(I-1,J,K)+XY(I-1,J,K)*ZY(I-1,J,K)+
     1 XZ(I-1,J,K)*ZZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*FX(I-1,J,K)
C
      SCX=0.5*VISE*Q12E/VOLE*(VTNE+VBNE-VTSE-VBSE)
     1   +0.5*VISE*Q13E/VOLE*(VTNE+VTSE-VBNE-VBSE)
     1   -0.5*VISW*Q12W/VOLW*(VTNW+VBNW-VTSW-VBSW)
     1   -0.5*VISW*Q13W/VOLW*(VTNW+VTSW-VBNW-VBSW)
C
C        NORTH-SOUTH FACE CROSS-DIFFUSION
C
      Q21N=(YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*XX(I,J+1,K)+YY(I,J+1,K)*XY(I,J+1,K)+
     1 YZ(I,J+1,K)*XZ(I,J+1,K))*FY(I,J,K)
C
      Q21S=(YX(I,J-1,K)*XX(I,J-1,K)+YY(I,J-1,K)*XY(I,J-1,K)+
     1 YZ(I,J-1,K)*XZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*FY(I,J-1,K)
C
      Q23N=(YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*ZX(I,J+1,K)+YY(I,J+1,K)*ZY(I,J+1,K)+
     1 YZ(I,J+1,K)*ZZ(I,J+1,K))*FY(I,J,K)
C
      Q23S=(YX(I,J-1,K)*ZX(I,J-1,K)+YY(I,J-1,K)*ZY(I,J-1,K)+
     1 YZ(I,J-1,K)*ZZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*FY(I,J-1,K)
C
      SCY=0.5*VISN*Q21N/VOLN*(VTNE+VBNE-VTNW-VBNW)
     1   +0.5*VISN*Q23N/VOLN*(VTNE+VTNW-VBNE-VBNW)
     1   -0.5*VISS*Q21S/VOLS*(VTSE+VBSE-VTSW-VBSW)
     1   -0.5*VISS*Q23S/VOLS*(VTSE+VTSW-VBSE-VBSW)
C
C     TOP-BOTTOM FACE CROSS-DIFUSION
C
      Q31T=(ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*XX(I,J,K+1)+ZY(I,J,K+1)*XY(I,J,K+1)+
     1 ZZ(I,J,K+1)*XZ(I,J,K+1))*FZ(I,J,K)
C
      Q31B=(ZX(I,J,K-1)*XX(I,J,K-1)+ZY(I,J,K-1)*XY(I,J,K-1)+
     1 ZZ(I,J,K-1)*XZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*FZ(I,J,K-1)
C
      Q32T=(ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*YX(I,J,K+1)+ZY(I,J,K+1)*YY(I,J,K+1)+
     1 ZZ(I,J,K+1)*YZ(I,J,K+1))*FZ(I,J,K)
C
      Q32B=(ZX(I,J,K-1)*YX(I,J,K-1)+ZY(I,J,K-1)*YY(I,J,K-1)+
     1 ZZ(I,J,K-1)*YZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*FZ(I,J,K-1)
C
      SCZ=0.5*VIST*Q31T/VOLT*(VTNE+VTSE-VTNW-VTSW)
     1   +0.5*VIST*Q32T/VOLT*(VTNE+VTNW-VTSE-VTSW)
     1   -0.5*VISB*Q31B/VOLB*(VBNE+VBSE-VBNW-VBSW)
     1   -0.5*VISB*Q32B/VOLB*(VBNE+VBNW-VBSE-VBSW)
C
      SU(I,J,K)=SU(I,J,K)+SCX+SCY+SCZ
C
C------ VARIABLE VISCOSITY SOURCE TERMS  --------------------------
C
C        EAST WEST U,V,W GRADIENT TERMS
C
      XXXYE=(1.-FX(I,J,K))*XX(I,J,K)*XY(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*XY(I+1,J,K)
      XYXYE=(1.-FX(I,J,K))*XY(I,J,K)*XY(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*XY(I+1,J,K)
      XZXYE=(1.-FX(I,J,K))*XZ(I,J,K)*XY(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*XY(I+1,J,K)
      XXYYE=(1.-FX(I,J,K))*XX(I,J,K)*YY(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*YY(I+1,J,K)
      XYYYE=(1.-FX(I,J,K))*XY(I,J,K)*YY(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*YY(I+1,J,K)
      XZYYE=(1.-FX(I,J,K))*XZ(I,J,K)*YY(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*YY(I+1,J,K)
      XXZYE=(1.-FX(I,J,K))*XX(I,J,K)*ZY(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*ZY(I+1,J,K)
      XYZYE=(1.-FX(I,J,K))*XY(I,J,K)*ZY(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*ZY(I+1,J,K)
      XZZYE=(1.-FX(I,J,K))*XZ(I,J,K)*ZY(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*ZY(I+1,J,K)
C
      XXXYW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*XY(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*XY(I,J,K)
      XYXYW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*XY(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*XY(I,J,K)
      XZXYW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*XY(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*XY(I,J,K)
      XXYYW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*YY(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*YY(I,J,K)
      XYYYW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*YY(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*YY(I,J,K)
      XZYYW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*YY(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*YY(I,J,K)
      XXZYW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*ZY(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*ZY(I,J,K)
      XYZYW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*ZY(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*ZY(I,J,K)
      XZZYW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*ZY(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*ZY(I,J,K)
C
      SC1=VISE/VOLE*( (U(I+1,J,K)-U(I,J,K))*XXXYE+(V(I+1,J,K)-V(I,J,K))
     1 *XYXYE+(W(I+1,J,K)-W(I,J,K))*XZXYE )
     1   -VISW/VOLW*( (U(I,J,K)-U(I-1,J,K))*XXXYW+(V(I,J,K)-V(I-1,J,K))
     1 *XYXYW+(W(I,J,K)-W(I-1,J,K))*XZXYW )
C
      SC2=VISE/VOLE*( (UTNE+UBNE-UTSE-UBSE)*XXYYE+(VTNE+VBNE-VTSE-VBSE)
     1 *XYYYE+(WTNE+WBNE-WTSE-WBSE)*XZYYE )
     1   -VISW/VOLW*( (UTNW+UBNW-UTSW-UBSW)*XXYYW+(VTNW+VBNW-VTSW-VBSW)
     1 *XYYYW+(WTNW+WBNW-WTSW-WBSW)*XZYYW )
      SC2=0.5*SC2
C
      SC3=VISE/VOLE*( (UTNE+UTSE-UBNE-UBSE)*XXZYE+(VTNE+VTSE-VBNE-VBSE)
     1 *XYZYE+(WTNE+WTSE-WBNE-WBSE)*XZZYE )
     1   -VISW/VOLW*( (UTNW+UTSW-UBNW-UBSW)*XXZYW+(VTNW+VTSW-VBNW-VBSW)
     1 *XYZYW+(WTNW+WTSW-WBNW-WBSW)*XZZYW )
      SC3=0.5*SC3
C
C      NORTH-SOUTH U,V,W GRADIENT TERMS
C
      YXXYN=(1.-FY(I,J,K))*YX(I,J,K)*XY(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*XY(I,J+1,K)
      YYXYN=(1.-FY(I,J,K))*YY(I,J,K)*XY(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*XY(I,J+1,K)
      YZXYN=(1.-FY(I,J,K))*YZ(I,J,K)*XY(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*XY(I,J+1,K)
      YXYYN=(1.-FY(I,J,K))*YX(I,J,K)*YY(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*YY(I,J+1,K)
      YYYYN=(1.-FY(I,J,K))*YY(I,J,K)*YY(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*YY(I,J+1,K)
      YZYYN=(1.-FY(I,J,K))*YZ(I,J,K)*YY(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*YY(I,J+1,K)
      YXZYN=(1.-FY(I,J,K))*YX(I,J,K)*ZY(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*ZY(I,J+1,K)
      YYZYN=(1.-FY(I,J,K))*YY(I,J,K)*ZY(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*ZY(I,J+1,K)
      YZZYN=(1.-FY(I,J,K))*YZ(I,J,K)*ZY(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*ZY(I,J+1,K)
C
      YXXYS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*XY(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*XY(I,J,K)
      YYXYS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*XY(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*XY(I,J,K)
      YZXYS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*XY(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*XY(I,J,K)
      YXYYS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*YY(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*YY(I,J,K)
      YYYYS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*YY(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*YY(I,J,K)
      YZYYS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*YY(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*YY(I,J,K)
      YXZYS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*ZY(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*ZY(I,J,K)
      YYZYS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*ZY(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*ZY(I,J,K)
      YZZYS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*ZY(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*ZY(I,J,K)
C
      SC4=VISN/VOLN*( (UTNE+UBNE-UTNW-UBNW)*YXXYN+(VTNE+VBNE-VTNW-VBNW)
     1 *YYXYN+(WTNE+WBNE-WTNW-WBNW)*YZXYN )
     1   -VISS/VOLS*( (UTSE+UBSE-UTSW-UBSW)*YXXYS+(VTSE+VBSE-VTSW-VBSW)
     1 *YYXYS+(WTSE+WBSE-WTSW-WBSW)*YZXYS )
      SC4=0.5*SC4
C
      SC5=VISN/VOLN*( (U(I,J+1,K)-U(I,J,K))*YXYYN+(V(I,J+1,K)-V(I,J,K))
     1 *YYYYN+(W(I,J+1,K)-W(I,J,K))*YZYYN )
     1   -VISS/VOLS*( (U(I,J,K)-U(I,J-1,K))*YXYYS+(V(I,J,K)-V(I,J-1,K))
     1 *YYYYS+(W(I,J,K)-W(I,J-1,K))*YZYYS )
C
      SC6=VISN/VOLN*( (UTNE+UTNW-UBNE-UBNW)*YXZYN+(VTNE+VTNW-VBNE-VBNW)
     1 *YYZYN+(WTNE+WTNW-WBNE-WBNW)*YZZYN )
     1   -VISS/VOLS*( (UTSE+UTSW-UBSE-UBSW)*YXZYS+(VTSE+VTSW-VBSE-VBSW)
     1 *YYZYS+(WTSE+WTSW-WBSE-WBSW)*YZZYS )
      SC6=0.5*SC6
C
C    TOP-BOTTOM U,V,W GRADIENT TERMS
C
      ZXXYT=(1.-FZ(I,J,K))*ZX(I,J,K)*XY(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*XY(I,J,K+1)
      ZYXYT=(1.-FZ(I,J,K))*ZY(I,J,K)*XY(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*XY(I,J,K+1)
      ZZXYT=(1.-FZ(I,J,K))*ZZ(I,J,K)*XY(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*XY(I,J,K+1)
      ZXYYT=(1.-FZ(I,J,K))*ZX(I,J,K)*YY(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*YY(I,J,K+1)
      ZYYYT=(1.-FZ(I,J,K))*ZY(I,J,K)*YY(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*YY(I,J,K+1)
      ZZYYT=(1.-FZ(I,J,K))*ZZ(I,J,K)*YY(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*YY(I,J,K+1)
      ZXZYT=(1.-FZ(I,J,K))*ZX(I,J,K)*ZY(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*ZY(I,J,K+1)
      ZYZYT=(1.-FZ(I,J,K))*ZY(I,J,K)*ZY(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*ZY(I,J,K+1)
      ZZZYT=(1.-FZ(I,J,K))*ZZ(I,J,K)*ZY(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*ZY(I,J,K+1)
C
      ZXXYB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*XY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*XY(I,J,K)
      ZYXYB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*XY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*XY(I,J,K)
      ZZXYB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*XY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*XY(I,J,K)
      ZXYYB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*YY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*YY(I,J,K)
      ZYYYB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*YY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*YY(I,J,K)
      ZZYYB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*YY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*YY(I,J,K)
      ZXZYB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*ZY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*ZY(I,J,K)
      ZYZYB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*ZY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*ZY(I,J,K)
      ZZZYB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*ZY(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*ZY(I,J,K)
C
      SC7=VIST/VOLT*( (UTNE+UTSE-UTNW-UTSW)*ZXXYT+(VTNE+VTSE-VTNW-VTSW)
     1 *ZYXYT+(WTNE+WTSE-WTNW-WTSW)*ZZXYT )
     1   -VISB/VOLB*( (UBNE+UBSE-UBNW-UBSW)*ZXXYB+(VBNE+VBSE-VBNW-VBSW)
     1 *ZYXYB+(WBNE+WBSE-WBNW-WBSW)*ZZXYB )
      SC7=0.5*SC7
C
      SC8=VIST/VOLT*( (UTNE+UTNW-UTSE-UTSW)*ZXYYT+(VTNE+VTNW-VTSE-VTSW)
     1 *ZYYYT+(WTNE+WTNW-WTSE-WTSW)*ZZYYT )
     1   -VISB/VOLB*( (UBNE+UBNW-UBSE-UBSW)*ZXYYB+(VBNE+VBNW-VBSE-VBSW)
     1 *ZYYYB+(WBNE+WBNW-WBSE-WBSW)*ZZYYB )
      SC8=0.5*SC8
C
      SC9=VIST/VOLT*( (U(I,J,K+1)-U(I,J,K))*ZXZYT+(V(I,J,K+1)-V(I,J,K))
     1 *ZYZYT+(W(I,J,K+1)-W(I,J,K))*ZZZYT )
     1   -VISB/VOLB*( (U(I,J,K)-U(I,J,K-1))*ZXZYB+(V(I,J,K)-V(I,J,K-1))
     1 *ZYZYB+(W(I,J,K)-W(I,J,K-1))*ZZZYB )
C
      SU(I,J,K)=SU(I,J,K)+SC1+SC2+SC3+SC4+SC5+SC6+SC7+SC8+SC9
C
   60 CONTINUE
C---------------------------------------------------------------------
C          ROTATION TERMS                                            *
C---------------------------------------------------------------------
      DO 70 K=2,NKM1
      DO 70 J=2,NJM1
      DO 70 I=2,NIM1
C----  CORIOLIS TERMS -----------------------------------------------
      SU(I,J,K)=SU(I,J,K)+2.*DEN(I,J,K)*VOL(I,J,K)*
     +                    (OMX*W(I,J,K)-OMZ*U(I,J,K))
C----  CENTRIFUGAL TERMS --------------------------------------------
      HXX=HX+XC(I,J,K)
      HYY=HY+YC(I,J,K)
      HZZ=HZ+ZC(I,J,K)
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*
     +                    (HYY*(OMX*OMX+OMZ*OMZ)-OMY*(HXX*OMX+HZZ*OMZ))
   70 CONTINUE
C------ BOUNDARY CONDITIONS   ---------------------------------------
      CALL MODV
C
      RESORV=0.
      DO 30 K=2,NKM1
      DO 30 J=2,NJM1
      DO 30 I=2,NIM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*V(I+1,J,K)+AW(I,J,K)*V(I-1,J,K)
     1     +AN(I,J,K)*V(I,J+1,K)+AS(I,J,K)*V(I,J-1,K)
     1     +AT(I,J,K)*V(I,J,K+1)+AB(I,J,K)*V(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*V(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
      RESORV=RESORV+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFV
      SU(I,J,K)=SU(I,J,K)+(1.-URFV)*AP(I,J,K)*V(I,J,K)
      DV1(I,J,K)=DV1(I,J,K)/AP(I,J,K)
      DV2(I,J,K)=DV2(I,J,K)/AP(I,J,K)
      DV3(I,J,K)=DV3(I,J,K)/AP(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPV
      CALL LISOLV(4,2,2,NI,NJ,NK,V)
C      CALL LISOLV(3,NJR+1,2,NI,NJRN+1,NK,V)
C      CALL LISOLV(NIR+1,2,2,NI-2,NJ,NK,V)
 40   CONTINUE
      RETURN
      END
C
      SUBROUTINE CALCW
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
      PR=1.
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*W(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFICIENTS  ----------------------------
C
      CALL ADVECT(1.,PR)
C
      CALL BQUICK(0,1.,W)
C
      DO 10 K=2,NKM1
      DO 10 J=2,NJM1
      DO 10 I=2,NIM1
C
C---------- CELL FACE DENSITIES ------------------------------------
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
      CE=CX(I  ,J,K)
      CW=CX(I-1,J,K)
      CN=CY(I,J  ,K)
      CS=CY(I,J-1,K)
      CT=CZ(I,J,K  )
      CB=CZ(I,J,K-1)
C
C     NON-MUSCL APPROACH
C
      CEPLUS=0.5*(1.+SIGN(1.,CE))
      CEMINS=0.5*(1.-SIGN(1.,CE))
      CWPLUS=0.5*(1.+SIGN(1.,CW))
      CWMINS=0.5*(1.-SIGN(1.,CW))
      CNPLUS=0.5*(1.+SIGN(1.,CN))
      CNMINS=0.5*(1.-SIGN(1.,CN))
      CSPLUS=0.5*(1.+SIGN(1.,CS))
      CSMINS=0.5*(1.-SIGN(1.,CS))
      CTPLUS=0.5*(1.+SIGN(1.,CT))
      CTMINS=0.5*(1.-SIGN(1.,CT))
      CBPLUS=0.5*(1.+SIGN(1.,CB))
      CBMINS=0.5*(1.-SIGN(1.,CB))
C
C      DIP3O2= W(  I+2,J,K  )-W(  I+1,J,K  )
C      DIP1O2= W(  I+1,J,K  )-W(  I  ,J,K  )
C      DIM1O2= W(  I  ,J,K  )-W(  I-1,J,K  )
C      DIM3O2= W(  I-1,J,K  )-W(  I-2,J,K  )
C
C      DJP3O2= W(  I,J+2,K  )-W(  I,J+1,K  )
C      DJP1O2= W(  I,J+1,K  )-W(  I,J  ,K  )
C      DJM1O2= W(  I,J  ,K  )-W(  I,J-1,K  )
C      DJM3O2= W(  I,J-1,K  )-W(  I,J-2,K  )
C
C      DKP3O2= W(  I,  J,K+2)-W(  I,  J,K+1)
C      DKP1O2= W(  I,  J,K+1)-W(  I,  J,K  )
C      DKM1O2= W(  I,  J,K  )-W(  I,  J,K-1)
C      DKM3O2= W(  I,  J,K-1)-W(  I,  J,K-2)
C
C     B.C.
C
C      IF(I.EQ.2) THEN
C      DIM1O2=2.*DIM1O2
C      DIM3O2=0.
C      END IF
C
C      IF(J.EQ.2) THEN
C      DJM1O2=2.*DJM1O2
C      DJM3O2=0.
C      END IF
C
C      IF(K.EQ.2) THEN
C      DKM1O2=2.*DKM1O2
C      DKM3O2=0.
C      END IF
C
C      IF(I.EQ.NIM1) THEN
C      DIP1O2=2.*DIP1O2
C      DIP3O2=0.
C      END IF
C
C      IF(J.EQ.NJM1) THEN
C      DJP1O2=2.*DJP1O2
C      DJP3O2=0.
C      END IF
C
C      IF(K.EQ.NKM1) THEN
C      DKP1O2=2.*DKP1O2
C      DKP3O2=0.
C      END IF
C
C       SUE=
C     1 CE*CEPLUS*
C     1 ( (1.-ETA)*DIM1O2+(1.+ETA)*DIP1O2 )
C     1-CE*CEMINS*
C     1 ( (1.-ETA)*DIP3O2+(1.+ETA)*DIP1O2 )
C
C       SUW=
C     1 CW*CWPLUS*
C     1 ( (1.-ETA)*DIM3O2+(1.+ETA)*DIM1O2 )
C     1-CW*CWMINS*
C     1 ( (1.-ETA)*DIP1O2+(1.+ETA)*DIM1O2 )
C
C       SUN=
C     1 CN*CNPLUS*
C     1 ( (1.-ETA)*DJM1O2+(1.+ETA)*DJP1O2 )
C     1-CN*CNMINS*
C     1 ( (1.-ETA)*DJP3O2+(1.+ETA)*DJP1O2 )
C
C       SUS=
C     1 CS*CSPLUS*
C     1 ( (1.-ETA)*DJM3O2+(1.+ETA)*DJM1O2 )
C     1-CS*CSMINS*
C     1 ( (1.-ETA)*DJP1O2+(1.+ETA)*DJM1O2 )
C
C       SUT=
C     1 CT*CTPLUS*
C     1 ( (1.-ETA)*DKM1O2+(1.+ETA)*DKP1O2 )
C     1-CT*CTMINS*
C     1 ( (1.-ETA)*DKP3O2+(1.+ETA)*DKP1O2 )
C
C       SUB=
C     1 CB*CBPLUS*
C     1 ( (1.-ETA)*DKM3O2+(1.+ETA)*DKM1O2 )
C     1-CB*CBMINS*
C     1 ( (1.-ETA)*DKP1O2+(1.+ETA)*DKM1O2 )
C
C      IF(I.EQ.NIM1) SUE=0.
C      IF(J.EQ.NJM1) SUN=0.
C      IF(K.EQ.NKM1) SUT=0.
C      IF(I.EQ.2   ) SUW=0.
C      IF(J.EQ.2   ) SUS=0.
C      IF(K.EQ.2   ) SUB=0.
C
C      SU(I,J,K)=SU(I,J,K)-0.25*( SUE-SUW+SUN-SUS+SUT-SUB )*RATIO
C
C------ REYNOLDS STRESS GRADIENT TERMS ------------------------------
C
C      EAST-WEST FACE AREAS
C
      AREXE=(1.-FX(I,J,K))*  XX(I,J,K)+  FX(I,J,K)*  XX(I+1,J,K)
      AREXW=(1.-FX(I-1,J,K))*XX(I-1,J,K)+FX(I-1,J,K)*XX(I,J,K)
      AREYE=(1.-FX(I,J,K))*  XY(I,J,K)+  FX(I,J,K)*  XY(I+1,J,K)
      AREYW=(1.-FX(I-1,J,K))*XY(I-1,J,K)+FX(I-1,J,K)*XY(I,J,K)
      AREZE=(1.-FX(I,J,K))*  XZ(I,J,K)+  FX(I,J,K)*  XZ(I+1,J,K)
      AREZW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)+FX(I-1,J,K)*XZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      WEST=(1-FX(I,J,K))*  W(I,J,K)+  FX(I,J,K)*  W(I+1,J,K)
      WWST=(1-FX(I-1,J,K))*W(I-1,J,K)+FX(I-1,J,K)*W(I,J,K)
      WEEST=WEST+(W(I+1,J,K)-W(I,J,K))
      IF(I.NE.NIM1)
     1 WEEST=(1-FX(I+1,J,K))*W(I+1,J,K)+FX(I+1,J,K)*W(I+2,J,K)
      IF(I.EQ.NIM1)
     1 WEEST=(1-FX(4,J,K))*W(4,J,K)+FX(4,J,K)*W(5,J,K)
      WWWST=WWST-(W(I,J,K)-W(I-1,J,K))
      IF(I.NE.2)
     1 WWWST=(1-FX(I-2,J,K))*W(I-2,J,K)+FX(I-2,J,K)*W(I-1,J,K)
      IF(I.EQ.2)
     1 WWWST=(1-FX(NI-4,J,K))*W(NI-4,J,K)+FX(NI-4,J,K)*W(NI-3,J,K)
C
C      EAST-WEST FACE STRESSES
C
C      FX(1,J,K)=.5
C      FX(NIM1,J,K)=0.5
C
      UWE=(1.-FX(I,J,K))*  UW(I,J,K)+  FX(I,J,K)*  UW(I+1,J,K)
      UWW=(1.-FX(I-1,J,K))*UW(I-1,J,K)+FX(I-1,J,K)*UW(I,J,K)
      VWE=(1.-FX(I,J,K))*  VW(I,J,K)+  FX(I,J,K)*  VW(I+1,J,K)
      VWW=(1.-FX(I-1,J,K))*VW(I-1,J,K)+FX(I-1,J,K)*VW(I,J,K)
      WWE=(1.-FX(I,J,K))*  WW(I,J,K)+  FX(I,J,K)*  WW(I+1,J,K)
      WWW=(1.-FX(I-1,J,K))*WW(I-1,J,K)+FX(I-1,J,K)*WW(I,J,K)
C
C        EAST-WEST CONTRIBUTION
C
      SUSTW=DENW*(AREXW*UWW+AREYW*VWW+AREZW*WWW)
      SUSTE=DENE*(AREXE*UWE+AREYE*VWE+AREZE*WWE)
C --- RHIE AND CHOW ------------------
C      if(j.eq.2.or.j.eq.njm1.or.k.eq.nkm1)go to 128
C      GO TO 128
      SUSTW=SUSTW+
     1            (1.-FX(I-1,J,K))*(WWST-WWWST)*(
     1     ((XX(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS)+
     1      (XY(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS)+
     1      (XZ(I-1,J,K)**2)*(VIS(I-1,J,K)-VISCOS))/VOL(I-1,J,K))
     1          + FX(I-1,J,K)*(WEST-WWST)*(
     1     ((XX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTE=SUSTE+
     1            (1.-FX(I,J,K))*(WEST-WWST)*(
     1     ((XX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (XZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FX(I,J,K)*(WEEST-WEST)*(
     1     ((XX(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS)+
     1      (XY(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS)+
     1      (XZ(I+1,J,K)**2)*(VIS(I+1,J,K)-VISCOS))/VOL(I+1,J,K))
C
  128 continue
C      FX(1,J,K)=0.
C      FX(NIM1,J,K)=1.
C
C      NORTH-SOUTH FACE AREAS
C
      AREXN=(1.-FY(I,J,K))*  YX(I,J,K)+  FY(I,J,K)*  YX(I,J+1,K)
      AREXS=(1.-FY(I,J-1,K))*YX(I,J-1,K)+FX(I,J-1,K)*YX(I,J,K)
      AREYN=(1.-FY(I,J,K))*  YY(I,J,K)+  FY(I,J,K)*  YY(I,J+1,K)
      AREYS=(1.-FY(I,J-1,K))*YY(I,J-1,K)+FY(I,J-1,K)*YY(I,J,K)
      AREZN=(1.-FY(I,J,K))*  YZ(I,J,K)+  FY(I,J,K)*  YZ(I,J+1,K)
      AREZS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)+FY(I,J-1,K)*YZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      WNRT=(1-FY(I,J,K))*  W(I,J,K)+  FY(I,J,K)*  W(I,J+1,K)
      WSTH=(1-FY(I,J-1,K))*W(I,J-1,K)+FY(I,J-1,K)*W(I,J,K)
      WNNRT=WNRT+(W(I,J+1,K)-W(I,J,K))
      IF(J.NE.NJM1)
     1 WNNRT=(1-FY(I,J+1,K))*W(I,J+1,K)+FY(I,J+1,K)*W(I,J+2,K)
      WSSTH=WSTH-(W(I,J,K)-W(I,J-1,K))
      IF(J.NE.2)
     1 WSSTH=(1-FY(I,J-2,K))*W(I,J-2,K)+FY(I,J-2,K)*W(I,J-1,K)
C
C      NORTH-SOUTH FACE STRESSES
C
C      FY(I,1,K)=.5
C      FY(I,NJM1,K)=0.5
C
      UWN=(1.-FY(I,J,K))*  UW(I,J,K)+  FY(I,J,K)*  UW(I,J+1,K)
      UWS=(1.-FY(I,J-1,K))*UW(I,J-1,K)+FY(I,J-1,K)*UW(I,J,K)
      VWN=(1.-FY(I,J,K))*  VW(I,J,K)+  FY(I,J,K)*  VW(I,J+1,K)
      VWS=(1.-FY(I,J-1,K))*VW(I,J-1,K)+FY(I,J-1,K)*VW(I,J,K)
      WWN=(1.-FY(I,J,K))*  WW(I,J,K)+  FY(I,J,K)*  WW(I,J+1,K)
      WWS=(1.-FY(I,J-1,K))*WW(I,J-1,K)+FY(I,J-1,K)*WW(I,J,K)
C
C        NORTH-SOUTH CONTRIBUTION
C
      SUSTS=DENS*(AREXS*UWS+AREYS*VWS+AREZS*WWS)
      SUSTN=DENN*(AREXN*UWN+AREYN*VWN+AREZN*WWN)
C --- RHIE AND CHOW ------------------
C      if(j.eq.2.or.j.eq.njm1.or.k.eq.nkm1)go to 129
C      GO TO 129
      SUSTS=SUSTS+
     1            (1.-FY(I,J-1,K))*(WSTH-WSSTH)*(
     1     ((YX(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS)+
     1      (YY(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS)+
     1      (YZ(I,J-1,K)**2)*(VIS(I,J-1,K)-VISCOS))/VOL(I,J-1,K))
     1          + FY(I,J-1,K)*(WNRT-WSTH)*(
     1     ((YX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTN=SUSTN+
     1            (1.-FY(I,J,K))*(WNRT-WSTH)*(
     1     ((YX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (YZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FY(I,J,K)*(WNNRT-WNRT)*(
     1     ((YX(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS)+
     1      (YY(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS)+
     1      (YZ(I,J+1,K)**2)*(VIS(I,J+1,K)-VISCOS))/VOL(I,J+1,K))
C
c      SUSTS=DENS*(AREXS*UWS+AREYS*VWS+AREZS*WWS)
c      SUSTN=DENN*(AREXN*UWN+AREYN*VWN+AREZN*WWN)
  129 continue
C      FY(I,1,K)=0.
C     FY(I,NJM1,K)=1.
C
C      TOP-BOTTOM FACE AREAS
C
      AREXT=(1.-FZ(I,J,K))*  ZX(I,J,K)+  FZ(I,J,K)*  ZX(I,J,K+1)
      AREXB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)+FZ(I,J,K-1)*ZX(I,J,K)
      AREYT=(1.-FZ(I,J,K))*  ZY(I,J,K)+  FZ(I,J,K)*  ZY(I,J,K+1)
      AREYB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)+FZ(I,J,K-1)*ZY(I,J,K)
      AREZT=(1.-FZ(I,J,K))*  ZZ(I,J,K)+  FZ(I,J,K)*  ZZ(I,J,K+1)
      AREZB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)+FZ(I,J,K-1)*ZZ(I,J,K)
C
C     CELL FACE VELOCITIES
C
      WTOP=(1-FZ(I,J,K))*  W(I,J,K)+  FZ(I,J,K)*  W(I,J,K+1)
      WBOT=(1-FZ(I,J,K-1))*W(I,J,K-1)+FZ(I,J,K-1)*W(I,J,K)
      WTTOP=WTOP+(W(I,J,K+1)-W(I,J,K))
      IF(K.NE.NKM1)
     1 WTTOP=(1-FZ(I,J,K+1))*W(I,J,K+1)+FZ(I,J,K+1)*W(I,J,K+2)
      WBBOT=WBOT-(W(I,J,K)-W(I,J,K-1))
      IF(K.NE.2)
     1 WBBOT=(1-FZ(I,J,K-2))*W(I,J,K-2)+FZ(I,J,K-2)*W(I,J,K-1)
C
C      TOP-BOTTOM FACE STRESSES
C
C      FZ(I,J,1)=0.5
c      FZ(I,J,NKM1)=0.5
C
      UWT=(1.-FZ(I,J,K))*  UW(I,J,K)+  FZ(I,J,K)*  UW(I,J,K+1)
      UWB=(1.-FZ(I,J,K-1))*UW(I,J,K-1)+FZ(I,J,K-1)*UW(I,J,K)
      VWT=(1.-FZ(I,J,K))  *VW(I,J,K)+  FZ(I,J,K)*  VW(I,J,K+1)
      VWB=(1.-FZ(I,J,K-1))*VW(I,J,K-1)+FZ(I,J,K-1)*VW(I,J,K)
      WWT=(1.-FZ(I,J,K))*  WW(I,J,K)+  FZ(I,J,K)*  WW(I,J,K+1)
      WWB=(1.-FZ(I,J,K-1))*WW(I,J,K-1)+FZ(I,J,K-1)*WW(I,J,K)
C
C        TOP-BOTTOM CONTRIBUTION
C
      SUSTB=DENB*(AREXB*UWB+AREYB*VWB+AREZB*WWB)
      SUSTT=DENT*(AREXT*UWT+AREYT*VWT+AREZT*WWT)
C --- RHIE AND CHOW ------------------
C      if(k.eq.2.or.k.eq.nkm1.or.j.eq.2.or.j.eq.njm1)go to 130
C      GO TO 130
      SUSTB=SUSTB+
     1            (1.-FZ(I,J,K-1))*(WBOT-WBBOT)*(
     1     ((ZX(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS)+
     1      (ZY(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS)+
     1      (ZZ(I,J,K-1)**2)*(VIS(I,J,K-1)-VISCOS))/VOL(I,J,K-1))
     1          + FZ(I,J,K-1)*(WTOP-WBOT)*(
     1     ((ZX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
C
      SUSTT=SUSTT+
     1            (1.-FZ(I,J,K))*(WTOP-WBOT)*(
     1     ((ZX(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZY(I,J,K)**2)*(VIS(I,J,K)-VISCOS)+
     1      (ZZ(I,J,K)**2)*(VIS(I,J,K)-VISCOS))/VOL(I,J,K))
     1          + FZ(I,J,K)*(WTTOP-WTOP)*(
     1     ((ZX(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS)+
     1      (ZY(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS)+
     1      (ZZ(I,J,K+1)**2)*(VIS(I,J,K+1)-VISCOS))/VOL(I,J,K+1))
      IF(K.EQ.2)SUSTB=SUSTT
c      IF(K.EQ.NKM1)SUSTT=SUSTB
C
  130 continue
C      FZ(I,J,1)=0.
c      FZ(I,J,NKM1)=1.
C
C   OVERALL STRESS GRADIENT TERM
C
      SU(I,J,K)=SU(I,J,K)+FLOAT(1-IEVM)*
     1          (SUSTW-SUSTE+SUSTS-SUSTN+SUSTB-SUSTT)
C
 10   CONTINUE
C
C -------- PRESSURE GRADIENT SOURCE TERM  ----------------------------
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C      CELL FACE PRESSURES
C
      PRESE=P(I,J,K)*(1.-FX(I,J,K))+P(I+1,J,K)*FX(I,J,K)
      PRESW=P(I-1,J,K)*(1.-FX(I-1,J,K))+P(I,J,K)*FX(I-1,J,K)
      PRESN=P(I,J,K)*(1.-FY(I,J,K))+P(I,J+1,K)*FY(I,J,K)
      PRESS=P(I,J-1,K)*(1.-FY(I,J-1,K))+P(I,J,K)*FY(I,J-1,K)
      PREST=P(I,J,K)*(1.-FZ(I,J,K))+P(I,J,K+1)*FZ(I,J,K)
      PRESB=P(I,J,K-1)*(1.-FZ(I,J,K-1))+P(I,J,K)*FZ(I,J,K-1)
C
C       AREA VECTOR AT CELL CENTER 
C
      DW1(I,J,K)=XZ(I,J,K)
      DW2(I,J,K)=YZ(I,J,K)
      DW3(I,J,K)=ZZ(I,J,K)
C
      SU(I,J,K)=SU(I,J,K)+DW1(I,J,K)*(PRESW-PRESE)
     1 +DW2(I,J,K)*(PRESS-PRESN)+DW3(I,J,K)*(PRESB-PREST)
C
   50 CONTINUE
C
C------ DIFFUSION  SOURCE TERMS ------------------------------------
C
      DO 60 K=2,NKM1
      DO 60 J=2,NJM1
      DO 60 I=2,NIM1
C
C ------ CONTROL VOLUMES FOR DUMMY (STAGGERED) VELOCITY NODES ----
C
       VOLE=VOL(I,J,K)*(1.-FX(I,J,K))+VOL(I+1,J,K)*FX(I,J,K)
       VOLW=VOL(I-1,J,K)*(1.-FX(I-1,J,K))+VOL(I,J,K)*FX(I-1,J,K)
       VOLN=VOL(I,J,K)*(1.-FY(I,J,K))+VOL(I,J+1,K)*FY(I,J,K)
       VOLS=VOL(I,J-1,K)*(1.-FY(I,J-1,K))+VOL(I,J,K)*FY(I,J-1,K)
       VOLT=VOL(I,J,K)*(1.-FZ(I,J,K))+VOL(I,J,K+1)*FZ(I,J,K)
       VOLB=VOL(I,J,K-1)*(1.-FZ(I,J,K-1))+VOL(I,J,K)*FZ(I,J,K-1)
C
C ------ CELL FACE VISCOCITIES ------------------------------------
C
       VISE=VIS(I,J,K)*(1.-FX(I,J,K))+VIS(I+1,J,K)*FX(I,J,K)
       VISW=VIS(I-1,J,K)*(1.-FX(I-1,J,K))+VIS(I,J,K)*FX(I-1,J,K)
       VISN=VIS(I,J,K)*(1.-FY(I,J,K))+VIS(I,J+1,K)*FY(I,J,K)
       VISS=VIS(I,J-1,K)*(1.-FY(I,J-1,K))+VIS(I,J,K)*FY(I,J-1,K)
       VIST=VIS(I,J,K)*(1.-FZ(I,J,K))+VIS(I,J,K+1)*FZ(I,J,K)
       VISB=VIS(I,J,K-1)*(1.-FZ(I,J,K-1))+VIS(I,J,K)*FZ(I,J,K-1)
C
       VISE=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISE
       VISW=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISW
       VISN=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISN
       VISS=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISS
       VIST=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VIST
       VISB=FLOAT(1-IEVM)*VISCOS+FLOAT(IEVM)*VISB
C
C     CALCULATE THE CORNER VALUES FOR U,V,W
C
      UTNE=
     1 U(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+U(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      UTNW=
     1 U(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+U(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+U(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+U(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+U(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      UBNE=
     1 U(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+U(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      UBNW=
     1 U(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+U(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+U(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+U(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      UTSE=
     1 U(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+U(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      UTSW=
     1 U(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+U(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+U(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+U(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+U(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      UBSE=
     1 U(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+U(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+U(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      UBSW=
     1 U(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+U(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+U(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+U(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+U(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+U(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      VTNE=
     1 V(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+V(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      VTNW=
     1 V(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+V(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+V(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+V(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+V(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      VBNE=
     1 V(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+V(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      VBNW=
     1 V(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+V(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+V(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+V(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      VTSE=
     1 V(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+V(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      VTSW=
     1 V(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+V(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+V(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+V(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+V(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      VBSE=
     1 V(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+V(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+V(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      VBSW=
     1 V(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+V(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+V(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+V(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+V(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+V(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      WTNE=
     1 W(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+W(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      WTNW=
     1 W(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+W(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+W(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+W(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+W(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      WBNE=
     1 W(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+W(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      WBNW=
     1 W(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+W(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+W(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+W(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      WTSE=
     1 W(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+W(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      WTSW=
     1 W(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+W(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+W(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+W(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+W(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      WBSE=
     1 W(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+W(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+W(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      WBSW=
     1 W(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+W(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+W(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+W(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+W(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+W(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
C ---- CROSS-DIFFUSION SOURCE TERMS  ------------------------------
C
C         EAST-WEST FACE CROSS-DIFFUSION
C
      Q12E=(XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*YX(I+1,J,K)+XY(I+1,J,K)*YY(I+1,J,K)+
     1 XZ(I+1,J,K)*YZ(I+1,J,K))*FX(I,J,K)
C
      Q12W=(XX(I-1,J,K)*YX(I-1,J,K)+XY(I-1,J,K)*YY(I-1,J,K)+
     1 XZ(I-1,J,K)*YZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*FX(I-1,J,K)
C
      Q13E=(XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*ZX(I+1,J,K)+XY(I+1,J,K)*ZY(I+1,J,K)+
     1 XZ(I+1,J,K)*ZZ(I+1,J,K))*FX(I,J,K)
C
      Q13W=(XX(I-1,J,K)*ZX(I-1,J,K)+XY(I-1,J,K)*ZY(I-1,J,K)+
     1 XZ(I-1,J,K)*ZZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*FX(I-1,J,K)
C
      SCX=0.5*VISE*Q12E/VOLE*(WTNE+WBNE-WTSE-WBSE)
     1   +0.5*VISE*Q13E/VOLE*(WTNE+WTSE-WBNE-WBSE)
     1   -0.5*VISW*Q12W/VOLW*(WTNW+WBNW-WTSW-WBSW)
     1   -0.5*VISW*Q13W/VOLW*(WTNW+WTSW-WBNW-WBSW)
C
C        NORTH-SOUTH FACE CROSS-DIFFUSION
C
      Q21N=(YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*XX(I,J+1,K)+YY(I,J+1,K)*XY(I,J+1,K)+
     1 YZ(I,J+1,K)*XZ(I,J+1,K))*FY(I,J,K)
C
      Q21S=(YX(I,J-1,K)*XX(I,J-1,K)+YY(I,J-1,K)*XY(I,J-1,K)+
     1 YZ(I,J-1,K)*XZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*FY(I,J-1,K)
C
      Q23N=(YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*ZX(I,J+1,K)+YY(I,J+1,K)*ZY(I,J+1,K)+
     1 YZ(I,J+1,K)*ZZ(I,J+1,K))*FY(I,J,K)
C
      Q23S=(YX(I,J-1,K)*ZX(I,J-1,K)+YY(I,J-1,K)*ZY(I,J-1,K)+
     1 YZ(I,J-1,K)*ZZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*FY(I,J-1,K)
C
      SCY=0.5*VISN*Q21N/VOLN*(WTNE+WBNE-WTNW-WBNW)
     1   +0.5*VISN*Q23N/VOLN*(WTNE+WTNW-WBNE-WBNW)
     1   -0.5*VISS*Q21S/VOLS*(WTSE+WBSE-WTSW-WBSW)
     1   -0.5*VISS*Q23S/VOLS*(WTSE+WTSW-WBSE-WBSW)
C
C     TOP-BOTTOM FACE CROSS-DIFUSION
C
      Q31T=(ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*XX(I,J,K+1)+ZY(I,J,K+1)*XY(I,J,K+1)+
     1 ZZ(I,J,K+1)*XZ(I,J,K+1))*FZ(I,J,K)
C
      Q31B=(ZX(I,J,K-1)*XX(I,J,K-1)+ZY(I,J,K-1)*XY(I,J,K-1)+
     1 ZZ(I,J,K-1)*XZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*FZ(I,J,K-1)
C
      Q32T=(ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*YX(I,J,K+1)+ZY(I,J,K+1)*YY(I,J,K+1)+
     1 ZZ(I,J,K+1)*YZ(I,J,K+1))*FZ(I,J,K)
C
      Q32B=(ZX(I,J,K-1)*YX(I,J,K-1)+ZY(I,J,K-1)*YY(I,J,K-1)+
     1 ZZ(I,J,K-1)*YZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*FZ(I,J,K-1)
C
      SCZ=0.5*VIST*Q31T/VOLT*(WTNE+WTSE-WTNW-WTSW)
     1   +0.5*VIST*Q32T/VOLT*(WTNE+WTNW-WTSE-WTSW)
     1   -0.5*VISB*Q31B/VOLB*(WBNE+WBSE-WBNW-WBSW)
     1   -0.5*VISB*Q32B/VOLB*(WBNE+WBNW-WBSE-WBSW)
C
      SU(I,J,K)=SU(I,J,K)+SCX+SCY+SCZ
C
C------ VARIABLE VISCOSITY SOURCE TERMS  --------------------------
C
C        EAST WEST U,V,W GRADIENT TERMS
C
      XXXZE=(1.-FX(I,J,K))*XX(I,J,K)*XZ(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*XZ(I+1,J,K)
      XYXZE=(1.-FX(I,J,K))*XY(I,J,K)*XZ(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*XZ(I+1,J,K)
      XZXZE=(1.-FX(I,J,K))*XZ(I,J,K)*XZ(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*XZ(I+1,J,K)
      XXYZE=(1.-FX(I,J,K))*XX(I,J,K)*YZ(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*YZ(I+1,J,K)
      XYYZE=(1.-FX(I,J,K))*XY(I,J,K)*YZ(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*YZ(I+1,J,K)
      XZYZE=(1.-FX(I,J,K))*XZ(I,J,K)*YZ(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*YZ(I+1,J,K)
      XXZZE=(1.-FX(I,J,K))*XX(I,J,K)*ZZ(I,J,K)
     1     +FX(I,J,K)*XX(I+1,J,K)*ZZ(I+1,J,K)
      XYZZE=(1.-FX(I,J,K))*XY(I,J,K)*ZZ(I,J,K)
     1     +FX(I,J,K)*XY(I+1,J,K)*ZZ(I+1,J,K)
      XZZZE=(1.-FX(I,J,K))*XZ(I,J,K)*ZZ(I,J,K)
     1     +FX(I,J,K)*XZ(I+1,J,K)*ZZ(I+1,J,K)
C
      XXXZW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*XZ(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*XZ(I,J,K)
      XYXZW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*XZ(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*XZ(I,J,K)
      XZXZW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*XZ(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*XZ(I,J,K)
      XXYZW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*YZ(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*YZ(I,J,K)
      XYYZW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*YZ(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*YZ(I,J,K)
      XZYZW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*YZ(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*YZ(I,J,K)
      XXZZW=(1.-FX(I-1,J,K))*XX(I-1,J,K)*ZZ(I-1,J,K)
     1     +FX(I-1,J,K)*XX(I,J,K)*ZZ(I,J,K)
      XYZZW=(1.-FX(I-1,J,K))*XY(I-1,J,K)*ZZ(I-1,J,K)
     1     +FX(I-1,J,K)*XY(I,J,K)*ZZ(I,J,K)
      XZZZW=(1.-FX(I-1,J,K))*XZ(I-1,J,K)*ZZ(I-1,J,K)
     1     +FX(I-1,J,K)*XZ(I,J,K)*ZZ(I,J,K)
C
      SC1=VISE/VOLE*( (U(I+1,J,K)-U(I,J,K))*XXXZE+(V(I+1,J,K)-V(I,J,K))
     1 *XYXZE+(W(I+1,J,K)-W(I,J,K))*XZXZE )
     1   -VISW/VOLW*( (U(I,J,K)-U(I-1,J,K))*XXXZW+(V(I,J,K)-V(I-1,J,K))
     1 *XYXZW+(W(I,J,K)-W(I-1,J,K))*XZXZW )
C
      SC2=VISE/VOLE*( (UTNE+UBNE-UTSE-UBSE)*XXYZE+(VTNE+VBNE-VTSE-VBSE)
     1 *XYYZE+(WTNE+WBNE-WTSE-WBSE)*XZYZE )
     1   -VISW/VOLW*( (UTNW+UBNW-UTSW-UBSW)*XXYZW+(VTNW+VBNW-VTSW-VBSW)
     1 *XYYZW+(WTNW+WBNW-WTSW-WBSW)*XZYZW )
      SC2=0.5*SC2
C
      SC3=VISE/VOLE*( (UTNE+UTSE-UBNE-UBSE)*XXZZE+(VTNE+VTSE-VBNE-VBSE)
     1 *XYZZE+(WTNE+WTSE-WBNE-WBSE)*XZZZE )
     1   -VISW/VOLW*( (UTNW+UTSW-UBNW-UBSW)*XXZZW+(VTNW+VTSW-VBNW-VBSW)
     1 *XYZZW+(WTNW+WTSW-WBNW-WBSW)*XZZZW )
      SC3=0.5*SC3
C
C      NORTH-SOUTH U,V,W GRADIENT TERMS
C
      YXXZN=(1.-FY(I,J,K))*YX(I,J,K)*XZ(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*XZ(I,J+1,K)
      YYXZN=(1.-FY(I,J,K))*YY(I,J,K)*XZ(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*XZ(I,J+1,K)
      YZXZN=(1.-FY(I,J,K))*YZ(I,J,K)*XZ(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*XZ(I,J+1,K)
      YXYZN=(1.-FY(I,J,K))*YX(I,J,K)*YZ(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*YZ(I,J+1,K)
      YYYZN=(1.-FY(I,J,K))*YY(I,J,K)*YZ(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*YZ(I,J+1,K)
      YZYZN=(1.-FY(I,J,K))*YZ(I,J,K)*YZ(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*YZ(I,J+1,K)
      YXZZN=(1.-FY(I,J,K))*YX(I,J,K)*ZZ(I,J,K)
     1     +FY(I,J,K)*YX(I,J+1,K)*ZZ(I,J+1,K)
      YYZZN=(1.-FY(I,J,K))*YY(I,J,K)*ZZ(I,J,K)
     1     +FY(I,J,K)*YY(I,J+1,K)*ZZ(I,J+1,K)
      YZZZN=(1.-FY(I,J,K))*YZ(I,J,K)*ZZ(I,J,K)
     1     +FY(I,J,K)*YZ(I,J+1,K)*ZZ(I,J+1,K)
C
      YXXZS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*XZ(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*XZ(I,J,K)
      YYXZS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*XZ(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*XZ(I,J,K)
      YZXZS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*XZ(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*XZ(I,J,K)
      YXYZS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*YZ(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*YZ(I,J,K)
      YYYZS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*YZ(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*YZ(I,J,K)
      YZYZS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*YZ(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*YZ(I,J,K)
      YXZZS=(1.-FY(I,J-1,K))*YX(I,J-1,K)*ZZ(I,J-1,K)
     1     +FY(I,J-1,K)*YX(I,J,K)*ZZ(I,J,K)
      YYZZS=(1.-FY(I,J-1,K))*YY(I,J-1,K)*ZZ(I,J-1,K)
     1     +FY(I,J-1,K)*YY(I,J,K)*ZZ(I,J,K)
      YZZZS=(1.-FY(I,J-1,K))*YZ(I,J-1,K)*ZZ(I,J-1,K)
     1     +FY(I,J-1,K)*YZ(I,J,K)*ZZ(I,J,K)
C
      SC4=VISN/VOLN*( (UTNE+UBNE-UTNW-UBNW)*YXXZN+(VTNE+VBNE-VTNW-VBNW)
     1 *YYXZN+(WTNE+WBNE-WTNW-WBNW)*YZXZN )
     1   -VISS/VOLS*( (UTSE+UBSE-UTSW-UBSW)*YXXZS+(VTSE+VBSE-VTSW-VBSW)
     1 *YYXZS+(WTSE+WBSE-WTSW-WBSW)*YZXZS )
      SC4=0.5*SC4
C
      SC5=VISN/VOLN*( (U(I,J+1,K)-U(I,J,K))*YXYZN+(V(I,J+1,K)-V(I,J,K))
     1 *YYYZN+(W(I,J+1,K)-W(I,J,K))*YZYZN )
     1   -VISS/VOLS*( (U(I,J,K)-U(I,J-1,K))*YXYZS+(V(I,J,K)-V(I,J-1,K))
     1 *YYYZS+(W(I,J,K)-W(I,J-1,K))*YZYZS )
C
      SC6=VISN/VOLN*( (UTNE+UTNW-UBNE-UBNW)*YXZZN+(VTNE+VTNW-VBNE-VBNW)
     1 *YYZZN+(WTNE+WTNW-WBNE-WBNW)*YZZZN )
     1   -VISS/VOLS*( (UTSE+UTSW-UBSE-UBSW)*YXZZS+(VTSE+VTSW-VBSE-VBSW)
     1 *YYZZS+(WTSE+WTSW-WBSE-WBSW)*YZZZS )
      SC6=0.5*SC6
C
C    TOP-BOTTOM U,V,W GRADIENT TERMS
C
      ZXXZT=(1.-FZ(I,J,K))*ZX(I,J,K)*XZ(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*XZ(I,J,K+1)
      ZYXZT=(1.-FZ(I,J,K))*ZY(I,J,K)*XZ(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*XZ(I,J,K+1)
      ZZXZT=(1.-FZ(I,J,K))*ZZ(I,J,K)*XZ(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*XZ(I,J,K+1)
      ZXYZT=(1.-FZ(I,J,K))*ZX(I,J,K)*YZ(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*YZ(I,J,K+1)
      ZYYZT=(1.-FZ(I,J,K))*ZY(I,J,K)*YZ(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*YZ(I,J,K+1)
      ZZYZT=(1.-FZ(I,J,K))*ZZ(I,J,K)*YZ(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*YZ(I,J,K+1)
      ZXZZT=(1.-FZ(I,J,K))*ZX(I,J,K)*ZZ(I,J,K)
     1     +FZ(I,J,K)*ZX(I,J,K+1)*ZZ(I,J,K+1)
      ZYZZT=(1.-FZ(I,J,K))*ZY(I,J,K)*ZZ(I,J,K)
     1     +FZ(I,J,K)*ZY(I,J,K+1)*ZZ(I,J,K+1)
      ZZZZT=(1.-FZ(I,J,K))*ZZ(I,J,K)*ZZ(I,J,K)
     1     +FZ(I,J,K)*ZZ(I,J,K+1)*ZZ(I,J,K+1)
C
      ZXXZB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*XZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*XZ(I,J,K)
      ZYXZB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*XZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*XZ(I,J,K)
      ZZXZB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*XZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*XZ(I,J,K)
      ZXYZB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*YZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*YZ(I,J,K)
      ZYYZB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*YZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*YZ(I,J,K)
      ZZYZB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*YZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*YZ(I,J,K)
      ZXZZB=(1.-FZ(I,J,K-1))*ZX(I,J,K-1)*ZZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZX(I,J,K)*ZZ(I,J,K)
      ZYZZB=(1.-FZ(I,J,K-1))*ZY(I,J,K-1)*ZZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZY(I,J,K)*ZZ(I,J,K)
      ZZZZB=(1.-FZ(I,J,K-1))*ZZ(I,J,K-1)*ZZ(I,J,K-1)
     1     +FZ(I,J,K-1)*ZZ(I,J,K)*ZZ(I,J,K)
C
      SC7=VIST/VOLT*( (UTNE+UTSE-UTNW-UTSW)*ZXXZT+(VTNE+VTSE-VTNW-VTSW)
     1 *ZYXZT+(WTNE+WTSE-WTNW-WTSW)*ZZXZT )
     1   -VISB/VOLB*( (UBNE+UBSE-UBNW-UBSW)*ZXXZB+(VBNE+VBSE-VBNW-VBSW)
     1 *ZYXZB+(WBNE+WBSE-WBNW-WBSW)*ZZXZB )
      SC7=0.5*SC7
C
      SC8=VIST/VOLT*( (UTNE+UTNW-UTSE-UTSW)*ZXYZT+(VTNE+VTNW-VTSE-VTSW)
     1 *ZYYZT+(WTNE+WTNW-WTSE-WTSW)*ZZYZT )
     1   -VISB/VOLB*( (UBNE+UBNW-UBSE-UBSW)*ZXYZB+(VBNE+VBNW-VBSE-VBSW)
     1 *ZYYZB+(WBNE+WBNW-WBSE-WBSW)*ZZYZB )
      SC8=0.5*SC8
C
      SC9=VIST/VOLT*( (U(I,J,K+1)-U(I,J,K))*ZXZZT+(V(I,J,K+1)-V(I,J,K))
     1 *ZYZZT+(W(I,J,K+1)-W(I,J,K))*ZZZZT )
     1   -VISB/VOLB*( (U(I,J,K)-U(I,J,K-1))*ZXZZB+(V(I,J,K)-V(I,J,K-1))
     1 *ZYZZB+(W(I,J,K)-W(I,J,K-1))*ZZZZB )
C
      SU(I,J,K)=SU(I,J,K)+SC1+SC2+SC3+SC4+SC5+SC6+SC7+SC8+SC9
C
   60 CONTINUE
C---------------------------------------------------------------------
C          ROTATION TERMS                                            *
C---------------------------------------------------------------------
      DO 70 K=2,NKM1
      DO 70 J=2,NJM1
      DO 70 I=2,NIM1
C----  CORIOLIS TERMS -----------------------------------------------
      SU(I,J,K)=SU(I,J,K)+2.*DEN(I,J,K)*VOL(I,J,K)*
     +                    (OMY*U(I,J,K)-OMX*V(I,J,K))
C----  CENTRIFUGAL TERMS --------------------------------------------
      HXX=HX+XC(I,J,K)
      HYY=HY+YC(I,J,K)
      HZZ=HZ+ZC(I,J,K)
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*
     +                    (HZZ*(OMX*OMX+OMY*OMY)-OMZ*(HXX*OMX+HYY*OMY))
   70 CONTINUE
C------ BOUNDARY CONDITIONS   ---------------------------------------
C
      CALL MODW
C
      RESORW=0.
      DO 30 K=2,NKM1
      DO 30 J=2,NJM1
      DO 30 I=2,NIM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*W(I+1,J,K)+AW(I,J,K)*W(I-1,J,K)
     1     +AN(I,J,K)*W(I,J+1,K)+AS(I,J,K)*W(I,J-1,K)
     1     +AT(I,J,K)*W(I,J,K+1)+AB(I,J,K)*W(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*W(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
      RESORW=RESORW+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFW
      SU(I,J,K)=SU(I,J,K)+(1.-URFW)*AP(I,J,K)*W(I,J,K)
      DW1(I,J,K)=DW1(I,J,K)/AP(I,J,K)
      DW2(I,J,K)=DW2(I,J,K)/AP(I,J,K)
      DW3(I,J,K)=DW3(I,J,K)/AP(I,J,K)
 30   CONTINUE
   36 CONTINUE
C
      DO 40 N=1,NSWPW
      CALL LISOLV(4,2,2,NI,NJ,NK,W)
C      CALL LISOLV(4,NJR+1,2,NI,NJRN+1,NK,W)
C      CALL LISOLV(NIR+1,2,2,NI-2,NJ,NK,W)
 40   CONTINUE
      RETURN
      END
C
      SUBROUTINE CALCKE
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      NJRSL=NJR+INWP
      NJRSL1=NJRSL+1
      NJRNL=NJRN-INWP+1
      NJRNL1=NJRNL-1
      NIEL=NIR+INWP
      NIEL1=NIEL+1
      NMWL=NIRMS-INWP+1
      NMWL1=NMWL-1
      NMEL=NIRMF+INWP
      NMEL=NI-1
      NMEL1=NMEL+1
      NIWL=NI-INWP-2
      NIWL1=NIWL-1
C
      PR=PRTKE
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*TKE(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(1.,PRTKE)
C
      CALL BQUICK(1,1.,TKE)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTKE,1,TKE)
C
C   NEAR-WALL REGIONS (TOP,BOTTOM,NORTH AND SOUTH) WHEN  THE 1-EQN
C   NEAR-WALL MODEL IS EMPLOYED
C
C-------------- SOUTH WALL ----------------------
      DO 60 K=2,NKM1
      DO 60 J=2,JED-1
      DO 60 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*PK1(I,J,K)*VOL(I,J,K)
C      SP(I,J,K)=-DEN(I,J,K)*SQRT(TKE(I,J,K))*VOL(I,J,K)/(ELM+TINY)
   60 CONTINUE
C------------- NORTH WALL ----------------------
      DO 70 K=2,NKM1
      DO 70 J=NJED+1,NJM1
      DO 70 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*PK1(I,J,K)*VOL(I,J,K)
C      SP(I,J,K)=-DEN(I,J,K)*SQRT(TKE(I,J,K))*VOL(I,J,K)/(ELM+TINY)
   70 CONTINUE
C------------- TOP WALL -----------------------
      DO 80 K=NKED+1,NKM1
      DO 80 J=2,NJM1
      DO 80 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*PK1(I,J,K)*VOL(I,J,K)
C      SP(I,J,K)=-DEN(I,J,K)*SQRT(TKE(I,J,K))*VOL(I,J,K)/(ELM+TINY)
   80 CONTINUE
C------------ BOTTOM WALL ---------------------
      DO 90 K=2,KED-1
      DO 90 J=2,NJM1
      DO 90 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*PK1(I,J,K)*VOL(I,J,K)
C      SP(I,J,K)=-DEN(I,J,K)*SQRT(TKE(I,J,K))*VOL(I,J,K)/(ELM+TINY)
   90 CONTINUE
C----------- RIBS  ---------------------------------------------
C      IF(INWM.EQ.1)THEN
C      DO 21 K=KED,NKED
C
C      DO 22 J=JED,NJRSL-1
C      DO 23 I=2,NIEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   23 CONTINUE
C      DO 24 I=NMWL+1,NMEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   24 CONTINUE
C      DO 25 I=NIWL+1,NI
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   25 CONTINUE
C   22 CONTINUE
C
C      DO 26 J=NJRNL+1,NJED
C      DO 27 I=2,NIEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   27 CONTINUE
C      DO 28 I=NMWL+1,NMEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   28 CONTINUE
C      DO 29 I=NIWL+1,NI
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   29 CONTINUE
C   26 CONTINUE
C
C   21 CONTINUE
C      ENDIF
C
C ----- DISIPPATION & GENERATION RATE -------------------------
C
C         MAIN FLOW REGION
C
C      DO 50 K=KED,NKED
C      DO 50 J=JED,NJED
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C-------  HIGH-RE DISSIPATION AND GENERATION ------------------
      SU1=DEN(I,J,K)*(PK1(I,J,K)-TED(I,J,K))*VOL(I,J,K)
C----  LOW-RE NEAR-WALL DISSIPATION RATE ---------------------
      SP(I,J,K)=SP(I,J,K)+FLOAT(1-INWM)*DEN(I,J,K)*VOL(I,J,K)*
     1           EDNW(I,J,K)/(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+MAX(SU1,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SU1,0.)/(TKE(I,J,K)+TINY)
 50   CONTINUE
C
      CALL MODTKE
C
      RESORK=0.
      DO 30 I=2,NIM1
      DO 30 J=2,NJM1
      DO 30 K=2,NKM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*TKE(I+1,J,K)+AW(I,J,K)*TKE(I-1,J,K)
     1     +AN(I,J,K)*TKE(I,J+1,K)+AS(I,J,K)*TKE(I,J-1,K)
     1     +AT(I,J,K)*TKE(I,J,K+1)+AB(I,J,K)*TKE(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*TKE(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
      RESORK=RESORK+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFK
      SU(I,J,K)=SU(I,J,K)+(1.-URFK)*AP(I,J,K)*TKE(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPK
      CALL LISOLV(3,2,2,NIM1,NJ,NK,TKE)
 40   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CALCUU
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      NJRSL=NJR+INWP
      NJRSL1=NJRSL+1
      NJRNL=NJRN-INWP+1
      NJRNL1=NJRNL-1
      NIEL=NIR+INWP
      NIEL1=NIEL+1
      NMWL=NIRMS-INWP+1
      NMWL1=NMWL-1
      NMEL=NIRMF+INWP
      NMEL=NI-1
      NMEL1=NMEL+1
      NIWL=NI-INWP-2
      NIWL1=NIWL-1
C
      PR=PRTKE
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*UU(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(1.,PRTKE)
C
      CALL BQUICK(1,1.,UU)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTKE,1,UU)
C
C   NEAR-WALL REGIONS (TOP,BOTTOM,NORTH AND SOUTH) WHEN  THE 1-EQN
C   NEAR-WALL MODEL IS EMPLOYED
C
C-------------- SOUTH WALL ----------------------
      DO 60 K=2,NKM1
      DO 60 J=2,JED-1
      DO 60 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
   60 CONTINUE
C------------- NORTH WALL ----------------------
      DO 70 K=2,NKM1
      DO 70 J=NJED+1,NJM1
      DO 70 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
   70 CONTINUE
C------------- TOP WALL -----------------------
      DO 80 K=NKED+1,NKM1
      DO 80 J=2,NJM1
      DO 80 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
   80 CONTINUE
C------------ BOTTOM WALL ---------------------
      DO 90 K=2,KED-1
      DO 90 J=2,NJM1
      DO 90 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS+TINY
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
   90 CONTINUE
C----------- RIBS  ---------------------------------------------
C      IF(INWM.EQ.1)THEN
C      DO 21 K=KED,NKED
C
C      DO 22 J=JED,NJRSL-1
C      DO 23 I=2,NIEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   23 CONTINUE
C      DO 24 I=NMWL+1,NMEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   24 CONTINUE
C      DO 25 I=NIWL+1,NI
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   25 CONTINUE
C   22 CONTINUE
C
C      DO 26 J=NJRNL+1,NJED
C      DO 27 I=2,NIEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   27 CONTINUE
C      DO 28 I=NMWL+1,NMEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   28 CONTINUE
C      DO 29 I=NIWL+1,NI
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-.263*XLST))
C      TED(I,J,K)=(TKE(I,J,K)**1.5)/(ELM+TINY)
C   29 CONTINUE
C   26 CONTINUE
C
C   21 CONTINUE
C      ENDIF
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C -------- GENERATION AND DISSIPATOIN SOURCE TERMS --------------------
C
      RT=DEN(I,J,K)*TKE(I,J,K)*TKE(I,J,K)/(VISCOS*(TED(I,J,K)+TINY))
      YSTR=DEN(I,J,K)*XL3(I,J,K)*(SQRT(TKE(I,J,K)))/VISCOS
      VISTUR=VIS(I,J,K)-VISCOS
      FE=FLOAT(INWM-1)*EXP(-RT/8.)+FLOAT(2-INWM)*EXP(-YSTR/3.)
C      FE=0.
      TEDT=TED(I,J,K)+FLOAT(INWM-1)*EDNW(I,J,K)
      GNUU=-2.*(UU(I,J,K)*DUDXC(I,J,K)+UV(I,J,K)*DUDYC(I,J,K)+
     1          UW(I,J,K)*DUDZC(I,J,K))
     1      +4.*OMZ*UV(I,J,K)-4.*OMY*UW(I,J,K)
C      GNUU=VISTUR*((DUDXC(I,J,K)**2)+(DUDYC(I,J,K)**2)+
C     1          (DUDZC(I,J,K)**2))
      SUGE=VOL(I,J,K)*DEN(I,J,K)*(GNUU-.667*(1.-FE)*TEDT)
      SU(I,J,K)=SU(I,J,K)+MAX(SUGE,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SUGE,0.)/(UU(I,J,K)+TINY)
      SP(I,J,K)=SP(I,J,K)-FE*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                    /(TKE(I,J,K)+TINY)
C
C---------   REDISTRIBUTION TERMS  (FIJ1 & FIJ2)  ---------------------
C            FIJ1
      SP(I,J,K)=SP(I,J,K)-CA1*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                     /(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+CA1*DEN(I,J,K)*VOL(I,J,K)*.667*TEDT
C            FIJ2
      FUU2=-CA2*(GNUU-.667*PK1(I,J,K))
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*FUU2
C
C --------- LOW-RE SOURCE TERMS  --------------------------------------
C
      FH=FLOAT(INWM-1)*(10.+2.6*RT)*EXP(-RT/20.)
     1  +FLOAT(2-INWM)*(10.2+7.5*YSTR)*EXP(-.05*YSTR)
      FJ=.06*(FLOAT(INWM-1)*EXP(-RT/8)+FLOAT(2-INWM)*EXP(-YSTR/3.))
C      FH=0.
C      FJ=0.
C
C          -(HIJ-2*HJJ/3)    TERM
C
      SP(I,J,K)=SP(I,J,K)-4.*FH*VISCOS*VOL(I,J,K)*DKDX(I,J,K)*
     1                    DKDX(I,J,K)/(3.*TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+.667*FH*VISCOS*VOL(I,J,K)*(
     1          DKDY(I,J,K)*DKDY(I,J,K)*VV(I,J,K)
     1         +DKDZ(I,J,K)*DKDZ(I,J,K)*WW(I,J,K)
     1         -DKDX(I,J,K)*DKDY(I,J,K)*UV(I,J,K)
     1         -DKDX(I,J,K)*DKDZ(I,J,K)*UW(I,J,K)
     1      +2.*DKDY(I,J,K)*DKDZ(I,J,K)*VW(I,J,K))/(TKE(I,J,K)+TINY)
C
C            J  TERM
C 
       SU(I,J,K)=SU(I,J,K)+FJ*DEN(I,J,K)*VOL(I,J,K)*TKE(I,J,K)*
     1                     2.*DUDXC(I,J,K)
C
C------ WALL REFLECTION TERMS  ----------------------------------------
C
      FW=FLOAT(INWM-1)*(1.-EXP(-RT/20.))*(1.+EXP(-RT/100.))
     +  +FLOAT(2-INWM)*(1.-EXP(-.12*YSTR))*(1.+EXP(-0.03*YSTR))
      FW=FW*VOL(I,J,K)*DEN(I,J,K)
C
      GNVV=-2.*(UV(I,J,K)*DVDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
     1          VW(I,J,K)*DVDZC(I,J,K))
     1      +4.*OMX*VW(I,J,K)-4*OMZ*UV(I,J,K)
      GNWW=-2.*(UW(I,J,K)*DWDXC(I,J,K)+VW(I,J,K)*DWDYC(I,J,K)+
     1          WW(I,J,K)*DWDZC(I,J,K))
     1      +4.*OMY*UW(I,J,K)-4.*OMX*VW(I,J,K)
      GNUV=-(VV(I,J,K)*DUDYC(I,J,K)+UU(I,J,K)*DVDXC(I,J,K)+
     1       UV(I,J,K)*(DUDXC(I,J,K)+DVDYC(I,J,K))+
     1       UW(I,J,K)*DVDZC(I,J,K)+VW(I,J,K)*DUDZC(I,J,K))
     1       +2.*OMZ*(VV(I,J,K)-UU(I,J,K))+2.*OMX*UW(I,J,K)
     1       -2.*OMY*VW(I,J,K)
      GNUW=-(WW(I,J,K)*DUDZC(I,J,K)+UU(I,J,K)*DWDXC(I,J,K)+
     1       UW(I,J,K)*(DUDXC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDYC(I,J,K)+VW(I,J,K)*DUDYC(I,J,K))
     1       +2.*OMY*(UU(I,J,K)-WW(I,J,K))+2.*OMZ*VW(I,J,K)
     1       -2.*OMX*UV(I,J,K) 
      GNVW=-(WW(I,J,K)*DVDZC(I,J,K)+VV(I,J,K)*DWDYC(I,J,K)+
     1       VW(I,J,K)*(DVDYC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDXC(I,J,K)+UW(I,J,K)*DVDXC(I,J,K))
     1       +2.*OMX*(WW(I,J,K)-UU(I,J,K))+2.*OMY*UV(I,J,K)
     1       -2.*OMZ*UW(I,J,K)
C
      FVV2=-CA2*(GNVV-.667*PK1(I,J,K))
      FWW2=-CA2*(GNWW-.667*PK1(I,J,K))
      FUV2=-CA2*GNUV
      FUW2=-CA2*GNUW
      FVW2=-CA2*GNVW
C
C---------- EAST WALL --------------------------------------------------
C
      RLE=ELE(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLE=MIN(1.,RLE)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLE*CA1W*TED(I,J,K)*2.*(COSXE(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLE*CA1W*TED(I,J,K)*(
     1    (COSYE(I,J,K)**2)*VV(I,J,K)+(COSZE(I,J,K)**2)*WW(I,J,K)
     2       -COSXE(I,J,K)*COSYE(I,J,K)*UV(I,J,K)
     3       -COSXE(I,J,K)*COSZE(I,J,K)*UW(I,J,K)
     4    +2.*COSYE(I,J,K)*COSZE(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLE*CA2W*(-2.*(COSXE(I,J,K)**2)*FUU2    
     1    +(COSYE(I,J,K)**2)*FVV2+(COSZE(I,J,K)**2)*FWW2
     2         -COSXE(I,J,K)*COSYE(I,J,K)*FUV2
     3         -COSXE(I,J,K)*COSZE(I,J,K)*FUW2
     4      +2.*COSYE(I,J,K)*COSZE(I,J,K)*FVW2           )
C
C---------- WEST WALL --------------------------------------------------
C
      RLW=ELW(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLW=MIN(1.,RLW)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLW*CA1W*TED(I,J,K)*2.*(COSXW(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLW*CA1W*TED(I,J,K)*(
     1    (COSYW(I,J,K)**2)*VV(I,J,K)+(COSZW(I,J,K)**2)*WW(I,J,K)
     2       -COSXW(I,J,K)*COSYW(I,J,K)*UV(I,J,K)
     3       -COSXW(I,J,K)*COSZW(I,J,K)*UW(I,J,K)
     4    +2.*COSYW(I,J,K)*COSZW(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLW*CA2W*(-2.*(COSXW(I,J,K)**2)*FUU2    
     1    +(COSYW(I,J,K)**2)*FVV2+(COSZW(I,J,K)**2)*FWW2
     2         -COSXW(I,J,K)*COSYW(I,J,K)*FUV2
     3         -COSXW(I,J,K)*COSZW(I,J,K)*FUW2
     4      +2.*COSYW(I,J,K)*COSZW(I,J,K)*FVW2           )
C
C---------- NORTH WALL -------------------------------------------------
C
      RLN=ELN(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLN=MIN(1.,RLN)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLN*CA1W*TED(I,J,K)*2.*(COSXN(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLN*CA1W*TED(I,J,K)*(
     1    (COSYN(I,J,K)**2)*VV(I,J,K)+(COSZN(I,J,K)**2)*WW(I,J,K)
     2       -COSXN(I,J,K)*COSYN(I,J,K)*UV(I,J,K)
     3       -COSXN(I,J,K)*COSZN(I,J,K)*UW(I,J,K)
     4    +2.*COSYN(I,J,K)*COSZN(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLN*CA2W*(-2.*(COSXN(I,J,K)**2)*FUU2    
     1    +(COSYN(I,J,K)**2)*FVV2+(COSZN(I,J,K)**2)*FWW2
     2         -COSXN(I,J,K)*COSYN(I,J,K)*FUV2
     3         -COSXN(I,J,K)*COSZN(I,J,K)*FUW2
     4      +2.*COSYN(I,J,K)*COSZN(I,J,K)*FVW2           )
C
C---------- SOUTH WALL -------------------------------------------------
C
      RLS=ELS(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLS=MIN(1.,RLS)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLS*CA1W*TED(I,J,K)*2.*(COSXS(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLS*CA1W*TED(I,J,K)*(
     1    (COSYS(I,J,K)**2)*VV(I,J,K)+(COSZS(I,J,K)**2)*WW(I,J,K)
     2       -COSXS(I,J,K)*COSYS(I,J,K)*UV(I,J,K)
     3       -COSXS(I,J,K)*COSZS(I,J,K)*UW(I,J,K)
     4    +2.*COSYS(I,J,K)*COSZS(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLS*CA2W*(-2.*(COSXS(I,J,K)**2)*FUU2    
     1    +(COSYS(I,J,K)**2)*FVV2+(COSZS(I,J,K)**2)*FWW2
     2         -COSXS(I,J,K)*COSYS(I,J,K)*FUV2
     3         -COSXS(I,J,K)*COSZS(I,J,K)*FUW2
     4      +2.*COSYS(I,J,K)*COSZS(I,J,K)*FVW2           )
C
C---------- TOP WALL -------------------------------------------------
C
      RLT=ELT(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLT=MIN(1.,RLT)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLT*CA1W*TED(I,J,K)*2.*(COSXT(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLT*CA1W*TED(I,J,K)*(
     1    (COSYT(I,J,K)**2)*VV(I,J,K)+(COSZT(I,J,K)**2)*WW(I,J,K)
     2       -COSXT(I,J,K)*COSYT(I,J,K)*UV(I,J,K)
     3       -COSXT(I,J,K)*COSZT(I,J,K)*UW(I,J,K)
     4    +2.*COSYT(I,J,K)*COSZT(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLT*CA2W*(-2.*(COSXT(I,J,K)**2)*FUU2    
     1    +(COSYT(I,J,K)**2)*FVV2+(COSZT(I,J,K)**2)*FWW2
     2         -COSXT(I,J,K)*COSYT(I,J,K)*FUV2
     3         -COSXT(I,J,K)*COSZT(I,J,K)*FUW2
     4      +2.*COSYT(I,J,K)*COSZT(I,J,K)*FVW2           )
C
C---------- BOTTOM WALL ----------------------------------------------
C
      RLB=ELB(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLB=MIN(1.,RLB)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLB*CA1W*TED(I,J,K)*2.*(COSXB(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLB*CA1W*TED(I,J,K)*(
     1    (COSYB(I,J,K)**2)*VV(I,J,K)+(COSZB(I,J,K)**2)*WW(I,J,K)
     2       -COSXB(I,J,K)*COSYB(I,J,K)*UV(I,J,K)
     3       -COSXB(I,J,K)*COSZB(I,J,K)*UW(I,J,K)
     4    +2.*COSYB(I,J,K)*COSZB(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLB*CA2W*(-2.*(COSXB(I,J,K)**2)*FUU2    
     1    +(COSYB(I,J,K)**2)*FVV2+(COSZB(I,J,K)**2)*FWW2
     2         -COSXB(I,J,K)*COSYB(I,J,K)*FUV2
     3         -COSXB(I,J,K)*COSZB(I,J,K)*FUW2
     4      +2.*COSYB(I,J,K)*COSZB(I,J,K)*FVW2           )
C
C ------ OVERAL SOURCE TREATMENT FOR POSITIVE ANSWERS ------------------
C
      SP(I,J,K)=SP(I,J,K)-MAX(-SU(I,J,K),0.)/(UU(I,J,K)+TINY)
      SU(I,J,K)=MAX(SU(I,J,K),0.)
 50   CONTINUE
C
      CALL MODTUU
C
      RESORK=0.
      DO 30 K=2,NKM1
      DO 30 J=2,NJM1
      DO 30 I=2,NIM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*UU(I+1,J,K)+AW(I,J,K)*UU(I-1,J,K)
     1     +AN(I,J,K)*UU(I,J+1,K)+AS(I,J,K)*UU(I,J-1,K)
     1     +AT(I,J,K)*UU(I,J,K+1)+AB(I,J,K)*UU(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*UU(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
      RESORK=RESORK+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFK
      SU(I,J,K)=SU(I,J,K)+(1.-URFK)*AP(I,J,K)*UU(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPK
      CALL LISOLV(3,2,2,NIM1,NJ,NK,UU)
 40   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CALCVV
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      PR=PRTKE
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*VV(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(1.,PRTKE)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTKE,1,VV)
C
      CALL BQUICK(1,1.,VV)
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C -------- GENERATION AND DISSIPATOIN SOURCE TERMS --------------------
C
      RT=DEN(I,J,K)*TKE(I,J,K)*TKE(I,J,K)/(VISCOS*(TED(I,J,K)+TINY))
      YSTR=DEN(I,J,K)*XL3(I,J,K)*(SQRT(TKE(I,J,K)))/VISCOS
      VISTUR=VIS(I,J,K)-VISCOS
      FE=FLOAT(INWM-1)*EXP(-RT/8.)+FLOAT(2-INWM)*EXP(-YSTR/3.)
C      FE=0.
      TEDT=TED(I,J,K)+FLOAT(INWM-1)*EDNW(I,J,K)
      GNVV=-2.*(UV(I,J,K)*DVDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
     1          VW(I,J,K)*DVDZC(I,J,K))
     1      +4.*OMX*VW(I,J,K)-4*OMZ*UV(I,J,K)
C      GNVV=VISTUR*((DVDXC(I,J,K)**2)+(DVDYC(I,J,K)**2)+
C     1          (DVDZC(I,J,K)**2))
      SUGE=VOL(I,J,K)*DEN(I,J,K)*(GNVV-.667*(1.-FE)*TEDT)
      SU(I,J,K)=SU(I,J,K)+MAX(SUGE,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SUGE,0.)/(VV(I,J,K)+TINY)
      SP(I,J,K)=SP(I,J,K)-FE*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                    /(TKE(I,J,K)+TINY)
C
C---------   REDISTRIBUTION TERMS  (FIJ1 & FIJ2)  ---------------------
C            FIJ1
      SP(I,J,K)=SP(I,J,K)-CA1*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                     /(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+CA1*DEN(I,J,K)*VOL(I,J,K)*.667*TEDT
C            FIJ2
      FVV2=-CA2*(GNVV-.667*PK1(I,J,K))
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*FVV2
C
C --------- LOW-RE SOURCE TERMS  --------------------------------------
C
      FH=FLOAT(INWM-1)*(10.+2.6*RT)*EXP(-RT/20.)
     1  +FLOAT(2-INWM)*(10.2+7.5*YSTR)*EXP(-.05*YSTR)
      FJ=.06*(FLOAT(INWM-1)*EXP(-RT/8)+FLOAT(2-INWM)*EXP(-YSTR/3.))
C      FH=0.
C      FJ=0.
C
C          -(HIJ-2*HJJ/3)    TERM
C
      SP(I,J,K)=SP(I,J,K)-4.*FH*VISCOS*VOL(I,J,K)*DKDY(I,J,K)*
     1                    DKDY(I,J,K)/(3.*TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+.667*FH*VISCOS*VOL(I,J,K)*(
     1          DKDX(I,J,K)*DKDX(I,J,K)*UU(I,J,K)
     1         +DKDZ(I,J,K)*DKDZ(I,J,K)*WW(I,J,K)
     1         -DKDX(I,J,K)*DKDY(I,J,K)*UV(I,J,K)
     1      +2.*DKDX(I,J,K)*DKDZ(I,J,K)*UW(I,J,K)
     1         -DKDY(I,J,K)*DKDZ(I,J,K)*VW(I,J,K))/(TKE(I,J,K)+TINY)
C
C            J  TERM
C 
       SU(I,J,K)=SU(I,J,K)+FJ*DEN(I,J,K)*VOL(I,J,K)*TKE(I,J,K)*
     1                     2.*DVDYC(I,J,K)
C
C------ WALL REFLECTION TERMS  ----------------------------------------
C
      FW=FLOAT(INWM-1)*(1.-EXP(-RT/20.))*(1.+EXP(-RT/100.))
     +  +FLOAT(2-INWM)*(1.-EXP(-.12*YSTR))*(1.+EXP(-0.03*YSTR))
      FW=FW*VOL(I,J,K)*DEN(I,J,K)
C
      GNUU=-2.*(UU(I,J,K)*DUDXC(I,J,K)+UV(I,J,K)*DUDYC(I,J,K)+
     1          UW(I,J,K)*DUDZC(I,J,K))
     1      +4.*OMZ*UV(I,J,K)-4.*OMY*UW(I,J,K)
C      GNVV=-2.*(UV(I,J,K)*DVDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
C     1          VW(I,J,K)*DVDZC(I,J,K))
      GNWW=-2.*(UW(I,J,K)*DWDXC(I,J,K)+VW(I,J,K)*DWDYC(I,J,K)+
     1          WW(I,J,K)*DWDZC(I,J,K))
     1      +4.*OMY*UW(I,J,K)-4.*OMX*VW(I,J,K)
      GNUV=-(VV(I,J,K)*DUDYC(I,J,K)+UU(I,J,K)*DVDXC(I,J,K)+
     1       UV(I,J,K)*(DUDXC(I,J,K)+DVDYC(I,J,K))+
     1       UW(I,J,K)*DVDZC(I,J,K)+VW(I,J,K)*DUDZC(I,J,K))
     1       +2.*OMZ*(VV(I,J,K)-UU(I,J,K))+2.*OMX*UW(I,J,K)
     1       -2.*OMY*VW(I,J,K)
      GNUW=-(WW(I,J,K)*DUDZC(I,J,K)+UU(I,J,K)*DWDXC(I,J,K)+
     1       UW(I,J,K)*(DUDXC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDYC(I,J,K)+VW(I,J,K)*DUDYC(I,J,K))
     1       +2.*OMY*(UU(I,J,K)-WW(I,J,K))+2.*OMZ*VW(I,J,K)
     1       -2.*OMX*UV(I,J,K) 
      GNVW=-(WW(I,J,K)*DVDZC(I,J,K)+VV(I,J,K)*DWDYC(I,J,K)+
     1       VW(I,J,K)*(DVDYC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDXC(I,J,K)+UW(I,J,K)*DVDXC(I,J,K))
     1       +2.*OMX*(WW(I,J,K)-UU(I,J,K))+2.*OMY*UV(I,J,K)
     1       -2.*OMZ*UW(I,J,K)
C
      FUU2=-CA2*(GNUU-.667*PK1(I,J,K))
C      FVV2=-CA2*(GNVV-.667*PK1(I,J,K))
      FWW2=-CA2*(GNWW-.667*PK1(I,J,K))
      FUV2=-CA2*GNUV
      FUW2=-CA2*GNUW
      FVW2=-CA2*GNVW
C
C---------- EAST WALL --------------------------------------------------
C
      RLE=ELE(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLE=MIN(1.,RLE)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLE*CA1W*TED(I,J,K)*2.*(COSYE(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLE*CA1W*TED(I,J,K)*(
     1    (COSXE(I,J,K)**2)*UU(I,J,K)+(COSZE(I,J,K)**2)*WW(I,J,K)
     2       -COSXE(I,J,K)*COSYE(I,J,K)*UV(I,J,K)
     3    +2.*COSXE(I,J,K)*COSZE(I,J,K)*UW(I,J,K)
     4       -COSYE(I,J,K)*COSZE(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLE*CA2W*(-2.*(COSYE(I,J,K)**2)*FVV2    
     1    +(COSXE(I,J,K)**2)*FUU2+(COSZE(I,J,K)**2)*FWW2
     2         -COSXE(I,J,K)*COSYE(I,J,K)*FUV2
     3      +2.*COSXE(I,J,K)*COSZE(I,J,K)*FUW2
     4         -COSYE(I,J,K)*COSZE(I,J,K)*FVW2           )
C
C---------- WEST WALL --------------------------------------------------
C
      RLW=ELW(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLW=MIN(1.,RLW)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLW*CA1W*TED(I,J,K)*2.*(COSYW(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLW*CA1W*TED(I,J,K)*(
     1    (COSXW(I,J,K)**2)*UU(I,J,K)+(COSZW(I,J,K)**2)*WW(I,J,K)
     2       -COSXW(I,J,K)*COSYW(I,J,K)*UV(I,J,K)
     3    +2.*COSXW(I,J,K)*COSZW(I,J,K)*UW(I,J,K)
     4       -COSYW(I,J,K)*COSZW(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLW*CA2W*(-2.*(COSYW(I,J,K)**2)*FVV2    
     1    +(COSXW(I,J,K)**2)*FUU2+(COSZW(I,J,K)**2)*FWW2
     2         -COSXW(I,J,K)*COSYW(I,J,K)*FUV2
     3      +2.*COSXW(I,J,K)*COSZW(I,J,K)*FUW2
     4         -COSYW(I,J,K)*COSZW(I,J,K)*FVW2           )
C
C---------- NORTH WALL -------------------------------------------------
C
      RLN=ELN(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLN=MIN(1.,RLN)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLN*CA1W*TED(I,J,K)*2.*(COSYN(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLN*CA1W*TED(I,J,K)*(
     1    (COSXN(I,J,K)**2)*UU(I,J,K)+(COSZN(I,J,K)**2)*WW(I,J,K)
     2       -COSXN(I,J,K)*COSYN(I,J,K)*UV(I,J,K)
     3    +2.*COSXN(I,J,K)*COSZN(I,J,K)*UW(I,J,K)
     4       -COSYN(I,J,K)*COSZN(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLN*CA2W*(-2.*(COSYN(I,J,K)**2)*FVV2    
     1    +(COSXN(I,J,K)**2)*FUU2+(COSZN(I,J,K)**2)*FWW2
     2         -COSXN(I,J,K)*COSYN(I,J,K)*FUV2
     3      +2.*COSXN(I,J,K)*COSZN(I,J,K)*FUW2
     4         -COSYN(I,J,K)*COSZN(I,J,K)*FVW2           )
C
C---------- SOUTH WALL -------------------------------------------------
C
      RLS=ELS(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLS=MIN(1.,RLS)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLS*CA1W*TED(I,J,K)*2.*(COSYS(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLS*CA1W*TED(I,J,K)*(
     1    (COSXS(I,J,K)**2)*UU(I,J,K)+(COSZS(I,J,K)**2)*WW(I,J,K)
     2       -COSXS(I,J,K)*COSYS(I,J,K)*UV(I,J,K)
     3    +2.*COSXS(I,J,K)*COSZS(I,J,K)*UW(I,J,K)
     4       -COSYS(I,J,K)*COSZS(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLS*CA2W*(-2.*(COSYS(I,J,K)**2)*FVV2    
     1    +(COSXS(I,J,K)**2)*FUU2+(COSZS(I,J,K)**2)*FWW2
     2         -COSXS(I,J,K)*COSYS(I,J,K)*FUV2
     3      +2.*COSXS(I,J,K)*COSZS(I,J,K)*FUW2
     4         -COSYS(I,J,K)*COSZS(I,J,K)*FVW2           )
C
C---------- TOP WALL -------------------------------------------------
C
      RLT=ELT(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLT=MIN(1.,RLT)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLT*CA1W*TED(I,J,K)*2.*(COSYT(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLT*CA1W*TED(I,J,K)*(
     1    (COSXT(I,J,K)**2)*UU(I,J,K)+(COSZT(I,J,K)**2)*WW(I,J,K)
     2       -COSXT(I,J,K)*COSYT(I,J,K)*UV(I,J,K)
     3    +2.*COSXT(I,J,K)*COSZT(I,J,K)*UW(I,J,K)
     4       -COSYT(I,J,K)*COSZT(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLT*CA2W*(-2.*(COSYT(I,J,K)**2)*FVV2    
     1    +(COSXT(I,J,K)**2)*FUU2+(COSZT(I,J,K)**2)*FWW2
     2         -COSXT(I,J,K)*COSYT(I,J,K)*FUV2
     3      +2.*COSXT(I,J,K)*COSZT(I,J,K)*FUW2
     4         -COSYT(I,J,K)*COSZT(I,J,K)*FVW2           )
C
C---------- BOTTOM WALL ----------------------------------------------
C
      RLB=ELB(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLB=MIN(1.,RLB)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLB*CA1W*TED(I,J,K)*2.*(COSYB(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLB*CA1W*TED(I,J,K)*(
     1    (COSXB(I,J,K)**2)*UU(I,J,K)+(COSZB(I,J,K)**2)*WW(I,J,K)
     2       -COSXB(I,J,K)*COSYB(I,J,K)*UV(I,J,K)
     3    +2.*COSXB(I,J,K)*COSZB(I,J,K)*UW(I,J,K)
     4       -COSYB(I,J,K)*COSZB(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLB*CA2W*(-2.*(COSYB(I,J,K)**2)*FVV2    
     1    +(COSXB(I,J,K)**2)*FUU2+(COSZB(I,J,K)**2)*FWW2
     2         -COSXB(I,J,K)*COSYB(I,J,K)*FUV2
     3      +2.*COSXB(I,J,K)*COSZB(I,J,K)*FUW2
     4         -COSYB(I,J,K)*COSZB(I,J,K)*FVW2           )
C
C ------ OVERAL SOURCE TREATMENT FOR POSITIVE ANSWERS ------------------
C
      SP(I,J,K)=SP(I,J,K)-MAX(-SU(I,J,K),0.)/(VV(I,J,K)+TINY)
      SU(I,J,K)=MAX(SU(I,J,K),0.)
 50   CONTINUE
C
      CALL MODTVV
C
C      RESORV=0.
      DO 30 K=2,NKM1
      DO 30 J=2,NJM1
      DO 30 I=2,NIM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*VV(I+1,J,K)+AW(I,J,K)*VV(I-1,J,K)
     1     +AN(I,J,K)*VV(I,J+1,K)+AS(I,J,K)*VV(I,J-1,K)
     1     +AT(I,J,K)*VV(I,J,K+1)+AB(I,J,K)*VV(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*VV(I,J,K)
C
      IF(SP(I,J,K).EQ.-GREAT) RESOR=0.
C      RESORV=RESORV+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFK
      SU(I,J,K)=SU(I,J,K)+(1.-URFK)*AP(I,J,K)*VV(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPK
      CALL LISOLV(3,2,2,NIM1,NJ,NK,VV)
 40   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CALCWW
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)

C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      PR=PRTKE
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*WW(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(1.,PRTKE)
C
      CALL BQUICK(1,1.,WW)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTKE,1,WW)
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C -------- GENERATION AND DISSIPATOIN SOURCE TERMS --------------------
C
      RT=DEN(I,J,K)*TKE(I,J,K)*TKE(I,J,K)/(VISCOS*(TED(I,J,K)+TINY))
      YSTR=DEN(I,J,K)*XL3(I,J,K)*(SQRT(TKE(I,J,K)))/VISCOS
      VISTUR=VIS(I,J,K)-VISCOS
      FE=FLOAT(INWM-1)*EXP(-RT/8.)+FLOAT(2-INWM)*EXP(-YSTR/3.)
C      FE=0.
      TEDT=TED(I,J,K)+FLOAT(INWM-1)*EDNW(I,J,K)
      GNWW=-2.*(UW(I,J,K)*DWDXC(I,J,K)+VW(I,J,K)*DWDYC(I,J,K)+
     1          WW(I,J,K)*DWDZC(I,J,K))
     1      +4.*OMY*UW(I,J,K)-4.*OMX*VW(I,J,K)
C      GNWW=VISTUR*((DWDXC(I,J,K)**2)+(DWDYC(I,J,K)**2)+
C     1          (DWDZC(I,J,K)**2))
      SUGE=VOL(I,J,K)*DEN(I,J,K)*(GNWW-.667*(1.-FE)*TEDT)
      SU(I,J,K)=SU(I,J,K)+MAX(SUGE,0.)
      SP(I,J,K)=SP(I,J,K)-MAX(-SUGE,0.)/(WW(I,J,K)+TINY)
      SP(I,J,K)=SP(I,J,K)-FE*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                    /(TKE(I,J,K)+TINY)
C
C---------   REDISTRIBUTION TERMS  (FIJ1 & FIJ2)  ---------------------
C            FIJ1
      SP(I,J,K)=SP(I,J,K)-CA1*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                     /(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+CA1*DEN(I,J,K)*VOL(I,J,K)*.667*TEDT
C            FIJ2
      FWW2=-CA2*(GNWW-.667*PK1(I,J,K))
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*FWW2
C
C --------- LOW-RE SOURCE TERMS  --------------------------------------
C
      FH=FLOAT(INWM-1)*(10.+2.6*RT)*EXP(-RT/20.)
     1  +FLOAT(2-INWM)*(10.2+7.5*YSTR)*EXP(-.05*YSTR)
      FJ=.06*(FLOAT(INWM-1)*EXP(-RT/8.)+FLOAT(2-INWM)*EXP(-YSTR/3.))
C      FH=0.
C      FJ=0.
C
C          -(HIJ-HJJ/3)    TERM
C
      SP(I,J,K)=SP(I,J,K)-4.*FH*VISCOS*VOL(I,J,K)*DKDZ(I,J,K)*
     1                    DKDZ(I,J,K)/(3.*TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+.667*FH*VISCOS*VOL(I,J,K)*(
     1          DKDX(I,J,K)*DKDX(I,J,K)*UU(I,J,K)
     1         +DKDY(I,J,K)*DKDY(I,J,K)*VV(I,J,K)
     1      +2.*DKDX(I,J,K)*DKDY(I,J,K)*UV(I,J,K)
     1         -DKDX(I,J,K)*DKDZ(I,J,K)*UW(I,J,K)
     1         -DKDY(I,J,K)*DKDZ(I,J,K)*VW(I,J,K))/(TKE(I,J,K)+TINY)
C
C            J  TERM
C 
       SU(I,J,K)=SU(I,J,K)+FJ*DEN(I,J,K)*VOL(I,J,K)*TKE(I,J,K)*
     1                     2.*DWDZC(I,J,K)
C
C------ WALL REFLECTION TERMS  ----------------------------------------
C
      FW=FLOAT(INWM-1)*(1.-EXP(-RT/20.))*(1.+EXP(-RT/100.))
     +  +FLOAT(2-INWM)*(1.-EXP(-.12*YSTR))*(1.+EXP(-0.03*YSTR))
      FW=FW*VOL(I,J,K)*DEN(I,J,K)
C
      GNUU=-2.*(UU(I,J,K)*DUDXC(I,J,K)+UV(I,J,K)*DUDYC(I,J,K)+
     1          UW(I,J,K)*DUDZC(I,J,K))
     1      +4.*OMZ*UV(I,J,K)-4.*OMY*UW(I,J,K)
      GNVV=-2.*(UV(I,J,K)*DVDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
     1          VW(I,J,K)*DVDZC(I,J,K))
     1      +4.*OMX*VW(I,J,K)-4*OMZ*UV(I,J,K)
C      GNWW=-2.*(UW(I,J,K)*DWDXC(I,J,K)+VW(I,J,K)*DWDYC(I,J,K)+
C     1          WW(I,J,K)*DWDZC(I,J,K))
      GNUV=-(VV(I,J,K)*DUDYC(I,J,K)+UU(I,J,K)*DVDXC(I,J,K)+
     1       UV(I,J,K)*(DUDXC(I,J,K)+DVDYC(I,J,K))+
     1       UW(I,J,K)*DVDZC(I,J,K)+VW(I,J,K)*DUDZC(I,J,K))
     1       +2.*OMZ*(VV(I,J,K)-UU(I,J,K))+2.*OMX*UW(I,J,K)
     1       -2.*OMY*VW(I,J,K)
      GNUW=-(WW(I,J,K)*DUDZC(I,J,K)+UU(I,J,K)*DWDXC(I,J,K)+
     1       UW(I,J,K)*(DUDXC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDYC(I,J,K)+VW(I,J,K)*DUDYC(I,J,K))
     1       +2.*OMY*(UU(I,J,K)-WW(I,J,K))+2.*OMZ*VW(I,J,K)
     1       -2.*OMX*UV(I,J,K) 
      GNVW=-(WW(I,J,K)*DVDZC(I,J,K)+VV(I,J,K)*DWDYC(I,J,K)+
     1       VW(I,J,K)*(DVDYC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDXC(I,J,K)+UW(I,J,K)*DVDXC(I,J,K))
     1       +2.*OMX*(WW(I,J,K)-UU(I,J,K))+2.*OMY*UV(I,J,K)
     1       -2.*OMZ*UW(I,J,K)
C
      FUU2=-CA2*(GNUU-.667*PK1(I,J,K))
      FVV2=-CA2*(GNVV-.667*PK1(I,J,K))
C      FWW2=-CA2*(GNWW-.667*PK1(I,J,K))
      FUV2=-CA2*GNUV
      FUW2=-CA2*GNUW
      FVW2=-CA2*GNVW
C
C---------- EAST WALL --------------------------------------------------
C
      RLE=ELE(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLE=MIN(1.,RLE)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLE*CA1W*TED(I,J,K)*2.*(COSZE(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLE*CA1W*TED(I,J,K)*(
     1    (COSXE(I,J,K)**2)*UU(I,J,K)+(COSYE(I,J,K)**2)*VV(I,J,K)
     2    +2.*COSXE(I,J,K)*COSYE(I,J,K)*UV(I,J,K)
     3       -COSXE(I,J,K)*COSZE(I,J,K)*UW(I,J,K)
     4       -COSYE(I,J,K)*COSZE(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLE*CA2W*(-2.*(COSZE(I,J,K)**2)*FWW2    
     1    +(COSXE(I,J,K)**2)*FUU2+(COSYE(I,J,K)**2)*FVV2
     2      +2.*COSXE(I,J,K)*COSYE(I,J,K)*FUV2
     3         -COSXE(I,J,K)*COSZE(I,J,K)*FUW2
     4         -COSYE(I,J,K)*COSZE(I,J,K)*FVW2           )
C
C---------- WEST WALL --------------------------------------------------
C
      RLW=ELW(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLW=MIN(1.,RLW)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLW*CA1W*TED(I,J,K)*2.*(COSZW(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLW*CA1W*TED(I,J,K)*(
     1    (COSXW(I,J,K)**2)*UU(I,J,K)+(COSYW(I,J,K)**2)*VV(I,J,K)
     2    +2.*COSXW(I,J,K)*COSYW(I,J,K)*UV(I,J,K)
     3       -COSXW(I,J,K)*COSZW(I,J,K)*UW(I,J,K)
     4       -COSYW(I,J,K)*COSZW(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLW*CA2W*(-2.*(COSZW(I,J,K)**2)*FWW2    
     1    +(COSXW(I,J,K)**2)*FUU2+(COSYW(I,J,K)**2)*FVV2
     2      +2.*COSXW(I,J,K)*COSYW(I,J,K)*FUV2
     3         -COSXW(I,J,K)*COSZW(I,J,K)*FUW2
     4         -COSYW(I,J,K)*COSZW(I,J,K)*FVW2           )
C
C---------- NORTH WALL -------------------------------------------------
C
      RLN=ELN(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLN=MIN(1.,RLN)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLN*CA1W*TED(I,J,K)*2.*(COSZN(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLN*CA1W*TED(I,J,K)*(
     1    (COSXN(I,J,K)**2)*UU(I,J,K)+(COSYN(I,J,K)**2)*VV(I,J,K)
     2    +2.*COSXN(I,J,K)*COSYN(I,J,K)*UV(I,J,K)
     3       -COSXN(I,J,K)*COSZN(I,J,K)*UW(I,J,K)
     4       -COSYN(I,J,K)*COSZN(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLN*CA2W*(-2.*(COSZN(I,J,K)**2)*FWW2    
     1    +(COSXN(I,J,K)**2)*FUU2+(COSYN(I,J,K)**2)*FVV2
     2      +2.*COSXN(I,J,K)*COSYN(I,J,K)*FUV2
     3         -COSXN(I,J,K)*COSZN(I,J,K)*FUW2
     4         -COSYN(I,J,K)*COSZN(I,J,K)*FVW2           )
C
C---------- SOUTH WALL -------------------------------------------------
C
      RLS=ELS(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLS=MIN(1.,RLS)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLS*CA1W*TED(I,J,K)*2.*(COSZS(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLS*CA1W*TED(I,J,K)*(
     1    (COSXS(I,J,K)**2)*UU(I,J,K)+(COSYS(I,J,K)**2)*VV(I,J,K)
     2    +2.*COSXS(I,J,K)*COSYS(I,J,K)*UV(I,J,K)
     3       -COSXS(I,J,K)*COSZS(I,J,K)*UW(I,J,K)
     4       -COSYS(I,J,K)*COSZS(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLS*CA2W*(-2.*(COSZS(I,J,K)**2)*FWW2    
     1    +(COSXS(I,J,K)**2)*FUU2+(COSYS(I,J,K)**2)*FVV2
     2      +2.*COSXS(I,J,K)*COSYS(I,J,K)*FUV2
     3         -COSXS(I,J,K)*COSZS(I,J,K)*FUW2
     4         -COSYS(I,J,K)*COSZS(I,J,K)*FVW2           )
C
C---------- TOP WALL -------------------------------------------------
C
      RLT=ELT(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLT=MIN(1.,RLT)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLT*CA1W*TED(I,J,K)*2.*(COSZT(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLT*CA1W*TED(I,J,K)*(
     1    (COSXT(I,J,K)**2)*UU(I,J,K)+(COSYT(I,J,K)**2)*VV(I,J,K)
     2    +2.*COSXT(I,J,K)*COSYT(I,J,K)*UV(I,J,K)
     3       -COSXT(I,J,K)*COSZT(I,J,K)*UW(I,J,K)
     4       -COSYT(I,J,K)*COSZT(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLT*CA2W*(-2.*(COSZT(I,J,K)**2)*FWW2    
     1    +(COSXT(I,J,K)**2)*FUU2+(COSYT(I,J,K)**2)*FVV2
     2      +2.*COSXT(I,J,K)*COSYT(I,J,K)*FUV2
     3         -COSXT(I,J,K)*COSZT(I,J,K)*FUW2
     4         -COSYT(I,J,K)*COSZT(I,J,K)*FVW2           )
C
C---------- BOTTOM WALL ----------------------------------------------
C
      RLB=ELB(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLB=MIN(1.,RLB)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLB*CA1W*TED(I,J,K)*2.*(COSZB(I,J,K)**2)
     1          /(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)+FW*RLB*CA1W*TED(I,J,K)*(
     1    (COSXB(I,J,K)**2)*UU(I,J,K)+(COSYB(I,J,K)**2)*VV(I,J,K)
     2    +2.*COSXB(I,J,K)*COSYB(I,J,K)*UV(I,J,K)
     3       -COSXB(I,J,K)*COSZB(I,J,K)*UW(I,J,K)
     4       -COSYB(I,J,K)*COSZB(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)+FW*RLB*CA2W*(-2.*(COSZB(I,J,K)**2)*FWW2    
     1    +(COSXB(I,J,K)**2)*FUU2+(COSYB(I,J,K)**2)*FVV2
     2      +2.*COSXB(I,J,K)*COSYB(I,J,K)*FUV2
     3         -COSXB(I,J,K)*COSZB(I,J,K)*FUW2
     4         -COSYB(I,J,K)*COSZB(I,J,K)*FVW2           )
C
C ------ OVERAL SOURCE TREATMENT FOR POSITIVE ANSWERS ------------------
C
      SP(I,J,K)=SP(I,J,K)-MAX(-SU(I,J,K),0.)/(WW(I,J,K)+TINY)
      SU(I,J,K)=MAX(SU(I,J,K),0.)
 50   CONTINUE
C
      CALL MODTWW
C
C      RESORW=0.
      DO 30 K=2,NKM1
      DO 30 J=2,NJM1
      DO 30 I=2,NIM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*WW(I+1,J,K)+AW(I,J,K)*WW(I-1,J,K)
     1     +AN(I,J,K)*WW(I,J+1,K)+AS(I,J,K)*WW(I,J-1,K)
     1     +AT(I,J,K)*WW(I,J,K+1)+AB(I,J,K)*WW(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*WW(I,J,K)
C
      IF(SP(I,J,K).EQ.-GREAT) RESOR=0.
C      RESORW=RESORW+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFK
      SU(I,J,K)=SU(I,J,K)+(1.-URFK)*AP(I,J,K)*WW(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPK
      CALL LISOLV(3,2,2,NIM1,NJ,NK,WW)
 40   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CALCUV
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      PR=PRTKE
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*UV(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(1.,PRTKE)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTKE,0,UV)
C
      CALL BQUICK(0,1.,UV)
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C -------- GENERATION AND DISSIPATOIN SOURCE TERMS --------------------
C
      RT=DEN(I,J,K)*TKE(I,J,K)*TKE(I,J,K)/(VISCOS*(TED(I,J,K)+TINY))
      YSTR=DEN(I,J,K)*XL3(I,J,K)*(SQRT(TKE(I,J,K)))/VISCOS
      FE=FLOAT(INWM-1)*EXP(-RT/8.)+FLOAT(2-INWM)*EXP(-YSTR/3.)
C      FE=0.
      TEDT=TED(I,J,K)+FLOAT(INWM-1)*EDNW(I,J,K)
      GNUV=-(VV(I,J,K)*DUDYC(I,J,K)+UU(I,J,K)*DVDXC(I,J,K)+
     1       UV(I,J,K)*(DUDXC(I,J,K)+DVDYC(I,J,K))+
     1       UW(I,J,K)*DVDZC(I,J,K)+VW(I,J,K)*DUDZC(I,J,K))
     1       +2.*OMZ*(VV(I,J,K)-UU(I,J,K))+2.*OMX*UW(I,J,K)
     1       -2.*OMY*VW(I,J,K)
      SU(I,J,K)=SU(I,J,K)+VOL(I,J,K)*DEN(I,J,K)*GNUV
      SP(I,J,K)=SP(I,J,K)-FE*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                    /(TKE(I,J,K)+TINY)
C
C---------   REDISTRIBUTION TERMS  (FIJ1 & FIJ2)  ---------------------
C            FIJ1
      SP(I,J,K)=SP(I,J,K)-CA1*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                     /(TKE(I,J,K)+TINY)
C            FIJ2
      FUV2=-CA2*GNUV
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*FUV2
C
C --------- LOW-RE SOURCE TERMS  --------------------------------------
C
      FH=FLOAT(INWM-1)*(10.+2.6*RT)*EXP(-RT/20.)
     1  +FLOAT(2-INWM)*(10.2+7.5*YSTR)*EXP(-.05*YSTR)
      FJ=.06*(FLOAT(INWM-1)*EXP(-RT/8)+FLOAT(2-INWM)*EXP(-YSTR/3.))
C
C      FH=0.
C      FJ=0.
C
C          -(HIJ-2*HJJ/3)    TERM
C
      SP(I,J,K)=SP(I,J,K)-FH*VISCOS*VOL(I,J,K)*(DKDX(I,J,K)*
     1    DKDX(I,J,K)+DKDY(I,J,K)*DKDY(I,J,K))/(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)-FH*VISCOS*VOL(I,J,K)*(
     1    DKDX(I,J,K)*DKDY(I,J,K)*(UU(I,J,K)+VV(I,J,K))
     1         +DKDX(I,J,K)*DKDZ(I,J,K)*VW(I,J,K)
     1         +DKDY(I,J,K)*DKDZ(I,J,K)*UW(I,J,K))/(TKE(I,J,K)+TINY)
C
C            J  TERM
C 
       SU(I,J,K)=SU(I,J,K)+FJ*DEN(I,J,K)*VOL(I,J,K)*TKE(I,J,K)*
     1                     (DUDYC(I,J,K)+DVDXC(I,J,K))
C
C------ WALL REFLECTION TERMS  ----------------------------------------
C
      FW=FLOAT(INWM-1)*(1.-EXP(-RT/20.))*(1.+EXP(-RT/100.))
     +  +FLOAT(2-INWM)*(1.-EXP(-.12*YSTR))*(1.+EXP(-0.03*YSTR))
      FW=FW*VOL(I,J,K)*DEN(I,J,K)
C
      GNUU=-2.*(UU(I,J,K)*DUDXC(I,J,K)+UV(I,J,K)*DUDYC(I,J,K)+
     1          UW(I,J,K)*DUDZC(I,J,K))
     1      +4.*OMZ*UV(I,J,K)-4.*OMY*UW(I,J,K)
      GNVV=-2.*(UV(I,J,K)*DVDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
     1          VW(I,J,K)*DVDZC(I,J,K))
     1      +4.*OMX*VW(I,J,K)-4*OMZ*UV(I,J,K)
      GNWW=-2.*(UW(I,J,K)*DWDXC(I,J,K)+VW(I,J,K)*DWDYC(I,J,K)+
     1          WW(I,J,K)*DWDZC(I,J,K))
     1      +4.*OMY*UW(I,J,K)-4.*OMX*VW(I,J,K)
C      GNUV=-(VV(I,J,K)*DUDYC(I,J,K)+UU(I,J,K)*DVDXC(I,J,K)+
C     1       UV(I,J,K)*(DUDXC(I,J,K)+DVDYC(I,J,K))+
C     1       UW(I,J,K)*DVDZC(I,J,K)+VW(I,J,K)*DUDZC(I,J,K))
      GNUW=-(WW(I,J,K)*DUDZC(I,J,K)+UU(I,J,K)*DWDXC(I,J,K)+
     1       UW(I,J,K)*(DUDXC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDYC(I,J,K)+VW(I,J,K)*DUDYC(I,J,K))
     1       +2.*OMY*(UU(I,J,K)-WW(I,J,K))+2.*OMZ*VW(I,J,K)
     1       -2.*OMX*UV(I,J,K) 
      GNVW=-(WW(I,J,K)*DVDZC(I,J,K)+VV(I,J,K)*DWDYC(I,J,K)+
     1       VW(I,J,K)*(DVDYC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDXC(I,J,K)+UW(I,J,K)*DVDXC(I,J,K))
     1       +2.*OMX*(WW(I,J,K)-UU(I,J,K))+2.*OMY*UV(I,J,K)
     1       -2.*OMZ*UW(I,J,K)
C
      FUU2=-CA2*(GNUU-.667*PK1(I,J,K))
      FVV2=-CA2*(GNVV-.667*PK1(I,J,K))
      FWW2=-CA2*(GNWW-.667*PK1(I,J,K))
C      FUV2=-CA2*GNUV
      FUW2=-CA2*GNUW
      FVW2=-CA2*GNVW
C
C---------- EAST WALL --------------------------------------------------
C
      RLE=ELE(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLE=MIN(1.,RLE)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLE*CA1W*TED(I,J,K)*1.5*((COSXE(I,J,K)**2)
     1          +(COSYE(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLE*CA1W*TED(I,J,K)*(
     1    COSXE(I,J,K)*COSYE(I,J,K)*(UU(I,J,K)+VV(I,J,K))
     3       +COSYE(I,J,K)*COSZE(I,J,K)*UW(I,J,K)
     4       +COSXE(I,J,K)*COSZE(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLE*CA2W*(    
     1        +COSXE(I,J,K)*COSYE(I,J,K)*(FUU2+FVV2)
     2      +((COSXE(I,J,K)**2)+(COSYE(I,J,K)**2))*FUV2
     3         +COSYE(I,J,K)*COSZE(I,J,K)*FUW2
     4         +COSXE(I,J,K)*COSZE(I,J,K)*FVW2           )
C
C---------- WEST WALL --------------------------------------------------
C
      RLW=ELW(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLW=MIN(1.,RLW)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLW*CA1W*TED(I,J,K)*1.5*((COSXW(I,J,K)**2)
     1          +(COSYW(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLW*CA1W*TED(I,J,K)*(
     1    COSXW(I,J,K)*COSYW(I,J,K)*(UU(I,J,K)+VV(I,J,K))
     3       +COSYW(I,J,K)*COSZW(I,J,K)*UW(I,J,K)
     4       +COSXW(I,J,K)*COSZW(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLW*CA2W*(    
     1        +COSXW(I,J,K)*COSYW(I,J,K)*(FUU2+FVV2)
     2      +((COSXW(I,J,K)**2)+(COSYW(I,J,K)**2))*FUV2
     3         +COSYW(I,J,K)*COSZW(I,J,K)*FUW2
     4         +COSXW(I,J,K)*COSZW(I,J,K)*FVW2           )
C
C---------- NORTH WALL -------------------------------------------------
C
      RLN=ELN(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLN=MIN(1.,RLN)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLN*CA1W*TED(I,J,K)*1.5*((COSXN(I,J,K)**2)
     1          +(COSYN(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLN*CA1W*TED(I,J,K)*(
     1    COSXN(I,J,K)*COSYN(I,J,K)*(UU(I,J,K)+VV(I,J,K))
     3       +COSYN(I,J,K)*COSZN(I,J,K)*UW(I,J,K)
     4       +COSXN(I,J,K)*COSZN(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLN*CA2W*(    
     1        +COSXN(I,J,K)*COSYN(I,J,K)*(FUU2+FVV2)
     2      +((COSXN(I,J,K)**2)+(COSYN(I,J,K)**2))*FUV2
     3         +COSYN(I,J,K)*COSZN(I,J,K)*FUW2
     4         +COSXN(I,J,K)*COSZN(I,J,K)*FVW2           )
C
C---------- SOUTH WALL -------------------------------------------------
C
      RLS=ELS(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLS=MIN(1.,RLS)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLS*CA1W*TED(I,J,K)*1.5*((COSXS(I,J,K)**2)
     1          +(COSYS(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLS*CA1W*TED(I,J,K)*(
     1    COSXS(I,J,K)*COSYS(I,J,K)*(UU(I,J,K)+VV(I,J,K))
     3       +COSYS(I,J,K)*COSZS(I,J,K)*UW(I,J,K)
     4       +COSXS(I,J,K)*COSZS(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLS*CA2W*(    
     1        +COSXS(I,J,K)*COSYS(I,J,K)*(FUU2+FVV2)
     2      +((COSXS(I,J,K)**2)+(COSYS(I,J,K)**2))*FUV2
     3         +COSYS(I,J,K)*COSZS(I,J,K)*FUW2
     4         +COSXS(I,J,K)*COSZS(I,J,K)*FVW2           )
C
C---------- TOP WALL -------------------------------------------------
C
      RLT=ELT(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLT=MIN(1.,RLT)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLT*CA1W*TED(I,J,K)*1.5*((COSXT(I,J,K)**2)
     1          +(COSYT(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLT*CA1W*TED(I,J,K)*(
     1    COSXT(I,J,K)*COSYT(I,J,K)*(UU(I,J,K)+VV(I,J,K))
     3       +COSYT(I,J,K)*COSZT(I,J,K)*UW(I,J,K)
     4       +COSXT(I,J,K)*COSZT(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLT*CA2W*(    
     1        +COSXT(I,J,K)*COSYT(I,J,K)*(FUU2+FVV2)
     2      +((COSXT(I,J,K)**2)+(COSYT(I,J,K)**2))*FUV2
     3         +COSYT(I,J,K)*COSZT(I,J,K)*FUW2
     4         +COSXT(I,J,K)*COSZT(I,J,K)*FVW2           )
C
C---------- BOTTOM WALL ----------------------------------------------
C
      RLB=ELB(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLB=MIN(1.,RLB)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLB*CA1W*TED(I,J,K)*1.5*((COSXB(I,J,K)**2)
     1          +(COSYB(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLB*CA1W*TED(I,J,K)*(
     1    COSXB(I,J,K)*COSYB(I,J,K)*(UU(I,J,K)+VV(I,J,K))
     3       +COSYB(I,J,K)*COSZB(I,J,K)*UW(I,J,K)
     4       +COSXB(I,J,K)*COSZB(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLB*CA2W*(    
     1        +COSXB(I,J,K)*COSYB(I,J,K)*(FUU2+FVV2)
     2      +((COSXB(I,J,K)**2)+(COSYB(I,J,K)**2))*FUV2
     3         +COSYB(I,J,K)*COSZB(I,J,K)*FUW2
     4         +COSXB(I,J,K)*COSZB(I,J,K)*FVW2           )
C
 50   CONTINUE
C
      CALL MODTUV
C
C      RESORK=0.
      DO 30 I=2,NIM1
      DO 30 J=2,NJM1
      DO 30 K=2,NKM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*UV(I+1,J,K)+AW(I,J,K)*UV(I-1,J,K)
     1     +AN(I,J,K)*UV(I,J+1,K)+AS(I,J,K)*UV(I,J-1,K)
     1     +AT(I,J,K)*UV(I,J,K+1)+AB(I,J,K)*UV(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*UV(I,J,K)
C
      IF(SP(I,J,K).EQ.-GREAT) RESOR=0.
C      RESORU=RESORU+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFK
      SU(I,J,K)=SU(I,J,K)+(1.-URFK)*AP(I,J,K)*UV(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPK
      CALL LISOLV(3,2,2,NIM1,NJ,NK,UV)
 40   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CALCUW
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      PR=PRTKE
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*UW(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(1.,PRTKE)
C
      CALL BQUICK(0,1.,UW)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTKE,0,UW)
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C -------- GENERATION AND DISSIPATOIN SOURCE TERMS --------------------
C
      RT=DEN(I,J,K)*TKE(I,J,K)*TKE(I,J,K)/(VISCOS*(TED(I,J,K)+TINY))
      YSTR=DEN(I,J,K)*XL3(I,J,K)*(SQRT(TKE(I,J,K)))/VISCOS
      FE=FLOAT(INWM-1)*EXP(-RT/8.)+FLOAT(2-INWM)*EXP(-YSTR/3.)
C      FE=0.
      TEDT=TED(I,J,K)+FLOAT(INWM-1)*EDNW(I,J,K)
      GNUW=-(WW(I,J,K)*DUDZC(I,J,K)+UU(I,J,K)*DWDXC(I,J,K)+
     1       UW(I,J,K)*(DUDXC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDYC(I,J,K)+VW(I,J,K)*DUDYC(I,J,K))
     1       +2.*OMY*(UU(I,J,K)-WW(I,J,K))+2.*OMZ*VW(I,J,K)
     1       -2.*OMX*UV(I,J,K) 
      SU(I,J,K)=SU(I,J,K)+VOL(I,J,K)*DEN(I,J,K)*GNUW
      SP(I,J,K)=SP(I,J,K)-FE*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                    /(TKE(I,J,K)+TINY)
C
C---------   REDISTRIBUTION TERMS  (FIJ1 & FIJ2)  ---------------------
C            FIJ1
      SP(I,J,K)=SP(I,J,K)-CA1*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                     /(TKE(I,J,K)+TINY)
C            FIJ2
      FUW2=-CA2*GNUW
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*FUW2
C
C --------- LOW-RE SOURCE TERMS  --------------------------------------
C
      FH=FLOAT(INWM-1)*(10.+2.6*RT)*EXP(-RT/20.)
     1  +FLOAT(2-INWM)*(10.2+7.5*YSTR)*EXP(-.05*YSTR)
      FJ=.06*(FLOAT(INWM-1)*EXP(-RT/8)+FLOAT(2-INWM)*EXP(-YSTR/3.))
C
C      FH=0.
C      FJ=0.
C
C          -(HIJ-2*HJJ/3)    TERM
C
      SP(I,J,K)=SP(I,J,K)-FH*VISCOS*VOL(I,J,K)*(DKDX(I,J,K)*
     1    DKDX(I,J,K)+DKDZ(I,J,K)*DKDZ(I,J,K))/(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)-FH*VISCOS*VOL(I,J,K)*(
     1    DKDX(I,J,K)*DKDZ(I,J,K)*(UU(I,J,K)+WW(I,J,K))
     1         +DKDX(I,J,K)*DKDY(I,J,K)*VW(I,J,K)
     1         +DKDY(I,J,K)*DKDZ(I,J,K)*UV(I,J,K))/(TKE(I,J,K)+TINY)
C
C            J  TERM
C 
       SU(I,J,K)=SU(I,J,K)+FJ*DEN(I,J,K)*VOL(I,J,K)*TKE(I,J,K)*
     1                     (DUDZC(I,J,K)+DWDXC(I,J,K))
C
C------ WALL REFLECTION TERMS  ----------------------------------------
C
      FW=FLOAT(INWM-1)*(1.-EXP(-RT/20.))*(1.+EXP(-RT/100.))
     +  +FLOAT(2-INWM)*(1.-EXP(-.12*YSTR))*(1.+EXP(-0.03*YSTR))
      FW=FW*VOL(I,J,K)*DEN(I,J,K)
C
      GNUU=-2.*(UU(I,J,K)*DUDXC(I,J,K)+UV(I,J,K)*DUDYC(I,J,K)+
     1          UW(I,J,K)*DUDZC(I,J,K))
     1      +4.*OMZ*UV(I,J,K)-4.*OMY*UW(I,J,K)
      GNVV=-2.*(UV(I,J,K)*DVDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
     1          VW(I,J,K)*DVDZC(I,J,K))
     1      +4.*OMX*VW(I,J,K)-4*OMZ*UV(I,J,K)
      GNWW=-2.*(UW(I,J,K)*DWDXC(I,J,K)+VW(I,J,K)*DWDYC(I,J,K)+
     1          WW(I,J,K)*DWDZC(I,J,K))
     1      +4.*OMY*UW(I,J,K)-4.*OMX*VW(I,J,K)
      GNUV=-(VV(I,J,K)*DUDYC(I,J,K)+UU(I,J,K)*DVDXC(I,J,K)+
     1       UV(I,J,K)*(DUDXC(I,J,K)+DVDYC(I,J,K))+
     1       UW(I,J,K)*DVDZC(I,J,K)+VW(I,J,K)*DUDZC(I,J,K))
     1       +2.*OMZ*(VV(I,J,K)-UU(I,J,K))+2.*OMX*UW(I,J,K)
     1       -2.*OMY*VW(I,J,K)
C      GNUW=-(WW(I,J,K)*DUDZC(I,J,K)+UU(I,J,K)*DWDXC(I,J,K)+
C     1       UW(I,J,K)*(DUDXC(I,J,K)+DWDZC(I,J,K))+
C     1       UV(I,J,K)*DWDYC(I,J,K)+VW(I,J,K)*DUDYC(I,J,K))
      GNVW=-(WW(I,J,K)*DVDZC(I,J,K)+VV(I,J,K)*DWDYC(I,J,K)+
     1       VW(I,J,K)*(DVDYC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDXC(I,J,K)+UW(I,J,K)*DVDXC(I,J,K))
     1       +2.*OMX*(WW(I,J,K)-UU(I,J,K))+2.*OMY*UV(I,J,K)
     1       -2.*OMZ*UW(I,J,K)
C
      FUU2=-CA2*(GNUU-.667*PK1(I,J,K))
      FVV2=-CA2*(GNVV-.667*PK1(I,J,K))
      FWW2=-CA2*(GNWW-.667*PK1(I,J,K))
      FUV2=-CA2*GNUV
C      FUW2=-CA2*GNUW
      FVW2=-CA2*GNVW
C
C---------- EAST WALL --------------------------------------------------
C
      RLE=ELE(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLE=MIN(1.,RLE)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLE*CA1W*TED(I,J,K)*1.5*((COSXE(I,J,K)**2)
     1          +(COSZE(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLE*CA1W*TED(I,J,K)*(
     1    COSXE(I,J,K)*COSZE(I,J,K)*(UU(I,J,K)+WW(I,J,K))
     3       +COSYE(I,J,K)*COSZE(I,J,K)*UV(I,J,K)
     4       +COSXE(I,J,K)*COSYE(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLE*CA2W*(    
     1        +COSXE(I,J,K)*COSZE(I,J,K)*(FUU2+FWW2)
     2      +((COSXE(I,J,K)**2)+(COSZE(I,J,K)**2))*FUW2
     3         +COSYE(I,J,K)*COSZE(I,J,K)*FUV2
     4         +COSXE(I,J,K)*COSYE(I,J,K)*FVW2           )
C
C---------- WEST WALL --------------------------------------------------
C
      RLW=ELW(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLW=MIN(1.,RLW)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLW*CA1W*TED(I,J,K)*1.5*((COSXW(I,J,K)**2)
     1          +(COSZW(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLW*CA1W*TED(I,J,K)*(
     1    COSXW(I,J,K)*COSZW(I,J,K)*(UU(I,J,K)+WW(I,J,K))
     3       +COSYW(I,J,K)*COSZW(I,J,K)*UV(I,J,K)
     4       +COSXW(I,J,K)*COSYW(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLW*CA2W*(    
     1        +COSXW(I,J,K)*COSZW(I,J,K)*(FUU2+FWW2)
     2      +((COSXW(I,J,K)**2)+(COSZW(I,J,K)**2))*FUW2
     3         +COSYW(I,J,K)*COSZW(I,J,K)*FUV2
     4         +COSXW(I,J,K)*COSYW(I,J,K)*FVW2           )
C
C---------- NORTH WALL -------------------------------------------------
C
      RLN=ELN(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLN=MIN(1.,RLN)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLN*CA1W*TED(I,J,K)*1.5*((COSXN(I,J,K)**2)
     1          +(COSZN(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLN*CA1W*TED(I,J,K)*(
     1    COSXN(I,J,K)*COSZN(I,J,K)*(UU(I,J,K)+WW(I,J,K))
     3       +COSYN(I,J,K)*COSZN(I,J,K)*UV(I,J,K)
     4       +COSXN(I,J,K)*COSYN(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLN*CA2W*(    
     1        +COSXN(I,J,K)*COSZN(I,J,K)*(FUU2+FWW2)
     2      +((COSXN(I,J,K)**2)+(COSZN(I,J,K)**2))*FUW2
     3         +COSYN(I,J,K)*COSZN(I,J,K)*FUV2
     4         +COSXN(I,J,K)*COSYN(I,J,K)*FVW2           )
C
C---------- SOUTH WALL -------------------------------------------------
C
      RLS=ELS(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLS=MIN(1.,RLS)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLS*CA1W*TED(I,J,K)*1.5*((COSXS(I,J,K)**2)
     1          +(COSZS(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLS*CA1W*TED(I,J,K)*(
     1    COSXS(I,J,K)*COSZS(I,J,K)*(UU(I,J,K)+WW(I,J,K))
     3       +COSYS(I,J,K)*COSZS(I,J,K)*UV(I,J,K)
     4       +COSXS(I,J,K)*COSYS(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLS*CA2W*(    
     1        +COSXS(I,J,K)*COSZS(I,J,K)*(FUU2+FWW2)
     2      +((COSXS(I,J,K)**2)+(COSZS(I,J,K)**2))*FUW2
     3         +COSYS(I,J,K)*COSZS(I,J,K)*FUV2
     4         +COSXS(I,J,K)*COSYS(I,J,K)*FVW2           )
C
C---------- TOP WALL -------------------------------------------------
C
      RLT=ELT(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLT=MIN(1.,RLT)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLT*CA1W*TED(I,J,K)*1.5*((COSXT(I,J,K)**2)
     1          +(COSZT(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLT*CA1W*TED(I,J,K)*(
     1    COSXT(I,J,K)*COSZT(I,J,K)*(UU(I,J,K)+WW(I,J,K))
     3       +COSYT(I,J,K)*COSZT(I,J,K)*UV(I,J,K)
     4       +COSXT(I,J,K)*COSYT(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLT*CA2W*(    
     1        +COSXT(I,J,K)*COSZT(I,J,K)*(FUU2+FWW2)
     2      +((COSXT(I,J,K)**2)+(COSZT(I,J,K)**2))*FUW2
     3         +COSYT(I,J,K)*COSZT(I,J,K)*FUV2
     4         +COSXT(I,J,K)*COSYT(I,J,K)*FVW2           )
C
C---------- BOTTOM WALL ----------------------------------------------
C
      RLB=ELB(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLB=MIN(1.,RLB)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLB*CA1W*TED(I,J,K)*1.5*((COSXB(I,J,K)**2)
     1          +(COSZB(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLB*CA1W*TED(I,J,K)*(
     1    COSXB(I,J,K)*COSZB(I,J,K)*(UU(I,J,K)+WW(I,J,K))
     3       +COSYB(I,J,K)*COSZB(I,J,K)*UV(I,J,K)
     4       +COSXB(I,J,K)*COSZB(I,J,K)*VW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLB*CA2W*(    
     1        +COSXB(I,J,K)*COSZB(I,J,K)*(FUU2+FWW2)
     2      +((COSXB(I,J,K)**2)+(COSZB(I,J,K)**2))*FUW2
     3         +COSYB(I,J,K)*COSZB(I,J,K)*FUV2
     4         +COSXB(I,J,K)*COSYB(I,J,K)*FVW2           )
C
 50   CONTINUE
C
      CALL MODTUW
C
C      RESORK=0.
      DO 30 I=2,NIM1
      DO 30 J=2,NJM1
      DO 30 K=2,NKM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*UW(I+1,J,K)+AW(I,J,K)*UW(I-1,J,K)
     1     +AN(I,J,K)*UW(I,J+1,K)+AS(I,J,K)*UW(I,J-1,K)
     1     +AT(I,J,K)*UW(I,J,K+1)+AB(I,J,K)*UW(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*UW(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
C      RESORV=RESORV+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFK
      SU(I,J,K)=SU(I,J,K)+(1.-URFK)*AP(I,J,K)*UW(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPK
      CALL LISOLV(3,2,2,NIM1,NJ,NK,UW)
 40   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CALCVW
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      PR=PRTKE
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*VW(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(1.,PRTKE)
C
      CALL BQUICK(0,1.,VW)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTKE,0,VW)
C
      DO 50 K=2,NKM1
      DO 50 J=2,NJM1
      DO 50 I=2,NIM1
C
C -------- GENERATION AND DISSIPATOIN SOURCE TERMS --------------------
C
      RT=DEN(I,J,K)*TKE(I,J,K)*TKE(I,J,K)/(VISCOS*(TED(I,J,K)+TINY))
      YSTR=DEN(I,J,K)*XL3(I,J,K)*(SQRT(TKE(I,J,K)))/VISCOS
      FE=FLOAT(INWM-1)*EXP(-RT/8.)+FLOAT(2-INWM)*EXP(-YSTR/3.)
C      FE=0.
      TEDT=TED(I,J,K)+FLOAT(INWM-1)*EDNW(I,J,K)
      GNVW=-(WW(I,J,K)*DVDZC(I,J,K)+VV(I,J,K)*DWDYC(I,J,K)+
     1       VW(I,J,K)*(DVDYC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDXC(I,J,K)+UW(I,J,K)*DVDXC(I,J,K))
     1       +2.*OMX*(WW(I,J,K)-UU(I,J,K))+2.*OMY*UV(I,J,K)
     1       -2.*OMZ*UW(I,J,K)
      SU(I,J,K)=SU(I,J,K)+VOL(I,J,K)*DEN(I,J,K)*GNVW
      SP(I,J,K)=SP(I,J,K)-FE*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                    /(TKE(I,J,K)+TINY)
C
C---------   REDISTRIBUTION TERMS  (FIJ1 & FIJ2)  ---------------------
C            FIJ1
      SP(I,J,K)=SP(I,J,K)-CA1*DEN(I,J,K)*VOL(I,J,K)*TEDT
     1                     /(TKE(I,J,K)+TINY)
C            FIJ2
      FVW2=-CA2*GNVW
      SU(I,J,K)=SU(I,J,K)+DEN(I,J,K)*VOL(I,J,K)*FVW2
C
C --------- LOW-RE SOURCE TERMS  --------------------------------------
C
      FH=FLOAT(INWM-1)*(10.+2.6*RT)*EXP(-RT/20.)
     1  +FLOAT(2-INWM)*(10.2+7.5*YSTR)*EXP(-.05*YSTR)
      FJ=.06*(FLOAT(INWM-1)*EXP(-RT/8)+FLOAT(2-INWM)*EXP(-YSTR/3.))
C
C      FH=0.
C      FJ=0.
C
C          -(HIJ-2*HJJ/3)    TERM
C
      SP(I,J,K)=SP(I,J,K)-FH*VISCOS*VOL(I,J,K)*(DKDY(I,J,K)*
     1    DKDY(I,J,K)+DKDZ(I,J,K)*DKDZ(I,J,K))/(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)-FH*VISCOS*VOL(I,J,K)*(
     1    DKDY(I,J,K)*DKDZ(I,J,K)*(VV(I,J,K)+WW(I,J,K))
     1         +DKDX(I,J,K)*DKDY(I,J,K)*UW(I,J,K)
     1         +DKDX(I,J,K)*DKDZ(I,J,K)*UV(I,J,K))/(TKE(I,J,K)+TINY)
C
C            J  TERM
C 
       SU(I,J,K)=SU(I,J,K)+FJ*DEN(I,J,K)*VOL(I,J,K)*TKE(I,J,K)*
     1                     (DVDZC(I,J,K)+DWDYC(I,J,K))
C
C------ WALL REFLECTION TERMS  ----------------------------------------
C
      FW=FLOAT(INWM-1)*(1.-EXP(-RT/20.))*(1.+EXP(-RT/100.))
     +  +FLOAT(2-INWM)*(1.-EXP(-.12*YSTR))*(1.+EXP(-0.03*YSTR))
      FW=FW*VOL(I,J,K)*DEN(I,J,K)
C
      GNUU=-2.*(UU(I,J,K)*DUDXC(I,J,K)+UV(I,J,K)*DUDYC(I,J,K)+
     1          UW(I,J,K)*DUDZC(I,J,K))
     1      +4.*OMZ*UV(I,J,K)-4.*OMY*UW(I,J,K)
      GNVV=-2.*(UV(I,J,K)*DVDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
     1          VW(I,J,K)*DVDZC(I,J,K))
     1      +4.*OMX*VW(I,J,K)-4*OMZ*UV(I,J,K)
      GNWW=-2.*(UW(I,J,K)*DWDXC(I,J,K)+VW(I,J,K)*DWDYC(I,J,K)+
     1          WW(I,J,K)*DWDZC(I,J,K))
     1      +4.*OMY*UW(I,J,K)-4.*OMX*VW(I,J,K)
      GNUV=-(VV(I,J,K)*DUDYC(I,J,K)+UU(I,J,K)*DVDXC(I,J,K)+
     1       UV(I,J,K)*(DUDXC(I,J,K)+DVDYC(I,J,K))+
     1       UW(I,J,K)*DVDZC(I,J,K)+VW(I,J,K)*DUDZC(I,J,K))
     1       +2.*OMZ*(VV(I,J,K)-UU(I,J,K))+2.*OMX*UW(I,J,K)
     1       -2.*OMY*VW(I,J,K)
      GNUW=-(WW(I,J,K)*DUDZC(I,J,K)+UU(I,J,K)*DWDXC(I,J,K)+
     1       UW(I,J,K)*(DUDXC(I,J,K)+DWDZC(I,J,K))+
     1       UV(I,J,K)*DWDYC(I,J,K)+VW(I,J,K)*DUDYC(I,J,K))
     1       +2.*OMY*(UU(I,J,K)-WW(I,J,K))+2.*OMZ*VW(I,J,K)
     1       -2.*OMX*UV(I,J,K) 
C      GNVW=-(WW(I,J,K)*DVDZC(I,J,K)+VV(I,J,K)*DWDYC(I,J,K)+
C     1       VW(I,J,K)*(DVDYC(I,J,K)+DWDZC(I,J,K))+
C     1       UV(I,J,K)*DWDXC(I,J,K)+UW(I,J,K)*DVDXC(I,J,K))
C
      FUU2=-CA2*(GNUU-.667*PK1(I,J,K))
      FVV2=-CA2*(GNVV-.667*PK1(I,J,K))
      FWW2=-CA2*(GNWW-.667*PK1(I,J,K))
      FUV2=-CA2*GNUV
      FUW2=-CA2*GNUW
      FVW2=-CA2*GNVW
C
C---------- EAST WALL --------------------------------------------------
C
      RLE=ELE(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLE=MIN(1.,RLE)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLE*CA1W*TED(I,J,K)*1.5*((COSYE(I,J,K)**2)
     1          +(COSZE(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLE*CA1W*TED(I,J,K)*(
     1    COSYE(I,J,K)*COSZE(I,J,K)*(VV(I,J,K)+WW(I,J,K))
     3       +COSXE(I,J,K)*COSZE(I,J,K)*UV(I,J,K)
     4       +COSXE(I,J,K)*COSYE(I,J,K)*UW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLE*CA2W*(    
     1        +COSYE(I,J,K)*COSZE(I,J,K)*(FVV2+FWW2)
     2      +((COSYE(I,J,K)**2)+(COSZE(I,J,K)**2))*FVW2
     3         +COSXE(I,J,K)*COSZE(I,J,K)*FUV2
     4         +COSXE(I,J,K)*COSYE(I,J,K)*FUW2           )
C
C---------- WEST WALL --------------------------------------------------
C
      RLW=ELW(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLW=MIN(1.,RLW)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLW*CA1W*TED(I,J,K)*1.5*((COSYW(I,J,K)**2)
     1          +(COSZW(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLW*CA1W*TED(I,J,K)*(
     1    COSYW(I,J,K)*COSZW(I,J,K)*(VV(I,J,K)+WW(I,J,K))
     3       +COSXW(I,J,K)*COSZW(I,J,K)*UV(I,J,K)
     4       +COSXW(I,J,K)*COSYW(I,J,K)*UW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLW*CA2W*(    
     1        +COSYW(I,J,K)*COSZW(I,J,K)*(FVV2+FWW2)
     2      +((COSYW(I,J,K)**2)+(COSZW(I,J,K)**2))*FVW2
     3         +COSXW(I,J,K)*COSZW(I,J,K)*FUV2
     4         +COSXW(I,J,K)*COSYW(I,J,K)*FUW2           )
C
C---------- NORTH WALL -------------------------------------------------
C
      RLN=ELN(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLN=MIN(1.,RLN)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLN*CA1W*TED(I,J,K)*1.5*((COSYN(I,J,K)**2)
     1          +(COSZN(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLN*CA1W*TED(I,J,K)*(
     1    COSYN(I,J,K)*COSZN(I,J,K)*(VV(I,J,K)+WW(I,J,K))
     3       +COSXN(I,J,K)*COSZN(I,J,K)*UV(I,J,K)
     4       +COSXN(I,J,K)*COSYN(I,J,K)*UW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLN*CA2W*(    
     1        +COSYN(I,J,K)*COSZN(I,J,K)*(FVV2+FWW2)
     2      +((COSYN(I,J,K)**2)+(COSZN(I,J,K)**2))*FVW2
     3         +COSXN(I,J,K)*COSZN(I,J,K)*FUV2
     4         +COSXN(I,J,K)*COSYN(I,J,K)*FUW2           )
C
C---------- SOUTH WALL -------------------------------------------------
C
      RLS=ELS(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLS=MIN(1.,RLS)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLS*CA1W*TED(I,J,K)*1.5*((COSYS(I,J,K)**2)
     1          +(COSZS(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLS*CA1W*TED(I,J,K)*(
     1    COSYS(I,J,K)*COSZS(I,J,K)*(VV(I,J,K)+WW(I,J,K))
     3       +COSXS(I,J,K)*COSZS(I,J,K)*UV(I,J,K)
     4       +COSXS(I,J,K)*COSYS(I,J,K)*UW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLS*CA2W*(    
     1        +COSYS(I,J,K)*COSZS(I,J,K)*(FVV2+FWW2)
     2      +((COSYS(I,J,K)**2)+(COSZS(I,J,K)**2))*FVW2
     3         +COSXS(I,J,K)*COSZS(I,J,K)*FUV2
     4         +COSXS(I,J,K)*COSYS(I,J,K)*FUW2           )
C
C---------- TOP WALL -------------------------------------------------
C
      RLT=ELT(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLT=MIN(1.,RLT)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLT*CA1W*TED(I,J,K)*1.5*((COSYT(I,J,K)**2)
     1          +(COSZT(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLT*CA1W*TED(I,J,K)*(
     1    COSYT(I,J,K)*COSZT(I,J,K)*(VV(I,J,K)+WW(I,J,K))
     3       +COSXT(I,J,K)*COSZT(I,J,K)*UV(I,J,K)
     4       +COSXT(I,J,K)*COSYT(I,J,K)*UW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLT*CA2W*(    
     1        +COSYT(I,J,K)*COSZT(I,J,K)*(FVV2+FWW2)
     2      +((COSYT(I,J,K)**2)+(COSZT(I,J,K)**2))*FVW2
     3         +COSXT(I,J,K)*COSZT(I,J,K)*FUV2
     4         +COSXT(I,J,K)*COSYT(I,J,K)*FUW2           )
C
C---------- BOTTOM WALL ----------------------------------------------
C
      RLB=ELB(I,J,K)*(TKE(I,J,K)**1.5)/(2.55*TED(I,J,K)+TINY)
      RLB=MIN(1.,RLB)
C---------- FIJ1W ----------------------------------------------------      
      SP(I,J,K)=SP(I,J,K)-FW*RLB*CA1W*TED(I,J,K)*1.5*((COSYB(I,J,K)**2)
     1          +(COSZB(I,J,K)**2))/(TKE(I,J,K)+TINY)
C
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLB*CA1W*TED(I,J,K)*(
     1    COSYB(I,J,K)*COSZB(I,J,K)*(VV(I,J,K)+WW(I,J,K))
     3       +COSXB(I,J,K)*COSZB(I,J,K)*UV(I,J,K)
     4       +COSXB(I,J,K)*COSZB(I,J,K)*UW(I,J,K) )/
     5            (TKE(I,J,K)+TINY)
C--------- FIJ2W -----------------------------------------------------
      SU(I,J,K)=SU(I,J,K)-1.5*FW*RLB*CA2W*(    
     1        +COSYB(I,J,K)*COSZB(I,J,K)*(FVV2+FWW2)
     2      +((COSYB(I,J,K)**2)+(COSZB(I,J,K)**2))*FVW2
     3         +COSXB(I,J,K)*COSZB(I,J,K)*FUV2
     4         +COSXB(I,J,K)*COSYB(I,J,K)*FUW2           )
C
 50   CONTINUE
C
      CALL MODTVW
C
      RESORK=0.
      DO 30 I=2,NIM1
      DO 30 J=2,NJM1
      DO 30 K=2,NKM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*VW(I+1,J,K)+AW(I,J,K)*VW(I-1,J,K)
     1     +AN(I,J,K)*VW(I,J+1,K)+AS(I,J,K)*VW(I,J-1,K)
     1     +AT(I,J,K)*VW(I,J,K+1)+AB(I,J,K)*VW(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*VW(I,J,K)
C
      IF(SP(I,J,K).EQ.-GREAT) RESOR=0.
      RESORK=RESORK+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFK
      SU(I,J,K)=SU(I,J,K)+(1.-URFK)*AP(I,J,K)*VW(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPK
      CALL LISOLV(3,2,2,NIM1,NJ,NK,VW)
 40   CONTINUE
C
      RETURN
      END
C
      SUBROUTINE CALCED
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
!	NYAP
	DIMENSION ED(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
!
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      NJRSL=NJR+INWP
      NJRSL1=NJRSL+1
      NJRNL=NJRN-INWP+1
      NJRNL1=NJRNL-1
      NIEL=NIR+INWP
      NIEL1=NIEL+1
      NMWL=NIRMS-INWP+1
      NMWL1=NMWL-1
      NMEL=NIRMF+INWP
      NMEL=NI-1
      NMEL1=NMEL+1
      NIWL=NI-INWP-2
      NIWL1=NIWL-1
C
      NJED1=NJED+1
      NKED1=NKED+1 
      PR=PRTED
C
      DO 100 K=2,NKM1
      DO 100 J=2,NJM1
      DO 100 I=2,NIM1
      SMP=CX(I-1,J,K)-CX(I,J,K)+CY(I,J-1,K)-CY(I,J,K)+
     1    CZ(I,J,K-1)-CZ(I,J,K)
      SU(I,J,K)=ABS(SMP)*TED(I,J,K)
      SP(I,J,K)=-ABS(SMP)
  100 CONTINUE
C
C ------- CALCULATE ADVECTION COEFFICIENTS  ---------------------
C
      CALL ADVECT(1.,PRTED)
C
      CALL BQUICK(1,1.,TED)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(1.,PRTED,1,TED)
C
 
      DO 10 I=2,NIM1
      DO 10 J=2,NJM1
      DO 10 K=2,NKM1
C      SU(I,J,K)=0.
C
C      XLST=DEN(I,J,K)*XL*ABS(TKE(I,J,K))**0.5/VISCOS
C      FMU=(1.-EXP(-0.0198*XLST))*(1.+5.29/XLST)
      RT=DEN(I,J,K)*(TKE(I,J,K)**2)/(TED(I,J,K)+TINY)/VISCOS
      rtl(i,j,k)=rt      
      F1=1.
      F2=1.
C      F2=1.-0.3*EXP(-RT**2)
C      XLD=2.50*XL*XLST/(XLST+5.29)+TINY
C      DAMP=EXP(-0.00057*XLST**2)
C      PK2(I,J,K)=C2*F2*ABS(TKE(I,J,K))**1.5/C1/XLD*DAMP
C
C     SU(I,J,K)=SU(I,J,K)
C     1+F1*C1*CMU*DEN(I,J,K)**2*(PK1(I,J,K)    )*TKE(I,J,K)
C     1*VOL(I,J,K)/(VIS(I,J,K)-VISCOS+TINY)
C
C      SP(I,J,K)=-F2*C2*CMU*DEN(I,J,K)**2*TKE(I,J,K)*VOL(I,J,K)
C      1/(VIS(I,J,K)-VISCOS+TINY)
C*	IF (INLEVM.EQ.0) C2L=C2*(1.-FLOAT(INWM-1)*0.3*EXP(-RT*RT))
C*	IF (INLEVM.EQ.1) C2L=C2
	C2L=C2*(1.-FLOAT(INWM-1)*0.3*EXP(-RT*RT))
      SP1=DEN(I,J,K)*VOL(I,J,K)*(
     1          C1*PK1(I,J,K)-C2L*TED(I,J,K))/(TKE(I,J,K)+TINY)
      SU(I,J,K)=SU(I,J,K)+MAX(SP1,0.)*TED(I,J,K)
      SP(I,J,K)=SP(I,J,K)-MAX(-SP1,0.)
C
C-------------- LOW-RE E-TERM --------------------------------
C
	SU(I,J,K)=SU(I,J,K)+FLOAT(INWM-1)*VOL(I,J,K)*DEN(I,J,K)
     1                    *ETRM(I,J,K)
C
C    YAP CORRECTION
C
      IF (IYAP.EQ.0) THEN
	XLTUR=ABS(TKE(I,J,K))**1.5/(TED(I,J,K)+TINY)
      XLWAL=2.55*XL3(I,J,K)
      VIST=VIS(I,J,K)-VISCOS
      SU(I,J,K)=SU(I,J,K)+
     1DEN(I,J,K)*AMAX1(0.83*(XLTUR/XLWAL-1.),0.)*((XLTUR/XLWAL)**2)
     1*TED(I,J,K)*TED(I,J,K)*VOL(I,J,K)/(TKE(I,J,K)+TINY)

	ELSE
!	NYAP CORRECTION TERM
	CL=2.55
	BE=0.1069
	SIJ2=(2.0*DUDXC(I,J,K))**2+(2.0*DVDYC(I,J,K))**2+
	1(2.0*DWDZC(I,J,K))**2+2.0*(DUDYC(I,J,K)+DVDXC(I,J,K))**2+
     12.0*(DUDZC(I,J,K)+DWDXC(I,J,K))**2+
     12.0*(DVDZC(I,J,K)+DWDYC(I,J,K))**2
	OMGIJ2=2.0*(DUDYC(I,J,K)-DVDXC(I,J,K))**2+
	12.0*(DUDZC(I,J,K)-DWDXC(I,J,K))**2+
     12.0*(DVDZC(I,J,K)-DWDYC(I,J,K))**2

      ED(I,J,K)=TED(I,J,K)+2.0*VISCOS*(DKDX(I,J,K)*DKDX(I,J,K)+
     1DKDY(I,J,K)*DKDY(I,J,K)+DKDZ(I,J,K)*DKDZ(I,J,K))/DEN(I,J,K)
	ST=AMAX1((TKE(I,J,K)/(TED(I,J,K)+TINY)),SQRT(VISCOS/ED(I,J,K)))*
	1SQRT(0.5*SIJ2)
	OMGT=AMAX1((TKE(I,J,K)/(TED(I,J,K)+TINY)),SQRT(VISCOS/ED(I,J,K)))*
	1SQRT(0.5*OMGIJ2)
	ETAP=AMAX1(ST,OMGT)
!	COMG=0.83*AMIN1(1.0,RT/5.0)/(0.8+0.7*((ETAP/3.33)**4)*EXP(-RT/12.5))


!******************	NLEVM3 REFINEMENT

      FS=(AMIN1(AMAX1(((AMIN1(1.0,RT/50.0))**2)*PK1(I,J,K)/(TED(I,J,K)
     1+TINY),0.0)/0.75,1.0))**2

	COMG=(1.85*(1-FS)*(AMIN1(1.0,RT/25.0))**2+0.75*
	1AMIN1(1.0,RT/40.0)*FS)/(0.8+0.7*((ETAP/3.33)**4)*EXP(-RT/12.5))

!*****************

	DLEDY=CL*(1-EXP(-BE*RT))+BE*CL*RT*EXP(-BE*RT)
	F=(SQRT(DLDXJ2(I,J,K))-DLEDY)/CL
	
      SU(I,J,K)=SU(I,J,K)+
     1DEN(I,J,K)*AMAX1(COMG*F*(F+1)**2,0.0)*TED(I,J,K)*TED(I,J,K)*
     1VOL(I,J,K)/(TKE(I,J,K)+TINY)
!
	ENDIF
C
 10   CONTINUE
C
C      RETURN
      CALL MODTED
C
      RESORE=0.
      DO 30 K=KED,NKED
      DO 30 J=JED,NJED
      DO 30 I=2,NIM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*TED(I+1,J,K)+AW(I,J,K)*TED(I-1,J,K)
     1     +AN(I,J,K)*TED(I,J+1,K)+AS(I,J,K)*TED(I,J-1,K)
     1     +AT(I,J,K)*TED(I,J,K+1)+AB(I,J,K)*TED(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*TED(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
      RESORE=RESORE+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFE
      SU(I,J,K)=SU(I,J,K)+(1.-URFE)*AP(I,J,K)*TED(I,J,K)
 30   CONTINUE
C
C      IF(INWM.EQ.2)THEN
      DO 40 N=1,NSWPE
      CALL LISOLV(3,JED,KED,NIM1,NJED1,NKED1,TED)
 40   CONTINUE
C      ENDIF
C
C      IF(INWM.EQ.1)THEN
C      DO 41 N=1,NSWPE
C      CALL LISOLV(3,NJRSL,KED,NIM1,NJRNL+1,NKED1,TED)
C      CALL LISOLV(NIEL,JED,KED,NMWL+1,NJED1,NKED1,TED)
C      CALL LISOLV(NMEL,JED,KED,NIWL+1,NJED1,NKED1,TED)
C   41 CONTINUE
C      ENDIF
C
      RETURN
      END
C
      SUBROUTINE CALCP
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C 
C*       DO 4 K=1,NK
C*       DO 4 J=1,NJ
C
C      DU1(1,J,K)=0.*DU1(2,J,K)
C      DU2(1,J,K)=0.*DU2(2,J,K)
C      DU3(1,J,K)=0.*DU3(2,J,K)
C      DU1(NI,J,K)=0.*DU1(NIM1,J,K)
C      DU2(NI,J,K)=0.*DU2(NIM1,J,K)
C      DU3(NI,J,K)=0.*DU3(NIM1,J,K)
C
C      DV1(1,J,K)=DV1(2,J,K)
C      DV2(1,J,K)=DV2(2,J,K)
C      DV3(1,J,K)=DV3(2,J,K)
C      DV1(NI,J,K)=DV1(NIM1,J,K)
C      DV2(NI,J,K)=DV2(NIM1,J,K)
C      DV3(NI,J,K)=DV3(NIM1,J,K)
C
C      DW1(1,J,K)=DW1(2,J,K)
C      DW2(1,J,K)=DW2(2,J,K)
C      DW3(1,J,K)=DW3(2,J,K)
C      DW1(NI,J,K)=DW1(NIM1,J,K)
C      DW2(NI,J,K)=DW2(NIM1,J,K)
C      DW3(NI,J,K)=DW3(NIM1,J,K)
C
C*    4 CONTINUE
C
C      DO 105 K=2,NKM1
C      P(NIR,NJR,K)=2.*P(NIR+1,NJR,K)-P(NIR+2,NJR,K)
C      P(NI-2,NJR,K)=2.*P(NI-3,NJR,K)-P(NI-4,NJR,K)
C      P(NIR,NJRN+1,K)=2.*P(NIR+1,NJRN+1,K)-P(NIR+2,NJRN+1,K)
C      P(NI-2,NJRN+1,K)=2.*P(NI-3,NJRN+1,K)-P(NI-4,NJRN+1,K)
C  105 CONTINUE
C
      DO 5 I=2,NI-2
      DO 5 J=2,NJ-1
      DO 5 K=2,NK-1
C
      PRWP=(1.-FX(I-1,J,K))*P(I-1,J,K)+FX(I-1,J,K)*P(I,J,K)
      PREP=(1.-FX(I,J,K))*P(I,J,K)+FX(I,J,K)*P(I+1,J,K)
      PRWE=PREP
      PREE=(1.-FX(I+1,J,K))*P(I+1,J,K)+FX(I+1,J,K)*P(I+2,J,K)
C
      UE=
     1 (1.-FX(I,J,K))*(U(I,J,K)-DU1(I,J,K)*(PRWP-PREP))
     1 +FX(I,J,K)*(U(I+1,J,K)-DU1(I+1,J,K)*(PRWE-PREE))
     1 +((1.-FX(I,J,K))*DU1(I,J,K)+FX(I,J,K)*DU1(I+1,J,K))
     1 *(P(I,J,K)-P(I+1,J,K))

C
      VE=
     1 (1.-FX(I,J,K))*(V(I,J,K)-DV1(I,J,K)*(PRWP-PREP))
     1 +FX(I,J,K)*(V(I+1,J,K)-DV1(I+1,J,K)*(PRWE-PREE))
     1 +((1.-FX(I,J,K))*DV1(I,J,K)+FX(I,J,K)*DV1(I+1,J,K))
     1 *(P(I,J,K)-P(I+1,J,K))
C
      WE=
     1 (1.-FX(I,J,K))*(W(I,J,K)-DW1(I,J,K)*(PRWP-PREP))
     1 +FX(I,J,K)*(W(I+1,J,K)-DW1(I+1,J,K)*(PRWE-PREE))
     1 +((1.-FX(I,J,K))*DW1(I,J,K)+FX(I,J,K)*DW1(I+1,J,K))
     1 *(P(I,J,K)-P(I+1,J,K))
C
      FLOW=UE*((Y(I,J-1,K)-Y(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K))
     1        -(Z(I,J-1,K)-Z(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K)))
     1    +VE*((Z(I,J-1,K)-Z(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K))
     1        -(X(I,J-1,K)-X(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K)))
     1    +WE*((X(I,J-1,K)-X(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K))
     1        -(Y(I,J-1,K)-Y(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K)))
      CX(I,J,K)=0.5*FLOW
C
 5    CONTINUE
C
C      DO 106 K=2,NKM1
C      P(NIR,NJR,K)=P(NIR,NJR+1,K)
C      P(NI-2,NJR,K)=P(NI-2,NJR+1,K)
C      P(NIR,NJRN+1,K)=P(NIR,NJRN,K)
C      P(NI-2,NJRN+1,K)=P(NI-2,NJRN,K)
  106 CONTINUE
C
      DO 6 I=2,NI-1
      DO 6 J=2,NJ-2
      DO 6 K=2,NK-1
C
      PRSP=(1.-FY(I,J-1,K))*P(I,J-1,K)+FY(I,J-1,K)*P(I,J,K)
      PRNP=(1.-FY(I,J,K))*P(I,J,K)+FY(I,J,K)*P(I,J+1,K)
      PRSN=PRNP
      PRNN=(1.-FY(I,J+1,K))*P(I,J+1,K)+FY(I,J+1,K)*P(I,J+2,K)
C
      UN=
     1 (1.-FY(I,J,K))*(U(I,J,K)-DU2(I,J,K)*(PRSP-PRNP))
     1 +FY(I,J,K)*(U(I,J+1,K)-DU2(I,J+1,K)*(PRSN-PRNN))
     1 +((1.-FY(I,J,K))*DU2(I,J,K)+FY(I,J,K)*DU2(I,J+1,K))
     1 *(P(I,J,K)-P(I,J+1,K))

C
      VN=
     1 (1.-FY(I,J,K))*(V(I,J,K)-DV2(I,J,K)*(PRSP-PRNP))
     1 +FY(I,J,K)*(V(I,J+1,K)-DV2(I,J+1,K)*(PRSN-PRNN))
     1 +((1.-FY(I,J,K))*DV2(I,J,K)+FY(I,J,K)*DV2(I,J+1,K))
     1 *(P(I,J,K)-P(I,J+1,K))
C
      WN=
     1 (1.-FY(I,J,K))*(W(I,J,K)-DW2(I,J,K)*(PRSP-PRNP))
     1 +FY(I,J,K)*(W(I,J+1,K)-DW2(I,J+1,K)*(PRSN-PRNN))
     1 +((1.-FY(I,J,K))*DW2(I,J,K)+FY(I,J,K)*DW2(I,J+1,K))
     1 *(P(I,J,K)-P(I,J+1,K))
C
      FLOW=UN*((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1        -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))
     1    +VN*((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1        -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))
     1    +WN*((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1        -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))
      CY(I,J,K)=0.5*FLOW
C
 6    CONTINUE
C
      DO 7 I=2,NI-1
      DO 7 J=2,NJ-1
      DO 7 K=2,NK-2
C
      PRBP=(1.-FZ(I,J,K-1))*P(I,J,K-1)+FZ(I,J,K-1)*P(I,J,K)
      PRTP=(1.-FZ(I,J,K))*P(I,J,K)+FZ(I,J,K)*P(I,J,K+1)
      PRBT=PRTP
      PRTT=(1.-FZ(I,J,K+1))*P(I,J,K+1)+FZ(I,J,K+1)*P(I,J,K+2)
C
      UT=
     1 (1.-FZ(I,J,K))*(U(I,J,K)-DU3(I,J,K)*(PRBP-PRTP))
     1 +FZ(I,J,K)*(U(I,J,K+1)-DU3(I,J,K+1)*(PRBT-PRTT))
     1 +((1.-FZ(I,J,K))*DU3(I,J,K)+FZ(I,J,K)*DU3(I,J,K+1))
     1 *(P(I,J,K)-P(I,J,K+1))
C
      VT=
     1 (1.-FZ(I,J,K))*(V(I,J,K)-DV3(I,J,K)*(PRBP-PRTP))
     1 +FZ(I,J,K)*(V(I,J,K+1)-DV3(I,J,K+1)*(PRBT-PRTT))
     1 +((1.-FZ(I,J,K))*DV3(I,J,K)+FZ(I,J,K)*DV3(I,J,K+1))
     1 *(P(I,J,K)-P(I,J,K+1))
C
      WT=
     1 (1.-FZ(I,J,K))*(W(I,J,K)-DW3(I,J,K)*(PRBP-PRTP))
     1 +FZ(I,J,K)*(W(I,J,K+1)-DW3(I,J,K+1)*(PRBT-PRTT))
     1 +((1.-FZ(I,J,K))*DW3(I,J,K)+FZ(I,J,K)*DW3(I,J,K+1))
     1 *(P(I,J,K)-P(I,J,K+1))
C
      FLOW=UT*((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1        -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))
     1    +VT*((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1        -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))
     1    +WT*((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1        -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))
      CZ(I,J,K)=0.5*FLOW
C
 7    CONTINUE
C
      I=1
      DO 8 J=2,NJ-1
      DO 8 K=2,NK-1
      FLOW=U(I,J,K)*((Y(I,J,K)-Y(I,J-1,K-1))*(Z(I,J-1,K)-Z(I,J,K-1))
     1              -(Z(I,J,K)-Z(I,J-1,K-1))*(Y(I,J-1,K)-Y(I,J,K-1)))
     1    +V(I,J,K)*((Z(I,J,K)-Z(I,J-1,K-1))*(X(I,J-1,K)-X(I,J,K-1))
     1              -(X(I,J,K)-X(I,J-1,K-1))*(Z(I,J-1,K)-Z(I,J,K-1)))
     1    +W(I,J,K)*((X(I,J,K)-X(I,J-1,K-1))*(Y(I,J-1,K)-Y(I,J,K-1))
     1              -(Y(I,J,K)-Y(I,J-1,K-1))*(X(I,J-1,K)-X(I,J,K-1)))
      CX(I,J,K)=0.5*FLOW
C      CX(I,J,K)=CX(NI-3,J,K)
       CX(I+1,J,K)=CX(NI-2,J,K)
C      CX(I-1,J,K)=CX(NI-3,J,K)
C
 8    CONTINUE
C
      I=NI-1
      DO 9 J=2,NJ-1
      DO 9 K=2,NK-1
      FLOW=U(I,J,K)*((Y(I,J,K)-Y(I,J-1,K-1))*(Z(I,J-1,K)-Z(I,J,K-1))
     1                -(Z(I,J,K)-Z(I,J-1,K-1))*(Y(I,J-1,K)-Y(I,J,K-1)))
     1  +V(I,J,K)*((Z(I,J,K)-Z(I,J-1,K-1))*(X(I,J-1,K)-X(I,J,K-1))
     1                -(X(I,J,K)-X(I,J-1,K-1))*(Z(I,J-1,K)-Z(I,J,K-1)))
     1   +W(I,J,K)*((X(I,J,K)-X(I,J-1,K-1))*(Y(I,J-1,K)-Y(I,J,K-1))
     1                -(Y(I,J,K)-Y(I,J-1,K-1))*(X(I,J-1,K)-X(I,J,K-1)))
      CX(I,J,K)=0.5*FLOW
      CX(3,J,K)=CX(I,J,K)
C      CX(I,J,K)=CX(3,J,K)
C      CX(1,J,K)=CX(I,J,K)
C
 9    CONTINUE
C
      J=1
      DO 10 K=2,NK-1
      DO 10 I=2,NI-1
      FLOW=U(I,J,K)*((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1              -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))
     1    +V(I,J,K)*((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1              -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))
     1    +W(I,J,K)*((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1              -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))
      CY(I,J,K)=0.*0.5*FLOW
C
 10   CONTINUE
C
      J=NJ-1
      DO 11 K=2,NK-1
      DO 11 I=2,NI-1
      FLOW=U(I,J+1,K)*((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1                -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))
     1    +V(I,J+1,K)*((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1                -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))
     1    +W(I,J+1,K)*((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1                -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))
      CY(I,J,K)=0.*0.5*FLOW
C
 11   CONTINUE
C
      K=1
      DO 12 I=2,NI-1
      DO 12 J=2,NJ-1
      FLOW=U(I,J,K)*((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1              -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))
     1    +V(I,J,K)*((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1              -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))
     1    +W(I,J,K)*((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1              -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))
      CZ(I,J,K)=0.5*FLOW
C
 12   CONTINUE
C
      K=NK-1
      DO 13 I=2,NI-1
      DO 13 J=2,NJ-1
      FLOW=U(I,J,K+1)*((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1                -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))
     1    +V(I,J,K+1)*((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1                -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))
     1    +W(I,J,K+1)*((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1                -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))
      CZ(I,J,K)=0.5*FLOW
C
 13   CONTINUE
C---------- RIBS ----------------------------------------------------
C      GO TO 96
C      DO 91 K=2,NKM1
C
C      DO 92 J=2,NJR
C      CX(NIR,J,K)=0.
C      CX(NIRMS,J,K)=0.
C      CX(NIRMF,J,K)=0.
C      CX(NI-3,J,K)=0.
C   92 CONTINUE
C      DO 93 J=NJRN+1,NJM1
C      CX(NIR,J,K)=0.
C      CX(NIRMS,J,K)=0.
C      CX(NIRMF,J,K)=0.
C      CX(NI-3,J,K)=0.
C   93 CONTINUE
C
C      DO 94 I=2,NIR
C      CY(I,NJR,K)=0.
C      CY(I,NJRN,K)=0.
C   94 CONTINUE
C      DO 95 I=NIRMS+1,NIRMF
C      CY(I,NJR,K)=0.
C      CY(I,NJRN,K)=0.
C   95 CONTINUE
C      DO 95 I=NI-2, NIM1
C      CY(I,NJR,K)=0.
C      CY(I,NJRN,K)=0.
C   95 CONTINUE
C   91 CONTINUE
C   96 CONTINUE
C---------------------------------------------------------------
      RESORM=0.
      DO 15 I=3,NI-1
      DO 15 J=2,NJ-1
      DO 15 K=2,NK-1
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
      SMP= DENE*CX(I,J,K)-DENW*CX(I-1,J,K)
     1    +DENN*CY(I,J,K)-DENS*CY(I,J-1,K)
     1    +DENT*CZ(I,J,K)-DENB*CZ(I,J,K-1)
C
      SU(I,J,K)=-SMP
      SP(I,J,K)=0.
      RESORM=RESORM+ABS(SMP)
 15   CONTINUE
C
      DO 16 I=2,NI-1
      DO 16 J=2,NJ-1
      DO 16 K=2,NK-1
      DUE=(1.-FX(I,J,K))*DU1(I,J,K)+FX(I,J,K)*DU1(I+1,J,K)
      DVE=(1.-FX(I,J,K))*DV1(I,J,K)+FX(I,J,K)*DV1(I+1,J,K)
      DWE=(1.-FX(I,J,K))*DW1(I,J,K)+FX(I,J,K)*DW1(I+1,J,K)
C
      DUN=(1.-FY(I,J,K))*DU2(I,J,K)+FY(I,J,K)*DU2(I,J+1,K)
      DVN=(1.-FY(I,J,K))*DV2(I,J,K)+FY(I,J,K)*DV2(I,J+1,K)
      DWN=(1.-FY(I,J,K))*DW2(I,J,K)+FY(I,J,K)*DW2(I,J+1,K)
C
      DUT=(1.-FZ(I,J,K))*DU3(I,J,K)+FZ(I,J,K)*DU3(I,J,K+1)
      DVT=(1.-FZ(I,J,K))*DV3(I,J,K)+FZ(I,J,K)*DV3(I,J,K+1)
      DWT=(1.-FZ(I,J,K))*DW3(I,J,K)+FZ(I,J,K)*DW3(I,J,K+1)
C
      FLOW=DUE*((Y(I,J-1,K)-Y(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K))
     1         -(Z(I,J-1,K)-Z(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K)))
     1    +DVE*((Z(I,J-1,K)-Z(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K))
     1         -(X(I,J-1,K)-X(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K)))
     1    +DWE*((X(I,J-1,K)-X(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K))
     1         -(Y(I,J-1,K)-Y(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K)))
      CXM(I,J,K)=0.5*FLOW
C
      FLOW=DUN*((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1         -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))
     1    +DVN*((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1         -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))
     1    +DWN*((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1         -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))
      CYM(I,J,K)=0.5*FLOW
C
      FLOW=DUT*((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1         -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))
     1    +DVT*((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1         -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))
     1    +DWT*((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1         -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))
      CZM(I,J,K)=0.5*FLOW
 16   CONTINUE
C
      DO 20 I=2,NI-1
      DO 20 J=2,NJ-1
      DO 20 K=2,NK-1
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
      FLOWE=CXM(I  ,J,K)
      FLOWW=CXM(I-1,J,K)
C
      AE(I,J,K)=DENE*FLOWE
      AW(I,J,K)=DENW*FLOWW
C
      FLOWN=CYM(I,J  ,K)
      FLOWS=CYM(I,J-1,K)
C
      AN(I,J,K)=DENN*FLOWN
      AS(I,J,K)=DENS*FLOWS
C
      FLOWT=CZM(I,J,K  )
      FLOWB=CZM(I,J,K-1)
C
      AT(I,J,K)=DENT*FLOWT
      AB(I,J,K)=DENB*FLOWB
C
 20   CONTINUE
C
      CALL MODP
C
      DO 30 I=3,NI-1
      DO 30 J=2,NJ-1
      DO 30 K=2,NK-1
      PP(I,J,K)=0.
      SUPP(I,J,K)=0.
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+
     1          AT(I,J,K)+AB(I,J,K)-SP(I,J,K)
 30   CONTINUE
      DO  I=1,NI
      DO  J=1,NJ
      DO  K=1,NK
      PP(I,J,K)=0.
      SUPP(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
*      DO 40 N=1,NSWPP
      go to 401
c
      DO  I=2,NI-1
      DO  J=2,NJ-1
      DO  K=2,NK-1
C
c      SU(I,J,K)=SU(I,J,K)-SUPP(I,J,K)
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
      DUE2=(1.-FX(I,J,K))*DU2(I,J,K)+FX(I,J,K)*DU2(I+1,J,K)
      DVE2=(1.-FX(I,J,K))*DV2(I,J,K)+FX(I,J,K)*DV2(I+1,J,K)
      DWE2=(1.-FX(I,J,K))*DW2(I,J,K)+FX(I,J,K)*DW2(I+1,J,K)
C
      DUE3=(1.-FX(I,J,K))*DU3(I,J,K)+FX(I,J,K)*DU3(I+1,J,K)
      DVE3=(1.-FX(I,J,K))*DV3(I,J,K)+FX(I,J,K)*DV3(I+1,J,K)
      DWE3=(1.-FX(I,J,K))*DW3(I,J,K)+FX(I,J,K)*DW3(I+1,J,K)
C
      DUW2=(1.-FX(I-1,J,K))*DU2(I-1,J,K)+FX(I-1,J,K)*DU2(I,J,K)
      DVW2=(1.-FX(I-1,J,K))*DV2(I-1,J,K)+FX(I-1,J,K)*DV2(I,J,K)
      DWW2=(1.-FX(I-1,J,K))*DW2(I-1,J,K)+FX(I-1,J,K)*DW2(I,J,K)
C
      DUW3=(1.-FX(I-1,J,K))*DU3(I-1,J,K)+FX(I-1,J,K)*DU3(I,J,K)
      DVW3=(1.-FX(I-1,J,K))*DV3(I-1,J,K)+FX(I-1,J,K)*DV3(I,J,K)
      DWW3=(1.-FX(I-1,J,K))*DW3(I-1,J,K)+FX(I-1,J,K)*DW3(I,J,K)
C
      DUN1=(1.-FY(I,J,K))*DU1(I,J,K)+FY(I,J,K)*DU1(I,J+1,K)
      DVN1=(1.-FY(I,J,K))*DV1(I,J,K)+FY(I,J,K)*DV1(I,J+1,K)
      DWN1=(1.-FY(I,J,K))*DW1(I,J,K)+FY(I,J,K)*DW1(I,J+1,K)
C
      DUN3=(1.-FY(I,J,K))*DU3(I,J,K)+FY(I,J,K)*DU3(I,J+1,K)
      DVN3=(1.-FY(I,J,K))*DV3(I,J,K)+FY(I,J,K)*DV3(I,J+1,K)
      DWN3=(1.-FY(I,J,K))*DW3(I,J,K)+FY(I,J,K)*DW3(I,J+1,K)
C
      DUS1=(1.-FY(I,J-1,K))*DU1(I,J-1,K)+FY(I,J-1,K)*DU1(I,J,K)
      DVS1=(1.-FY(I,J-1,K))*DV1(I,J-1,K)+FY(I,J-1,K)*DV1(I,J,K)
      DWS1=(1.-FY(I,J-1,K))*DW1(I,J-1,K)+FY(I,J-1,K)*DW1(I,J,K)
C
      DUS3=(1.-FY(I,J-1,K))*DU3(I,J-1,K)+FY(I,J-1,K)*DU3(I,J,K)
      DVS3=(1.-FY(I,J-1,K))*DV3(I,J-1,K)+FY(I,J-1,K)*DV3(I,J,K)
      DWS3=(1.-FY(I,J-1,K))*DW3(I,J-1,K)+FY(I,J-1,K)*DW3(I,J,K)
C
      DUT1=(1.-FZ(I,J,K))*DU1(I,J,K)+FZ(I,J,K)*DU1(I,J,K+1)
      DVT1=(1.-FZ(I,J,K))*DV1(I,J,K)+FZ(I,J,K)*DV1(I,J,K+1)
      DWT1=(1.-FZ(I,J,K))*DW1(I,J,K)+FZ(I,J,K)*DW1(I,J,K+1)
C
      DUT2=(1.-FZ(I,J,K))*DU2(I,J,K)+FZ(I,J,K)*DU2(I,J,K+1)
      DVT2=(1.-FZ(I,J,K))*DV2(I,J,K)+FZ(I,J,K)*DV2(I,J,K+1)
      DWT2=(1.-FZ(I,J,K))*DW2(I,J,K)+FZ(I,J,K)*DW2(I,J,K+1)
C
      DUB1=(1.-FZ(I,J,K-1))*DU1(I,J,K-1)+FZ(I,J,K-1)*DU1(I,J,K)
      DVB1=(1.-FZ(I,J,K-1))*DV1(I,J,K-1)+FZ(I,J,K-1)*DV1(I,J,K)
      DWB1=(1.-FZ(I,J,K-1))*DW1(I,J,K-1)+FZ(I,J,K-1)*DW1(I,J,K)
C
      DUB2=(1.-FZ(I,J,K-1))*DU2(I,J,K-1)+FZ(I,J,K-1)*DU2(I,J,K)
      DVB2=(1.-FZ(I,J,K-1))*DV2(I,J,K-1)+FZ(I,J,K-1)*DV2(I,J,K)
      DWB2=(1.-FZ(I,J,K-1))*DW2(I,J,K-1)+FZ(I,J,K-1)*DW2(I,J,K)
C
      FLOW2=DUE2*((Y(I,J-1,K)-Y(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K))
     1         -(Z(I,J-1,K)-Z(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K)))
     1    +DVE2*((Z(I,J-1,K)-Z(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K))
     1         -(X(I,J-1,K)-X(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K)))
     1    +DWE2*((X(I,J-1,K)-X(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K))
     1         -(Y(I,J-1,K)-Y(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K)))
C
      FLOW3=DUE3*((Y(I,J-1,K)-Y(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K))
     1         -(Z(I,J-1,K)-Z(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K)))
     1    +DVE3*((Z(I,J-1,K)-Z(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K))
     1         -(X(I,J-1,K)-X(I,J,K-1))*(Z(I,J-1,K-1)-Z(I,J,K)))
     1    +DWE3*((X(I,J-1,K)-X(I,J,K-1))*(Y(I,J-1,K-1)-Y(I,J,K))
     1         -(Y(I,J-1,K)-Y(I,J,K-1))*(X(I,J-1,K-1)-X(I,J,K)))
C
      AE2=0.5*FLOW2*DENE
      AE3=0.5*FLOW3*DENE
C
      FLOW2=
     1  DUW2*((Y(I-1,J-1,K)-Y(I-1,J,K-1))*(Z(I-1,J-1,K-1)-Z(I-1,J,K))
     1      -(Z(I-1,J-1,K)-Z(I-1,J,K-1))*(Y(I-1,J-1,K-1)-Y(I-1,J,K)))
     1 +DVW2*((Z(I-1,J-1,K)-Z(I-1,J,K-1))*(X(I-1,J-1,K-1)-X(I-1,J,K))
     1      -(X(I-1,J-1,K)-X(I-1,J,K-1))*(Z(I-1,J-1,K-1)-Z(I-1,J,K)))
     1 +DWW2*((X(I-1,J-1,K)-X(I-1,J,K-1))*(Y(I-1,J-1,K-1)-Y(I-1,J,K))
     1      -(Y(I-1,J-1,K)-Y(I-1,J,K-1))*(X(I-1,J-1,K-1)-X(I-1,J,K)))
C
      FLOW3=
     1  DUW3*((Y(I-1,J-1,K)-Y(I-1,J,K-1))*(Z(I-1,J-1,K-1)-Z(I-1,J,K))
     1      -(Z(I-1,J-1,K)-Z(I-1,J,K-1))*(Y(I-1,J-1,K-1)-Y(I-1,J,K)))
     1 +DVW3*((Z(I-1,J-1,K)-Z(I-1,J,K-1))*(X(I-1,J-1,K-1)-X(I-1,J,K))
     1      -(X(I-1,J-1,K)-X(I-1,J,K-1))*(Z(I-1,J-1,K-1)-Z(I-1,J,K)))
     1 +DWW3*((X(I-1,J-1,K)-X(I-1,J,K-1))*(Y(I-1,J-1,K-1)-Y(I-1,J,K))
     1      -(Y(I-1,J-1,K)-Y(I-1,J,K-1))*(X(I-1,J-1,K-1)-X(I-1,J,K)))
C
      AW2=0.5*FLOW2*DENW
      AW3=0.5*FLOW3*DENW
C
      FLOW1=DUN1*((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1         -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))
     1    +DVN1*((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1         -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))
     1    +DWN1*((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1         -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))
C
      FLOW3=DUN3*((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1         -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))
     1    +DVN3*((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1         -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))
     1    +DWN3*((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1         -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))
      AN1=0.5*FLOW1*DENN
      AN3=0.5*FLOW3*DENN
C
      FLOW1=
     1  DUS1*((Y(I,J-1,K)-Y(I-1,J-1,K-1))*(Z(I,J-1,K-1)-Z(I-1,J-1,K))
     1      -(Z(I,J-1,K)-Z(I-1,J-1,K-1))*(Y(I,J-1,K-1)-Y(I-1,J-1,K)))
     1 +DVS1*((Z(I,J-1,K)-Z(I-1,J-1,K-1))*(X(I,J-1,K-1)-X(I-1,J-1,K))
     1      -(X(I,J-1,K)-X(I-1,J-1,K-1))*(Z(I,J-1,K-1)-Z(I-1,J-1,K)))
     1 +DWS1*((X(I,J-1,K)-X(I-1,J-1,K-1))*(Y(I,J-1,K-1)-Y(I-1,J-1,K))
     1      -(Y(I,J-1,K)-Y(I-1,J-1,K-1))*(X(I,J-1,K-1)-X(I-1,J-1,K)))
C
      FLOW3=
     1  DUS3*((Y(I,J-1,K)-Y(I-1,J-1,K-1))*(Z(I,J-1,K-1)-Z(I-1,J-1,K))
     1      -(Z(I,J-1,K)-Z(I-1,J-1,K-1))*(Y(I,J-1,K-1)-Y(I-1,J-1,K)))
     1 +DVS3*((Z(I,J-1,K)-Z(I-1,J-1,K-1))*(X(I,J-1,K-1)-X(I-1,J-1,K))
     1      -(X(I,J-1,K)-X(I-1,J-1,K-1))*(Z(I,J-1,K-1)-Z(I-1,J-1,K)))
     1 +DWS3*((X(I,J-1,K)-X(I-1,J-1,K-1))*(Y(I,J-1,K-1)-Y(I-1,J-1,K))
     1      -(Y(I,J-1,K)-Y(I-1,J-1,K-1))*(X(I,J-1,K-1)-X(I-1,J-1,K)))
C
      AS1=0.5*FLOW1*DENS
      AS3=0.5*FLOW3*DENS
C
      FLOW1=DUT1*((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1         -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))
     1    +DVT1*((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1         -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))
     1    +DWT1*((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1         -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))
C
      FLOW2=DUT2*((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1         -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))
     1    +DVT2*((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1         -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))
     1    +DWT2*((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1         -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))
C
      AT1=0.5*FLOW1*DENT
      AT2=0.5*FLOW2*DENT
C
      FLOW1=
     1  DUB1*((Y(I,J,K-1)-Y(I-1,J-1,K-1))*(Z(I-1,J,K-1)-Z(I,J-1,K-1))
     1      -(Z(I,J,K-1)-Z(I-1,J-1,K-1))*(Y(I-1,J,K-1)-Y(I,J-1,K-1)))
     1 +DVB1*((Z(I,J,K-1)-Z(I-1,J-1,K-1))*(X(I-1,J,K-1)-X(I,J-1,K-1))
     1      -(X(I,J,K-1)-X(I-1,J-1,K-1))*(Z(I-1,J,K-1)-Z(I,J-1,K-1)))
     1 +DWB1*((X(I,J,K-1)-X(I-1,J-1,K-1))*(Y(I-1,J,K-1)-Y(I,J-1,K-1))
     1      -(Y(I,J,K-1)-Y(I-1,J-1,K-1))*(X(I-1,J,K-1)-X(I,J-1,K-1)))
C
      FLOW2=
     1  DUB2*((Y(I,J,K-1)-Y(I-1,J-1,K-1))*(Z(I-1,J,K-1)-Z(I,J-1,K-1))
     1      -(Z(I,J,K-1)-Z(I-1,J-1,K-1))*(Y(I-1,J,K-1)-Y(I,J-1,K-1)))
     1 +DVB2*((Z(I,J,K-1)-Z(I-1,J-1,K-1))*(X(I-1,J,K-1)-X(I,J-1,K-1))
     1      -(X(I,J,K-1)-X(I-1,J-1,K-1))*(Z(I-1,J,K-1)-Z(I,J-1,K-1)))
     1 +DWB2*((X(I,J,K-1)-X(I-1,J-1,K-1))*(Y(I-1,J,K-1)-Y(I,J-1,K-1))
     1      -(Y(I,J,K-1)-Y(I-1,J-1,K-1))*(X(I-1,J,K-1)-X(I,J-1,K-1)))
C
      AB1=0.5*FLOW1*DENB
      AB2=0.5*FLOW2*DENB
C 
      PPNE=(1.-FY(I,J,K))*
     1         ((1.-FX(I,J,K))*PP(I,J,K)+FX(I,J,K)*PP(I+1,J,K)) +
     1   FY(I,J,K)*
     1     ((1-FX(I,J+1,K))*PP(I,J+1,K)+FX(I,J+1,K)*PP(I+1,J+1,K))
C 
      PPSE=(1.-FY(I,J-1,K))*
cc     1         ((1.-FX(I,J-1,K))*PP(I,J-1,K)+FX(I,J-1,K)*PP(I,J,K)) +
     1         ((1.-FX(I,J-1,K))*PP(I,J-1,K)+
     1               FX(I,J-1,K)*PP(I+1,J-1,K)) +
     1   FY(I,J-1,K)*
     1     ((1-FX(I,J,K))*PP(I,J,K)+FX(I,J,K)*PP(I+1,J,K))
C
      PPTE=(1.-FZ(I,J,K))*
     1      ((1.-FX(I,J,K))*PP(I,J,K)+FX(I,J,K)*PP(I+1,J,K)) +
     1   FZ(I,J,K)*
     1     ((1.-FX(I,J,K+1))*PP(I,J,K+1)+FX(I,J,K+1)*PP(I+1,J,K+1))
C
      PPBE=(1.-FZ(I,J,K-1))*
     1      ((1.-FX(I,J,K-1))*PP(I,J,K-1)+FX(I,J,K-1)*PP(I+1,J,K-1)) +
     1   FZ(I,J,K-1)*
     1     ((1.-FX(I,J,K))*PP(I,J,K)+FX(I,J,K)*PP(I+1,J,K))
C
      PPNW=(1.-FY(I-1,J,K))*
     1         ((1.-FX(I-1,J,K))*PP(I-1,J,K)+FX(I-1,J,K)*PP(I,J,K)) +
     1   FY(I-1,J,K)*
     1     ((1-FX(I-1,J+1,K))*PP(I-1,J+1,K)+FX(I-1,J+1,K)*PP(I,J+1,K))
C 
      PPSW=(1.-FY(I-1,J-1,K))*
cc     1   ((1.-FX(I-1,J-1,K))*PP(I-1,J-1,K)+FX(I-1,J-1,K)*PP(I-1,J,K))+
     1   ((1.-FX(I-1,J-1,K))*PP(I-1,J-1,K)+FX(I-1,J-1,K)*PP(I,J-1,K))+
     1   FY(I-1,J-1,K)*
     1     ((1-FX(I-1,J,K))*PP(I-1,J,K)+FX(I-1,J,K)*PP(I,J,K))
C
      PPTW=(1.-FZ(I-1,J,K))*
     1      ((1.-FX(I-1,J,K))*PP(I-1,J,K)+FX(I-1,J,K)*PP(I,J,K)) +
     1   FZ(I-1,J,K)*
     1    ((1.-FX(I-1,J,K+1))*PP(I-1,J,K+1)+FX(I-1,J,K+1)*PP(I,J,K+1))
C
      PPBW=(1.-FZ(I-1,J,K-1))*
     1    ((1.-FX(I-1,J,K-1))*PP(I-1,J,K-1)+FX(I-1,J,K-1)*PP(I,J,K-1))+
     1   FZ(I-1,J,K-1)*
     1     ((1.-FX(I-1,J,K))*PP(I-1,J,K)+FX(I-1,J,K)*PP(I,J,K))
C
      PPTN=(1.-FZ(I,J,K))*
     1   ((1.-FY(I,J,K))*PP(I,J,K)+FY(I,J,K)*PP(I,J+1,K))+
     1      FZ(I,J,K)*
     1   ((1-FY(I,J,K+1))*PP(I,J,K+1)+FY(I,J,K+1)*PP(I,J+1,K+1))
C
      PPBN=(1.-FZ(I,J,K-1))*
     1   ((1.-FY(I,J,K-1))*PP(I,J,K-1)+FY(I,J,K-1)*PP(I,J+1,K-1))+
     1      FZ(I,J,K-1)*
     1   ((1-FY(I,J,K))*PP(I,J,K)+FY(I,J,K)*PP(I,J+1,K))
C
      PPTS=(1.-FZ(I,J-1,K))*
     1   ((1.-FY(I,J-1,K))*PP(I,J-1,K)+FY(I,J-1,K)*PP(I,J,K))+
     1      FZ(I,J-1,K)*
     1   ((1-FY(I,J-1,K+1))*PP(I,J-1,K+1)+FY(I,J-1,K+1)*PP(I,J,K+1))
C
      PPBS=(1.-FZ(I,J-1,K-1))*
     1 ((1.-FY(I,J-1,K-1))*PP(I,J-1,K-1)+FY(I,J-1,K-1)*PP(I,J,K-1))+
     1      FZ(I,J-1,K-1)*
     1   ((1-FY(I,J-1,K))*PP(I,J-1,K)+FY(I,J-1,K)*PP(I,J,K))

c      IF(J.EQ.2) THEN
c      PPSE=
c     1         ((1.-FX(I,J-1,K))*PP(I,J-1,K)+FX(I,J-1,K)*PP(I,J,K)) 
c      PPSW=
c     1   ((1.-FX(I-1,J-1,K))*PP(I-1,J-1,K)+FX(I-1,J-1,K)*PP(I-1,J,K))
c      PPSE=
c     1         ((1.-FX(I,J,K))*PP(I,J-1,K)+FX(I,J,K)*PP(I,J,K)) 
c      PPSW=
c     1   ((1.-FX(I,J,K))*PP(I-1,J-1,K)+FX(I,J,K)*PP(I-1,J,K))
c      ENDIF
c      IF(J.EQ.NJM1) THEN
c*      PPNW=
c*     1     ((1-FX(I-1,J+1,K))*PP(I-1,J+1,K)+FX(I-1,J+1,K)*PP(I,J+1,K))
c*      PPNE=
c*     1     ((1-FX(I,J+1,K))*PP(I,J+1,K)+FX(I,J+1,K)*PP(I+1,J+1,K))
c      PPNW=
c     1     ((1-FX(I,J,K))*PP(I-1,J+1,K)+FX(I,J,K)*PP(I,J+1,K))
c      PPNE=
c     1     ((1-FX(I,J,K))*PP(I,J+1,K)+FX(I,J,K)*PP(I+1,J+1,K))
c      ENDIF
c      IF(K.EQ.NKM1) THEN
c*      PPTW=
c*     1    ((1.-FX(I-1,J,K+1))*PP(I-1,J,K+1)+FX(I-1,J,K+1)*PP(I,J,K+1))
c*      PPTE=
c*     1     ((1.-FX(I,J,K+1))*PP(I,J,K+1)+FX(I,J,K+1)*PP(I+1,J,K+1))
c      PPTW=
c     1    ((1.-FX(I,J,K))*PP(I-1,J,K+1)+FX(I,J,K)*PP(I,J,K+1))
c      PPTE=
c     1     ((1.-FX(I,J,K))*PP(I,J,K+1)+FX(I,J,K)*PP(I+1,J,K+1))
c      ENDIF
c      IF(K.EQ.NKM1) THEN
c*      PPBE=
c*     1      ((1.-FX(I,J,K-1))*PP(I,J,K-1)+FX(I,J,K-1)*PP(I+1,J,K-1)) 
c*      PPBW=
c*     1    ((1.-FX(I-1,J,K-1))*PP(I-1,J,K-1)+FX(I-1,J,K-1)*PP(I,J,K-1))
c      PPBE=
c     1      ((1.-FX(I,J,K))*PP(I,J,K-1)+FX(I,J,K)*PP(I+1,J,K-1)) 
c      PPBW=
c     1    ((1.-FX(I,J,K))*PP(I-1,J,K-1)+FX(I,J,K)*PP(I,J,K-1))
c      ENDIF
C
      DPNSE = PPNE - PPSE
      DPTBE = PPTE - PPBE 
      DPNSW = PPNW - PPSW
      DPTBW = PPTW - PPBW
      DPEWN = PPNE - PPNW
      DPTBN = PPTN - PPBN
      DPEWS = PPSE - PPSW
      DPTBS = PPTS - PPBS
      DPEWT = PPTE - PPTW
      DPNST = PPTN - PPTS
      DPEWB = PPBE - PPBW
      DPNSB = PPBN - PPBS
C
      SUPP(I,J,K)= 
     1             AE2*DPNSE + AE3*DPTBE - AW2*DPNSW - AW3*DPTBW
     1           + AN1*DPEWN + AN3*DPTBN - AS1*DPEWS - AS3*DPTBS
     1           + AT1*DPEWT + AT2*DPNST - AB1*DPEWB - AB2*DPNSB
C
      SU(I,J,K)=SU(I,J,K)+SUPP(I,J,K)*0.0
C
	ENDDO
	ENDDO
	ENDDO
  401   continue
C
C
      DO 40 N=1,NSWPP
      CALL LISOLV(4,2,2,NI,NJ,NK,PP)
C      CALL LISOLV(4,NJR+1,2,NI,NJRN+1,NK,PP)
C      CALL LISOLV(NIR+1,2,2,NI-2,NJ,NK,PP)
 40   CONTINUE
C
      PPCN=PP(3,3,3)
      DO 41 J=1,NJ
      DO 41 K=1,NK
C      PP(NI,J,K)=PP(NI-1,J,K)+0.5*(PP(NI-1,J,K)-PP(NI-2,J,K))
C      PP(NI,J,K)=(1.-FX(NI-2,J,K))*PP(NI-2,J,K)+FX(NI-2,J,K)
C     +           *PP(NIM1,J,K)
C      PP(1,J,K)=PP(2,J,K)-0.5*(PP(3,J,K)-PP(2,J,K))
      PP(2,J,K)=PP(3,J,K)+PP(NI-2,J,K)-PP(NI-1,J,K)
C     ++(PP(NI-2,J,K)-PP(NI-1,J,K))
      PP(3,J,K)=PP(4,J,K)
      PP(NI,J,K)=PP(NI-1,J,K)
C(1.-FX(NI-1,J,K))*PP(NI-2,J,K)+
C     1                FX(NI-1,J,K)*PP(NI-1,J,K)
 41   CONTINUE
C
      DO 42 K=1,NK
      DO 42 I=1,NI
      PP(I,1,K)=PP(I,2,K)-0.5*(PP(I,3,K)-PP(I,2,K))
      PP(I,NJ,K)=PP(I,NJ-1,K)+0.5*(PP(I,NJ-1,K)-PP(I,NJ-2,K))
 42   CONTINUE
C
      DO 43 I=1,NI
      DO 43 J=1,NJ
      PP(I,J,1)=PP(I,J,2)-.5*(PP(I,J,3)-PP(I,J,2))
      PP(I,J,NK)=PP(I,J,NK-1)+.5*(PP(I,J,NK-1)-PP(I,J,NK-2))
 43   CONTINUE
C     CORRECT CELL CENTER VELOCITY
C
      DO 50 I=4,NI-1
      DO 50 J=2,NJ-1
      DO 50 K=2,NK-1
      PPE=PP(I,J,K)*(1.-FX(I,J,K))+PP(I+1,J,K)*FX(I,J,K)
      PPW=PP(I-1,J,K)*(1.-FX(I-1,J,K))+PP(I,J,K)*FX(I-1,J,K)
      PPT=PP(I,J,K)*(1.-FZ(I,J,K))+PP(I,J,K+1)*FZ(I,J,K)
      PPB=PP(I,J,K-1)*(1.-FZ(I,J,K-1))+PP(I,J,K)*FZ(I,J,K-1)
C
      U(I,J,K)=U(I,J,K)+DU1(I,J,K)*(PPW-PPE)+DU3(I,J,K)*(PPB-PPT)
      V(I,J,K)=V(I,J,K)+DV1(I,J,K)*(PPW-PPE)+DV3(I,J,K)*(PPB-PPT)
      W(I,J,K)=W(I,J,K)+DW1(I,J,K)*(PPW-PPE)+DW3(I,J,K)*(PPB-PPT)
 50   CONTINUE
      DO 51 I=4,NI-1
      DO 51 J=2,NJ-1
      DO 51 K=2,NK-1
      PPN=PP(I,J,K)*(1.-FY(I,J,K))+PP(I,J+1,K)*FY(I,J,K)
      PPS=PP(I,J-1,K)*(1.-FY(I,J-1,K))+PP(I,J,K)*FY(I,J-1,K)
C
      U(I,J,K)=U(I,J,K)+DU2(I,J,K)*(PPS-PPN)
      V(I,J,K)=V(I,J,K)+DV2(I,J,K)*(PPS-PPN)
      W(I,J,K)=W(I,J,K)+DW2(I,J,K)*(PPS-PPN)
 51   CONTINUE
C
C     CORRECT CELL SURFACE VELOCITY
C
      DO 60 I=4,NI-1
      DO 60 J=2,NJ-1
      DO 60 K=2,NK-1
C
      UEC=((1.-FX(I,J,K))*DU1(I,J,K)+FX(I,J,K)*DU1(I+1,J,K))*
     1 (PP(I,J,K)-PP(I+1,J,K))
      VEC=((1.-FX(I,J,K))*DV1(I,J,K)+FX(I,J,K)*DV1(I+1,J,K))*
     1 (PP(I,J,K)-PP(I+1,J,K))
      WEC=((1.-FX(I,J,K))*DW1(I,J,K)+FX(I,J,K)*DW1(I+1,J,K))*
     1 (PP(I,J,K)-PP(I+1,J,K))
C
      FLOW=UEC*((Y(I,J,K)-Y(I,J-1,K-1))*(Z(I,J-1,K)-Z(I,J,K-1))
     1         -(Z(I,J,K)-Z(I,J-1,K-1))*(Y(I,J-1,K)-Y(I,J,K-1)))
     1    +VEC*((Z(I,J,K)-Z(I,J-1,K-1))*(X(I,J-1,K)-X(I,J,K-1))
     1         -(X(I,J,K)-X(I,J-1,K-1))*(Z(I,J-1,K)-Z(I,J,K-1)))
     1    +WEC*((X(I,J,K)-X(I,J-1,K-1))*(Y(I,J-1,K)-Y(I,J,K-1))
     1         -(Y(I,J,K)-Y(I,J-1,K-1))*(X(I,J-1,K)-X(I,J,K-1)))
C
       CX(I,J,K)=0.5*FLOW+CX(I,J,K)
C
      UNC=((1.-FY(I,J,K))*DU2(I,J,K)+FY(I,J,K)*DU2(I,J+1,K))*
     1 (PP(I,J,K)-PP(I,J+1,K))
      VNC=((1.-FY(I,J,K))*DV2(I,J,K)+FY(I,J,K)*DV2(I,J+1,K))*
     1 (PP(I,J,K)-PP(I,J+1,K))
      WNC=((1.-FY(I,J,K))*DW2(I,J,K)+FY(I,J,K)*DW2(I,J+1,K))*
     1 (PP(I,J,K)-PP(I,J+1,K))
C
      FLOW=UNC*((Y(I,J,K)-Y(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K))
     1         -(Z(I,J,K)-Z(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K)))
     1    +VNC*((Z(I,J,K)-Z(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K))
     1         -(X(I,J,K)-X(I-1,J,K-1))*(Z(I,J,K-1)-Z(I-1,J,K)))
     1    +WNC*((X(I,J,K)-X(I-1,J,K-1))*(Y(I,J,K-1)-Y(I-1,J,K))
     1         -(Y(I,J,K)-Y(I-1,J,K-1))*(X(I,J,K-1)-X(I-1,J,K)))
      CY(I,J,K)=0.5*FLOW+CY(I,J,K)
C
      UTC=((1.-FZ(I,J,K))*DU3(I,J,K)+FZ(I,J,K)*DU3(I,J,K+1))*
     1 (PP(I,J,K)-PP(I,J,K+1))
      VTC=((1.-FZ(I,J,K))*DV3(I,J,K)+FZ(I,J,K)*DV3(I,J,K+1))*
     1 (PP(I,J,K)-PP(I,J,K+1))
      WTC=((1.-FZ(I,J,K))*DW3(I,J,K)+FZ(I,J,K)*DW3(I,J,K+1))*
     1 (PP(I,J,K)-PP(I,J,K+1))
C
      FLOW=UTC*((Y(I,J,K)-Y(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K))
     1         -(Z(I,J,K)-Z(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K)))
     1    +VTC*((Z(I,J,K)-Z(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K))
     1         -(X(I,J,K)-X(I-1,J-1,K))*(Z(I-1,J,K)-Z(I,J-1,K)))
     1    +WTC*((X(I,J,K)-X(I-1,J-1,K))*(Y(I-1,J,K)-Y(I,J-1,K))
     1         -(Y(I,J,K)-Y(I-1,J-1,K))*(X(I-1,J,K)-X(I,J-1,K)))
       CZ(I,J,K)=0.5*FLOW+CZ(I,J,K)
 60   CONTINUE
C
      PPREF=PP(IPREF,JPREF,KPREF)
      DO 70 I=4,NI-1
      DO 70 J=2,NJ-1
      DO 70 K=2,NK-1
      P(I,J,K)=P(I,J,K)+URFP*(PP(I,J,K)-PPREF)
 70   CONTINUE
C
      PCN=P(2,3,3)
      DO 71 J=1,NJ
      DO 71 K=1,NK
C      P(1,J,K)=P(2,J,K)
C      P(NI,J,K)=P(NI-1,J,K)
C      P(NI,J,K)=P(NI-1,J,K)+0.5*(P(NI-1,J,K)-P(NI-2,J,K))
C      P(1,J,K)=(1.-FX(2,J,K))*P(2,J,K)+FX(2,J,K)*P(3,J,K)-DPEX
 71   CONTINUE
C
      DO 72 K=1,NK
      DO 72 I=1,NI
      P(I,1,K)=P(I,2,K)
      P(I,NJ,K)=P(I,NJ-1,K)
C*      TOARA=sqrt(yx(i,2,k)**2.+yy(i,2,k)**2.+yz(i,2,k)**2.)
C*      UN2=(u(i,2,k)*yx(i,2,k)+v(i,2,k)*yy(i,2,k)+
C*     &w(i,2,k)*yz(i,2,k))/TOARA
C*      P(I,1,K)=P(I,2,K)+DEN(I,2,K)*ABS(UN2)*UN2
C*      TOARA=sqrt(yx(i,NJM1,k)**2.+yy(i,NJM1,k)**2.+yz(i,NJM1,k)**2.)
C*      UNNJM1=(u(i,NJM1,k)*yx(i,NJM1,k)+v(i,NJM1,k)*yy(i,NJM1,k)+
C*     &w(i,NJM1,k)*yz(i,NJM1,k))/TOARA
C*      P(I,NJ,K)=P(I,NJM1,K)+DEN(I,NJM1,K)*ABS(UNNJM1)*UNNJM1
*      P(I,1,K)=P(I,2,K)-0.5*(P(I,3,K)-P(I,2,K))
*      P(I,NJ,K)=P(I,NJ-1,K)+0.5*(P(I,NJ-1,K)-P(I,NJ-2,K))
 72   CONTINUE
C
      DO 73 I=1,NI
      DO 73 J=1,NJ
      P(I,J,1)=P(I,J,2)-.5*(P(I,J,3)-P(I,J,2))
      P(I,J,NK)=P(I,J,NK-1)+.5*(P(I,J,NK-1)-P(I,J,NK-2))
      P(I,J,1)=P(I,J,2)-.5*DEN(I,J,2)*W(I,J,2)*ABS(W(I,J,2))
      P(I,J,NK)=P(I,J,NKM1)+.5*DEN(I,J,NKM1)*
     &W(I,J,NKM1)*ABS(W(I,J,NKM1))
 73   CONTINUE
      CALL MODU
      CALL MODP
      CALL MODXU
C
      RETURN
      END
C
      SUBROUTINE CALCT
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      PR=PRTTM
C
      DO 2 K=2,NKM1
      DO 2 J=2,NJM1
      DO 2 I=2,NIM1
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
      SMP= DENE*CX(I,J,K)-DENW*CX(I-1,J,K)
     1    +DENN*CY(I,J,K)-DENS*CY(I,J-1,K)
     1    +DENT*CZ(I,J,K)-DENB*CZ(I,J,K-1)
      CP=ABS(SMP)
      SU(I,J,K)=CP*T(I,J,K)
      SP(I,J,K)=-CP
    2 CONTINUE
C
C ------ CALCULATE ADVECTION COEFFICIENTS -----------------------
C
      CALL ADVECT(PRLTM,PRTTM)
C
      CALL BQUICK(0,1.,T)
C
C ------ CALCULATE CROSS-DIFFUSION SOURCE TERMS -----------------
C
      CALL CRDIFF(PRLTM,PRTTM,0,T)
C
      CALL MODTT
C
      RESORT=0.
      DO 30 K=2,NKM1
      DO 30 J=2,NJM1
      DO 30 I=2,NIM1
      AP(I,J,K)=AE(I,J,K)+AW(I,J,K)+AN(I,J,K)+AS(I,J,K)+AT(I,J,K)
     1 +AB(I,J,K)-SP(I,J,K)
      RESOR=AE(I,J,K)*T(I+1,J,K)+AW(I,J,K)*T(I-1,J,K)
     1     +AN(I,J,K)*T(I,J+1,K)+AS(I,J,K)*T(I,J-1,K)
     1     +AT(I,J,K)*T(I,J,K+1)+AB(I,J,K)*T(I,J,K-1)+SU(I,J,K)
     1     -AP(I,J,K)*T(I,J,K)
C
      IF(AP(I,J,K).GE.GREAT) RESOR=0.
      RESORT=RESORT+ABS(RESOR)
C
C     UNDERRELAX
C
      AP(I,J,K)=AP(I,J,K)/URFT
      SU(I,J,K)=SU(I,J,K)+(1.-URFT)*AP(I,J,K)*T(I,J,K)
 30   CONTINUE
C
      DO 40 N=1,NSWPT
      CALL LISOLV(3,2,2,NI,NJ,NK,T)
 40   CONTINUE
C      CALL MODTT
C      TCOR=300.-T(2,NK,2)
C      DO 41 K=1,NK
C      DO 41 J=1,NJ
C      DO 41 I=1,NI
C   41 T(I,J,K)=T(I,J,K)+TCOR
C
      RETURN
      END
C
      SUBROUTINE GRID1
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
C      OPEN(13,FILE='XGRID')
C      OPEN(14,FILE='YGRID')
C      OPEN(15,FILE='ZGRID')
      NIM2=NIM1-1
      NJM2=NJM1-1
      NKM2=NKM1-1
C
      XLREF=100.
C     DX=3.33333 /NIM2
      DX=XLREF/3./NIM2
      DY=1./NJM2
      DZ=1./NKM2
C
C     INITIALIZATION
C
      DO 5 I=1,NI
      DO 5 J=1,NJ
      DO 5 K=1,NK
      X(I,J,K)=0.
      Y(I,J,K)=0.
      Z(I,J,K)=0.
      XC(I,J,K)=0.
      YC(I,J,K)=0.
      ZC(I,J,K)=0.
C
      XX(I,J,K)=0.
      XY(I,J,K)=0.
      XZ(I,J,K)=0.
      YX(I,J,K)=0.
      YY(I,J,K)=0.
      YZ(I,J,K)=0.
      ZX(I,J,K)=0.
      ZY(I,J,K)=0.
      ZZ(I,J,K)=0.
      XJACOB(I,J,K)=0.
      VOL(I,J,K)=0.
C
      FX(I,J,K)=0.
      FY(I,J,K)=0.
      FZ(I,J,K)=0.
 5    CONTINUE
C
C      READ(13,*)(((X(I,J,K),I=1,NIM1),J=1,NJM1),K=1,NKM1)
C      READ(14,*)(((Y(I,J,K),I=1,NIM1),J=1,NJM1),K=1,NKM1)
C      READ(15,*)(((Z(I,J,K),I=1,NIM1),J=1,NJM1),K=1,NKM1)
      READ(13,*)X
      READ(14,*)Y
      READ(15,*)Z
C
C     CALCULATE CELL CENTER COORDINATES
C
      DO 10 K=2,NKM1
      DO 10 J=2,NJM1
      DO 10 I=2,NI
C      DO 10 I=2,NIM1
      XC(I,J,K)=0.125*(X(I-1,J-1,K-1)+X(I-1,J,K-1)+X(I,J,K-1)
     1 +X(I,J-1,K-1)+X(I-1,J-1,K)+X(I-1,J,K)+X(I,J,K)+X(I,J-1,K))
      YC(I,J,K)=0.125*(Y(I-1,J-1,K-1)+Y(I-1,J,K-1)+Y(I,J,K-1)
     1 +Y(I,J-1,K-1)+Y(I-1,J-1,K)+Y(I-1,J,K)+Y(I,J,K)+Y(I,J-1,K))
      ZC(I,J,K)=0.125*(Z(I-1,J-1,K-1)+Z(I-1,J,K-1)+Z(I,J,K-1)
     1 +Z(I,J-1,K-1)+Z(I-1,J-1,K)+Z(I-1,J,K)+Z(I,J,K)+Z(I,J-1,K))
 10   CONTINUE
C
C   BOUNDARY VALUES FOR CELL CENTER COORDINATES
C
      DO 910 J=2,NJM1
      DO 910 I=2,NI
C      DO 910 I=2,NIM1
      XC(I,J,1)=.25*(X(I,J,1)+X(I-1,J,1)+X(I,J-1,1)+X(I-1,J-1,1))
      YC(I,J,1)=.25*(Y(I,J,1)+Y(I-1,J,1)+Y(I,J-1,1)+Y(I-1,J-1,1))
      ZC(I,J,1)=.25*(Z(I,J,1)+Z(I-1,J,1)+Z(I,J-1,1)+Z(I-1,J-1,1))
C
      XC(I,J,NK)=.25*(X(I,J,NKM1)+X(I-1,J,NKM1)+X(I,J-1,NKM1)+
     +                X(I-1,J-1,NKM1))
      YC(I,J,NK)=.25*(Y(I,J,NKM1)+Y(I-1,J,NKM1)+Y(I,J-1,NKM1)+
     +                Y(I-1,J-1,NKM1))
      ZC(I,J,NK)=.25*(Z(I,J,NKM1)+Z(I-1,J,NKM1)+Z(I,J-1,NKM1)+
     +                Z(I-1,J-1,NKM1))
  910 CONTINUE
C
      DO 911 K=2,NKM1
C      DO 911 I=2,NIM1
      DO 911 I=2,NI
      XC(I,1,K)=.25*(X(I,1,K)+X(I-1,1,K)+X(I,1,K-1)+X(I-1,1,K-1))
      YC(I,1,K)=.25*(Y(I,1,K)+Y(I-1,1,K)+Y(I,1,K-1)+Y(I-1,1,K-1))
      ZC(I,1,K)=.25*(Z(I,1,K)+Z(I-1,1,K)+Z(I,1,K-1)+Z(I-1,1,K-1))
C
      XC(I,NJ,K)=.25*(X(I,NJM1,K)+X(I-1,NJM1,K)+X(I,NJM1,K-1)+
     +                X(I-1,NJM1,K-1))
      YC(I,NJ,K)=.25*(Y(I,NJM1,K)+Y(I-1,NJM1,K)+Y(I,NJM1,K-1)+
     +                Y(I-1,NJM1,K-1))
      ZC(I,NJ,K)=.25*(Z(I,NJM1,K)+Z(I-1,NJM1,K)+Z(I,NJM1,K-1)+
     +                Z(I-1,NJM1,K-1))
  911 CONTINUE
C
      DO 912 K=2,NKM1
      DO 912 J=2,NJM1
      XC(1,J,K)=.25*(X(1,J,K)+X(1,J-1,K)+X(1,J,K-1)+X(1,J-1,K-1))
      YC(1,J,K)=.25*(Y(1,J,K)+Y(1,J-1,K)+Y(1,J,K-1)+Y(1,J-1,K-1))
      ZC(1,J,K)=.25*(Z(1,J,K)+Z(1,J-1,K)+Z(1,J,K-1)+Z(1,J-1,K-1))
C
C      XC(NI,J,K)=.25*(X(NIM1,J,K)+X(NIM1,J-1,K)+X(NIM1,J,K-1)+
C     +                X(NIM1,J-1,K-1))
C      YC(NI,J,K)=.25*(Y(NIM1,J,K)+Y(NIM1,J-1,K)+Y(NIM1,J,K-1)+
C     +                Y(NIM1,J-1,K-1))
C      ZC(NI,J,K)=.25*(Z(NIM1,J,K)+Z(NIM1,J-1,K)+Z(NIM1,J,K-1)+
C     +                Z(NIM1,J-1,K-1))
  912 CONTINUE
C
      XC(1,1,1)=X(1,1,1)
      YC(1,1,1)=Y(1,1,1)
      ZC(1,1,1)=Z(1,1,1)
      XC(1,NJ,1)=X(1,NJM1,1)
      YC(1,NJ,1)=Y(1,NJM1,1)
      ZC(1,NJ,1)=Z(1,NJM1,1)
C      XC(NI,1,1)=X(NIM1,1,1)
C      YC(NI,1,1)=Y(NIM1,1,1)
C      ZC(NI,1,1)=Z(NIM1,1,1)
C      XC(NI,NJ,1)=X(NIM1,NJM1,1)
C      YC(NI,NJ,1)=Y(NIM1,NJM1,1)
C      ZC(NI,NJ,1)=Z(NIM1,NJM1,1)
      DO 113 K=2,NKM1
      XC(1,1,K)=.5*(X(1,1,K)+X(1,1,K-1))
      YC(1,1,K)=.5*(Y(1,1,K)+Y(1,1,K-1))
      ZC(1,1,K)=.5*(Z(1,1,K)+Z(1,1,K-1))
      XC(1,NJ,K)=.5*(X(1,NJM1,K)+X(1,NJM1,K-1))
      YC(1,NJ,K)=.5*(Y(1,NJM1,K)+Y(1,NJM1,K-1))
      ZC(1,NJ,K)=.5*(Z(1,NJM1,K)+Z(1,NJM1,K-1))
      XC(NI,1,K)=.5*(X(NIM1,1,K)+X(NIM1,1,K-1))
      YC(NI,1,K)=.5*(Y(NIM1,1,K)+Y(NIM1,1,K-1))
      ZC(NI,1,K)=.5*(Z(NIM1,1,K)+Z(NIM1,1,K-1))
C      XC(NI,NJ,K)=.5*(X(NIM1,NJM1,K)+X(NIM1,NJM1,K-1))
C      YC(NI,NJ,K)=.5*(Y(NIM1,NJM1,K)+Y(NIM1,NJM1,K-1))
C      ZC(NI,NJ,K)=.5*(Z(NIM1,NJM1,K)+Z(NIM1,NJM1,K-1))
  113 CONTINUE
      XC(1,1,NK)=X(1,1,NKM1)
      YC(1,1,NK)=Y(1,1,NKM1)
      ZC(1,1,NK)=Z(1,1,NKM1)
      XC(1,NJ,NK)=X(1,NJM1,NKM1)
      YC(1,NJ,NK)=Y(1,NJM1,NKM1)
      ZC(1,NJ,NK)=Z(1,NJM1,NKM1)
C      XC(NI,1,NK)=X(NIM1,1,NKM1)
C      YC(NI,1,NK)=Y(NIM1,1,NKM1)
C     ZC(NI,1,NK)=Z(NIM1,1,NKM1)
C      XC(NI,NJ,NK)=X(NIM1,NJM1,NKM1)
C      YC(NI,NJ,NK)=Y(NIM1,NJM1,NKM1)
C      ZC(NI,NJ,NK)=Z(NIM1,NJM1,NKM1)
C
      DO 114 J=2,NJM1
      XC(1,J,1)=.5*(X(1,J,1)+X(1,J-1,1))
      YC(1,J,1)=.5*(Y(1,J,1)+Y(1,J-1,1))
      ZC(1,J,1)=.5*(Z(1,J,1)+Z(1,J-1,1))
      XC(NI,J,1)=.5*(X(NIM1,J,1)+X(NIM1,J-1,1))
      YC(NI,J,1)=.5*(Y(NIM1,J,1)+Y(NIM1,J-1,1))
      ZC(NI,J,1)=.5*(Z(NIM1,J,1)+Z(NIM1,J-1,1))
      XC(1,J,NK)=.5*(X(1,J,NKM1)+X(1,J-1,NKM1))
      YC(1,J,NK)=.5*(Y(1,J,NKM1)+Y(1,J-1,NKM1))
      ZC(1,J,NK)=.5*(Z(1,J,NKM1)+Z(1,J-1,NKM1))
C      XC(NI,J,NK)=.5*(X(NIM1,J,NKM1)+X(NIM1,J-1,NKM1))
C      YC(NI,J,NK)=.5*(Y(NIM1,J,NKM1)+Y(NIM1,J-1,NKM1))
C      ZC(NI,J,NK)=.5*(Z(NIM1,J,NKM1)+Z(NIM1,J-1,NKM1))
  114 CONTINUE
C
      DO 115 I=2,NIM1
      XC(I,1,1)=0.5*(X(I,1,1)+X(I-1,1,1))
      YC(I,1,1)=0.5*(Y(I,1,1)+Y(I-1,1,1))
      ZC(I,1,1)=0.5*(Z(I,1,1)+Z(I-1,1,1))
      XC(I,NJ,1)=0.5*(X(I,NJM1,1)+X(I-1,NJM1,1))
      YC(I,NJ,1)=0.5*(Y(I,NJM1,1)+Y(I-1,NJM1,1))
      ZC(I,NJ,1)=0.5*(Z(I,NJM1,1)+Z(I-1,NJM1,1))
      XC(I,1,NK)=0.5*(X(I,1,NKM1)+X(I-1,1,NKM1))
      YC(I,1,NK)=0.5*(Y(I,1,NKM1)+Y(I-1,1,NKM1))
      ZC(I,1,NK)=0.5*(Z(I,1,NKM1)+Z(I-1,1,NKM1))
      XC(I,NJ,NK)=0.5*(X(I,NJM1,NKM1)+X(I-1,NJM1,NKM1))
      YC(I,NJ,NK)=0.5*(Y(I,NJM1,NKM1)+Y(I-1,NJM1,NKM1))
      ZC(I,NJ,NK)=0.5*(Z(I,NJM1,NKM1)+Z(I-1,NJM1,NKM1))
  115 CONTINUE
C
C     CALCULATE METRIC AND JACOBIAN
C
      DO 20 I=2,NIM1
      DO 20 J=2,NJM1
      DO 20 K=2,NKM1
C
      XE=0.25*(X(I,J,K-1)+X(I,J-1,K-1)+X(I,J,K)+X(I,J-1,K))
      YE=0.25*(Y(I,J,K-1)+Y(I,J-1,K-1)+Y(I,J,K)+Y(I,J-1,K))
      ZE=0.25*(Z(I,J,K-1)+Z(I,J-1,K-1)+Z(I,J,K)+Z(I,J-1,K))
C
      XW=0.25*(X(I-1,J-1,K-1)+X(I-1,J,K-1)+X(I-1,J-1,K)+X(I-1,J,K))
      YW=0.25*(Y(I-1,J-1,K-1)+Y(I-1,J,K-1)+Y(I-1,J-1,K)+Y(I-1,J,K))
      ZW=0.25*(Z(I-1,J-1,K-1)+Z(I-1,J,K-1)+Z(I-1,J-1,K)+Z(I-1,J,K))
C
      XN=0.25*(X(I-1,J,K-1)+X(I,J,K-1)+X(I-1,J,K)+X(I,J,K))
      YN=0.25*(Y(I-1,J,K-1)+Y(I,J,K-1)+Y(I-1,J,K)+Y(I,J,K))
      ZN=0.25*(Z(I-1,J,K-1)+Z(I,J,K-1)+Z(I-1,J,K)+Z(I,J,K))
C
      XS=0.25*(X(I-1,J-1,K-1)+X(I,J-1,K-1)+X(I-1,J-1,K)+X(I,J-1,K))
      YS=0.25*(Y(I-1,J-1,K-1)+Y(I,J-1,K-1)+Y(I-1,J-1,K)+Y(I,J-1,K))
      ZS=0.25*(Z(I-1,J-1,K-1)+Z(I,J-1,K-1)+Z(I-1,J-1,K)+Z(I,J-1,K))
C
      XT=0.25*(X(I-1,J-1,K)+X(I-1,J,K)+X(I,J,K)+X(I,J-1,K))
      YT=0.25*(Y(I-1,J-1,K)+Y(I-1,J,K)+Y(I,J,K)+Y(I,J-1,K))
      ZT=0.25*(Z(I-1,J-1,K)+Z(I-1,J,K)+Z(I,J,K)+Z(I,J-1,K))
C
      XB=0.25*(X(I-1,J-1,K-1)+X(I-1,J,K-1)+X(I,J,K-1)+X(I,J-1,K-1))
      YB=0.25*(Y(I-1,J-1,K-1)+Y(I-1,J,K-1)+Y(I,J,K-1)+Y(I,J-1,K-1))
      ZB=0.25*(Z(I-1,J-1,K-1)+Z(I-1,J,K-1)+Z(I,J,K-1)+Z(I,J-1,K-1))
C
      XX(I,J,K)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
      XY(I,J,K)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
      XZ(I,J,K)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
      YX(I,J,K)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
      YY(I,J,K)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
      YZ(I,J,K)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
      ZX(I,J,K)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
      ZY(I,J,K)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
      ZZ(I,J,K)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
      XJACOB(I,J,K)=(XE-XW)*(YN-YS)*(ZT-ZB)+(XT-XB)*(YE-YW)*(ZN-ZS)
     1 +(XN-XS)*(YT-YB)*(ZE-ZW)
     1             -(XE-XW)*(YT-YB)*(ZN-ZS)-(XT-XB)*(YN-YS)*(ZE-ZW)
     1 -(XN-XS)*(YE-YW)*(ZT-ZB)
C
 20   CONTINUE
C
C     INTERPOLATION FACTOR IN X,Y,Z DIRECTION
C
      DO 30 I=2,NIM1
      DO 30 J=2,NJM1
      DO 30 K=2,NKM1
C
      XE=0.25*(X(I,J,K-1)+X(I,J-1,K-1)+X(I,J,K)+X(I,J-1,K))
      YE=0.25*(Y(I,J,K-1)+Y(I,J-1,K-1)+Y(I,J,K)+Y(I,J-1,K))
      ZE=0.25*(Z(I,J,K-1)+Z(I,J-1,K-1)+Z(I,J,K)+Z(I,J-1,K))
C
      FXBOT=SQRT( (XC(I+1,J,K)-XC(I,J,K))**2+(YC(I+1,J,K)-YC(I,J,K))**2
     1+(ZC(I+1,J,K)-ZC(I,J,K))**2 )
      FXTOP=SQRT( (XE-XC(I,J,K))**2+(YE-YC(I,J,K))**2
     1+(ZE-ZC(I,J,K))**2 )
      FXBOT1=SQRT( (XC(I+1,J,K)-XE)**2+(YC(I+1,J,K)-YE)**2
     1            +(ZC(I+1,J,K)-ZE)**2)
C      FX(I,J,K)=FXTOP/FXBOT
      FX(I,J,K)=FXTOP/(FXTOP+FXBOT1)
 30   CONTINUE
C
      DO 40 I=2,NIM1
      DO 40 J=2,NJ-2
      DO 40 K=2,NKM1
C
      XN=0.25*(X(I-1,J,K-1)+X(I,J,K-1)+X(I-1,J,K)+X(I,J,K))
      YN=0.25*(Y(I-1,J,K-1)+Y(I,J,K-1)+Y(I-1,J,K)+Y(I,J,K))
      ZN=0.25*(Z(I-1,J,K-1)+Z(I,J,K-1)+Z(I-1,J,K)+Z(I,J,K))
C
      FYBOT=SQRT( (XC(I,J+1,K)-XC(I,J,K))**2+(YC(I,J+1,K)-YC(I,J,K))**2
     1+(ZC(I,J+1,K)-ZC(I,J,K))**2 )
      FYTOP=SQRT( (XN-XC(I,J,K))**2+(YN-YC(I,J,K))**2
     1+(ZN-ZC(I,J,K))**2 )
      FYBOT1=SQRT((XC(I,J+1,K)-XN)**2+(YC(I,J+1,K)-YN)**2+
     1            (ZC(I,J+1,K)-ZN)**2)
C      FY(I,J,K)=FYTOP/FYBOT
      FY(I,J,K)=FYTOP/(FYTOP+FYBOT1)
 40   CONTINUE
C
      DO 50 I=2,NIM1
      DO 50 J=2,NJM1
      DO 50 K=2,NK-2
C
      XT=0.25*(X(I-1,J-1,K)+X(I-1,J,K)+X(I,J,K)+X(I,J-1,K))
      YT=0.25*(Y(I-1,J-1,K)+Y(I-1,J,K)+Y(I,J,K)+Y(I,J-1,K))
      ZT=0.25*(Z(I-1,J-1,K)+Z(I-1,J,K)+Z(I,J,K)+Z(I,J-1,K))
C
      FZBOT=SQRT( (XC(I,J,K+1)-XC(I,J,K))**2+(YC(I,J,K+1)-YC(I,J,K))**2
     1+(ZC(I,J,K+1)-ZC(I,J,K))**2 )
      FZTOP=SQRT( (XT-XC(I,J,K))**2+(YT-YC(I,J,K))**2
     1+(ZT-ZC(I,J,K))**2 )
      FZ(I,J,K)=FZTOP/FZBOT
 50   CONTINUE
C
C       BOUNDARY FOR FX,FY,FZ
C
C******************************************************************
C
C*****FX
C
C NS
C0000000
      J=1
CCCC
C      DO 60 I=2,NIM2
      DO 60 I=2,NIM1
      DO 60 K=2,NKM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J,K-1)+Z(I,J,K-1))
C
      XE=0.25*(X(I+1,J,K)+X(I,J,K)+X(I,J,K-1)+X(I+1,J,K-1))
      YE=0.25*(Y(I+1,J,K)+Y(I,J,K)+Y(I,J,K-1)+Y(I+1,J,K-1))
      ZE=0.25*(Z(I+1,J,K)+Z(I,J,K)+Z(I,J,K-1)+Z(I+1,J,K-1))
C
      XM=0.5*(X(I,J,K)+X(I,J,K-1))
      YM=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZM=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J  ,K)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J,K)=FXTOP/(FXTOP+FXBOT1)
60    CONTINUE
C
      K=1
      DO 61 I=2,NIM1
C      DO 61 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J  ,K)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J,K)=FXTOP/(FXTOP+FXBOT1)
 61   CONTINUE
C
C
      K=NKM1
      DO 62 I=2,NIM1
C      DO 62 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J  ,K+1)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J,K+1)=FXTOP/(FXTOP+FXBOT1)
 62   CONTINUE
C
C
C00000000
C
C0000000
      J=NJM1
CCCC
      DO 70 I=2,NIM1
C      DO 70 I=2,NIM2
      DO 70 K=2,NKM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J,K-1)+Z(I,J,K-1))
C
      XE=0.25*(X(I+1,J,K)+X(I,J,K)+X(I,J,K-1)+X(I+1,J,K-1))
      YE=0.25*(Y(I+1,J,K)+Y(I,J,K)+Y(I,J,K-1)+Y(I+1,J,K-1))
      ZE=0.25*(Z(I+1,J,K)+Z(I,J,K)+Z(I,J,K-1)+Z(I+1,J,K-1))
C
      XM=0.5*(X(I,J,K)+X(I,J,K-1))
      YM=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZM=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J+1,K)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J+1,K)=FXTOP/(FXTOP+FXBOT1)
70    CONTINUE
C
      K=1
      DO 71 I=2,NIM1
C      DO 71 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J+1,K)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J+1,K)=FXTOP/(FXTOP+FXBOT1)
 71   CONTINUE
C
C
      K=NKM1
      DO 72 I=2,NIM1
C      DO 72 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J+1,K+1)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J+1,K+1)=FXTOP/(FXTOP+FXBOT1)
 72   CONTINUE
C
C00000000
C
C BOTTOM
C0000000
      K=1
CCCC
      DO 80 I=2,NIM1
C      DO 80 I=2,NIM2
      DO 80 J=2,NJM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J-1,K)+X(I,J-1,K))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J-1,K)+Y(I,J-1,K))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J-1,K)+Z(I,J-1,K))
C
      XE=0.25*(X(I+1,J,K)+X(I,J,K)+X(I,J-1,K)+X(I+1,J-1,K))
      YE=0.25*(Y(I+1,J,K)+Y(I,J,K)+Y(I,J-1,K)+Y(I+1,J-1,K))
      ZE=0.25*(Z(I+1,J,K)+Z(I,J,K)+Z(I,J-1,K)+Z(I+1,J-1,K))
C
      XM=0.5*(X(I,J,K)+X(I,J-1,K))
      YM=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZM=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J,K  )=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J,K)=FXTOP/(FXTOP+FXBOT1)
80    CONTINUE
C
      J=1
      DO 81 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J,K  )=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J,K)=FXTOP/(FXTOP+FXBOT1)
 81   CONTINUE
C
C
      J=NJM1
      DO 82 I=2,NIM1
C      DO 82 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J+1,K  )=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J+1,K)=FXTOP/(FXTOP+FXBOT1)
 82   CONTINUE
C
C00000000
C
C STEP
C0000000
C      GO TO 841
C      K=NKR
CCCC
C      DO 380 I=2,NIR-1
C      DO 380 J=2,NJM1
C      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J-1,K)+X(I,J-1,K))
C      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J-1,K)+Y(I,J-1,K))
C      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J-1,K)+Z(I,J-1,K))
C      XE=0.25*(X(I+1,J,K)+X(I,J,K)+X(I,J-1,K)+X(I+1,J-1,K))
C      YE=0.25*(Y(I+1,J,K)+Y(I,J,K)+Y(I,J-1,K)+Y(I+1,J-1,K))
C      ZE=0.25*(Z(I+1,J,K)+Z(I,J,K)+Z(I,J-1,K)+Z(I+1,J-1,K))
C      XM=0.5*(X(I,J,K)+X(I,J-1,K))
C      YM=0.5*(Y(I,J,K)+Y(I,J-1,K))
C      ZM=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
C      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J,K  )=FXTOP/FXBOT
C 380    CONTINUE
C
C      J=1
C      DO 381 I=2,NIR-1
C      XP=0.5*(X(I,J,K)+X(I-1,J,K))
C      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
C      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
C      XE=0.5*(X(I,J,K)+X(I+1,J,K))
C      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
C      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J,K  )=FXTOP/FXBOT
C 381   CONTINUE
C
C
C      J=NJM1
C      DO 382 I=2,NIR-1
C      XP=0.5*(X(I,J,K)+X(I-1,J,K))
C      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
C      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
C      XE=0.5*(X(I,J,K)+X(I+1,J,K))
C      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
C      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J+1,K  )=FXTOP/FXBOT
C  382 CONTINUE
C  841 CONTINUE
C00000000
C
C  TOP
C00000000
      K=NKM1
CCCC
      DO 90 I=2,NIM1
C      DO 90 I=2,NIM2
      DO 90 J=2,NJM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J-1,K)+X(I,J-1,K))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J-1,K)+Y(I,J-1,K))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J-1,K)+Z(I,J-1,K))
C
      XE=0.25*(X(I+1,J,K)+X(I,J,K)+X(I,J-1,K)+X(I+1,J-1,K))
      YE=0.25*(Y(I+1,J,K)+Y(I,J,K)+Y(I,J-1,K)+Y(I+1,J-1,K))
      ZE=0.25*(Z(I+1,J,K)+Z(I,J,K)+Z(I,J-1,K)+Z(I+1,J-1,K))
C
      XM=0.5*(X(I,J,K)+X(I,J-1,K))
      YM=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZM=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J,K+1)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J,K+1)=FXTOP/(FXTOP+FXBOT1)
90    CONTINUE
C
      J=1
      DO 91 I=2,NIM1
C      DO 91 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J,K+1)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J,K+1)=FXTOP/(FXTOP+FXBOT1)
 91   CONTINUE
C
C
      J=NJM1
      DO 92 I=2,NIM1
C      DO 92 I=2,NIM2
      XP=0.5*(X(I,J,K)+X(I-1,J,K))
      YP=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZP=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      XE=0.5*(X(I,J,K)+X(I+1,J,K))
      YE=0.5*(Y(I,J,K)+Y(I+1,J,K))
      ZE=0.5*(Z(I,J,K)+Z(I+1,J,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FXTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FXBOT=SQRT( (XE-XP)**2+(YE-YP)**2+(ZE-ZP)**2 )
C      FX(I,J+1,K+1)=FXTOP/FXBOT
      FXBOT1=SQRT( (XE-XM)**2+(YE-YM)**2+(ZE-ZM)**2)
      FX(I,J+1,K+1)=FXTOP/(FXTOP+FXBOT1)
 92   CONTINUE
C
C00000000
C
C EW
C
      DO 100 J=1,NJ
      DO 100 K=1,NK
      FX(1,J,K)=0.
c      FX(NIM1,J,K)=1.
 100  CONTINUE
C      DO 93 J=1,NJ
C      DO 93 K=1,NKR
C      FX(NIR,J,K)=0.
C   93 CONTINUE
C******************************************************************
C
C******************************************************************
C
C*****FY
C
C BOTTOM
C0000000
      K=1
CCCC
      DO 110 J=2,NJM2
      DO 110 I=2,NIM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J-1,K)+X(I,J-1,K))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J-1,K)+Y(I,J-1,K))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J-1,K)+Z(I,J-1,K))
C
      XN=0.25*(X(I,J+1,K)+X(I,J,K)+X(I-1,J,K)+X(I-1,J+1,K))
      YN=0.25*(Y(I,J+1,K)+Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J+1,K))
      ZN=0.25*(Z(I,J+1,K)+Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J+1,K))
C
      XM=0.5*(X(I,J,K)+X(I-1,J,K))
      YM=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZM=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I,J,K  )=FYTOP/FYBOT
      FY(I,J,K)=FYTOP/(FYTOP+FYBOT1)
110   CONTINUE
C
      I=1
      DO 111 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I,J,K  )=FYTOP/FYBOT
      FY(I,J,K)=FYTOP/(FYTOP+FYBOT1)
 111  CONTINUE
C
      I=NIM1
      DO 112 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I+1,J,K  )=FYTOP/FYBOT
      FY(I+1,J,K)=FYTOP/(FYTOP+FYBOT1)
 112  CONTINUE
C
C00000000
C STEP
C0000000
C      GO TO 842
C      K=NKR
CCCC
C      DO 410 J=2,NJM2
C      DO 410 I=2,NIR
C      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J-1,K)+X(I,J-1,K))
C      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J-1,K)+Y(I,J-1,K))
C      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J-1,K)+Z(I,J-1,K))
C
C      XN=0.25*(X(I,J+1,K)+X(I,J,K)+X(I-1,J,K)+X(I-1,J+1,K))
C      YN=0.25*(Y(I,J+1,K)+Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J+1,K))
C      ZN=0.25*(Z(I,J+1,K)+Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J+1,K))
C
C      XM=0.5*(X(I,J,K)+X(I-1,J,K))
C      YM=0.5*(Y(I,J,K)+Y(I-1,J,K))
C      ZM=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
C      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I,J,K  )=FYTOP/FYBOT
 410  CONTINUE
C
C      I=1
C      DO 411 J=2,NJM2
C      XP=0.5*(X(I,J,K)+X(I,J-1,K))
C      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
C      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
C      XN=0.5*(X(I,J,K)+X(I,J+1,K))
C      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
C      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I,J,K  )=FYTOP/FYBOT
C 411  CONTINUE
C
C      I=NIR
C      DO 412 J=2,NJM2
C      XP=0.5*(X(I,J,K)+X(I,J-1,K))
C      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
C      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
C      XN=0.5*(X(I,J,K)+X(I,J+1,K))
C     YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
C      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I+1,J,K  )=FYTOP/FYBOT
C 412  CONTINUE
C 842  CONTINUE
C
C0000000
C TOP
C0000000
      K=NKM1
CCCC
      DO 120 J=2,NJM2
      DO 120 I=2,NIM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J-1,K)+X(I,J-1,K))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J-1,K)+Y(I,J-1,K))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J-1,K)+Z(I,J-1,K))
C
      XN=0.25*(X(I,J+1,K)+X(I,J,K)+X(I-1,J,K)+X(I-1,J+1,K))
      YN=0.25*(Y(I,J+1,K)+Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J+1,K))
      ZN=0.25*(Z(I,J+1,K)+Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J+1,K))
C
      XM=0.5*(X(I,J,K)+X(I-1,J,K))
      YM=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZM=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I,J,K+1)=FYTOP/FYBOT
      FY(I,J,K+1)=FYTOP/(FYTOP+FYBOT1)
120   CONTINUE
C
      I=1
      DO 121 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
      FY(I,J,K+1)=FYTOP/FYBOT
      FY(I,J,K+1)=FYTOP/(FYTOP+FYBOT1)
 121  CONTINUE
C
      I=NIM1
      DO 122 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I+1,J,K+1)=FYTOP/FYBOT
      FY(I+1,J,K+1)=FYTOP/(FYTOP+FYBOT1)
 122  CONTINUE
C
C00000000
C
C EW
C0000000
      I=1
CCCC
      DO 130 J=2,NJM2
      DO 130 K=2,NKM1
      XP=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J-1,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K-1)+Z(I,J,K-1))
C
      XN=0.25*(X(I,J+1,K)+X(I,J,K)+X(I,J,K-1)+X(I,J+1,K-1))
      YN=0.25*(Y(I,J+1,K)+Y(I,J,K)+Y(I,J,K-1)+Y(I,J+1,K-1))
      ZN=0.25*(Z(I,J+1,K)+Z(I,J,K)+Z(I,J,K-1)+Z(I,J+1,K-1))
C
      XM=0.5*(X(I,J,K)+X(I,J,K-1))
      YM=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZM=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I  ,J,K)=FYTOP/FYBOT
      FY(I,J,K)=FYTOP/(FYTOP+FYBOT1)
130   CONTINUE
C
      K=1
      DO 131 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
C      FY(I  ,J,K)=FYTOP/FYBOT
      FY(I,J,K)=FYTOP/(FYTOP+FYBOT1)
 131  CONTINUE
C
      K=NKM1
      DO 132 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I  ,J,K+1)=FYTOP/FYBOT
      FY(I,J,K+1)=FYTOP/(FYTOP+FYBOT1)
 132  CONTINUE
C
C      GO TO 845
C      I=NIR
CCCC
C      DO 330 J=2,NJM1
C      DO 330 K=2,NKR-1
C      XP=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J-1,K-1)+X(I,J,K-1))
C      YP=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K-1)+Y(I,J,K-1))
C      ZP=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K-1)+Z(I,J,K-1))
C
C      XN=0.25*(X(I,J+1,K)+X(I,J,K)+X(I,J,K-1)+X(I,J+1,K-1))
C      YN=0.25*(Y(I,J+1,K)+Y(I,J,K)+Y(I,J,K-1)+Y(I,J+1,K-1))
C      ZN=0.25*(Z(I,J+1,K)+Z(I,J,K)+Z(I,J,K-1)+Z(I,J+1,K-1))
C
C      XM=0.5*(X(I,J,K)+X(I,J,K-1))
C      YM=0.5*(Y(I,J,K)+Y(I,J,K-1))
C      ZM=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
C      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I  ,J,K)=FYTOP/FYBOT
C330   CONTINUE
C
C      K=1
C      DO 331 J=2,NJM1
C      XP=0.5*(X(I,J,K)+X(I,J-1,K))
C      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
C      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
C      XN=0.5*(X(I,J,K)+X(I,J+1,K))
C      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
C      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I  ,J,K)=FYTOP/FYBOT
C 331  CONTINUE
C
C     K=NKR
C      DO 332 J=2,NJM1
C      XP=0.5*(X(I,J,K)+X(I,J-1,K))
C      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
C      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
C      XN=0.5*(X(I,J,K)+X(I,J+1,K))
C      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
C      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I  ,J,K+1)=FYTOP/FYBOT
C 332  CONTINUE
C 845  CONTINUE
C
C00000000
C**********************************************************
C0000000
      I=NIM1
CCCC
      DO 140 J=2,NJM2
      DO 140 K=2,NKM1
      XP=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J-1,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K-1)+Z(I,J,K-1))
C
      XN=0.25*(X(I,J+1,K)+X(I,J,K)+X(I,J,K-1)+X(I,J+1,K-1))
      YN=0.25*(Y(I,J+1,K)+Y(I,J,K)+Y(I,J,K-1)+Y(I,J+1,K-1))
      ZN=0.25*(Z(I,J+1,K)+Z(I,J,K)+Z(I,J,K-1)+Z(I,J+1,K-1))
C
      XM=0.5*(X(I,J,K)+X(I,J,K-1))
      YM=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZM=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I+1,J,K)=FYTOP/FYBOT
      FY(I+1,J,K)=FYTOP/(FYTOP+FYBOT1)
140   CONTINUE
C
      K=1
      DO 141 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I+1,J,K)=FYTOP/FYBOT
      FY(I+1,J,K)=FYTOP/(FYTOP+FYBOT1)
 141  CONTINUE
C
      K=NKM1
      DO 142 J=2,NJM2
      XP=0.5*(X(I,J,K)+X(I,J-1,K))
      YP=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZP=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      XN=0.5*(X(I,J,K)+X(I,J+1,K))
      YN=0.5*(Y(I,J,K)+Y(I,J+1,K))
      ZN=0.5*(Z(I,J,K)+Z(I,J+1,K))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FYTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FYBOT1=SQRT( (XN-XM)**2+(YN-YM)**2+(ZN-ZM)**2)
      FYBOT=SQRT( (XN-XP)**2+(YN-YP)**2+(ZN-ZP)**2 )
C      FY(I+1,J,K+1)=FYTOP/FYBOT
      FY(I+1,J,K+1)=FYTOP/(FYTOP+FYBOT1)
 142  CONTINUE
C
C00000000
C
C NS
C
      DO 150 K=1,NK
      DO 150 I=1,NI
      FY(I,1,K)=0.
      FY(I,NJM1,K)=1.
 150  CONTINUE
C
      RETURN
      END
      SUBROUTINE GRID2
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
      NIM2=NIM1-1
      NJM2=NJM1-1
      NKM2=NKM1-1
C
C*****FZ
C
C WEST
C0000000
      I=1
CCCC
      DO 160 K=2,NKM2
      DO 160 J=2,NJM1
      XP=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J-1,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K-1)+Z(I,J,K-1))
C
      XT=0.25*(X(I,J,K+1)+X(I,J,K)+X(I,J-1,K)+X(I,J-1,K+1))
      YT=0.25*(Y(I,J,K+1)+Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K+1))
      ZT=0.25*(Z(I,J,K+1)+Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K+1))
C
      XM=0.5*(X(I,J,K)+X(I,J-1,K))
      YM=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZM=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I  ,J,K)=FZTOP/FZBOT
160   CONTINUE
C
      J=1
      DO 161 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I  ,J,K)=FZTOP/FZBOT
 161  CONTINUE
C
      J=NJM1
      DO 162 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I  ,J+1,K)=FZTOP/FZBOT
 162  CONTINUE
C
C0000000
C
C STEP
C0000000
C      GO TO 843
C      I=NIR
CCCC
C      DO 360 K=2,NKR-1
C      DO 360 J=2,NJM1
C      XP=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J-1,K-1)+X(I,J,K-1))
C      YP=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K-1)+Y(I,J,K-1))
C      ZP=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K-1)+Z(I,J,K-1))
C
C      XT=0.25*(X(I,J,K+1)+X(I,J,K)+X(I,J-1,K)+X(I,J-1,K+1))
C      YT=0.25*(Y(I,J,K+1)+Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K+1))
C      ZT=0.25*(Z(I,J,K+1)+Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K+1))
C
C      XM=0.5*(X(I,J,K)+X(I,J-1,K))
C      YM=0.5*(Y(I,J,K)+Y(I,J-1,K))
C      ZM=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
C      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
C      FZ(I  ,J,K)=FZTOP/FZBOT
C 360  CONTINUE
C
C      J=1
C      DO 361 K=2,NKR-1
C      XP=0.5*(X(I,J,K)+X(I,J,K-1))
CC      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
C      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
C      XT=0.5*(X(I,J,K)+X(I,J,K+1))
C      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
C      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
C      FZ(I  ,J,K)=FZTOP/FZBOT
C  361  CONTINUE
C
C      J=NJM1
C      DO 362 K=2,NKR-1
C      XP=0.5*(X(I,J,K)+X(I,J,K-1))
C      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
C      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
C      XT=0.5*(X(I,J,K)+X(I,J,K+1))
C      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
C      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
C      XM=X(I,J,K)
C      YM=Y(I,J,K)
C      ZM=Z(I,J,K)
C
C      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
C      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
C      FZ(I  ,J+1,K)=FZTOP/FZBOT
C 362  CONTINUE
C  843 CONTINUE
C
C0000000000
C EAST
C0000000000
      I=NIM1
CCCC
      DO 170 K=2,NKM2
      DO 170 J=2,NJM1
      XP=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J-1,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K-1)+Z(I,J,K-1))
C
      XT=0.25*(X(I,J,K+1)+X(I,J,K)+X(I,J-1,K)+X(I,J-1,K+1))
      YT=0.25*(Y(I,J,K+1)+Y(I,J,K)+Y(I,J-1,K)+Y(I,J-1,K+1))
      ZT=0.25*(Z(I,J,K+1)+Z(I,J,K)+Z(I,J-1,K)+Z(I,J-1,K+1))
C
      XM=0.5*(X(I,J,K)+X(I,J-1,K))
      YM=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZM=0.5*(Z(I,J,K)+Z(I,J-1,K))
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I+1,J,K)=FZTOP/FZBOT
170   CONTINUE
C
      J=1
      DO 171 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I+1,J,K)=FZTOP/FZBOT
 171  CONTINUE
C
      J=NJM1
      DO 172 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I+1,J+1,K)=FZTOP/FZBOT
 172  CONTINUE
C
C00000000
C
C NS
C0000000
      J=1
CCCC
      DO 180 K=2,NKM2
      DO 180 I=2,NIM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J,K-1)+Z(I,J,K-1))
C
      XT=0.25*(X(I,J,K+1)+X(I,J,K)+X(I-1,J,K)+X(I-1,J,K+1))
      YT=0.25*(Y(I,J,K+1)+Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J,K+1))
      ZT=0.25*(Z(I,J,K+1)+Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J,K+1))
C
      XM=0.5*(X(I,J,K)+X(I-1,J,K))
      YM=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZM=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I,J  ,K)=FZTOP/FZBOT
180   CONTINUE
C
      I=1
      DO 181 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I,J  ,K)=FZTOP/FZBOT
 181  CONTINUE
C
      I=NIM1
      DO 182 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I+1,J  ,K)=FZTOP/FZBOT
 182  CONTINUE
C
C0000000
C0000000
      J=NJM1
CCCC
      DO 190 K=2,NKM2
      DO 190 I=2,NIM1
      XP=0.25*(X(I,J,K)+X(I-1,J,K)+X(I-1,J,K-1)+X(I,J,K-1))
      YP=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J,K-1)+Y(I,J,K-1))
      ZP=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J,K-1)+Z(I,J,K-1))
C
      XT=0.25*(X(I,J,K+1)+X(I,J,K)+X(I-1,J,K)+X(I-1,J,K+1))
      YT=0.25*(Y(I,J,K+1)+Y(I,J,K)+Y(I-1,J,K)+Y(I-1,J,K+1))
      ZT=0.25*(Z(I,J,K+1)+Z(I,J,K)+Z(I-1,J,K)+Z(I-1,J,K+1))
C
      XM=0.5*(X(I,J,K)+X(I-1,J,K))
      YM=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZM=0.5*(Z(I,J,K)+Z(I-1,J,K))
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I,J+1,K)=FZTOP/FZBOT
190   CONTINUE
C
      I=1
      DO 191 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I,J+1,K)=FZTOP/FZBOT
 191  CONTINUE
C
      I=NIM1
      DO 192 K=2,NKM2
      XP=0.5*(X(I,J,K)+X(I,J,K-1))
      YP=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZP=0.5*(Z(I,J,K)+Z(I,J,K-1))
C
      XT=0.5*(X(I,J,K)+X(I,J,K+1))
      YT=0.5*(Y(I,J,K)+Y(I,J,K+1))
      ZT=0.5*(Z(I,J,K)+Z(I,J,K+1))
C
      XM=X(I,J,K)
      YM=Y(I,J,K)
      ZM=Z(I,J,K)
C
      FZTOP=SQRT( (XM-XP)**2+(YM-YP)**2+(ZM-ZP)**2 )
      FZBOT=SQRT( (XT-XP)**2+(YT-YP)**2+(ZT-ZP)**2 )
      FZ(I+1,J+1,K)=FZTOP/FZBOT
 192  CONTINUE
C
C0000000
C TB
C
      DO 200 I=1,NI
      DO 200 J=1,NJ
      FZ(I,J,1)=0.
      FZ(I,J,NKM1)=1.
 200  CONTINUE
C      DO 400 I=1,NIR
C     DO 400 J=1,NJ
C      FZ(I,J,NKR)=0.
C 400  CONTINUE
C
C     CALCULATE THE VOLUME
C
      DO 210 I=2,NIM1
      DO 210 J=2,NJM1
      DO 210 K=2,NKM1
      X1=X(I-1,J-1,K-1)
      X2=X(I-1,J,K-1)
      X3=X(I,J,K-1)
      X4=X(I,J-1,K-1)
      X5=X(I-1,J-1,K)
      X6=X(I-1,J,K)
      X7=X(I,J,K)
      X8=X(I,J-1,K)
C
      Y1=Y(I-1,J-1,K-1)
      Y2=Y(I-1,J,K-1)
      Y3=Y(I,J,K-1)
      Y4=Y(I,J-1,K-1)
      Y5=Y(I-1,J-1,K)
      Y6=Y(I-1,J,K)
      Y7=Y(I,J,K)
      Y8=Y(I,J-1,K)
C
      Z1=Z(I-1,J-1,K-1)
      Z2=Z(I-1,J,K-1)
      Z3=Z(I,J,K-1)
      Z4=Z(I,J-1,K-1)
      Z5=Z(I-1,J-1,K)
      Z6=Z(I-1,J,K)
      Z7=Z(I,J,K)
      Z8=Z(I,J-1,K)
C
C
      CALL CVOL(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     1          X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,VOLUME)
      VOL(I,J,K)=VOLUME
 210  CONTINUE
C
C     MODIFY VOLUME AND METRIC ON BOUNDARY
C
C  WEST
C
      I=1
      DO 220 J=2,NJM1
      DO 220 K=2,NKM1
      VOL(I,J,K)=0.5*VOL(I+1,J,K)
C
      XE=XC(I+1,J,K)
      YE=YC(I+1,J,K)
      ZE=ZC(I+1,J,K)
      XW=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J,K-1)+X(I,J-1,K-1))
      YW=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J,K-1)+Y(I,J-1,K-1))
      ZW=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J,K-1)+Z(I,J-1,K-1))
      XN=0.5*(X(I,J,K)+X(I,J,K-1))
      YN=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZN=0.5*(Z(I,J,K)+Z(I,J,K-1))
      XS=0.5*(X(I,J-1,K)+X(I,J-1,K-1))
      YS=0.5*(Y(I,J-1,K)+Y(I,J-1,K-1))
      ZS=0.5*(Z(I,J-1,K)+Z(I,J-1,K-1))
      XT=0.5*(X(I,J,K)+X(I,J-1,K))
      YT=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZT=0.5*(Z(I,J,K)+Z(I,J-1,K))
      XB=0.5*(X(I,J,K-1)+X(I,J-1,K-1))
      YB=0.5*(Y(I,J,K-1)+Y(I,J-1,K-1))
      ZB=0.5*(Z(I,J,K-1)+Z(I,J-1,K-1))
C
      XX(I,J,K)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
      XY(I,J,K)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
      XZ(I,J,K)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
      YX(I,J,K)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
      YY(I,J,K)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
      YZ(I,J,K)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
      ZX(I,J,K)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
      ZY(I,J,K)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
      ZZ(I,J,K)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
 220  CONTINUE
C
C
C      GO TO 424
C      I=NIR
C      DO 420 J=2,NJM1
C      DO 420 K=2,NKR
C      VOL(I,J,K)=0.5*VOL(I+1,J,K)
C
C      XE=XC(I+1,J,K)
C      YE=YC(I+1,J,K)
C      ZE=ZC(I+1,J,K)
C      XW=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J,K-1)+X(I,J-1,K-1))
C      YW=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J,K-1)+Y(I,J-1,K-1))
CC      ZW=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J,K-1)+Z(I,J-1,K-1))
C      XN=0.5*(X(I,J,K)+X(I,J,K-1))
C      YN=0.5*(Y(I,J,K)+Y(I,J,K-1))
C      ZN=0.5*(Z(I,J,K)+Z(I,J,K-1))
C      XS=0.5*(X(I,J-1,K)+X(I,J-1,K-1))
C      YS=0.5*(Y(I,J-1,K)+Y(I,J-1,K-1))
C      ZS=0.5*(Z(I,J-1,K)+Z(I,J-1,K-1))
C      XT=0.5*(X(I,J,K)+X(I,J-1,K))
C      YT=0.5*(Y(I,J,K)+Y(I,J-1,K))
C      ZT=0.5*(Z(I,J,K)+Z(I,J-1,K))
C      XB=0.5*(X(I,J,K-1)+X(I,J-1,K-1))
C      YB=0.5*(Y(I,J,K-1)+Y(I,J-1,K-1))
C      ZB=0.5*(Z(I,J,K-1)+Z(I,J-1,K-1))
C
C      XX(I,J,K)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
C      XY(I,J,K)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
C      XZ(I,J,K)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
C      YX(I,J,K)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
C      YY(I,J,K)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
C      YZ(I,J,K)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
C      ZX(I,J,K)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
C      ZY(I,J,K)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
C      ZZ(I,J,K)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
C 420  CONTINUE
C 424  CONTINUE
C
C EAST
C
      I=NIM1
      DO 221 J=2,NJM1
      DO 221 K=2,NKM1
      VOL(I+1,J,K)=0.5*VOL(I,J,K)
C
      XW=XC(I,J,K)
      YW=YC(I,J,K)
      ZW=ZC(I,J,K)
      XE=0.25*(X(I,J,K)+X(I,J-1,K)+X(I,J,K-1)+X(I,J-1,K-1))
      YE=0.25*(Y(I,J,K)+Y(I,J-1,K)+Y(I,J,K-1)+Y(I,J-1,K-1))
      ZE=0.25*(Z(I,J,K)+Z(I,J-1,K)+Z(I,J,K-1)+Z(I,J-1,K-1))
      XN=0.5*(X(I,J,K)+X(I,J,K-1))
      YN=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZN=0.5*(Z(I,J,K)+Z(I,J,K-1))
      XS=0.5*(X(I,J-1,K)+X(I,J-1,K-1))
      YS=0.5*(Y(I,J-1,K)+Y(I,J-1,K-1))
      ZS=0.5*(Z(I,J-1,K)+Z(I,J-1,K-1))
      XT=0.5*(X(I,J,K)+X(I,J-1,K))
      YT=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZT=0.5*(Z(I,J,K)+Z(I,J-1,K))
      XB=0.5*(X(I,J,K-1)+X(I,J-1,K-1))
      YB=0.5*(Y(I,J,K-1)+Y(I,J-1,K-1))
      ZB=0.5*(Z(I,J,K-1)+Z(I,J-1,K-1))
C
      XX(I+1,J,K)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
      XY(I+1,J,K)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
      XZ(I+1,J,K)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
      YX(I+1,J,K)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
      YY(I+1,J,K)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
      YZ(I+1,J,K)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
      ZX(I+1,J,K)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
      ZY(I+1,J,K)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
      ZZ(I+1,J,K)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
 221  CONTINUE
C
C SOUTH
C
      J=1
      DO 230 K=2,NKM1
      DO 230 I=2,NIM1
      VOL(I,J,K)=0.5*VOL(I,J+1,K)
C
      XE=0.5*(X(I,J,K)+X(I,J,K-1))
      YE=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZE=0.5*(Z(I,J,K)+Z(I,J,K-1))
      XW=0.5*(X(I-1,J,K)+X(I-1,J,K-1))
      YW=0.5*(Y(I-1,J,K)+Y(I-1,J,K-1))
      ZW=0.5*(Z(I-1,J,K)+Z(I-1,J,K-1))
      XN=XC(I,J+1,K)
      YN=YC(I,J+1,K)
      ZN=ZC(I,J+1,K)
      XS=0.25*(X(I,J,K)+X(I-1,J,K)+X(I,J,K-1)+X(I-1,J,K-1))
      YS=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I,J,K-1)+Y(I-1,J,K-1))
      ZS=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I,J,K-1)+Z(I-1,J,K-1))
      XT=0.5*(X(I,J,K)+X(I-1,J,K))
      YT=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZT=0.5*(Z(I,J,K)+Z(I-1,J,K))
      XB=0.5*(X(I,J,K-1)+X(I-1,J,K-1))
      YB=0.5*(Y(I,J,K-1)+Y(I-1,J,K-1))
      ZB=0.5*(Z(I,J,K-1)+Z(I-1,J,K-1))
C
      XX(I,J,K)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
      XY(I,J,K)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
      XZ(I,J,K)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
      YX(I,J,K)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
      YY(I,J,K)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
      YZ(I,J,K)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
      ZX(I,J,K)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
      ZY(I,J,K)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
      ZZ(I,J,K)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
 230  CONTINUE
C
C NORTH
C
      J=NJM1
      DO 231 K=2,NKM1
      DO 231 I=2,NIM1
      VOL(I,J+1,K)=0.5*VOL(I,J,K)
C
      XE=0.5*(X(I,J,K)+X(I,J,K-1))
      YE=0.5*(Y(I,J,K)+Y(I,J,K-1))
      ZE=0.5*(Z(I,J,K)+Z(I,J,K-1))
      XW=0.5*(X(I-1,J,K)+X(I-1,J,K-1))
      YW=0.5*(Y(I-1,J,K)+Y(I-1,J,K-1))
      ZW=0.5*(Z(I-1,J,K)+Z(I-1,J,K-1))
      XS=XC(I,J,K)
      YS=YC(I,J,K)
      ZS=ZC(I,J,K)
      XN=0.25*(X(I,J,K)+X(I-1,J,K)+X(I,J,K-1)+X(I-1,J,K-1))
      YN=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I,J,K-1)+Y(I-1,J,K-1))
      ZN=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I,J,K-1)+Z(I-1,J,K-1))
      XT=0.5*(X(I,J,K)+X(I-1,J,K))
      YT=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZT=0.5*(Z(I,J,K)+Z(I-1,J,K))
      XB=0.5*(X(I,J,K-1)+X(I-1,J,K-1))
      YB=0.5*(Y(I,J,K-1)+Y(I-1,J,K-1))
      ZB=0.5*(Z(I,J,K-1)+Z(I-1,J,K-1))
C
      XX(I,J+1,K)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
      XY(I,J+1,K)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
      XZ(I,J+1,K)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
      YX(I,J+1,K)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
      YY(I,J+1,K)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
      YZ(I,J+1,K)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
      ZX(I,J+1,K)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
      ZY(I,J+1,K)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
      ZZ(I,J+1,K)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
 231  CONTINUE
C
C BOTTOM
C
      K=1
      DO 240 I=2,NIM1
      DO 240 J=2,NJM1
      VOL(I,J,K)=0.5*VOL(I,J,K+1)
C
      XE=0.5*(X(I,J,K)+X(I,J-1,K))
      YE=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZE=0.5*(Z(I,J,K)+Z(I,J-1,K))
      XW=0.5*(X(I-1,J,K)+X(I-1,J-1,K))
      YW=0.5*(Y(I-1,J,K)+Y(I-1,J-1,K))
      ZW=0.5*(Z(I-1,J,K)+Z(I-1,J-1,K))
      XN=0.5*(X(I,J,K)+X(I-1,J,K))
      YN=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZN=0.5*(Z(I,J,K)+Z(I-1,J,K))
      XS=0.5*(X(I,J-1,K)+X(I-1,J-1,K))
      YS=0.5*(Y(I,J-1,K)+Y(I-1,J-1,K))
      ZS=0.5*(Z(I,J-1,K)+Z(I-1,J-1,K))
      XT=XC(I,J,K+1)
      YT=YC(I,J,K+1)
      ZT=ZC(I,J,K+1)
      XB=0.25*(X(I,J,K)+X(I-1,J,K)+X(I,J-1,K)+X(I-1,J-1,K))
      YB=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I,J-1,K)+Y(I-1,J-1,K))
      ZB=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I,J-1,K)+Z(I-1,J-1,K))
C
      XX(I,J,K)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
      XY(I,J,K)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
      XZ(I,J,K)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
      YX(I,J,K)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
      YY(I,J,K)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
      YZ(I,J,K)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
      ZX(I,J,K)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
      ZY(I,J,K)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
      ZZ(I,J,K)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
 240  CONTINUE
C
C TOP
C
      K=NKM1
      DO 241 I=2,NIM1
      DO 241 J=2,NJM1
      VOL(I,J,K+1)=0.5*VOL(I,J,K)
C
      XE=0.5*(X(I,J,K)+X(I,J-1,K))
      YE=0.5*(Y(I,J,K)+Y(I,J-1,K))
      ZE=0.5*(Z(I,J,K)+Z(I,J-1,K))
      XW=0.5*(X(I-1,J,K)+X(I-1,J-1,K))
      YW=0.5*(Y(I-1,J,K)+Y(I-1,J-1,K))
      ZW=0.5*(Z(I-1,J,K)+Z(I-1,J-1,K))
      XN=0.5*(X(I,J,K)+X(I-1,J,K))
      YN=0.5*(Y(I,J,K)+Y(I-1,J,K))
      ZN=0.5*(Z(I,J,K)+Z(I-1,J,K))
      XS=0.5*(X(I,J-1,K)+X(I-1,J-1,K))
      YS=0.5*(Y(I,J-1,K)+Y(I-1,J-1,K))
      ZS=0.5*(Z(I,J-1,K)+Z(I-1,J-1,K))
      XB=XC(I,J,K)
      YB=YC(I,J,K)
      ZB=ZC(I,J,K)
      XT=0.25*(X(I,J,K)+X(I-1,J,K)+X(I,J-1,K)+X(I-1,J-1,K))
      YT=0.25*(Y(I,J,K)+Y(I-1,J,K)+Y(I,J-1,K)+Y(I-1,J-1,K))
      ZT=0.25*(Z(I,J,K)+Z(I-1,J,K)+Z(I,J-1,K)+Z(I-1,J-1,K))
C
      XX(I,J,K+1)=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
      XY(I,J,K+1)=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
      XZ(I,J,K+1)=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
      YX(I,J,K+1)=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
      YY(I,J,K+1)=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
      YZ(I,J,K+1)=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
      ZX(I,J,K+1)=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
      ZY(I,J,K+1)=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
      ZZ(I,J,K+1)=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
 241  CONTINUE
C
C ----------------INTERNODAL DISTANCES --------------------
      DO 351 K=1,NK
      DO 351 J=1,NJ
      DO 351 I=1,NIM1
      DSI(I,J,K)=SQRT( (XC(I+1,J,K)-XC(I,J,K))**2+
     1  (YC(I+1,J,K)-YC(I,J,K))**2+(ZC(I+1,J,K)-ZC(I,J,K))**2 )
  351 CONTINUE
C
      DO 352 K=1,NK
      DO 352 J=1,NJM1
      DO 352 I=1,NI
      DSJ(I,J,K)=SQRT( (XC(I,J+1,K)-XC(I,J,K))**2+
     1     (YC(I,J+1,K)-YC(I,J,K))**2+(ZC(I,J+1,K)-ZC(I,J,K))**2 )
  352 CONTINUE
C
      DO 353 K=1,NKM1
      DO 353 J=1,NJ
      DO 353 I=1,NI
      DSK(I,J,K)=SQRT( (XC(I,J,K+1)-XC(I,J,K))**2+
     1   (YC(I,J,K+1)-YC(I,J,K))**2+(ZC(I,J,K+1)-ZC(I,J,K))**2 )
  353 CONTINUE     
C
      RETURN
      END
      SUBROUTINE FLORAT(UREF,VREF,WREF,XA,YA,ZA,XB,YB,ZB,
     1 XC,YC,ZC,XD,YD,ZD,FLOW)
C
C      CALCULATE THE MASS FLOW RATE
C
      FLOW= UREF*( (YC-YA)*(ZD-ZB)-(ZC-ZA)*(YD-YB) )
     1     +VREF*( (ZC-ZA)*(XD-XB)-(XC-XA)*(ZD-ZB) )
     1     +WREF*( (XC-XA)*(YD-YB)-(YC-YA)*(XD-XB) )
      FLOW=0.5*FLOW
C
      RETURN
      END
C
      SUBROUTINE AREA(XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     1 AREAX,AREAY,AREAZ)
C
C     CALCULATE AREA PROGECTION IN X,Y,Z DIRECTION
C
      AREAX=0.5*( (YC-YA)*(ZD-ZB)-(ZC-ZA)*(YD-YB) )
      AREAY=0.5*( (ZC-ZA)*(XD-XB)-(XC-XA)*(ZD-ZB) )
      AREAZ=0.5*( (XC-XA)*(YD-YB)-(YC-YA)*(XD-XB) )
C
      RETURN
      END
C
C
C
      SUBROUTINE CVOL(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     1               X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8,Z8,VOLUME)
C
C     CALCULATE THE VOLUME OF A CELL
C
      CALL AREA(X1,Y1,Z1,X5,Y5,Z5,X8,Y8,Z8,X4,Y4,Z4,
     1 A1584X,A1584Y,A1584Z)
      CALL AREA(X1,Y1,Z1,X2,Y2,Z2,X6,Y6,Z6,X5,Y5,Z5,
     1 A1265X,A1265Y,A1265Z)
      CALL AREA(X1,Y1,Z1,X4,Y4,Z4,X3,Y3,Z3,X2,Y2,Z2,
     1 A1432X,A1432Y,A1432Z)
C
      VOLUME=(X7-X1)*(A1584X+A1265X+A1432X)
     1      +(Y7-Y1)*(A1584Y+A1265Y+A1432Y)
     1      +(Z7-Z1)*(A1584Z+A1265Z+A1432Z)
      VOLUME=VOLUME/3.
      RETURN
      END
C
      FUNCTION APFUN(ISCHEM,PEC)
      GO TO (10,20,30,40),ISCHEM
 10   APFUN=1.-0.5*ABS(PEC)
      GO TO 100
 20   APFUN=1.
      GO TO 100
 30   APFUN=AMAX1(0.,1.-0.5*ABS(PEC))
      GO TO 100
 40   APFUN=AMAX1(0.,(1.-0.1*ABS(PEC))**5)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE LISOLV(ISTART,JSTART,KSTART,NI,NJ,NK,PHI)
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 ,NM=150)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      DIMENSION  PHI(NX,NY,NZ),A(NM,NM),C(NM,NM)
CSEZ***   RED-BLACK ADI SOLUTION FORMULATION   ***
CSEZ***     SERGE ZARANTONELLO  7/18/88        ***
CSEZ***         AMDAHL CORPORATION             ***
C
      NIM1 = NI-1
      NJM1 = NJ-1
      NKM1 = NK-1
      JSTM1= JSTART-1
      ISTM1= ISTART-1
      KSTM1= KSTART-1
C
CSEZ***  ITERATIVE ADI LOOP.  N-S/E-W SWEEPS  ***
C
      DO 1000 K=KSTART,NKM1
C
CSEZ***  INITIALIZE FOR COLUMN SWEEPS  ***
C
      DO 100 I=ISTART,NIM1
       C(I,JSTM1)=PHI(I,JSTM1,K)
       A(I,JSTM1)=0.0
100   CONTINUE
C
      IBEG = ISTART
C
CSEZ***  BLACK N-S COLUMNS FIRST,  RED N-S COLUMNS NEXT  ***
C
150   CONTINUE
C
CSEZ***  SOLVE TRIDIAGONAL SYSTEMS IN PARALLEL   ***
C
      DO 120 J=JSTART,NJM1
      DO 110 I=IBEG,NIM1,2
       A(I,J)  =  AN(I,J,K)
       C(I,J)  =  AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)
     1  +AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)+SU(I,J,K)
       TERM    =  1./(AP(I,J,K)-AS(I,J,K)*A(I,J-1))
       A(I,J)  =  A(I,J)*TERM
       C(I,J)  = (C(I,J)+AS(I,J,K)*C(I,J-1))*TERM
  110 CONTINUE
  120 CONTINUE
C
      DO 140 JJ= JSTART,NJM1
       J=NJ+JSTM1-JJ
       DO 130 I = IBEG,NIM1,2
        PHI(I,J,K)=A(I,J)*PHI(I,J+1,K)+C(I,J)
  130 CONTINUE
  140 CONTINUE
C
      IF(IBEG.EQ.ISTART) THEN
       IBEG=ISTART+1
       GO TO 150
      ENDIF
C
      JBEG = JSTART
C
CSEZ***  INITIALIZE FOR ROW SWEEPS  ***
C
      DO 200 J=JSTART,NJM1
       C(ISTM1,J)=PHI(ISTM1,J,K)
       A(ISTM1,J)=0.0
200   CONTINUE
C
C
CSEZ***  BLACK W-E ROWS FIRST, RED W-E ROWS NEXT  ***
C
250   CONTINUE
C
CSEZ***  SOLVE TRIDIAGONAL SYSTEMS IN PARALLEL  ***
C
      DO 220 I=ISTART,NIM1
      DO 210 J=JBEG,NJM1,2
       A(I,J)  =  AE(I,J,K)
       C(I,J)  =  AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)
     1  +AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)+SU(I,J,K)
       TERM    =  1./(AP(I,J,K)-AW(I,J,K)*A(I-1,J))
       A(I,J)  =  A(I,J)*TERM
       C(I,J)  = (C(I,J)+AW(I,J,K)*C(I-1,J))*TERM
  210 CONTINUE
  220 CONTINUE
C
      DO 240 II= ISTART,NIM1
       I=NI+ISTM1-II
       DO 230 J = JBEG,NJM1,2
        PHI(I,J,K)=A(I,J)*PHI(I+1,J,K)+C(I,J)
  230 CONTINUE
  240 CONTINUE
C
      IF(JBEG.EQ.JSTART) THEN
       JBEG=JSTART+1
       GO TO 250
      ENDIF
 1000 CONTINUE
C
CSEZ***  END OF SWEEP  ***
C
C
C
CSEZ***  ITERATIVE ADI LOOP.  T-B/E-W SWEEPS  ***
C
C      GO TO 2001
      DO 2000 J=JSTART,NJM1
C
CSEZ***  INITIALIZE FOR COLUMN SWEEPS  ***
C
      DO 300 I=ISTART,NIM1
       C(I,KSTM1)=PHI(I,J,KSTM1)
       A(I,KSTM1)=0.0
 300   CONTINUE
C
      IBEG = ISTART
C
CSEZ***  BLACK T-B COLUMNS FIRST,  RED T-B COLUMNS NEXT  ***
C
 350   CONTINUE
C
CSEZ***  SOLVE TRIDIAGONAL SYSTEMS IN PARALLEL   ***
C
      DO 320 K=KSTART,NKM1
      DO 310 I=IBEG,NIM1,2
       A(I,K)  =  AT(I,J,K)
       C(I,K)  =  AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)
     1  +AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)+SU(I,J,K)
       TERM    =  1./(AP(I,J,K)-AB(I,J,K)*A(I,K-1))
       A(I,K)  =  A(I,K)*TERM
       C(I,K)  = (C(I,K)+AB(I,J,K)*C(I,K-1))*TERM
  310 CONTINUE
  320 CONTINUE
C
      DO 340 KK= KSTART,NKM1
       K=NK+KSTM1-KK
       DO 330 I = IBEG,NIM1,2
        PHI(I,J,K)=A(I,K)*PHI(I,J,K+1)+C(I,K)
  330 CONTINUE
  340 CONTINUE
C
      IF(IBEG.EQ.ISTART) THEN
       IBEG=ISTART+1
       GO TO 350
      ENDIF
C
      KBEG = KSTART
C
CSEZ***  INITIALIZE FOR ROW SWEEPS  ***
C
      DO 400 K=KSTART,NKM1
       C(ISTM1,K)=PHI(ISTM1,J,K)
       A(ISTM1,K)=0.0
 400   CONTINUE
C
C
CSEZ***  BLACK W-E ROWS FIRST, RED W-E ROWS NEXT  ***
C
 450   CONTINUE
C
CSEZ***  SOLVE TRIDIAGONAL SYSTEMS IN PARALLEL  ***
C
      DO 420 I=ISTART,NIM1
      DO 410 K=KBEG,NKM1,2
       A(I,K)  =  AE(I,J,K)
       C(I,K)  =  AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)
     1  +AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)+SU(I,J,K)
       TERM    =  1./(AP(I,J,K)-AW(I,J,K)*A(I-1,K))
       A(I,K)  =  A(I,K)*TERM
       C(I,K)  = (C(I,K)+AW(I,J,K)*C(I-1,K))*TERM
  410 CONTINUE
  420 CONTINUE
C
      DO 440 II= ISTART,NIM1
       I=NI+ISTM1-II
       DO 430 K = KBEG,NKM1,2
        PHI(I,J,K)=A(I,K)*PHI(I+1,J,K)+C(I,K)
  430 CONTINUE
  440 CONTINUE
C
      IF(KBEG.EQ.KSTART) THEN
       KBEG=KSTART+1
       GO TO 450
      ENDIF
 2000 CONTINUE
 2001 CONTINUE
C
CSEZ***  END OF SWEEP  ***
C
C      GO TO 4000
C
CSEZ***  ITERATIVE ADI LOOP.  N-S/T-B SWEEPS  ***
C
      DO 3000 I=ISTART,NIM1
C
CSEZ***  INITIALIZE FOR COLUMN SWEEPS  ***
C
      DO 500 K=KSTART,NKM1
       C(JSTM1,K)=PHI(I,JSTM1,K)
       A(JSTM1,K)=0.0
 500   CONTINUE
C
      KBEG = KSTART
C
CSEZ***  BLACK N-S COLUMNS FIRST,  RED N-S COLUMNS NEXT  ***
C
 550   CONTINUE
C
CSEZ***  SOLVE TRIDIAGONAL SYSTEMS IN PARALLEL   ***
C
      DO 520 J=JSTART,NJM1
      DO 510 K=KBEG,NKM1,2
       A(J,K)  =  AN(I,J,K)
       C(J,K)  =  AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)
     1  +AT(I,J,K)*PHI(I,J,K+1)+AB(I,J,K)*PHI(I,J,K-1)+SU(I,J,K)
       TERM    =  1./(AP(I,J,K)-AS(I,J,K)*A(J-1,K))
       A(J,K)  =  A(J,K)*TERM
       C(J,K)  = (C(J,K)+AS(I,J,K)*C(J-1,K))*TERM
  510 CONTINUE
  520 CONTINUE
C
      DO 540 JJ= JSTART,NJM1
       J=NJ+JSTM1-JJ
       DO 530 K = KBEG,NKM1,2
        PHI(I,J,K)=A(J,K)*PHI(I,J+1,K)+C(J,K)
  530 CONTINUE
  540 CONTINUE
C
      IF(KBEG.EQ.KSTART) THEN
       KBEG=KSTART+1
       GO TO 550
      ENDIF
C
      JBEG = JSTART
C
CSEZ***  INITIALIZE FOR ROW SWEEPS  ***
C
      DO 600 J=JSTART,NJM1
       C(J,KSTM1)=PHI(I,J,KSTM1)
       A(J,KSTM1)=0.0
 600   CONTINUE
C
C
CSEZ***  BLACK T-B ROWS FIRST, RED T-B ROWS NEXT  ***
C
 650   CONTINUE
C
CSEZ***  SOLVE TRIDIAGONAL SYSTEMS IN PARALLEL  ***
C
      DO 620 K=KSTART,NKM1
      DO 610 J=JBEG,NJM1,2
       A(J,K)  =  AT(I,J,K)
       C(J,K)  =  AN(I,J,K)*PHI(I,J+1,K)+AS(I,J,K)*PHI(I,J-1,K)
     1  +AE(I,J,K)*PHI(I+1,J,K)+AW(I,J,K)*PHI(I-1,J,K)+SU(I,J,K)
       TERM    =  1./(AP(I,J,K)-AB(I,J,K)*A(J,K-1))
       A(J,K)  =  A(J,K)*TERM
       C(J,K)  = (C(J,K)+AB(I,J,K)*C(J,K-1))*TERM
  610 CONTINUE
  620 CONTINUE
C
      DO 640 KK= KSTART,NKM1
       K=NK+KSTM1-KK
      DO 630 J = JBEG,NJM1,2
        PHI(I,J,K)=A(J,K)*PHI(I,J,K+1)+C(J,K)
  630 CONTINUE
  640 CONTINUE
C
      IF(JBEG.EQ.JSTART) THEN
       JBEG=JSTART+1
       GO TO 650
      ENDIF
 3000 CONTINUE
C
CSEZ***  END OF SWEEP  ***
C
 4000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE FLORAX(UREF,VREF,WREF,X1,Y1,Z1,X2,Y2,Z2,
     1 X3,Y3,Z3,X4,Y4,Z4,FLOW)
C
C      CALCULATE THE MASS FLOW RATE IN X-DIRECTION
C
       XN=0.5*(X1+X4)
       XS=0.5*(X2+X3)
       XT=0.5*(X1+X2)
       XB=0.5*(X3+X4)
C
       YN=0.5*(Y1+Y4)
       YS=0.5*(Y2+Y3)
       YT=0.5*(Y1+Y2)
       YB=0.5*(Y3+Y4)
C
       ZN=0.5*(Z1+Z4)
       ZS=0.5*(Z2+Z3)
       ZT=0.5*(Z1+Z2)
       ZB=0.5*(Z3+Z4)
C
       AREAX=(YN-YS)*(ZT-ZB)-(YT-YB)*(ZN-ZS)
       AREAY=(ZN-ZS)*(XT-XB)-(ZT-ZB)*(XN-XS)
       AREAZ=(XN-XS)*(YT-YB)-(XT-XB)*(YN-YS)
C
       FLOW=UREF*AREAX+VREF*AREAY+WREF*AREAZ
C
       RETURN
       END
C
      SUBROUTINE FLORAY(UREF,VREF,WREF,X1,Y1,Z1,X2,Y2,Z2,
     1 X3,Y3,Z3,X4,Y4,Z4,FLOW)
C
C      CALCULATE THE MASS FLOW RATE IN X-DIRECTION
C
       XT=0.5*(X1+X2)
       XB=0.5*(X3+X4)
       XE=0.5*(X1+X4)
       XW=0.5*(X2+X3)
C
       YT=0.5*(Y1+Y2)
       YB=0.5*(Y3+Y4)
       YE=0.5*(Y1+Y4)
       YW=0.5*(Y2+Y3)
C
       ZT=0.5*(Z1+Z2)
       ZB=0.5*(Z3+Z4)
       ZE=0.5*(Z1+Z4)
       ZW=0.5*(Z2+Z3)
C
       AREAX=(YT-YB)*(ZE-ZW)-(YE-YW)*(ZT-ZB)
       AREAY=(ZT-ZB)*(XE-XW)-(ZE-ZW)*(XT-XB)
       AREAZ=(XT-XB)*(YE-YW)-(XE-XW)*(YT-YB)
C
       FLOW=UREF*AREAX+VREF*AREAY+WREF*AREAZ
C
       RETURN
       END
      SUBROUTINE FLORAZ(UREF,VREF,WREF,X1,Y1,Z1,X2,Y2,Z2,
     1 X3,Y3,Z3,X4,Y4,Z4,FLOW)
C
C      CALCULATE THE MASS FLOW RATE IN X-DIRECTION
C
       XE=0.5*(X1+X4)
       XW=0.5*(X2+X3)
       XN=0.5*(X1+X2)
       XS=0.5*(X3+X4)
C
       YE=0.5*(Y1+Y4)
       YW=0.5*(Y2+Y3)
       YN=0.5*(Y1+Y2)
       YS=0.5*(Y3+Y4)
C
       ZE=0.5*(Z1+Z4)
       ZW=0.5*(Z2+Z3)
       ZN=0.5*(Z1+Z2)
       ZS=0.5*(Z3+Z4)
C
       AREAX=(YE-YW)*(ZN-ZS)-(YN-YS)*(ZE-ZW)
       AREAY=(ZE-ZW)*(XN-XS)-(ZN-ZS)*(XE-XW)
       AREAZ=(XE-XW)*(YN-YS)-(XN-XS)*(YE-YW)
C
       FLOW=UREF*AREAX+VREF*AREAY+WREF*AREAZ
C
       RETURN
       END
      SUBROUTINE PROPS
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
       JED= (LWS*NJSL + (1-LWS)*2)*(2-INWM)    + (INWM-1)*2
       NJED=(LWN*NJNH + (1-LWN)*NJM1)*(2-INWM) + (INWM-1)*NJM1
       KED=(LWB*NKBL + (1-LWB)*2)*(2-INWM)     + (INWM-1)*2
       NKED=(LWT*NKTH + (1-LWT)*NKM1)*(2-INWM) + (INWM-1)*NKM1
C
      NJRSL=NJR+INWP
      NJRSL1=NJRSL+1
      NJRNL=NJRN-INWP+1
      NJRNL1=NJRNL-1
      NIEL=NIR+INWP
      NIEL1=NIEL+1
      NMWL=NIRMS-INWP+1
      NMWL1=NMWL-1
      NMEL=NIRMF+INWP
      NMEL=NI-1
      NMEL1=NMEL+1
      NIWL=NI-INWP-2
      NIWL1=NIWL-1
C
C         MAIN FLOW REGION 
C

      DO 10 K=KED,NKED
      DO 10 J=JED,NJED
      DO 10 I=2,NIM1
      RT=DEN(I,J,K)*(TKE(I,J,K)**2)/(TED(I,J,K)+TINY)/VISCOS
C*-------- NLEVM ------------
      IF (INLEVM.EQ.0) THEN
	DMU=(FLOAT(IEVM)*3.4/((1+.02*RT)**2))
     1   +(FLOAT(1-IEVM)*4./((1.+0.01*RT)*(1.+0.01*RT)))
      CMUL=FLOAT(2-INWM)*CMU+FLOAT(INWM-1)*CMU*EXP(-DMU)
	ELSE


!******************	NLEVM3 REFINEMENT

      FS=(AMIN1(AMAX1(((AMIN1(1.0,RT/50.0))**2)*PK1(I,J,K)/(TED(I,J,K)
     1+TINY),0.0)/0.75,1.0))**2
      FMUI=1.0-EXP(-SQRT(RT/90.0)-(RT/400.0)**2)
	CMUL=CMU1(I,J,K)*(FMUI*FS+(1.0-FS))	

!******************	

	ENDIF

C*---------------------------
C
      VISNEW=VISCOS+CMUL*DEN(I,J,K)*TKE(I,J,K)*TKE(I,J,K)/
     1           (TED(I,J,K)+TINY)
      VIS(I,J,K)=URFVIS*VISNEW+(1.-URFVIS)*VIS(I,J,K)
 10   CONTINUE

C
C   NEAR-WALL REGIONS WHEN THE 1-EQN NEAR-WALL MODEL IS EMPLOYED
C
C----------------SOUTH WALL --------------------------------------
      DO 50 K=2,NKM1
      DO 50 J=2,JED-1
      DO 50 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      FMU=(1.-EXP(-0.0198*XLST))*(1.+5.29/XLST)
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
      VISNEW=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
      VIS(I,J,K)=URFVIS*VISNEW+(1.-URFVIS)*VIS(I,J,K)
   50 CONTINUE
C--------------- NORTH WALL -------------------------------------
      DO 60 K=2,NKM1
      DO 60 J=NJED+1,NJM1
      DO 60 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
      VISNEW=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
      VIS(I,J,K)=URFVIS*VISNEW+(1.-URFVIS)*VIS(I,J,K)
   60 CONTINUE
C--------------- TOP WALL --------------------------------------
      DO 70 K=NKED+1, NKM1
      DO 70 J=2,NJM1
      DO 70 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
      VISNEW=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
      VIS(I,J,K)=URFVIS*VISNEW+(1.-URFVIS)*VIS(I,J,K)
   70 CONTINUE
C--------------- BOTTOM WALL -----------------------------------
      DO 80 K=2,KED-1
      DO 80 J=2,NJM1
      DO 80 I=2,NIM1
      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
      VISNEW=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
      VIS(I,J,K)=URFVIS*VISNEW+(1.-URFVIS)*VIS(I,J,K)
   80 CONTINUE
C
      DO 20 K=1,NK
      DO 20 I=1,NI
      VIS(I,1 ,K)=VISCOS
      VIS(I,NJ,K)=VISCOS
 20   CONTINUE
C
      DO 30 I=1,NI
      DO 30 J=1,NJ
      VIS(I,J,1 )=VIS(I,J,2)
C*      VIS(I,J,1 )=VISCOS
      VIS(I,J,NK)=VISCOS
 30   CONTINUE
C      IF(INWM.EQ.1)THEN
C----------- RIBS  ---------------------------------------------
C      DO 21 K=KED,NKED
C
C      DO 22 J=JED,NJRSL-1
C      DO 23 I=2,NIEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
C      VISTEP=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
C      VIS(I,J,K)=VISTEP
C   23 CONTINUE
C      DO 24 I=NMWL+1,NMEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
C      VISTEP=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
C      VIS(I,J,K)=VISTEP
C   24 CONTINUE
C      DO 25 I=NIWL+1,NI
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
C      VISTEP=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
C      VIS(I,J,K)=VISTEP
C   25 CONTINUE
C   22 CONTINUE
C
C      DO 26 J=NJRNL+1,NJED
C      DO 27 I=2,NIEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
C      VISTEP=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
C      VIS(I,J,K)=VISTEP
C   27 CONTINUE
C      DO 28 I=NMWL+1,NMEL-1
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
C      VISTEP=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
C      VIS(I,J,K)=VISTEP
C   28 CONTINUE
C      DO 29 I=NIWL+1,NI
C      XLST=DEN(I,J,K)*XL3(I,J,K)*(ABS(TKE(I,J,K))**0.5)/VISCOS
C      ELM=2.5*XL3(I,J,K)*(1.-EXP(-0.016*XLST))
C      VISTEP=VISCOS+DEN(I,J,K)*CMU*SQRT(TKE(I,J,K))*ELM
C      VIS(I,J,K)=VISTEP
C   29 CONTINUE
C   26 CONTINUE
C
C   21 CONTINUE
C      ENDIF
C   
C      DO 11 K=1,NK
C      DO 12 J=1,NJR
C      DO 13 I=1,NIR
C      VIS(I,J,K)=VISCOS
C   13 CONTINUE
C      DO 14 I=NIRMS+1,NIRMF
C      VIS(I,J,K)=VISCOS
C   14 CONTINUE
C      DO 15 I=NI-2,NI
C      VIS(I,J,K)=VISCOS
C   15 CONTINUE
C   12 CONTINUE
C
C      DO 16 J=NJRN+1,NJ
C      DO 17 I=1,NIR
C      VIS(I,J,K)=VISCOS
C   17 CONTINUE
C      DO 18 I=NIRMS+1,NIRMF
C      VIS(I,J,K)=VISCOS
C   18 CONTINUE
C      DO 19 I=NI-2,NI
C      VIS(I,J,K)=VISCOS
C   19 CONTINUE
C   16 CONTINUE
C
C   11 CONTINUE
C
      DO 40 J=1,NJ
      DO 40 K=1,NK
      VIS(NI,J,K)=VIS(4,J,K)
      VIS(1,J,K)=VIS(NI-3,J,K)
 40   CONTINUE
C
      RETURN
      END
      SUBROUTINE STRAIN
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
C--------CALCULATION OF STRAIN RATES ----------------------------------
C
! NYAP  
      DIMENSION TL(NX,NY,NZ) 
      DO  K=2,NKM1
      DO  J=2,NJM1
      DO  I=2,NIM1	  
      TL(I,J,K)=TKE(I,J,K)**1.5/(TED(I,J,K)+TINY)
      ENDDO
      ENDDO
      ENDDO

!

      DO 1 K=2,NKM1
      DO 1 J=2,NJM1
      DO 1 I=2,NIM1
! NYAP      
      DLX=(1.-FX(I  ,J,K))*TL(I  ,J,K)+FX(I  ,J,K)*TL(I+1,J,K)
     1   -(1.-FX(I-1,J,K))*TL(I-1,J,K)-FX(I-1,J,K)*TL(I  ,J,K)
      DLY=(1.-FY(I  ,J,K))*TL(I  ,J,K)+FY(I,J  ,K)*TL(I,J+1,K)
     1   -(1.-FY(I,J-1,K))*TL(I,J-1,K)-FY(I,J-1,K)*TL(I  ,J,K)
      DLZ=(1.-FZ(I  ,J,K))*TL(I  ,J,K)+FZ(I,J,K  )*TL(I,J,K+1)
     1   -(1.-FZ(I,J,K-1))*TL(I,J,K-1)-FZ(I,J,K-1)*TL(I  ,J,K)
      DLDX=(DLX*XX(I,J,K)+DLY*YX(I,J,K)+DLZ*ZX(I,J,K))/VOL(I,J,K)
      DLDY=(DLX*XY(I,J,K)+DLY*YY(I,J,K)+DLZ*ZY(I,J,K))/VOL(I,J,K)
      DLDZ=(DLX*XZ(I,J,K)+DLY*YZ(I,J,K)+DLZ*ZZ(I,J,K))/VOL(I,J,K)
	DLDXJ2(I,J,K)=DLDX*DLDX+DLDY*DLDY+DLDZ*DLDZ
!
      DUX=(1.-FX(I  ,J,K))*U(I  ,J,K)+FX(I  ,J,K)*U(I+1,J,K)
     1   -(1.-FX(I-1,J,K))*U(I-1,J,K)-FX(I-1,J,K)*U(I  ,J,K)
      DVX=(1.-FX(I  ,J,K))*V(I  ,J,K)+FX(I  ,J,K)*V(I+1,J,K)
     1   -(1.-FX(I-1,J,K))*V(I-1,J,K)-FX(I-1,J,K)*V(I  ,J,K)
      DWX=(1.-FX(I  ,J,K))*W(I  ,J,K)+FX(I  ,J,K)*W(I+1,J,K)
     1   -(1.-FX(I-1,J,K))*W(I-1,J,K)-FX(I-1,J,K)*W(I  ,J,K)

      DUY=(1.-FY(I  ,J,K))*U(I  ,J,K)+FY(I,J  ,K)*U(I,J+1,K)
     1   -(1.-FY(I,J-1,K))*U(I,J-1,K)-FY(I,J-1,K)*U(I  ,J,K)

      DVY=(1.-FY(I  ,J,K))*V(I  ,J,K)+FY(I,J  ,K)*V(I,J+1,K)
     1   -(1.-FY(I,J-1,K))*V(I,J-1,K)-FY(I,J-1,K)*V(I  ,J,K)
      DWY=(1.-FY(I  ,J,K))*W(I  ,J,K)+FY(I,J  ,K)*W(I,J+1,K)
     1   -(1.-FY(I,J-1,K))*W(I,J-1,K)-FY(I,J-1,K)*W(I  ,J,K)
      DUZ=(1.-FZ(I  ,J,K))*U(I  ,J,K)+FZ(I,J,K  )*U(I,J,K+1)
     1   -(1.-FZ(I,J,K-1))*U(I,J,K-1)-FZ(I,J,K-1)*U(I  ,J,K)
      DVZ=(1.-FZ(I  ,J,K))*V(I  ,J,K)+FZ(I,J,K  )*V(I,J,K+1)
     1   -(1.-FZ(I,J,K-1))*V(I,J,K-1)-FZ(I,J,K-1)*V(I  ,J,K)
      DWZ=(1.-FZ(I  ,J,K))*W(I  ,J,K)+FZ(I,J,K  )*W(I,J,K+1)
     1   -(1.-FZ(I,J,K-1))*W(I,J,K-1)-FZ(I,J,K-1)*W(I  ,J,K)
C
      DUDXC(I,J,K)=(DUX*XX(I,J,K)
     1            +DUY*YX(I,J,K)
     1            +DUZ*ZX(I,J,K))/VOL(I,J,K)
C
      DUDYC(I,J,K)=(DUX*XY(I,J,K)
     1            +DUY*YY(I,J,K)
     1            +DUZ*ZY(I,J,K))/VOL(I,J,K)
C
      DUDZC(I,J,K)=(DUX*XZ(I,J,K)
     1            +DUY*YZ(I,J,K)
     1            +DUZ*ZZ(I,J,K))/VOL(I,J,K)
C
      DVDXC(I,J,K)=(DVX*XX(I,J,K)
     1            +DVY*YX(I,J,K)
     1            +DVZ*ZX(I,J,K))/VOL(I,J,K)
C
      DVDYC(I,J,K)=(DVX*XY(I,J,K)
     1            +DVY*YY(I,J,K)
     1            +DVZ*ZY(I,J,K))/VOL(I,J,K)
C
      DVDZC(I,J,K)=(DVX*XZ(I,J,K)
     1            +DVY*YZ(I,J,K)
     1            +DVZ*ZZ(I,J,K))/VOL(I,J,K)
C
      DWDXC(I,J,K)=(DWX*XX(I,J,K)
     1            +DWY*YX(I,J,K)
     1            +DWZ*ZX(I,J,K))/VOL(I,J,K)
C
      DWDYC(I,J,K)=(DWX*XY(I,J,K)
     1            +DWY*YY(I,J,K)
     1            +DWZ*ZY(I,J,K))/VOL(I,J,K)
C
      DWDZC(I,J,K)=(DWX*XZ(I,J,K)
     1            +DWY*YZ(I,J,K)
     1            +DWZ*ZZ(I,J,K))/VOL(I,J,K)
    1 CONTINUE
C
C----- CALCULATION OF EVM GENERATION TERM ---------------
      DO 2 K=2,NKM1
      DO 2 J=2,NJM1
      DO 2 I=2,NIM1
      VISTUR=VIS(I,J,K)-VISCOS
      IF (INLEVM.EQ.0) THEN
	PK1(I,J,K)=FLOAT(IEVM)*(VISTUR*(
     1 2.*(DUDXC(I,J,K))**2+2.*(DVDYC(I,J,K))**2+2.*(DWDZC(I,J,K))**2
     1+(DUDYC(I,J,K)+DVDXC(I,J,K))**2+(DVDZC(I,J,K)+DWDYC(I,J,K))**2
     1+(DWDXC(I,J,K)+DUDZC(I,J,K))**2 )/DEN(I,J,K))
     2          -FLOAT(1-IEVM)*(
     2     UU(I,J,K)*DUDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)
     2               +WW(I,J,K)*DWDZC(I,J,K)
     2    +UV(I,J,K)*(DUDYC(I,J,K)+DVDXC(I,J,K))
     2    +UW(I,J,K)*(DUDZC(I,J,K)+DWDXC(I,J,K))
     2    +VW(I,J,K)*(DVDZC(I,J,K)+DWDYC(I,J,K)))
	ELSE
C*-------- NLEVM ------------
	PK1(I,J,K)=-(UU(I,J,K)*DUDXC(I,J,K)+VV(I,J,K)*DVDYC(I,J,K)+
	1WW(I,J,K)*DWDZC(I,J,K)+UV(I,J,K)*(DUDYC(I,J,K)+DVDXC(I,J,K))
     1+UW(I,J,K)*(DUDZC(I,J,K)+DWDXC(I,J,K))+
     1VW(I,J,K)*(DVDZC(I,J,K)+DWDYC(I,J,K)))
	ENDIF
C*---------------------------
    2 CONTINUE
C
C-------- STRAIN RATE BOUNDARY CONDITIONS
C
      DO 3 K=2,NKM1
      DO 3 J=2,NJM1
C---------WEST BOUNDARY (ENTRY)
      DUDXC(1,J,K)=DUDXC(NI-3,J,K)
      DUDYC(1,J,K)=DUDYC(NI-3,J,K)
      DUDZC(1,J,K)=DUDZC(NI-3,J,K)      
      DVDXC(1,J,K)=DVDXC(NI-3,J,K)
      DVDYC(1,J,K)=DVDYC(NI-3,J,K)
      DVDZC(1,J,K)=DVDZC(NI-3,J,K)         
      DWDXC(1,J,K)=DWDXC(NI-3,J,K)
      DWDYC(1,J,K)=DWDYC(NI-3,J,K)
      DWDZC(1,J,K)=DWDZC(NI-3,J,K) 
C--------- EAST BOUNDARY    (EXIT)
      DUDXC(NI,J,K)=DUDXC(4,J,K)
      DUDYC(NI,J,K)=DUDYC(4,J,K)
      DUDZC(NI,J,K)=DUDZC(4,J,K)      
      DVDXC(NI,J,K)=DVDXC(4,J,K)
      DVDYC(NI,J,K)=DVDYC(4,J,K)
      DVDZC(NI,J,K)=DVDZC(4,J,K)         
      DWDXC(NI,J,K)=DWDXC(4,J,K)
      DWDYC(NI,J,K)=DWDYC(4,J,K)
      DWDZC(NI,J,K)=DWDZC(4,J,K) 
    3 CONTINUE
C
      DO 4 J=2,NJM1
      DO 4 I=2,NIM1
C------- BOTTOM BOUNDARY  (WALL OR SYMMETRY)
      DUDXC(I,J,1)=FLOAT(1-LWB)*DUDXC(I,J,2)
      DUDYC(I,J,1)=FLOAT(1-LWB)*DUDYC(I,J,2)
      DUDZC(I,J,1)=FLOAT(LWB)*DUDZC(I,J,2)
      DVDXC(I,J,1)=FLOAT(1-LWB)*DVDXC(I,J,2)
      DVDYC(I,J,1)=FLOAT(1-LWB)*DVDYC(I,J,2)
      DVDZC(I,J,1)=FLOAT(LWB)*DVDZC(I,J,2)
      DWDXC(I,J,1)=0.
      DWDYC(I,J,1)=0.
      DWDZC(I,J,1)=0.
C------- TOP BOUNDARY   (WALL OR SYMMETRY)
      DUDXC(I,J,NK)=FLOAT(1-LWT)*DUDXC(I,J,NKM1)
      DUDYC(I,J,NK)=FLOAT(1-LWT)*DUDYC(I,J,NKM1)
      DUDZC(I,J,NK)=FLOAT(LWT)*DUDZC(I,J,NKM1)
      DVDXC(I,J,NK)=FLOAT(1-LWT)*DVDXC(I,J,NKM1)
      DVDYC(I,J,NK)=FLOAT(1-LWT)*DVDYC(I,J,NKM1)
      DVDZC(I,J,NK)=FLOAT(LWT)*DVDZC(I,J,NKM1)
      DWDXC(I,J,NK)=0.
      DWDYC(I,J,NK)=0.
      DWDZC(I,J,NK)=0.
    4 CONTINUE
C
      DO 5 K=2,NKM1
      DO 5 I=2,NIM1
C------- SOUTH BOUNDARY (WALL OR SYMMETRY)
      DUDXC(I,1,K)=FLOAT(LWS)*DUDXC(I,2,K)
      DUDYC(I,1,K)=FLOAT(LWS)*DUDYC(I,2,K)
      DUDZC(I,1,K)=FLOAT(LWS)*DUDZC(I,2,K)
      DVDXC(I,1,K)=FLOAT(LWS)*DVDXC(I,2,K)
      DVDYC(I,1,K)=FLOAT(LWS)*DVDYC(I,2,K)
      DVDZC(I,1,K)=FLOAT(LWS)*DVDZC(I,2,K)
      DWDXC(I,1,K)=FLOAT(1-LWS)*DWDXC(I,2,K)
      DWDYC(I,1,K)=FLOAT(LWS)*DWDYC(I,2,K)
      DWDZC(I,1,K)=FLOAT(1-LWS)*DWDZC(I,2,K)
C-------  NORTH BOUNDARY (WALL OR SYMMETRY)
      DUDXC(I,NJ,K)=FLOAT(LWN)*DUDXC(I,NJM1,K)
      DUDYC(I,NJ,K)=FLOAT(LWN)*DUDYC(I,NJM1,K)
      DUDZC(I,NJ,K)=FLOAT(LWN)*DUDZC(I,NJM1,K)
      DVDXC(I,NJ,K)=FLOAT(LWN)*DVDXC(I,NJM1,K)
      DVDYC(I,NJ,K)=FLOAT(LWN)*DVDYC(I,NJM1,K)
      DVDZC(I,NJ,K)=FLOAT(LWN)*DVDZC(I,NJM1,K)
      DWDXC(I,NJ,K)=FLOAT(1-LWN)*DWDXC(I,NJM1,K)
      DWDYC(I,NJ,K)=FLOAT(LWN)*DWDYC(I,NJM1,K)
      DWDZC(I,NJ,K)=FLOAT(1-LWN)*DWDZC(I,NJM1,K)
    5 CONTINUE   
C
C-------- NEAR-WALL DISSIPATION RATE AND ROOT-K DERIVATIVES --------
C
      DO 7 K=2,NKM1
      DO 7 J=2,NJM1
      DO 7 I=2,NIM1
C ------ K DERIVATIVES ----------------------------------------------
      DKX=(1.-FX(I  ,J,K))*SQRT(TKE(I  ,J,K))
     1    +FX(I  ,J,K)*SQRT(TKE(I+1,J,K))
     1    -(1.-FX(I-1,J,K))*SQRT(TKE(I-1,J,K))
     1    -FX(I-1,J,K)*SQRT(TKE(I  ,J,K))
      DKY=(1.-FY(I  ,J,K))*SQRT(TKE(I  ,J,K))
     1    +FY(I,J  ,K)*SQRT(TKE(I,J+1,K))
     1    -(1.-FY(I,J-1,K))*SQRT(TKE(I,J-1,K))
     1    -FY(I,J-1,K)*SQRT(TKE(I  ,J,K))
      DKZ=(1.-FZ(I  ,J,K))*SQRT(TKE(I  ,J,K))
     1    +FZ(I,J,K  )*SQRT(TKE(I,J,K+1))
     1    -(1.-FZ(I,J,K-1))*SQRT(TKE(I,J,K-1))
     1    -FZ(I,J,K-1)*SQRT(TKE(I  ,J,K))
C
      DKDX(I,J,K)=(DKX*XX(I,J,K)
     1            +DKY*YX(I,J,K)
     1            +DKZ*ZX(I,J,K))/VOL(I,J,K)
C
      DKDY(I,J,K)=(DKX*XY(I,J,K)
     1            +DKY*YY(I,J,K)
     1            +DKZ*ZY(I,J,K))/VOL(I,J,K)
C
      DKDZ(I,J,K)=(DKX*XZ(I,J,K)
     1            +DKY*YZ(I,J,K)
     1            +DKZ*ZZ(I,J,K))/VOL(I,J,K)
C----- NEAR-WALL EPSILON --------------------------------
      EDNW(I,J,K)=2.*VISCOS*(DKDX(I,J,K)*DKDX(I,J,K)+
     1                       DKDY(I,J,K)*DKDY(I,J,K)+
     1                       DKDZ(I,J,K)*DKDZ(I,J,K))/DEN(I,J,K)
    7 CONTINUE
C*-------- NLEVM ------------------		
	IF (INLEVM.EQ.1) CALL NLEVM1
C*---------------------------------

c      IF(INWM.EQ.1)RETURN
C
C------ DOUBLE DERIVATIVES OF VELOCITY AND E-TERM ----------------
C
      DO 6 K=2,NKM1
      DO 6 J=2,NJM1
      DO 6 I=2,NIM1
C
C---------- U DERIVATIVES
C
      DUX1=(1.-FX(I  ,J,K))*DUDXC(I  ,J,K)+FX(I  ,J,K)*DUDXC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DUDXC(I-1,J,K)-FX(I-1,J,K)*DUDXC(I  ,J,K)
      DUX2=(1.-FY(I  ,J,K))*DUDXC(I  ,J,K)+FY(I,J  ,K)*DUDXC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DUDXC(I,J-1,K)-FY(I,J-1,K)*DUDXC(I  ,J,K)
      DUX3=(1.-FZ(I  ,J,K))*DUDXC(I  ,J,K)+FZ(I,J,K  )*DUDXC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DUDXC(I,J,K-1)-FZ(I,J,K-1)*DUDXC(I  ,J,K)
C
      DUXX=(DUX1*XX(I,J,K)
     1    +DUX2*YX(I,J,K)
     1    +DUX3*ZX(I,J,K))/VOL(I,J,K)
C
      DUXY=(DUX1*XY(I,J,K)
     1    +DUX2*YY(I,J,K)
     1    +DUX3*ZY(I,J,K))/VOL(I,J,K)
C
      DUXZ=(DUX1*XZ(I,J,K)
     1    +DUX2*YZ(I,J,K)
     1    +DUX3*ZZ(I,J,K))/VOL(I,J,K)
C
      DUY1=(1.-FX(I  ,J,K))*DUDYC(I  ,J,K)+FX(I  ,J,K)*DUDYC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DUDYC(I-1,J,K)-FX(I-1,J,K)*DUDYC(I  ,J,K)
      DUY2=(1.-FY(I  ,J,K))*DUDYC(I  ,J,K)+FY(I,J  ,K)*DUDYC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DUDYC(I,J-1,K)-FY(I,J-1,K)*DUDYC(I  ,J,K)
      DUY3=(1.-FZ(I  ,J,K))*DUDYC(I  ,J,K)+FZ(I,J,K  )*DUDYC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DUDYC(I,J,K-1)-FZ(I,J,K-1)*DUDYC(I  ,J,K)
C
      DUYY=(DUY1*XY(I,J,K)
     1    +DUY2*YY(I,J,K)
     1    +DUY3*ZY(I,J,K))/VOL(I,J,K)
C
      DUYZ=(DUY1*XZ(I,J,K)
     1    +DUY2*YZ(I,J,K)
     1    +DUY3*ZZ(I,J,K))/VOL(I,J,K)
C
      DUZ1=(1.-FX(I  ,J,K))*DUDZC(I  ,J,K)+FX(I  ,J,K)*DUDZC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DUDZC(I-1,J,K)-FX(I-1,J,K)*DUDZC(I  ,J,K)
      DUZ2=(1.-FY(I  ,J,K))*DUDZC(I  ,J,K)+FY(I,J  ,K)*DUDZC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DUDZC(I,J-1,K)-FY(I,J-1,K)*DUDZC(I  ,J,K)
      DUZ3=(1.-FZ(I  ,J,K))*DUDZC(I  ,J,K)+FZ(I,J,K  )*DUDZC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DUDZC(I,J,K-1)-FZ(I,J,K-1)*DUDZC(I  ,J,K)
C
      DUZZ=(DUZ1*XZ(I,J,K)
     1    +DUZ2*YZ(I,J,K)
     1    +DUZ3*ZZ(I,J,K))/VOL(I,J,K)
C
C---------- V DERIVATIVES
C
      DVX1=(1.-FX(I  ,J,K))*DVDXC(I  ,J,K)+FX(I  ,J,K)*DVDXC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DVDXC(I-1,J,K)-FX(I-1,J,K)*DVDXC(I  ,J,K)
      DVX2=(1.-FY(I  ,J,K))*DVDXC(I  ,J,K)+FY(I,J  ,K)*DVDXC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DVDXC(I,J-1,K)-FY(I,J-1,K)*DVDXC(I  ,J,K)
      DVX3=(1.-FZ(I  ,J,K))*DVDXC(I  ,J,K)+FZ(I,J,K  )*DVDXC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DVDXC(I,J,K-1)-FZ(I,J,K-1)*DVDXC(I  ,J,K)
C
      DVXX=(DVX1*XX(I,J,K)
     1     +DVX2*YX(I,J,K)
     1     +DVX3*ZX(I,J,K))/VOL(I,J,K)
C
      DVXY=(DVX1*XY(I,J,K)
     1     +DVX2*YY(I,J,K)
     1     +DVX3*ZY(I,J,K))/VOL(I,J,K)
C
      DVXZ=(DVX1*XZ(I,J,K)
     1     +DVX2*YZ(I,J,K)
     1     +DVX3*ZZ(I,J,K))/VOL(I,J,K)
C
      DVY1=(1.-FX(I  ,J,K))*DVDYC(I  ,J,K)+FX(I  ,J,K)*DVDYC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DVDYC(I-1,J,K)-FX(I-1,J,K)*DVDYC(I  ,J,K)
      DVY2=(1.-FY(I  ,J,K))*DVDYC(I  ,J,K)+FY(I,J  ,K)*DVDYC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DVDYC(I,J-1,K)-FY(I,J-1,K)*DVDYC(I  ,J,K)
      DVY3=(1.-FZ(I  ,J,K))*DVDYC(I  ,J,K)+FZ(I,J,K  )*DVDYC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DVDYC(I,J,K-1)-FZ(I,J,K-1)*DVDYC(I  ,J,K)
C
      DVYY=(DVY1*XY(I,J,K)
     1     +DVY2*YY(I,J,K)
     1     +DVY3*ZY(I,J,K))/VOL(I,J,K)
C
      DVYZ=(DVY1*XZ(I,J,K)
     1     +DVY2*YZ(I,J,K)
     1     +DVY3*ZZ(I,J,K))/VOL(I,J,K)
C
      DVZ1=(1.-FX(I  ,J,K))*DVDZC(I  ,J,K)+FX(I  ,J,K)*DVDZC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DVDZC(I-1,J,K)-FX(I-1,J,K)*DVDZC(I  ,J,K)
      DVZ2=(1.-FY(I  ,J,K))*DVDZC(I  ,J,K)+FY(I,J  ,K)*DVDZC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DVDZC(I,J-1,K)-FY(I,J-1,K)*DVDZC(I  ,J,K)
      DVZ3=(1.-FZ(I  ,J,K))*DVDZC(I  ,J,K)+FZ(I,J,K  )*DVDZC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DVDZC(I,J,K-1)-FZ(I,J,K-1)*DVDZC(I  ,J,K)
C
      DVZZ=(DVZ1*XZ(I,J,K)
     1     +DVZ2*YZ(I,J,K)
     1     +DVZ3*ZZ(I,J,K))/VOL(I,J,K)
C
C---------- W DERIVATIVES
C
      DWX1=(1.-FX(I  ,J,K))*DWDXC(I  ,J,K)+FX(I  ,J,K)*DWDXC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DWDXC(I-1,J,K)-FX(I-1,J,K)*DWDXC(I  ,J,K)
      DWX2=(1.-FY(I  ,J,K))*DWDXC(I  ,J,K)+FY(I,J  ,K)*DWDXC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DWDXC(I,J-1,K)-FY(I,J-1,K)*DWDXC(I  ,J,K)
      DWX3=(1.-FZ(I  ,J,K))*DWDXC(I  ,J,K)+FZ(I,J,K  )*DWDXC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DWDXC(I,J,K-1)-FZ(I,J,K-1)*DWDXC(I  ,J,K)
C
      DWXX=(DWX1*XX(I,J,K)
     1     +DWX2*YX(I,J,K)
     1     +DWX3*ZX(I,J,K))/VOL(I,J,K)
C
      DWXY=(DWX1*XY(I,J,K)
     1     +DWX2*YY(I,J,K)
     1     +DWX3*ZY(I,J,K))/VOL(I,J,K)
C
      DWXZ=(DWX1*XZ(I,J,K)
     1     +DWX2*YZ(I,J,K)
     1     +DWX3*ZZ(I,J,K))/VOL(I,J,K)
C
      DWY1=(1.-FX(I  ,J,K))*DWDYC(I  ,J,K)+FX(I  ,J,K)*DWDYC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DWDYC(I-1,J,K)-FX(I-1,J,K)*DWDYC(I  ,J,K)
      DWY2=(1.-FY(I  ,J,K))*DWDYC(I  ,J,K)+FY(I,J  ,K)*DWDYC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DWDYC(I,J-1,K)-FY(I,J-1,K)*DWDYC(I  ,J,K)
      DWY3=(1.-FZ(I  ,J,K))*DWDYC(I  ,J,K)+FZ(I,J,K  )*DWDYC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DWDYC(I,J,K-1)-FZ(I,J,K-1)*DWDYC(I  ,J,K)
C
      DWYY=(DWY1*XY(I,J,K)
     1     +DWY2*YY(I,J,K)
     1     +DWY3*ZY(I,J,K))/VOL(I,J,K)
C
      DWYZ=(DWY1*XZ(I,J,K)
     1     +DWY2*YZ(I,J,K)
     1     +DWY3*ZZ(I,J,K))/VOL(I,J,K)
C
      DWZ1=(1.-FX(I  ,J,K))*DWDZC(I  ,J,K)+FX(I  ,J,K)*DWDZC(I+1,J,K)
     1    -(1.-FX(I-1,J,K))*DWDZC(I-1,J,K)-FX(I-1,J,K)*DWDZC(I  ,J,K)
      DWZ2=(1.-FY(I  ,J,K))*DWDZC(I  ,J,K)+FY(I,J  ,K)*DWDZC(I,J+1,K)
     1    -(1.-FY(I,J-1,K))*DWDZC(I,J-1,K)-FY(I,J-1,K)*DWDZC(I  ,J,K)
      DWZ3=(1.-FZ(I  ,J,K))*DWDZC(I  ,J,K)+FZ(I,J,K  )*DWDZC(I,J,K+1)
     1    -(1.-FZ(I,J,K-1))*DWDZC(I,J,K-1)-FZ(I,J,K-1)*DWDZC(I  ,J,K)
C
      DWZZ=(DWZ1*XZ(I,J,K)
     1     +DWZ2*YZ(I,J,K)
     1     +DWZ3*ZZ(I,J,K))/VOL(I,J,K)
C
C------------   E  TERM    ---------------------
C
      VISTUR=VIS(I,J,K)-VISCOS
	ETRM(I,J,K)=2.*VISCOS*VISTUR*(
     1            DUXX*DUXX + DUYY*DUYY + DUZZ*DUZZ +
     1            DVXX*DVXX + DVYY*DVYY + DVZZ*DVZZ +
     1            DWXX*DWXX + DWYY*DWYY + DWZZ*DWZZ +
     1        2.*(DUXY*DUXY + DUXZ*DUXZ + DUYZ*DUYZ +
     1            DVXY*DVXY + DVXZ*DVXZ + DVYZ*DVYZ +
     1            DWXY*DWXY + DWXZ*DWXZ + DWYZ*DWYZ ) )
     1            /( DEN(I,J,K)*DEN(I,J,K))
C*-------- NLEVM ------------
	IF (INLEVM.EQ.1) THEN
	RTLDT =DEN(I,J,K)*(TKE(I,J,K)**2)/(TED(I,J,K)+TINY)/VISCOS
	IF (RTLDT.LE.250.0) THEN
	ETRM(I,J,K)=0.0022*DEN(I,J,K)*STLD(I,J,K)*(TKE(I,J,K)**2)*
	1ETRM(I,J,K)/(2.0*VISCOS*TED(I,J,K)+tiny)
	ENDIF
	IF (RTLDT.GT.250.0) then 
	ETRM(I,J,K)=0.0
	endif
	ENDIF
C*---------------------------
    6 CONTINUE
C
C ------------ EVM TURBULENT STRESSES ------------
C
      IF(IEVM.EQ.1.OR.(NREAD.EQ.0.AND.IT.EQ.1))THEN
      RIT=FLOAT(IT-1)
      CRIT=MIN(RIT,1.)
      CRIT=MAX(CRIT,FLOAT(NREAD))
      CRIT=1.
	IF (INLEVM.EQ.0) THEN
      DO 11 K=2,NKM1
      DO 11 J=2,NJM1
      DO 11 I=2,NIM1
      VISTUR=VIS(I,J,K)-VISCOS
      UU(I,J,K)=2.0*TKE(I,J,K)/3.0-VISTUR*2.*DUDXC(I,J,K)
      VV(I,J,K)=2.0*TKE(I,J,K)/3.0-VISTUR*2.*DVDYC(I,J,K)
      WW(I,J,K)=2.0*TKE(I,J,K)/3.0-VISTUR*2.*DWDXC(I,J,K)
      UV(I,J,K)=-CRIT*VISTUR*(DUDYC(I,J,K)+DVDXC(I,J,K))
      UW(I,J,K)=-CRIT*VISTUR*(DWDXC(I,J,K)+DUDZC(I,J,K))
      VW(I,J,K)=-CRIT*VISTUR*(DWDYC(I,J,K)+DVDZC(I,J,K))
   11 CONTINUE
	ENDIF
      DO 12 J=2,NJM1
      DO 12 I=2,NIM1
      UU(I,J,1)=UU(I,J,2)
      VV(I,J,1)=VV(I,J,2)
      WW(I,J,1)=WW(I,J,2)
      UV(I,J,1)=UV(I,J,2)
C*      UU(I,J,1)=0.
C*      VV(I,J,1)=0.
C*      WW(I,J,1)=0.
C*      UV(I,J,1)=0.
      UW(I,J,1)=0.
      VW(I,J,1)=0.
      UU(I,J,NK)=0.
      VV(I,J,NK)=0.
      WW(I,J,NK)=0.
      UV(I,J,NK)=0.
      UW(I,J,NK)=0.
      VW(I,J,NK)=0.
   12 CONTINUE
      CALL MODTUU
      CALL MODTVV
      CALL MODTWW
      CALL MODTUV
      CALL MODTUW
      CALL MODTVW
      GO TO 15
      DO 13 K=2,NKM1
      DO 13 J=2,NJM1
      UU(NI,J,K)=UU(NIM1,J,K)
      VV(NI,J,K)=VV(NIM1,J,K)
      WW(NI,J,K)=WW(NIM1,J,K)
      UV(NI,J,K)=UV(NIM1,J,K)
      UW(NI,J,K)=UW(NIM1,J,K)
      VW(NI,J,K)=VW(NIM1,J,K)
      UU(1,J,K)=UU(NIM1,J,K)
      VV(1,J,K)=VV(NIM1,J,K)
      WW(1,J,K)=WW(NIM1,J,K)
      UV(1,J,K)=UV(NIM1,J,K)
      UW(1,J,K)=UW(NIM1,J,K)
      VW(1,J,K)=VW(NIM1,J,K)
   13 CONTINUE
      DO 14 K=2,NKM1
      DO 14 I=2,NIM1
      UU(I,1,K)=0.
      VV(I,1,K)=0.
      WW(I,1,K)=0.
      UV(I,1,K)=0.
      VW(I,1,K)=0.
      UW(I,1,K)=0.
      UU(I,NJ,K)=0.
      VV(I,NJ,K)=0.
      WW(I,NJ,K)=0.
      UV(I,NJ,K)=0.
      VW(I,NJ,K)=0.
      UW(I,NJ,K)=0.
   14 CONTINUE
   15 CONTINUE
      ENDIF
      RETURN
      END
C
      SUBROUTINE ADVECT(PRL,PRT)
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
C--------CALCULATION OF ADVECTION COEFFICIENTS ------------------------
C
      DO 1 K=2,NKM1
      DO 1 J=2,NJM1
      DO 1 I=2,NIM1
C
C---------- CELL FACE VISCOCITIES ------------------------------------
C
C      FX(1,J,K)=0.5
C      FX(NIM1,J,K)=0.5
C      FY(I,1,K)=0.5
C      FY(I,NJM1,K)=0.5
C      FZ(I,J,1)=0.5
      FZ(I,J,NKM1)=0.5
C
       VISE=VIS(I,J,K)*(1.-FX(I,J,K))+VIS(I+1,J,K)*FX(I,J,K)
       VISW=VIS(I-1,J,K)*(1.-FX(I-1,J,K))+VIS(I,J,K)*FX(I-1,J,K)
       VISN=VIS(I,J,K)*(1.-FY(I,J,K))+VIS(I,J+1,K)*FY(I,J,K)
       VISS=VIS(I,J-1,K)*(1.-FY(I,J-1,K))+VIS(I,J,K)*FY(I,J-1,K)
       VIST=VIS(I,J,K)*(1.-FZ(I,J,K))+VIS(I,J,K+1)*FZ(I,J,K)
       VISB=VIS(I,J,K-1)*(1.-FZ(I,J,K-1))+VIS(I,J,K)*FZ(I,J,K-1)
C
C      FX(1,J,K)=0.
C      FX(NIM1,J,K)=1.
C      FY(I,1,K)=0.
C      FY(I,NJM1,K)=1.
C      FZ(I,J,1)=0.
      FZ(I,J,NKM1)=1.
C
      VISE=(VISCOS/PRL)+(VISE-VISCOS)/PRT
      VISW=(VISCOS/PRL)+(VISW-VISCOS)/PRT
      VISN=(VISCOS/PRL)+(VISN-VISCOS)/PRT
      VISS=(VISCOS/PRL)+(VISS-VISCOS)/PRT
      VIST=(VISCOS/PRL)+(VIST-VISCOS)/PRT
      VISB=(VISCOS/PRL)+(VISB-VISCOS)/PRT
C
C ------ CONTROL VOLUMES FOR DUMMY (STAGGERED) VELOCITY NODES ----
C
       VOLE=VOL(I,J,K)*(1.-FX(I,J,K))+VOL(I+1,J,K)*FX(I,J,K)
       VOLW=VOL(I-1,J,K)*(1.-FX(I-1,J,K))+VOL(I,J,K)*FX(I-1,J,K)
       VOLN=VOL(I,J,K)*(1.-FY(I,J,K))+VOL(I,J+1,K)*FY(I,J,K)
       VOLS=VOL(I,J-1,K)*(1.-FY(I,J-1,K))+VOL(I,J,K)*FY(I,J-1,K)
       VOLT=VOL(I,J,K)*(1.-FZ(I,J,K))+VOL(I,J,K+1)*FZ(I,J,K)
       VOLB=VOL(I,J,K-1)*(1.-FZ(I,J,K-1))+VOL(I,J,K)*FZ(I,J,K-1)
C
C------ EAST-WEST ADVECTION --------------------------------------
C
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
C
      Q11E=(1.-FX(I,J,K))*(XX(I,J,K)**2+XY(I,J,K)**2+XZ(I,J,K)**2)
     1+FX(I,J,K)*(XX(I+1,J,K)**2+XY(I+1,J,K)**2+XZ(I+1,J,K)**2)
      Q11W=(1.-FX(I-1,J,K))*(XX(I-1,J,K)**2+XY(I-1,J,K)**2+
     1                       XZ(I-1,J,K)**2)
     1+FX(I-1,J,K)*(XX(I,J,K)**2+XY(I,J,K)**2+XZ(I,J,K)**2)
C
      DE=VISE*Q11E/VOLE
      DW=VISW*Q11W/VOLW
C
      FLOWE=CX(I  ,J,K)
      FLOWW=CX(I-1,J,K)
C
      CE=DENE*FLOWE
      CW=DENW*FLOWW
      PECE=CE/DE
      PECW=CW/DW
C
C        UPWIND COEFFICIENTS
C
      AE(I,J,K)=DE+AMAX1(-CE,0.)
      AW(I,J,K)=DW+AMAX1( CW,0.)
C
C -------- NORTH-SOUTH ADVECTION  -------------------------------------
C
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
C
      Q22N=(1.-FY(I,J,K))*(YX(I,J,K)**2+YY(I,J,K)**2+YZ(I,J,K)**2)
     1+FY(I,J,K)*(YX(I,J+1,K)**2+YY(I,J+1,K)**2+YZ(I,J+1,K)**2)
      Q22S=(1.-FY(I,J-1,K))*(YX(I,J-1,K)**2+YY(I,J-1,K)**2+
     1                       YZ(I,J-1,K)**2)
     1+FY(I,J-1,K)*(YX(I,J,K)**2+YY(I,J,K)**2+YZ(I,J,K)**2)
C
      DN=VISN*Q22N/VOLN
      DS=VISS*Q22S/VOLS
C
      FLOWN=CY(I,J  ,K)
      FLOWS=CY(I,J-1,K)
C
      CN=DENN*FLOWN
      CS=DENS*FLOWS
      PECN=CN/DN
      PECS=CS/DS
C
C      UPWIND COEFFICIENTS  
C
      AN(I,J,K)=DN+AMAX1(-CN,0.)
      AS(I,J,K)=DS+AMAX1( CS,0.)
C
C------ TOP-BOTTOM ADVECTION ----------------------------------------
C
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C
      Q33T=(1.-FZ(I,J,K))*(ZX(I,J,K)**2+ZY(I,J,K)**2+ZZ(I,J,K)**2)
     1+FZ(I,J,K)*(ZX(I,J,K+1)**2+ZY(I,J,K+1)**2+ZZ(I,J,K+1)**2)
      Q33B=(1.-FZ(I,J,K-1))*(ZX(I,J,K-1)**2+ZY(I,J,K-1)**2+
     1                       ZZ(I,J,K-1)**2)
     1+FZ(I,J,K-1)*(ZX(I,J,K)**2+ZY(I,J,K)**2+ZZ(I,J,K)**2)
C
      DT=VIST*Q33T/VOLT
      DB=VISB*Q33B/VOLB
C
      FLOWT=CZ(I,J,K  )
      FLOWB=CZ(I,J,K-1)
C
      CT=DENT*FLOWT
      CB=DENB*FLOWB
      PECT=CT/DT
      PECB=CB/DB
C
C    UPWIND COEFFICIENTS
C
      AT(I,J,K)=DT+AMAX1(-CT,0.)
      AB(I,J,K)=DB+AMAX1( CB,0.)
C
    1 CONTINUE
C      
      RETURN
      END
      SUBROUTINE CRDIFF(PRL,PRT,ITRB,PHI)
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
C
      DIMENSION PHI(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
C--------CALCULATION OF CROSS-DIFFUSION TERMS ------------------------
C
C      IF(IEVM.EQ.1)ITRB=1
C
      DO 1 K=2,NKM1
      DO 1 J=2,NJM1
      DO 1 I=2,NIM1
C
C---------- CELL FACE VISCOCITIES ------------------------------------
C
       VISE=VIS(I,J,K)*(1.-FX(I,J,K))+VIS(I+1,J,K)*FX(I,J,K)
       VISW=VIS(I-1,J,K)*(1.-FX(I-1,J,K))+VIS(I,J,K)*FX(I-1,J,K)
       VISN=VIS(I,J,K)*(1.-FY(I,J,K))+VIS(I,J+1,K)*FY(I,J,K)
       VISS=VIS(I,J-1,K)*(1.-FY(I,J-1,K))+VIS(I,J,K)*FY(I,J-1,K)
       VIST=VIS(I,J,K)*(1.-FZ(I,J,K))+VIS(I,J,K+1)*FZ(I,J,K)
       VISB=VIS(I,J,K-1)*(1.-FZ(I,J,K-1))+VIS(I,J,K)*FZ(I,J,K-1)
C
      VISE=(VISCOS/PRL)+(VISE-VISCOS)/PRT
      VISW=(VISCOS/PRL)+(VISW-VISCOS)/PRT
      VISN=(VISCOS/PRL)+(VISN-VISCOS)/PRT
      VISS=(VISCOS/PRL)+(VISS-VISCOS)/PRT
      VIST=(VISCOS/PRL)+(VIST-VISCOS)/PRT
      VISB=(VISCOS/PRL)+(VISB-VISCOS)/PRT
C
C ------ CONTROL VOLUMES FOR DUMMY (STAGGERED) VELOCITY NODES ----
C
       VOLE=VOL(I,J,K)*(1.-FX(I,J,K))+VOL(I+1,J,K)*FX(I,J,K)
       VOLW=VOL(I-1,J,K)*(1.-FX(I-1,J,K))+VOL(I,J,K)*FX(I-1,J,K)
       VOLN=VOL(I,J,K)*(1.-FY(I,J,K))+VOL(I,J+1,K)*FY(I,J,K)
       VOLS=VOL(I,J-1,K)*(1.-FY(I,J-1,K))+VOL(I,J,K)*FY(I,J-1,K)
       VOLT=VOL(I,J,K)*(1.-FZ(I,J,K))+VOL(I,J,K+1)*FZ(I,J,K)
       VOLB=VOL(I,J,K-1)*(1.-FZ(I,J,K-1))+VOL(I,J,K)*FZ(I,J,K-1)
C
C      CALCULATE THE CORNER VALUES OF PHI
C
      PHITNE=
     1 PHI(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+PHI(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+PHI(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+PHI(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+PHI(I,J,K+1)*(1.-FX(I,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+PHI(I,J+1,K+1)*(1.-FX(I,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+PHI(I+1,J+1,K+1)*FX(I,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+PHI(I+1,J,K+1)*FX(I,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      PHITNW=
     1 PHI(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+PHI(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*(1.-FZ(I,J,K))
     1+PHI(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*(1.-FZ(I,J,K))
     1+PHI(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*(1.-FZ(I,J,K))
     1+PHI(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*(1.-FY(I,J,K+1))*FZ(I,J,K)
     1+PHI(I-1,J+1,K+1)*(1.-FX(I-1,J+1,K+1))*FY(I,J,K+1)*FZ(I,J,K)
     1+PHI(I,J+1,K+1)*FX(I-1,J+1,K+1)*FY(I,J,K+1)*FZ(I,J,K)
     1+PHI(I,J,K+1)*FX(I-1,J,K+1)*(1.-FY(I,J,K+1))*FZ(I,J,K)
C
      PHIBNE=
     1 PHI(I,J,K-1)*(1.-FX(I,J,K-1))*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+PHI(I,J+1,K-1)*(1.-FX(I,J+1,K-1))*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+PHI(I+1,J+1,K-1)*FX(I,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+PHI(I+1,J,K-1)*FX(I,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+PHI(I,J,K)*(1.-FX(I,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+PHI(I,J+1,K)*(1.-FX(I,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+PHI(I+1,J+1,K)*FX(I,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+PHI(I+1,J,K)*FX(I,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      PHIBNW=
     1 PHI(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*(1.-FY(I,J,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+PHI(I-1,J+1,K-1)*(1.-FX(I-1,J+1,K-1))*FY(I,J,K-1)*
     1                (1.-FZ(I,J,K-1))
     1+PHI(I,J+1,K-1)*FX(I-1,J+1,K-1)*FY(I,J,K-1)*(1.-FZ(I,J,K-1))
     1+PHI(I,J,K-1)*FX(I-1,J,K-1)*(1.-FY(I,J,K-1))*(1.-FZ(I,J,K-1))
     1+PHI(I-1,J,K)*(1.-FX(I-1,J,K))*(1.-FY(I,J,K))*FZ(I,J,K-1)
     1+PHI(I-1,J+1,K)*(1.-FX(I-1,J+1,K))*FY(I,J,K)*FZ(I,J,K-1)
     1+PHI(I,J+1,K)*FX(I-1,J+1,K)*FY(I,J,K)*FZ(I,J,K-1)
     1+PHI(I,J,K)*FX(I-1,J,K)*(1.-FY(I,J,K))*FZ(I,J,K-1)
C
      PHITSE=
     1 PHI(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+PHI(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+PHI(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+PHI(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+PHI(I,J-1,K+1)*(1.-FX(I,J-1,K+1))*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
     1+PHI(I,J,K+1)*(1.-FX(I,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+PHI(I+1,J,K+1)*FX(I,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+PHI(I+1,J-1,K+1)*FX(I,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      PHITSW=
     1 PHI(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+PHI(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+PHI(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*(1.-FZ(I,J,K))
     1+PHI(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*(1.-FZ(I,J,K))
     1+PHI(I-1,J-1,K+1)*(1.-FX(I-1,J-1,K+1))*(1.-FY(I,J-1,K+1))*
     1                  FZ(I,J,K)
     1+PHI(I-1,J,K+1)*(1.-FX(I-1,J,K+1))*FY(I,J-1,K+1)*FZ(I,J,K)
     1+PHI(I,J,K+1)*FX(I-1,J,K+1)*FY(I,J-1,K+1)*FZ(I,J,K)
     1+PHI(I,J-1,K+1)*FX(I-1,J-1,K+1)*(1.-FY(I,J-1,K+1))*FZ(I,J,K)
C
      PHIBSE=
     1 PHI(I,J-1,K-1)*(1.-FX(I,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+PHI(I,J,K-1)*(1.-FX(I,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+PHI(I+1,J,K-1)*FX(I,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+PHI(I+1,J-1,K-1)*FX(I,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                 (1.-FZ(I,J,K-1))
     1+PHI(I,J-1,K)*(1.-FX(I,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+PHI(I,J,K)*(1.-FX(I,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+PHI(I+1,J,K)*FX(I,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+PHI(I+1,J-1,K)*FX(I,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
      PHIBSW=
     1 PHI(I-1,J-1,K-1)*(1.-FX(I-1,J-1,K-1))*(1.-FY(I,J-1,K-1))*
     1                  (1.-FZ(I,J,K-1))
     1+PHI(I-1,J,K-1)*(1.-FX(I-1,J,K-1))*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+PHI(I,J,K-1)*FX(I-1,J,K-1)*FY(I,J-1,K-1)*(1.-FZ(I,J,K-1))
     1+PHI(I,J-1,K-1)*FX(I-1,J-1,K-1)*(1.-FY(I,J-1,K-1))*
     1                (1.-FZ(I,J,K-1))
     1+PHI(I-1,J-1,K)*(1.-FX(I-1,J-1,K))*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
     1+PHI(I-1,J,K)*(1.-FX(I-1,J,K))*FY(I,J-1,K)*FZ(I,J,K-1)
     1+PHI(I,J,K)*FX(I-1,J,K)*FY(I,J-1,K)*FZ(I,J,K-1)
     1+PHI(I,J-1,K)*FX(I-1,J-1,K)*(1.-FY(I,J-1,K))*FZ(I,J,K-1)
C
C------- EAST-WEST PROJECTED AREAS -----------------------------------
C
      Q12E=(XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*YX(I+1,J,K)+XY(I+1,J,K)*YY(I+1,J,K)+
     1 XZ(I+1,J,K)*YZ(I+1,J,K))*FX(I,J,K)
C
      Q12W=(XX(I-1,J,K)*YX(I-1,J,K)+XY(I-1,J,K)*YY(I-1,J,K)+
     1 XZ(I-1,J,K)*YZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*YX(I,J,K)+XY(I,J,K)*YY(I,J,K)+
     1 XZ(I,J,K)*YZ(I,J,K))*FX(I-1,J,K)
C
      Q13E=(XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*(1.-FX(I,J,K))+
     1     (XX(I+1,J,K)*ZX(I+1,J,K)+XY(I+1,J,K)*ZY(I+1,J,K)+
     1 XZ(I+1,J,K)*ZZ(I+1,J,K))*FX(I,J,K)
C
      Q13W=(XX(I-1,J,K)*ZX(I-1,J,K)+XY(I-1,J,K)*ZY(I-1,J,K)+
     1 XZ(I-1,J,K)*ZZ(I-1,J,K))*(1.-FX(I-1,J,K))+
     1     (XX(I,J,K)*ZX(I,J,K)+XY(I,J,K)*ZY(I,J,K)+
     1 XZ(I,J,K)*ZZ(I,J,K))*FX(I-1,J,K)
C
C ------ EAST-WEST CONTRIBUTION ---------------------------------------
C
      SCX=0.5*VISE*Q12E/VOLE*(PHITNE+PHIBNE-PHITSE-PHIBSE)
     1   +0.5*VISE*Q13E/VOLE*(PHITNE+PHITSE-PHIBNE-PHIBSE)
     1   -0.5*VISW*Q12W/VOLW*(PHITNW+PHIBNW-PHITSW-PHIBSW)
     1   -0.5*VISW*Q13W/VOLW*(PHITNW+PHITSW-PHIBNW-PHIBSW)
C
C -------- NORTH-SOUTH PROJECTED AREAS ------------------------------
C
      Q21N=(YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*XX(I,J+1,K)+YY(I,J+1,K)*XY(I,J+1,K)+
     1 YZ(I,J+1,K)*XZ(I,J+1,K))*FY(I,J,K)
C
      Q21S=(YX(I,J-1,K)*XX(I,J-1,K)+YY(I,J-1,K)*XY(I,J-1,K)+
     1 YZ(I,J-1,K)*XZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*XX(I,J,K)+YY(I,J,K)*XY(I,J,K)+
     1 YZ(I,J,K)*XZ(I,J,K))*FY(I,J-1,K)
C
      Q23N=(YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*(1.-FY(I,J,K))+
     1     (YX(I,J+1,K)*ZX(I,J+1,K)+YY(I,J+1,K)*ZY(I,J+1,K)+
     1 YZ(I,J+1,K)*ZZ(I,J+1,K))*FY(I,J,K)
C
      Q23S=(YX(I,J-1,K)*ZX(I,J-1,K)+YY(I,J-1,K)*ZY(I,J-1,K)+
     1 YZ(I,J-1,K)*ZZ(I,J-1,K))*(1.-FY(I,J-1,K))+
     1     (YX(I,J,K)*ZX(I,J,K)+YY(I,J,K)*ZY(I,J,K)+
     1 YZ(I,J,K)*ZZ(I,J,K))*FY(I,J-1,K)
C
C ------------ NORTH-SOUTH CONTRIBUTION -----------------------------
C
      SCY=0.5*VISN*Q21N/VOLN*(PHITNE+PHIBNE-PHITNW-PHIBNW)
     1   +0.5*VISN*Q23N/VOLN*(PHITNE+PHITNW-PHIBNE-PHIBNW)
     1   -0.5*VISS*Q21S/VOLS*(PHITSE+PHIBSE-PHITSW-PHIBSW)
     1   -0.5*VISS*Q23S/VOLS*(PHITSE+PHITSW-PHIBSE-PHIBSW)
C
C ------- TOP-BOTTOM PROJECTED AREAS --------------------------------
C
      Q31T=(ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*XX(I,J,K+1)+ZY(I,J,K+1)*XY(I,J,K+1)+
     1 ZZ(I,J,K+1)*XZ(I,J,K+1))*FZ(I,J,K)
C
      Q31B=(ZX(I,J,K-1)*XX(I,J,K-1)+ZY(I,J,K-1)*XY(I,J,K-1)+
     1 ZZ(I,J,K-1)*XZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*XX(I,J,K)+ZY(I,J,K)*XY(I,J,K)+
     1 ZZ(I,J,K)*XZ(I,J,K))*FZ(I,J,K-1)
C
      Q32T=(ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*(1.-FZ(I,J,K))+
     1     (ZX(I,J,K+1)*YX(I,J,K+1)+ZY(I,J,K+1)*YY(I,J,K+1)+
     1 ZZ(I,J,K+1)*YZ(I,J,K+1))*FZ(I,J,K)
C
      Q32B=(ZX(I,J,K-1)*YX(I,J,K-1)+ZY(I,J,K-1)*YY(I,J,K-1)+
     1 ZZ(I,J,K-1)*YZ(I,J,K-1))*(1.-FZ(I,J,K-1))+
     1     (ZX(I,J,K)*YX(I,J,K)+ZY(I,J,K)*YY(I,J,K)+
     1 ZZ(I,J,K)*YZ(I,J,K))*FZ(I,J,K-1)
C
C -------- TOP-BOTTOM CONTRIBUTION ----------------------------------
C
      SCZ=0.5*VIST*Q31T/VOLT*(PHITNE+PHITSE-PHITNW-PHITSW)
     1   +0.5*VIST*Q32T/VOLT*(PHITNE+PHITNW-PHITSE-PHITSW)
     1   -0.5*VISB*Q31B/VOLB*(PHIBNE+PHIBSE-PHIBNW-PHIBSW)
     1   -0.5*VISB*Q32B/VOLB*(PHIBNE+PHIBNW-PHIBSE-PHIBSW)
C
      SUCD=SCX+SCY+SCZ
      SU(I,J,K)=SU(I,J,K)+FLOAT(1-ITRB)*SUCD
      SU(I,J,K)=SU(I,J,K)+FLOAT(ITRB)*MAX(SUCD,0.)
      SP(I,J,K)=SP(I,J,K)-FLOAT(ITRB)*MAX(-SUCD,0.)/(PHI(I,J,K)+TINY)
C      SU(I,J,K)=SU(I,J,K)+SCX+SCY+SCZ
C
    1 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE BQUICK(IPOS,URFQ,PHI)
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
C
      DIMENSION PHI(NX,NY,NZ),RINW(NX,NY),RINE(NX,NY),RINS(NX,NY),
     1 RINN(NX,NY),RINZ(NX,NY),RLMX(NX),RLMY(NY),RLMZ(NZ),SUQ(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
C
      DO 110 K=1,NK
      DO 110 J=1,NJ
      DO 110 I=1,NI
      SUQ(I,J,K)=0.
      RINW(I,J)=1.
      RINE(I,J)=1.
      RINS(I,J)=1.
      RINN(I,J)=1.
      RINZ(I,J)=1.
  110 CONTINUE
C
C      DO 111 J=2,NJR+2
C      DO 112 I=2,NIR+2
C      RINW(I,J)=0.
C      RINS(I,J)=0.
C      RINN(I,J-1)=0.
C      RINE(I-1,J)=0.
C      RINZ(I,J)=0.
C  112 CONTINUE
C      DO 113 I=NI-3,NI
C      RINS(I-1,J)=0.
C      RINN(I-1,J-1)=0.
C      RINZ(I-1,J)=0.
C      RINE(I-1,J)=0.
C  113 RINW(I,J)=0.
C  111 CONTINUE
C      DO 114 J=NJRN-1,NJM1
C      DO 115 I=2,NIR+2
C      RINS(I,J+1)=0.
C      RINN(I,J)=0.
C      RINZ(I,J)=0.
C      RINW(I,J)=0.
C      RINE(I-1,J)=0.
C  115 CONTINUE
C      DO 116 I=NI-3,NI
C      RINS(I-1,J+1)=0.
C      RINN(I-1,J)=0.
C      RINZ(I-1,J)=0.
C      RINE(I-1,J)=0.
C  116 RINW(I,J)=0.
C  114 CONTINUE
C
      RONE=1.
c      GO TO 8
C---- WEST  FACE ----------------------------------------
      DO 1 K=3,NK-2
      DO 1 J=3,NJ-2
C
      DO 11 I=2,NIM1
C ---------- BOUNDING CRITERIA ----------------------------
      GRMNP=(DSI(I-1,J,K)/(DSI(I-1,J,K)+DSI(I,J,K)))**2
      GRMXP=DSI(I-1,J,K)*(2.*DSI(I,J,K)+DSI(I-1,J,K))/
     1      ((DSI(I-1,J,K)+DSI(I,J,K))**2)
      DNUMP=PHI(I,J,K)-PHI(I-1,J,K)
      DDENP=PHI(I+1,J,K)-PHI(I-1,J,K)
      DDENP=DDENP+SIGN(TINY,DDENP)
      GRDP=DNUMP/DDENP
C   LIMITER 1
      RLIM1P=SIGN(RONE,GRDP-GRMNP)
      RLIM2P=SIGN(RONE,GRMXP-GRDP)
      RLIMP=MAX(RLIM1P,0.)*MAX(RLIM2P,0.)
C  LIMITER 2
      RLIM1P=GRDP/GRMNP
      RLIM2P=(1.-GRDP)/(1.-GRMXP)
      RLIM1P=MAX(RLIM1P,0.)
      RLIM1P=MIN(RLIM1P,1.)
      RLIM2P=MAX(RLIM2P,0.)
      RLIM2P=MIN(RLIM2P,1.)
      RLIMP=RLIM1P*RLIM2P
      RLMX(I)=RLIMP
c      RLIMP=1.
   11 CONTINUE
C
      DO 12 I=4,NI-2
C---------- CELL FACE DENSITY-----------------------
      DENW=(1.-FX(I-1,J,K))*DEN(I-1,J,K)+FX(I-1,J,K)*DEN(I,J,K)
C  ------------------   MASS FLUX  ----------------
      CW=DENW*CX(I-1,J,K)
C
      CWPL=MAX(CW,0.)
      CWMN=-MAX(-CW,0.)
C
C -------- QUICK SOURCE  TERMS ----------------
C 
      ALFP=FX(I-1,J,K)*(DSI(I-2,J,K)+FX(I-1,J,K)*DSI(I-1,J,K))
     1     /(DSI(I-2,J,K)+DSI(I-1,J,K))
      BETP=FX(I-1,J,K)*(DSI(I-1,J,K)**2)*(1.-FX(I-1,J,K))
     1     /(DSI(I-2,J,K)*(DSI(I-2,J,K)+DSI(I-1,J,K)))
      ALFM=(1.-FX(I-1,J,K))*((1.-FX(I-1,J,K))*DSI(I-1,J,K)+
     1     DSI(I,J,K))/((DSI(I-1,J,K)+DSI(I,J,K)))
      BETM=(1.-FX(I-1,J,K))*(DSI(I-1,J,K)**2)*FX(I-1,J,K)
     1     /(DSI(I,J,K)*(DSI(I,J,K)+DSI(I-1,J,K)))
C
C      ALFP=3./8.
C      ALFM=3./8.
C      BETP=1./8.
C      BETM=1./8.
C
      RLIM=MIN(RLMX(I),RLMX(I-1),URFQ)
C
      SUQW=RLIM*CWPL*(ALFP*(PHI(I,J,K)-PHI(I-1,J,K))+BETP*
     1           (PHI(I-1,J,K)-PHI(I-2,J,K)))
     1    +RLIM*CWMN*(ALFM*(PHI(I-1,J,K)-PHI(I,J,K))+BETM*
     1           (PHI(I,J,K)-PHI(I+1,J,K)))
C
      SUQ(I,J,K)=SUQ(I,J,K)+RINW(I,J)*SUQW
C
c      AW(I,J,K)=AW(I,J,K)+RINW(I,J)*RLIM*(-CWPL*ALFP+CWMN*ALFM)
C
   12 CONTINUE
    1 CONTINUE
C
      DO 2 K=3,NK-2
      DO 2 J=3,NJ-2
C ---------- BOUNDING CRITERIA ----------------------------
      DO 21 I=2,NIM1
      GRMNP=(DSI(I-1,J,K)/(DSI(I-1,J,K)+DSI(I,J,K)))**2
      GRMXP=DSI(I-1,J,K)*(2.*DSI(I,J,K)+DSI(I-1,J,K))/
     1      ((DSI(I-1,J,K)+DSI(I,J,K))**2)
      DNUMP=PHI(I,J,K)-PHI(I-1,J,K)
      DDENP=PHI(I+1,J,K)-PHI(I-1,J,K)
      DDENP=DDENP+SIGN(TINY,DDENP)
      GRDP=DNUMP/DDENP
C   LIMITER 1
      RLIM1P=SIGN(RONE,GRDP-GRMNP)
      RLIM2P=SIGN(RONE,GRMXP-GRDP)
      RLIMP=MAX(RLIM1P,0.)*MAX(RLIM2P,0.)
C  LIMITER 2
      RLIM1P=GRDP/GRMNP
      RLIM2P=(1.-GRDP)/(1.-GRMXP)
      RLIM1P=MAX(RLIM1P,0.)
      RLIM1P=MIN(RLIM1P,1.)
      RLIM2P=MAX(RLIM2P,0.)
      RLIM2P=MIN(RLIM2P,1.)
      RLIMP=RLIM1P*RLIM2P
c      RLIMP=1.
      RLMX(I)=RLIMP
   21 CONTINUE
C
      DO 22 I=3,NI-3
C---- CELL FACE DENSITY
      DENE=(1.-FX(I,J,K))*DEN(I,J,K)+FX(I,J,K)*DEN(I+1,J,K)
C----  MASS FLUX
      CE=DENE*CX(I  ,J,K)
C
      CEPL=MAX(CE,0.)
      CEMN=-MAX(-CE,0.)
C
C -------- QUICK SOURCE  TERMS ----------------
C
      ALFP=FX(I,J,K)*(DSI(I-1,J,K)+FX(I,J,K)*DSI(I,J,K))
     1     /(DSI(I-1,J,K)+DSI(I,J,K))
      BETP=FX(I,J,K)*(DSI(I,J,K)**2)*(1.-FX(I,J,K))
     1     /(DSI(I-1,J,K)*(DSI(I-1,J,K)+DSI(I,J,K)))
      ALFM=(1.-FX(I,J,K))*((1.-FX(I,J,K))*DSI(I,J,K)+
     1   DSI(I+1,J,K))/((DSI(I,J,K)+DSI(I+1,J,K)))
      BETM=(1.-FX(I,J,K))*(DSI(I,J,K)**2)*FX(I,J,K)
     1     /(DSI(I+1,J,K)*(DSI(I,J,K)+DSI(I+1,J,K)))
C
C      ALFP=3./8.
C      ALFM=3./8.
C      BETP=1./8.
C      BETM=1./8.
C
      RLIM=MIN(RLMX(I),RLMX(I+1),URFQ)
C
      SUQE=RLIM*CEPL*(ALFP*(PHI(I+1,J,K)-PHI(I,J,K))+BETP*
     1           (PHI(I,J,K)-PHI(I-1,J,K)))
     1    +RLIM*CEMN*(ALFM*(PHI(I,J,K)-PHI(I+1,J,K))+BETM*
     1           (PHI(I+1,J,K)-PHI(I+2,J,K)))
C
      SUQ(I,J,K)=SUQ(I,J,K)-RINE(I,J)*SUQE
C
c      AE(I,J,K)=AE(I,J,K)+RINE(I,J)*RLIM*(-CEPL*ALFP+CEMN*ALFM)
C
   22 CONTINUE
    2 CONTINUE
    8 CONTINUE
c      go to 7
C
      DO 3 K=3,NK-2
      DO 3 I=3,NI-2
C ---------- BOUNDING CRITERIA ----------------------------
      DO 31 J=2,NJM1
      GRMNP=(DSJ(I,J-1,K)/(DSJ(I,J-1,K)+DSJ(I,J,K)))**2
      GRMXP=DSJ(I,J-1,K)*(2.*DSJ(I,J,K)+DSJ(I,J-1,K))/
     1      ((DSJ(I,J-1,K)+DSJ(I,J,K))**2)
      DNUMP=PHI(I,J,K)-PHI(I,J-1,K)
      DDENP=PHI(I,J+1,K)-PHI(I,J-1,K)
      DDENP=DDENP+SIGN(TINY,DDENP)
      GRDP=DNUMP/DDENP
C  LIMITER 1
      RLIM1P=SIGN(RONE,GRDP-GRMNP)
      RLIM2P=SIGN(RONE,GRMXP-GRDP)
      RLIMP=MAX(RLIM1P,0.)*MAX(RLIM2P,0.)
c      grmxp=0.75
c      grmnp=0.25
C  LIMITER 2
      RLIM1P=GRDP/GRMNP
      RLIM2P=(1.-GRDP)/(1.-GRMXP)
      RLIM1P=MAX(RLIM1P,0.)
      RLIM1P=MIN(RLIM1P,1.)
      RLIM2P=MAX(RLIM2P,0.)
      RLIM2P=MIN(RLIM2P,1.)
      RLIMP=RLIM1P*RLIM2P
      RLMY(J)=RLIMP
c      RLIMP=1.
   31 CONTINUE
      DO 32 J=3,NJM1
C ------ CELL FACE DENSITY --------------
      DENS=(1.-FY(I,J-1,K))*DEN(I,J-1,K)+FY(I,J-1,K)*DEN(I,J,K)
C ------- MASS FLUX  --------------------
      CS=DENS*CY(I,J-1,K)
C
      CSPL=MAX(CS,0.)
      CSMN=-MAX(-CS,0.)
C
C -------- QUICK SOURCE  TERMS ----------------
C
      GAMP=FY(I,J-1,K)*(DSJ(I,J-2,K)+FY(I,J-1,K)*DSJ(I,J-1,K))
     1      /(DSJ(I,J-2,K)+DSJ(I,J-1,K))
      DELP=FY(I,J-1,K)*(DSJ(I,J-1,K)**2)*(1.-FY(I,J-1,K))
     1      /(DSJ(I,J-2,K)*(DSJ(I,J-2,K)+DSJ(I,J-1,K)))
      GAMM=(1.-FY(I,J-1,K))*((1.-FY(I,J-1,K))*DSJ(I,J-1,K)+DSJ(I,J,K))
     1      /((DSJ(I,J-1,K)+DSJ(I,J,K)))
      DELM=(1.-FY(I,J-1,K))*(DSJ(I,J-1,K)**2)*FY(I,J-1,K)
     1      /(DSJ(I,J,K)*(DSJ(I,J-1,K)+DSJ(I,J,K)))
C      gamp=3./8.
C      gamm=3./8.
C      delp=1./8.
C      delm=1./8.
C
      RLIM=MIN(RLMY(J),RLMY(J-1),URFQ)
C
      SUQS=RLIM*CSPL*(GAMP*(PHI(I,J,K)-PHI(I,J-1,K))+DELP*
     1           (PHI(I,J-1,K)-PHI(I,J-2,K)))
     1    +RLIM*CSMN*(GAMM*(PHI(I,J-1,K)-PHI(I,J,K))+DELM*
     1           (PHI(I,J,K)-PHI(I,J+1,K)))
C
      SUQ(I,J,K)=SUQ(I,J,K)+RINS(I,J)*SUQS
C
c      AS(I,J,K)=AS(I,J,K)+RINS(I,J)*RLIM*(-CSPL*GAMP+CSMN*GAMM)
C
   32 CONTINUE
    3 CONTINUE
C
      DO 4 K=3,NK-2
      DO 4 I=3,NI-2
C ---------- BOUNDING CRITERIA ----------------------
      DO 41 J=2,NJM1
      GRMNP=(DSJ(I,J-1,K)/(DSJ(I,J-1,K)+DSJ(I,J,K)))**2
      GRMXP=DSJ(I,J-1,K)*(2.*DSJ(I,J,K)+DSJ(I,J-1,K))/
     1      ((DSJ(I,J-1,K)+DSJ(I,J,K))**2)
      DNUMP=PHI(I,J,K)-PHI(I,J-1,K)
      DDENP=PHI(I,J+1,K)-PHI(I,J-1,K)
      DDENP=DDENP+SIGN(TINY,DDENP)
      GRDP=DNUMP/DDENP
C  LIMITER 1
      RLIM1P=SIGN(RONE,GRDP-GRMNP)
      RLIM2P=SIGN(RONE,GRMXP-GRDP)
      RLIMP=MAX(RLIM1P,0.)*MAX(RLIM2P,0.)
c      grmxp=0.75
c      grmnp=0.25
C  LIMITER 2
      RLIM1P=GRDP/GRMNP
      RLIM2P=(1.-GRDP)/(1.-GRMXP)
      RLIM1P=MAX(RLIM1P,0.)
      RLIM1P=MIN(RLIM1P,1.)
      RLIM2P=MAX(RLIM2P,0.)
      RLIM2P=MIN(RLIM2P,1.)
      RLIMP=RLIM1P*RLIM2P
c      RLIMP=1.
      RLMY(J)=RLIMP
   41 CONTINUE
C
      DO 42 J=2,NJ-2
C -----  CELL  FACE  DENSITY ----------------
      DENN=(1.-FY(I,J,K))*DEN(I,J,K)+FY(I,J,K)*DEN(I,J+1,K)
C  ------  MASS FLUX  -----------------------
      CN=DENN*CY(I,J  ,K)
C
      CNPL=MAX(CN,0.)
      CNMN=-MAX(-CN,0.)
C -------- QUICK SOURCE  TERMS ----------------
C
      GAMP=FY(I,J,K)*(DSJ(I,J-1,K)+FY(I,J,K)*DSJ(I,J,K))
     1      /(DSJ(I,J-1,K)+DSJ(I,J,K))
      DELP=FY(I,J,K)*(DSJ(I,J,K)**2)*(1.-FY(I,J,K))
     1      /(DSJ(I,J-1,K)*(DSJ(I,J-1,K)+DSJ(I,J,K)))
      GAMM=(1.-FY(I,J,K))*((1.-FY(I,J,K))*DSJ(I,J,K)+DSJ(I,J+1,K))
     1      /((DSJ(I,J,K)+DSJ(I,J+1,K)))
      DELM=(1.-FY(I,J,K))*(DSJ(I,J,K)**2)*FY(I,J,K)
     1      /(DSJ(I,J+1,K)*(DSJ(I,J,K)+DSJ(I,J+1,K)))
C      gamp=3./8.
C      gamm=3./8.
C      delp=1./8.
C      delm=1./8.
C
      RLIM=MIN(RLMY(J),RLMY(J+1),URFQ)
C
      SUQN=RLIM*CNPL*(GAMP*(PHI(I,J+1,K)-PHI(I,J,K))+DELP*
     1           (PHI(I,J,K)-PHI(I,J-1,K)))
     1    +RLIM*CNMN*(GAMM*(PHI(I,J,K)-PHI(I,J+1,K))+DELM*
     1           (PHI(I,J+1,K)-PHI(I,J+2,K)))
C
      SUQ(I,J,K)=SUQ(I,J,K)-RINN(I,J)*SUQN
C
c      AN(I,J,K)=AN(I,J,K)+RINN(I,J)*RLIM*(-CNPL*GAMP+CNMN*GAMM)
C
   42 CONTINUE
    4 CONTINUE
C
C      GO TO 7
      DO 5 J=3,NJ-2
      DO 5 I=3,NI-2
      DO 51 K=2,NKM1
C ---------- BOUNDING CRITERIA ----------------------------
      GRMNP=(DSK(I,J,K-1)/(DSK(I,J,K-1)+DSK(I,J,K)))**2
      GRMXP=DSK(I,J,K-1)*(2.*DSK(I,J,K)+DSK(I,J,K-1))/
     1      ((DSK(I,J,K-1)+DSK(I,J,K))**2)
      DNUMP=PHI(I,J,K)-PHI(I,J,K-1)
      DDENP=PHI(I,J,K+1)-PHI(I,J,K-1)
      DDENP=DDENP+SIGN(TINY,DDENP)
      GRDP=DNUMP/DDENP
C   LIMITER 1
      RLIM1P=SIGN(RONE,GRDP-GRMNP)
      RLIM2P=SIGN(RONE,GRMXP-GRDP)
      RLIMP=MAX(RLIM1P,0.)*MAX(RLIM2P,0.)
C  LIMITER 2
      RLIM1P=GRDP/GRMNP
      RLIM2P=(1.-GRDP)/(1.-GRMXP)
      RLIM1P=MAX(RLIM1P,0.)
      RLIM1P=MIN(RLIM1P,1.)
      RLIM2P=MAX(RLIM2P,0.)
      RLIM2P=MIN(RLIM2P,1.)
      RLIMP=RLIM1P*RLIM2P
c      RLIMP=1.
      RLMZ(K)=RLIMP
   51 CONTINUE
C
      DO 52 K=3,NKM1
C ----- CELL  FACE   DENSITY   -------------------
      DENB=(1.-FZ(I,J,K-1))*DEN(I,J,K-1)+FZ(I,J,K-1)*DEN(I,J,K)
C  ------  MASS FLUX  ----------------------------
      CB=DENB*CZ(I,J,K-1)
C
      CBPL=MAX(CB,0.)
      CBMN=-MAX(-CB,0.)
C
C -------- QUICK SOURCE  TERMS ----------------
C
      EPSP=FZ(I,J,K-1)*(DSK(I,J,K-2)+FZ(I,J,K-1)*DSK(I,J,K-1))
     1      /(DSK(I,J,K-2)+DSK(I,J,K-1))
      ZETP=FZ(I,J,K-1)*(DSK(I,J,K-1)**2)*(1.-FZ(I,J,K-1))
     1      /(DSK(I,J,K-2)*(DSK(I,J,K-2)+DSK(I,J,K-1)))
      EPSM=(1.-FZ(I,J,K-1))*((1.-FZ(I,J,K-1))*DSK(I,J,K-1)+DSK(I,J,K))
     1      /((DSK(I,J,K-1)+DSK(I,J,K)))
      ZETM=(1.-FZ(I,J,K-1))*(DSK(I,J,K-1)**2)*FZ(I,J,K-1)
     1      /(DSK(I,J,K)*(DSK(I,J,K-1)+DSK(I,J,K)))
C      epsp=3./8.
C      epsm=3./8.
C      zetp=1./8.
C      zetm=1./8.
C
      RLIM=MIN(RLMZ(K),RLMZ(K-1),URFQ)
C
      SUQB=RLIM*CBPL*(EPSP*(PHI(I,J,K)-PHI(I,J,K-1))+ZETP*
     1           (PHI(I,J,K-1)-PHI(I,J,K-2)))
     1    +RLIM*CBMN*(EPSM*(PHI(I,J,K-1)-PHI(I,J,K))+ZETM*
     1           (PHI(I,J,K)-PHI(I,J,K+1)))
C
      SUQ(I,J,K)=SUQ(I,J,K)+RINZ(I,J)*SUQB
C
c      AB(I,J,K)=AB(I,J,K)+RINZ(I,J)*RLIM*(-CBPL*EPSP+CBMN*EPSM)
C
   52 CONTINUE
    5 CONTINUE
C
      DO 6 J=3,NJ-2
      DO 6 I=3,NI-2
C ---------- BOUNDING CRITERIA ----------------------------
      DO 61 K=2,NKM1
      GRMNP=(DSK(I,J,K-1)/(DSK(I,J,K-1)+DSK(I,J,K)))**2
      GRMXP=DSK(I,J,K-1)*(2.*DSK(I,J,K)+DSK(I,J,K-1))/
     1      ((DSK(I,J,K-1)+DSK(I,J,K))**2)
      DNUMP=PHI(I,J,K)-PHI(I,J,K-1)
      DDENP=PHI(I,J,K+1)-PHI(I,J,K-1)
      DDENP=DDENP+SIGN(TINY,DDENP)
      GRDP=DNUMP/DDENP
C   LIMITER 1
      RLIM1P=SIGN(RONE,GRDP-GRMNP)
      RLIM2P=SIGN(RONE,GRMXP-GRDP)
      RLIMP=MAX(RLIM1P,0.)*MAX(RLIM2P,0.)
C  LIMITER 2
      RLIM1P=GRDP/GRMNP
      RLIM2P=(1.-GRDP)/(1.-GRMXP)
      RLIM1P=MAX(RLIM1P,0.)
      RLIM1P=MIN(RLIM1P,1.)
      RLIM2P=MAX(RLIM2P,0.)
      RLIM2P=MIN(RLIM2P,1.)
      RLIMP=RLIM1P*RLIM2P
      RLMZ(K)=RLIMP
   61 CONTINUE
C
c      RLIMP=1.
      DO 62 K=2,NK-2
C  ------  CELL FACE DENSITY   ---------------------
      DENT=(1.-FZ(I,J,K))*DEN(I,J,K)+FZ(I,J,K)*DEN(I,J,K+1)
C
      CT=DENT*CZ(I,J,K  )
C
      CTPL=MAX(CT,0.)
      CTMN=-MAX(-CT,0.)
C
C -------- QUICK SOURCE  TERMS ----------------
C
      EPSP=FZ(I,J,K)*(DSK(I,J,K-1)+FZ(I,J,K)*DSK(I,J,K))
     1      /(DSK(I,J,K-1)+DSK(I,J,K))
      ZETP=FZ(I,J,K)*(DSK(I,J,K)**2)*(1.-FZ(I,J,K))
     1      /(DSK(I,J,K-1)*(DSK(I,J,K-1)+DSK(I,J,K)))
      EPSM=(1.-FZ(I,J,K))*((1.-FZ(I,J,K))*DSK(I,J,K)+DSK(I,J,K+1))
     1      /((DSK(I,J,K)+DSK(I,J,K+1)))
      ZETM=(1.-FZ(I,J,K))*(DSK(I,J,K)**2)*FZ(I,J,K)
     1      /(DSK(I,J,K+1)*(DSK(I,J,K)+DSK(I,J,K+1)))
C      epsp=3./8.
C      epsm=3./8.
C      zetp=1./8.
C      zetm=1./8.
C
      RLIM=MIN(RLMZ(K),RLMZ(K+1),URFQ)
C
      SUQT=RLIM*CTPL*(EPSP*(PHI(I,J,K+1)-PHI(I,J,K))+ZETP*
     1           (PHI(I,J,K)-PHI(I,J,K-1)))
     1    +RLIM*CTMN*(EPSM*(PHI(I,J,K)-PHI(I,J,K+1))+ZETM*
     1           (PHI(I,J,K+1)-PHI(I,J,K+2)))
C
      SUQ(I,J,K)=SUQ(I,J,K)-RINZ(I,J)*SUQT
C
c      AT(I,J,K)=AT(I,J,K)+RINZ(I,J)*RLIM*(-CTPL*EPSP+CTMN*EPSM)
C
   62 CONTINUE
    6 CONTINUE
C
    7 continue
C
      IF(IPOS.EQ.0)THEN
      DO K=2,NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+SUQ(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
      IF(IPOS.EQ.1)THEN
      DO K=2,NKM1
      DO J=2,NJM1
      DO I=2,NIM1
      SU(I,J,K)=SU(I,J,K)+MAX(0.,SUQ(I,J,K))
      SP(I,J,K)=SP(I,J,K)-MAX(0.,-SUQ(I,J,K))/(PHI(I,J,K)+TINY)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
      RETURN
      END
C
	SUBROUTINE NLEVM1
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
C
      COMMON /COEF/ISS(NX),INS(NX),AP(NX,NY,NZ),
     1 AE(NX,NY,NZ),AW(NX,NY,NZ),AN(NX,NY,NZ),AS(NX,NY,NZ),AT(NX,NY,NZ),
     1 AB(NX,NY,NZ),SU(NX,NY,NZ),SP(NX,NY,NZ),
     1 CX(NX,NY,NZ),CY(NX,NY,NZ),CZ(NX,NY,NZ),
     1 CXM(NX,NY,NZ),CYM(NX,NY,NZ),CZM(NX,NY,NZ)
C
      COMMON /VARS/
     1 U(NX,NY,NZ),V(NX,NY,NZ),W(NX,NY,NZ),P(NX,NY,NZ),PP(NX,NY,NZ),
     1 TKE(NX,NY,NZ),TED(NX,NY,NZ),DEN(NX,NY,NZ),VIS(NX,NY,NZ),
     1 PK1(NX,NY,NZ),PK2(NX,NY,NZ),T(NX,NY,NZ),
     1 DU1(NX,NY,NZ),DU2(NX,NY,NZ),DU3(NX,NY,NZ),
     1 DV1(NX,NY,NZ),DV2(NX,NY,NZ),DV3(NX,NY,NZ),
     1 DW1(NX,NY,NZ),DW2(NX,NY,NZ),DW3(NX,NY,NZ),SUPP(NX,NY,NZ)
C
      COMMON /GEOM/
     1 X(NX,NY,NZ),Y(NX,NY,NZ),Z(NX,NY,NZ),XC(NX,NY,NZ),YC(NX,NY,NZ),
     1 ZC(NX,NY,NZ),FX(NX,NY,NZ),FY(NX,NY,NZ),FZ(NX,NY,NZ),
     1 XX(NX,NY,NZ),XY(NX,NY,NZ),XZ(NX,NY,NZ),YX(NX,NY,NZ),
     1 YY(NX,NY,NZ),YZ(NX,NY,NZ),ZX(NX,NY,NZ),ZY(NX,NY,NZ),
     1 ZZ(NX,NY,NZ),XJACOB(NX,NY,NZ),VOL(NX,NY,NZ),XL3(NX,NY,NZ)
     1 ,DSI(NX,NY,NZ),DSJ(NX,NY,NZ),DSK(NX,NY,NZ)
C
      COMMON /STRN/
     1 DUDXC(NX,NY,NZ),DUDYC(NX,NY,NZ),DUDZC(NX,NY,NZ),DVDXC(NX,NY,NZ),
     1 DVDYC(NX,NY,NZ),DVDZC(NX,NY,NZ),DWDXC(NX,NY,NZ),DWDYC(NX,NY,NZ),
     1 DWDZC(NX,NY,NZ),ETRM(NX,NY,NZ),DKDX(NX,NY,NZ),DKDY(NX,NY,NZ),
     1 DKDZ(NX,NY,NZ),EDNW(NX,NY,NZ),DLDXJ2(NX,NY,NZ),ETAL(NX,NY,NZ)
C
      COMMON /TSTRS/
     1 UU(NX,NY,NZ),VV(NX,NY,NZ),WW(NX,NY,NZ),UV(NX,NY,NZ),
     1 UW(NX,NY,NZ),VW(NX,NY,NZ),YPLUS(NX-2,NY-2,NZ-2),
     1 NU(NX,NY,NZ),TBL(NX)
C
      COMMON /WLRFL/
     1 ELE(NX,NY,NZ),COSXE(NX,NY,NZ),COSYE(NX,NY,NZ),COSZE(NX,NY,NZ),
     1 ELW(NX,NY,NZ),COSXW(NX,NY,NZ),COSYW(NX,NY,NZ),COSZW(NX,NY,NZ),
     1 ELN(NX,NY,NZ),COSXN(NX,NY,NZ),COSYN(NX,NY,NZ),COSZN(NX,NY,NZ),
     1 ELS(NX,NY,NZ),COSXS(NX,NY,NZ),COSYS(NX,NY,NZ),COSZS(NX,NY,NZ),
     1 ELT(NX,NY,NZ),COSXT(NX,NY,NZ),COSYT(NX,NY,NZ),COSZT(NX,NY,NZ),
     1 ELB(NX,NY,NZ),COSXB(NX,NY,NZ),COSYB(NX,NY,NZ),COSZB(NX,NY,NZ)
C
      COMMON /LBTIM/
     1    U0(NX,NY,NZ),V0(NX,NY,NZ),T0(NX,NY,NZ),W0(NX,NY,NZ),
     2    TKE0(NX,NY,NZ),TED0(NX,NY,NZ),SUTED0(NX,NY,NZ),
     3    UU0(NX,NY,NZ),VV0(NX,NY,NZ),WW0(NX,NY,NZ),UV0(NX,NY,NZ),
     4    UW0(NX,NY,NZ),VW0(NX,NY,NZ),SUU0(NX,NY,NZ),SUV0(NX,NY,NZ),
     5    SUW0(NX,NY,NZ),SUT0(NX,NY,NZ),SUTKE0(NX,NY,NZ),
     6    SUUU0(NX,NY,NZ),SUVV0(NX,NY,NZ),SUWW0(NX,NY,NZ),
     7    SUUV0(NX,NY,NZ),SUUW0(NX,NY,NZ),SUVW0(NX,NY,NZ)
	 
       COMMON /BSECS/ LSTPS(NSEC),LSTPN(NSEC),LSTPE(NSEC)
     1       ,LSTPW(NSEC)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
	DIMENSION DUIDXJ(3,3),S(3,3),OMG(3,3),UIUJ(3,3),DEL(3,3)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)

	C_1=-0.1
	C_2= 0.1
	C_3= 0.26
	C_5= 0.0
	ETA1=0.0

	DO 301 M=1,3
	DO 301 N=1,3
	DEL(M,N)=0.0
	IF (M.EQ.N) DEL(M,N)=1.0	!KRONECKER DELTA
 301	CONTINUE


 	DO 302 I=2,NI-1
	DO 302 J=2,NJ-1
	DO 302 K=2,NK-1
	VISTUR=VIS(I,J,K)-VISCOS
	DUIDXJ(1,1)=DUDXC(I,J,K)
	DUIDXJ(1,2)=DUDYC(I,J,K)
	DUIDXJ(1,3)=DUDZC(I,J,K)
	DUIDXJ(2,1)=DVDXC(I,J,K)
	DUIDXJ(2,2)=DVDYC(I,J,K)
	DUIDXJ(2,3)=DVDZC(I,J,K)
	DUIDXJ(3,1)=DWDXC(I,J,K)
	DUIDXJ(3,2)=DWDYC(I,J,K)
	DUIDXJ(3,3)=DWDZC(I,J,K)

C	S   =0.0				!STRAIN RATE TENSOR
C	OMG =0.0				!VORTICITY RATE TENSOR
	S2  =0.0
	OMG2=0.0
	DO 303 M=1,3
	DO 303 N=1,3
	S(M,  N)=DUIDXJ(M,N)+DUIDXJ(N,M)
	OMG(M,N)=DUIDXJ(M,N)-DUIDXJ(N,M)
	S2  =S2+S(M,N)*S(M,N)
	OMG2=OMG2+OMG(M,N)*OMG(M,N)
 303	CONTINUE
	STLD(I,J,K)  =TKE(I,J,K)*SQRT(0.5*S2  )/(TED(I,J,K)+tiny)
	OMGTLD(I,J,K)=TKE(I,J,K)*SQRT(0.5*OMG2)/(TED(I,J,K)+tiny)
	ETA1  =AMAX1(STLD(I,J,K),OMGTLD(I,J,K))
	ETAL(I,J,K)=ETA1
	RTLDT =DEN(I,J,K)*(TKE(I,J,K)*TKE(I,J,K))/(TED(I,J,K)+TINY)/VISCOS
	FRS   =0.235*EXP(-RTLDT/400.0)*(AMAX1(0.0,ETA1-3.333))**2
	IF (ICRAFT.EQ.0) THEN
	CMU1(I,J,K)  =0.3*(1.0-EXP(-0.36*EXP(0.75*ETA1)))/
	1(1.0+0.35*ETA1**1.5)
	ELSE
      
!******************	NLEVM3 REFINEMENT	
	
	BETA=0.5
      ALPHA=3.5-BETA*(1.0-EXP(-(RTLDT/400.)**2))
!	CMU1(I,J,K)  =AMIN1(0.09,1.2/(1.0+3.5*ETA1+FRS))
	CMU1(I,J,K)  =AMIN1(0.09,1.2/(1.0+ALPHA*ETA1+FRS))

!******************
	ENDIF
	C_4=-10.0*CMU1(I,J,K)*CMU1(I,J,K)
	C_6=-5.0 *CMU1(I,J,K)*CMU1(I,J,K)
	C_7= 5.0 *CMU1(I,J,K)*CMU1(I,J,K)

	DO 304 M=1,3
	DO 304 N=1,3
	IF (M.LE.N) THEN
	SS    =0.0
	OMGG  =0.0
	TRM_C2=0.0

	DO 305 L=1,3
	SS    =SS+S(M,L)*S(L,N)
	OMGG  =OMGG+OMG(M,L)*OMG(N,L)
	TRM_C2=TRM_C2+OMG(M,L)*S(L,N)+OMG(N,L)*S(L,M)
 305	CONTINUE

	TRM_C4=0.0
	TRM_C5=0.0
	DO 306 L=1,3
	DO 306 Q=1,3
	TRM_C4=TRM_C4+(S(L,M)*OMG(Q,N)+S(L,N)*OMG(Q,M))*S(L,Q)
	TRM_C5=TRM_C5+OMG(M,L)*OMG(L,Q)*S(Q,N)+S(M,L)*OMG(L,Q)*OMG(Q,N)

	DO R=1,3
	TRM_C5=TRM_C5-2.0*S(L,Q)*OMG(Q,R)*OMG(R,L)*DEL(M,N)/3.0
	ENDDO
 306	CONTINUE

	UIUJ(M,N)=2.0*TKE(I,J,K)*DEL(M,N)/3.0-VISTUR*S(M,N)/DEN(I,J,K)+
	1(VISTUR*TKE(I,J,K)/(DEN(I,J,K)*TED(I,J,K)+tiny))*
	1(C_1*(SS-S2*DEL(M,N)/3.0)+C_2*TRM_C2+C_3*(OMGG-OMG2*DEL(M,N)/3.0))
	1+(VISTUR*(TKE(I,J,K)*TKE(I,J,K))/
	1((TED(I,J,K)*TED(I,J,K))*DEN(I,J,K)+tiny))*
	1(C_4*TRM_C4+C_5*TRM_C5+C_6*S(M,N)*S2+C_7*S(M,N)*OMG2)
	ENDIF
 304	CONTINUE
	UU(I,J,K)=UIUJ(1,1)
	UV(I,J,K)=UIUJ(1,2)
	UW(I,J,K)=UIUJ(1,3)
	VV(I,J,K)=UIUJ(2,2)
	VW(I,J,K)=UIUJ(2,3)
	WW(I,J,K)=UIUJ(3,3)
 302	CONTINUE
	RETURN
	END

      SUBROUTINE CHECK_SUB
C      USE DFLIB
      PARAMETER(NX=103,NZ=35,NY=67,NSEC=10 )
      COMMON /PARA/
     1 NI,NJ,NK,NIM1,NJM1,NKM1,GREAT,TINY,IPREF,JPREF,KPREF,
     1 FLOWIN,XMONIN,SORMAX,AU,AD,MAXIT,I90,I180,RATIO,ETA,UIN,
     1 RESORU,RESORV,RESORW,RESORM,RESORK,RESORE,RESORT,
     1 NSWPU,NSWPV,NSWPW,NSWPP,NSWPK,NSWPE,NSWPT,IT,NIRMS,NIRMF,
     1 URFU,URFV,URFW,URFP,URFK,URFE,URFVIS,CA1,CA2,CA1W,CA2W,DPEX,
     1 ISCMU,ISCMK,ISCME,INWP,LWS,LWN,LWT,LWB,INWM,IEVM,IYAP,NREAD,
     1 CMU,C1,C2,CAPPA,ELOG,VISCOS,DENSIT,PRTKE,PRTED,TKEIN,TEDIN,
     1 PRLTM,PRTTM,CT,QWOCP,URFT,TTBLK,FLOWEN,DTR1,DTR2,DTR3,
     1 NJSL,NJSL1,NJNH,NJNH1,NKBL,NKBL1,NKTH,NKTH1,NIR,NJR,NJRN
     1 ,OMX,OMY,OMZ,HX,HY,HZ,LSCN,LSCS,LSCE,LSCW,QS(NX,NZ),QN(NX,NZ)
      COMMON /NLEVM/
     1   STLD(NX,NY,NZ),OMGTLD(NX,NY,NZ),INLEVM,
     1   CMU1(NX,NY,NZ),ICRAFT,RTL(NX,NY,NZ),uvtermL(NX,NY,NZ),
	 1   uvterm1(NX,NY,NZ),uvterm2(NX,NY,NZ),uvterm3(NX,NY,NZ),
	 1   uvterm4(NX,NY,NZ),uvterm6(NX,NY,NZ),uvterm7(NX,NY,NZ),
	 1   uutermL(NX,NY,NZ),uuterm1(NX,NY,NZ),uuterm2(NX,NY,NZ),
	 1   uuterm3(NX,NY,NZ),uuterm4(NX,NY,NZ),uuterm6(NX,NY,NZ),
	 1   uuterm7(NX,NY,NZ),vvtermL(NX,NY,NZ),vvterm1(NX,NY,NZ),
	 1   vvterm2(NX,NY,NZ),vvterm3(NX,NY,NZ),vvterm4(NX,NY,NZ),
	 1   vvterm6(NX,NY,NZ),vvterm7(NX,NY,NZ)
      DIMENSION RCNT(NY),RCNB(NY),RCTN(NZ),RCTS(NZ)
	     
      CHARACTER KEY1
	WRITE(*,   *) "		CHECKING OF THE INITIAL PARAMETERS"
      WRITE(*,   *)
	WRITE(*,1100) " MESH IS:",NI,"*",NJ,"*",NK
	WRITE(*,1101) " IEVM  =",IEVM
	WRITE(*,1101) " INWM  =",INWM
	WRITE(*,1101) " IYAP  =",IYAP
	WRITE(*,1101) " ICRAFT  =",ICRAFT
	WRITE(*,1101) " INLEVM=",INLEVM
	WRITE(*,1101) " INWP  =",INWP
	WRITE(*,1101) " NREAD =",NREAD
	WRITE(*,   *) "UNDER-RELAXATION FACTORS"
	WRITE(*,1102) " URFU  =",URFU
	WRITE(*,1102) " URFV  =",URFV
	WRITE(*,1102) " URFW  =",URFW
	WRITE(*,1102) " URFK  =",URFK
	WRITE(*,1102) " URFE  =",URFE
	WRITE(*,1102) " URFT  =",URFT
	WRITE(*,1102) " URFVIS=",URFVIS
1103	CONTINUE
C	WRITE(*,   *) "DO YOU WANT TO CONTINUE?(Y/N)"
C      KEY1=GETCHARQQ()
C	IF(KEY1.EQ.'y' .OR. KEY1.EQ.'Y') THEN 
C	GOTO 1104
C	ELSE
C	IF(KEY1.EQ.'n' .OR. KEY1.EQ.'N') STOP
C	GOTO 1103
C	ENDIF

 1104	CONTINUE
 1100 FORMAT(A9,I3,A2,I3,A2,I3)
 1101 FORMAT(A8,I2)
 1102 FORMAT(A8,F4.2)

	RETURN
	END