      module msise00_gemini
C     module to avoid link conflicts
      use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
      use msise00_data_gemini, only : parm7g, ptm, pdm, pavgm, imr

      private
      public :: gtd7, tselec, meters

      contains

      SUBROUTINE GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C
C     NRLMSISE-00
C     -----------
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere
C
C        NEW FEATURES:
C          *Extensive satellite drag database used in model generation
C          *Revised O2 (and O) in lower thermosphere
C          *Additional nonlinear solar activity term
C          *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
C           At high altitudes (> 500 km), hot atomic oxygen or ionized
C           oxygen can become appreciable for some ranges of subroutine
C           inputs, thereby affecting drag on satellites and debris. We
C           group these species under the term "anomalous oxygen," since
C           their individual variations are not presently separable with
C           the drag data used to define this model component.
C
C        SUBROUTINES FOR SPECIAL OUTPUTS:
C
C        HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY
C        (SUBROUTINE GTD7D, OUTPUT D(6))
C           For atmospheric drag calculations at altitudes above 500 km,
C           call SUBROUTINE GTD7D to compute the "effective total mass
C           density" by including contributions from "anomalous oxygen."
C           See "NOTES ON OUTPUT VARIABLES" below on D(6).
C
C        PRESSURE GRID (SUBROUTINE GHP7)
C          See subroutine GHP7 to specify outputs at a pressure level
C          rather than at an altitude.
C
C        OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES:
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
C     NOTES ON OUTPUT VARIABLES:
C        TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
C
C        O, H, and N are set to zero below 72.5 km
C
C        T(1), Exospheric temperature, is set to global average for
C        altitudes below 120 km. The 120 km gradient is left at global
C        average value for altitudes below 72 km.
C
C        D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7
C        and GTD7D
C
C          SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
C          species labeled by indices 1-5 and 7-8 in output variable D.
C          This includes He, O, N2, O2, Ar, H, and N but does NOT include
C          anomalous oxygen (species index 9).
C
C          SUBROUTINE GTD7D -- D(6) is the "effective total mass density
C          for drag" and is the sum of the mass densities of all species
C          in this model, INCLUDING anomalous oxygen.
C
C     SWITCHES: The following is for test and special purposes:
C
C        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW),
C        WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
C        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C        FOR THE FOLLOWING VARIATIONS
C               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
C               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
C               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
C               7 - DIURNAL               8 - SEMIDIURNAL
C               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
C              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
C              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
C              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
C              16 - ALL TINF VAR         17 - ALL TLB VAR
C              18 - ALL TN1 VAR           19 - ALL S VAR
C              20 - ALL TN2 VAR           21 - ALL NLB VAR
C              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
C
C        To get current values of SW: CALL TRETRV(SW)
C
      Real,Intent(OUT):: D(9), T(2)
      Real,Intent(In) :: SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP(7)
      Integer,Intent(IN)::IYD,MASS

      Real DS(9),TS(2)
      real ZN3(5),ZN2(4),SV(25)

      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO7/TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)

      COMMON/PARM7g/PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),
     $ PMA(100,10),SAM(100)

      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      SAVE

      DATA MN3/5/,ZN3/32.5,20.,15.,10.,0./
      DATA MN2/4/,ZN2/72.5,55.,45.,32.5/
      DATA ZMIX/62.5/,ALAST/99999./,MSSL/-999/
      DATA SV/25*1./
      IF(ISW.NE.64999) CALL TSELEC(SV)
C
C        Test for changed input
      V1=VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,1)
C       Latitude variation of gravity (none for SW(2)=0)
      XLAT=GLAT
      IF(SW(2).EQ.0) XLAT=45.
      CALL GLATF(XLAT,GSURF,RE)
C
      XMM=PDM(5,3)
C
C       THERMOSPHERE/MESOSPHERE (above ZN2(1))
      ALTT=MAX(ALT,ZN2(1))
      MSS=MASS
C       Only calculate N2 in thermosphere if alt in mixed region
      IF(ALT.LT.ZMIX.AND.MASS.GT.0) MSS=28
C       Only calculate thermosphere if input parameters changed
C         or altitude above ZN2(1) in mesosphere
      IF(V1.EQ.1..OR.ALT.GT.ZN2(1).OR.ALAST.GT.ZN2(1).OR.MSS.NE.MSSL)
     $ THEN
        CALL GTS7(IYD,SEC,ALTT,GLAT,GLONG,STL,F107A,F107,AP,MSS,DS,TS)
        DM28M=DM28
C         metric adjustment
        IF(IMR.EQ.1) DM28M=DM28*1.E6
        MSSL=MSS
      ENDIF
      T(1)=TS(1)
      T(2)=TS(2)
      IF(ALT.GE.ZN2(1)) THEN
        DO 5 J=1,9
          D(J)=DS(J)
    5   CONTINUE
        GOTO 10
      ENDIF
C
C       LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3(1) and ZN2(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
C         Only calculate nodes if input changed
       IF(V1.EQ.1..OR.ALAST.GE.ZN2(1)) THEN
        TGN2(1)=TGN1(2)
        TN2(1)=TN1(5)
        TN2(2)=PMA(1,1)*PAVGM(1)/(1.-SW(20)*GLOB7S(PMA(1,1)))
        TN2(3)=PMA(1,2)*PAVGM(2)/(1.-SW(20)*GLOB7S(PMA(1,2)))
        TN2(4)=PMA(1,3)*PAVGM(3)/(1.-SW(20)*SW(22)*GLOB7S(PMA(1,3)))
        TGN2(2)=PAVGM(9)*PMA(1,10)*(1.+SW(20)*SW(22)*GLOB7S(PMA(1,10)))
     $  *TN2(4)*TN2(4)/(PMA(1,3)*PAVGM(3))**2
        TN3(1)=TN2(4)
       ENDIF
C Including ZN3(1) in the jump condition creates a model coverage gap at that exact altitude
C       IF(ALT.GE.ZN3(1)) GOTO 6
       IF(ALT.GT.ZN3(1)) GOTO 6
C
C       LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
C         Temperature at nodes and gradients at end nodes
C         Inverse temperature a linear function of spherical harmonics
C         Only calculate nodes if input changed
        IF(V1.EQ.1..OR.ALAST.GE.ZN3(1)) THEN
         TGN3(1)=TGN2(2)
         TN3(2)=PMA(1,4)*PAVGM(4)/(1.-SW(22)*GLOB7S(PMA(1,4)))
         TN3(3)=PMA(1,5)*PAVGM(5)/(1.-SW(22)*GLOB7S(PMA(1,5)))
         TN3(4)=PMA(1,6)*PAVGM(6)/(1.-SW(22)*GLOB7S(PMA(1,6)))
         TN3(5)=PMA(1,7)*PAVGM(7)/(1.-SW(22)*GLOB7S(PMA(1,7)))
         TGN3(2)=PMA(1,8)*PAVGM(8)*(1.+SW(22)*GLOB7S(PMA(1,8)))
     $   *TN3(5)*TN3(5)/(PMA(1,7)*PAVGM(7))**2
        ENDIF
    6   CONTINUE
        IF(MASS.EQ.0) GOTO 50
C          LINEAR TRANSITION TO FULL MIXING BELOW ZN2(1)
        DMC=0
        IF(ALT.GT.ZMIX) DMC=1.-(ZN2(1)-ALT)/(ZN2(1)-ZMIX)
        DZ28=DS(3)
C      ***** N2 DENSITY ****
        DMR=DS(3)/DM28M-1.
        D(3)=DENSM(ALT,DM28M,XMM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
        D(3)=D(3)*(1.+DMR*DMC)
C      ***** HE DENSITY ****
        D(1)=0
        IF(MASS.NE.4.AND.MASS.NE.48) GOTO 204
          DMR=DS(1)/(DZ28*PDM(2,1))-1.
          D(1)=D(3)*PDM(2,1)*(1.+DMR*DMC)
  204   CONTINUE
C      **** O DENSITY ****
        D(2)=0
        D(9)=0
  216   CONTINUE
C      ***** O2 DENSITY ****
        D(4)=0
        IF(MASS.NE.32.AND.MASS.NE.48) GOTO 232
          DMR=DS(4)/(DZ28*PDM(2,4))-1.
          D(4)=D(3)*PDM(2,4)*(1.+DMR*DMC)
  232   CONTINUE
C      ***** AR DENSITY ****
        D(5)=0
        IF(MASS.NE.40.AND.MASS.NE.48) GOTO 240
          DMR=DS(5)/(DZ28*PDM(2,5))-1.
          D(5)=D(3)*PDM(2,5)*(1.+DMR*DMC)
  240   CONTINUE
C      ***** HYDROGEN DENSITY ****
        D(7)=0
C      ***** ATOMIC NITROGEN DENSITY ****
        D(8)=0
C
C       TOTAL MASS DENSITY
C
        IF(MASS.EQ.48) THEN
         D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8))
         IF(IMR.EQ.1) D(6)=D(6)/1000.
         ENDIF
         T(2)=TZ
   10 CONTINUE
      GOTO 90
   50 CONTINUE
      DD=DENSM(ALT,1.,0.,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
      T(2)=TZ
   90 CONTINUE
      ALAST=ALT
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTD7D(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,
     $ D,T)
C
C     NRLMSISE-00
C     -----------
C        This subroutine provides Effective Total Mass Density for
C        output D(6) which includes contributions from "anomalous
C        oxygen" which can affect satellite drag above 500 km.  This
C        subroutine is part of the distribution package for the
C        Neutral Atmosphere Empirical Model from the surface to lower
C        exosphere.  See subroutine GTD7 for more extensive comments.
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES:
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3) [includes anomalous oxygen]
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      DIMENSION D(9),T(2),AP(7),DS(9),TS(2)
      COMMON/METSEL/IMR
      CALL GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C       TOTAL MASS DENSITY
C
        IF(MASS.EQ.48) THEN
         D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8)+16.*D(9))
         IF(IMR.EQ.1) D(6)=D(6)/1000.
         ENDIF
      END SUBROUTINE GTD7D
C-----------------------------------------------------------------------
      SUBROUTINE GHP7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,
     $  D,T,PRESS)
C       FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD7
C     INPUT:
C        IYD - YEAR AND DAY AS YYDDD
C        SEC - UT(SEC)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS)
C        F107A - 3 MONTH AVERAGE OF F10.7 FLUX
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
C                    TO CURRENT TIME
C        PRESS - PRESSURE LEVEL(MB)
C     OUTPUT:
C        ALT - ALTITUDE(KM)
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3)
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - HOT O NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      COMMON/PARMB/GSURF,RE
      COMMON/METSEL/IMR
      DIMENSION D(9),T(2),AP(7)
      SAVE
      DATA BM/1.3806E-19/,RGAS/831.4/
      DATA TEST/.00043/,LTEST/12/
      PL=ALOG10(PRESS)
C      Initial altitude estimate
      IF(PL.GE.-5.) THEN
         IF(PL.GT.2.5) ZI=18.06*(3.00-PL)
         IF(PL.GT..75.AND.PL.LE.2.5) ZI=14.98*(3.08-PL)
         IF(PL.GT.-1..AND.PL.LE..75) ZI=17.8*(2.72-PL)
         IF(PL.GT.-2..AND.PL.LE.-1.) ZI=14.28*(3.64-PL)
         IF(PL.GT.-4..AND.PL.LE.-2.) ZI=12.72*(4.32-PL)
         IF(PL.LE.-4.) ZI=25.3*(.11-PL)
         IDAY=MOD(IYD,1000)
         CL=GLAT/90.
         CL2=CL*CL
         IF(IDAY.LT.182) CD=1.-IDAY/91.25
         IF(IDAY.GE.182) CD=IDAY/91.25-3.
         CA=0
         IF(PL.GT.-1.11.AND.PL.LE.-.23) CA=1.0
         IF(PL.GT.-.23) CA=(2.79-PL)/(2.79+.23)
         IF(PL.LE.-1.11.AND.PL.GT.-3.) CA=(-2.93-PL)/(-2.93+1.11)
         Z=ZI-4.87*CL*CD*CA-1.64*CL2*CA+.31*CA*CL
      ENDIF
      IF(PL.LT.-5.) Z=22.*(PL+4.)**2+110
C      ITERATION LOOP
      L=0
   10 CONTINUE
        L=L+1
        CALL GTD7(IYD,SEC,Z,GLAT,GLONG,STL,F107A,F107,AP,48,D,T)
        XN=D(1)+D(2)+D(3)+D(4)+D(5)+D(7)+D(8)
        P=BM*XN*T(2)
        IF(IMR.EQ.1) P=P*1.E-6
        DIFF=PL-ALOG10(P)
        IF(ABS(DIFF).LT.TEST .OR. L.EQ.LTEST) GOTO 20
        XM=D(6)/XN/1.66E-24
        IF(IMR.EQ.1) XM = XM*1.E3
        G=GSURF/(1.+Z/RE)**2
        SH=RGAS*T(2)/(XM*G)
C         New altitude estimate using scale height
        IF(L.LT.6) THEN
          Z=Z-SH*DIFF*2.302
        ELSE
          Z=Z-SH*DIFF
        ENDIF
        GOTO 10
   20 CONTINUE
      IF(L.EQ.LTEST) write(stderr,100) PRESS,DIFF
  100 FORMAT(1X,29HGHP7 NOT CONVERGING FOR PRESS, 1PE12.2,E12.2)
      ALT=Z

      END SUBROUTINE GHP7
C-----------------------------------------------------------------------
      SUBROUTINE GLATF(LAT,GV,REFF)
C      CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
C      RADIUS (REFF)
      REAL LAT
      SAVE
      DATA DGTR/1.74533E-2/
      C2 = COS(2.*DGTR*LAT)
      GV = 980.616*(1.-.0026373*C2)
      REFF = 2.*GV/(3.085462E-6 + 2.27E-9*C2)*1.E-5
      END SUBROUTINE GLATF
C-----------------------------------------------------------------------
      FUNCTION VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,IC)
C       Test if geophysical variables or switches changed and save
C       Return 0 if unchanged and 1 if changed
      DIMENSION AP(7),IYDL(2),SECL(2),GLATL(2),GLL(2),STLL(2)
      DIMENSION FAL(2),FL(2),APL(7,2),SWL(25,2),SWCL(25,2)
      COMMON/CSW/SW(25),ISW,SWC(25)
      SAVE
      DATA IYDL/2*-999/,SECL/2*-999./,GLATL/2*-999./,GLL/2*-999./
      DATA STLL/2*-999./,FAL/2*-999./,FL/2*-999./,APL/14*-999./
      DATA SWL/50*-999./,SWCL/50*-999./
      VTST7=0
      IF(IYD.NE.IYDL(IC)) GOTO 10
      IF(SEC.NE.SECL(IC)) GOTO 10
      IF(GLAT.NE.GLATL(IC)) GOTO 10
      IF(GLONG.NE.GLL(IC)) GOTO 10
      IF(STL.NE.STLL(IC)) GOTO 10
      IF(F107A.NE.FAL(IC)) GOTO 10
      IF(F107.NE.FL(IC)) GOTO 10
      DO 5 I=1,7
        IF(AP(I).NE.APL(I,IC)) GOTO 10
    5 CONTINUE
      DO 7 I=1,25
        IF(SW(I).NE.SWL(I,IC)) GOTO 10
        IF(SWC(I).NE.SWCL(I,IC)) GOTO 10
    7 CONTINUE
      GOTO 20
   10 CONTINUE
      VTST7=1
      IYDL(IC)=IYD
      SECL(IC)=SEC
      GLATL(IC)=GLAT
      GLL(IC)=GLONG
      STLL(IC)=STL
      FAL(IC)=F107A
      FL(IC)=F107
      DO 15 I=1,7
        APL(I,IC)=AP(I)
   15 CONTINUE
      DO 16 I=1,25
        SWL(I,IC)=SW(I)
        SWCL(I,IC)=SWC(I)
   16 CONTINUE
   20 CONTINUE
      END FUNCTION VTST7
C-----------------------------------------------------------------------
      SUBROUTINE GTS7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
C
C     Thermospheric portion of NRLMSISE-00
C     See GTD7 for more extensive comments
C
C        OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
C
C     INPUT VARIABLES:
C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
C              (Year ignored in current model)
C        SEC - UT(SEC)
C        ALT - ALTITUDE(KM) (>72.5 km)
C        GLAT - GEODETIC LATITUDE(DEG)
C        GLONG - GEODETIC LONGITUDE(DEG)
C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
C           - ARRAY CONTAINING:
C             (1) DAILY AP
C             (2) 3 HR AP INDEX FOR CURRENT TIME
C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
C                    TO CURRENT TIME
C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
C                    TO CURRENT TIME
C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                 MASS 17 IS Anomalous O ONLY.)
C
C     NOTES ON INPUT VARIABLES:
C        UT, Local Time, and Longitude are used independently in the
C        model and are not of equal importance for every situation.
C        For the most physically realistic calculation these three
C        variables should be consistent (STL=SEC/3600+GLONG/15).
C        The Equation of Time departures from the above formula
C        for apparent local time can be included if available but
C        are of minor importance.
c
C        F107 and F107A values used to generate the model correspond
C        to the 10.7 cm radio flux at the actual distance of the Earth
C        from the Sun rather than the radio flux at 1 AU. The following
C        site provides both classes of values:
C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
C
C        F107, F107A, and AP effects are neither large nor well
C        established below 80 km and these parameters should be set to
C        150., 150., and 4. respectively.
C
C     OUTPUT VARIABLES:
C        D(1) - HE NUMBER DENSITY(CM-3)
C        D(2) - O NUMBER DENSITY(CM-3)
C        D(3) - N2 NUMBER DENSITY(CM-3)
C        D(4) - O2 NUMBER DENSITY(CM-3)
C        D(5) - AR NUMBER DENSITY(CM-3)
C        D(6) - TOTAL MASS DENSITY(GM/CM3) [Anomalous O NOT included]
C        D(7) - H NUMBER DENSITY(CM-3)
C        D(8) - N NUMBER DENSITY(CM-3)
C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
C        T(1) - EXOSPHERIC TEMPERATURE
C        T(2) - TEMPERATURE AT ALT
C
      DIMENSION ZN1(5),ALPHA(9)
      COMMON/GTS3C/TLB,S,DB04,DB16,DB28,DB32,DB40,DB48,DB01,ZA,T0,Z0
     & ,G0,RL,DD,DB14,TR12
      COMMON/MESO7/TN1(5),TN2(4),TN3(5),TGN1(2),TGN2(2),TGN3(2)
      DIMENSION D(9),T(2),MT(11),AP(7),ALTL(8)

      COMMON/PARM7g/PT(150),PD(150,9),PS(150),PDL(25,2),PTL(100,4),
     $ PMA(100,10),SAM(100)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/TTEST/TINFG,GB,ROUT,TT(15)
      COMMON/DMIX/DM04,DM16,DM28,DM32,DM40,DM01,DM14
      COMMON/METSEL/IMR
      SAVE
      DATA MT/48,0,4,16,28,32,40,1,49,14,17/
      DATA ALTL/200.,300.,160.,250.,240.,450.,320.,450./
      DATA MN1/5/,ZN1/120.,110.,100.,90.,72.5/
      DATA DGTR/1.74533E-2/,DR/1.72142E-2/,ALAST/-999./
      DATA ALPHA/-0.38,0.,0.,0.,0.17,0.,-0.38,0.,0./
C        Test for changed input
      V2=VTST7(IYD,SEC,GLAT,GLONG,STL,F107A,F107,AP,2)
C
      YRD=IYD
      ZA=PDL(16,2)
      ZN1(1)=ZA
      DO 2 J=1,9
        D(J)=0.
    2 CONTINUE
C        TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
      IF(ALT.GT.ZN1(1)) THEN
        IF(V2.EQ.1..OR.ALAST.LE.ZN1(1)) TINF=PTM(1)*PT(1)
     $  *(1.+SW(16)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PT))
      ELSE
        TINF=PTM(1)*PT(1)
      ENDIF
      T(1)=TINF
C          GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
      IF(ALT.GT.ZN1(5)) THEN
        IF(V2.EQ.1.OR.ALAST.LE.ZN1(5)) G0=PTM(4)*PS(1)
     $   *(1.+SW(19)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PS))
      ELSE
        G0=PTM(4)*PS(1)
      ENDIF
C      Calculate these temperatures only if input changed
      IF(V2.EQ.1. .OR. ALT.LT.300.)
     $  TLB=PTM(2)*(1.+SW(17)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,
     $  F107A,F107,AP,PD(1,4)))*PD(1,4)
       S=G0/(TINF-TLB)
C       Lower thermosphere temp variations not significant for
C        density above 300 km
       IF(ALT.LT.300.) THEN
        IF(V2.EQ.1..OR.ALAST.GE.300.) THEN
         TN1(2)=PTM(7)*PTL(1,1)/(1.-SW(18)*GLOB7S(PTL(1,1)))
         TN1(3)=PTM(3)*PTL(1,2)/(1.-SW(18)*GLOB7S(PTL(1,2)))
         TN1(4)=PTM(8)*PTL(1,3)/(1.-SW(18)*GLOB7S(PTL(1,3)))
         TN1(5)=PTM(5)*PTL(1,4)/(1.-SW(18)*SW(20)*GLOB7S(PTL(1,4)))
         TGN1(2)=PTM(9)*PMA(1,9)*(1.+SW(18)*SW(20)*GLOB7S(PMA(1,9)))
     $   *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
        ENDIF
       ELSE
        TN1(2)=PTM(7)*PTL(1,1)
        TN1(3)=PTM(3)*PTL(1,2)
        TN1(4)=PTM(8)*PTL(1,3)
        TN1(5)=PTM(5)*PTL(1,4)
        TGN1(2)=PTM(9)*PMA(1,9)
     $  *TN1(5)*TN1(5)/(PTM(5)*PTL(1,4))**2
       ENDIF
C
      Z0=ZN1(4)
      T0=TN1(4)
      TR12=1.
C
      IF(MASS.EQ.0) GO TO 50
C       N2 variation factor at Zlb
      G28=SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,
     & AP,PD(1,3))
      DAY=AMOD(YRD,1000.)
C        VARIATION OF TURBOPAUSE HEIGHT
      ZHF=PDL(25,2)
     $    *(1.+SW(5)*PDL(25,1)*SIN(DGTR*GLAT)*COS(DR*(DAY-PT(14))))
      YRD=IYD
      T(1)=TINF
      XMM=PDM(5,3)
      Z=ALT
C
      DO 10 J = 1,11
      IF(MASS.EQ.MT(J))   GO TO 15
   10 CONTINUE
      WRITE(6,100) MASS
      GO TO 90
   15 IF(Z.GT.ALTL(6).AND.MASS.NE.28.AND.MASS.NE.48) GO TO 17
C
C       **** N2 DENSITY ****
C
C      Diffusive density at Zlb
      DB28 = PDM(1,3)*EXP(G28)*PD(1,3)
C      Diffusive density at Alt
      D(3)=DENSU(Z,DB28,TINF,TLB, 28.,ALPHA(3),T(2),PTM(6),S,MN1,ZN1,
     & TN1,TGN1)
      DD=D(3)
C      Turbopause
      ZH28=PDM(3,3)*ZHF
      ZHM28=PDM(4,3)*PDL(6,2)
      XMD=28.-XMM
C      Mixed density at Zlb
      B28=DENSU(ZH28,DB28,TINF,TLB,XMD,ALPHA(3)-1.,TZ,PTM(6),S,MN1,
     & ZN1,TN1,TGN1)
      IF(Z.GT.ALTL(3).OR.SW(15).EQ.0.) GO TO 17
C      Mixed density at Alt
      DM28=DENSU(Z,B28,TINF,TLB,XMM,ALPHA(3),TZ,PTM(6),S,MN1,
     & ZN1,TN1,TGN1)
C      Net density at Alt
      D(3)=DNET(D(3),DM28,ZHM28,XMM,28.)
   17 CONTINUE
      GO TO (20,50,20,25,90,35,40,45,25,48,46),  J
   20 CONTINUE
C
C       **** HE DENSITY ****
C
C       Density variation factor at Zlb
      G4 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,1))
C      Diffusive density at Zlb
      DB04 = PDM(1,1)*EXP(G4)*PD(1,1)
C      Diffusive density at Alt
      D(1)=DENSU(Z,DB04,TINF,TLB, 4.,ALPHA(1),T(2),PTM(6),S,MN1,ZN1,
     & TN1,TGN1)
      DD=D(1)
      IF(Z.GT.ALTL(1).OR.SW(15).EQ.0.) GO TO 24
C      Turbopause
      ZH04=PDM(3,1)
C      Mixed density at Zlb
      B04=DENSU(ZH04,DB04,TINF,TLB,4.-XMM,ALPHA(1)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM04=DENSU(Z,B04,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM04=ZHM28
C      Net density at Alt
      D(1)=DNET(D(1),DM04,ZHM04,XMM,4.)
C      Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,1)/B04)
      ZC04=PDM(5,1)*PDL(1,2)
      HC04=PDM(6,1)*PDL(2,2)
C      Net density corrected at Alt
      D(1)=D(1)*CCOR(Z,RL,HC04,ZC04)
   24 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   25 CONTINUE
C
C      **** O DENSITY ****
C
C       Density variation factor at Zlb
      G16= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,2))
C      Diffusive density at Zlb
      DB16 =  PDM(1,2)*EXP(G16)*PD(1,2)
C       Diffusive density at Alt
      D(2)=DENSU(Z,DB16,TINF,TLB, 16.,ALPHA(2),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(2)
      IF(Z.GT.ALTL(2).OR.SW(15).EQ.0.) GO TO 34
C  Corrected from PDM(3,1) to PDM(3,2)  12/2/85
C       Turbopause
      ZH16=PDM(3,2)
C      Mixed density at Zlb
      B16=DENSU(ZH16,DB16,TINF,TLB,16-XMM,ALPHA(2)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM16=DENSU(Z,B16,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM16=ZHM28
C      Net density at Alt
      D(2)=DNET(D(2),DM16,ZHM16,XMM,16.)
C   3/16/99 Change form to match O2 departure from diff equil near 150
C   km and add dependence on F10.7
C      RL=ALOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
      RL=PDM(2,2)*PDL(17,2)*(1.+SW(1)*PDL(24,1)*(F107A-150.))
      HC16=PDM(6,2)*PDL(4,2)
      ZC16=PDM(5,2)*PDL(3,2)
      HC216=PDM(6,2)*PDL(5,2)
      D(2)=D(2)*CCOR2(Z,RL,HC16,ZC16,HC216)
C       Chemistry correction
      HCC16=PDM(8,2)*PDL(14,2)
      ZCC16=PDM(7,2)*PDL(13,2)
      RC16=PDM(4,2)*PDL(15,2)
C      Net density corrected at Alt
      D(2)=D(2)*CCOR(Z,RC16,HCC16,ZCC16)
   34 CONTINUE
      IF(MASS.NE.48.AND.MASS.NE.49) GO TO 90
   35 CONTINUE
C
C       **** O2 DENSITY ****
C
C       Density variation factor at Zlb
      G32= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,5))
C      Diffusive density at Zlb
      DB32 = PDM(1,4)*EXP(G32)*PD(1,5)
C       Diffusive density at Alt
      D(4)=DENSU(Z,DB32,TINF,TLB, 32.,ALPHA(4),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      IF(MASS.EQ.49) THEN
         DD=DD+2.*D(4)
      ELSE
         DD=D(4)
      ENDIF
      IF(SW(15).EQ.0.) GO TO 39
      IF(Z.GT.ALTL(4)) GO TO 38
C       Turbopause
      ZH32=PDM(3,4)
C      Mixed density at Zlb
      B32=DENSU(ZH32,DB32,TINF,TLB,32.-XMM,ALPHA(4)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM32=DENSU(Z,B32,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM32=ZHM28
C      Net density at Alt
      D(4)=DNET(D(4),DM32,ZHM32,XMM,32.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,4)/B32)
      HC32=PDM(6,4)*PDL(8,2)
      ZC32=PDM(5,4)*PDL(7,2)
      D(4)=D(4)*CCOR(Z,RL,HC32,ZC32)
   38 CONTINUE
C      Correction for general departure from diffusive equilibrium above Zlb
      HCC32=PDM(8,4)*PDL(23,2)
      HCC232=PDM(8,4)*PDL(23,1)
      ZCC32=PDM(7,4)*PDL(22,2)
      RC32=PDM(4,4)*PDL(24,2)*(1.+SW(1)*PDL(24,1)*(F107A-150.))
C      Net density corrected at Alt
      D(4)=D(4)*CCOR2(Z,RC32,HCC32,ZCC32,HCC232)
   39 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   40 CONTINUE
C
C       **** AR DENSITY ****
C
C       Density variation factor at Zlb
      G40= SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,6))
C      Diffusive density at Zlb
      DB40 = PDM(1,5)*EXP(G40)*PD(1,6)
C       Diffusive density at Alt
      D(5)=DENSU(Z,DB40,TINF,TLB, 40.,ALPHA(5),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(5)
      IF(Z.GT.ALTL(5).OR.SW(15).EQ.0.) GO TO 44
C       Turbopause
      ZH40=PDM(3,5)
C      Mixed density at Zlb
      B40=DENSU(ZH40,DB40,TINF,TLB,40.-XMM,ALPHA(5)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM40=DENSU(Z,B40,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM40=ZHM28
C      Net density at Alt
      D(5)=DNET(D(5),DM40,ZHM40,XMM,40.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,5)/B40)
      HC40=PDM(6,5)*PDL(10,2)
      ZC40=PDM(5,5)*PDL(9,2)
C      Net density corrected at Alt
      D(5)=D(5)*CCOR(Z,RL,HC40,ZC40)
   44 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   45 CONTINUE
C
C        **** HYDROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G1 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,7))
C      Diffusive density at Zlb
      DB01 = PDM(1,6)*EXP(G1)*PD(1,7)
C       Diffusive density at Alt
      D(7)=DENSU(Z,DB01,TINF,TLB,1.,ALPHA(7),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(7)
      IF(Z.GT.ALTL(7).OR.SW(15).EQ.0.) GO TO 47
C       Turbopause
      ZH01=PDM(3,6)
C      Mixed density at Zlb
      B01=DENSU(ZH01,DB01,TINF,TLB,1.-XMM,ALPHA(7)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM01=DENSU(Z,B01,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM01=ZHM28
C      Net density at Alt
      D(7)=DNET(D(7),DM01,ZHM01,XMM,1.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,6)*ABS(PDL(18,2))/B01)
      HC01=PDM(6,6)*PDL(12,2)
      ZC01=PDM(5,6)*PDL(11,2)
      D(7)=D(7)*CCOR(Z,RL,HC01,ZC01)
C       Chemistry correction
      HCC01=PDM(8,6)*PDL(20,2)
      ZCC01=PDM(7,6)*PDL(19,2)
      RC01=PDM(4,6)*PDL(21,2)
C      Net density corrected at Alt
      D(7)=D(7)*CCOR(Z,RC01,HCC01,ZCC01)
   47 CONTINUE
      IF(MASS.NE.48)   GO TO 90
   48 CONTINUE
C
C        **** ATOMIC NITROGEN DENSITY ****
C
C       Density variation factor at Zlb
      G14 = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,8))
C      Diffusive density at Zlb
      DB14 = PDM(1,7)*EXP(G14)*PD(1,8)
C       Diffusive density at Alt
      D(8)=DENSU(Z,DB14,TINF,TLB,14.,ALPHA(8),T(2),PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      DD=D(8)
      IF(Z.GT.ALTL(8).OR.SW(15).EQ.0.) GO TO 49
C       Turbopause
      ZH14=PDM(3,7)
C      Mixed density at Zlb
      B14=DENSU(ZH14,DB14,TINF,TLB,14.-XMM,ALPHA(8)-1.,
     $  T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
C      Mixed density at Alt
      DM14=DENSU(Z,B14,TINF,TLB,XMM,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
      ZHM14=ZHM28
C      Net density at Alt
      D(8)=DNET(D(8),DM14,ZHM14,XMM,14.)
C       Correction to specified mixing ratio at ground
      RL=ALOG(B28*PDM(2,7)*ABS(PDL(3,1))/B14)
      HC14=PDM(6,7)*PDL(2,1)
      ZC14=PDM(5,7)*PDL(1,1)
      D(8)=D(8)*CCOR(Z,RL,HC14,ZC14)
C       Chemistry correction
      HCC14=PDM(8,7)*PDL(5,1)
      ZCC14=PDM(7,7)*PDL(4,1)
      RC14=PDM(4,7)*PDL(6,1)
C      Net density corrected at Alt
      D(8)=D(8)*CCOR(Z,RC14,HCC14,ZCC14)
   49 CONTINUE
      IF(MASS.NE.48) GO TO 90
   46 CONTINUE
C
C        **** Anomalous OXYGEN DENSITY ****
C
      G16H = SW(21)*GLOBE7(YRD,SEC,GLAT,GLONG,STL,F107A,F107,AP,PD(1,9))
      DB16H = PDM(1,8)*EXP(G16H)*PD(1,9)
      THO=PDM(10,8)*PDL(7,1)
      DD=DENSU(Z,DB16H,THO,THO,16.,ALPHA(9),T2,PTM(6),S,MN1,
     $ ZN1,TN1,TGN1)
      ZSHT=PDM(6,8)
      ZMHO=PDM(5,8)
      ZSHO=SCALH(ZMHO,16.,THO)
      D(9)=DD*EXP(-ZSHT/ZSHO*(EXP(-(Z-ZMHO)/ZSHT)-1.))
      IF(MASS.NE.48) GO TO 90
C
C       TOTAL MASS DENSITY
C
      D(6) = 1.66E-24*(4.*D(1)+16.*D(2)+28.*D(3)+32.*D(4)+40.*D(5)+
     &       D(7)+14.*D(8))
      DB48=1.66E-24*(4.*DB04+16.*DB16+28.*DB28+32.*DB32+40.*DB40+DB01+
     &        14.*DB14)
      GO TO 90
C       TEMPERATURE AT ALTITUDE
   50 CONTINUE
      Z=ABS(ALT)
      DDUM  = DENSU(Z,1., TINF,TLB,0.,0.,T(2),PTM(6),S,MN1,ZN1,TN1,TGN1)
   90 CONTINUE
C       ADJUST DENSITIES FROM CGS TO KGM
      IF(IMR.EQ.1) THEN
        DO 95 I=1,9
          D(I)=D(I)*1.E6
   95   CONTINUE
        D(6)=D(6)/1000.
      ENDIF
      ALAST=ALT
      RETURN
  100 FORMAT(1X,'MASS', I5, '  NOT VALID')
      END Subroutine GTS7
C-----------------------------------------------------------------------
      SUBROUTINE METERS(METER)
      implicit none (external)
C      Convert outputs to Kg & Meters if METER true
      logical,Intent(In) :: METER
      Integer IMR
      COMMON/METSEL/IMR
      SAVE
      IMR=0
      IF(METER) IMR=1
      END SUBROUTINE METERS
C-----------------------------------------------------------------------
      FUNCTION SCALH(ALT,XM,TEMP)
C      Calculate scale height (km)
      COMMON/PARMB/GSURF,RE
      SAVE
      DATA RGAS/831.4/
      G=GSURF/(1.+ALT/RE)**2
      SCALH=RGAS*TEMP/(G*XM)
      END FUNCTION SCALH
C-----------------------------------------------------------------------
      FUNCTION GLOBE7(YRD,SEC,LAT,LONG,TLOC,F107A,F107,AP,P)
C       CALCULATE G(L) FUNCTION
C       Upper Thermosphere Parameters
      REAL LAT, LONG
      DIMENSION P(150),SV(25),AP(7)
      COMMON/TTEST/TINF,GB,ROUT,T(15)
      COMMON/CSW/SW(25),ISW,SWC(25)
      COMMON/LPOLY/PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ IYR,DAY,DF,DFA,APD,APDF,APT(4),XLONG
      SAVE
      DATA DGTR/1.74533E-2/,DR/1.72142E-2/, XL/1000./,TLL/1000./
      DATA SW9/1./,DAYL/-1./,P14/-1000./,P18/-1000./,P32/-1000./
      DATA HR/.2618/,SR/7.2722E-5/,SV/25*1./,NSW/14/,P39/-1000./
C       3hr Magnetic activity functions
C      Eq. A24d
      G0(A)=(A-4.+(P(26)-1.)*(A-4.+(EXP(-ABS(P(25))*(A-4.))-1.)/ABS(P(25
     *))))
C       Eq. A24c
       SUMEX(EX)=1.+(1.-EX**19)/(1.-EX)*EX**(.5)
C       Eq. A24a
      SG0(EX)=(G0(AP(2))+(G0(AP(3))*EX+G0(AP(4))*EX*EX+G0(AP(5))*EX**3
     $ +(G0(AP(6))*EX**4+G0(AP(7))*EX**12)*(1.-EX**8)/(1.-EX))
     $ )/SUMEX(EX)
      IF(ISW.NE.64999) CALL TSELEC(SV)
      DO 10 J=1,14
       T(J)=0
   10 CONTINUE
      IF(SW(9).GT.0) SW9=1.
      IF(SW(9).LT.0) SW9=-1.
      IYR = YRD/1000.
      DAY = YRD - IYR*1000.
      XLONG=LONG
C      Eq. A22 (remainder of code)
      IF(XL.EQ.LAT)   GO TO 15
C          CALCULATE LEGENDRE POLYNOMIALS
      C = SIN(LAT*DGTR)
      S = COS(LAT*DGTR)
      C2 = C*C
      C4 = C2*C2
      S2 = S*S
      PLG(2,1) = C
      PLG(3,1) = 0.5*(3.*C2 -1.)
      PLG(4,1) = 0.5*(5.*C*C2-3.*C)
      PLG(5,1) = (35.*C4 - 30.*C2 + 3.)/8.
      PLG(6,1) = (63.*C2*C2*C - 70.*C2*C + 15.*C)/8.
      PLG(7,1) = (11.*C*PLG(6,1) - 5.*PLG(5,1))/6.
C     PLG(8,1) = (13.*C*PLG(7,1) - 6.*PLG(6,1))/7.
      PLG(2,2) = S
      PLG(3,2) = 3.*C*S
      PLG(4,2) = 1.5*(5.*C2-1.)*S
      PLG(5,2) = 2.5*(7.*C2*C-3.*C)*S
      PLG(6,2) = 1.875*(21.*C4 - 14.*C2 +1.)*S
      PLG(7,2) = (11.*C*PLG(6,2)-6.*PLG(5,2))/5.
C     PLG(8,2) = (13.*C*PLG(7,2)-7.*PLG(6,2))/6.
C     PLG(9,2) = (15.*C*PLG(8,2)-8.*PLG(7,2))/7.
      PLG(3,3) = 3.*S2
      PLG(4,3) = 15.*S2*C
      PLG(5,3) = 7.5*(7.*C2 -1.)*S2
      PLG(6,3) = 3.*C*PLG(5,3)-2.*PLG(4,3)
      PLG(7,3)=(11.*C*PLG(6,3)-7.*PLG(5,3))/4.
      PLG(8,3)=(13.*C*PLG(7,3)-8.*PLG(6,3))/5.
      PLG(4,4) = 15.*S2*S
      PLG(5,4) = 105.*S2*S*C
      PLG(6,4)=(9.*C*PLG(5,4)-7.*PLG(4,4))/2.
      PLG(7,4)=(11.*C*PLG(6,4)-8.*PLG(5,4))/3.
      XL=LAT
   15 CONTINUE
      IF(TLL.EQ.TLOC)   GO TO 16
      IF(SW(7).EQ.0.AND.SW(8).EQ.0.AND.SW(14).EQ.0) GOTO 16
      STLOC = SIN(HR*TLOC)
      CTLOC = COS(HR*TLOC)
      S2TLOC = SIN(2.*HR*TLOC)
      C2TLOC = COS(2.*HR*TLOC)
      S3TLOC = SIN(3.*HR*TLOC)
      C3TLOC = COS(3.*HR*TLOC)
      TLL = TLOC
   16 CONTINUE
      IF(DAY.NE.DAYL.OR.P(14).NE.P14) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P(18).NE.P18) CD18=COS(2.*DR*(DAY-P(18)))
      IF(DAY.NE.DAYL.OR.P(32).NE.P32) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P(39).NE.P39) CD39=COS(2.*DR*(DAY-P(39)))
      DAYL = DAY
      P14 = P(14)
      P18 = P(18)
      P32 = P(32)
      P39 = P(39)
C         F10.7 EFFECT
      DF = F107 - F107A
      DFA=F107A-150.
      T(1) =  P(20)*DF*(1.+P(60)*DFA) + P(21)*DF*DF + P(22)*DFA
     & + P(30)*DFA**2
      F1 = 1. + (P(48)*DFA +P(20)*DF+P(21)*DF*DF)*SWC(1)
      F2 = 1. + (P(50)*DFA+P(20)*DF+P(21)*DF*DF)*SWC(1)
C        TIME INDEPENDENT
      T(2) =
     1  (P(2)*PLG(3,1) + P(3)*PLG(5,1)+P(23)*PLG(7,1))
     & +(P(15)*PLG(3,1))*DFA*SWC(1)
     2 +P(27)*PLG(2,1)
C        SYMMETRICAL ANNUAL
      T(3) =
     1 (P(19) )*CD32
C        SYMMETRICAL SEMIANNUAL
      T(4) =
     1 (P(16)+P(17)*PLG(3,1))*CD18
C        ASYMMETRICAL ANNUAL
      T(5) =  F1*
     1  (P(10)*PLG(2,1)+P(11)*PLG(4,1))*CD14
C         ASYMMETRICAL SEMIANNUAL
      T(6) =    P(38)*PLG(2,1)*CD39
C        DIURNAL
      IF(SW(7).EQ.0) GOTO 200
      T71 = (P(12)*PLG(3,2))*CD14*SWC(5)
      T72 = (P(13)*PLG(3,2))*CD14*SWC(5)
      T(7) = F2*
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2) + P(28)*PLG(6,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2) +P(29)*PLG(6,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5)
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
      T(8) = F2*
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0) GOTO 220
      T(14) = F2*
     1 ((P(40)*PLG(4,4)+(P(94)*PLG(5,4)+P(47)*PLG(7,4))*CD14*SWC(5))*
     $ S3TLOC
     2 +(P(41)*PLG(4,4)+(P(95)*PLG(5,4)+P(49)*PLG(7,4))*CD14*SWC(5))*
     $ C3TLOC)
  220 CONTINUE
C          MAGNETIC ACTIVITY BASED ON DAILY AP

      IF(SW9.EQ.-1.) GO TO 30
      APD=(AP(1)-4.)
      P44=P(44)
      P45=P(45)
      IF(P44.LT.0) P44=1.E-5
      APDF = APD+(P45-1.)*(APD+(EXP(-P44  *APD)-1.)/P44)
      IF(SW(9).EQ.0) GOTO 40
      T(9)=APDF*(P(33)+P(46)*PLG(3,1)+P(35)*PLG(5,1)+
     $ (P(101)*PLG(2,1)+P(102)*PLG(4,1)+P(103)*PLG(6,1))*CD14*SWC(5)+
     $ (P(122)*PLG(2,2)+P(123)*PLG(4,2)+P(124)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(125))))
      GO TO 40
   30 CONTINUE
      IF(P(52).EQ.0) GO TO 40
      EXP1 = EXP(-10800.*ABS(P(52))/(1.+P(139)*(45.-ABS(LAT))))
      IF(EXP1.GT..99999) EXP1=.99999
      IF(P(25).LT.1.E-4) P(25)=1.E-4
      APT(1)=SG0(EXP1)
C      APT(2)=SG2(EXP1)
c      APT(3)=SG0(EXP2)
C      APT(4)=SG2(EXP2)
      IF(SW(9).EQ.0) GOTO 40
      T(9) = APT(1)*(P(51)+P(97)*PLG(3,1)+P(55)*PLG(5,1)+
     $ (P(126)*PLG(2,1)+P(127)*PLG(4,1)+P(128)*PLG(6,1))*CD14*SWC(5)+
     $ (P(129)*PLG(2,2)+P(130)*PLG(4,2)+P(131)*PLG(6,2))*SWC(7)*
     $ COS(HR*(TLOC-P(132))))
  40  CONTINUE
      IF(SW(10).EQ.0.OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      IF(SW(11).EQ.0) GOTO 230
      T(11)= (1.+P(81)*DFA*SWC(1))*
     $((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $ +P(104)*PLG(2,2)+P(105)*PLG(4,2)+P(106)*PLG(6,2)
     $ +SWC(5)*(P(110)*PLG(2,2)+P(111)*PLG(4,2)+P(112)*PLG(6,2))*CD14)*
     $     COS(DGTR*LONG)
     $ +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $ +P(107)*PLG(2,2)+P(108)*PLG(4,2)+P(109)*PLG(6,2)
     $ +SWC(5)*(P(113)*PLG(2,2)+P(114)*PLG(4,2)+P(115)*PLG(6,2))*CD14)*
     $  SIN(DGTR*LONG))
  230 CONTINUE
C        UT AND MIXED UT,LONGITUDE
      IF(SW(12).EQ.0) GOTO 240
      T(12)=(1.+P(96)*PLG(2,1))*(1.+P(82)*DFA*SWC(1))*
     $(1.+P(120)*PLG(2,1)*SWC(5)*CD14)*
     $((P(69)*PLG(2,1)+P(70)*PLG(4,1)+P(71)*PLG(6,1))*
     $     COS(SR*(SEC-P(72))))
      T(12)=T(12)+SWC(11)*
     $ (P(77)*PLG(4,3)+P(78)*PLG(6,3)+P(79)*PLG(8,3))*
     $     COS(SR*(SEC-P(80))+2.*DGTR*LONG)*(1.+P(138)*DFA*SWC(1))
  240 CONTINUE
C        UT,LONGITUDE MAGNETIC ACTIVITY
      IF(SW(13).EQ.0) GOTO 48
      IF(SW9.EQ.-1.) GO TO 45
      T(13)= APDF*SWC(11)*(1.+P(121)*PLG(2,1))*
     $((P( 61)*PLG(3,2)+P( 62)*PLG(5,2)+P( 63)*PLG(7,2))*
     $     COS(DGTR*(LONG-P( 64))))
     $ +APDF*SWC(11)*SWC(5)*
     $ (P(116)*PLG(2,2)+P(117)*PLG(4,2)+P(118)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(119)))
     $ + APDF*SWC(12)*
     $ (P( 84)*PLG(2,1)+P( 85)*PLG(4,1)+P( 86)*PLG(6,1))*
     $     COS(SR*(SEC-P( 76)))
      GOTO 48
   45 CONTINUE
      IF(P(52).EQ.0) GOTO 48
      T(13)=APT(1)*SWC(11)*(1.+P(133)*PLG(2,1))*
     $((P(53)*PLG(3,2)+P(99)*PLG(5,2)+P(68)*PLG(7,2))*
     $     COS(DGTR*(LONG-P(98))))
     $ +APT(1)*SWC(11)*SWC(5)*
     $ (P(134)*PLG(2,2)+P(135)*PLG(4,2)+P(136)*PLG(6,2))*
     $     CD14*COS(DGTR*(LONG-P(137)))
     $ +APT(1)*SWC(12)*
     $ (P(56)*PLG(2,1)+P(57)*PLG(4,1)+P(58)*PLG(6,1))*
     $     COS(SR*(SEC-P(59)))
   48 CONTINUE
C  PARMS NOT USED: 83, 90,100,140-150
   49 CONTINUE
      TINF=P(31)
      DO 50 I = 1,NSW
   50 TINF = TINF + ABS(SW(I))*T(I)
      GLOBE7 = TINF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TSELEC(SV)
C        SET SWITCHES
C        Output in  COMMON/CSW/SW(25),ISW,SWC(25)
C        SW FOR MAIN TERMS, SWC FOR CROSS TERMS
C
C        TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SV),
C        WHERE SV IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
C        FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
C
C        To get current values of SW: CALL TRETRV(SW)
C
      DIMENSION SV(25),SAV(25),SVV(25)
      COMMON/CSW/SW(25),ISW,SWC(25)
      SAVE
      DO 100 I = 1,25
        SAV(I)=SV(I)
        SW(I)=AMOD(SV(I),2.)
        IF(ABS(SV(I)).EQ.1.OR.ABS(SV(I)).EQ.2.) THEN
          SWC(I)=1.
        ELSE
          SWC(I)=0.
        ENDIF
  100 CONTINUE
      ISW=64999
      RETURN
!      ENTRY TRETRV(SVV)
      DO 200 I=1,25
        SVV(I)=SAV(I)
  200 CONTINUE
      END
C-----------------------------------------------------------------------
      FUNCTION GLOB7S(P)
C      VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
      REAL LONG
      COMMON/LPOLY/PLG(9,4),CTLOC,STLOC,C2TLOC,S2TLOC,C3TLOC,S3TLOC,
     $ IYR,DAY,DF,DFA,APD,APDF,APT(4),LONG
      COMMON/CSW/SW(25),ISW,SWC(25)
      DIMENSION P(*),T(14)
      SAVE
      DATA DR/1.72142E-2/,DGTR/1.74533E-2/,PSET/2./
      DATA DAYL/-1./,P32,P18,P14,P39/4*-1000./
C       CONFIRM PARAMETER SET
      IF(P(100).EQ.0) P(100)=PSET
      IF(P(100).NE.PSET) THEN
        write(stderr,900) PSET,P(100)
  900   FORMAT(1X,'WRONG PARAMETER SET FOR GLOB7S',3F10.1)
        error stop
      ENDIF
      DO 10 J=1,14
        T(J)=0.
   10 CONTINUE
      IF(DAY.NE.DAYL.OR.P32.NE.P(32)) CD32=COS(DR*(DAY-P(32)))
      IF(DAY.NE.DAYL.OR.P18.NE.P(18)) CD18=COS(2.*DR*(DAY-P(18)))
      IF(DAY.NE.DAYL.OR.P14.NE.P(14)) CD14=COS(DR*(DAY-P(14)))
      IF(DAY.NE.DAYL.OR.P39.NE.P(39)) CD39=COS(2.*DR*(DAY-P(39)))
      DAYL=DAY
      P32=P(32)
      P18=P(18)
      P14=P(14)
      P39=P(39)
C
C       F10.7
      T(1)=P(22)*DFA
C       TIME INDEPENDENT
      T(2)=P(2)*PLG(3,1)+P(3)*PLG(5,1)+P(23)*PLG(7,1)
     $     +P(27)*PLG(2,1)+P(15)*PLG(4,1)+P(60)*PLG(6,1)
C       SYMMETRICAL ANNUAL
      T(3)=(P(19)+P(48)*PLG(3,1)+P(30)*PLG(5,1))*CD32
C       SYMMETRICAL SEMIANNUAL
      T(4)=(P(16)+P(17)*PLG(3,1)+P(31)*PLG(5,1))*CD18
C       ASYMMETRICAL ANNUAL
      T(5)=(P(10)*PLG(2,1)+P(11)*PLG(4,1)+P(21)*PLG(6,1))*CD14
C       ASYMMETRICAL SEMIANNUAL
      T(6)=(P(38)*PLG(2,1))*CD39
C        DIURNAL
      IF(SW(7).EQ.0) GOTO 200
      T71 = P(12)*PLG(3,2)*CD14*SWC(5)
      T72 = P(13)*PLG(3,2)*CD14*SWC(5)
      T(7) =
     1 ((P(4)*PLG(2,2) + P(5)*PLG(4,2)
     2 + T71)*CTLOC
     4 + (P(7)*PLG(2,2) + P(8)*PLG(4,2)
     5 + T72)*STLOC)
  200 CONTINUE
C        SEMIDIURNAL
      IF(SW(8).EQ.0) GOTO 210
      T81 = (P(24)*PLG(4,3)+P(36)*PLG(6,3))*CD14*SWC(5)
      T82 = (P(34)*PLG(4,3)+P(37)*PLG(6,3))*CD14*SWC(5)
      T(8) =
     1 ((P(6)*PLG(3,3) + P(42)*PLG(5,3) + T81)*C2TLOC
     3 +(P(9)*PLG(3,3) + P(43)*PLG(5,3) + T82)*S2TLOC)
  210 CONTINUE
C        TERDIURNAL
      IF(SW(14).EQ.0) GOTO 220
      T(14) = P(40)*PLG(4,4)*S3TLOC
     $ +P(41)*PLG(4,4)*C3TLOC
  220 CONTINUE
C       MAGNETIC ACTIVITY
      IF(SW(9).EQ.0) GOTO 40
      IF(SW(9).EQ.1)
     $ T(9)=APDF*(P(33)+P(46)*PLG(3,1)*SWC(2))
      IF(SW(9).EQ.-1)
     $ T(9)=(P(51)*APT(1)+P(97)*PLG(3,1)*APT(1)*SWC(2))
   40 CONTINUE
      IF(SW(10).EQ.0.OR.SW(11).EQ.0.OR.LONG.LE.-1000.) GO TO 49
C        LONGITUDINAL
      T(11)= (1.+PLG(2,1)*(P(81)*SWC(5)*COS(DR*(DAY-P(82)))
     $           +P(86)*SWC(6)*COS(2.*DR*(DAY-P(87))))
     $        +P(84)*SWC(3)*COS(DR*(DAY-P(85)))
     $           +P(88)*SWC(4)*COS(2.*DR*(DAY-P(89))))
     $ *((P(65)*PLG(3,2)+P(66)*PLG(5,2)+P(67)*PLG(7,2)
     $   +P(75)*PLG(2,2)+P(76)*PLG(4,2)+P(77)*PLG(6,2)
     $    )*COS(DGTR*LONG)
     $  +(P(91)*PLG(3,2)+P(92)*PLG(5,2)+P(93)*PLG(7,2)
     $   +P(78)*PLG(2,2)+P(79)*PLG(4,2)+P(80)*PLG(6,2)
     $    )*SIN(DGTR*LONG))
   49 CONTINUE
      TT=0.
      DO 50 I=1,14
   50 TT=TT+ABS(SW(I))*T(I)
      GLOB7S=TT
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSU(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,
     $  MN1,ZN1,TN1,TGN1)
C       Calculate Temperature and Density Profiles for MSIS models
C       New lower thermo polynomial 10/30/89
      DIMENSION ZN1(MN1),TN1(MN1),TGN1(2),XS(5),YS(5),Y2OUT(5)
      COMMON/PARMB/GSURF,RE
      COMMON/LSQV/MP,II,JG,LT,QPB(50),IERR,IFUN,N,J,DV(60)
      SAVE
      DATA RGAS/831.4/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
CCCCCCWRITE(6,*) 'DB',ALT,DLB,TINF,TLB,XM,ALPHA,ZLB,S2,MN1,ZN1,TN1
      DENSU=1.
C        Joining altitude of Bates and spline
      ZA=ZN1(1)
      Z=AMAX1(ALT,ZA)
C      Geopotential altitude difference from ZLB
      ZG2=ZETA(Z,ZLB)
C      Bates temperature
      TT=TINF-(TINF-TLB)*EXP(-S2*ZG2)
      TA=TT
      TZ=TT
      DENSU=TZ
      IF(ALT.GE.ZA) GO TO 10
C
C       CALCULATE TEMPERATURE BELOW ZA
C      Temperature gradient at ZA from Bates profile
      DTA=(TINF-TA)*S2*((RE+ZLB)/(RE+ZA))**2
      TGN1(1)=DTA
      TN1(1)=TA
      Z=AMAX1(ALT,ZN1(MN1))
      MN=MN1
      Z1=ZN1(1)
      Z2=ZN1(MN)
      T1=TN1(1)
      T2=TN1(MN)
C      Geopotental difference from Z1
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 20 K=1,MN
        XS(K)=ZETA(ZN1(K),Z1)/ZGDIF
        YS(K)=1./TN1(K)
   20 CONTINUE
C        End node derivatives
      YD1=-TGN1(1)/(T1*T1)*ZGDIF
      YD2=-TGN1(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1./Y
      DENSU=TZ
   10 IF(XM.EQ.0.) GO TO 50
C
C      CALCULATE DENSITY ABOVE ZA
      GLB=GSURF/(1.+ZLB/RE)**2
      GAMMA=XM*GLB/(S2*RGAS*TINF)
      EXPL=EXP(-S2*GAMMA*ZG2)
      IF(EXPL.GT.50.OR.TT.LE.0.) THEN
        EXPL=50.
      ENDIF
C       Density at altitude
      DENSA=DLB*(TLB/TT)**(1.+ALPHA+GAMMA)*EXPL
      DENSU=DENSA
      IF(ALT.GE.ZA) GO TO 50
C
C      CALCULATE DENSITY BELOW ZA
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       integrate spline temperatures
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50..OR.TZ.LE.0.) THEN
        EXPL=50.
      ENDIF
C       Density at altitude
      DENSU=DENSU*(T1/TZ)**(1.+ALPHA)*EXP(-EXPL)
   50 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
      FUNCTION DENSM(ALT,D0,XM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,TGN2)
C       Calculate Temperature and Density Profiles for lower atmos.
      DIMENSION ZN3(MN3),TN3(MN3),TGN3(2),XS(10),YS(10),Y2OUT(10)
      DIMENSION ZN2(MN2),TN2(MN2),TGN2(2)
      COMMON/PARMB/GSURF,RE
      COMMON/FIT/TAF
      COMMON/LSQV/MP,II,JG,LT,QPB(50),IERR,IFUN,N,J,DV(60)
      SAVE
      DATA RGAS/831.4/
      ZETA(ZZ,ZL)=(ZZ-ZL)*(RE+ZL)/(RE+ZZ)
      DENSM=D0
      IF(ALT.GT.ZN2(1)) GOTO 50
C      STRATOSPHERE/MESOSPHERE TEMPERATURE
      Z=AMAX1(ALT,ZN2(MN2))
      MN=MN2
      Z1=ZN2(1)
      Z2=ZN2(MN)
      T1=TN2(1)
      T2=TN2(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 210 K=1,MN
        XS(K)=ZETA(ZN2(K),Z1)/ZGDIF
        YS(K)=1./TN2(K)
  210 CONTINUE
      YD1=-TGN2(1)/(T1*T1)*ZGDIF
      YD2=-TGN2(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       Temperature at altitude
      TZ=1./Y
      IF(XM.EQ.0.) GO TO 20
C
C      CALCULATE STRATOSPHERE/MESOSPHERE DENSITY
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C       Integrate temperature profile
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.) EXPL=50.
C       Density at altitude
      DENSM=DENSM*(T1/TZ)*EXP(-EXPL)
   20 CONTINUE
      IF(ALT.GT.ZN3(1)) GOTO 50
C
C      TROPOSPHERE/STRATOSPHERE TEMPERATURE
      Z=ALT
      MN=MN3
      Z1=ZN3(1)
      Z2=ZN3(MN)
      T1=TN3(1)
      T2=TN3(MN)
      ZG=ZETA(Z,Z1)
      ZGDIF=ZETA(Z2,Z1)
C       Set up spline nodes
      DO 220 K=1,MN
        XS(K)=ZETA(ZN3(K),Z1)/ZGDIF
        YS(K)=1./TN3(K)
  220 CONTINUE
      YD1=-TGN3(1)/(T1*T1)*ZGDIF
      YD2=-TGN3(2)/(T2*T2)*ZGDIF*((RE+Z2)/(RE+Z1))**2
C       Calculate spline coefficients
      CALL SPLINE(XS,YS,MN,YD1,YD2,Y2OUT)
      X=ZG/ZGDIF
      CALL SPLINT(XS,YS,Y2OUT,MN,X,Y)
C       temperature at altitude
      TZ=1./Y
      IF(XM.EQ.0.) GO TO 30
C
C      CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY
C
      GLB=GSURF/(1.+Z1/RE)**2
      GAMM=XM*GLB*ZGDIF/RGAS
C        Integrate temperature profile
      CALL SPLINI(XS,YS,Y2OUT,MN,X,YI)
      EXPL=GAMM*YI
      IF(EXPL.GT.50.) EXPL=50.
C        Density at altitude
      DENSM=DENSM*(T1/TZ)*EXP(-EXPL)
   30 CONTINUE
   50 CONTINUE
      IF(XM.EQ.0) DENSM=TZ
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
C        CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
C        ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
C        X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        N: SIZE OF ARRAYS X,Y
C        YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
C                 >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
C        Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      SAVE
      IF(YP1.GT..99E30) THEN
        Y2(1)=0
        U(1)=0
      ELSE
        Y2(1)=-.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     $    /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
   11 CONTINUE
      IF(YPN.GT..99E30) THEN
        QN=0
        UN=0
      ELSE
        QN=.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
   12 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
C        CALCULATE CUBIC SPLINE INTERP VALUE
C        ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA FOR INTERPOLATION
C        Y: OUTPUT VALUE
      DIMENSION XA(N),YA(N),Y2A(N)
      SAVE
      KLO=1
      KHI=N
    1 CONTINUE
      IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF(H.EQ.0) write(stderr,*) 'BAD XA INPUT TO SPLINT'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     $  ((A*A*A-A)*Y2A(KLO)+(B*B*B-B)*Y2A(KHI))*H*H/6.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPLINI(XA,YA,Y2A,N,X,YI)
C       INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
C        XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
C        Y2A: ARRAY OF SECOND DERIVATIVES
C        N: SIZE OF ARRAYS XA,YA,Y2A
C        X: ABSCISSA ENDPOINT FOR INTEGRATION
C        Y: OUTPUT VALUE
      DIMENSION XA(N),YA(N),Y2A(N)
      SAVE
      YI=0
      KLO=1
      KHI=2
    1 CONTINUE
      IF(X.GT.XA(KLO).AND.KHI.LE.N) THEN
        XX=X
        IF(KHI.LT.N) XX=AMIN1(X,XA(KHI))
        H=XA(KHI)-XA(KLO)
        A=(XA(KHI)-XX)/H
        B=(XX-XA(KLO))/H
        A2=A*A
        B2=B*B
        YI=YI+((1.-A2)*YA(KLO)/2.+B2*YA(KHI)/2.+
     $     ((-(1.+A2*A2)/4.+A2/2.)*Y2A(KLO)+
     $     (B2*B2/4.-B2/2.)*Y2A(KHI))*H*H/6.)*H
        KLO=KLO+1
        KHI=KHI+1
        GOTO 1
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION DNET(DD,DM,ZHM,XMM,XM)
C       TURBOPAUSE CORRECTION FOR MSIS MODELS
C         Root mean density
C       8/20/80
C          DD - diffusive density
C          DM - full mixed density
C          ZHM - transition scale length
C          XMM - full mixed molecular weight
C          XM  - species molecular weight
C          DNET - combined density
      SAVE
      A=ZHM/(XMM-XM)
      IF(DM.GT.0.AND.DD.GT.0) GOTO 5
        WRITE(stderr,*) 'DNET LOG ERROR',DM,DD,XM
        error stop
        IF(DD.EQ.0.AND.DM.EQ.0) DD=1.
        IF(DM.EQ.0) GOTO 10
        IF(DD.EQ.0) GOTO 20
    5 CONTINUE
      YLOG=A*ALOG(DM/DD)
      IF(YLOG.LT.-10.) GO TO 10
      IF(YLOG.GT.10.)  GO TO 20
        DNET=DD*(1.+EXP(YLOG))**(1/A)
        GO TO 50
   10 CONTINUE
        DNET=DD
        GO TO 50
   20 CONTINUE
        DNET=DM
        GO TO 50
   50 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR(ALT, R,H1,ZH)
C        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
C        ALT - altitude
C        R - target ratio
C        H1 - transition scale length
C        ZH - altitude of 1/2 R
      SAVE
      E=(ALT-ZH)/H1
      IF(E.GT.70.) GO TO 20
      IF(E.LT.-70.) GO TO 10
        EX=EXP(E)
        CCOR=R/(1.+EX)
        GO TO 50
   10   CCOR=R
        GO TO 50
   20   CCOR=0.
        GO TO 50
   50 CONTINUE
      CCOR=EXP(CCOR)
       RETURN
      END
C-----------------------------------------------------------------------
      FUNCTION  CCOR2(ALT, R,H1,ZH,H2)
C       O&O2 CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
      E1=(ALT-ZH)/H1
      E2=(ALT-ZH)/H2
      IF(E1.GT.70. .OR. E2.GT.70.) GO TO 20
      IF(E1.LT.-70. .AND. E2.LT.-70) GO TO 10
        EX1=EXP(E1)
        EX2=EXP(E2)
        CCOR2=R/(1.+.5*(EX1+EX2))
        GO TO 50
   10   CCOR2=R
        GO TO 50
   20   CCOR2=0.
        GO TO 50
   50 CONTINUE
      CCOR2=EXP(CCOR2)
       RETURN
      END
C-----------------------------------------------------------------------
      end module msise00_gemini
