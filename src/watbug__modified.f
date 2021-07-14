      program watbug

      real T(481),P(481),PE(481),APE(481),D(481),
     1          DAYS(13),H(481),AE(481),ST(481),DST(481),DEF(481),
     2          SUR(481),DY(481)

      real MON(481)

      character*72 LABEL

      data DAYS/0.0,31.0,28.0,31.0,30.0,31.0,30.0,31.0,31.0,30.0,31.0,
     1          30.0,31.0/

      real LAT,INDEX
      integer TUNIT, PUNIT

c
c     ********************************************************************
c     *                                                                  *
c     *  This algorithm was developed by C. Willmott at the Department   *
c     *  of Geography, University of Delaware in 1978 in order to        *
c     *  facilitate the calculation of climatic water budgets.  A        *
c     *  minimum amount of data (i.e. air temperature, precipitation     *
c     *  and a few initial parameters) and no "look up" tables are       *
c     *  required as all relationships are explicitly specified.  The    *
c     *  program was refined on a Burroughs' B7700 although standard     *
c     *  (ANSI compatible) Fortran was used. It should, therefore, run   *
c     *  with few or no modifications on most moderate to large sized    *
c     *  machines.  If problems are encountered, however, users are      *
c     *  urged to contact the author.                                    *
c     *                                                                  *
c     ********************************************************************
c
c     Initial parameters:
c
c     'LABEL'  72 character alphanumeric problem title.
c
c     'N'      Total number of months or days over which soil moisture
c              balancing is to occur.  If N equals '1', balancing
c              does not occur ans ST(1) must be specified.
c
c     'NT'     Total number of months or days over which the water budget 
c              is to be calculated.
c
c     'KD'     The day of the month where the first calculations
c              are to begin. KD must be less than or equal
c              to the number of days in month KM.
c
c     'KM'     The first month of calculations.
c              KM must be between 0 and 12.
c
c     'KY'     The first year of calculations. Last two digits only(!).
c
c     'FC'     Soil water holding or field capacity ot the top (only)
c              soil layer in millimeters.
c
c     'SM'     Determines a resistance function of soil water to removal
c              by evapotranspiration.  Blank or 0 indicates that the
c              availability of soil moisture to evapotranspiration will
c              decline linearly with the ratio of actual to potential
c              maximum soil moisture. Any other numeric designation
c              will result in an alternative procedure where moisture is
c              withdrawn at the maximum rate until the actual/potential
c              ratio drops below 0.7 at which time a linear decline in 
c              availability is assumed (see Mather, 1974: 106 - curves
c              C and G).
c
c     'LAT'    The latitude in degrees.
c
c     'DT'     Time differential. Blank or 0 indicates monthly
c              calculations.  Any other number causes daily calculations.
c
c     'TUNIT'  Designates the units of air temperature. 1 means
c              the raw temperature data are in degrees Fahrenheit.
c              2 means degrees Kelvin. Any other numeric
c              designation or a blank means degrees Celsius.
c
c     'PUNIT'  Units of precipitation. 1 means the raw data
c              are in centimeters. 2 means inches. 3 means hundreths
c              of an inch. Other designations or blanks means millimeters.
c
c     'ST'     Estimated soil moisture content of the top soil
c              layer just prior to beginning of calculations.
c              ST only needs to be specified when balancing 
c              is not to be done (see note below). ST(1) is in mm.
c
c     'HEAT'   Estimated heat index. It needs to be specified
c              only when soil moisture balancing does not occur.
c              (Note: balancing should only be specified for
c              periods containing one or more complete years of data).
c     'INDEX'  Should be set greater than zero when calculations 
c              for a subsequent station are to follow these.  (Note:
c              control parameters and data must be included
c              sequentially in the input data set for each station
c              that is to be evaluated).
c
c     Read initial parameters:
c

   10 read(5,1000,end=20,err=20) LABEL,N,NT,KD,KM,KY,FC,SM,LAT,DT,
     1                           TUNIT,PUNIT,ST(1),HEAT,INDEX

      NNN = NNN + 1

c
c     Assumed parameter values (i.e. when they are not specified).
c

      if (N .eq. 0)  N = 1
      if (NT .eq. 0)  NT = N
      if ((KD .eq. 0) .and. (DT .ne. 0.0))  KD = 1
      if ((KD .eq. 0) .and. (DT .eq. 0.0))  KD = 15
      if (KM .eq. 0)  KM = 1

c
c     Set array sizes for calculating a soil water balance
c

      M = N + 1

c
c     Call the main subprogram which controls all calculations
c 

      call main(N,NT,M,FC,LAT,KD,KM,KY,DT,DY,HEAT,SM,T,P,PE,APE,D,
     1          DAYS,H,AE,ST,DST,DEF,SUR,MON,LABEL,TUNIT,PUNIT)

c
c     Test to see if subsequent stations are to be evaluated.
c

      if (INDEX .GT. 0.0)  go to 10

   20 continue
      stop 
 1000 format(A72,/,5I5,4F5.0,2I5,3F5.0)

      end program watbug

c*************************************************************************

      subroutine main(N,NT,M,FC,LAT,KD,KM,KY,DT,DY,HEAT,SM,T,P,PE,APE,D,
     1                DAYS,H,AE,ST,DST,DEF,SUR,MON,LABEL,TUNIT,PUNIT)
      
      real LAT,T(M),P(M),PE(M),APE(M),D(M),DAYS(13),H(M),DY(M),AE(M),
     1     ST(M),DST(M),DEF(M),SUR(M),MON(M)

      character*8 FMT
      character*72 LABEL
      dimension OUT(10),SUMM(5),SUMY(5)
      integer TUNIT, PUNIT
      integer IND(5) /3,4,8,9,10/

      H = -9999.


c
c     Read the data format (FMT).
c

      read(5,1000)  FMT
      NNN = 0

c 
c     Unit correction factors
c

      C1 = 1.0
      C2 = 1.0
      FK = 0.0
      if (TUNIT .eq. 1)  FK = 32.0
      if (TUNIT .eq. 1)  C1 = 5.0/9.0
      if (TUNIT .eq. 2)  FK = 273.16
      if (PUNIT .eq. 1)  C2 = 10.0
      if (PUNIT .eq. 2)  C2 = 25.4
      if (PUNIT .eq. 3)  C2 = 0.254

c 
c     Read air temperature and precipitation data
c
      T = -999.
      P = -999.

      do 10 I=1,N
        read(5,FMT,end=290,err=280)  T(I),P(I)
        NNN = NNN + 1
c
c     Unit translations
c
        T(I) = C1 * (T(I) - FK)
        P(I) = C2 * P(I)

   10 continue

c
c     Original code has a Y2K problem 
c
      KY = KY + 1900

c
c     Test for daily, monthly, day by day or month by month budgeting
c

      if ((DT .ne. 0.0) .and. (N .eq. 1))  go to 70
      if ((DT .ne. 0.0) .and. (N .eq. 1)) go to 20
      if (DT .ne. 0.0)  go to 60

c***
c***  Here for monthly balancing
c***

      call date(N,M,KD,KM,DY,MON,DT,DAYS)
      call mather(N,M,H,T,HEAT,A,PE,APE,DAYS,LAT,DL,KD,KM,DT,MON,DY)
      call diff(N,M,P,APE,D,DEF)
      call bal(N,M,ST,D,FC,SM,SUR,DST,DT,KM)
      call evapo(N,M,D,AE,APE,P,DST,DEF)

c
c     Write monthly input data and results
c
   20 continue

      write(6,1010) LABEL
      write(6,1020) N,NT,FC,LAT

      if (N .eq. 1) go to 130   

      I = 0
      go to 40
   30 continue
      call conv(SUMY,5,1,5)
      write(6,1030)  SUMY
   40 call init(SUMY,5)
      write(6,1040)  KY
      
      KY = KY + 1
      write(6,1050)
   50 I = I + 1

c
c     Round off to nearest whole number and get totals before writing
c

      call output(PE,APE,P,D,ST,DST,AE,DEF,SUR,M,OUT,I)
      call toty(OUT,10,IND,SUMY,5)
      call conv(OUT,10,2,10)
      OUT(1)=T(I)

      write(6,1060)  int(MON(I)),OUT
      if ((I .lt. N) .and. (MON(I) .eq. 12)) go to 30
      if ((I .eq. N) .and. (NT .gt. N)) go to 130
      if (I .eq. N) go to 290
      go to 50

   60 continue

c***
c*** Here for daily balancing
c***
      call date(N,M,KD,KM,DY,MON,DT,DAYS)
      call mather(N,M,H,T,HEAT,A,PE,APE,DAYS,LAT,DL,KD,KM,DT,MON,DY)
      call diff(N,M,P,APE,D,DEF)
      call bal(N,M,ST,D,FC,SM,SUR,DST,DT,KM)
      call evapo(N,M,D,AE,APE,P,DST,DEF)

c
c     Write daily input data and results.
c
   70 continue

      write(6,1010) LABEL
      write(6,1070) N,NT,FC,LAT

      if (N .eq. 1) go to 190

      I = 0
      go to 90
   80 continue
      call conv(SUMM,5,1,5)
      call conv(SUMY,5,1,5)
      
      write(6,1080) SUMM
      write(6,1030) SUMY

   90 call init(SUMY,5)
      write(6,1040)  KY
      KY = KY + 1
      go to 110
  100 continue
      call conv(SUMM,5,1,5)
      write(6,1080)  SUMM
  110 call init(SUMM,5)
      write(6,1090)  MON(I+1)
      write(6,1100)
  120 I = I + 1
      KD = DY(I)
      KM = MON(I)

c
c     Round off to nearest whole numbers and get totals before writing.
c

      call output(PE,APE,P,D,ST,DST,AE,DEF,SUR,M,OUT,I)
      call totm(OUT,10,IND,SUMM,5)
      call toty(OUT,10,IND,SUMY,5)
      call conv(OUT,10,2,10)

      OUT(1) = T(I)

      write(6,1060) KD,OUT
      if ((I .lt. N) .and. (MON(I) .eq.12) .and. (DY(I) .eq. DAYS(KM + 1))) go to 80
      if ((I .lt. N) .and. (DY(I) .eq. DAYS(KM + 1))) go to 100
      if ((I .eq. N) .and. (NT .le. N)) go to 290
      if (I .eq. N) go to 190
      go to 120

c***
c***  Here for month-by-month calculations.
c***

  130 continue
      N1 = 0
      if (N .gt. 1)  N1 = 1

c
c     Get the initial soil moisture, month and day.
c

      if (N .gt. 1)  ST(1) = ST(N)
      KD = 15
      if (N .gt. 1)  KM = MON(N) + 1
      if (KM .ge. 13)  KM = 1
      DY(1) = KD
      MON(1) = KM

c
c     Set initial parameters.
c 

      NN = NT - N
      if (N .eq. 1)  NN = NT
      if (NNN .eq. 1)  NNN = 0
      N = 1
      LL = 0
      M = N + 1

c
c     Test for appropriate labels.
c

      if (NNN .eq. 0)  go to 150
      if ((KM .gt. 1) .and. (N1 .eq. 1)) go to 170
      if (KM .gt. 1)  go to 160

  140 continue

c
c     Write labels, year, and last years totals.
c

      call conv(SUMY,5,1,5)
      write(6,1030)  SUMY
  150 call init(SUMY,5)
      write(6,1040)  KY
      KY = KY + 1
  160 continue
      write(6,1050)

c
c     Read input data and call budget subroutines.
c

  170 LL = LL + 1
      if (NNN .eq. 0)  go to 180
      read(5,FMT,end=290,err=280)  T(N),P(N)
      T(N) = C1 * (T(N) - FK)
      P(N) = C2 * P(N)
  180 NNN = NNN + 1

      call mather(N,M,H,T,HEAT,A,PE,APE,DAYS,LAT,DL,KD,KM,DT,MON,DY)
      call diff(N,M,P,APE,D,DEF)
      D(N + 1) = D(N)
      call bal(N,M,ST,D,FC,SM,SUR,DST,DT,KM)
      ST(N) = ST(N + 1)
      SUR(N) = SUR(N + 1)
      DST(N) = DST(N + 1)
      call evapo(N,M,D,AE,APE,P,DST,DEF)

c
c     Round off to nearest whole number and get totals before writing.
c

      call output(PE,APE,P,D,ST,DST,AE,DEF,SUR,M,OUT,N)
      call toty(OUT,10,IND,SUMY,5)
      call conv(OUT,10,2,10)
      OUT(1) = T(N)

c
c     Write results and get next month.
c

      write(6,1060)  KM,OUT
      call date(N,M,KD,KM,DY,MON,DT,DAYS)

      if (NN .eq. LL)  go to 290
      if (KM .eq. 1)  go to 140
      go to 170

c***
c***  Here for day-by-day calculations.
c***

  190 continue
      N1 = 0
      if (N .gt. 1)  N1 = 1

c
c     Get the initial soil moisture, month, and day.
c

      if (N .gt. 1)  ST(1) = ST(N)
      if (N .gt. 1)  KD = DY(N) + 1
      if (N .gt. 1)  KM = MON(N)
      if (KD .gt. DAYS(KM + 1))  go to 200
      DY(1) = KD
      MON(1) = KM
      go to 210
  200 continue
      KM = KM + 1
      KD = 1
      if (KM .ge. 13)  KM = 1
      DY(1) = KD
      MON(1) = KM
  
  210 continue

c
c     Initialize parameters.
c

      NN = NT - N
      if (N .eq. 1)  NN = NT
      if (NNN .eq. 1)   NNN = 0
      N = 1
      L = 0
      M = N + 1

c
c     Test for appropriate labels.
c
      if (NNN .eq. 0) go to 230
      if ((KD .ne. 1) .and. (KM .ne. 1) .and. (N1 .eq.1))  go to 260
      if ((KD .ne. 1) .or. (KM .ne. 1)) go to 240

  220 continue

c
c     Write labels, the year, month, and last years or months totals.
c

      call conv(SUMM,5,1,5)
      call conv(SUMY,5,1,5)
      write(6,1080) SUMM
      write(6,1030)  SUMY
  230 call init(SUMY,5)
      write(6,1040)  KY
      KY = KY + 1
      go to 250
  240 continue
      call conv(SUMM,5,1,5)
      write(6,1080)  SUMM
  250 call init(SUMM,5)
      write(6,1090)  KM
      write(6,1100)
      
c
c     Read input data and call budget subroutines.
c

  260 L = L + 1

      if (NNN .eq. 0)  go to 270
      read(5,FMT,end=290,err=280)  T(N),P(N)
      T(N) = C1 * (T(N) - FK)
      P(N) = C2 * P(N)
  270 NNN = NNN + 1

      call mather(N,M,H,T,HEAT,A,PE,APE,DAYS,LAT,DL,KD,KM,DT,MON,DY)
      call diff(N,M,P,APE,D,DEF)
      D(N + 1) = D(N)
      call bal(N,M,ST,D,FC,SM,SUR,DST,DT,KM)
      ST(N) = ST(N+1)
      SUR(N) = SUR(N+1)
      DST(N) = DST(N+1)
      call evapo(N,M,D,AE,APE,P,DST,DEF)

c
c     Round off to nearest whole number and get totals before writing.
c

      call output(PE,APE,P,D,ST,DST,AE,DEF,SUR,M,OUT,N)
      call totm(OUT,10,IND,SUMM,5)
      call toty(OUT,10,IND,SUMY,5)
      call conv(OUT,10,2,10)
      OUT(1) = T(N)

c
c     Write results and get a new date.
c

      write(6,1060)  KD,OUT
      call date(N,M,KD,KM,DY,MON,DT,DAYS)

      if (NN .eq. L)  go to 290
      if ((KM .eq. 1) .and. (KD .eq. 1)) go to 220
      if (KD .eq. 1) go to 240
      go to 260

c
c     Write final messages and totals.
c

  280 continue
      write(6,1110)  NNN + 1
      go to 300

  290 continue
      if (DT .ne. 0.0)  write(6,1080)  SUMM
      write(6,1030)  SUMY
      write(6,1120)  NNN
      
  300 continue
      write(6,1130)

      return

 1000 format (A)
 1010 format (//,A)
 1020 format (//,'  No. of months over which balancing occurs is ',I5,//,
     1 '  Total No. of months evaluated is ',I5,//,
     2 '  Soil moisture capacity is ',F5.1,' mm ',//,
     3 '  Latitude is ',F4.1)
     
 1030 format (/,'  Yearly totals',5X,2F7.1,21X,3F7.1)
 1040 format (///,'  Year is ', I4)
 1050 format (//,'  MO     TEMP    UPE    APE   PREC   DIFF     ST    DST     AE',
     1        '    DEF   SURP',/)
 1060 format (I4,2X,F7.1,9F7.1)
 1070 format (//,'   No. of days over which balancing occurs is ',I5,//,
     1        '   Total No. of days evaluated is ',I5,//,
     2        '   Soil moisture capacity is ',F5.1,' mm ',//,
     3        '   Latitude is ',F4.1)
 1080 format (/,'  Monthly totals',4X,2F7.1,21X,3F7.1)
 1090 format (//,'   Month is ',I2)
 1100 format (//,'  DY     TEMP    UPE    APE   PREC   DIFF     ST    DST     AE',
     1        '    DEF   SURP',/)
 1110 format (///,'   Error encountered in the data at record ',I5)
 1120 format (///,'   Processing terminated after record ',I5)
 1130 format ('1')
      end subroutine main
      
c*************************************************************************      
      subroutine date(N,M,KD,KM,DY,MON,DT,DAYS)

        real DY(M),MON(M),DAYS(12)

c
c       Generate day and month designations.
c
c       Test for monthly, daily, month-by-month, or day-by-day calculations.
c

        if ((DT .ne. 0.0) .and. (N .eq. 1))  go to 60
        if ((DT .eq. 0.0) .and. (N .eq. 1))  go to 50
        if (DT .ne. 0.0)  go to 20

c
c       Monthly calculations.
c

        KD = 15
        KM = KM - 1
        if (KM .le. 0)  KM = 0
        do 10 I=1,N
            KM = KM + 1
            if (KM .ge. 13)  KM = 1
            MON(I) = KM
            DY(I) = KD
  10    continue
        go to 80

c
c       Daily calculations.
c

  20    continue
        K = 1
        KD = KD - 1
        if (KD .le. 0)  KD = 0
        if (KD .gt. 0)  K = KD
        KM = KM - 1
        if (KM .le. 0)  KM = 0
        J = 0
  30    KM = KM + 1
        if (KD .gt. DAYS(KM+1))  KD = 1
        if (KM .ge. 13)  KM = 1
        do 40 I=K,DAYS(KM+1)
            J = J + 1
            MON(J) = KM
            DY(J) = I
  40    continue
        if (J .ge. N)  go to 80
        K = 1
        go to 30

c
c       Month-by-month calculations.
c

  50    continue

        KD = 15
        KM = KM + 1
        if (KM .ge. 13)  KM = 1
        DY(N) = KD
        MON(N) = KM
        go to 80

c
c       Day-by-day calculations.
c

  60    continue
        KD = KD + 1
        DY(N) = KD
        MON(N) = KM
        if (KD .gt. DAYS(KM+1))  go to 70
        go to 80

  70    continue
        KM = KM + 1
        KD = 1
        if (KM .ge. 13)  KM = 1
        DY(N) = KD
        MON(N) = KM

  80   continue
       return

      end subroutine date

c*************************************************************************
      subroutine mather(N,M,H,T,HEAT,A,PE,APE,DAYS,LAT,DL,KD,KM,DT,MON,DY)

        real LAT,H(M),T(M),PE(M),APE(M),DAYS(12),MON(M),DY(M)

c
c       Calculate potential evapotranspiration.
c
c       When 'LAT' is greater than 50 degrees, the daylength correction
c       remains equal to that for 50 degs. 'ALAT' is, therefore,
c       used as the argument for subrouting 'day'.
c

        ALAT = LAT
        if (ALAT .ge. 50.0)  ALAT = 50.0

c
c       Calculate the heat index during balancing on the first call
c       of 'mather'. On the second call, go directly to 'PE' calculations.
c

        if (N .lt. 12)  go to 40
        if ((N .lt. 365) .and. (DT .ne. 0.0)) go to 40

        XN = N
        HEAT = 0.0
        do 30 I=1,N
            if (T(I) .le. 0.0)  go to 10
            H(I) = (T(I) / 5.0) ** 1.514
            go to 20
  10        H(I) = 0.0
  20        continue
            HEAT = HEAT + H(I)
  30    continue

c
c       Adjust 'HEAT' for budgets greater than a year.
c

        HEAT = HEAT * 12.0 / XN

c
c       Note: 'A' is an empirically derived exponent based upon 'HEAT'
c

  40    continue
        A = 6.75 / 10.0**7.0 * HEAT**3.0
     1      - 7.71 / 10.0**5.0 * HEAT**2.0 
     2      + 1.79 / 10.0**2.0 * HEAT
     3      + 0.49

c
c       Get initial monthly PE, I.E. based upon 30 days in
c       a month and 12 hours in a day
c       Note: PE(I) and APE(I) are calculated in mm / month or day.
c

        do 70 I=1,N
            if (T(I) .le. 0.0)  go to 50
            PE(I) = 16.0 * (10.0 * T(I) / HEAT)**A

c
c           Correct for temperatures greater than 26.5 deg C.
c           See Thornthwaite (1948) for explanation.
c

            if (T(I) .ge. 26.5)  PE(I) = (-41.58547 + 3.22441 * T(I)
     1                                   -0.04325 * T(I)**2.0) * 10.0

c
c           Correct PE(I) for daily calculations.
c

            if (DT .ne. 0.0)  PE(I) = PE(I) / 30.0
            go to 60
   50       PE(I) = 0.0
   60       continue
            KD = DY(I)
            KM = MON(I)

c
c           Adjust PE for daylength and the number of days in a month.
c

            call day(DAYS,ALAT,KD,KM,DT,DECD,DL)
            APE(I) = PE(I) * (DAYS(KM+1) / 30.0) * (DL / 12.0)
            
  70        continue
            return
      end subroutine mather

c*************************************************************************
      subroutine day(DAYS,LAT,KD,KM,DT,DECD,DL)

        real DAYS(12)
        real LAT

c
c       Calculate the number of hours in a 
c       day and the solar declination associated with that
c       day.  The input required includes: the month (KM), 
c       the dat (KD) and the latitude (LAT).
c

        X = 0.0
        do 10  I=1,KM
            X = X + DAYS(I)
   10   continue
        SUM = X + KD

c
c       Get the number of days since the vernal equinox (March 21).
c

        DAYL = SUM - 80.0
        if (DAYL .le. 0.0)  DAYL = 285.0 + SUM

c
c       Calculate the declination.
c

        DECD = 23.45 * sin(DAYL / 365.0 * 6.2832)
        DECR = DECD * 0.017453

c
c       Calculate number of hours of daylight corresponding to
c       day (KD) and month (KM) (see Sellers, 1965).
c

        CZ = cos(1.5708 + 0.01745 * (100.0 / 60.0))
        ALAT = LAT * 0.017453
        XX = cos(DECR) * cos(ALAT)
        if (XX .le. 0.0)  go to 20

        CSH = (CZ - sin(DECR) * sin(ALAT)) / XX
        H = acos(CSH)
        DL = 24.0 * H / 3.1416
        go to 30

c
c       Error message - divide by zero or less.
c

   20   write(6,1000)

   30   continue
        return
 1000   format ('0', ' Error - divide by zero or less - lat. ',//,
     1          '   or the declination is probably incorrect ')
      
      end subroutine day
        
c*************************************************************************
      subroutine diff(N,M,P,APE,D,DEF)

        real P(M),APE(M),D(M),DEF(M)

c
c       Compare APE(I) with precipitation (P(I)).
c

        do 10 I=1,N
            D(I) = P(I) - APE(I)
            DEF(I) = 0.0
            if (D(I) .lt. 0.0)  DEF(I) = D(I)
  10    continue
        return
      
      end subroutine diff

c*************************************************************************
      subroutine bal(N,M,ST,D,FC,SM,SUR,DST,DT,KM)

        real ST(M),D(M),SUR(M),DST(M)

c
c       Iterate for soil moisture terms that balance
c       the water budget on the first call of 'bal'. On a second call
c       and/or when N is '1', 'bal' does month-by-month or day-by-day
c       soil moisture calculations.
c

        if (N .eq. 1)  go to 10

        ST(N+1) = 0.0
c***
c***    Shouldn't ST(1) = FC?
c***
        ST(1) = 300.0
        DST(1) = 0.0
        K = 0
        Z = 0.0
   10   continue

        do 80 I=2,N+1

            if (D(I) .ge. 1.e-5)  go to 50

c
c           Test for monthly or daily withdrawal. 
c

            if (DT .ne. 0.0)  go to 30

c
c           Withdrawal for monthly budgets (Note: this is done on
c           an approximate day-by-day basis).
c

            X1 = ST(I - 1)
            do 20 J=1,30
                RATIO = ST(I-1) / FC
                ST(I) = ST(I-1) + D(I) * RATIO / 30.0
                if ((RATIO .ge. 0.7) .and. (SM .gt. 0.0))  ST(I) = ST(I-1) + D(I) / 30.0
                ST(I-1) = ST(I)
   20       continue
            ST(I-1) = X1
            go to 40
   30       continue

c
c           Withdrawal for daily budgets.
c

            ST(I) = ST(I-1) + ST(I-1) / FC * D(I)
            if ((ST(I-1) / FC .ge. 0.7) .and. (SM .gt. 0.0))  ST(I) = ST(I-1) + D(I)

  40        continue
            if (ST(I) .le. 1.0)  ST(I) = 1.0
            go to 70

  50        continue

            ST(I) = ST(I-1) + D(I)
            if (ST(I) .ge. FC)  go to 60
            SUR(I) = 0.0
            go to 70

  60        continue
            SUR(I) = ST(I) - FC
            ST(I) = FC

  70        continue
            DST(I) = ST(I) - ST(I-1)
            if (D(I) .le. 0.0)  SUR(I) = 0.0

  80        continue

            if (N .eq. 1)  go to 160

            K = K + 1

c
c           Tests for balances.
c            

            if (K .gt. 50) go to 160
            XX = abs(ST(N+1) + DST(1) - ST(1))
            if ((XX .lt. 1.0) .and. (Z .eq. 1.0))  go to 160
            if (XX .lt. 1.0)  go to 90
            ST(1) = ST(N+1) + DST(1)
            go to 10

  90        continue
            if (D(1) .ge. 0.0)  go to 130
            if (DT .ne. 0.0)  go to 110

c
c           Balance for the first month
c

            X2 = ST(N+1)
            do 100 L=1,30
                RATIO = ST(N+1) / FC
                ST(1) = ST(N+1) + D(1) * RATIO / 30.0
                if ((RATIO .ge. 0.7) .and. (SM .gt. 0.0))  ST(1) = ST(N+1) + D(1) / 30.0
                ST(N+1) = ST(1)
  100       continue
            ST(N+1) = X2
            go to 120
  110       continue

c
c           Balance for the first day.
c

            ST(1) = ST(N+1) + ST(N+1) / FC * D(1)
            if ((ST(N+1) / FC .ge. 0.7) .and. (SM .gt. 0.0))  ST(1) = ST(N+1) + D(1)

  120       continue
            if (ST(1) .le. 1.0)  ST(1) = 1.0
            go to 150

  130       continue
            ST(1) = ST(N+1) + D(1)
            if (ST(1) .ge. FC)  go to 140
            SUR(1) = 0.0
            go to 150

  140       SUR(1) = ST(1) - FC
            ST(1) = FC

  150       continue
            DST(1) = ST(1) - ST(N+1)
            if (D(1) .le. 0.0)  SUR(1) = 0.0
            Z = 1.0
            go to 10

  160       continue

            return

      end subroutine bal

c*************************************************************************
      subroutine evapo(N,M,D,AE,APE,P,DST,DEF)

        real D(M),AE(M),APE(M),P(M),DST(M),DEF(M)

c
c       Calculate actual evapotranspiration and deficit.
c

        do 10 I=1,N
            AE(I) = APE(I)
            if (D(I) .lt. 0.0)  AE(I) = P(I) + abs(DST(I))
            DEF(I) = APE(I) - AE(I)

  10    continue
        
        return

      end subroutine evapo

c*************************************************************************
      subroutine init(SUM,N)

        dimension SUM(N)

c
c       Initialize array 'SUM' with zeroes.
c

        I = 0
   10   I = I + 1
        SUM(I) = 0.0
        if (I .lt. N)  go to 10

        return

      end subroutine init

c*************************************************************************
      subroutine output(PE,APE,P,D,ST,DST,AE,DEF,SUR,M,OUT,L)

        dimension PE(M),APE(M),P(M),D(M),ST(M),DST(M),AE(M),DEF(M),SUR(M),
     1            OUT(10)

c
c       Fill the output array 'OUT'.
c

        OUT(2) = PE(L)
        OUT(3) = APE(L)
        OUT(4) = P(L)
        OUT(5) = D(L)
        OUT(6) = ST(L)
        OUT(7) = DST(L)
        OUT(8) = AE(L)
        OUT(9) = DEF(L)
        OUT(10) = SUR(L)

        return

      end subroutine output

c*************************************************************************
      subroutine totm(X,N,IND,SUM,NN)

        dimension SUM(NN),IND(NN),X(N)

c
c       Sum values of 'X' specified by 'IND' over the month.
c

        I = 0
   10   I = I + 1
        J = IND(I)
        SUM(I) = SUM(I) + X(J)
        if (I .lt. NN)  go to 10

        return

      end subroutine totm

c*************************************************************************
      subroutine toty(X,N,IND,SUM,NN)

        dimension SUM(NN),IND(NN),X(N)

c
c       Sum values of 'X' specified by 'IND' over the year.
c

        I = 0
 10     I = I + 1
        J = IND(I)
        SUM(I) = SUM(I) + X(J)
        if (I .lt. NN)  go to 10

        return

      end subroutine toty

c*************************************************************************
      subroutine conv(X,NUM,MIN,MAX)

        dimension X(NUM)

c
c       Round X(I) to nearest whole number.
c

        I = MIN - 1
   10   I = I + 1
        if (X(I) .eq. 0.0)  go to 20
        Y = abs(X(I))
        J = X(I) / Y
        K = Y + 0.5
        X(I) = K * J
   20   if (I .lt. MAX) go to 10
        
        return
      end subroutine conv