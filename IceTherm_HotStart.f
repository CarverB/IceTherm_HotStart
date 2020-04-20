C      IceTherm_HotStart is a modifed version of the thermal model 
C      from Bierson et al. 2018 https://doi.org/10.1016/j.icarus.2018.03.007
C      This version has been modified to ignore porosity in the ice shell
C      and include the option to set the initial condition to that of a 
C      liquid ocean insted of an ice shell.
C      
C      You are welcome to use or adapt this code providing you cite
C      C. J. Bierson, F. Nimmo, and S. A Stern, Evidence for a hot start 
C      and early ocean on Pluto, Nature Astronomy, 2020

            implicit none
          ! Accocate arrarys
          REAL*8, ALLOCATABLE :: T(:),DTT(:),R(:),TP(:),
     >     DRI(:),PSI(:),CPI(:),
     >     DTDTA(:),FMEL(:),HEAT(:),
     >     DUMMY1D(:),
     >     PP(:),GRV(:), MASS_SUM(:),
     >     DENI(:),CONDI(:), flux(:), TMELI(:)
            
          

            
          REAL*8 ::RHO, RHOC, CP, CPC, DR, DRC,
     >      COND, CONDC, ammoniaX, TMEL, DTMEL, XLH, RAD, RADC,
     >      TSTOP, TIM, DTS, TINIT, TINIT2, TCUTF, 
     >      XMASSC, FSRAD, CPM, TS, HEAT0,
     >      H1, H2, FBR, FOUT, EMEL, HFAC,
     >      FS, DELTARHO, HINT,
     >      MAXDTS, MelDRHO, Layer
     
         REAL*8 :: HPROD !heat production function
     
          !Physical constants
          
          REAL*8 :: PI
          
          
          INTEGER :: IMAX, IMAXC, IMAXT, I, J, K, KK, ILOG
          INTEGER :: IFRZ, IMEL
          INTEGER*8 :: KMAX
          
          ! coeffs for themral dependance of conducivity
          REAL*8 :: CONDa0, CONDa1
          
          ! Coeffs for calculating the melting temperature
          REAL*8 :: MEL1,MEL2,MEL3,MEL4,MEL5,MEL6
          
          LOGICAL :: CONDT, FixedTmel, CONSERVEMASS, HotStart
          CHARACTER(LEN=40) OUTFILE, INFILE

            ! This sets the variables that can be contained within the input file.
          NAMELIST/INP/COND,RHO,CP,CONDC,RHOC,CPC,
     C      TMEL,DTMEL,XLH,RAD,RADC,TSTOP,HFAC,TS,ILOG,
     C      TINIT,TINIT2,TCUTF,IMAX,IMAXC,KMAX,CONDT,
     C      FixedTmel, TMEL, ammoniaX, CONSERVEMASS, HotStart
          
C
C         Read optional output filename        
          IF (IARGC() >0) then
             call GETARG(1,OUTFILE)
          ELSE
             OUTFILE='IceTherm.out'
          ENDIF

          ! Open output file
          OPEN(UNIT=24,FILE=OUTFILE)
          
          INFILE='IceTherm.inp'
C

C
          ! Read input file
          CALL READINPUT(INFILE, RAD, RADC, HFAC, TSTOP, 
     >     CP, CPC, COND, CONDC, RHO, RHOC,
     >     TMEL, TS, ammoniaX, TINIT, TINIT2, IMAX, IMAXC,
     >     CONDT, FixedTmel, CONSERVEMASS, MAXDTS, HotStart)
          

        
          
          ! if T dependance of conductivity, k(T)=conda0+conda1/T
          CONDa0=0.4685
          CONDa1=488.12 
            
C
C         TMEL IS MELTING POINT, DTMEL IS FINITE INTERVAL OVER
C         WHICH MELTING TAKES PLACE, XLH IS LATENT HEAT
  
          
          !Coeffecents for pressure and amonia melt temperature
          Mel1=273.1
            Mel2=7.95E-8
            Mel3=9.6E-17
            mel4=53.8
            mel5=650.0
            mel6=4.0E-8
            
          
c          tmel=176.
          DTMEL=3.
          XLH=3.33E5
c          xlh=1.3e5
C

C
          ILOG=0
c          HFAC=0.001
c
          MelDRHO=1000.0-RHO ! Density change when ice melts

C

c          TINIT=170.
c          TINIT2=1000.
C
C         CAN ARTIFICIALLY STOP HEAT FLUX OUT OF CORE
C         TO VERIFY STEFAN SOLUTION 

          !TCUTF=30.E6
          TCUTF=1.E30
C
                    
C
C         NUMBER OF NODES IN OUTER AND INNER LAYERS 
C         The sum of these can not exceed the memory allocated above


          IMAXT=IMAX+IMAXC
            ! Accocate arrays for each layer
           ALLOCATE(T(IMAXT),DTT(IMAXT),R(IMAXT),TP(IMAXT),
     >       DRI(IMAXT),PSI(IMAXT),CPI(IMAXT),
     >       DTDTA(IMAXT),FMEL(IMAXT),HEAT(IMAXT),DUMMY1D(IMAXT),
     >       MASS_SUM(IMAXT),
     >       PP(IMAXT),GRV(IMAXT),
     >       DENI(IMAXT),CONDI(IMAXT), flux(IMAXT), TMELI(IMAXT))
C
C        Max NUMBER OF TIMESTEPS

          KMAX=200000000



          ! Start output file
          ! Write parameter values to output file
          WRITE(24,INP)

          WRITE(24,*) "Log10 time (yrs),    R (m),    T (K),    ",
     >     "Melt Fraction,    Density (kg/m^3),    ",
     >     "Flux (W/m^2)"
          



          CPM=CP+(XLH/DTMEL)

C
C         Setup grid

          DR=(RAD-RADC)/(FLOAT(IMAX)-0.5)
          DRC=RADC/(0.5+FLOAT(IMAXC-1))
          PI=4.*ATAN(1.)

          XMASSC=(4.*PI*(RADC**3.)*RHOC/3.)

C
C         INITIALIZE arrays

          HEAT0=RHO*((CP*DTMEL)+XLH)
          FMEL=0.0
          HEAT=HEAT0
          
          DO I=1,IMAXT
           IF(I.LE.IMAXC)THEN
               R(I)=FLOAT(I-1)*DRC
                  DRI(I)=DRC
           ELSE
               R(I)=RADC+(DR/2.)+(FLOAT(I-IMAXC-1)*DR)
                  DRI(I)=DR
               
           ENDIF
           T(I)=TINIT
           IF(TINIT2.NE.0.AND.I.LE.IMAXC)T(I)=TINIT2
           IF (HOTSTART) THEN
               IF(I.GE.IMAXC) FMEL(I)=1.0 ! Hot Start
               IF(I.GE.IMAXC) HEAT(I)=0.0
               !IF(I.GE.IMAXC)HEAT(I)=HEAT0
           ELSE
               IF(I.GE.IMAXC)HEAT(I)=HEAT0
           ENDIF
          END DO
          
  
C         INITIALIZE density
!$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(I)                  
           DO I=1,IMAXT
             DENI(I)=RHO

             IF (FMEL(I) .GT. 0.0) THEN
                 DENI(I)=RHO+MelDRHO*FMEL(I)
             ENDIF
             
             IF (CONDT) THEN
                   COND=CONDa0 + CONDa1/T(I)
             ENDIF
             
             CONDI(I)=COND
             CPI(I)=CP

             IF(I.LE.IMAXC)THEN
                  DENI(I)=RHOC
                  CONDI(I)=CONDC
                  CPI(I)=CPC
             ENDIF
             
            
             
           END DO
!$OMP END PARALLEL DO          

           MASS_SUM(1)=DENI(1)*4.0/3.0*PI*(DRI(1)/2.0)**3
           
           Layer=(DENI(1)*4.0/3.0*PI* 
     >          ((R(2)-DRI(2)/2.0)**3-(R(1)+DRI(1)/2.0)**3))
           
           DO I=2,IMAXT
               ! this isn't quite right as the mass is not evenly distributed above and below the center point
           
               MASS_SUM(I)=MASS_SUM(I-1) + Layer/2 !add upper half of previous layer
           
               Layer=(DENI(I)*4.0/3.0*PI* 
     >          ((R(I)+DRI(I)/2.0)**3-(R(I)-DRI(I)/2.0)**3))
               
               MASS_SUM(I)=MASS_SUM(I) + Layer/2 !add lower half of this layer
               
               
           ENDDO


      !PRECALCULATE GRAVITY AND PRESSURE
      CALL PRESSURE_GRAV(R, DRI, DENI, PP, MASS_SUM, GRV, IMAXT)
            
       !write(*,*) PP
       !Initalize Melting Temperature
      IF (FixedTmel) THEN
           TMELI=TMEL
      ELSE
           TMELI=Mel1-Mel2*PP - Mel3*PP*PP - 
     >     Mel4*ammoniaX - Mel5*ammoniaX*ammoniaX - Mel6*ammoniaX*PP
      ENDIF


      IF (HotStart) THEN
          DO I=IMAXC+1,IMAXT
              T(I)=TMELI(I)
          ENDDO

          ! Make Cold lid
          DO I=IMAXT-5,IMAXT
                T(I)=TS
                FMEL(I)=0.
                HEAT(I)=HEAT0
                DENI(I)=RHO
          ENDDO
      endif
  
      ! WRITE(*,*) R, DRI, DENI, PP, GRV,TMELI
     

      !START OF MAIN TIME LOOP

      KK=0
      IFRZ=0
      TIM=0
       
      PSI=1.0 ! 1=no density change
      
      DO K=1,KMAX
       
            !Reset Variables
            DTT=0.0
            TP=T
        
            ! calculate timestep CFL 
            ! DTS = timestep in seconds
            DTS=MIN(MINVAL(0.3*DRI**2*(DENI*CPI)/CONDI), MAXDTS)
            
         IF (.NOT. DTS==DTS) THEN ! Check for NAN
             WRITE(*,*) "DTS =NaN. Stopping"
             STOP
         ELSEIF (DTS<0.0) THEN
              WRITE(*,*) "DTS <0. Stopping"
              !WRITE(*,*) K,DRI, DENI, CONDI
             STOP                
         ENDIF
         
            !write(*,*) DTS
            
            TIM=TIM+DTS/(3.15E7) !time in years
            !TIM=FLOAT(K-1)*DTS/(3.15E7)
            IF(TIM.GE.TSTOP)GOTO 109
            
            HINT=HFAC*HPROD(TIM)
            FSRAD=HINT*XMASSC/(4.*PI*RAD*RAD)
            
          !Update pressure and melting temperature
          
          ! Note that this melt temperature change could cause phase changes
          ! This is currently not included in the model (because this change is tiny)
            IF (.NOT. FixedTmel .AND. CONSERVEMASS) THEN
                 TMELI=Mel1-Mel2*PP - Mel3*PP*PP - 
     >             Mel4*ammoniaX - Mel5*ammoniaX*ammoniaX - 
     >             Mel6*ammoniaX*PP
            ENDIF
               
C           NOW INVESTIGATE WHETHER MELTING STARTS
C           OR CONTINUES THIS TIMESTEP and
C           where the ice-ocean interface is



         IMEL=0
         DO I=IMAXC+1,IMAXT
              
            IF(FMEL(I).GE.1.0 .AND. FMEL(I+1).LE.0.0 .AND.IFRZ.EQ.1)THEN
                   IMEL=I
                   EXIT
            ELSEIF(FMEL(I).LT.1.0)THEN
                    IMEL=I
                    EXIT
             ENDIF
         
         END DO
            
            
          ! Temperature dependant thermal conductivity  
          IF (CONDT) THEN
               CONDI(IMAXC+1:IMAXT)=CONDa0 + CONDa1/T(IMAXC+1:IMAXT)
          ELSE
               CONDI(IMAXC+1:IMAXT)=COND
          ENDIF

        

C
C           Solve thermal diffusion equation
!$OMP PARALLEL DO PRIVATE(H1, H2, I), DEFAULT(SHARED)
            DO I=2,IMAXT-1
            
            
            IF (I.GE.(IMAXC).AND.I.LE.(IMEL+1)) THEN
               CYCLE
            ENDIF

              
             H1=(R(I)-DRI(I)/2)**2*(T(I)-T(I-1))/
     >         (DRI(I)/CONDI(I)+DRI(I-1)/CONDI(I-1))
     
             H2=(R(I)+DRI(I)/2)**2*(T(I+1)-T(I))
     >         /(DRI(I)/CONDI(I)+DRI(I+1)/CONDI(I+1))

             DTT(I)=2/(DENI(I)*CPI(I)*DRI(I)*R(I)*R(I))*(H2-H1)
            
            IF (I.LE.IMAXC) THEN
               DTT(I)=DTT(I)+(HINT/CPI(I))
            ENDIF
              
            END DO
!$OMP END PARALLEL DO         
     
               
C

         
C
!         Solve for regions at the ice-ocean interface
         
         I=IMAXC
         
         H1=(R(I)-DRI(I)/2)**2*(T(I)-T(I-1))/
     >          (DRI(I)/CONDI(I)+DRI(I-1)/CONDI(I-1))
     
            IF (SUM(FMEL) .LE. 0.01) THEN ! No ocean
         H2=(R(I)+DRI(I)/2.0)**2*(T(I+1)-T(I))
     >         /(DRI(I)/CONDI(I)+DRI(I+1)/CONDI(I+1))
            ELSE
             H2=(R(I)+DRI(I)/2.0)**2*(T(I+1)-T(I))*CONDI(I) !T(I+1) should be the melt temperature
     >         /(DRI(I)/2.0)
         ENDIF

         DTT(I)=2.0/(DENI(I)*CPI(I)*DRI(I)*R(I)*R(I))*(H2-H1)
                                
         DTT(I)=DTT(I)+(HINT/CPI(I))
         
         FBR=2.0*H2/(R(IMEL)**2)
         
         ! IMEL+1
         I=IMEL+1
         IF (SUM(FMEL) .LE. 0.01) THEN ! No ocean
               H1=(R(I)-DRI(I)/2.0)**2*(T(I)-T(I-1))/
     >        (DRI(I)/CONDI(I)+DRI(I-1)/CONDI(I-1))
            ELSE
                H1=(R(I)-DRI(I)/2.0)**2*(T(I)-T(I-1))*CONDI(I)/
     >       (DRI(I)/2.0)
         ENDIF
       
          H2=(R(I)+DRI(I)/2.0)**2*(T(I+1)-T(I))
     >         /(DRI(I)/CONDI(I)+DRI(I+1)/CONDI(I+1))
       
          DTT(I)=2.0/(DENI(I)*CPI(I)*DRI(I)*R(I)*R(I))*(H2-H1)
           
          FOUT=2.0*H1/(R(IMEL)**2)
            ! IMEL

          I=IMEL
         
          IF(TIM.GE.TCUTF)FBR=0.

         
          EMEL=(FOUT-FBR)/DRI(IMEL)
          DTT(IMEL)=EMEL/(DENI(IMEL)*CPI(IMEL))
         

          DO I=2,IMAXT-1
              TP(I)=T(I)+(DTT(I)*DTS)
          END DO

C
C         CHECK TO SEE IF MELTING STARTS
         
         ! New Melting/Freezing
         
         IF(TP(IMEL).GE.(TMELI(IMEL)-DTMEL).OR. ANY(FMEL.GT.0))THEN
           I=IMEL
           ! Normal case where the energy avaiable for melt
            ! comes from the flux difference between bottom and top
            ! of the ocean
           EMEL=-(FBR-FOUT)*DTS/DRI(IMEL)
           IF(FMEL(IMEL) .LE. 0.0 .AND. 
     >           TP(IMEL).GT.(TMELI(IMEL)-DTMEL)) THEN
               ! Case when the energy was already sent into warming a 
               ! layer that should now melt
               EMEL=(TP(IMEL)-(TMELI(IMEL)-DTMEL))*DENI(IMEL)*CPM 
           ENDIF
           
           IF (EMEL.LT.0.0) THEN !check if freezing for ocean top determination
            IFRZ=1
           ELSEIF (EMEL.GT.0.0) THEN
            IFRZ=0
           ENDIF
         

           HEAT(IMEL)=HEAT(IMEL)-EMEL
           
            FMEL(IMEL)=1.-(HEAT(IMEL)/HEAT0) !!!! melting
              

           
            IF(FMEL(IMEL).GT.1.0) THEN !all melted
              !Set state
              FMEL(IMEL)=1.0
              HEAT(IMEL)=0.0
             ELSEIF (FMEL(IMEL).LT.0.0) THEN !all frozen
             
                !set state
                FMEL(IMEL)=0.0
                HEAT(IMEL)=HEAT0
             ENDIF
              
              TP(IMEL)=(TMELI(IMEL)-DTMEL)+(FMEL(IMEL)*DTMEL)
              ! Calculate values for mass conservation
              DELTARHO=(RHO+MelDRHO*FMEL(IMEL))-DENI(IMEL)            
              PSI(IMEL)=(DENI(IMEL))/((RHO+MelDRHO*FMEL(IMEL)))
              DENI(IMEL)=DENI(IMEL)+DELTARHO
         
         ENDIF
         
            
        
        IF (CONSERVEMASS) THEN 
            !recalculate layer thickness
            CALL VOLUMECHANGE(R, DRI, PSI, IMAXT)
            PSI=1.0 !reset PSI
            !Recalculate Pressure
            CALL PRESSURE_GRAV(R, DRI, DENI, PP, MASS_SUM, GRV, IMAXT)
        ENDIF


C
C           UPDATE ALL POINTS

! Enforce boundary conditions
            TP(IMAXT)=TS
            TP(1)=TP(2)

            IF (FMEL(IMAXT) >0.0) THEN
                FMEL(IMAXT) =0.0 
                IFRZ = 1.0

            ENDIF

            ! Update temperature
            T=TP
C
C           Output VALUES FOR PLOTTING
C           Output every 1000 timesteps
            IF(MOD(K-1,1000).EQ.0)THEN
C
C             OUTPUT TO FILE

              DO I=1,IMAXT
C                Calculate flux at I-1/2
                 if (I.GT.1) then
                 flux(I)=0.5 *(CONDI(I)+CONDI(I-1))*
     >            (T(I-1)-T(I))/(R(I)-R(I-1))
      !flux(i)=2*(T(I)-T(I-1))/
       !>           (DRI(I)/CONDI(I)+DRI(I-1)/CONDI(I-1))
                 else
                 flux(i)=0.0
                 endif

                 WRITE(24,*)LOG10(TIM),R(I),T(I),FMEL(I),DENI(I), 
     >            flux(I)
              END DO
            ENDIF
          END DO

109       CONTINUE


C         Output at code completion
!          FS=CONDI(IMAXT)*(T(IMAXT-1)-T(IMAXT))/DRI(IMAX)
!          WRITE(6,*) 'SURFACE HEAT FLUX ',FS
C
C         THICKNESS

!          THICK=RAD-(R(IMEL)+(DR*FMEL(IMEL)))
!          WRITE(6,*) 'SHELL THICKNESS ',THICK/1.E3,
!     >          R(IMEL),IMEL,FMEL(IMEL)




          END
C
C**********************************************************************
C
          FUNCTION HPROD(TIM)
        
          IMPLICIT NONE
          REAL*8:: HPROD 
          REAL*8, INTENT(IN) :: TIM 
          
          REAL*8 ::CU238, CU235, CTH, CK, HU238, HU235, HTH, HK,
     >       TU238, TU235, TTH, TK, CU238T, CU235T, CTHT, CKT    
C
C         CALCULATES HEAT PRODUCTION (W/KG) AT TIME T (YRS)
C
C         original value          CU238=16.9E-9
          CU238=19.9E-9
          CU235=5.4E-9
          CTH=38.7E-9
          CK=737.9E-9
c          ck=50.e-9

 
        
          HU238=94.65E-6
          HU235=568.7E-6
          HTH=26.38E-6
          HK=29.17E-6

          TU238=4.47E9
          TU235=7.04E8
          TTH=1.4E10
          TK=1.28E9

c          tu238=tu238*1.e4
c          tu235=tu235*1.e4
c          tth=tth*1.e4
c          tk=tk*1.e4

          CU238T=CU238*EXP(-TIM*LOG(2.)/TU238)
          CU235T=CU235*EXP(-TIM*LOG(2.)/TU235)
          CTHT=CTH*EXP(-TIM*LOG(2.)/TTH)
          CKT=CK*EXP(-TIM*LOG(2.)/TK)

          HPROD=(CU238T*HU238)+(CU235T*HU235)+
     C             (CTHT*HTH)+(CKT*HK)
          RETURN
          END
C
C*************************************************************************
C


          SUBROUTINE  VOLUMECHANGE(R, DRI, PSI, N)
          
          IMPLICIT NONE
            
            INTEGER, INTENT(IN) :: N
            
            REAL*8, INTENT(INOUT) :: R(N), DRI(N), PSI(N)         
            REAL*8 :: DELTAR, Rup
            
            !Real*8, DIMENSION(N) :: Da,Db, rlow
            Real*8 :: Da,Db, rlow
          
          INTEGER :: I, J
            
                 
       DO I=1,N
       
            IF (PSI(I).EQ.1.0) CYCLE ! No density change
    
            IF (I.LT.N) THEN
                Rup=R(I)+DRI(I)/2.0;
            ELSE
                Rup=R(I); !Surface
            ENDIF
                        
            
            
            !Rf=(Rlow*(1.0-PSI(I))+PSI(I)*Rup**3.0)**(1.0/3.0)
            !DELTAR=Rf-Rup
            
            DELTAR=(((1-DRI(I)/Rup)**3*(1.0-PSI(I))
     >      +PSI(I))**(1.0/3.0)-1.0)*Rup
 
            !Rf=((Rup**3-Rlow**3)*PSI(I)+Rlow**3)**(1.0/3.0)
            !DELTAR=Rf-Rup
            
            !write(*,*) Rf, Rup, PSI(I)
            
            DRI(I)=DRI(I)+DELTAR
            
            R(I)=R(I)+DELTAR/2.0;
            
            Da=DELTAR
            
            DO J=I+1,N

                Rlow=R(J)-DRI(J)/2.0
                
                
                Db=(((1+(Da-DRI(J))/Rlow)**3-
     >               (1-DRI(J)/Rlow)**3+1.0)**(1.0/3.0)-1.0)*Rlow
     
                R(J)=R(J)+(Da+Db)/2
                DRI(J)=DRI(J)+(Db-Da) !DRI(J) also changes
                
               Da=Db !for next layer
                
 
           ENDDO
            
            ! Vector form
            
            !This works but is marginally slower.
            
            
               ! Rup refers to layer where density change happened
               ! Rlow is an array of layer bottoms
            
c                Rlow(I+1:N)=R(I+1:N)-DRI(I+1:N)/2.0
            
c                Da(I+1:N)= (((1+(DeltaR-(Rlow(I+1:N)-Rup))/Rup)**3.0-
c     >           (1-(Rlow(I+1:N)-Rup)/Rup)**3.0+1.0)**(1.0/3.0)-1.0)*Rup
                
c                Db(I+1:N)=(((1+(Da(I+1:N)-DRI(I+1:N))/Rlow(I+1:N))**3.0-
c     >                    (1-DRI(I+1:N)/Rlow(I+1:N))**3.0+
c     >                    1.0)**(1.0/3.0)-1.0)*Rlow(I+1:N)
     
c                R(I+1:N)=R(I+1:N)+(Da(I+1:N)+Db(I+1:N))/2.0
c                DRI(I+1:N)=DRI(I+1:N)+(Db(I+1:N)-Da(I+1:N)) !DRI(J) also changes
                
            

            
            
            !IF (I.LT.N .AND. I.GT.1) THEN
            !    R(I)=R(I)+DELTAR/2.0;
            !    R(I+1:N)=R(I+1:N)+DELTAR;
            !ELSE
            !    R(I)=R(I)+DELTAR;   ! Boundary layers  
            !ENDIF

          ENDDO
            
       RETURN                   
       END 

C***********************************************************************

          SUBROUTINE PRESSURE_GRAV(R, DR, DENSITY, PRESSURE, MASS_SUM,
     >      GRAV, N)
          
          IMPLICIT NONE
          
          INTEGER, INTENT(IN) :: N
          REAL*8, INTENT(IN) :: R(N), DENSITY(N), DR(N), MASS_SUM(N)
          REAL*8, INTENT(OUT) :: PRESSURE(N), GRAV(N)
          
          INTEGER :: I
          REAL*8 :: BIGG, PI
          
          PI=4*ATAN(1.0)
          BIGG=6.6742E-11
          

C          CALCULATE PRESSURE
           PRESSURE(N)=0
           DO I=N-1,1,-1
           
             GRAV(I)=BIGG*MASS_SUM(I)/(R(I)*R(I))
             
             PRESSURE(I)=PRESSURE(I+1)+DENSITY(I)*GRAV(I)*
     C        (DR(I+1)/2.0+DR(I)/2.0)
            
            
           END DO
          
          
          RETURN 
          END

C***********************************************************************

          SUBROUTINE READINPUT(INFILE, RAD, RADC, HFAC, TSTOP, 
     >     CP, CPC, COND, CONDC, RHO, RHOC,
     >     TMEL, TS, ammoniaX, TINIT, TINIT2, IMAX, IMAXC,
     >     CONDT, FixedTmel, CONSERVEMASS, MAXDTS, HotStart)
          
          IMPLICIT NONE
          
          REAL*8, INTENT(OUT) :: RAD, RADC, HFAC, TSTOP, 
     >     CP, CPC, COND, CONDC, RHO, RHOC,
     >     TMEL, TS, ammoniaX, TINIT, TINIT2, MAXDTS
     
          INTEGER, INTENT(OUT) :: IMAX, IMAXC
          
          LOGICAL, INTENT(OUT) :: CONDT, FixedTmel, CONSERVEMASS, 
     >                            HotStart

          CHARACTER(LEN=40), INTENT(IN) :: INFILE
          
          !WRITE(*,*) INFILE
          
          OPEN(UNIT=22,FILE=INFILE)
          
          READ(22,*) !Text line
          READ(22,*) !Text line
          READ(22,*) !Text line
          READ(22,*) RAD, RADC, IMAX, IMAXC, HFAC, TSTOP, MAXDTS
          READ(22,*) !Text line
          READ(22,*) !Text line
          READ(22,*) CP, CPC, COND, CONDC, RHO, RHOC
          READ(22,*) !Text line
          READ(22,*) TMEL, TS, ammoniaX, TINIT, TINIT2
          READ(22,*) !Text line
          READ(22,*) !Text line
          READ(22,*) CONDT, FixedTmel, CONSERVEMASS, HotStart
          
          CLOSE(22)
                
          
          
          RETURN
          END









