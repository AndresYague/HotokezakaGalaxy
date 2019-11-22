!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! Calculate radioactive nuclei and abundances for a given timeDistance.in  !!!
!!! files. Also takes a list of taus for the calculation.                    !!!
!!!                                                                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM HotoCalc
    USE MPI
    IMPLICIT NONE
    
    ! MPI variables
    INTEGER::rank, nProc, ierror
    
    ! Program variables
    DOUBLE PRECISION, ALLOCATABLE::taus(:), tEvents(:, :), dEvents(:, :)
    DOUBLE PRECISION, ALLOCATABLE::tArray(:), abundArray(:, :), prodFactor(:, :)
    DOUBLE PRECISION::valFactor, alpha, vt, diffCoef, hscale
    INTEGER::uni, nTau, nTimes, lenEvent, nEvents, ii, jj, kk, redNEvents
    INTEGER::lenFactor, nFactor
    CHARACTER(5)::sRank
    LOGICAL::isMaster
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Initialize MPI
    CALL MPI_INIT(ierror)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    
    ! Identify master
    isMaster = (rank.EQ.0)
    
    ! Read tau parameters
    uni = 16
    OPEN(UNIT = uni, FILE = "input/tauList.in")
    
    READ(uni, *) nTau
    ALLOCATE(taus(nTau))
    READ(uni, *) taus
    
    CLOSE(UNIT = uni)
    
    ! Read alpha, vt, and hscale
    OPEN(UNIT = uni, FILE = "input/inputHotoGalaxy.in")
    
    READ(uni, *) alpha
    READ(uni, *) vt
    
    CLOSE(UNIT = uni)
    
    ! Now read the events.in (which contains hscale)
    OPEN(UNIT = uni, FILE = "input/events.in")
    
    READ(uni, *) nEvents, lenEvent, hscale
    ! Distribute the events among the processes to conserve memory
    redNEvents = 0
    DO ii = 1, nEvents
        IF (rank.NE.MOD(ii - 1, nProc)) CYCLE
        
        redNEvents = redNEvents + 1
    END DO
    
    ! Calculate the diffusion coefficient multiplied by pi (3.14)
    diffCoef = alpha*vt/7*hscale/0.2*1e-3*3.14
    
    ! Allocate tEvents and dEvents with the reduced size redNEvents
    ALLOCATE(tEvents(lenEvent, redNEvents))
    ALLOCATE(dEvents(lenEvent, redNEvents))
    
    ! Now read the information
    jj = 1
    DO ii = 1, nEvents
        IF (rank.EQ.MOD(ii - 1, nProc)) THEN
            READ(uni, *) tEvents(:, jj)
            jj = jj + 1
        ELSE
            READ(uni, *)
        END IF
    END DO
    
    jj = 1
    DO ii = 1, nEvents
        IF (rank.EQ.MOD(ii - 1, nProc)) THEN
            READ(uni, *) dEvents(:, jj)
            jj = jj + 1
        ELSE
            READ(uni, *)
        END IF
    END DO
    
    CLOSE(UNIT = uni)
    
    ! Now read the prodFactor
    OPEN(UNIT = uni, FILE = "input/prodFactor.in")
    
    READ(uni, *) lenFactor, nFactor
    
    ! If there is only one production factor, then we fill the array with that
    ! value. We know is only one if lenFactor = nFactor = 1
    IF ((lenFactor.NE.1).AND.(lenFactor.NE.lenEvent)) THEN
        PRINT*, "Error! There should be the same number of events"
        PRINT*, " and production factors!"
        STOP
    END IF
    IF ((nFactor.NE.1).AND.(nFactor.NE.nEvents)) THEN
        PRINT*, "Error! There should be the same number of events"
        PRINT*, " and production factors!"
        STOP
    END IF
    IF (lenFactor.EQ.1) THEN
        READ(uni, *) valFactor
    END IF
    
    ! Allocate prodFactor with the reduced size redNEvents
    ALLOCATE(prodFactor(lenEvent, redNEvents))
    
    ! Now read the information
    IF (lenFactor.EQ.1) THEN
        prodFactor = valFactor
    ELSE
        jj = 1
        DO ii = 1, nEvents
            IF (rank.EQ.MOD(ii - 1, nProc)) THEN
                READ(uni, *) prodFactor(:, jj)
                jj = jj + 1
            ELSE
                READ(uni, *)
            END IF
        END DO
    END IF
    
    CLOSE(UNIT = uni)
    
    ! Allocate and read tArray, allocate abundArray
    OPEN(UNIT = uni, FILE = "input/tArray.in")
    
    READ(uni, *) nTimes
    ALLOCATE(tArray(nTimes), abundArray(nTimes, redNEvents))
    READ(uni, *) tArray
    
    CLOSE(UNIT = uni)
    
    ! Open output files. Write the rank to the filename.
    WRITE(sRank, '(I5)') rank
    OPEN(UNIT = uni, FILE = "output/Output"//TRIM(ADJUSTL(sRank))//".out")
    
    ! Write temporal array first and leave a space
    WRITE(uni, *) tArray
    WRITE(uni, *)
    
    ! Do one separate run for each tau
    DO ii = 1, nTau
        WRITE(uni, *) "#", taus(ii)
        
        ! Calculate
        kk = 1
        DO jj = 1, nEvents
            ! Divide events evenly
            IF (rank.NE.MOD(jj - 1, nProc)) CYCLE
            
            CALL decayingAbund(abundArray(:, kk), tEvents(:, kk), &
                               dEvents(:, kk), tArray, taus(ii), &
                               prodFactor(:, kk), diffCoef, hscale)
            
            kk = kk + 1
        END DO
        
        ! Write
        kk = 1
        DO jj = 1, nEvents
            ! Divide events evenly
            IF (rank.NE.MOD(jj - 1, nProc)) CYCLE
            
            WRITE(uni, *) jj, abundArray(:, kk)
            kk = kk + 1
        END DO
    END DO
    
    CLOSE(uni)
    
    DEALLOCATE(taus, tEvents, dEvents, tArray, abundArray, prodFactor)
    CALL MPI_FINALIZE(ierror)
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This function generates the decaying abundances array.                   !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -abundArray, the empty abundances array.                                 !!!
!!! -tEvents, the array with the polluting events times.                     !!!
!!! -dEvents, the array with the polluting events distances.                 !!!
!!! -tArray, the array with the sampling temporal points.                    !!!
!!! -tau, the tau value.                                                     !!!
!!! -prodFactor, the array with the production factor.                       !!!
!!! -diffCoef, the diffusion coefficient times pi.                           !!!
!!! -hscale, the scale height for the galaxy model.                          !!!
!!!                                                                          !!!
!!! At the output, abundArray will be updated.                               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE decayingAbund(abundArray, tEvents, dEvents, tArray, tau, &
                         prodFactor, diffCoef, hscale)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::abundArray(:), tEvents(:), dEvents(:), tArray(:), tau
    DOUBLE PRECISION::prodFactor(:), diffCoef, hscale
    
    ! Local
    DOUBLE PRECISION, ALLOCATABLE::dt(:), r2(:)
    DOUBLE PRECISION::invTau, totSum, factor
    INTEGER::ii, jj, tLen, iiEvent, lenEvent
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Store inverted tau
    invTau = 1.0/tau
    
    ! Store lengths
    tLen = SIZE(tArray)
    lenEvent = SIZE(tEvents)
    
    ! Allocate arrays
    ALLOCATE(dt(lenEvent), r2(lenEvent))
    
    ! Now calculate the values for all points
    iiEvent = 1
    DO ii = 1, tLen
        ! Get latest event
        DO WHILE (tEvents(iiEvent + 1).LT.tArray(ii))
            iiEvent = iiEvent + 1
            
            ! Check that we don't go over the array size
            IF (iiEvent.GT.lenEvent - 1) THEN
                iiEvent = lenEvent
                EXIT
            END IF
        END DO
        
        ! Calculate dt and r2
        dt(1:iiEvent) = tArray(ii) - tEvents(1:iiEvent)
        r2(1:iiEvent) = diffCoef*dt(1:iiEvent)
        
        factor = 1.D0; totSum = 0.D0
        DO jj = 1, iiEvent
            factor = DEXP(-dEvents(jj)**2/(4*r2(jj)) - dt(jj)*invTau)
            factor = factor/MIN((4*r2(jj))**1.5, 8*hscale*r2(jj))
            
            totSum = totSum + prodFactor(jj)*factor
        END DO
        abundArray(ii) = totSum
    END DO
    
    DEALLOCATE(dt, r2)
    
END SUBROUTINE decayingAbund

END PROGRAM HotoCalc
