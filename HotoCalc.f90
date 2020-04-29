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
    DOUBLE PRECISION, ALLOCATABLE::taus(:), tEvents(:, :), dEvents(:, :), nothing(:)
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
    OPEN(UNIT = uni, FILE = "input/events.in", FORM = "UNFORMATTED", &
        ACCESS = "STREAM")
    
    READ(uni) nEvents, lenEvent, hscale
    
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
    ALLOCATE(nothing(lenEvent))
    
    ! Read times
    jj = 1
    DO ii = 1, nEvents
        IF (rank.EQ.MOD(ii - 1, nProc)) THEN
            READ(uni) tEvents(:, jj)
            jj = jj + 1
        ELSE
            READ(uni) nothing
        END IF
    END DO
    
    ! Read distances
    jj = 1
    DO ii = 1, nEvents
        IF (rank.EQ.MOD(ii - 1, nProc)) THEN
            READ(uni) dEvents(:, jj)
            jj = jj + 1
        ELSE
            READ(uni) nothing
        END IF
    END DO
    
    ! Now read the prodFactor
    READ(uni) lenFactor, nFactor
    
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
        READ(uni) valFactor
    END IF
    
    ! Allocate prodFactor with the reduced size redNEvents
    ALLOCATE(prodFactor(lenEvent, redNEvents))
    
    ! Read the information on the prodFactor
    IF (lenFactor.EQ.1) THEN
        prodFactor = valFactor
    ELSE
        jj = 1
        DO ii = 1, nEvents
            IF (rank.EQ.MOD(ii - 1, nProc)) THEN
                READ(uni) prodFactor(:, jj)
                jj = jj + 1
            ELSE
                READ(uni) nothing
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
        abundArray = 0.D0
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
    
    DEALLOCATE(taus, tEvents, dEvents, tArray, abundArray, prodFactor, nothing)
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
    DOUBLE PRECISION::invTau, totSum, factor, dt, r2
    INTEGER::ii, jj, tLen, iiEvent, lenEvent
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Store inverted tau
    invTau = 1.0/tau
    
    ! Store lengths
    tLen = SIZE(tArray)
    lenEvent = SIZE(tEvents)
    
    ! Now calculate the values for all points
    ii = 1
    DO iiEvent = 1, lenEvent
        ! Get the next measuring point
        DO WHILE ((ii.LT.tLen).AND.(tArray(ii).LT.tEvents(iiEvent)))
            ii = ii + 1
        END DO
        
        ! Now apply to every point in the future
        DO jj = ii, tLen
            ! Calculate dt and r2
            dt = tArray(jj) - tEvents(iiEvent)
            r2 = diffCoef*dt
            
            ! Get the factor
            factor = DEXP(-dEvents(iiEvent)**2/(4*r2) - dt*invTau)
            factor = factor/MIN((4*r2)**1.5, 8*hscale*r2)
            
            ! Check if the factor is too small
            IF (factor.LT.1.D-100) EXIT
            
            abundArray(jj) = abundArray(jj) + prodFactor(iiEvent)*factor
        END DO
    END DO
    
END SUBROUTINE decayingAbund

END PROGRAM HotoCalc
