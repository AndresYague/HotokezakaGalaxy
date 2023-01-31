!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! Calculate radioactive nuclei and abundances for a given timeDistance.in  !!!
!!! files. Also takes a list of taus for the calculation.                    !!!
!!!                                                                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM HotokezakaCalc
    USE MPI
    IMPLICIT NONE

    ! MPI variables
    INTEGER::rank, nProc, ierror

    ! Program variables
    DOUBLE PRECISION, ALLOCATABLE::taus(:), tEvents(:), dEvents(:)
    DOUBLE PRECISION, ALLOCATABLE::tArray(:), abundArray(:), prodFactor(:)
    DOUBLE PRECISION::alpha, vt, diffCoef, hscale
    INTEGER::uni, uni2, nTau, nTimes, lenEvent, nEvents, ii, jj
    INTEGER::lenFactor
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
    uni = 16; uni2 = 17
    OPEN(UNIT = uni, FILE = "input/tauList.in")

    READ(uni, *) nTau
    ALLOCATE(taus(nTau))
    READ(uni, *) taus

    CLOSE(UNIT = uni)

    ! Read alpha and vt
    OPEN(UNIT = uni, FILE = "input/inputHotokezakaGalaxy.in")

    READ(uni, *) alpha
    READ(uni, *) vt

    CLOSE(UNIT = uni)

    ! Allocate and read tArray, allocate abundArray
    OPEN(UNIT = uni, FILE = "input/tArray.in")

    READ(uni, *) nTimes
    ALLOCATE(tArray(nTimes), abundArray(nTimes))
    READ(uni, *) tArray

    CLOSE(UNIT = uni)

    ! Open output files. Write the rank to the filename.
    WRITE(sRank, '(I5)') rank
    OPEN(UNIT = uni2, FILE = "output/Output"//TRIM(ADJUSTL(sRank))//".out")

    ! Write temporal array first and leave a space
    WRITE(uni2, *) tArray
    WRITE(uni2, *)

    ! Read the events.in (which contains hscale)
    OPEN(UNIT = uni, FILE = "input/events.in", FORM = "UNFORMATTED", &
        ACCESS = "STREAM")

    READ(uni) nEvents, lenEvent, hscale

    ! Calculate the diffusion coefficient multiplied by pi (3.14)
    diffCoef = alpha*vt/7*hscale/0.2*1e-3*3.14

    ! Allocate tEvents and dEvents
    ALLOCATE(tEvents(lenEvent), dEvents(lenEvent))

    ! Read times, distances, production factors, and integrate
    DO ii = 1, nEvents
        ! For one run, read tEvents and dEvents
        READ(uni) tEvents
        READ(uni) dEvents

        ! Read the production factor
        READ(uni) lenFactor
        IF (.NOT.ALLOCATED(prodFactor)) THEN
            ALLOCATE(prodFactor(lenFactor))
        END IF
        READ(uni) prodFactor

        IF (rank.EQ.MOD(ii - 1, nProc)) THEN
            ! Integrate
            WRITE(uni2, *) "#", ii

            ! Do one separate run for each tau
            DO jj = 1, nTau
                abundArray = 0.D0

                ! Calculate
                CALL decayingAbund(abundArray, tEvents, dEvents, tArray, &
                                   taus(jj), prodFactor, diffCoef, hscale)

                ! Write
                WRITE(uni2, *) taus(jj), abundArray
            END DO
        END IF
    END DO

    CLOSE(UNIT = uni)
    CLOSE(UNIT = uni2)

    ! Clean
    DEALLOCATE(taus, tEvents, dEvents, tArray, abundArray)
    IF (ALLOCATED(prodFactor)) DEALLOCATE(prodFactor)
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
    DOUBLE PRECISION::invTau, factor, dt, r2
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

            ! Stop calculating if it decayed too much
            IF ((dt * invTau).GT.100) EXIT

            IF (SIZE(prodFactor).GT.1) THEN
                abundArray(jj) = abundArray(jj) + prodFactor(iiEvent)*factor
            ELSE
                abundArray(jj) = abundArray(jj) + prodFactor(1)*factor
            END IF
        END DO
    END DO

END SUBROUTINE decayingAbund

END PROGRAM HotokezakaCalc
