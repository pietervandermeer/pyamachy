PROGRAM test_hdf5
    USE noise_reader
    IMPLICIT NONE

    REAL(8), DIMENSION(1024) :: noise
    INTEGER :: error

    CALL read_noise(3301, 1, noise, error)
    WRITE(*,*) noise

END PROGRAM test_hdf5
