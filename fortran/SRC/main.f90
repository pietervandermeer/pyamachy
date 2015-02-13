PROGRAM test_hdf5
    USE noise_reader
    USE mask_reader
    IMPLICIT NONE

    REAL(8), DIMENSION(1024) :: noise
    INTEGER, DIMENSION(1024) :: mask
    INTEGER :: error

    CALL read_noise(3301, 1, noise, error)
    WRITE(*,*) noise

    CALL read_mask(5400, mask, error)
    WRITE(*,*) mask

END PROGRAM test_hdf5
