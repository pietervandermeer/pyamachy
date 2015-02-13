PROGRAM test_hdf5
    USE noise_reader
    USE mask_reader
    USE smr_reader
    IMPLICIT NONE

    !
    ! variables
    !

    ! SMR database with various calibration options
    CHARACTER(LEN=256) :: smr_name = "/SCIA/SDMF31/pieter/sdmf_smr_123_.h5"
    !CHARACTER(LEN=256) :: smr_name = "/SCIA/SDMF31/pieter/sdmf_smr_full.h5"
    !CHARACTER(LEN=256) :: smr_name = "/SCIA/SDMF31/pieter/sdmf_smr_raw_.h5"
    !CHARACTER(LEN=256) :: smr_name = "/SCIA/SDMF31/pieter/sdmf_smr_12_.h5"
    !CHARACTER(LEN=256) :: smr_name = "/SCIA/SDMF31/pieter/sdmf_smr_3.h5"
    !CHARACTER(LEN=256) :: smr_name = "/SCIA/SDMF31/pieter/sdmf_smr_.h5"

    ! noise array
    REAL(8), DIMENSION(1024) :: noise
    ! mask array
    INTEGER, DIMENSION(1024) :: mask

    ! smr array, errors, variance
    REAL(8), DIMENSION(8192) :: smr, smr_errors, smr_var
    ! error code
    INTEGER :: error

    !
    ! code
    !

    !CALL read_noise(3301, 1, noise, error)
    !WRITE(*,*) noise

    !CALL read_mask(4165, mask, error)
    !CALL read_mask(4152, mask, error)
    !WRITE(*,*) mask

    CALL read_smr(4204, smr_name, smr, smr_errors, smr_var, error)
    WRITE(*,*) smr

END PROGRAM test_hdf5
