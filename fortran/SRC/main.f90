PROGRAM test_hdf5
    USE noise_reader
    USE mask_reader
    USE smr_reader
    USE dark_reader
    USE analogoffset_reader

    IMPLICIT NONE

    !
    ! variables
    !

    CHARACTER(LEN=256) :: intmonthlies_name = "/SCIA/SDMF31/pieter/interpolated_monthlies_long__.h5"
    !CHARACTER(LEN=256) :: intmonthlies = "/SCIA/SDMF31/pieter/interpolated_monthlies_short__.h5"

    ! dark database for short or long exposure times, (long for earthshine, short for smr)
    CHARACTER(LEN=256) :: dark_name = "/SCIA/SDMF31/pieter/vardark_long__.h5"
    !CHARACTER(LEN=256) :: dark_name = "/SCIA/SDMF31/pieter/vardark_short__.h5"

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
    ! dark array, errors
    REAL(8), DIMENSION(1024) :: dark, dark_errors, analog_offset, ao_errors
    REAL(8) :: orbit_phase
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

    CALL read_dark(9000, orbit_phase, dark_name, dark, dark_errors, error)
    WRITE(*,*) dark
    WRITE(*,*) dark_errors

    CALL read_analogoffset(9000, intmonthlies_name, analog_offset, ao_errors, error)
    WRITE(*,*) analog_offset
    WRITE(*,*) ao_errors

END PROGRAM test_hdf5
