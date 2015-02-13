MODULE smr_reader
  USE hdf5

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE read_from_dset(group_id, dset_name, index, buf, error)
    INTEGER(HID_T), INTENT(IN) :: group_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(HSIZE_T), INTENT(IN) :: index
    REAL(8), DIMENSION(8192), INTENT(OUT) :: buf
    INTEGER, INTENT(OUT) :: error

    INTEGER(HSIZE_T), DIMENSION(1:2) :: dim, dim_, dim_out, maxdim_out
    INTEGER(HSIZE_T), DIMENSION(1:2) :: file_slab_start, file_slab_count, file_slab_stride, file_slab_block, mem_2d_start
    INTEGER(HSIZE_T), DIMENSION(1) :: mem_slab_start, mem_slab_count
    INTEGER(HID_T) :: dset_id, memspace_id, filespace_id, xfer_prp

    dim=(/8192,1/)
    dim_=(/1,8192/)

    CALL h5dopen_f(group_id, TRIM(dset_name), dset_id, error)
    IF (error < 0) THEN
        WRITE(*,*) "failed to open dataset:",dset_name, "error=",error
        RETURN
    ENDIF
    WRITE(*,*) "dataset opened", dset_id, "error=", error

    CALL h5screate_simple_f(2, dim, memspace_id, error) 
    WRITE(*,*) "memspace created", memspace_id, "error=", error

    CALL h5dget_space_f(dset_id, filespace_id, error)
    WRITE(*,*) "got filespace", filespace_id, "error=", error

    !
    ! Select hyperslab in dataset (a single orbit x 8192 pixels).
    !

    file_slab_start = (/0,INT(index)/)
    file_slab_count = (/1,1/)
    file_slab_stride = (/1,1/)
    file_slab_block = dim
    CALL h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, file_slab_start, file_slab_count, error, file_slab_stride, file_slab_block)
    WRITE(*,*) "file hyperslab selected, error=", error

    !
    ! Select hyperslab in memory (the entire buffer).
    !

!    mem_2d_start = (/0,0/)
!    mem_slab_start = 0
!    mem_slab_count = 8192
!    CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, mem_2d_start, file_slab_count, error, file_slab_stride, file_slab_block) 
!    WRITE(*,*) "mem hyperslab selected, error=", error

    CALL h5sget_simple_extent_dims_f(memspace_id, dim_out, maxdim_out, error)
    WRITE(*,*) "dim_out=", dim_out, "maxdim_out=", maxdim_out, "error=", error

    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, dim, error, file_space_id=filespace_id, mem_space_id=memspace_id)
    WRITE(*,*) "hyperslab read", dset_id, "error=", error

    IF (error < 0) THEN
        WRITE(*,*) "failed to read hyperslab!"
        RETURN
    ENDIF

    CALL h5sclose_f(memspace_id, error)
    CALL h5sclose_f(filespace_id, error)
    CALL h5dclose_f(dset_id, error)

  END SUBROUTINE read_from_dset

  ! read sun-mean-reference spectrum for a single orbit (channel 8)
  ! orbit: absolute orbit [1..53000]
  ! file_name: name of existing smr database
  ! smr: the smr data
  ! errors: errors on the smr data
  ! variances: variances on the smr data
  ! error: error return code
  SUBROUTINE read_smr(orbit, file_name, smr, errors, variances, error) 
    INTEGER, INTENT(IN) :: orbit
    CHARACTER(LEN=256), INTENT(IN) :: file_name
    REAL(8), DIMENSION(8192), INTENT(OUT) :: smr, errors, variances
    INTEGER, INTENT(OUT) :: error

    INTEGER(HID_T) :: file_id
    CHARACTER(LEN=16) :: group_name, dataset_name
    INTEGER(HSIZE_T) :: index

    CALL h5open_f(error)
    IF (error < 0) THEN
        WRITE(*,*) "unable to open hdf5 lib!"
        RETURN
    ENDIF

    CALL h5fopen_f(TRIM(file_name), H5F_ACC_RDONLY_F, file_id, error)
    IF (error < 0) THEN
        WRITE(*,*) "unable to open file:",file_name,"error=",error
        RETURN
    ENDIF
    WRITE(*,*) "opened", file_id, "error=", error

    !
    ! find index in the orbit list
    !

    CALL get_orbit_index(file_id, orbit, index)

    WRITE(*,*) "index= ", index

    IF (index < 0) THEN
        WRITE(*,*) "orbit ", orbit," not found!"
        RETURN
    ENDIF

    dataset_name = "smr"
    CALL read_from_dset(file_id, dataset_name, index, smr, error)
    IF (error < 0) THEN
        WRITE(*,*) "failed to read from:", dataset_name
        RETURN
    ENDIF
    dataset_name = "smrError"
    CALL read_from_dset(file_id, dataset_name, index, errors, error)
    IF (error < 0) THEN
        WRITE(*,*) "failed to read from:", dataset_name
        RETURN
    ENDIF
    dataset_name = "smrVariance"
    CALL read_from_dset(file_id, dataset_name, index, variances, error)
    IF (error < 0) THEN
        WRITE(*,*) "failed to read from:", dataset_name
        RETURN
    ENDIF

    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

  END SUBROUTINE read_smr

  ! orbit: absolute orbit [1..53000]
  ! index: index in orbit list (and other datasets)
  SUBROUTINE get_orbit_index(group_id, orbit, index)
    INTEGER, INTENT(IN) :: orbit
    INTEGER(HSIZE_T), INTENT(OUT) :: index
    INTEGER, DIMENSION(:), ALLOCATABLE :: orbit_list

    INTEGER(HID_T) :: group_id, dset_id, space_id, filespace_id, xfer_prp
    INTEGER :: error, ndims
    INTEGER(HSIZE_T), DIMENSION(1) :: dim, maxdim
    LOGICAL :: orbit_found = .FALSE.

    CALL h5dopen_f(group_id, "orbitList", dset_id, error)

    CALL h5dget_space_f(dset_id, space_id, error)

    CALL h5sget_simple_extent_dims_f(space_id, dim, maxdim, error)

    ALLOCATE(orbit_list(dim(1)))

    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, orbit_list, dim, error)

    IF (error < 0) THEN
        WRITE(*,*) "failed to read orbit list!"
        RETURN
    ENDIF

    CALL h5sclose_f(space_id, error)

    CALL h5dclose_f(dset_id, error)

    !
    ! search the orbit list
    !

    DO index=1, dim(1)
        !WRITE(*,*) orbit_list(index)
        IF (orbit_list(index) .EQ. orbit) THEN
            orbit_found = .TRUE.
            EXIT
        ENDIF 
    ENDDO

    DEALLOCATE(orbit_list)

    IF (orbit_found) THEN
        ! adjust for c-style indexing..
        index = index -1    
    ELSE
        index = -1
    ENDIF

  END SUBROUTINE get_orbit_index

END MODULE smr_reader
