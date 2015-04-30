MODULE noise_reader
  USE hdf5

  IMPLICIT NONE

!  CHARACTER(LEN=256) :: file_name = "/SCIA/SDMF31/pieter/noise.h5"
  CHARACTER(LEN=256) :: file_name = "/SCIA/SDMF31/pieter/noise30.h5"

  CONTAINS

  ! threshold noise numbers to flags
  SUBROUTINE threshold_noise(data, threshold, flags)
    REAL(8), DIMENSION(1024), INTENT(IN) :: data
    REAL(8), INTENT(IN) :: threshold
    LOGICAL, DIMENSION(1024), INTENT(OUT) :: flags

    flags = data < threshold
  END SUBROUTINE threshold_noise

  ! read a single orbit of the noise product (channel 8 only)
  ! orbit: absolute orbit [1..53000]
  ! pet_id: 0,1,2 -> (0.125s, 0.5s, 1.0s)
  SUBROUTINE read_noise(orbit, pet_id, data, error) 
    INTEGER, INTENT(IN) :: orbit, pet_id
    REAL(8), DIMENSION(1024), INTENT(INOUT) :: data
    INTEGER, INTENT(OUT) :: error

    INTEGER(HID_T) :: file_id, group_id, dset_id, memspace_id, filespace_id, xfer_prp
    CHARACTER(LEN=16) :: group_name, dataset_name
    INTEGER(HSIZE_T), DIMENSION(1:2) :: dim, dim_out, maxdim_out
    INTEGER(HSIZE_T), DIMENSION(1:2) :: file_slab_start, file_slab_count, file_slab_stride, file_slab_block, mem_2d_start
    INTEGER(HSIZE_T), DIMENSION(1) :: mem_slab_start, mem_slab_count
    INTEGER(HSIZE_T) :: index

    WRITE(*,*) "read(", orbit, pet_id, ")"

    CALL h5open_f(error)
    CALL h5fopen_f(TRIM(file_name), H5F_ACC_RDONLY_F, file_id, error)

    WRITE(*,*) "opened", file_id, "error=", error

    if (pet_id .eq. 0) then
        group_name = 'pet0.125'
    else if (pet_id .eq. 1) then
        group_name = 'pet0.5'
    else if (pet_id .eq. 2) then
        group_name = 'pet1.0'
    end if
    WRITE(*,*) "group=", group_name

    call h5gopen_f(file_id, group_name, group_id, error)
    WRITE(*,*) "group opened", group_id, "error=", error

    dataset_name = "noise"
    CALL h5dopen_f(group_id, TRIM(dataset_name), dset_id, error)
    WRITE(*,*) "dataset opened", dset_id, "error=", error

    dim=(/1024,1/)

    CALL h5screate_simple_f(2, dim, memspace_id, error) 
    WRITE(*,*) "memspace created", memspace_id, "error=", error

    CALL h5dget_space_f(dset_id, filespace_id, error)
    WRITE(*,*) "got filespace", filespace_id, "error=", error

    !
    ! find index in the orbit list
    !

    CALL get_orbit_index(group_id, orbit, index)

    WRITE(*,*) "index ", index

    IF (index < 0) THEN
        WRITE(*,*) "orbit ", orbit," not found!"
        RETURN
    ENDIF

    !
    ! Select hyperslab in dataset (a single orbit x 1024 pixels).
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

    mem_2d_start = (/0,0/)
    mem_slab_start = 0
    mem_slab_count = 1024
    CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, mem_2d_start, file_slab_count, error, file_slab_stride, file_slab_block) 
    WRITE(*,*) "mem hyperslab selected, error=", error

    CALL h5sget_simple_extent_dims_f(memspace_id, dim_out, maxdim_out, error)
    WRITE(*,*) "dim_out=", dim_out, "maxdim_out=", maxdim_out, "error=", error

    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dim, error, file_space_id=filespace_id, mem_space_id=memspace_id)
    WRITE(*,*) "hyperslab read", dset_id, "error=", error

    IF (error < 0) THEN
        WRITE(*,*) "failed to read hyperslab!"
        RETURN
    ENDIF

    CALL h5sclose_f(memspace_id, error)
    CALL h5sclose_f(filespace_id, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(group_id, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)

  END SUBROUTINE read_noise

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

    CALL h5dopen_f(group_id, "orbits", dset_id, error)

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

END MODULE noise_reader
