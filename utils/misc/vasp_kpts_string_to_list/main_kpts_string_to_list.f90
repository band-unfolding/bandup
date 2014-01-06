program generate_suitable_SCKPTS_from_pckpts
use io_routines
use math
implicit none
character(len=127) :: input_file_prim_cell, input_file_supercell, input_file_pc_kpts, outpt_file_SC_kpts
real(kind=dp), dimension(1:3,1:3) :: dir_latt_pc, dir_latt_SC, b_matrix_pc, B_matrix_SC
real(kind=dp), dimension(1:3) :: point, previous_point, point_reduced_to_SCBZ
real(kind=dp), parameter :: tol_for_vec_equality = 1d-5
integer, dimension(:), allocatable :: n_kpts_dirs
type(vec3d), dimension(:), allocatable :: list, list_of_necessary_SCKPTS
integer :: n_nec_SCKPTS,ikpt,ikpt2,i,ios
logical :: point_already_used
character(len=254) :: file_line
    

    input_file_prim_cell = 'prim_cell_lattice.in'
    input_file_supercell = 'supercell_lattice.in'
    input_file_pc_kpts = 'KPOINTS_prim_cell.in'
    outpt_file_SC_kpts = 'KPOINTS_supercell.out'
    call read_unit_cell(input_file=input_file_prim_cell,latt=dir_latt_pc(:,:))
    call read_unit_cell(input_file=input_file_supercell,latt=dir_latt_SC(:,:))
    call get_rec_latt(latt=dir_latt_pc,rec_latt=b_matrix_pc)
    call get_rec_latt(latt=dir_latt_SC,rec_latt=B_matrix_SC)
    call define_pckpts_to_be_checked(input_file_pc_kpts, b_matrix_pc(:,:), list,n_kpts_dirs)

    n_nec_SCKPTS = 0
    allocate(list_of_necessary_SCKPTS(1:size(list)))
    do ikpt = 1, size(list_of_necessary_SCKPTS)
        list_of_necessary_SCKPTS(ikpt)%coord(:) = 0.0d0
    enddo

    do ikpt=1, size(list)

        point(:) = list(ikpt)%coord(:)
        call reduce_point_to_bz(point,B_matrix_SC,point_reduced_to_SCBZ)
        point(:) = coords_cart_vec_in_new_basis(cart_vec=point_reduced_to_SCBZ,new_basis=B_matrix_SC)

        point_already_used = .FALSE.
        if(ikpt>1)then
            do ikpt2 = ikpt-1,1,-1
                previous_point(:) = list_of_necessary_SCKPTS(ikpt2)%coord(:)
                if(same_vector(point,previous_point,tol_for_vec_equality))then
                    point_already_used = .TRUE.
                    exit
                endif
            enddo
        endif

        if(.not.point_already_used)then
            n_nec_SCKPTS = n_nec_SCKPTS + 1
            list_of_necessary_SCKPTS(n_nec_SCKPTS)%coord(:) = point(:)
        endif

    enddo

    open(unit=03, file=outpt_file_SC_kpts)
        open(unit=04, file=input_file_pc_kpts)
            read(04,'(A)')file_line
            write(03,'(A)')trim(file_line)//' (header of the input kpts file)'
        close(04)
        write(03,'(I0)')n_nec_SCKPTS
        write(03,'(A)')'Reciprocal' 
        do ikpt=1, n_nec_SCKPTS
            point(:) = list_of_necessary_SCKPTS(ikpt)%coord(:)
            write(03,'(3(2X,f10.7),2X,I1)')(point(i), i=1,3),1
        enddo

        write(03,'(A)')''
        write(03,'(A)')''
        write(03,'(A)')'! These SC-KPTS unfold into the pc-kpts listed below: '
        write(03,'(A)')'! (Fractional coords. wrt. the pc rec. latt. vectors) '
        open(unit=04, file=input_file_pc_kpts)
            ios=0
            do while(ios==0)
                read(04,'(A)',iostat=ios)file_line
                if(ios==0) write(03,'(A)')trim(file_line)
            enddo
        close(04)

    close(03)


end program generate_suitable_SCKPTS_from_pckpts
