!    Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!
!    This file is part of BandUP: Band Unfolding code for Plane-wave based calculations.
!
!    BandUP is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BandUP is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BandUP.  If not, see <http://www.gnu.org/licenses/>.
program get_wavefunc_SCKPTS_needed_for_EBS
use io_routines
use math
implicit none
character(len=127) :: input_file_prim_cell, input_file_supercell, input_file_pc_kpts, outpt_file_SC_kpts, outpt_file_SC_kpts_not_reduced
real(kind=dp), dimension(:,:), allocatable :: k_starts, k_ends
real(kind=dp), dimension(1:3,1:3) :: dir_latt_pc, dir_latt_SC, b_matrix_pc, B_matrix_SC
real(kind=dp), dimension(1:3) :: point
integer, dimension(:), allocatable :: n_pckpts_dirs
integer :: n_irr_unfolding_SCKPTS,ikpt,i,ios,nkpts,start_k,end_k,n_selec_pcbz_dirs,idir
type(vec3d), dimension(:), allocatable :: list, irr_unfolding_SCKPTS
character(len=254) :: file_line
! Variables for the symmetry analysis
type(symmetry_operation), dimension(:), allocatable :: symops_pc, symops_SC
type(vec3d), dimension(:), allocatable :: considered_kpts_list
character(len=11) :: int_symb_pcbz, int_symb_SCBZ
character(len=10) :: schoenflies_notation_pc, schoenflies_notation_SC
real(kind=dp), parameter :: symprec = 1d-6
integer :: num_symm_op_pcbz, num_symm_op_SCBZ
logical :: get_all_kpts_needed_for_EBS_averaging
!!!**********************************************************************************************

    get_all_kpts_needed_for_EBS_averaging = .TRUE.
    input_file_prim_cell = 'prim_cell_lattice.in'
    input_file_supercell = 'supercell_lattice.in'
    input_file_pc_kpts = 'KPOINTS_prim_cell.in'
    outpt_file_SC_kpts = 'KPOINTS_supercell.out'
    outpt_file_SC_kpts_not_reduced = 'KPOINTS_supercell_not_reduced.out'

    call print_welcome_messages(package_version)
    write(*,'(A)')'               Pre-processing utility "get_SCKPTS_pre_BandUP"'
    write(*,'(A)')'   >>> Getting the SC-KPTS you will need for your plane-wave calculation <<<'
    write(*,*)
    call read_unit_cell(input_file=input_file_prim_cell,latt=dir_latt_pc(:,:))
    call get_rec_latt(latt=dir_latt_pc,rec_latt=b_matrix_pc)
    ! Symmetry analysis for the pcbz
    call get_symm(nsym=num_symm_op_pcbz, symops=symops_pc, & 
                  schoenflies=schoenflies_notation_pc, international_symb=int_symb_pcbz, &
                  symprec=symprec, lattice=b_matrix_pc(:,:))
    write(*,'(A,I0,A)')'Found ',num_symm_op_pcbz,' symmetry operations for the pcbz.' 
    write(*,'(5A)')'    * Space group of the pcbz: ',trim(adjustl(schoenflies_notation_pc)),'(',trim(adjustl(int_symb_pcbz)),').'
    ! End of symmetry analysis for the pcbz
    call read_unit_cell(input_file=input_file_supercell,latt=dir_latt_SC(:,:))
    call get_rec_latt(latt=dir_latt_SC,rec_latt=B_matrix_SC)
    ! Symmetry analysis for the SCBZ
    call get_symm(nsym=num_symm_op_SCBZ, symops=symops_SC, & 
                  schoenflies=schoenflies_notation_SC, international_symb=int_symb_SCBZ, &
                  symprec=symprec, lattice=B_matrix_SC(:,:))
    write(*,'(A,I0,A)')'Found ',num_symm_op_SCBZ,' symmetry operations for the SCBZ.' 
    write(*,'(5A)')'    * Space group of the SCBZ: ',trim(adjustl(schoenflies_notation_SC)),'(',trim(adjustl(int_symb_SCBZ)),').'
    ! End of symmetry analysis for the SCBZ
    !! Reading selected pckpts from the input file
    call read_pckpts_selected_by_user(k_starts=k_starts, k_ends=k_ends, ndirs=n_selec_pcbz_dirs, n_kpts_dirs=n_pckpts_dirs, &
                                      input_file=input_file_pc_kpts,b_matrix_pc=b_matrix_pc(:,:))
    allocate(list(1:sum(n_pckpts_dirs(:))))
    start_k=1
    do idir=1,n_selec_pcbz_dirs
        end_k = start_k + n_pckpts_dirs(idir) - 1
        list(start_k:end_k) = kpts_line(k_starts(idir,:),k_ends(idir,:),n_pckpts_dirs(idir))
        start_k = end_k + 1
    enddo


    if(get_all_kpts_needed_for_EBS_averaging  .and. (.not. int_symb_SCBZ == int_symb_pcbz))then
        !! Obtaining all pckpts that are equivalent to the selected ones by symm. ops. of the pcbz
        !! This is needed when the SC and the pc do not possess the same symmetry, in which cases
        !! they will not be equivalent w.r.t. symm. ops. of the  SCBZ
        call get_star(list_of_all_generated_points=considered_kpts_list, points=list, lattice=b_matrix_pc(:,:))
    else
        allocate(considered_kpts_list(1:size(list)))
        considered_kpts_list = list
    endif
    deallocate(list)

    !! Getting the smallest number of SCKPTS possible for the calculation of the
    !! EBS for the selected directions in the pcbz
    call get_irr_kpts(n_irr_kpts=n_irr_unfolding_SCKPTS,irr_kpts_list=irr_unfolding_SCKPTS, &
                      kpts_list=considered_kpts_list,rec_latt=B_matrix_SC,reduce_to_bz=.TRUE.)
     
    write(*,'(A,I0,A)')'A total of ',n_irr_unfolding_SCKPTS," irreducible SCKPTS are needed to obtain an EBS along the selected direction(s) of the reference pcbz."
    write(*,"(3A)")'>>> The SCBZ-Kpoints you will need to run your plane-wave calculation have been stored in the file "',trim(adjustl(outpt_file_SC_kpts)),'".'
    if(get_all_kpts_needed_for_EBS_averaging .and. (.not. int_symb_SCBZ == int_symb_pcbz))then
        write(*,"(A)")                "===================================================================================================="
        write(*,"(A)")                "NOTICE:"
        write(*,'(A,/,A,/,A,/,A,/,A)')"       We have considered more pcbz directions than what you asked for. The reason is that, since", &
                                      "       the symmetry groups of the SCBZ and the pcbz are not equal, the pc-kpts in such directions", &
                                      "       will not necessarily be equivalent by symmetry operations of the SCBZ. This should not be", &
                                      "       forgotten when calculating the EBS.", &
                                      "       Don't worry, though. Only irreducible complementary directions (if any) have been considered."
        write(*,"(A)")                "===================================================================================================="
    endif
    !!! Writing results to the output file
    open(unit=03, file=outpt_file_SC_kpts)
        open(unit=04, file=input_file_pc_kpts)
            read(04,'(A)')file_line
            write(03,'(A)')trim(file_line)//' (this is exactly the header of the input kpts file)'
        close(04)
        write(03,'(I0)')n_irr_unfolding_SCKPTS
        write(03,'(A)')'Reciprocal (fractional) coords. w.r.t. the SCRL vectors:' 
        do ikpt=1, n_irr_unfolding_SCKPTS
            point(:) = irr_unfolding_SCKPTS(ikpt)%coord(:)
            write(03,'(3(2X,f10.7),2X,I1)')(point(i), i=1,3),1
        enddo
        write(03,'(A)')''
        write(03,'(A)')''
        write(03,'(A)')'! The above SCKPTS (and/or some other SCKPTS related to them by symm. ops. of the SCBZ)' 
        write(03,'(A)')'! unfold onto the pckpts listed below (selected by you) (and/or some other pckpts related to them by symm. ops. of the pcbz): '
        write(03,'(A)')'! (Fractional coords. w.r.t. the pcrl vectors) '
        open(unit=04, file=input_file_pc_kpts)
            ios=0
            do while(ios==0)
                read(04,'(A)',iostat=ios)file_line
                if(ios==0) write(03,'(A)')trim(file_line)
            enddo
        close(04)
    close(03)

end program get_wavefunc_SCKPTS_needed_for_EBS
