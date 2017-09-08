!! Copyright (C) 2013-2016 Paulo V. C. Medeiros
!!
!! This file is part of BandUP:
!! Band Unfolding code for Plane-wave based calculations.
!!
!! BandUP is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BandUP is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BandUP.  If not, see <http://www.gnu.org/licenses/>.
program get_wavefunc_SCKPTS_needed_for_EBS
use constants_and_types
use cla_wrappers
use general_io
use io_routines
use read_vasp_files
use write_vasp_files
use symmetry, only: get_prim_cell, get_irr_sc_kpts, analise_symm_pc_sc, &
                    get_pcbz_dirs_2b_used_for_ebs
use band_unfolding
implicit none
type(vec3d), dimension(:), allocatable :: irr_unfolding_SCKPTS, &
                                          irr_unfolding_SCKPTS_frac_coords, &
                                          aux_irr_unfolding_SCKPTS, &
                                          aux_irr_unfolding_SCKPTS_frac_coords
type(crystal_3D) :: crystal_pc, crystal_SC
real(kind=dp), dimension(:,:), allocatable :: k_starts, k_ends
real(kind=dp), dimension(1:3) :: point
integer :: n_irr_unfolding_SCKPTS,ikpt,i,ios,n_selec_pcbz_dirs,idir, &
           n_decimals, total_size_float_format, icoord
character(len=str_len) :: file_line, str_n_decimals, &
                          str_total_size_float_format, float_format
! Variables for the symmetry analysis
type(irr_bz_directions), dimension(:), allocatable :: &
    all_dirs_used_for_EBS_along_pcbz_dir
type(vec3d), dimension(:), allocatable :: considered_kpts_list
type(selected_pcbz_directions) :: pckpts_to_be_checked
! Geometric Unfolding Relations
type(geom_unfolding_relations_for_each_SCKPT) :: GUR 
integer, dimension(:), allocatable :: n_pckpts_dirs, &
                                      n_dirs_for_EBS_along_pcbz_dir
integer :: i_req_dir, ikpt2, i_irr_kpt, aux_n_irr_unfolding_SCKPTS
!!*****************************************************************************

    call print_welcome_messages(package_version)
    write(*,'(A)')&
        '               Pre-processing utility "get_SCKPTS_pre_BandUP"'
    write(*,'(A)')&
        "   >>> Getting the SC-KPTS you'll need for your plane-wave calc. <<<"
    write(*,*)
    call get_commline_args(args, running_from_main_code=.FALSE.)
    call get_crystal_from_file(&
             crystal_pc,input_file=args%input_file_prim_cell, &
             stop_if_file_not_found=.TRUE. &
         )
    call get_prim_cell(crystal_pc, symprec=default_symprec)
    call get_crystal_from_file(&
             crystal_SC,input_file=args%input_file_supercell, &
             stop_if_file_not_found=.TRUE. &
         )
    call get_prim_cell(crystal_SC, symprec=default_symprec)
    call write_attempted_pc_assoc_with_input_unit_cell_and_SC(&
             crystal_pc, crystal_SC &
         )
    call verify_commens(crystal_pc, crystal_SC, args)
    call analise_symm_pc_SC(crystal_pc, crystal_SC)

    call read_pckpts_selected_by_user(&
             k_starts=k_starts, k_ends=k_ends, &
             ndirs=n_selec_pcbz_dirs, n_kpts_dirs=n_pckpts_dirs, &
             input_file=args%input_file_pc_kpts, &
             b_matrix_pc=crystal_pc%rec_latt_vecs &
         )
    call get_pcbz_dirs_2b_used_for_EBS(&
             all_dirs_used_for_EBS_along_pcbz_dir, crystal_pc, crystal_SC, &
             k_starts, k_ends, args &
         )
    call print_symm_analysis_for_selected_pcbz_dirs(&
             all_dirs_used_for_EBS_along_pcbz_dir &
         )
    call define_pckpts_to_be_checked(&
             pckpts_to_be_checked, all_dirs_used_for_EBS_along_pcbz_dir, &
             n_pckpts_dirs(:) &
         )

    ! Getting a list of all considered pc-kpts, including the ones chosen 
    ! by the user and those obtained by symm. for the complementary directions
    call define_pckpts_to_be_checked(&
             pckpts_to_be_checked, all_dirs_used_for_EBS_along_pcbz_dir, &
             n_pckpts_dirs(:) &
         )
    allocate(n_dirs_for_EBS_along_pcbz_dir(1:n_selec_pcbz_dirs))
    do idir=1,n_selec_pcbz_dirs
        n_dirs_for_EBS_along_pcbz_dir(idir) = &
            size(all_dirs_used_for_EBS_along_pcbz_dir(idir)%irr_dir(:))
    enddo
    allocate(&
        considered_kpts_list(1:dot_product(n_pckpts_dirs(:), &
                             n_dirs_for_EBS_along_pcbz_dir(:))) &
    )
    ikpt = 0
    do idir=1,n_selec_pcbz_dirs
        do i_req_dir=1,n_dirs_for_EBS_along_pcbz_dir(idir)
            do ikpt2=1, n_pckpts_dirs(idir)
                ikpt = ikpt + 1
                considered_kpts_list(ikpt)%coord(:) = &
                    pckpts_to_be_checked%selec_pcbz_dir(idir)%&
                        needed_dir(i_req_dir)%pckpt(ikpt2)%coords(:)
            enddo
        enddo
    enddo
    !! Getting the smallest possible number of SCBZ KPTS for the calculation
    !! of the EBS along the selected direction(s) of the pcbz
    write(*,'(A)')'Getting the needed SCBZ-Kpoints...'
    call get_irr_SC_kpts(&
             n_irr_kpts=aux_n_irr_unfolding_SCKPTS, &
             irr_kpts_list=aux_irr_unfolding_SCKPTS, &
             irr_kpts_list_frac_coords=aux_irr_unfolding_SCKPTS_frac_coords, &
             kpts_list=considered_kpts_list, &
             crystal=crystal_SC, args=args, reduce_to_bz=.TRUE. &
         )
    call get_geom_unfolding_relations(&
             GUR, aux_irr_unfolding_SCKPTS, &
             pckpts_to_be_checked, crystal_pc, crystal_SC &
         )
    
    n_irr_unfolding_SCKPTS = count(GUR%SCKPT_used_for_unfolding(:))
    allocate(irr_unfolding_SCKPTS(1:n_irr_unfolding_SCKPTS), &
             irr_unfolding_SCKPTS_frac_coords(1:n_irr_unfolding_SCKPTS))
    i_irr_kpt = 0
    do ikpt=1,size(GUR%SCKPT_used_for_unfolding(:))
        if(GUR%SCKPT_used_for_unfolding(ikpt))then
            i_irr_kpt = i_irr_kpt + 1
            irr_unfolding_SCKPTS(i_irr_kpt) = aux_irr_unfolding_SCKPTS(ikpt)
            irr_unfolding_SCKPTS_frac_coords(i_irr_kpt) = &
                aux_irr_unfolding_SCKPTS_frac_coords(ikpt)
        endif
    enddo

    call get_geom_unfolding_relations(&
             GUR, irr_unfolding_SCKPTS, &
             pckpts_to_be_checked, crystal_pc, crystal_SC &
         )
    call print_message_success_determining_GUR(&
             GUR, stop_if_GUR_fails, is_main_code=.FALSE. &
         ) 
    if((GUR%n_pckpts /= GUR%n_folding_pckpts) .and. stop_if_GUR_fails) stop

    n_irr_unfolding_SCKPTS = count(GUR%SCKPT_used_for_unfolding(:))
    if(count(.not. GUR%SCKPT_used_for_unfolding) > 0)then
        write(*,'(A,I0,A)')&
            'WARNING: There seems to be more SC-KPTS than the necessarry (', &
            count(.not. GUR%SCKPT_used_for_unfolding),' too many).'
    endif

    if(print_GUR_pre_unfolding_utility)then 
        call print_geom_unfolding_relations(&
                 GUR, irr_unfolding_SCKPTS, crystal_pc, crystal_SC &
             )
    endif

    if(args%no_symm_avg)then
        write(*,'(A,I0,A)')&
            'A total of ',n_irr_unfolding_SCKPTS, &
            " SCBZ-Kpoints will be needed to perform the unfolding along the &
            selected direction(s) of the reference pcbz."
    else
        write(*,'(A,I0,A)')&
            'A total of ',n_irr_unfolding_SCKPTS, &
            " SCBZ-Kpoints will be needed to obtain a symmetry-averaged EBS &
            along the selected direction(s) of the reference pcbz."
    endif
    write(*,"(3A)")&
        '>>> The SCBZ-Kpts you will need to run your plane-wave calc have &
        been stored in the file "',trim(adjustl(args%out_file_SC_kpts)),'".'

    do idir=1, size(all_dirs_used_for_EBS_along_pcbz_dir(:))
        if(all_dirs_used_for_EBS_along_pcbz_dir(idir)%ncompl_dirs > 0)then
            write(*,"(10(A,/),A)")&
"===========================================================================",&
"NOTICE:                                                                    ",&
"    We have considered more pcbz directions than you asked for.            ",&
"    This is because the SC and pc belong to different symmetry groups, and,",&
"    therefore, some pcbz kpts that are equivalent by symmetry ops of the pc",&
"    might not be equivalent by symmetry ops of the SC.                     ",&
"    Don't worry, though: We only kept irreducible complementary directions.",&
"    >> If you don't want BandUP to consider these extra pcbz direcs., then ",&
'       run both the "kpts-sc-get" AND "unfold" tasks with the option       ',&
"                                '-no_symm_avg'.                            ",&
"============================================================================="
            exit
        endif
    enddo
    !!! Writing results to the output file
    n_decimals = nint(abs((log10(default_symprec))))
    ! The numbers will be between -1 and 1
    total_size_float_format = 3 + n_decimals
    write(str_n_decimals,'(I0)') n_decimals
    write(str_total_size_float_format,'(I0)') total_size_float_format
    float_format = 'f' // trim(adjustl(str_total_size_float_format)) // &
                   '.' // trim(adjustl(str_n_decimals))
    open(unit=03, file=args%out_file_SC_kpts)
        open(unit=04, file=args%input_file_pc_kpts)
            read(04,'(A)')file_line
            write(03,'(A)')trim(adjustl(file_line))//' (this is &
                           exactly the header of the input kpts &
                           file) ' // file_header_BandUP()
        close(04)
        write(03,'(I0)')n_irr_unfolding_SCKPTS
        write(03,'(A)')&
            'Reciprocal (fractional) coords. w.r.t. the SCRL vectors:' 
        do ikpt=1, n_irr_unfolding_SCKPTS
            point(:) = irr_unfolding_SCKPTS_frac_coords(ikpt)%coord(:)
            do icoord=1,3
                point(icoord) = &
                    real(nint(point(icoord)/default_symprec), kind=dp) * &
                    default_symprec
            enddo
            if(trim(adjustl(args%pw_code))=='abinit')then
                write(03,'(3(2X,'//trim(adjustl(float_format))//'),2X,I1)')&
                    (point(i), i=1,3)
            else
                write(03,'(3(2X,'//trim(adjustl(float_format))//'),2X,I1)')&
                    (point(i), i=1,3), 1
            endif
        enddo
        write(03,'(A)')''
        write(03,'(A)')''
        write(03,'(A)')'! The above SCKPTS (and/or some other SCKPTS &
                        related to them by symm. ops. of the SCBZ)' 
        write(03,'(A)')'! unfold onto the pckpts listed below &
                        (selected by you) (and/or some other pckpts &
                        related to them by symm. ops. of the pcbz): '
        write(03,'(A)')'! (Fractional coords. w.r.t. the pcrl vectors) '
        open(unit=04, file=args%input_file_pc_kpts)
            ios=0
            do while(ios==0)
                read(04,'(A)',iostat=ios)file_line
                if(ios==0) write(03,'(A)')trim(file_line)
            enddo
        close(04)
    close(03)


end program get_wavefunc_SCKPTS_needed_for_EBS
