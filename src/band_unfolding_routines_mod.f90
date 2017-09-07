!! Copyright (C) 2013-2017 Paulo V. C. Medeiros
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

!==============================================================================
! MODULE: band_unfolding 
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION:
!> Contains all routines employed to perform unfolding in BandUP.
!==============================================================================

module band_unfolding
!$ use omp_lib
use constants_and_types
use cla_wrappers
use lists_and_seqs, only: append, list_index, kpts_line
use math, only: integral_delta_x_minus_x0, delta, trace_ab, versor, cross 
use symmetry, only: get_symm, get_star, pt_eqv_by_point_group_symop
use crystals, only: vec_in_latt, check_if_pc_and_SC_are_commensurate
use time, only: time_now
use general_io
use io_routines
implicit none
PRIVATE
PUBLIC :: get_geom_unfolding_relations, define_pckpts_to_be_checked, &
          select_coeffs_to_calc_spectral_weights, get_delta_Ns_for_output, &
          update_GUR_indices, perform_unfolding, verify_commens, &
          allocate_UnfoldedQuantities, allocate_UnfoldedQuantitiesForOutput

CONTAINS


subroutine allocate_UnfoldedQuantities(delta_N,pckpts_to_be_checked)
implicit none
type(UnfoldedQuantities), intent(out) :: delta_N
!! Geometric Unfolding Relations
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked 
integer :: n_selec_pcbz_dir,i_selc_pcbz_dir,n_needed_dirs,i_needed_dirs,nkpts

    n_selec_pcbz_dir = size(pckpts_to_be_checked%selec_pcbz_dir(:))
    allocate(delta_N%selec_pcbz_dir(1:n_selec_pcbz_dir))
    do i_selc_pcbz_dir=1,n_selec_pcbz_dir
        n_needed_dirs = size(&
                            pckpts_to_be_checked%&
                                selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(:) &
                            )
        allocate(&
            delta_N%selec_pcbz_dir(i_selc_pcbz_dir)%&
                needed_dir(1:n_needed_dirs) &
        )
        do i_needed_dirs=1,n_needed_dirs
            nkpts = size(&
                        pckpts_to_be_checked%selec_pcbz_dir(i_selc_pcbz_dir)%&
                            needed_dir(i_needed_dirs)%pckpt(:) &
                    )
            allocate(&
                delta_N%selec_pcbz_dir(i_selc_pcbz_dir)%&
                    needed_dir(i_needed_dirs)%pckpt(1:nkpts) &
            )
        enddo
    enddo

end subroutine allocate_UnfoldedQuantities

subroutine allocate_UnfoldedQuantitiesForOutput(delta_N,pckpts_to_be_checked)
implicit none
type(UnfoldedQuantitiesForOutput), intent(out) :: delta_N
!! Geometric Unfolding Relations
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked 
integer :: n_selec_pcbz_dir,i_selc_pcbz_dir,nkpts

    n_selec_pcbz_dir = size(pckpts_to_be_checked%selec_pcbz_dir(:))
    allocate(delta_N%pcbz_dir(1:n_selec_pcbz_dir))
    do i_selc_pcbz_dir=1,n_selec_pcbz_dir
        nkpts = size(&
                    pckpts_to_be_checked%selec_pcbz_dir(i_selc_pcbz_dir)%&
                        needed_dir(1)%pckpt(:) &
                )
        allocate(delta_N%pcbz_dir(i_selc_pcbz_dir)%pckpt(1:nkpts))
    enddo

end subroutine allocate_UnfoldedQuantitiesForOutput


subroutine update_GUR_indices(&
               GUR, i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt&
           )
implicit none
type(geom_unfolding_relations_for_each_SCKPT), intent(inout) :: GUR
integer, intent(in) :: i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt

    GUR%current_index%i_SCKPT = i_SCKPT
    GUR%current_index%i_selec_pcbz_dir = i_selec_pcbz_dir
    GUR%current_index%i_needed_dirs = i_needed_dirs
    GUR%current_index%ipc_kpt = ipc_kpt

end subroutine update_GUR_indices


subroutine verify_commens(crystal_pc, crystal_SC, args)
implicit none
type(crystal_3D), intent(in) :: crystal_pc, crystal_SC
type(comm_line_args), intent(in) :: args
real(kind=dp), dimension(1:3, 1:3) :: matrix_M
logical :: are_commens

    call check_if_pc_and_SC_are_commensurate(&
             are_commens, matrix_M, crystal_pc, crystal_SC, &
             tol=default_tol_for_int_commens_test &
         )
    call print_message_commens_test(&
             commensurate=are_commens, M=matrix_M,&
             stop_if_not_commens=args%stop_if_not_commensurate &
         ) 
    if(args%stop_if_not_commensurate .and. .not. are_commens) stop

end subroutine verify_commens


subroutine get_GUR_not_public(GUR,list_of_SCKPTS, pckpts_to_be_checked, &
                              input_crystal_pc, input_crystal_SC, &
                              vec_in_latt_tol_for_vec_eq, verbose)
!! Copyright (C) 2013-2017 Paulo V. C. Medeiros
implicit none
!! GUR = Geometric Unfolding Relations
type(geom_unfolding_relations_for_each_SCKPT), intent(out) :: GUR 
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(vec3d), dimension(:), intent(in) :: list_of_SCKPTS
type(crystal_3D), intent(in) :: input_crystal_pc, input_crystal_SC
real(kind=dp), intent(in), optional :: vec_in_latt_tol_for_vec_eq
logical, intent(in), optional :: verbose
integer :: nkpts, n_selec_pcbz_dirs, i_SCKPT,ipc_kpt,ieqv_SCKPT,isym, &
           i_selec_pcbz_dir,i_needed_dirs,alloc_stat, i
integer, dimension(:), allocatable :: n_dirs_for_EBS_along_pcbz_dir, &
                                      n_pckpts_dirs
logical, dimension(:,:,:), allocatable :: pc_kpt_already_folded
real(kind=dp) :: vec_in_latt_tolerance_for_vec_eq
real(kind=dp), dimension(1:3) :: pc_kpt, current_SCKPT, &
                                 SCKPT_eqv_to_current_SCKPT, &
                                 trial_folding_G, origin_for_spin_proj
real(kind=dp), dimension(1:3,1:3) :: B_matrix_SC
type(star), dimension(:), allocatable :: SKPTS_eqv_to_SKPT
logical :: print_stuff
type(crystal_3D) :: crystal_pc, crystal_SC

    crystal_SC = input_crystal_SC
    crystal_pc = input_crystal_pc
    print_stuff = .FALSE.
    if(present(verbose))then
        print_stuff = verbose
    endif

    vec_in_latt_tolerance_for_vec_eq = default_tol_for_vec_equality 
    if(present(vec_in_latt_tol_for_vec_eq))then
        vec_in_latt_tolerance_for_vec_eq = abs(vec_in_latt_tol_for_vec_eq)
    endif

    if(print_stuff)then
        write(*,"(A)")"Verifying geometric unfolding relations &
                       between pcbz and SCBZ wave-vectors..."
    endif

    B_matrix_SC = crystal_SC%rec_latt_vecs
    GUR%B_matrix_SC = B_matrix_SC
    GUR%b_matrix_pc = crystal_pc%rec_latt_vecs
    if(args%origin_for_spin_proj_passed_in_rec)then
        args%origin_for_spin_proj_cartesian(:) = 0.0_dp
        do i=1,3
            args%origin_for_spin_proj_cartesian = &
                args%origin_for_spin_proj_cartesian + &
                args%origin_for_spin_proj_rec(i) * B_matrix_SC(i,:)
        enddo
    endif

    if(args%no_symm_sckpts)then
        if(print_stuff)then
            write(*,"(A)")'    * Running with "-no_symm_sckpts".'
            write(*,"(A)")'    * The SC symmetry will NOT be used'
            write(*,"(A)")'      to reduce the number of SC K-points needed.'
        endif
        crystal_SC%nsym = 1
        deallocate(crystal_SC%symops, stat=alloc_stat)
        allocate(crystal_SC%symops(1:1))
        crystal_SC%symops(1)%translation_fractional_coords(:) = 0.0_dp
        crystal_SC%symops(1)%translation_cartesian_coords(:) = 0.0_dp
        crystal_SC%symops(1)%rotation_fractional_coords = identity_3D
        crystal_SC%symops(1)%rotation_cartesian_coords = identity_3D
    else
        ! This fails if use_pc_to_get_symm=.TRUE.
        call get_symm(&
                 crystal=crystal_SC, use_pc_to_get_symm=.FALSE., &
                 symprec=default_symprec &
             ) 
    endif
    call get_star(star_of_pt=SKPTS_eqv_to_SKPT, points=list_of_SCKPTS, &
                  crystal=crystal_SC, &
                  tol_for_vec_equality=default_tol_for_vec_equality, &
                  symprec=default_symprec, reduce_to_bz=.TRUE.)
    !! Allocating and initializing table
    nkpts = size(list_of_SCKPTS)
    allocate(GUR%list_of_SCKPTS(1:nkpts))
    GUR%list_of_SCKPTS = list_of_SCKPTS
    n_selec_pcbz_dirs = size(pckpts_to_be_checked%selec_pcbz_dir(:))
    allocate(n_dirs_for_EBS_along_pcbz_dir(1:n_selec_pcbz_dirs))
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir) = &
            size(&
                pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%&
                    needed_dir(:) &
            )
    enddo
    allocate(n_pckpts_dirs(1:n_selec_pcbz_dirs))
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        n_pckpts_dirs(i_selec_pcbz_dir) = &
            size(&
                pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%&
                    needed_dir(1)%pckpt(:) &
            ) 
    enddo
    GUR%n_folding_pckpts = 0
    GUR%n_pckpts = 0
    deallocate(GUR%SCKPT, stat=alloc_stat)
    allocate(GUR%SCKPT(1:nkpts))
    deallocate(GUR%SCKPT_used_for_unfolding, stat=alloc_stat)
    allocate(GUR%SCKPT_used_for_unfolding(1:nkpts))
    do i_SCKPT=1,nkpts
        deallocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir, stat=alloc_stat)
        allocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(1:n_selec_pcbz_dirs))
        GUR%SCKPT_used_for_unfolding(i_SCKPT) = .FALSE.
        do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
            deallocate(&
                GUR%SCKPT(i_SCKPT)%&
                    selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir, &
                stat=alloc_stat &
            )
            allocate(&
                GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)% &
                    needed_dir(&
                        1:n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir) &
                              ) &
            )
            do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)
                deallocate(&
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)% &
                        needed_dir(i_needed_dirs)%pckpt, &
                    stat=alloc_stat &
                )
                allocate(&
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)% &
                        needed_dir(i_needed_dirs)%&
                        pckpt(1:n_pckpts_dirs(i_selec_pcbz_dir)) &
                )
                do ipc_kpt=1,n_pckpts_dirs(i_selec_pcbz_dir)
                    GUR%n_pckpts = GUR%n_pckpts + 1
                    pc_kpt(:) = pckpts_to_be_checked%&
                                    selec_pcbz_dir(i_selec_pcbz_dir)% &
                                    needed_dir(i_needed_dirs)%&
                                    pckpt(ipc_kpt)% &
                                    coords(:)
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%&
                        needed_dir(i_needed_dirs)%&
                        pckpt(ipc_kpt)%coords(:) &
                        = pc_kpt(:)
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%&
                        needed_dir(i_needed_dirs)%&
                        pckpt(ipc_kpt)%Scoords(:) &
                        = pc_kpt(:) 
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%&
                        needed_dir(i_needed_dirs)%&
                        pckpt(ipc_kpt)%folds &
                        = .FALSE.
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%&
                        needed_dir(i_needed_dirs)%&
                        pckpt(ipc_kpt)%Sfolding_vec(:) &
                        = 0.0_dp 
                enddo
            enddo
        enddo
    enddo
    GUR%n_pckpts = (GUR%n_pckpts)/nkpts
    deallocate(pc_kpt_already_folded,stat=alloc_stat)
    allocate(&
        pc_kpt_already_folded(&
            1:n_selec_pcbz_dirs, &
            1:maxval(n_dirs_for_EBS_along_pcbz_dir(:)), &
            1:maxval(n_pckpts_dirs(:)) &
        ) &
    )
    pc_kpt_already_folded = .FALSE.
    !! Obtaining geometric unfolding relations
    do i_SCKPT=1,nkpts
        ! Ignoring the 1st "args%n_sckpts_to_skip" SCKPTS
        ! This might be used, e.g., together with hybrid functional calcs 
        ! when using VASP.
        ! Since everything is initialized to "not folding", this is safe.
        if(i_SCKPT <= args%n_sckpts_to_skip) cycle 
        do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
            do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)
                do ipc_kpt=1, n_pckpts_dirs(i_selec_pcbz_dir)
                    if(pc_kpt_already_folded(&
                           i_selec_pcbz_dir,i_needed_dirs,ipc_kpt&
                       ))then
                        ! It has already folded before.
                        ! Since a pckpt can only fold into one SCKPT,
                        ! let's move on to the next pckpt.
                        cycle 
                    endif
                    pc_kpt(:) = pckpts_to_be_checked%&
                                    selec_pcbz_dir(i_selec_pcbz_dir)% &
                                    needed_dir(i_needed_dirs)% &
                                    pckpt(ipc_kpt)%coords(:)
                    current_SCKPT = list_of_SCKPTS(i_SCKPT)%coord(:)
                    do ieqv_SCKPT=1,SKPTS_eqv_to_SKPT(i_SCKPT) % neqv
                        SCKPT_eqv_to_current_SCKPT(:) = &
                            SKPTS_eqv_to_SKPT(i_SCKPT)%&
                                eqv_pt(ieqv_SCKPT) % coord(:)
                        trial_folding_G(:) = pc_kpt(:) - &
                                             SCKPT_eqv_to_current_SCKPT(:)
                        if(vec_in_latt(&
                               vec=trial_folding_G, latt=B_matrix_SC, &
                               tolerance=vec_in_latt_tolerance_for_vec_eq &
                           ))then
                            GUR%SCKPT_used_for_unfolding(i_SCKPT) = .TRUE.
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selec_pcbz_dir)% &
                                needed_dir(i_needed_dirs)%&
                                pckpt(ipc_kpt)%folds &
                                = .TRUE.
                            GUR%n_folding_pckpts = GUR%n_folding_pckpts + 1
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selec_pcbz_dir)% &
                                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)% &
                                coords_actual_unfolding_K &
                                = SCKPT_eqv_to_current_SCKPT(:)
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selec_pcbz_dir)% &
                                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)% &
                                coords_SCKPT_used_for_coeffs &
                                = current_SCKPT(:)
                            isym = SKPTS_eqv_to_SKPT(i_SCKPT)%&
                                       eqv_pt(ieqv_SCKPT)%symop
                            pc_kpt = pt_eqv_by_point_group_symop(&
                                         point=pc_kpt, &
                                         symops=crystal_SC%symops, &
                                         isym=isym, &
                                         fractional_coords=.FALSE., &
                                         invert_symop=.TRUE. &
                                     )
                            origin_for_spin_proj = &
                                pt_eqv_by_point_group_symop(&
                                    point=args%origin_for_spin_proj_cartesian,&
                                    symops=crystal_SC%symops, isym=isym, &
                                    fractional_coords=.FALSE., &
                                    invert_symop=.TRUE. &
                            )
                            ! Message from Paulo:
                            ! The prefix "S" means "symmetrized". 
                            ! This is a little trick that allows me to use the 
                            ! coefficients of a SC wavefunction
                            ! psi(K',n) to calculate the spectral weights 
                            ! associated with a SC wf
                            ! psi(K,n), where K' = SK and S is a symme op
                            ! of the crystal's point group.
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selec_pcbz_dir)%&
                                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                                Scoords(:) &
                                = pc_kpt(:)
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selec_pcbz_dir)%&
                                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                                Sorigin_for_spin_proj(:) &
                                = origin_for_spin_proj(:)
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selec_pcbz_dir)%&
                                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)% &
                                Sfolding_vec(:) &
                                = pc_kpt(:) - current_SCKPT(:)
                            pc_kpt_already_folded(&
                                i_selec_pcbz_dir, &
                                i_needed_dirs,ipc_kpt &
                                ) &
                                = .TRUE.
                            ! It has already folded:
                            ! No need to keep looking for this pckpt.
                            exit
                        else 
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selec_pcbz_dir)%&
                                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                                folds &
                                = .FALSE.
                        endif
                    ! Loop over all the SCKPTS eqv. to the one on the WF file 
                    ! by SC symmops
                    enddo 
                ! Loop over the pckpts along the dirs. of the pcbz that are 
                ! eqv. to the selec. ones by symm. ops. of the pcbz and 
                ! NOT by symm. ops. of the SCBZ
                enddo 
            ! Loop over all dirs. of the pcbz that are eqv. to the selec. ones 
            ! by symm. ops. of the pcbz and NOT by symm. ops. of the SCBZ
            enddo 
        enddo ! Loop over the selected pcbz directions
    enddo ! Loop over the SCKPTS found on the wavefunction file

    if(print_stuff)then
        write(*,"(A,/)")"    * Done."
    endif


end subroutine get_GUR_not_public


subroutine get_geom_unfolding_relations(&
               GUR,list_of_SCKPTS, pckpts_to_be_checked, &
               input_crystal_pc, input_crystal_SC, verbose &
           )
!! Copyright (C) 2013-2017 Paulo V. C. Medeiros
!! This routine is a wrapper for the internal routine get_GUR_not_public, which
!! determines the GUR between PC and SC. If get_GUR_not_public fails in the
!! first attempt, than the present routine calls it again using a sligtly
!! increased value for the tolerance for vector equality to be used in the
!! routine that checks if vectors belong to a given Bravais lattice. 
!! The present routine keeps doing that until the GUR are successfully 
!! determined or until a maximum value for the tolerance parameter is reached.
implicit none
!! Geometric Unfolding Relations
type(geom_unfolding_relations_for_each_SCKPT), intent(out) :: GUR 
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(vec3d), dimension(:), intent(in) :: list_of_SCKPTS
type(crystal_3D), intent(in) :: input_crystal_pc, input_crystal_SC
logical, intent(in), optional :: verbose
! Local vars
integer :: n_attempts
real(kind=dp) :: vec_in_latt_tol_for_vec_eq
logical :: GUR_successfully_determined


    GUR_successfully_determined = .FALSE.
    vec_in_latt_tol_for_vec_eq = abs(default_tol_for_vec_equality)
    n_attempts = 0
    do while((.not. GUR_successfully_determined) .and. &
              vec_in_latt_tol_for_vec_eq <= abs(max_tol_for_vec_equality))
        call get_GUR_not_public(GUR,list_of_SCKPTS, pckpts_to_be_checked, &
                                input_crystal_pc, input_crystal_SC, &
                                vec_in_latt_tol_for_vec_eq, verbose)
        GUR_successfully_determined = GUR%n_pckpts == GUR%n_folding_pckpts
        vec_in_latt_tol_for_vec_eq = vec_in_latt_tol_for_vec_eq + &
                                     0.05_dp*abs(default_tol_for_vec_equality)
        n_attempts = n_attempts + 1
        if((.not. GUR_successfully_determined) .and. &
           vec_in_latt_tol_for_vec_eq <= abs(max_tol_for_vec_equality))then
            write(*,'(A,I0,4(A,/))') &
'WARNING (get_geom_unfolding_relations): Failed to determine GURs (attempt #',&
                                         n_attempts, ').', &
'                                        Increasing, by &
                                         0.05*default_tol_for_vec_equality,', &
'                                        the tol on the check if vectors', &
'                                        belong to lattice and trying again.'
        endif
    enddo
    vec_in_latt_tol_for_vec_eq = vec_in_latt_tol_for_vec_eq - &
                                 0.1_dp * abs(default_tol_for_vec_equality)

end subroutine get_geom_unfolding_relations


subroutine define_pckpts_to_be_checked(&
               pckpts_to_be_checked, &
               all_dirs_used_for_EBS_along_pcbz_dir, &
               nkpts_selected_dirs &
           )
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(selected_pcbz_directions), intent(out) :: pckpts_to_be_checked
type(irr_bz_directions), dimension(:), intent(in) :: &
    all_dirs_used_for_EBS_along_pcbz_dir
integer, dimension(:), intent(in) :: nkpts_selected_dirs
integer :: n_selec_pcbz_dirs, i_selec_pcbz_dir, nk, ikpt, i_needed_dirs, &
           n_dirs_for_EBS_along_pcbz_dir, alloc_stat
real(kind=dp), dimension(1:3) :: kstart,kend
type(vec3d), dimension(:), allocatable :: kline

    n_selec_pcbz_dirs = size(all_dirs_used_for_EBS_along_pcbz_dir)
    deallocate(pckpts_to_be_checked%selec_pcbz_dir,stat=alloc_stat)
    allocate(pckpts_to_be_checked%selec_pcbz_dir(1:n_selec_pcbz_dirs))
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        n_dirs_for_EBS_along_pcbz_dir = &
            size(&
                all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%&
                    irr_dir(:)&
            )
        deallocate(&
            pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir, &
            stat=alloc_stat&
        )
        allocate(&
            pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(1:n_dirs_for_EBS_along_pcbz_dir) &
        )
        nk = nkpts_selected_dirs(i_selec_pcbz_dir)
        do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir
            deallocate(&
                pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%&
                    needed_dir(i_needed_dirs)%pckpt, &
                stat=alloc_stat&
            )
            allocate(&
                pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%&
                    needed_dir(i_needed_dirs)%pckpt(1:nk), &
                stat=alloc_stat&
            )
            kstart = all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%&
                         irr_dir(i_needed_dirs)%kstart  
            kend = all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%&
                       irr_dir(i_needed_dirs)%kend
            deallocate(kline,stat=alloc_stat)
            allocate(kline(1:nk))
            kline = kpts_line(kstart,kend,nk)
            do ikpt=1,nk
                pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%&
                    needed_dir(i_needed_dirs)%pckpt(ikpt)%coords(:) &
                    = kline(ikpt)%coord(:)
            enddo
        enddo
    enddo


end subroutine define_pckpts_to_be_checked


subroutine select_coeffs_to_calc_spectral_weights(&
               selected_coeff_indices, wf, crystal_pc, GUR &
           )
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
integer, dimension(:), allocatable, intent(out) :: selected_coeff_indices
type(pw_wavefunction), intent(in) :: wf
type(crystal_3D), intent(in) :: crystal_pc
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR
! Local variables
integer :: alloc_stat, ig, i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt
real(kind=dp), dimension(1:3) :: trial_pc_g, folding_G, current_SCKPT, &
                                 symmetrized_unf_pc_kpt
real(kind=dp) :: tol, tol_for_vec_equality


    i_SCKPT = GUR%current_index%i_SCKPT
    i_selec_pcbz_dir = GUR%current_index%i_selec_pcbz_dir
    i_needed_dirs = GUR%current_index%i_needed_dirs
    ipc_kpt = GUR%current_index%ipc_kpt

    current_SCKPT = GUR%list_of_SCKPTS(i_SCKPT)%coord(:)
    symmetrized_unf_pc_kpt(:) = GUR%SCKPT(i_SCKPT)%&
                                    selec_pcbz_dir(i_selec_pcbz_dir)%&
                                    needed_dir(i_needed_dirs)%&
                                    pckpt(ipc_kpt)%Scoords(:)
    folding_G(:) = symmetrized_unf_pc_kpt(:) - current_SCKPT(:)

    tol_for_vec_equality = default_tol_for_vec_equality
    ! tol changed from 0.1 to 0.9 on 2016/01/15
    tol = 0.9_dp * tol_for_vec_equality
    deallocate(selected_coeff_indices, stat=alloc_stat)
    do while((.not. allocated(selected_coeff_indices)) .and. &
              (tol <= max_tol_for_vec_equality))
        ! changed from tol=10_dp*tol, 2016/01/15
        tol = tol + 0.1_dp*tol_for_vec_equality
        !$omp parallel do &
        !$omp schedule(guided) default(none) &
        !$omp shared(wf, folding_G, crystal_pc, tol, selected_coeff_indices) &
        !$omp private(ig, trial_pc_g) 
        do ig=1, wf%n_pw
            trial_pc_g = wf%G_cart(ig)%coord(:) - folding_G(:)
            if(vec_in_latt(&
                   vec=trial_pc_g, latt=crystal_pc%rec_latt_vecs, &
                   tolerance=tol &
               ))then
                 !$omp critical
                 call append(item=ig,list=selected_coeff_indices)
                 !$omp end critical
             endif
        enddo
    enddo
    if(abs(tol - tol_for_vec_equality) > epsilon(1.0_dp) .and. &
       allocated(selected_coeff_indices))then
        write(*,'(3(A,/),2(A,f0.6),A,/,A)') &
'    WARNING (select_coeffs_to_calc_spectral_weights): Problems selecting ',& 
'            the coeffs to calc. the spec. weights for the current pc-kpt.', &
'            The tolerace for testing vector equality had to be increased ', &
'            from ', tol_for_vec_equality,' to ',tol,'.', &
'            Please remember to double-check that your results are consistent.'
    endif

end subroutine select_coeffs_to_calc_spectral_weights


function spectral_weight_for_coeff(&
             coeff, selected_coeff_indices, add_elapsed_time_to &
         ) &
result(spectral_weight)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
! coeff = coeff(ikpt, iband)
complex(kind=kind_cplx_coeffs), dimension(:,:), intent(in) :: coeff
real(kind=dp), dimension(1:size(coeff,dim=2)) :: spectral_weight
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: n_bands_SC_calculation, iband, i, icoeff
real(kind=dp) :: stime, ftime

    stime = time_now()

    n_bands_SC_calculation = size(coeff,dim=2)
    spectral_weight(:) = 0.0_dp
    !$omp parallel do &
    !$omp schedule(guided) default(none) &
    !$omp private(iband,i,icoeff) &
    !$omp shared(n_bands_SC_calculation,selected_coeff_indices, &
    !$omp        spectral_weight,coeff)
    do iband = 1, n_bands_SC_calculation
        do i=1, size(selected_coeff_indices)
            icoeff = selected_coeff_indices(i)
            spectral_weight(iband) = spectral_weight(iband) + &
                                     abs(coeff(icoeff,iband))**2.0_dp
        enddo
    enddo

    if(present(add_elapsed_time_to))then
        ftime = time_now()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end function spectral_weight_for_coeff


subroutine calc_spectral_function(&
               SF_at_pckpt,energies,SC_calc_ener,spectral_weight,std_dev, &
               add_elapsed_time_to &
           )
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: SF_at_pckpt
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight, &
                                           SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in), optional :: std_dev
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, stime, ftime, sigma

    stime = time_now()

    sigma = 0.025_dp
    if(present(std_dev)) sigma = std_dev
    ! Calculating the spectral function A(pc_kpt,E) 
    n_SC_bands = size(spectral_weight)
    deallocate(SF_at_pckpt, stat=alloc_stat)
    allocate(SF_at_pckpt(1:size(energies)))
    SF_at_pckpt(:) = 0.0_dp
    !$omp parallel do &
    !$omp schedule(guided) default(none) &
    !$omp private(iener,E,i_SC_band,E_SC_band)  &
    !$omp shared(energies,n_SC_bands,SC_calc_ener,SF_at_pckpt, &
    !$omp        spectral_weight,sigma)
    do iener=1,size(energies)  
        ! For E in the list of energy values for the E(k) plot
        E = energies(iener)
        do i_SC_band=1,n_SC_bands
            E_SC_band = SC_calc_ener(i_SC_band)
            SF_at_pckpt(iener) = SF_at_pckpt(iener) + &
                                 spectral_weight(i_SC_band) * &
                                 delta(E_SC_band - E, std_dev=sigma)
        enddo
    enddo

    if(present(add_elapsed_time_to))then
        ftime = time_now()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine calc_spectral_function


function lambda(pc_ener, SC_ener, delta_e, std_dev) result(rtn)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp) :: rtn
real(kind=dp), intent(in) :: pc_ener, SC_ener, delta_e
real(kind=dp), intent(in), optional :: std_dev
real(kind=dp) :: sigma

    sigma = epsilon(1.0_dp)
    if(present(std_dev)) sigma = std_dev
    rtn = integral_delta_x_minus_x0(x0=SC_ener, &
                                    lower_lim=pc_ener - 0.5_dp*delta_e, &
                                    upper_lim=pc_ener + 0.5_dp*delta_e, &
                                    std_dev=sigma)

end function lambda


subroutine calc_delta_N_pckpt(&
               delta_N_pckpt, energies, SC_calc_ener, &
               spectral_weight, std_dev &
           )
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: delta_N_pckpt
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in), optional :: std_dev
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, delta_e, sigma

    
    sigma = epsilon(1.0_dp)
    if(present(std_dev)) sigma = std_dev
    delta_e = energies(2)  - energies(1) !! Uniform energy grid 
    n_SC_bands = size(spectral_weight)
    deallocate(delta_N_pckpt, stat=alloc_stat)
    allocate(delta_N_pckpt(1:size(energies)))
    delta_N_pckpt(:) = 0.0_dp
    !$omp parallel do &
    !$omp schedule(guided) default(none) &
    !$omp private(iener,E,i_SC_band,E_SC_band)  &
    !$omp shared(energies,n_SC_bands,SC_calc_ener,delta_N_pckpt, &
    !$omp        spectral_weight, delta_e, sigma)
    do iener=1,size(energies)  
        ! For E in the list of energy values for the E(k) plot
        E = energies(iener)
        do i_SC_band=1,n_SC_bands
            E_SC_band = SC_calc_ener(i_SC_band)
            delta_N_pckpt(iener) = delta_N_pckpt(iener) + &
                                   spectral_weight(i_SC_band) * &
                                   lambda(pc_ener=E, SC_ener=E_SC_band, &
                                          delta_e=delta_e, std_dev=sigma)
        enddo
    enddo

end subroutine calc_delta_N_pckpt


subroutine get_delta_Ns_for_EBS(&
               delta_N_pckpt, energy_grid, SC_calc_ener, spectral_weight, &
               smearing_for_delta_function, add_elapsed_time_to &
           )
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: delta_N_pckpt
real(kind=dp), dimension(:), intent(in) :: energy_grid, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in), optional :: smearing_for_delta_function
real(kind=dp) :: sigma, stime, ftime
real(kind=dp), intent(inout), optional :: add_elapsed_time_to

    stime = time_now()
    sigma = epsilon(1.0_dp)
    if(present(smearing_for_delta_function))then
        sigma = max(sigma,abs(smearing_for_delta_function))
    endif
    call calc_delta_N_pckpt(&
             delta_N_pckpt, energy_grid, SC_calc_ener, spectral_weight, sigma &
         )
    if(present(add_elapsed_time_to))then
        ftime = time_now()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine get_delta_Ns_for_EBS


subroutine get_delta_Ns_for_output(&
               delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS, &
               delta_N, all_dirs_used_for_EBS_along_pcbz_dir, &
               pckpts_to_be_checked, times &
           )
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(UnfoldedQuantitiesForOutput), target, intent(out) :: &
    delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
type(UnfoldedQuantities), target, intent(in) :: delta_N
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(irr_bz_directions), dimension(:), intent(in) :: &
    all_dirs_used_for_EBS_along_pcbz_dir
type(timekeeping), intent(inout), optional :: times
integer :: nener, i_selec_pcbz_dir, ipc_kpt, i_needed_dirs, iener, i_rho, &
           n_rhos, m1, m2, not_avgd_unf_dens_el_1D_index, &
           avgd_unf_dens_el_1D_index
real(kind=dp), dimension(:), allocatable :: avrgd_dNs, avrgd_spin_proj, &
                                            avrgd_parallel_proj
real(kind=dp), dimension(:,:), allocatable :: avrgd_sigma
real(kind=dp) :: weight, stime, ftime
logical :: output_spin_info
real(kind=dp), dimension(:), pointer :: symm_avg_dNs
type(UnfoldDensityOpContainer), pointer :: avgd_rho, not_avgd_rho
type(UnfoldDensityOpContainer), dimension(:), pointer :: symm_avg_rhos


    stime = time_now()
    nener = size(delta_N%selec_pcbz_dir(1)%needed_dir(1)%pckpt(1)%dN(:))
    allocate(avrgd_dNs(1:nener))
    call allocate_UnfoldedQuantitiesForOutput(&
             delta_N_only_selected_dirs, pckpts_to_be_checked &
         )
    call allocate_UnfoldedQuantitiesForOutput(&
             delta_N_symm_avrgd_for_EBS, pckpts_to_be_checked &
         )

    delta_N_only_selected_dirs%n_SC_bands = delta_N%n_SC_bands
    delta_N_symm_avrgd_for_EBS%n_SC_bands = delta_N%n_SC_bands
    do i_selec_pcbz_dir=1,size(delta_N%selec_pcbz_dir(:))
        delta_N_only_selected_dirs%pcbz_dir(i_selec_pcbz_dir) = &
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)
        do ipc_kpt=1, size(&
                          delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                              needed_dir(1)%pckpt(:) &
                      )
            avrgd_dNs(:) = 0.0_dp
            allocate(&
                delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                    pckpt(ipc_kpt)%dN(1:nener) &
            )
            do i_needed_dirs=1, size(&
                                    delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                                        needed_dir(:) &
                                )
                weight=all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%&
                           irr_dir(i_needed_dirs)%weight
                avrgd_dNs = avrgd_dNs + &
                            weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                                       needed_dir(i_needed_dirs)%&
                                       pckpt(ipc_kpt)%dN(:)
            enddo
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                pckpt(ipc_kpt)%dN(:) &
                = avrgd_dNs
       enddo     
    enddo
    ftime = time_now()
    if(present(times))then
        times%calc_dN = times%calc_dN + (ftime - stime)
    endif

    if(args%write_unf_dens_op)then
        stime = time_now()
        ! Symmetry-averaging the unfolding density operator
        do i_selec_pcbz_dir=1, size(delta_N%selec_pcbz_dir(:))
            do ipc_kpt=1, size(&
                              delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                                  needed_dir(1)%pckpt(:) &
                          )
                symm_avg_dNs => delta_N_symm_avrgd_for_EBS%&
                                    pcbz_dir(i_selec_pcbz_dir)%&
                                    pckpt(ipc_kpt)%dN
                allocate(&
                    delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                        pckpt(ipc_kpt)%rhos(1:nener) &
                )
                symm_avg_rhos => delta_N_symm_avrgd_for_EBS%&
                                     pcbz_dir(i_selec_pcbz_dir)%&
                                     pckpt(ipc_kpt)%rhos

                do i_needed_dirs=1, size(&
                                        delta_N%&
                                            selec_pcbz_dir(i_selec_pcbz_dir)%&
                                            needed_dir(:) &
                                    )
                    weight = &
                        all_dirs_used_for_EBS_along_pcbz_dir(&
                            i_selec_pcbz_dir &
                        )%irr_dir(i_needed_dirs)%weight
                    n_rhos = size(&
                                 delta_N%selec_pcbz_dir(i_selec_pcbz_dir)% &
                                     needed_dir(i_needed_dirs)% &
                                     pckpt(ipc_kpt)%rhos &
                             )
                    do i_rho=1,n_rhos
                        ! The unf. dens. ops. are nbandsXnbands matrices, 
                        ! so I have not stored (not even calculated) them at 
                        ! PC energy grid points at which delta_N(k,E) is 
                        ! near zero. This is why I have this
                        ! i_rho --> iener map here.
                        iener = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)% &
                                    needed_dir(i_needed_dirs)% &
                                    pckpt(ipc_kpt)%rhos(i_rho)%&
                                    iener_in_full_pc_egrid
                        if(symm_avg_dNs(iener) < 1E-3) cycle
                        not_avgd_rho => delta_N%&
                                            selec_pcbz_dir(i_selec_pcbz_dir)% &
                                            needed_dir(i_needed_dirs)% &
                                            pckpt(ipc_kpt)%rhos(i_rho)
                        if(.not. allocated(not_avgd_rho%band_indices)) cycle 
                        avgd_rho => symm_avg_rhos(iener)
                        avgd_rho%iener_in_full_pc_egrid = iener
                        avgd_rho%nbands = delta_N%n_SC_bands
                        do m1=1,delta_N%n_SC_bands
                            do m2=1,delta_N%n_SC_bands
                                not_avgd_unf_dens_el_1D_index = &
                                    list_index(&
                                        [m1,m2], not_avgd_rho%band_indices &
                                    )
                                ! Nothing to add if the [m1,m2] matrix element 
                                ! of the "non-avgd" operator does not exist 
                                ! This means it's zero. The op is very sparse.
                                if(not_avgd_unf_dens_el_1D_index<1) cycle
                                if(.not. allocated(avgd_rho%band_indices))then
                                    avgd_unf_dens_el_1D_index = 0
                                else
                                    avgd_unf_dens_el_1D_index = &
                                        list_index(&
                                            item=[m1,m2], &
                                            list=avgd_rho%band_indices &
                                        )
                                endif
                                if(avgd_unf_dens_el_1D_index<1)then
                                    ! Non-avgd op exists, but no corresponding
                                    ! entry (matrix element) for the averaged 
                                    ! op has been calculated.
                                    ! The "append" subroutine allocates 
                                    ! argument arrays if not already allocated
                                    call append(item=&
                                           cmplx(&
                                             weight * &
                                             not_avgd_rho%rho(&
                                               not_avgd_unf_dens_el_1D_index),&
                                             kind=kind_cplx_coeffs &
                                           ), &
                                           list=avgd_rho%rho &
                                         )
                                    call append(&
                                             item=[m1,m2], &
                                             list=avgd_rho%band_indices &
                                         )
                                else
                                    ! Avgd op matrix el. [m1, m2] already 
                                    ! started to be calcd, but will generally
                                    ! have a different linear index in the 
                                    ! symm-avgd and no symm-avgd ops
                                    avgd_rho%rho(avgd_unf_dens_el_1D_index) = &
                                        avgd_rho%&
                                            rho(avgd_unf_dens_el_1D_index) + &
                                        weight* &
                                        not_avgd_rho%&
                                            rho(not_avgd_unf_dens_el_1D_index)
                                endif
                            enddo ! SC band index m2
                        enddo ! SC band index m2
                    enddo ! i_rho index
                enddo ! i_needed_dirs index
            enddo ! i_pc_kpt
        enddo ! i_selec_pcbz_dir index
        ftime = time_now()
        if(present(times))then
            times%calc_rho = times%calc_rho + (ftime - stime)
        endif
    endif

    output_spin_info = allocated(&
                           delta_N%selec_pcbz_dir(1)%needed_dir(1)%&
                               pckpt(1)%spin_proj_perp &
                       )
    if(.not. output_spin_info) return
    stime = time_now()
    allocate(avrgd_spin_proj(1:nener))
    allocate(avrgd_parallel_proj(1:nener))
    allocate(avrgd_sigma(1:nener,1:3))
    do i_selec_pcbz_dir=1,size(delta_N%selec_pcbz_dir(:))
        delta_N_only_selected_dirs%pcbz_dir(i_selec_pcbz_dir) = &
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)
        do ipc_kpt=1, size(&
                          delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                              needed_dir(1)%pckpt(:) &
                      )
            avrgd_spin_proj(:) = 0.0_dp
            avrgd_parallel_proj(:) = 0.0_dp
            avrgd_sigma(:,:) = 0.0_dp
            allocate(&
                delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                    pckpt(ipc_kpt)%spin_proj_perp(1:nener) &
            )
            allocate(&
                delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                    pckpt(ipc_kpt)%spin_proj_para(1:nener) &
            )
            allocate(&
                delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                    pckpt(ipc_kpt)%sigma(1:nener,1:3) &
            )
            do i_needed_dirs=1, size(&
                                    delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                                        needed_dir(:) &
                                )
                weight=all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%&
                           irr_dir(i_needed_dirs)%weight
                avrgd_spin_proj = avrgd_spin_proj + weight * &
                                  delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                                      needed_dir(i_needed_dirs)%&
                                      pckpt(ipc_kpt)%spin_proj_perp(:)
                avrgd_parallel_proj = avrgd_parallel_proj + weight * &
                                      delta_N%&
                                          selec_pcbz_dir(i_selec_pcbz_dir)%&
                                          needed_dir(i_needed_dirs)%&
                                          pckpt(ipc_kpt)%spin_proj_para(:)
                avrgd_sigma = avrgd_sigma + weight * &
                              delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                                  needed_dir(i_needed_dirs)%&
                                  pckpt(ipc_kpt)%sigma(:,:)
            enddo
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                pckpt(ipc_kpt)%spin_proj_perp(:) &
                = avrgd_spin_proj(:)
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                pckpt(ipc_kpt)%spin_proj_para(:) &
                = avrgd_parallel_proj(:)
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%&
                pckpt(ipc_kpt)%sigma(:,:) &
                = avrgd_sigma(:,:)
       enddo     
    enddo
    ftime = time_now()
    if(present(times))then
        times%calc_pauli_vec_projs = times%calc_pauli_vec_projs + &
                                     (ftime - stime)
    endif

end subroutine get_delta_Ns_for_output


subroutine calc_rho(rho, delta_N, pc_ener, delta_e, wf, &
                    selected_coeff_indices, std_dev, add_elapsed_time_to)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! rho is the 'unfolding-density operator', defined in
!! Phys. Rev. B. 91, 041116(R) (2015)
!! Please make sure that delta_N /= 0 before calling this subroutine
implicit none
complex(kind=kind_cplx_coeffs), dimension(:,:), allocatable, intent(out) :: rho
real(kind=dp), intent(in) :: delta_N, pc_ener, delta_e
type(pw_wavefunction), intent(in) :: wf
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), intent(in), optional :: std_dev
real(kind=dp) :: min_prefactor, lamb_ee1, lamb_ee2, prefactor, sigma, &
                 stime, ftime
! m1, m2: SC band ind
integer :: m1, m2, n_SC_bands, pw_index, alloc_stat, n_coeffs, n_coeffs_pf, & 
           n_spinor, i_spinor, pf_pw_index, n_picked_bands
! pf = partial function
complex(kind=kind_cplx_coeffs), dimension(:), allocatable :: pf_coeffs_m1, &
                                                             pf_coeffs_m2 
real(kind=dp), intent(inout), optional :: add_elapsed_time_to


    stime = time_now()
    if(delta_N < epsilon(1.0_dp))then
        write(*,'(A)')'ERROR (calc_rho): delta_N = 0.'
        write(*,'(A)')'Stopping now.'
        stop
    endif

    min_prefactor = 1E-5_dp
    sigma = epsilon(1.0_dp)
    if(present(std_dev)) sigma = std_dev

    n_spinor = size(wf%pw_coeffs, dim=1)
    n_coeffs = size(wf%pw_coeffs, dim=2)
    n_SC_bands = size(wf%pw_coeffs, dim=3)
    n_coeffs_pf = size(selected_coeff_indices)
    deallocate(rho, stat=alloc_stat)
    allocate(rho(1:n_SC_bands, 1:n_SC_bands))

    rho(:,:) = 0.0_dp
    if(args%perform_unfold)then
        !$omp parallel do &
        !$omp schedule(guided) default(none) &
        !$omp shared(n_SC_bands, pc_ener, wf, delta_e, sigma, &
        !$omp        min_prefactor, delta_N, n_coeffs_pf, &
        !$omp        n_spinor, selected_coeff_indices, rho) &
        !$omp private(m1, m2, lamb_ee1, lamb_ee2, prefactor, pf_coeffs_m1, &
        !$omp         pf_coeffs_m2, alloc_stat, i_spinor, pf_pw_index)
        do m1=1, n_SC_bands - 1
            lamb_ee1 = lambda(pc_ener=pc_ener, SC_ener=wf%band_energies(m1), &
                              delta_e=delta_e, std_dev=sigma)
            if(lamb_ee1 < min_prefactor) cycle
            do m2=m1, n_SC_bands
                lamb_ee2 = lambda(&
                               pc_ener=pc_ener, SC_ener=wf%band_energies(m2), &
                               delta_e=delta_e, std_dev=sigma &
                           )
                prefactor = (lamb_ee1 * lamb_ee2) / delta_N
                if(prefactor < min_prefactor) cycle
            
                allocate(&
                    pf_coeffs_m1(1:n_coeffs_pf), pf_coeffs_m2(1:n_coeffs_pf) &
                )
                do i_spinor=1, n_spinor
                    pf_pw_index = 0
                    do pw_index = 1, n_coeffs_pf
                        pf_pw_index = pf_pw_index + 1
                        pf_coeffs_m1(pf_pw_index) = &
                            wf%pw_coeffs(&
                                   i_spinor, &
                                   selected_coeff_indices(pw_index), &
                                   m1 &
                               )
                        pf_coeffs_m2(pf_pw_index) = &
                            wf%pw_coeffs(&
                                   i_spinor, &
                                   selected_coeff_indices(pw_index), &
                                   m2 &
                               )
                    enddo
                    rho(m1,m2) = rho(m1,m2) + &
                                 dot_product(pf_coeffs_m1, pf_coeffs_m2)
                enddo
                deallocate(pf_coeffs_m1, pf_coeffs_m2, stat=alloc_stat)
                rho(m1,m2) = prefactor * rho(m1,m2)
                rho(m2,m1) = conjg(rho(m1,m2))
            enddo
        enddo
    else
        n_picked_bands = 0
        do m1=1, n_SC_bands
            lamb_ee1 = lambda(pc_ener=pc_ener, SC_ener=wf%band_energies(m1), &
                              delta_e=delta_e, std_dev=sigma)
            if(lamb_ee1 < min_prefactor) cycle
            n_picked_bands = n_picked_bands + 1
            do m2=1, n_SC_bands
                lamb_ee2 = lambda(&
                               pc_ener=pc_ener, SC_ener=wf%band_energies(m2), &
                               delta_e=delta_e, std_dev=sigma &
                           )
                prefactor = (lamb_ee1 * lamb_ee2) / delta_N
                if(prefactor < min_prefactor) cycle
                rho(m1,m2) = 1.0_dp
            enddo
        enddo
        if(n_picked_bands > 0) rho = rho / real(n_picked_bands, kind=dp)
    endif

    if(present(add_elapsed_time_to))then
        ftime = time_now()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine calc_rho


subroutine get_SC_pauli_matrix_elmnts(&
               SC_pauli_matrix_elmnts, pauli_mtx_elmts_already_calc, &
               wf, selected_pc_ener, delta_e, add_elapsed_time_to &
           )
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
complex(kind=kind_cplx_coeffs), dimension(:,:,:), &
    allocatable, intent(inout) :: SC_pauli_matrix_elmnts
logical, dimension(:,:), allocatable, intent(inout) :: &
    pauli_mtx_elmts_already_calc
type(pw_wavefunction), intent(in) :: wf
real(kind=dp), intent(in) :: selected_pc_ener, delta_e
complex(kind=kind_cplx_coeffs) :: inner_prod_Km2Alpha_Km1Beta, &
                                  inner_prod_Km2Beta_Km1Alpha, &
                                  inner_prod_Km2Alpha_Km1Alpha, &
                                  inner_prod_Km2Beta_Km1Beta
complex(kind=kind_cplx_coeffs), parameter :: &
    J = (0.0_kind_cplx_coeffs, 1.0_kind_cplx_coeffs), &
    cplx_zero = (0.0_kind_cplx_coeffs, 0.0_kind_cplx_coeffs)
integer :: n_SC_bands, alloc_stat, m1, m2, m3, m4
! list of SC bands to enter the calc of the matrix elements
integer, dimension(:), allocatable :: picked_SC_bands 
real(kind=dp) :: lambd, stime, ftime
real(kind=dp), intent(inout), optional :: add_elapsed_time_to


    stime = time_now()
    n_SC_bands = size(wf%pw_coeffs, dim=3)
    if(.not. allocated(pauli_mtx_elmts_already_calc))then
        allocate(pauli_mtx_elmts_already_calc(1:n_SC_bands, 1:n_SC_bands))
        pauli_mtx_elmts_already_calc(:,:) = .FALSE.
    endif
    if(.not. any(pauli_mtx_elmts_already_calc))then
        deallocate(SC_pauli_matrix_elmnts, stat=alloc_stat)
        allocate(SC_pauli_matrix_elmnts(1:3, 1:n_SC_bands, 1:n_SC_bands))
        SC_pauli_matrix_elmnts(:,:,:) = cplx_zero
    endif

    do m1=1, n_SC_bands
        lambd = lambda(&
                    pc_ener=selected_pc_ener, SC_ener=wf%band_energies(m1), &
                    delta_e=delta_e &
                )
        if(lambd > 1E-1)then
            call append(list=picked_SC_bands, item=m1)
        endif
    enddo
    if(allocated(picked_SC_bands))then
        !$omp parallel do collapse(2) &
        !$omp schedule(guided) default(none) &
        !$omp shared(wf, SC_pauli_matrix_elmnts, picked_SC_bands, &
        !$omp        pauli_mtx_elmts_already_calc) &
        !$omp private(m1, m2, m3, m4, inner_prod_Km2Alpha_Km1Beta, &
        !$omp         inner_prod_Km2Beta_Km1Alpha, &
        !$omp         inner_prod_Km2Alpha_Km1Alpha, inner_prod_Km2Beta_Km1Beta)
        do m3=1,size(picked_SC_bands)
            do m4=1,size(picked_SC_bands)
                m1 = picked_SC_bands(m3)
                m2 = picked_SC_bands(m4)
                if(pauli_mtx_elmts_already_calc(m2,m1)) cycle
                inner_prod_Km2Alpha_Km1Beta = dot_product(&
                                                  wf%pw_coeffs(1,:,m2), &
                                                  wf%pw_coeffs(2,:,m1) &
                                              )
                inner_prod_Km2Beta_Km1Alpha = dot_product(&
                                                  wf%pw_coeffs(2,:,m2), &
                                                  wf%pw_coeffs(1,:,m1) &
                                              )
                inner_prod_Km2Alpha_Km1Alpha = dot_product(&
                                                   wf%pw_coeffs(1,:,m2), &
                                                   wf%pw_coeffs(1,:,m1) &
                                               )
                inner_prod_Km2Beta_Km1Beta = dot_product(&
                                                 wf%pw_coeffs(2,:,m2), &
                                                 wf%pw_coeffs(2,:,m1) &
                                             )

                SC_pauli_matrix_elmnts(1,m2,m1) = inner_prod_Km2Alpha_Km1Beta+&
                                                  inner_prod_Km2Beta_Km1Alpha
                SC_pauli_matrix_elmnts(2,m2,m1) = &
                    J*(inner_prod_Km2Beta_Km1Alpha-inner_prod_Km2Alpha_Km1Beta)
                SC_pauli_matrix_elmnts(3,m2,m1) = &
                    inner_prod_Km2Alpha_Km1Alpha - inner_prod_Km2Beta_Km1Beta
                pauli_mtx_elmts_already_calc(m2,m1) = .TRUE.

                SC_pauli_matrix_elmnts(:,m1,m2) = &
                    conjg(SC_pauli_matrix_elmnts(:,m2,m1))
                pauli_mtx_elmts_already_calc(m1,m2) = .TRUE.
            enddo
        enddo
        deallocate(picked_SC_bands)
    endif

    if(present(add_elapsed_time_to))then
        ftime = time_now()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif


end subroutine get_SC_pauli_matrix_elmnts


subroutine calc_spin_projections(spin_proj_perp, spin_proj_para, &
                                 pauli_vector, cart_coords_kpt, origin)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! This is a BETA feature, and should be used with caution
implicit none
real(kind=dp), intent(out) :: spin_proj_perp, spin_proj_para
real(kind=dp), dimension(1:3), intent(in) :: pauli_vector, cart_coords_kpt
! In cartesian coordinates
real(kind=dp), dimension(1:3), intent(in), optional :: origin 
real(kind=dp), dimension(1:3) :: versor_para, versor_perp, &
                                 normal2proj_plane, cart_coords_origin

    cart_coords_origin = 0.0_dp
    if(present(origin)) cart_coords_origin = origin
    normal2proj_plane = versor(args%normal_to_proj_plane)
    versor_para = versor(cart_coords_kpt - cart_coords_origin)
    versor_perp = versor(cross(normal2proj_plane, versor_para)) 
    spin_proj_perp = dot_product(pauli_vector, versor_perp)
    spin_proj_para = dot_product(pauli_vector, versor_para)

end subroutine calc_spin_projections


subroutine perform_unfolding(&
               delta_N, times, GUR, wf, selected_coeff_indices, energy_grid &
           )
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none 
type(UnfoldedQuantities), target, intent(inout) :: delta_N
type(timekeeping), intent(inout) :: times
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR
type(pw_wavefunction), intent(in) :: wf
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), dimension(:), intent(in) :: energy_grid

integer :: n_WF_components, n_bands_SC_calculation, i_spinor, iener, nener, &
           iener2, alloc_stat, i, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt, &
           i_SCKPT, m1, m2, unf_dens_file_unit
real(kind=dp), dimension(:), allocatable :: spectral_weight
real(kind=dp) :: selected_pc_ener, delta_e, dN, unf_spin_proj_perp, &
                 unf_spin_proj_para
real(kind=dp), dimension(1:3) :: symmetrized_unf_pc_kpt, true_unf_pc_kpt, &
                                 unf_pauli_vec, projection_origin
complex(kind=kind_cplx_coeffs), dimension(:,:), allocatable :: temp_rho_matrix
type(UnfoldDensityOpContainer), dimension(:), pointer :: rhos
complex(kind=kind_cplx_coeffs), dimension(:,:,:), allocatable :: &
    SC_pauli_matrix_elmnts
integer, dimension(:), allocatable :: pc_energies_to_calc_eigenvalues
logical, dimension(:,:), allocatable :: pauli_mtx_elmts_already_calc


    i_SCKPT = GUR%current_index%i_SCKPT
    i_selec_pcbz_dir = GUR%current_index%i_selec_pcbz_dir
    i_needed_dirs = GUR%current_index%i_needed_dirs
    ipc_kpt = GUR%current_index%ipc_kpt
    ! In cartesion coors!
    symmetrized_unf_pc_kpt(:) = GUR%SCKPT(i_SCKPT)%&
                                    selec_pcbz_dir(i_selec_pcbz_dir)%&
                                    needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                                    Scoords(:)
    ! In cartesion coors!
    true_unf_pc_kpt(:) = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%&
                             needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords(:)

    ! 1 for a regular WF, and 2 for a spinor WF
    n_WF_components = size(wf%pw_coeffs, dim=1)
    n_bands_SC_calculation = size(wf%pw_coeffs, dim=3)
    allocate(spectral_weight(1:n_bands_SC_calculation))

    if(args%perform_unfold)then
        spectral_weight = 0.0_dp
        do i_spinor=1,n_WF_components
           ! Calculating spectral_weights
            spectral_weight = spectral_weight + &
                              spectral_weight_for_coeff(&
                                  wf%pw_coeffs(i_spinor,:,:), &
                                  selected_coeff_indices, &
                                  add_elapsed_time_to=times%calc_spec_weights &
                              )
        enddo
    else
        write(*,'(A)')&
            '    **** WARNING: You are using the option "-dont_unfold". &
                               You are NOT performing unfolding!! ****'
        spectral_weight = 1.0_dp
    endif
    
    ! Calculating the delta_Ns
    delta_N%n_SC_bands = n_bands_SC_calculation
    call get_delta_Ns_for_EBS(&
             delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                 needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN, &
             energy_grid, wf%band_energies, spectral_weight, &
             add_elapsed_time_to=times%calc_dN &
         )

    if(calc_spec_func_explicitly)then
        ! Calculating the spectral_function (optional)
        call calc_spectral_function(&
                 delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                     needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%SF, &
                 energy_grid, wf%band_energies, spectral_weight, &
                 add_elapsed_time_to=times%calc_SF &
             )
    endif

    if(args%write_unf_dens_op .or. n_WF_components==2)then
        ! Points in the PC energy grid at which the 
        ! unfolding-density op will be calcd.
        nener = size(&
                    delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                        needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN &
                )
        delta_e = energy_grid(2) - energy_grid(1) ! uniform energy grid
        !$omp parallel do &
        !$omp schedule(guided) default(none) &
        !$omp shared(nener, delta_N, i_selec_pcbz_dir, i_needed_dirs, &
        !$omp        ipc_kpt, pc_energies_to_calc_eigenvalues) &
        !$omp private(iener, dN)
        do iener=1, nener
            dN = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                     needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(iener)
            if(dN >= 1E-3_dp)then 
                ! Can hardly be considered a band otherwise;
                ! no point calculating rho
                !$omp critical
                call append(list=pc_energies_to_calc_eigenvalues, item=iener)
                !$omp end critical
            endif
        enddo

        ! Calculating the unfolding-density operator at the selected energies
        allocate(&
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                rhos(1:size(pc_energies_to_calc_eigenvalues)) &
        )
        rhos => delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                    needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%rhos
        do iener2=1, size(pc_energies_to_calc_eigenvalues)
            iener = pc_energies_to_calc_eigenvalues(iener2)
            selected_pc_ener = energy_grid(iener)
            dN = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                     needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(iener)
            call calc_rho(temp_rho_matrix, &
                          dN, energy_grid(iener), delta_e, wf, &
                          selected_coeff_indices, &
                          add_elapsed_time_to=times%calc_rho)
            rhos(iener2)%iener_in_full_pc_egrid = iener
            rhos(iener2)%nbands = n_bands_SC_calculation
            do m1=1,n_bands_SC_calculation
                do m2=1,n_bands_SC_calculation
                    if(abs(temp_rho_matrix(m1,m2)) < 1E-5) cycle
                    call append(&
                             item=temp_rho_matrix(m1,m2), &
                             list=rhos(iener2)%rho &
                         )
                    call append(item=[m1,m2], list=rhos(iener2)%band_indices)
                enddo
            enddo
            deallocate(temp_rho_matrix)
        enddo
    endif


    if(n_WF_components==2)then
        allocate(&
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                spin_proj_perp(1:nener) &
        )
        allocate(&
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                spin_proj_para(1:nener) &
        )
        allocate(&
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%sigma(1:nener,1:3) &
        )
        delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
            needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp(:) &
            = 0.0_dp 
        delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
            needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para(:) &
            = 0.0_dp 
        delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
            needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%sigma(:,:) &
            = 0.0_dp

        deallocate(&
            pauli_mtx_elmts_already_calc, SC_pauli_matrix_elmnts, &
            stat=alloc_stat &
        )
        do iener2=1, size(pc_energies_to_calc_eigenvalues)
            iener = pc_energies_to_calc_eigenvalues(iener2)
            selected_pc_ener = energy_grid(iener)
            dN = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                     needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(iener)

            call get_SC_pauli_matrix_elmnts(&
                     SC_pauli_matrix_elmnts, pauli_mtx_elmts_already_calc, &
                     wf, selected_pc_ener, delta_e, &
                     add_elapsed_time_to=times%calc_pauli_vec &
                 )
            do i=1,3
                unf_pauli_vec(i) = real(&
                                       trace_AB(&
                                           A=rhos(iener2), &
                                           B=SC_pauli_matrix_elmnts(i,:,:) &
                                       ), &
                                       kind=dp &
                                   )
            enddo
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%sigma(iener,:) &
                = unf_pauli_vec(:)

            projection_origin = GUR%SCKPT(i_SCKPT)%&
                                    selec_pcbz_dir(i_selec_pcbz_dir)%&
                                    needed_dir(i_needed_dirs)%&
                                    pckpt(ipc_kpt)%Sorigin_for_spin_proj(:)
            call calc_spin_projections(&
                     unf_spin_proj_perp, unf_spin_proj_para, unf_pauli_vec, &
                     symmetrized_unf_pc_kpt, origin=projection_origin &
                 )
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                spin_proj_perp(iener) &
                = unf_spin_proj_perp 
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%&
                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                spin_proj_para(iener) &
                = unf_spin_proj_para
            !! rho is allocated in the calc_rho routine
            !! rho(iener2) is no longer needed
            !!deallocate(rho(iener2)%rho, stat=alloc_stat) 
        enddo
        deallocate(&
            pauli_mtx_elmts_already_calc, SC_pauli_matrix_elmnts, &
            stat=alloc_stat &
        )
    endif

end subroutine perform_unfolding


end module band_unfolding
