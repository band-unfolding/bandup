!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!!
!! This file is part of BandUP: Band Unfolding code for Plane-wave based calculations.
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

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module band_unfolding
! ===================
! This module contains the routines that are directly related to the unfolding
! of the bands.
!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!===============================================================================
! MODULE: band_unfolding 
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION:
!> Contains all routines employed to perform unfolding in BandUP.
!===============================================================================

module band_unfolding
!$ use omp_lib
use constants_and_types
use cla_wrappers
use math
use general_io
use io_routines
implicit none
PRIVATE
PUBLIC :: get_geom_unfolding_relations, define_pckpts_to_be_checked, &
          select_coeffs_to_calc_spectral_weights, get_delta_Ns_for_output, &
          update_GUR_indices, perform_unfolding, verify_commens

CONTAINS


subroutine update_GUR_indices(GUR, i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt)
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

    call check_if_pc_and_SC_are_commensurate(are_commens, matrix_M, crystal_pc, crystal_SC, &
                                             tol=default_tol_for_int_commens_test)
    call print_message_commens_test(commensurate=are_commens,M=matrix_M,&
                                    stop_if_not_commens=args%stop_if_not_commensurate) 
    if(args%stop_if_not_commensurate .and. .not. are_commens) stop

end subroutine verify_commens


subroutine get_geom_unfolding_relations(GUR,list_of_SCKPTS, pckpts_to_be_checked, input_crystal_SC, verbose)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(geom_unfolding_relations_for_each_SCKPT), intent(out) :: GUR !! Geometric Unfolding Relations
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(vec3d), dimension(:), intent(in) :: list_of_SCKPTS
type(crystal_3D), intent(in) :: input_crystal_SC
logical, intent(in), optional :: verbose
integer :: nkpts, n_selec_pcbz_dirs, i_SCKPT,ipc_kpt,ieqv_SCKPT,isym, &
           i_selec_pcbz_dir,i_needed_dirs,alloc_stat, i
integer, dimension(:), allocatable :: n_dirs_for_EBS_along_pcbz_dir,n_pckpts_dirs
logical, dimension(:,:,:), allocatable :: pc_kpt_already_folded
real(kind=dp), dimension(1:3) :: pc_kpt, current_SCKPT, SCKPT_eqv_to_current_SCKPT, &
                                 trial_folding_G, origin_for_spin_proj
real(kind=dp), dimension(1:3,1:3) :: B_matrix_SC
type(star), dimension(:), allocatable :: SKPTS_eqv_to_SKPT
logical :: print_stuff
type(crystal_3D) :: crystal_SC

    crystal_SC = input_crystal_SC
    print_stuff = .FALSE.
    if(present(verbose))then
        print_stuff = verbose
    endif

    if(print_stuff)then
        write(*,"(A)")"Verifying geometric unfolding relations between pcbz and SCBZ wave-vectors... "
    endif

    B_matrix_SC = crystal_SC%rec_latt_vecs
    if(args%origin_for_spin_proj_passed_in_rec)then
        args%origin_for_spin_proj_cartesian(:) = 0.0_dp
        do i=1,3
            args%origin_for_spin_proj_cartesian = args%origin_for_spin_proj_cartesian + &
                args%origin_for_spin_proj_rec(i) * B_matrix_SC(i,:)
        enddo
    endif

    if(args%no_symm_sckpts)then
        if(print_stuff)then
            write(*,"(A)")'    * Running with the flag "-no_symm_sckpts". The symmetry of the SC'
            write(*,"(A)")'      will NOT be used to reduce the number of SC K-points needed.'
        endif
        crystal_SC%nsym = 1
        deallocate(crystal_SC%symops, stat=alloc_stat)
        allocate(crystal_SC%symops(1:1))
        crystal_SC%symops(1)%translation_fractional_coords(:) = 0.0_dp
        crystal_SC%symops(1)%translation_cartesian_coords(:) = 0.0_dp
        crystal_SC%symops(1)%rotation_fractional_coords = identity_3D
        crystal_SC%symops(1)%rotation_cartesian_coords = identity_3D
    else
        call get_symm(crystal=crystal_SC, use_pc_to_get_symm=.FALSE., symprec=default_symprec) ! fails if use_pc_to_get_symm=.TRUE.
    endif
    call get_star(star_of_pt=SKPTS_eqv_to_SKPT, points=list_of_SCKPTS, crystal=crystal_SC, &
                  tol_for_vec_equality=default_tol_for_vec_equality, &
                  symprec=default_symprec, reduce_to_bz=.TRUE.)
    !! Allocating and initializing table
    nkpts = size(list_of_SCKPTS)
    allocate(GUR%list_of_SCKPTS(1:nkpts))
    GUR%list_of_SCKPTS = list_of_SCKPTS
    n_selec_pcbz_dirs = size(pckpts_to_be_checked%selec_pcbz_dir(:))
    allocate(n_dirs_for_EBS_along_pcbz_dir(1:n_selec_pcbz_dirs))
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir) = size(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(:))
    enddo
    allocate(n_pckpts_dirs(1:n_selec_pcbz_dirs))
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        n_pckpts_dirs(i_selec_pcbz_dir) = size(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)%pckpt(:)) 
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
            deallocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir, stat=alloc_stat)
            allocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1:n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)))
            do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)
                deallocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt, stat=alloc_stat)
                allocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(1:n_pckpts_dirs(i_selec_pcbz_dir)))
                do ipc_kpt=1,n_pckpts_dirs(i_selec_pcbz_dir)
                    GUR%n_pckpts = GUR%n_pckpts + 1
                    pc_kpt(:) = pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords(:)
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords(:) = pc_kpt(:) 
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:) = pc_kpt(:) 
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds = .FALSE.
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Sfolding_vec(:) = 0.0_dp 
                enddo
            enddo
        enddo
    enddo
    GUR%n_pckpts = (GUR%n_pckpts)/nkpts
    deallocate(pc_kpt_already_folded,stat=alloc_stat)
    allocate(pc_kpt_already_folded(1:n_selec_pcbz_dirs,1:maxval(n_dirs_for_EBS_along_pcbz_dir(:)),1:maxval(n_pckpts_dirs(:))))
    pc_kpt_already_folded = .FALSE.
    !! Obtaining geometric unfolding relations
    do i_SCKPT=1,nkpts
        do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
            do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)
                do ipc_kpt=1, n_pckpts_dirs(i_selec_pcbz_dir)
                    if(pc_kpt_already_folded(i_selec_pcbz_dir,i_needed_dirs,ipc_kpt))then
                        cycle ! It has already folded before. Since a pckpt can only fold into one SCKPT, let's move on to the next pckpt.
                    endif
                    pc_kpt(:) = pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords(:)
                    current_SCKPT = list_of_SCKPTS(i_SCKPT)%coord(:)
                    do ieqv_SCKPT=1,SKPTS_eqv_to_SKPT(i_SCKPT) % neqv
                        SCKPT_eqv_to_current_SCKPT(:) = SKPTS_eqv_to_SKPT(i_SCKPT) % eqv_pt(ieqv_SCKPT) % coord(:)
                        trial_folding_G(:) = pc_kpt(:) - SCKPT_eqv_to_current_SCKPT(:)
                        if(vec_in_latt(vec=trial_folding_G, latt=B_matrix_SC,tolerance=default_tol_for_vec_equality))then
                            GUR%SCKPT_used_for_unfolding(i_SCKPT) = .TRUE.
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds = .TRUE.
                            GUR%n_folding_pckpts = GUR%n_folding_pckpts + 1
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords_actual_unfolding_K = &
                                SCKPT_eqv_to_current_SCKPT(:)
                            isym = SKPTS_eqv_to_SKPT(i_SCKPT) % eqv_pt(ieqv_SCKPT) % symop
                            pc_kpt =  pt_eqv_by_point_group_symop(point=pc_kpt,symops=crystal_SC%symops,isym=isym, fractional_coords=.FALSE.,invert_symop=.TRUE.)
                            origin_for_spin_proj =  pt_eqv_by_point_group_symop(point=args%origin_for_spin_proj_cartesian, &
                                                                                symops=crystal_SC%symops,isym=isym, &
                                                                                fractional_coords=.FALSE.,invert_symop=.TRUE.)
                            ! Message from Paulo:
                            ! The prefix "S" means "symmetrized". This is a little trick I came up with
                            ! that allows me to use the coefficients of a SC wavefunction psi(K',n) to
                            ! calculate the spectral weights associated with a SC wavefunction psi(K,n),
                            ! where K' = SK and S is a symmetry operation of the crystal's point group.
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:) = pc_kpt(:)
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Sorigin_for_spin_proj(:) = origin_for_spin_proj(:)
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Sfolding_vec(:) = &
                                pc_kpt(:) - current_SCKPT(:)
                            pc_kpt_already_folded(i_selec_pcbz_dir,i_needed_dirs,ipc_kpt) = .TRUE.
                            exit ! It has already folded. No need to keep looking for this pckpt.
                        else 
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds = .FALSE.
                        endif
                    enddo ! Loop over all the SCKPTS eqv. to the one on the WF file by symm. ops. of the SC
                enddo ! Loop over the pckpts along the dirs. of the pcbz that are eqv. to the selec. ones by 
                      ! symm. ops. of the pcbz and NOT by symm. ops. of the SCBZ
            enddo ! Loop over all dirs. of the pcbz that are eqv. to the selec. ones by symm. ops. of the pcbz and NOT by symm. ops. of the SCBZ
        enddo ! Loop over the selected pcbz directions
    enddo ! Loop over the SCKPTS found on the wavefunction file

    if(print_stuff)then
        write(*,"(A,/)")"    * Done."
    endif


end subroutine get_geom_unfolding_relations


subroutine define_pckpts_to_be_checked(pckpts_to_be_checked,all_dirs_used_for_EBS_along_pcbz_dir,nkpts_selected_dirs)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(selected_pcbz_directions), intent(out) :: pckpts_to_be_checked
type(irr_bz_directions), dimension(:), intent(in) :: all_dirs_used_for_EBS_along_pcbz_dir
integer, dimension(:), intent(in) :: nkpts_selected_dirs
integer :: n_selec_pcbz_dirs,i_selec_pcbz_dir,nk,ikpt,i_needed_dirs,n_dirs_for_EBS_along_pcbz_dir,alloc_stat
real(kind=dp), dimension(1:3) :: kstart,kend
type(vec3d), dimension(:), allocatable :: kline

    n_selec_pcbz_dirs = size(all_dirs_used_for_EBS_along_pcbz_dir)
    deallocate(pckpts_to_be_checked%selec_pcbz_dir,stat=alloc_stat)
    allocate(pckpts_to_be_checked%selec_pcbz_dir(1:n_selec_pcbz_dirs))
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        n_dirs_for_EBS_along_pcbz_dir = size(all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(:))
        deallocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir,stat=alloc_stat)
        allocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1:n_dirs_for_EBS_along_pcbz_dir))
        nk = nkpts_selected_dirs(i_selec_pcbz_dir)
        do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir
            deallocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt,stat=alloc_stat)
            allocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(1:nk),stat=alloc_stat)
            kstart = all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%kstart  
            kend = all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%kend
            deallocate(kline,stat=alloc_stat)
            allocate(kline(1:nk))
            kline = kpts_line(kstart,kend,nk)
            do ikpt=1,nk
                pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ikpt)%coords(:) = kline(ikpt)%coord(:)
            enddo
        enddo
    enddo


end subroutine define_pckpts_to_be_checked


subroutine select_coeffs_to_calc_spectral_weights(selected_coeff_indices, iall_G, crystal_pc, &
                                                  crystal_SC, GUR)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
integer, dimension(:), allocatable, intent(out) :: selected_coeff_indices
integer, dimension(:,:), intent(in) :: iall_G !Coordinates of the G points in the [b1, b2, b3] basis.
type(crystal_3D), intent(in) :: crystal_pc, crystal_SC
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR

integer :: alloc_stat, ig, nplane, i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt
real(kind=dp), dimension(1:3) :: SC_G, trial_pc_g, g, folding_G, &
                                 current_SCKPT, symmetrized_unf_pc_kpt
real(kind=dp), dimension(1:3,1:3) :: b_matrix_pc, B_matrix_SC
real(kind=dp) :: tol, tol_for_vec_equality


    b_matrix_pc = crystal_pc%rec_latt_vecs
    B_matrix_SC = crystal_SC%rec_latt_vecs

    i_SCKPT = GUR%current_index%i_SCKPT
    i_selec_pcbz_dir = GUR%current_index%i_selec_pcbz_dir
    i_needed_dirs = GUR%current_index%i_needed_dirs
    ipc_kpt = GUR%current_index%ipc_kpt

    current_SCKPT = GUR%list_of_SCKPTS(i_SCKPT)%coord(:)
    symmetrized_unf_pc_kpt(:) = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)% &
                                    needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:)
    folding_G(:) = symmetrized_unf_pc_kpt(:) - current_SCKPT(:)

    tol_for_vec_equality = default_tol_for_vec_equality
    tol = 0.1_dp * tol_for_vec_equality
    nplane = size(iall_G, dim=2)
    deallocate(selected_coeff_indices, stat=alloc_stat)
    do while((.not. allocated(selected_coeff_indices)) .and. &
              (tol <= max_tol_for_vec_equality))
        tol = 10_dp * tol
        !$omp parallel do &
        !$    schedule(guided) default(none) &
        !$    shared(nplane,iall_G,B_matrix_SC,folding_G,b_matrix_pc,tol,selected_coeff_indices) &
        !$    private(ig,g,SC_G,trial_pc_g) 
        do ig=1, nplane
            SC_G = iall_G(1,ig)*B_matrix_SC(1,:) + iall_G(2,ig)*B_matrix_SC(2,:) + iall_G(3,ig)*B_matrix_SC(3,:)
            trial_pc_g = SC_G(:) - folding_G(:)
            if(vec_in_latt(vec=trial_pc_g, latt=b_matrix_pc, tolerance=tol))then
                 !$omp critical
                 call append(item=ig,list=selected_coeff_indices) ! pc_g + folding_G has to belong to the 
                 !$omp end critical
             endif
        enddo
    enddo
    if(abs(tol - tol_for_vec_equality) > epsilon(1.0_dp) .and. allocated(selected_coeff_indices))then
        write(*,'(A)') &
        '    WARNING (select_coeffs_to_calc_spectral_weights): Problems selecting the coeffs to calculate the spec. weights for the current pc-kpt.'
        write(*,'(2(A,f0.6),A)') &
        '            The tolerace for testing vector equality had to be increased from ',tol_for_vec_equality,' to ',tol,'.'
        write(*,'(A)')'            The results might not be entirely correct.'
    endif

end subroutine select_coeffs_to_calc_spectral_weights


function spectral_weight_for_coeff(coeff, selected_coeff_indices, add_elapsed_time_to) result(spectral_weight)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
complex(kind=kind_cplx_coeffs), dimension(:,:), intent(in) :: coeff !Coefficients, coeff(ikpt, iband)
real(kind=dp), dimension(1:size(coeff,dim=2)) :: spectral_weight
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: n_bands_SC_calculation, iband, i, icoeff
real(kind=dp) :: stime, ftime

    stime = time()

    n_bands_SC_calculation = size(coeff,dim=2)
    spectral_weight(:) = 0.0_dp
    !$omp parallel do &
    !$    schedule(guided) default(none) &
    !$    private(iband,i,icoeff) &
    !$    shared(n_bands_SC_calculation,selected_coeff_indices,spectral_weight,coeff)
    do iband = 1, n_bands_SC_calculation
        do i=1, size(selected_coeff_indices)
            icoeff = selected_coeff_indices(i)
            spectral_weight(iband) = spectral_weight(iband) + abs(coeff(icoeff,iband))**2.0_dp
        enddo
    enddo

    if(present(add_elapsed_time_to))then
        ftime = time()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end function spectral_weight_for_coeff


subroutine calc_spectral_function(SF_at_pckpt,energies,SC_calc_ener,spectral_weight,std_dev,add_elapsed_time_to)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: SF_at_pckpt
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight, &
                                           SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in), optional :: std_dev
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, stime, ftime, sigma

    stime = time()

    sigma = 0.025_dp
    if(present(std_dev)) sigma = std_dev
    ! Calculating the spectral function A(pc_kpt,E) 
    n_SC_bands = size(spectral_weight)
    deallocate(SF_at_pckpt, stat=alloc_stat)
    allocate(SF_at_pckpt(1:size(energies)))
    SF_at_pckpt(:) = 0.0_dp
    !$omp parallel do &
    !$    schedule(guided) default(none) &
    !$    private(iener,E,i_SC_band,E_SC_band)  &
    !$    shared(energies,n_SC_bands,SC_calc_ener,SF_at_pckpt,spectral_weight,sigma)
    do iener=1,size(energies)  ! For E in the list of energy values for the E(k) plot
        E = energies(iener)
        do i_SC_band=1,n_SC_bands
            E_SC_band = SC_calc_ener(i_SC_band)
            SF_at_pckpt(iener) = SF_at_pckpt(iener) + &
                                 spectral_weight(i_SC_band)*delta(E_SC_band - E, std_dev=sigma)
        enddo
    enddo

    if(present(add_elapsed_time_to))then
        ftime = time()
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


subroutine calc_delta_N_pckpt(delta_N_pckpt, energies, SC_calc_ener, spectral_weight, std_dev)
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
    !$    schedule(guided) default(none) &
    !$    private(iener,E,i_SC_band,E_SC_band)  &
    !$    shared(energies,n_SC_bands,SC_calc_ener,delta_N_pckpt,spectral_weight, delta_e,sigma)
    do iener=1,size(energies)  ! For E in the list of energy values for the E(k) plot
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


subroutine get_delta_Ns_for_EBS(delta_N_pckpt, energy_grid, SC_calc_ener, spectral_weight, &
                                smearing_for_delta_function, add_elapsed_time_to)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: delta_N_pckpt
real(kind=dp), dimension(:), intent(in) :: energy_grid, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in), optional :: smearing_for_delta_function
real(kind=dp) :: sigma, stime, ftime
real(kind=dp), intent(inout), optional :: add_elapsed_time_to

    stime = time()
    sigma = epsilon(1.0_dp)
    if(present(smearing_for_delta_function))then
        sigma = max(sigma,abs(smearing_for_delta_function))
    endif
    call calc_delta_N_pckpt(delta_N_pckpt,energy_grid,SC_calc_ener,spectral_weight,sigma)
    if(present(add_elapsed_time_to))then
        ftime = time()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine get_delta_Ns_for_EBS


subroutine get_delta_Ns_for_output(delta_N_only_selected_dirs,delta_N_symm_avrgd_for_EBS, &
                                   delta_N,all_dirs_used_for_EBS_along_pcbz_dir,pckpts_to_be_checked)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(UnfoldedQuantitiesForOutput), intent(out) :: delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
type(UnfoldedQuantities), intent(in) :: delta_N
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(irr_bz_directions), dimension(:), intent(in) :: all_dirs_used_for_EBS_along_pcbz_dir
integer :: nener, i_selec_pcbz_dir, ipc_kpt, i_needed_dirs
real(kind=dp), dimension(:), allocatable :: avrgd_dNs, avrgd_spin_proj, avrgd_parallel_proj
real(kind=dp), dimension(:,:), allocatable :: avrgd_sigma
real(kind=dp) :: weight
logical :: output_spin_info


    nener = size(delta_N%selec_pcbz_dir(1)%needed_dir(1)%pckpt(1)%dN(:))
    allocate(avrgd_dNs(1:nener))
    call allocate_UnfoldedQuantitiesForOutput(delta_N_only_selected_dirs, pckpts_to_be_checked)
    call allocate_UnfoldedQuantitiesForOutput(delta_N_symm_avrgd_for_EBS, pckpts_to_be_checked)
    do i_selec_pcbz_dir=1,size(delta_N%selec_pcbz_dir(:))
        delta_N_only_selected_dirs%pcbz_dir(i_selec_pcbz_dir) = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)
        do ipc_kpt=1, size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)%pckpt(:))
            avrgd_dNs(:) = 0.0_dp
            allocate(delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%dN(1:nener))
            do i_needed_dirs=1,size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(:))
                weight = all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%weight
                avrgd_dNs = avrgd_dNs + weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(:)
            enddo
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%dN(:) = avrgd_dNs(:)
       enddo     
    enddo

    output_spin_info = allocated(delta_N%selec_pcbz_dir(1)%needed_dir(1)%pckpt(1)%spin_proj_perp)
    if(.not. output_spin_info) return
    allocate(avrgd_spin_proj(1:nener))
    allocate(avrgd_parallel_proj(1:nener))
    allocate(avrgd_sigma(1:nener,1:3))
    do i_selec_pcbz_dir=1,size(delta_N%selec_pcbz_dir(:))
        delta_N_only_selected_dirs%pcbz_dir(i_selec_pcbz_dir) = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)
        do ipc_kpt=1, size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)%pckpt(:))
            avrgd_spin_proj(:) = 0.0_dp
            avrgd_parallel_proj(:) = 0.0_dp
            avrgd_sigma(:,:) = 0.0_dp
            allocate(delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_perp(1:nener))
            allocate(delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_para(1:nener))
            allocate(delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%sigma(1:nener,1:3))
            do i_needed_dirs=1,size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(:))
                weight = all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%weight
                avrgd_spin_proj = avrgd_spin_proj + weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp(:)
                avrgd_parallel_proj = avrgd_parallel_proj + weight * &
                                      delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para(:)
                avrgd_sigma = avrgd_sigma + weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%sigma(:,:)
            enddo
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_perp(:) = avrgd_spin_proj(:)
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_para(:) = avrgd_parallel_proj(:)
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%sigma(:,:) = avrgd_sigma(:,:)
       enddo     
    enddo
  

end subroutine get_delta_Ns_for_output



subroutine calc_rho(rho, delta_N, pc_ener, SC_ener_bands, delta_e, coeffs, &
                                      selected_coeff_indices, std_dev, add_elapsed_time_to)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! Please make sure that delta_N /= 0 before calling this subroutine
implicit none
complex(kind=kind_cplx_coeffs), dimension(:,:), allocatable, intent(out) :: rho
real(kind=dp), intent(in) :: delta_N, pc_ener, delta_e
real(kind=dp), intent(in), dimension(:) :: SC_ener_bands
complex(kind=kind_cplx_coeffs), dimension(:,:,:), intent(in) :: coeffs
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), intent(in), optional :: std_dev
real(kind=dp) :: min_prefactor, lamb_ee1, lamb_ee2, prefactor, sigma, stime, ftime
integer :: m1, m2, n_SC_bands, pw_index, alloc_stat, n_coeffs, n_coeffs_pf, & ! m1 and m2 are SC band indices
           n_spinor, i_spinor, pf_pw_index, n_picked_bands
complex(kind=kind_cplx_coeffs), dimension(:), allocatable :: pf_coeffs_m1, pf_coeffs_m2 ! pf = partial function
real(kind=dp), intent(inout), optional :: add_elapsed_time_to


    stime = time()
    if(delta_N < epsilon(1.0_dp))then
        write(*,'(A)')'ERROR (calc_rho): delta_N = 0.'
        write(*,'(A)')'Stopping now.'
        stop
    endif

    min_prefactor = 1E-5_dp
    sigma = epsilon(1.0_dp)
    if(present(std_dev)) sigma = std_dev

    n_spinor = size(coeffs, dim=1)
    n_coeffs = size(coeffs, dim=2)
    n_SC_bands = size(coeffs, dim=3)
    n_coeffs_pf = size(selected_coeff_indices)
    deallocate(rho, stat=alloc_stat)
    allocate(rho(1:n_SC_bands, 1:n_SC_bands))

    rho(:,:) = 0.0_dp
    if(args.perform_unfold)then
        !$omp parallel do &
        !$    schedule(guided) default(none) &
        !$    shared(n_SC_bands, pc_ener, SC_ener_bands, delta_e, sigma, min_prefactor, delta_N, &
        !$           n_coeffs_pf, n_spinor, coeffs, selected_coeff_indices, rho) &
        !$    private(m1, m2, lamb_ee1, lamb_ee2, prefactor, pf_coeffs_m1, pf_coeffs_m2, alloc_stat, &
        !$            i_spinor, pf_pw_index)
        do m1=1, n_SC_bands - 1
            lamb_ee1 = lambda(pc_ener=pc_ener, SC_ener=SC_ener_bands(m1), delta_e=delta_e, std_dev=sigma)
            if(lamb_ee1 < min_prefactor) cycle
            do m2=m1, n_SC_bands
                lamb_ee2 = lambda(pc_ener=pc_ener, SC_ener=SC_ener_bands(m2), delta_e=delta_e, std_dev=sigma)
                prefactor = (lamb_ee1 * lamb_ee2) / delta_N
                if(prefactor < min_prefactor) cycle
            
                allocate(pf_coeffs_m1(1:n_coeffs_pf), pf_coeffs_m2(1:n_coeffs_pf))
                do i_spinor=1, n_spinor
                    pf_pw_index = 0
                    do pw_index = 1, n_coeffs_pf
                        pf_pw_index = pf_pw_index + 1
                        pf_coeffs_m1(pf_pw_index) = coeffs(i_spinor, selected_coeff_indices(pw_index), m1)
                        pf_coeffs_m2(pf_pw_index) = coeffs(i_spinor, selected_coeff_indices(pw_index), m2)
                    enddo
                    rho(m1,m2) = rho(m1,m2) + dot_product(pf_coeffs_m1, pf_coeffs_m2)
                enddo
                deallocate(pf_coeffs_m1, pf_coeffs_m2, stat=alloc_stat)
                rho(m1,m2) = prefactor * rho(m1,m2)
                rho(m2,m1) = conjg(rho(m1,m2))
            enddo
        enddo
    else
        n_picked_bands = 0
        do m1=1, n_SC_bands
            lamb_ee1 = lambda(pc_ener=pc_ener, SC_ener=SC_ener_bands(m1), delta_e=delta_e, std_dev=sigma)
            if(lamb_ee1 < min_prefactor) cycle
            n_picked_bands = n_picked_bands + 1
            do m2=1, n_SC_bands
                lamb_ee2 = lambda(pc_ener=pc_ener, SC_ener=SC_ener_bands(m2), delta_e=delta_e, std_dev=sigma)
                prefactor = (lamb_ee1 * lamb_ee2) / delta_N
                if(prefactor < min_prefactor) cycle
                rho(m1,m2) = 1.0_dp
            enddo
        enddo
        if(n_picked_bands > 0) rho = rho / real(n_picked_bands, kind=dp)
    endif

    if(present(add_elapsed_time_to))then
        ftime = time()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine calc_rho


subroutine get_SC_pauli_matrix_elmnts(SC_pauli_matrix_elmnts, pauli_mtx_elmts_already_calc, &
                                      coeff, SC_band_energies, selected_pc_ener, delta_e, &
                                      add_elapsed_time_to)
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
complex(kind=kind_cplx_coeffs), dimension(:,:,:), allocatable, intent(inout) :: SC_pauli_matrix_elmnts
logical, dimension(:,:), allocatable, intent(inout) :: pauli_mtx_elmts_already_calc
real(kind=dp), dimension(:), intent(in) :: SC_band_energies
real(kind=dp), intent(in) :: selected_pc_ener, delta_e
complex(kind=kind_cplx_coeffs), dimension(:,:,:), intent(in) :: coeff ! coefficients(i_spinor_comp, iplane, iband)
complex(kind=kind_cplx_coeffs) :: inner_prod_Km2Alpha_Km1Beta, inner_prod_Km2Beta_Km1Alpha, &
                                  inner_prod_Km2Alpha_Km1Alpha, inner_prod_Km2Beta_Km1Beta
complex(kind=kind_cplx_coeffs), parameter :: J = (0.0_kind_cplx_coeffs, 1.0_kind_cplx_coeffs), &
                                             cplx_zero = (0.0_kind_cplx_coeffs, 0.0_kind_cplx_coeffs)
integer :: n_SC_bands, alloc_stat, m1, m2, m3, m4
integer, dimension(:), allocatable :: picked_SC_bands ! list of SC bands to enter the calc of the matrix elements
real(kind=dp) :: lambd, stime, ftime
real(kind=dp), intent(inout), optional :: add_elapsed_time_to


    stime = time()
    n_SC_bands = size(coeff, dim=3)
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
        lambd = lambda(pc_ener=selected_pc_ener, SC_ener=SC_band_energies(m1), delta_e=delta_e)
        if(lambd > 1E-1)then
            call append(list=picked_SC_bands, item=m1)
        endif
    enddo
    if(allocated(picked_SC_bands))then
        !$omp parallel do collapse(2) &
        !$    schedule(guided) default(none) &
        !$    shared(coeff, SC_pauli_matrix_elmnts, picked_SC_bands, pauli_mtx_elmts_already_calc) &
        !$    private(m1, m2, m3, m4, inner_prod_Km2Alpha_Km1Beta, inner_prod_Km2Beta_Km1Alpha, &
        !$            inner_prod_Km2Alpha_Km1Alpha, inner_prod_Km2Beta_Km1Beta) 
        do m3=1,size(picked_SC_bands)
            do m4=1,size(picked_SC_bands)
                m1 = picked_SC_bands(m3)
                m2 = picked_SC_bands(m4)
                if(pauli_mtx_elmts_already_calc(m2,m1)) cycle
                inner_prod_Km2Alpha_Km1Beta = dot_product(coeff(1,:,m2), coeff(2,:,m1))
                inner_prod_Km2Beta_Km1Alpha = dot_product(coeff(2,:,m2), coeff(1,:,m1))
                inner_prod_Km2Alpha_Km1Alpha = dot_product(coeff(1,:,m2), coeff(1,:,m1))
                inner_prod_Km2Beta_Km1Beta = dot_product(coeff(2,:,m2), coeff(2,:,m1))

                SC_pauli_matrix_elmnts(1,m2,m1) = inner_prod_Km2Alpha_Km1Beta + inner_prod_Km2Beta_Km1Alpha
                SC_pauli_matrix_elmnts(2,m2,m1) = J * (inner_prod_Km2Beta_Km1Alpha - inner_prod_Km2Alpha_Km1Beta) 
                SC_pauli_matrix_elmnts(3,m2,m1) = inner_prod_Km2Alpha_Km1Alpha - inner_prod_Km2Beta_Km1Beta
                pauli_mtx_elmts_already_calc(m2,m1) = .TRUE.

                SC_pauli_matrix_elmnts(:,m1,m2) = conjg(SC_pauli_matrix_elmnts(:,m2,m1))
                pauli_mtx_elmts_already_calc(m1,m2) = .TRUE.
            enddo
        enddo
        deallocate(picked_SC_bands)
    endif

    if(present(add_elapsed_time_to))then
        ftime = time()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif


end subroutine get_SC_pauli_matrix_elmnts


subroutine calc_spin_projections(spin_proj_perp, spin_proj_para, pauli_vector, cart_coords_kpt, origin)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! This is a BETA feature, and should be used with caution
implicit none
real(kind=dp), intent(out) :: spin_proj_perp, spin_proj_para
real(kind=dp), dimension(1:3), intent(in) :: pauli_vector, cart_coords_kpt
real(kind=dp), dimension(1:3), intent(in), optional :: origin ! In cartesian coordinates
real(kind=dp), dimension(1:3) :: versor_para, versor_perp, normal2proj_plane, cart_coords_origin

    cart_coords_origin = 0.0_dp
    if(present(origin)) cart_coords_origin = origin
    normal2proj_plane = versor(args%normal_to_proj_plane)
    versor_para = versor(cart_coords_kpt - cart_coords_origin)
    versor_perp = versor(cross(normal2proj_plane, versor_para)) 
    spin_proj_perp = dot_product(pauli_vector, versor_perp)
    spin_proj_para = dot_product(pauli_vector, versor_para)

end subroutine calc_spin_projections


subroutine perform_unfolding(delta_N, times, GUR, coefficients, selected_coeff_indices, &
                             energy_grid, ener_SC_bands)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none 
type(UnfoldedQuantities), intent(inout) :: delta_N
type(timekeeping), intent(inout) :: times
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR
complex(kind=kind_cplx_coeffs), dimension(:,:,:), intent(in) :: coefficients ! Coefficients, coefficients(i_spinor_comp, ikpt, iband)
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), dimension(:), intent(in) :: energy_grid, ener_SC_bands

integer :: n_WF_components, n_bands_SC_calculation, i_spinor, iener, nener, iener2, alloc_stat, i, &
           i_selec_pcbz_dir, i_needed_dirs, ipc_kpt, i_SCKPT
real(kind=dp), dimension(:), allocatable :: spectral_weight
real(kind=dp) :: selected_pc_ener, delta_e, dN, unf_spin_proj_perp, unf_spin_proj_para
real(kind=dp), dimension(1:3) :: symmetrized_unf_pc_kpt, unf_pauli_vec, projection_origin
complex(kind=kind_cplx_coeffs), dimension(:,:), allocatable :: rho
complex(kind=kind_cplx_coeffs), dimension(:,:,:), allocatable :: SC_pauli_matrix_elmnts
integer, dimension(:), allocatable :: pc_energies_to_calc_eigenvalues
logical, dimension(:,:), allocatable :: pauli_mtx_elmts_already_calc


    i_SCKPT = GUR%current_index%i_SCKPT
    i_selec_pcbz_dir = GUR%current_index%i_selec_pcbz_dir
    i_needed_dirs = GUR%current_index%i_needed_dirs
    ipc_kpt = GUR%current_index%ipc_kpt
    ! In cartesion coors!
    symmetrized_unf_pc_kpt(:) = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)% &
                                    needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:)

    n_WF_components = size(coefficients, dim=1) ! 1 for a regular WF, and 2 for a spinor WF
    n_bands_SC_calculation = size(coefficients, dim=3)
    allocate(spectral_weight(1:n_bands_SC_calculation))

    if(args.perform_unfold)then
        spectral_weight = 0.0_dp
        do i_spinor=1,n_WF_components
           ! Calculating spectral_weights
            spectral_weight = spectral_weight + spectral_weight_for_coeff(coefficients(i_spinor,:,:), &
                                                                          selected_coeff_indices, &
                                                                          add_elapsed_time_to=times%calc_spec_weights)
        enddo
    else
        write(*,'(A)')'    **** WARNING: You are using the option "-dont_unfold". You are NOT performing unfolding!! ****'
        spectral_weight = 1.0_dp
    endif
    
    ! Calculating the delta_Ns
    call get_delta_Ns_for_EBS(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN, &
                              energy_grid, ener_SC_bands, spectral_weight, &
                              add_elapsed_time_to=times%calc_dN)

    if(calc_spec_func_explicitly)then
        ! Calculating the spectral_function (optional)
        call calc_spectral_function(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%SF, &
                                    energy_grid, ener_SC_bands, spectral_weight, &
                                    add_elapsed_time_to=times%calc_SF)
    endif

    if(n_WF_components==2)then
        nener = size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN)
        allocate(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp(1:nener))
        allocate(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para(1:nener))
        allocate(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%sigma(1:nener,1:3))
        delta_e = energy_grid(2) - energy_grid(1) ! uniform energy grid

        !$omp parallel do &
        !$    schedule(guided) default(none) &
        !$    shared(nener, delta_N, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt, &
        !$           pc_energies_to_calc_eigenvalues) &
        !$    private(iener, dN)
        do iener=1, nener
            dN = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(iener)
            if(dN < 1E-3_dp)then ! Can hardly be considered a band; no point calculating rho
                delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp(iener) = 0.0_dp 
                delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para(iener) = 0.0_dp 
                delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%sigma(iener,:) = 0.0_dp
            else
                !$omp critical
                call append(list=pc_energies_to_calc_eigenvalues, item=iener)
                !$omp end critical
            endif
        enddo

        deallocate(pauli_mtx_elmts_already_calc, SC_pauli_matrix_elmnts, stat=alloc_stat)
        do iener2=1, size(pc_energies_to_calc_eigenvalues)
            iener = pc_energies_to_calc_eigenvalues(iener2)
            selected_pc_ener = energy_grid(iener)
            dN = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(iener)
            call calc_rho(rho, dN, energy_grid(iener), ener_SC_bands, delta_e, &
                                            coefficients, selected_coeff_indices, &
                                            add_elapsed_time_to=times%calc_rho)

            call get_SC_pauli_matrix_elmnts(SC_pauli_matrix_elmnts, pauli_mtx_elmts_already_calc, &
                                            coefficients, ener_SC_bands, selected_pc_ener, delta_e, &
                                            add_elapsed_time_to=times%calc_pauli_vec)
            do i=1,3
                unf_pauli_vec(i) = real(trace_AB(A=rho, B=SC_pauli_matrix_elmnts(i,:,:)), kind=dp)
            enddo
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%sigma(iener,:) = unf_pauli_vec(:)

            projection_origin = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Sorigin_for_spin_proj(:)
            call calc_spin_projections(unf_spin_proj_perp, unf_spin_proj_para, unf_pauli_vec, &
                                       symmetrized_unf_pc_kpt, origin=projection_origin)
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp(iener) = unf_spin_proj_perp 
            delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para(iener) = unf_spin_proj_para
            deallocate(rho, stat=alloc_stat) !! rho is allocated in the calc_rho routine
        enddo
        deallocate(pauli_mtx_elmts_already_calc, SC_pauli_matrix_elmnts, stat=alloc_stat)
    endif

end subroutine perform_unfolding


end module band_unfolding


