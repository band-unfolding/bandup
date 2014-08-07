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

module band_unfolding
!$ use omp_lib
use math
use general_io
implicit none
PRIVATE
PUBLIC :: get_geom_unfolding_relations, define_pckpts_to_be_checked, &
          select_coeffs_to_calc_spectral_weights, get_delta_Ns_for_output, &
          perform_unfolding

CONTAINS

subroutine get_geom_unfolding_relations(GUR,list_of_SCKPTS,pckpts_to_be_checked,crystal_SC,verbose)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(geom_unfolding_relations_for_each_SCKPT), intent(out) :: GUR !! Geometric Unfolding Relations
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(vec3d), dimension(:), intent(in) :: list_of_SCKPTS
type(crystal_3D), intent(in) :: crystal_SC
logical, intent(in), optional :: verbose
integer :: nkpts, n_selec_pcbz_dirs, i_SCKPT,ipc_kpt,ieqv_SCKPT,isym, &
           i_selec_pcbz_dir,i_needed_dirs,alloc_stat
integer, dimension(:), allocatable :: n_dirs_for_EBS_along_pcbz_dir,n_pckpts_dirs
logical, dimension(:,:,:), allocatable :: pc_kpt_already_folded
real(kind=dp), dimension(1:3) :: pc_kpt, current_SCKPT,SCKPT_eqv_to_current_SCKPT,trial_folding_G
real(kind=dp), dimension(1:3,1:3) :: B_matrix_SC
type(star), dimension(:), allocatable :: SKPTS_eqv_to_SKPT
type(symmetry_operation), dimension(:), allocatable :: symops
logical :: print_stuff

    print_stuff = .FALSE.
    if(present(verbose))then
        print_stuff = verbose
    endif

    if(print_stuff)then
        write(*,"(A)",advance="no")"Verifying geometric unfolding relations between pcbz and SCBZ wave-vectors... "
    endif
    B_matrix_SC = crystal_SC%rec_latt_vecs
    call get_star(star_of_pt=SKPTS_eqv_to_SKPT, points=list_of_SCKPTS, crystal=crystal_SC, symm_ops=symops, &
                  tol_for_vec_equality=default_tol_for_vec_equality, &
                  symprec=default_symprec, reduce_to_bz=.TRUE., &
                  try_to_reduce_to_a_pc=.FALSE.) ! fails if try_to_reduce_to_a_pc=.TRUE.
    !! Allocating and initializing table
    nkpts = size(list_of_SCKPTS)
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
                            pc_kpt =  pt_eqv_by_point_group_symop(point=pc_kpt,symops=symops,isym=isym, fractional_coords=.FALSE.,invert_symop=.TRUE.)
                            ! Message from Paulo:
                            ! The prefix "S" means "symmetrized". This is a little trick I came up with
                            ! that allows me to use the coefficients of a SC wavefunction psi(K',n) to
                            ! calculate the spectral weights associated with a SC wavefunction psi(K,n),
                            ! where K' = SK and S is a symmetry operation of the crystal's point group.
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:) = pc_kpt(:)
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
        write(*,"(A,/)")"Done."
    endif


end subroutine get_geom_unfolding_relations


subroutine define_pckpts_to_be_checked(pckpts_to_be_checked,dirs_req_for_symmavgd_EBS_along_pcbz_dir,nkpts_selected_dirs)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(selected_pcbz_directions), intent(out) :: pckpts_to_be_checked
type(irr_bz_directions), dimension(:), intent(in) :: dirs_req_for_symmavgd_EBS_along_pcbz_dir
integer, dimension(:), intent(in) :: nkpts_selected_dirs
integer :: n_selec_pcbz_dirs,i_selec_pcbz_dir,nk,ikpt,i_needed_dirs,n_dirs_for_EBS_along_pcbz_dir,alloc_stat
real(kind=dp), dimension(1:3) :: kstart,kend
type(vec3d), dimension(:), allocatable :: kline

    n_selec_pcbz_dirs = size(dirs_req_for_symmavgd_EBS_along_pcbz_dir)
    deallocate(pckpts_to_be_checked%selec_pcbz_dir,stat=alloc_stat)
    allocate(pckpts_to_be_checked%selec_pcbz_dir(1:n_selec_pcbz_dirs))
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        n_dirs_for_EBS_along_pcbz_dir = size(dirs_req_for_symmavgd_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(:))
        deallocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir,stat=alloc_stat)
        allocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1:n_dirs_for_EBS_along_pcbz_dir))
        nk = nkpts_selected_dirs(i_selec_pcbz_dir)
        do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir
            deallocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt,stat=alloc_stat)
            allocate(pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(1:nk),stat=alloc_stat)
            kstart = dirs_req_for_symmavgd_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%kstart  
            kend = dirs_req_for_symmavgd_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%kend
            deallocate(kline,stat=alloc_stat)
            allocate(kline(1:nk))
            kline = kpts_line(kstart,kend,nk)
            do ikpt=1,nk
                pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ikpt)%coords(:) = kline(ikpt)%coord(:)
            enddo
        enddo
    enddo


end subroutine define_pckpts_to_be_checked


subroutine select_coeffs_to_calc_spectral_weights(selected_coeff_indices, &
                                                  iall_G,b_matrix_pc,B_matrix_SC, & 
                                                  folding_G, &
                                                  tol_for_vec_equality)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
integer, dimension(:), allocatable, intent(out) :: selected_coeff_indices
integer, dimension(:,:), intent(in) :: iall_G !Coordinates of the G points in the [b1, b2, b3] basis.
real(kind=dp), dimension(1:3,1:3), intent(in) :: b_matrix_pc,B_matrix_SC
real(kind=dp), dimension(1:3), intent(in) :: folding_G
real(kind=dp), intent(in) :: tol_for_vec_equality
integer :: alloc_stat, ig, nplane
real(kind=dp), dimension(1:3) :: SC_G, trial_pc_g,g
real(kind=dp) :: tol

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
complex(kind=sp), dimension(:,:), intent(in) :: coeff !Coefficients, coeff(ikpt, iband)
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



subroutine calc_unfolded_delta_Q_pckpt(delta_Q_pckpt, Q_SC, energies,SC_calc_ener,spectral_weight,sigma)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
! Q = Q[current_SCKPT](1:n_SC_bands) is any property calculated, using the SC,
! at the current SC-KPT and for all SC bands.
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: delta_Q_pckpt
! Make sure size(Q_SC) == n_SC_bands, because I don't check this!!
real(kind=dp), dimension(:), intent(in) :: Q_SC, energies, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in), optional :: sigma
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, delta_e, std_dev

    
    std_dev = epsilon(1D0)
    if(present(sigma)) std_dev = sigma
    delta_e = energies(2)  - energies(1) !! Uniform energy grid 
    n_SC_bands = size(spectral_weight)
    deallocate(delta_Q_pckpt, stat=alloc_stat)
    allocate(delta_Q_pckpt(1:size(energies)))
    delta_Q_pckpt(:) = 0.0_dp
    !$omp parallel do &
    !$    schedule(guided) default(none) &
    !$    private(iener,E,i_SC_band,E_SC_band)  &
    !$    shared(energies,n_SC_bands,SC_calc_ener,delta_Q_pckpt,spectral_weight, Q_SC, delta_e,std_dev)
    do iener=1,size(energies)  ! For E in the list of energy values for the E(k) plot
        E = energies(iener)
        do i_SC_band=1,n_SC_bands
            E_SC_band = SC_calc_ener(i_SC_band)
            delta_Q_pckpt(iener) = delta_Q_pckpt(iener) + &
                                      spectral_weight(i_SC_band) * Q_SC(i_SC_band) * &
                                      integral_delta_x_minus_x0(x0=E_SC_band, &
                                                                lower_lim=E-0.5_dp*delta_e, &
                                                                upper_lim=E+0.5_dp*delta_e, &
                                                                std_dev=std_dev)
        enddo
    enddo

end subroutine calc_unfolded_delta_Q_pckpt



subroutine calc_delta_N_pckpt(delta_N_pckpt,energies,SC_calc_ener,spectral_weight,sigma)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: delta_N_pckpt
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in) :: sigma
integer :: alloc_stat
real(kind=dp), dimension(:), allocatable :: N_SC

    deallocate(N_SC, stat=alloc_stat)
    allocate(N_SC(1:size(SC_calc_ener)))
    N_SC = 1.0_dp ! By construction
    call calc_unfolded_delta_Q_pckpt(delta_N_pckpt, N_SC, energies, SC_calc_ener, spectral_weight, sigma)


end subroutine calc_delta_N_pckpt


subroutine calc_SC_spin_info_from_spinor_WF(expec_val_pauli_M, spin_proj_perp, spin_proj_para, &
                                            cart_coords_kpt, coeff)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! This is a BETA feature, and should be used with caution
implicit none
real(kind=dp), dimension(:,:), allocatable, intent(out) :: expec_val_pauli_M
real(kind=dp), dimension(:), allocatable, intent(out) :: spin_proj_perp, spin_proj_para
real(kind=dp), dimension(1:3), intent(in) :: cart_coords_kpt
complex(kind=sp), dimension(:,:,:), intent(in) :: coeff ! coefficients(i_spinor_comp, ikpt, iband)
real(kind=dp), dimension(1:3) :: quantization_axis, versor_k, versor_normal_to_k
complex(kind=sp) :: inner_prod
integer :: iband, nbands

    quantization_axis = saxis ! saxis is defined in the general_io and options_setup modules
    versor_k = versor(cart_coords_kpt)                                       
    versor_normal_to_k = versor(cross(quantization_axis, versor_k)) 
 
    nbands = size(coeff, dim=3)
    if(.not. allocated(expec_val_pauli_M)) allocate(expec_val_pauli_M(1:nbands, 1:3))
    if(.not. allocated(spin_proj_perp)) allocate(spin_proj_perp(1:nbands))
    if(.not. allocated(spin_proj_para)) allocate(spin_proj_para(1:nbands))
    do iband=1, nbands
        inner_prod = inner_product(coeff(1,:,iband), coeff(2,:,iband))
        !! expec_val_pauli_M = <psi_spinor| (sigma_x, sigma_y, sigma_z) |psi_spinor>
        expec_val_pauli_M(iband,1) = 2.0_dp * real(inner_prod)
        expec_val_pauli_M(iband,2) = 2.0_dp * aimag(inner_prod)
        expec_val_pauli_M(iband,3) = real(inner_product(coeff(1,:,iband), coeff(1,:,iband)) - & 
                                          inner_product(coeff(2,:,iband), coeff(2,:,iband)))    

        spin_proj_perp(iband) = dot_product(expec_val_pauli_M(iband,:), versor_normal_to_k)
        spin_proj_para(iband) = dot_product(expec_val_pauli_M(iband,:), versor_k)
    enddo

end subroutine calc_SC_spin_info_from_spinor_WF



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
                                   delta_N,dirs_req_for_symmavgd_EBS_along_pcbz_dir,pckpts_to_be_checked)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(UnfoldedQuantitiesForOutput), intent(out) :: delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
type(UnfoldedQuantities), intent(in) :: delta_N
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(irr_bz_directions), dimension(:), intent(in) :: dirs_req_for_symmavgd_EBS_along_pcbz_dir
integer :: nener, i_selec_pcbz_dir, ipc_kpt, i_needed_dirs
real(kind=dp), dimension(:), allocatable :: avrgd_dNs, avrgd_spin_proj, avrgd_parallel_proj
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
                weight = dirs_req_for_symmavgd_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%weight
                avrgd_dNs = avrgd_dNs + weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(:)
            enddo
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%dN(:) = avrgd_dNs(:)
       enddo     
    enddo

    output_spin_info = allocated(delta_N%selec_pcbz_dir(1)%needed_dir(1)%pckpt(1)%spin_proj_perp)
    if(.not. output_spin_info) return
    allocate(avrgd_spin_proj(1:nener))
    allocate(avrgd_parallel_proj(1:nener))
    do i_selec_pcbz_dir=1,size(delta_N%selec_pcbz_dir(:))
        delta_N_only_selected_dirs%pcbz_dir(i_selec_pcbz_dir) = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)
        do ipc_kpt=1, size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)%pckpt(:))
            avrgd_spin_proj(:) = 0.0_dp
            avrgd_parallel_proj(:) = 0.0_dp
            allocate(delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_perp(1:nener))
            allocate(delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_para(1:nener))
            do i_needed_dirs=1,size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(:))
                weight = dirs_req_for_symmavgd_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%weight
                avrgd_spin_proj = avrgd_spin_proj + weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp(:)
                avrgd_parallel_proj = avrgd_parallel_proj + weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para(:)
            enddo
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_perp(:) = avrgd_spin_proj(:)
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%spin_proj_para(:) = avrgd_parallel_proj(:)
       enddo     
    enddo
  

end subroutine get_delta_Ns_for_output


subroutine perform_unfolding(delta_N, time_calc_spectral_weights, time_spent_calculating_delta_Ns, time_calc_spectral_function, & 
                             i_selec_pcbz_dir, i_needed_dirs, ipc_kpt, cart_coords_pckpt, & 
                             coefficients, selected_coeff_indices, energy_grid, ener_SC_bands, calc_spec_func_explicitly)
implicit none 
type(UnfoldedQuantities), intent(inout) :: delta_N
real(kind=dp), intent(out) :: time_calc_spectral_weights, time_spent_calculating_delta_Ns, time_calc_spectral_function
integer, intent(in) :: i_selec_pcbz_dir, i_needed_dirs, ipc_kpt
real(kind=dp), dimension(1:3), intent(in) :: cart_coords_pckpt
complex(kind=sp), dimension(:,:,:), intent(in) :: coefficients ! Coefficients, coefficients(i_spinor_comp, ikpt, iband)
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), dimension(:), intent(in) :: energy_grid, ener_SC_bands
logical, intent(in) :: calc_spec_func_explicitly
integer :: n_WF_components, n_bands_SC_calculation, i_spinor
real(kind=dp), dimension(:), allocatable :: spectral_weight
real(kind=dp), dimension(:), allocatable :: spin_proj_perp, spin_proj_para, & 
                                            unf_spin_proj_perp, unf_spin_proj_para
real(kind=dp), dimension(:,:), allocatable :: expec_val_pauli_M

    n_WF_components = size(coefficients, dim=1) ! 1 for a regular WF, and 2 for a spinor WF
    n_bands_SC_calculation = size(coefficients, dim=3)
    allocate(spectral_weight(1:n_bands_SC_calculation))

    spectral_weight = 0.0_dp
    do i_spinor=1,n_WF_components
        ! Calculating spectral_weights
        spectral_weight = spectral_weight + spectral_weight_for_coeff(coefficients(i_spinor,:,:), &
                                                                      selected_coeff_indices, &
                                                                      add_elapsed_time_to=time_calc_spectral_weights)

    enddo
    ! Calculating the delta_Ns
    call get_delta_Ns_for_EBS(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN, &
                              energy_grid, ener_SC_bands, spectral_weight, &
                              add_elapsed_time_to=time_spent_calculating_delta_Ns)
    if(calc_spec_func_explicitly)then
        ! Calculating the spectral_function (optional)
        call calc_spectral_function(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%SF, &
                                    energy_grid, ener_SC_bands, spectral_weight, &
                                    add_elapsed_time_to=time_calc_spectral_function)
    endif

    if(n_WF_components==2)then
        call calc_SC_spin_info_from_spinor_WF(expec_val_pauli_M, spin_proj_perp, spin_proj_para, &
                                              cart_coords_pckpt, coefficients)
        call calc_unfolded_delta_Q_pckpt(unf_spin_proj_perp, spin_proj_perp, energy_grid, ener_SC_bands, spectral_weight)
        call calc_unfolded_delta_Q_pckpt(unf_spin_proj_para, spin_proj_para, energy_grid, ener_SC_bands, spectral_weight)
        allocate(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp(1:size(unf_spin_proj_perp)))
        allocate(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para(1:size(unf_spin_proj_para)))
        delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_perp = unf_spin_proj_perp 
        delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%spin_proj_para = unf_spin_proj_para
    endif

end subroutine perform_unfolding


end module band_unfolding


