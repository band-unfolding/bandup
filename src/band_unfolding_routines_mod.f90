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
implicit none
PRIVATE
PUBLIC :: get_geom_unfolding_relations, define_pckpts_to_be_checked, &
          select_coeffs_to_calc_spectral_weights,calc_spectral_weights, &
          calc_spectral_function, get_delta_Ns_for_EBS, get_delta_Ns_for_output

CONTAINS

subroutine get_geom_unfolding_relations(GUR,list_of_SCKPTS,pckpts_to_be_checked,B_matrix_SC,verbose)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
type(geom_unfolding_relations_for_each_SCKPT), intent(out) :: GUR !! Geometric Unfolding Relations
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(vec3d), dimension(:), intent(in) :: list_of_SCKPTS
real(kind=dp), dimension(1:3,1:3), intent(in) :: B_matrix_SC
logical, intent(in), optional :: verbose
integer :: nkpts, n_selec_pcbz_dirs, i_SCKPT,ipc_kpt,ieqv_SCKPT,nsym,isym, &
           i_selec_pcbz_dir,i_needed_dirs,alloc_stat
integer, dimension(:), allocatable :: n_dirs_for_EBS_along_pcbz_dir,n_pckpts_dirs
logical, dimension(:,:,:), allocatable :: pc_kpt_already_folded
real(kind=dp), dimension(1:3) :: pc_kpt, current_SCKPT,SCKPT_eqv_to_current_SCKPT,trial_folding_G
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
    call get_symm(nsym=nsym, symops=symops, symprec=1D-6, lattice=B_matrix_SC) 
    call get_star(star_of_pt=SKPTS_eqv_to_SKPT,points=list_of_SCKPTS,lattice=B_matrix_SC, &
                        tol_for_vec_equality=tol_for_vec_equality,reduce_to_bz=.TRUE.)
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
    allocate(GUR%SCKPT(1:nkpts))
    do i_SCKPT=1,nkpts
        allocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(1:n_selec_pcbz_dirs))
        do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
            allocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1:n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)))
            do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)
                allocate(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(1:n_pckpts_dirs(i_selec_pcbz_dir)))
                do ipc_kpt=1,n_pckpts_dirs(i_selec_pcbz_dir)
                    GUR%n_pckpts = GUR%n_pckpts + 1
                    pc_kpt(:) = pckpts_to_be_checked%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords(:)
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords(:) = pc_kpt(:) 
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:) = pc_kpt(:) 
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds = .FALSE.
                    GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Sfolding_vec(:) = 0.0d0 
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
                        if(vec_in_latt(vec=trial_folding_G, latt=B_matrix_SC,tolerance=tol_for_vec_equality))then
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds = .TRUE.
                            GUR%n_folding_pckpts = GUR%n_folding_pckpts + 1
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords_actual_unfolding_K = &
                                SCKPT_eqv_to_current_SCKPT(:)
                            isym = SKPTS_eqv_to_SKPT(i_SCKPT) % eqv_pt(ieqv_SCKPT) % symop
                            pc_kpt =  pt_eqv_by_point_group_symop(point=pc_kpt,symops=symops,isym=isym, fractional_coords=.FALSE.,invert_symop=.TRUE.)
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:) = pc_kpt(:)
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Sfolding_vec(:) = pc_kpt(:) - current_SCKPT(:)
                            pc_kpt_already_folded(i_selec_pcbz_dir,i_needed_dirs,ipc_kpt) = .TRUE.
                            exit ! It has already folded. No need to keep looking for this pckpt.
                        else 
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds = .FALSE.
                        endif
                    enddo ! Loop over all the SCKPTS eqv. to the one on the WF file by symm. ops. of the SC
                enddo ! Loop over the pckpts along the dirs. of the pcbz that are eqv. to the selec. ones by symm. ops. of the pcbz and NOT by symm. ops. of the SCBZ
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
                                                  iall_G,b_matrix_pc,B_matrix_SC,nplane, & 
                                                  folding_G, &
                                                  tol_for_vec_equality)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
integer, dimension(:), allocatable, intent(out) :: selected_coeff_indices
integer, dimension(:,:), intent(in) :: iall_G !Coordinates of the G points in the [b1, b2, b3] basis.
real(kind=dp), dimension(1:3,1:3), intent(in) :: b_matrix_pc,B_matrix_SC
integer, intent(in) :: nplane
real(kind=dp), dimension(1:3), intent(in) :: folding_G
real(kind=dp), intent(in) :: tol_for_vec_equality
integer :: alloc_stat, ig
real(kind=dp), dimension(1:3) :: SC_G, trial_pc_g,g

    deallocate(selected_coeff_indices, stat=alloc_stat)
    !$omp parallel do &
    !$    schedule(guided) default(none) &
    !$    shared(nplane,iall_G,B_matrix_SC,folding_G,b_matrix_pc,tol_for_vec_equality,selected_coeff_indices) &
    !$    private(ig,g,SC_G,trial_pc_g)
    do ig=1, nplane
        SC_G = iall_G(1,ig)*B_matrix_SC(1,:) + iall_G(2,ig)*B_matrix_SC(2,:) + iall_G(3,ig)*B_matrix_SC(3,:)
        trial_pc_g = SC_G(:) - folding_G(:)
        if(vec_in_latt(vec=trial_pc_g, latt=b_matrix_pc, tolerance=tol_for_vec_equality))then
             !$omp critical
             call append(item=ig,list=selected_coeff_indices) ! pc_g + folding_G has to belong to the 
             !$omp end critical
         endif
    enddo

end subroutine select_coeffs_to_calc_spectral_weights


subroutine calc_spectral_weights(spectral_weight,coeff,selected_coeff_indices,add_elapsed_time_to)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: spectral_weight
complex(kind=sp), dimension(:,:), intent(in) :: coeff !Coefficients, coeff(ikpt, iband)
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: n_bands_SC_calculation, alloc_stat, iband, i, icoeff
real(kind=dp) :: stime, ftime

    stime = time()

    n_bands_SC_calculation = size(coeff,dim=2)
    deallocate(spectral_weight, stat=alloc_stat)
    allocate(spectral_weight(1:n_bands_SC_calculation))
    spectral_weight(:) = 0.0d0
    !$omp parallel do &
    !$    schedule(guided) default(none) &
    !$    private(iband,i,icoeff) &
    !$    shared(n_bands_SC_calculation,selected_coeff_indices,spectral_weight,coeff)
    do iband = 1, n_bands_SC_calculation
        do i=1, size(selected_coeff_indices)
            icoeff = selected_coeff_indices(i)
            spectral_weight(iband) = spectral_weight(iband) + abs(coeff(icoeff,iband))**2.0d0
        enddo
    enddo

    if(present(add_elapsed_time_to))then
        ftime = time()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine calc_spectral_weights

subroutine calc_spectral_function(SF_at_pckpt,energies,SC_calc_ener,spectral_weight,sigma,add_elapsed_time_to)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: SF_at_pckpt
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight, &
                                           SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in) :: sigma
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, stime, ftime

    stime = time()

    ! Calculating the spectral function A(pc_kpt,E) 
    n_SC_bands = size(spectral_weight)
    deallocate(SF_at_pckpt, stat=alloc_stat)
    allocate(SF_at_pckpt(1:size(energies)))
    SF_at_pckpt(:) = 0.0d0
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


subroutine calc_delta_N_pckpt(delta_N_pckpt,energies,SC_calc_ener,spectral_weight,sigma)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: delta_N_pckpt
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in) :: sigma
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, delta_e


    delta_e = energies(2)  - energies(1) !! Uniform energy grid 
    n_SC_bands = size(spectral_weight)
    deallocate(delta_N_pckpt, stat=alloc_stat)
    allocate(delta_N_pckpt(1:size(energies)))
    delta_N_pckpt(:) = 0.0d0
    !$omp parallel do &
    !$    schedule(guided) default(none) &
    !$    private(iener,E,i_SC_band,E_SC_band)  &
    !$    shared(energies,n_SC_bands,SC_calc_ener,delta_N_pckpt,spectral_weight,delta_e,sigma)
    do iener=1,size(energies)  ! For E in the list of energy values for the E(k) plot
        E = energies(iener)
        do i_SC_band=1,n_SC_bands
            E_SC_band = SC_calc_ener(i_SC_band)
            delta_N_pckpt(iener) = delta_N_pckpt(iener) + &
                                   spectral_weight(i_SC_band)*integral_delta_x_minus_x0(x0=E_SC_band, &
                                                                                        lower_lim=E-0.5d0*delta_e, &
                                                                                        upper_lim=E+0.5d0*delta_e, &
                                                                                        std_dev=sigma)
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
    sigma = epsilon(1D0)
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
type(delta_Ns_for_output), intent(out) :: delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
type(delta_Ns), intent(in) :: delta_N
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
type(irr_bz_directions), dimension(:), intent(in) :: dirs_req_for_symmavgd_EBS_along_pcbz_dir
integer :: nener,i_selec_pcbz_dir,ipc_kpt,i_needed_dirs
real(kind=dp), dimension(:), allocatable :: avrgd_dNs
real(kind=dp) :: weight

    call allocate_delta_Ns_for_output_type(delta_N_only_selected_dirs,pckpts_to_be_checked)
    call allocate_delta_Ns_for_output_type(delta_N_symm_avrgd_for_EBS,pckpts_to_be_checked)
    nener = size(delta_N%selec_pcbz_dir(1)%needed_dir(1)%pckpt(1)%dN(:))
    allocate(avrgd_dNs(1:nener))
    do i_selec_pcbz_dir=1,size(delta_N%selec_pcbz_dir(:))
        delta_N_only_selected_dirs%pcbz_dir(i_selec_pcbz_dir) = delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)
        do ipc_kpt=1, size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(1)%pckpt(:))
            avrgd_dNs(:) = 0.0d0
            allocate(delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%dN(1:nener))
            do i_needed_dirs=1,size(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(:))
                weight = dirs_req_for_symmavgd_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(i_needed_dirs)%weight
                avrgd_dNs = avrgd_dNs + weight*delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN(:)
            enddo
            delta_N_symm_avrgd_for_EBS%pcbz_dir(i_selec_pcbz_dir)%pckpt(ipc_kpt)%dN(:) = avrgd_dNs(:)
       enddo     
    enddo

end subroutine get_delta_Ns_for_output


subroutine append(item, list)
implicit none
integer, intent(in) :: item
integer, dimension(:), allocatable, intent(inout) :: list
integer, dimension(:), allocatable :: new_list

if(.not.allocated(list))then
    allocate(list(1:1))
    list(1) = item
else
    allocate(new_list(1:size(list)+1))
    new_list(1:size(list)) = list
    new_list(size(list)+1) = item
    deallocate(list)
    allocate(list(1:size(new_list)))
    list = new_list
    deallocate(new_list)
endif

end subroutine append


end module band_unfolding


