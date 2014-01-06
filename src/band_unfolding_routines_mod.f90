!! Copyright (C) 2013 Paulo V. C. Medeiros
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
PUBLIC :: generate_folding_table, select_coeffs_to_calc_spectral_weights,calc_spectral_weights, &
          calc_spectral_function, get_delta_Ns_for_EBS

CONTAINS


subroutine generate_folding_table(fold_into_each_other,list_SC_kpts_in_wavecar,pckpts_to_be_checked,B_matrix_SC,tol_for_vec_equality)
implicit none
logical, dimension(:,:), allocatable, intent(out) :: fold_into_each_other
type(vec3d), dimension(:), allocatable, intent(in) :: list_SC_kpts_in_wavecar, pckpts_to_be_checked
real(kind=dp), dimension(1:3,1:3), intent(in) :: B_matrix_SC
real(kind=dp), intent(in) :: tol_for_vec_equality
integer :: nkpts, n_input_pc_kpts, i_SCKPT,ipc_kpt
logical, dimension(:), allocatable :: pc_kpt_already_folded
real(kind=dp), dimension(1:3) :: pc_kpt, current_SCKPT,trial_folding_G

    nkpts = size(list_SC_kpts_in_wavecar)
    n_input_pc_kpts = size(pckpts_to_be_checked)
    allocate(pc_kpt_already_folded(1:n_input_pc_kpts))
    pc_kpt_already_folded(:) = .FALSE.
    allocate(fold_into_each_other(1:nkpts,1:n_input_pc_kpts))
    fold_into_each_other(:,:) = .FALSE.
    
    do i_SCKPT=1,nkpts
        do ipc_kpt=1, n_input_pc_kpts
            if(pc_kpt_already_folded(ipc_kpt))then
                cycle ! It has already folded before. Since a pc-kpt can fold into only one SC-KPT, let's move on to the next pc-kpt.
            endif
            pc_kpt(:) = pckpts_to_be_checked(ipc_kpt)%coord(:)
            current_SCKPT = list_SC_kpts_in_wavecar(i_SCKPT)%coord(:)
            trial_folding_G(:) = pc_kpt(:) - current_SCKPT(:)
            if(vec_in_latt(vec=trial_folding_G, latt=B_matrix_SC,tolerance=tol_for_vec_equality))then
                fold_into_each_other(i_SCKPT,ipc_kpt) = .TRUE.
                pc_kpt_already_folded(ipc_kpt) = .TRUE. 
            endif
        enddo
    enddo

end subroutine generate_folding_table


subroutine select_coeffs_to_calc_spectral_weights(selected_coeff_indices, &
                                                  iall_G,b_matrix_pc,B_matrix_SC,nplane,folding_G,&
                                                  tol_for_vec_equality)
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
             call append(item=ig,list=selected_coeff_indices)
             !$omp end critical
         endif
    enddo

end subroutine select_coeffs_to_calc_spectral_weights


subroutine calc_spectral_weights(coeff,selected_coeff_indices,spectral_weight,add_elapsed_time_to)
implicit none
complex(kind=sp), dimension(:,:), intent(in) :: coeff !Coefficients, coeff(ikpt, iband)
integer, dimension(:), intent(in) :: selected_coeff_indices
real(kind=dp), dimension(:), allocatable, intent(out) :: spectral_weight
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

subroutine calc_spectral_function(energies,sigma,spectral_weight,SC_calc_ener,spectral_function_pckpt,add_elapsed_time_to)
implicit none
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight
real(kind=dp), intent(in) :: sigma
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), dimension(:), allocatable, intent(out) :: spectral_function_pckpt
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, stime, ftime

    stime = time()

    ! Calculating the spectral function A(pc_kpt,E) 
    n_SC_bands = size(spectral_weight)
    deallocate(spectral_function_pckpt, stat=alloc_stat)
    allocate(spectral_function_pckpt(1:size(energies)))
    spectral_function_pckpt(:) = 0.0d0
    !$omp parallel do &
    !$    schedule(guided) default(none) &
    !$    private(iener,E,i_SC_band,E_SC_band)  &
    !$    shared(energies,n_SC_bands,SC_calc_ener,spectral_function_pckpt,spectral_weight,sigma)
    do iener=1,size(energies)  ! For E in the list of energy values for the E(k) plot
        E = energies(iener)
        do i_SC_band=1,n_SC_bands
            E_SC_band = SC_calc_ener(i_SC_band)
            spectral_function_pckpt(iener) = spectral_function_pckpt(iener) + &
                                             spectral_weight(i_SC_band)*delta(E_SC_band - E, std_dev=sigma)
        enddo
    enddo

    if(present(add_elapsed_time_to))then
        ftime = time()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine calc_spectral_function


function integral_delta_x_minus_x0(x0,lower_lim,upper_lim,std_dev) result(rtn)
implicit none
real(kind=dp) :: rtn
real(kind=dp), intent(in) :: x0,lower_lim,upper_lim,std_dev
real(kind=dp) :: x1, x2


    x1 = (lower_lim - x0)/(sqrt(2.0d0)*std_dev)
    x2 = (upper_lim - x0)/(sqrt(2.0d0)*std_dev)
    rtn = 0.5d0*(erf(x2) - erf(x1))

end function integral_delta_x_minus_x0


subroutine calc_delta_N_pckpt(delta_N_pckpt,energies,spectral_weight,SC_calc_ener,sigma)
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: delta_N_pckpt
real(kind=dp), dimension(:), intent(in) :: energies, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in) :: sigma
integer :: alloc_stat,n_SC_bands,iener,i_SC_band
real(kind=dp) :: E, E_SC_band, delta_e


    delta_e = energies(2)  - energies(1)   
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


subroutine get_delta_Ns_for_EBS(pc_E_vs_k, pc_E_vs_k_intensity, &
                                                ipc_kpt, n_input_pc_kpts, energy_grid,spectral_weight,SC_calc_ener, &
                                                smearing_for_delta_function, add_elapsed_time_to)
real(kind=dp), dimension(:,:), allocatable, intent(inout) :: pc_E_vs_k,pc_E_vs_k_intensity
integer, intent(in) :: ipc_kpt, n_input_pc_kpts
real(kind=dp), dimension(:), intent(in) :: energy_grid, spectral_weight
real(kind=dp), dimension(:), intent(in) :: SC_calc_ener ! E(K;i_SC_band)
real(kind=dp), intent(in), optional :: smearing_for_delta_function
integer :: n_ener, iener, ipc_kpt_aux
real(kind=dp) :: sigma, stime, ftime
real(kind=dp), dimension(:), allocatable :: delta_N_pckpt
real(kind=dp), intent(inout), optional :: add_elapsed_time_to

    stime = time()

    sigma = epsilon(1D0)
    if(present(smearing_for_delta_function))then
        sigma = max(sigma,abs(smearing_for_delta_function))
    endif
    n_ener = size(energy_grid)
    if(.not.allocated(pc_E_vs_k))then
        allocate(pc_E_vs_k(1:n_ener,1:n_input_pc_kpts))
        do ipc_kpt_aux = 1, n_input_pc_kpts
            pc_E_vs_k(:,ipc_kpt_aux) = energy_grid(:)
        enddo
    endif
    if(.not.allocated(pc_E_vs_k_intensity))then
        allocate(pc_E_vs_k_intensity(1:n_ener,1:n_input_pc_kpts))
        pc_E_vs_k_intensity(:,:) = 0.0d0
    endif
    call calc_delta_N_pckpt(delta_N_pckpt,energy_grid,spectral_weight,SC_calc_ener,sigma)
    do iener=1,n_ener
        pc_E_vs_k_intensity(iener,ipc_kpt) = delta_N_pckpt(iener)
    enddo
    
    if(present(add_elapsed_time_to))then
        ftime = time()
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine get_delta_Ns_for_EBS


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


