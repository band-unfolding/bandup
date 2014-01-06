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

program BandUP_main
use math
use io_routines
use read_wavecar
use band_unfolding
!$ use omp_lib
implicit none
character(len=127), parameter :: input_file_prim_cell = 'prim_cell_lattice.in', &
                                 input_file_pc_kpts = 'KPOINTS_prim_cell.in', &
                                 input_file_energies = 'energy_info.in', &
                                 output_file = 'unfolded_band_structure.dat'
type(vec3d), dimension(:), allocatable :: pckpts_to_be_checked, list_SC_kpts_in_wavecar
real(kind=dp), dimension(:), allocatable :: ener_SC_bands 
complex(kind=sp), dimension(:,:), allocatable :: coefficients ! Coefficients, coefficients(ikpt, iband)
integer(8) :: file_size
integer :: nplane, nkpts, n_input_pc_kpts, ipc_kpt, i_SCKPT,     & ! nkpts is the total number of k-points in the sampling of the SCBZ
           n_bands_SC_calculation, alloc_stat, n_folding_pckpts, n_folding_pckpts_parsed
integer, dimension(:), allocatable :: nkpts_dirs, selected_coeff_indices, iselec_searched_energies
integer, dimension(:,:), allocatable :: iall_G ! Coordinates of the G points in the [b1, b2, b3] basis.
real(kind=dp) :: e_fermi,ecut, vbz, dist_btween_neigh_searched_pckpts, tol_for_vec_equality, &
                 E_start, E_end, delta_e, stime, ftime, nearest_neigh_dist_SC_kpts,                &
                 total_time_reading_wavecar, elapsed, zero_of_kpts_scale,          &
                 time_calc_spectral_weights,time_calc_spectral_function,time_spent_calculating_delta_Ns
real(kind=dp), dimension(1:3) :: pc_kpt, current_SCKPT_coords, current_SCKPT, folding_G
real(kind=dp), dimension(:), allocatable :: spectral_weight, list_of_energy_values, spectral_function_pckpt
real(kind=dp), dimension(1:3,1:3) :: dir_latt_pc, b_matrix_pc, B_matrix_SC
real(kind=dp), dimension(:,:), allocatable :: pc_E_vs_k, EBS_delta_Ns
logical :: calc_spec_func_explicitly, coeffs_read_once_for_current_SCKPT, wavecar_contains_needed_coeffs
logical, dimension(:,:), allocatable :: fold_into_each_other
!!*****************************************************************************************************************************************
calc_spec_func_explicitly = .FALSE.

!$ call omp_set_dynamic(.FALSE.)
stime = time()
total_time_reading_wavecar = 0.0d0
time_calc_spectral_weights = 0.0d0
time_calc_spectral_function = 0.0d0
time_spent_calculating_delta_Ns = 0.0d0

call print_welcome_messages()
! Getting info about the SC rec. latt. from WAVECAR
call read_from_wavefunc_file(rec_latt_B_matrix=B_matrix_SC, total_nkpts=nkpts, ENCUT=ecut,file_size_in_bytes=file_size)

call read_unit_cell(input_file=input_file_prim_cell, latt=dir_latt_pc(:,:))
call get_rec_latt(latt=dir_latt_pc,rec_latt=b_matrix_pc,rec_latt_vol=vbz)

call get_SCKPTS_contained_in_wavecar(nkpts,B_matrix_SC,list_SC_kpts_in_wavecar,nearest_neigh_dist_SC_kpts)
dist_btween_neigh_searched_pckpts = max(min(0.5d0*nearest_neigh_dist_SC_kpts,1E-4),1E-5)
tol_for_vec_equality = 0.5d0*dist_btween_neigh_searched_pckpts

call define_pckpts_to_be_checked(input_file_pc_kpts, b_matrix_pc(:,:), &
                                 list=pckpts_to_be_checked, n_kpts_dirs=nkpts_dirs, &
                                 min_dk=dist_btween_neigh_searched_pckpts,zero_of_kpts_scale=zero_of_kpts_scale)
n_input_pc_kpts = size(pckpts_to_be_checked(:))

call generate_folding_table(fold_into_each_other,list_SC_kpts_in_wavecar,pckpts_to_be_checked,B_matrix_SC,tol_for_vec_equality)
n_folding_pckpts = count(fold_into_each_other==.TRUE.)
call print_folding_table_results(fold_into_each_other(:,:),list_SC_kpts_in_wavecar(:),pckpts_to_be_checked(:),b_matrix_pc,B_matrix_SC)

call read_energy_info_for_band_search(input_file_energies,e_fermi,E_start,E_end,delta_e)
call real_seq(first_term=E_start,last_term=E_end,increment=delta_e, return_list=list_of_energy_values)

call print_last_messages_before_unfolding(file_size,B_matrix_SC,vbz,E_start,E_end,delta_e,e_fermi)

n_folding_pckpts_parsed = 0
! By looping first over the SCKPTS I avoid reading the wavecar more than once for
!     each SCKPT, saving thus a lot of time when parsing big files 
do i_SCKPT=1,nkpts
    coeffs_read_once_for_current_SCKPT = .FALSE.
    current_SCKPT = list_SC_kpts_in_wavecar(i_SCKPT)%coord(:)
    do ipc_kpt=1, n_input_pc_kpts
        if(fold_into_each_other(i_SCKPT,ipc_kpt))then
            pc_kpt(:) = pckpts_to_be_checked(ipc_kpt)%coord(:)
            folding_G(:) = pc_kpt(:) - current_SCKPT(:)
            ! Reading WAVECAR file
            if(.not.coeffs_read_once_for_current_SCKPT)then
                deallocate(ener_SC_bands,coefficients,iall_G,stat=alloc_stat)
                call read_from_wavefunc_file(spin_channel=1,i_selected_kpt=i_SCKPT, kpt_frac_coords=current_SCKPT_coords, &
                                             energies_bands=ener_SC_bands, n_plane_waves=nplane, &
                                             coeff=coefficients,i_allg_in_file=iall_G, &
                                             elapsed_time=elapsed,add_elapsed_time_to=total_time_reading_wavecar)
                write(*,'(A,I0,A,f0.1,A)')'Info for SC-Kpoint #',i_SCKPT,' read from WAVECAR in ',elapsed,'s.'
                n_bands_SC_calculation = size(coefficients,dim=2)
                coeffs_read_once_for_current_SCKPT = .TRUE.
            endif
            call select_coeffs_to_calc_spectral_weights(selected_coeff_indices,iall_G,b_matrix_pc,B_matrix_SC,nplane,folding_G,tol_for_vec_equality)

            wavecar_contains_needed_coeffs = allocated(selected_coeff_indices)
            if(wavecar_contains_needed_coeffs)then
                ! Calculating spectral_weights
                n_folding_pckpts_parsed = n_folding_pckpts_parsed + 1
                call calc_spectral_weights(coefficients,selected_coeff_indices,spectral_weight,add_elapsed_time_to=time_calc_spectral_weights)
                deallocate(selected_coeff_indices,stat=alloc_stat)

                if(calc_spec_func_explicitly)then
                    ! Calculating the spectral_function (optional)
                    call calc_spectral_function(list_of_energy_values,0.025d0,spectral_weight,ener_SC_bands,spectral_function_pckpt, &
                                                add_elapsed_time_to=time_calc_spectral_function)
                endif

                ! Calculating the delta_Ns
                call get_delta_Ns_for_EBS(pc_E_vs_k, EBS_delta_Ns, &
                                          ipc_kpt,n_input_pc_kpts,list_of_energy_values,spectral_weight,ener_SC_bands, &
                                          add_elapsed_time_to=time_spent_calculating_delta_Ns)
                ! Done.

                if(coeffs_read_once_for_current_SCKPT.and.(ipc_kpt==n_input_pc_kpts))then
                    deallocate(iall_G,stat=alloc_stat)
                    deallocate(coefficients,stat=alloc_stat)
                    deallocate(ener_SC_bands,stat=alloc_stat)
                endif
                deallocate(spectral_weight,stat=alloc_stat)
                deallocate(spectral_function_pckpt, iselec_searched_energies, stat=alloc_stat)
            else
                write(*,'(A)')"WARNING: Not enough coefficients in the WAVECAR file to parse this wave-vector."
            endif ! If wavecar_contains_needed_coeffs
        endif ! If pc k folds into SC K
    enddo ! Loop over pckpts
enddo ! Loop over SC KPTS

call say_goodbye_and_save_results(n_input_pc_kpts,n_folding_pckpts,n_folding_pckpts_parsed, &
                                  pc_E_vs_k, EBS_delta_Ns, pckpts_to_be_checked, nkpts_dirs, &
                                  zero_of_kpts_scale, e_fermi, output_file)

ftime = time()
call print_final_times(stime,ftime, &
                 total_time_reading_wavecar,time_calc_spectral_weights,&
                 calc_spec_func_explicitly,time_calc_spectral_function,&
                 time_spent_calculating_delta_Ns)

!!********************************************************************************************************************************
end program BandUP_main
