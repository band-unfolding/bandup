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
                                 output_file_only_user_selec_direcs = 'unfolded_EBS_not-symmetry_averaged.dat', &
                                 output_file_symm_averaged_EBS = 'unfolded_EBS_symmetry-averaged.dat'
type(vec3d), dimension(:), allocatable :: list_of_SCKPTS
real(kind=dp), dimension(:), allocatable :: ener_SC_bands 
complex(kind=sp), dimension(:,:), allocatable :: coefficients ! Coefficients, coefficients(ikpt, iband)
integer(8) :: file_size
integer :: nplane, nkpts, ipc_kpt, i_SCKPT, n_selec_pcbz_dirs,   & ! nkpts is the total number of k-points in the sampling of the SCBZ
           n_bands_SC_calculation, alloc_stat, n_folding_pckpts_parsed, &
           i_selec_pcbz_dir, i_needed_dirs 
integer, dimension(:), allocatable :: selected_coeff_indices,n_pckpts_dirs,n_dirs_for_EBS_along_pcbz_dir, &
                                      neqv_dirs_pcbz, neqv_dirs_SCBZ, ncompl_dirs, n_irr_compl_dirs
integer, dimension(:,:), allocatable :: iall_G ! Coordinates of the G points in the [b1, b2, b3] basis.
real(kind=dp) :: e_fermi,ecut, vbz, &
                 E_start, E_end, delta_e, stime, ftime,      &
                 total_time_reading_wavecar, elapsed, zero_of_kpts_scale,          &
                 time_calc_spectral_weights,time_calc_spectral_function,time_spent_calculating_delta_Ns
real(kind=dp), dimension(1:3) :: symmetrized_unf_pc_kpt, current_SCKPT_coords, current_SCKPT, folding_G
real(kind=dp), dimension(:), allocatable :: spectral_weight, energy_grid
real(kind=dp), dimension(1:3,1:3) :: dir_latt_pc, b_matrix_pc,B_matrix_SC
real(kind=dp), dimension(:,:), allocatable :: k_starts,k_ends
logical :: calc_spec_func_explicitly,stop_when_a_pckpt_cannot_be_parsed,pckpt_folds, &
           coeffs_read_once_for_current_SCKPT, wavecar_contains_needed_coeffs
type(irr_bz_directions), dimension(:), allocatable :: dirs_req_for_symmavgd_EBS_along_pcbz_dir
type(geom_unfolding_relations_for_each_SCKPT) :: GUR !! Geometric Unfolding Relations
type(selected_pcbz_directions) :: pckpts_to_be_checked
type(delta_Ns) :: delta_N
type(delta_Ns_for_output) :: delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
!!*****************************************************************************************************************************************
calc_spec_func_explicitly = .FALSE.
stop_when_a_pckpt_cannot_be_parsed = .FALSE.
!$ call omp_set_dynamic(.FALSE.)
stime = time()
total_time_reading_wavecar = 0.0d0
time_calc_spectral_weights = 0.0d0
time_calc_spectral_function = 0.0d0
time_spent_calculating_delta_Ns = 0.0d0

call print_welcome_messages(package_version)
call read_from_wavefunc_file(rec_latt_B_matrix=B_matrix_SC, total_nkpts=nkpts, ENCUT=ecut,file_size_in_bytes=file_size)
call read_unit_cell(input_file=input_file_prim_cell, latt=dir_latt_pc(:,:))
call get_rec_latt(latt=dir_latt_pc,rec_latt=b_matrix_pc,rec_latt_vol=vbz)
call get_SCKPTS_contained_in_wavecar(nkpts,B_matrix_SC,list_of_SCKPTS)
call read_pckpts_selected_by_user(k_starts=k_starts, k_ends=k_ends, ndirs=n_selec_pcbz_dirs, n_kpts_dirs=n_pckpts_dirs, &
                                  input_file=input_file_pc_kpts,b_matrix_pc=b_matrix_pc(:,:), &
                                  zero_of_kpts_scale=zero_of_kpts_scale)
call get_all_irr_dirs_req_for_symmavgd_EBS(dirs_req_for_symmavgd_EBS_along_pcbz_dir, n_dirs_for_EBS_along_pcbz_dir, &
                                           neqv_dirs_pcbz, neqv_dirs_SCBZ, ncompl_dirs, n_irr_compl_dirs,&
                                           b_matrix_pc=b_matrix_pc(:,:),B_matrix_SC=B_matrix_SC(:,:), &
                                           k_starts=k_starts(:,:),k_ends=k_ends(:,:))
call print_symm_analysis_for_selected_pcbz_dirs(dirs_req_for_symmavgd_EBS_along_pcbz_dir, &
                                                neqv_dirs_pcbz, neqv_dirs_SCBZ, ncompl_dirs, n_irr_compl_dirs)
call define_pckpts_to_be_checked(pckpts_to_be_checked,dirs_req_for_symmavgd_EBS_along_pcbz_dir,n_pckpts_dirs(:))
call get_geom_unfolding_relations(GUR,list_of_SCKPTS,pckpts_to_be_checked,B_matrix_SC)
call print_geom_unfolding_relations(GUR,list_of_SCKPTS,b_matrix_pc,B_matrix_SC)
call read_energy_info_for_band_search(input_file_energies,e_fermi,E_start,E_end,delta_e)
call real_seq(first_term=E_start,last_term=E_end,increment=delta_e, return_list=energy_grid)
call print_last_messages_before_unfolding(file_size,B_matrix_SC,vbz,E_start,E_end,delta_e,e_fermi)

call allocate_delta_Ns_type(delta_N,pckpts_to_be_checked)
n_folding_pckpts_parsed = 0
! By looping first over the SCKPTS I avoid reading the wavecar more than once for
!     each SCKPT, saving thus a lot of time when parsing big files 
do i_SCKPT=1,nkpts
    coeffs_read_once_for_current_SCKPT = .FALSE.
    current_SCKPT = list_of_SCKPTS(i_SCKPT)%coord(:)
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        do i_needed_dirs=1,n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)
            do ipc_kpt=1, n_pckpts_dirs(i_selec_pcbz_dir)
                pckpt_folds = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds
                if(pckpt_folds)then
                    symmetrized_unf_pc_kpt(:) = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:)
                    folding_G(:) = symmetrized_unf_pc_kpt(:) - current_SCKPT(:)
                    ! Reading WAVECAR file
                    if(.not.coeffs_read_once_for_current_SCKPT)then
                        deallocate(ener_SC_bands,coefficients,iall_G,stat=alloc_stat)
                        call read_from_wavefunc_file(spin_channel=1,i_selected_kpt=i_SCKPT, kpt_frac_coords=current_SCKPT_coords, &
                                                     energies_bands=ener_SC_bands, n_plane_waves=nplane, &
                                                     coeff=coefficients,i_allg_in_file=iall_G, &
                                                     elapsed_time=elapsed,add_elapsed_time_to=total_time_reading_wavecar)
                        write(*,'(A,I0,A,f0.1,A)')'Info for SC-Kpoint K(',i_SCKPT,') read from the wavefunctions file in ',elapsed,'s.'
                        n_bands_SC_calculation = size(coefficients,dim=2)
                        coeffs_read_once_for_current_SCKPT = .TRUE.
                    endif
                    call select_coeffs_to_calc_spectral_weights(selected_coeff_indices,iall_G,b_matrix_pc,B_matrix_SC,nplane, &
                                                                folding_G, tol_for_vec_equality)
                    wavecar_contains_needed_coeffs = allocated(selected_coeff_indices)
                    if(wavecar_contains_needed_coeffs)then
                        ! Calculating spectral_weights
                        n_folding_pckpts_parsed = n_folding_pckpts_parsed + 1
                        call calc_spectral_weights(spectral_weight,coefficients,selected_coeff_indices,add_elapsed_time_to=time_calc_spectral_weights)
                        deallocate(selected_coeff_indices,stat=alloc_stat)
                        if(calc_spec_func_explicitly)then
                            ! Calculating the spectral_function (optional)
                            call calc_spectral_function(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%SF, &
                                                        energy_grid,ener_SC_bands,spectral_weight,sigma=0.025d0, &
                                                        add_elapsed_time_to=time_calc_spectral_function)
                        endif
                        ! Calculating the delta_Ns
                        call get_delta_Ns_for_EBS(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%dN, &
                                                  energy_grid,ener_SC_bands,spectral_weight, &
                                                  add_elapsed_time_to=time_spent_calculating_delta_Ns)
                        ! Done.
                        if(coeffs_read_once_for_current_SCKPT  .and. &
                           i_selec_pcbz_dir==n_selec_pcbz_dirs .and. &
                           i_needed_dirs==n_dirs_for_EBS_along_pcbz_dir(i_selec_pcbz_dir) .and. &
                           ipc_kpt==n_pckpts_dirs(i_selec_pcbz_dir))then
                            !! We're done with the current SCKPS when we enter here
                            deallocate(iall_G,coefficients,ener_SC_bands,stat=alloc_stat)
                        endif
                        deallocate(spectral_weight,stat=alloc_stat) !! We don't need to keep the spectral weights on the memory
                    else
                        if(stop_when_a_pckpt_cannot_be_parsed)then
                            write(*,'(A)')"    ERROR: Not enough coefficients in the WAVECAR file to unfold one of the selected pc wave-vectors."
                            write(*,'(A)')"           Stopping now."
                            stop
                        else
                            write(*,'(A)')"    WARNING: Not enough coefficients in the WAVECAR file to unfold this pc wave-vector."
                            write(*,'(A)')"             Be careful with your results: they might be incomplete or even wrong."
                        endif
                    endif ! If wavecar_contains_needed_coeffs
                endif ! If pc k folds into SC K
            enddo ! Loop over pckpts
        enddo ! Loop over all directions actually needed to get a proper EBS along the selected directions
    enddo ! Loop over selected directions
enddo ! Loop over SC KPTS

!! Preparing output
!! >> delta_N_only_selected_dirs:  Unfolded EBS strictly along the directions requested by the user. Can be good, for instance, to investigate anisotropy effects
!! >> delta_N_symm_avrgd_for_EBS: Symmetry-averaged unfolded EBS. This is the kind of EBS you'll see in my paper [Phys. Rev. B 89, 041407(R) (2014)].
call get_delta_Ns_for_output(delta_N_only_selected_dirs,delta_N_symm_avrgd_for_EBS,delta_N,dirs_req_for_symmavgd_EBS_along_pcbz_dir,pckpts_to_be_checked)

call say_goodbye_and_save_results(output_file_only_user_selec_direcs,output_file_symm_averaged_EBS, &
                                  delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS, &
                                  pckpts_to_be_checked,energy_grid,e_fermi,zero_of_kpts_scale, &
                                  GUR%n_pckpts,GUR%n_folding_pckpts,n_folding_pckpts_parsed)

ftime = time()
call print_final_times(stime,ftime, &
                 total_time_reading_wavecar,time_calc_spectral_weights,&
                 calc_spec_func_explicitly,time_calc_spectral_function,&
                 time_spent_calculating_delta_Ns)

!!********************************************************************************************************************************
end program BandUP_main
