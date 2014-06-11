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
use general_io
use io_routines
use read_wavecar
use read_vasp_files
use write_vasp_files
use band_unfolding
!$ use omp_lib
implicit none
character(len=127), parameter :: input_file_prim_cell = 'prim_cell_lattice.in', &
                                 input_file_supercell = 'supercell_lattice.in', &
                                 input_file_pc_kpts = 'KPOINTS_prim_cell.in', &
                                 input_file_energies = 'energy_info.in', &
                                 output_file_only_user_selec_direcs = 'unfolded_EBS_not-symmetry_averaged.dat', &
                                 output_file_symm_averaged_EBS = 'unfolded_EBS_symmetry-averaged.dat'
type(vec3d), dimension(:), allocatable :: list_of_SCKPTS
type(crystal_3D), allocatable :: crystal_pc, crystal_SC
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
real(kind=dp), dimension(1:3) :: symmetrized_unf_pc_kpt, current_SCKPT, folding_G
real(kind=dp), dimension(:), allocatable :: spectral_weight, energy_grid
real(kind=dp), dimension(1:3,1:3) :: b_matrix_pc,B_matrix_SC,SC_latt_vecs,matrix_M
real(kind=dp), dimension(:,:), allocatable :: k_starts,k_ends
logical :: calc_spec_func_explicitly,stop_when_a_pckpt_cannot_be_parsed,stop_if_not_commensurate,pckpt_folds, &
           coeffs_read_once_for_current_SCKPT, wavecar_contains_needed_coeffs, are_commens, crystal_SC_read_from_file, &
           pc_is_prim_cell, write_attempted_pc_corresp_to_input_pc, SC_is_prim_cell, write_attempted_pc_corresp_to_SC, &
           stop_if_GUR_fails
type(irr_bz_directions), dimension(:), allocatable :: dirs_req_for_symmavgd_EBS_along_pcbz_dir
type(geom_unfolding_relations_for_each_SCKPT) :: GUR !! Geometric Unfolding Relations
type(selected_pcbz_directions) :: pckpts_to_be_checked
type(delta_Ns) :: delta_N
type(delta_Ns_for_output) :: delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
type(crystal_3D), allocatable :: crystal_SC_reduced_to_prim_cell, crystal_pc_reduced_to_prim_cell
!!*****************************************************************************************************************************************
calc_spec_func_explicitly = .FALSE.
stop_when_a_pckpt_cannot_be_parsed = .TRUE.
stop_if_not_commensurate = .FALSE.
stop_if_GUR_fails = .TRUE.
write_attempted_pc_corresp_to_input_pc = .TRUE.
write_attempted_pc_corresp_to_SC = .TRUE.
!$ call omp_set_dynamic(.FALSE.)
stime = time()
total_time_reading_wavecar = 0.0_dp
time_calc_spectral_weights = 0.0_dp
time_calc_spectral_function = 0.0_dp
time_spent_calculating_delta_Ns = 0.0_dp

call print_welcome_messages(package_version)
call read_from_wavefunc_file(rec_latt_B_matrix=B_matrix_SC, latt_vecs=SC_latt_vecs, &
                             total_nkpts=nkpts, ENCUT=ecut,file_size_in_bytes=file_size)
call get_crystal_from_file(crystal_pc,input_file=input_file_prim_cell, &
                           stop_if_file_not_found=.TRUE.)
call get_crystal_from_file(crystal_SC,input_file=input_file_supercell, &
                           stop_if_file_not_found=.FALSE.,success=crystal_SC_read_from_file)
if(.not. crystal_SC_read_from_file)then
    call create_crystal(crystal_SC, latt_vecs=SC_latt_vecs)
endif
call get_rec_latt(latt=crystal_pc%latt_vecs(:,:),rec_latt=b_matrix_pc,rec_latt_vol=vbz)
call check_if_pc_and_SC_are_commensurate(commensurate=are_commens, M=matrix_M, &
                                         b_matrix_pc=b_matrix_pc, B_matrix_SC=B_matrix_SC, &
                                         tol=default_tol_for_int_commens_test)
call print_message_commens_test(commensurate=are_commens,M=matrix_M,stop_if_not_commens=stop_if_not_commensurate) 
if(stop_if_not_commensurate .and. .not. are_commens) stop
call get_SCKPTS_contained_in_wavecar(nkpts,B_matrix_SC,list_of_SCKPTS)
call read_pckpts_selected_by_user(k_starts=k_starts, k_ends=k_ends, ndirs=n_selec_pcbz_dirs, n_kpts_dirs=n_pckpts_dirs, &
                                  input_file=input_file_pc_kpts,b_matrix_pc=b_matrix_pc(:,:), &
                                  zero_of_kpts_scale=zero_of_kpts_scale)
call get_all_irr_dirs_req_for_symmavgd_EBS(dirs_req_for_symmavgd_EBS_along_pcbz_dir, &
                                           n_dirs_for_EBS_along_pcbz_dir, &
                                           neqv_dirs_pcbz, neqv_dirs_SCBZ, &
                                           ncompl_dirs, n_irr_compl_dirs,&
                                           crystal_pc_reduced_to_prim_cell, &
                                           pc_is_prim_cell, &
                                           crystal_SC_reduced_to_prim_cell, &
                                           SC_is_prim_cell, &
                                           crystal_pc=crystal_pc, crystal_SC=crystal_SC, &
                                           k_starts=k_starts(:,:),k_ends=k_ends(:,:))
call write_attempted_pc_assoc_with_input_unit_cell_and_SC(crystal_pc_reduced_to_prim_cell, &
                                                                crystal_SC_reduced_to_prim_cell, &
                                                                pc_is_prim_cell,SC_is_prim_cell, &
                                                                write_attempted_pc_corresp_to_input_pc, &
                                                                write_attempted_pc_corresp_to_SC)
call print_symm_analysis_for_selected_pcbz_dirs(dirs_req_for_symmavgd_EBS_along_pcbz_dir, &
                                                neqv_dirs_pcbz, neqv_dirs_SCBZ, ncompl_dirs, n_irr_compl_dirs)
call define_pckpts_to_be_checked(pckpts_to_be_checked,dirs_req_for_symmavgd_EBS_along_pcbz_dir,n_pckpts_dirs(:))
call get_geom_unfolding_relations(GUR,list_of_SCKPTS,pckpts_to_be_checked,crystal_SC)
call print_message_success_determining_GUR(GUR, stop_if_GUR_fails, is_main_code=.TRUE.) 
if((GUR%n_pckpts /= GUR%n_folding_pckpts) .and. stop_if_GUR_fails) stop
call print_geom_unfolding_relations(GUR,list_of_SCKPTS,b_matrix_pc,B_matrix_SC)
call read_energy_info_for_band_search(input_file_energies,e_fermi,E_start,E_end,delta_e)
call real_seq(first_term=E_start,last_term=E_end,increment=delta_e, return_list=energy_grid)
call print_last_messages_before_unfolding(file_size,nkpts,B_matrix_SC,vbz,E_start,E_end,delta_e,e_fermi)
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
                    ! The "symmetrized" unfolding pc_kpt is a little trick I came up with
                    ! that allows me to use the coefficients of a SC wavefunction psi(K',n) to
                    ! calculate the spectral weights associated with a SC wavefunction psi(K,n),
                    ! where K' = SK and S is a symmetry operation of the crystal's point group.
                    symmetrized_unf_pc_kpt(:) = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%Scoords(:)
                    folding_G(:) = symmetrized_unf_pc_kpt(:) - current_SCKPT(:)
                    ! Reading WAVECAR file
                    if(.not.coeffs_read_once_for_current_SCKPT)then
                        write(*,"(A,I0,A)")'Reading plane-wave coefficients for SC-Kpoint K(',i_SCKPT,') from the wavefunctions file... '
                        deallocate(ener_SC_bands,coefficients,iall_G,stat=alloc_stat)
                        call read_from_wavefunc_file(spin_channel=1,i_selected_kpt=i_SCKPT, &
                                                     energies_bands=ener_SC_bands, n_plane_waves=nplane, &
                                                     coeff=coefficients,i_allg_in_file=iall_G, &
                                                     elapsed_time=elapsed,add_elapsed_time_to=total_time_reading_wavecar)
                        write(*,'(A,f0.1,A)')'    * Done in ',elapsed,'s.'
                        n_bands_SC_calculation = size(coefficients,dim=2)
                        coeffs_read_once_for_current_SCKPT = .TRUE.
                    endif
                    call select_coeffs_to_calc_spectral_weights(selected_coeff_indices,iall_G,b_matrix_pc,B_matrix_SC,nplane, &
                                                                folding_G, default_tol_for_vec_equality)
                    wavecar_contains_needed_coeffs = allocated(selected_coeff_indices)
                    if(wavecar_contains_needed_coeffs)then
                        ! Calculating spectral_weights
                        n_folding_pckpts_parsed = n_folding_pckpts_parsed + 1
                        call calc_spectral_weights(spectral_weight,coefficients,selected_coeff_indices,add_elapsed_time_to=time_calc_spectral_weights)
                        deallocate(selected_coeff_indices,stat=alloc_stat)
                        if(calc_spec_func_explicitly)then
                            ! Calculating the spectral_function (optional)
                            call calc_spectral_function(delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%SF, &
                                                        energy_grid,ener_SC_bands,spectral_weight,sigma=0.025_dp, &
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
                        call print_message_pckpt_cannot_be_parsed(stop_when_a_pckpt_cannot_be_parsed)
                        if(stop_when_a_pckpt_cannot_be_parsed) stop
                    endif ! If wavecar_contains_needed_coeffs
                endif ! If pc k folds into SC K
            enddo ! Loop over pckpts
        enddo ! Loop over all directions actually needed to get a proper EBS along the selected directions
    enddo ! Loop over selected directions
enddo ! Loop over SC KPTS

!! Preparing output
!! >> delta_N_only_selected_dirs:  Unfolded EBS strictly along the directions requested by the user.
!!                                 Can be good, for instance, to investigate anisotropy effects.
!! >> delta_N_symm_avrgd_for_EBS: Symmetry-averaged unfolded EBS. 
!!                                This is the kind of EBS you'll see in my paper [Phys. Rev. B 89, 041407(R) (2014)].
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
