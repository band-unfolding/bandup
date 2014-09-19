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
!$ use omp_lib
use constants_and_types
use cla_wrappers
use math
use general_io
use io_routines
use read_wavecar
use read_vasp_files
use write_vasp_files
use band_unfolding
implicit none
type(UnfoldedQuantities) :: delta_N
type(UnfoldedQuantitiesForOutput) :: delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
type(geom_unfolding_relations_for_each_SCKPT) :: GUR
type(irr_bz_directions), dimension(:), allocatable :: all_dirs_used_for_EBS_along_pcbz_dir
type(selected_pcbz_directions) :: pckpts_to_be_checked
type(crystal_3D) :: crystal_SC, crystal_pc
type(vec3d), dimension(:), allocatable :: list_of_SCKPTS
type(timekeeping) :: times
complex(kind=kind_cplx_coeffs), dimension(:,:,:), allocatable :: coefficients
real(kind=dp), dimension(:,:), allocatable :: k_starts,k_ends
real(kind=dp), dimension(:), allocatable :: energy_grid, ener_SC_bands
real(kind=dp) :: e_fermi, vbz, E_start, E_end, delta_e, elapsed, zero_of_kpts_scale
integer, dimension(:,:), allocatable :: iall_G ! Coordinates of the G points in the SCRL.
integer, dimension(:), allocatable :: selected_coeff_indices,n_pckpts_dirs
integer(8) :: file_size
integer :: nkpts, ipc_kpt, i_SCKPT, n_selec_pcbz_dirs, & ! nkpts = #SCkpts
           alloc_stat, n_folding_pckpts_parsed, i_selec_pcbz_dir, i_needed_dirs
logical :: pckpt_folds, crystal_SC_read_from_file, spinor_wf
!!*************************************************************************************************
call initialize(times)
!$ call omp_set_dynamic(.FALSE.)
call print_welcome_messages(package_version)
call get_commline_args(args)

call get_crystal_from_file(crystal_SC,input_file=args%input_file_supercell, &
                           stop_if_file_not_found=.FALSE., success=crystal_SC_read_from_file)
call read_from_wavefunc_file(rec_latt_B_matrix=crystal_SC%rec_latt_vecs, &
                             latt_vecs=crystal_SC%latt_vecs, total_nkpts=nkpts, &
                             file_size_in_bytes=file_size, spinor_wf=spinor_wf)
if(.not. crystal_SC_read_from_file) call create_crystal(crystal_SC, latt_vecs=crystal_SC%latt_vecs)
call get_prim_cell(crystal_SC, symprec=default_symprec)
call get_crystal_from_file(crystal_pc,input_file=args%input_file_prim_cell, stop_if_file_not_found=.TRUE.)
call get_prim_cell(crystal_pc, symprec=default_symprec)
call get_rec_latt(latt=crystal_pc%latt_vecs(:,:),rec_latt=crystal_pc%rec_latt_vecs,rec_latt_vol=vbz)
call write_attempted_pc_assoc_with_input_unit_cell_and_SC(crystal_pc, crystal_SC)
call verify_commens(crystal_pc, crystal_SC, args)
call analise_symm_pc_SC(crystal_pc, crystal_SC)

call get_SCKPTS_contained_in_wavecar(list_of_SCKPTS, crystal_SC%rec_latt_vecs)
call read_pckpts_selected_by_user(k_starts=k_starts, k_ends=k_ends, ndirs=n_selec_pcbz_dirs, n_kpts_dirs=n_pckpts_dirs, &
                                  input_file=args%input_file_pc_kpts,b_matrix_pc=crystal_pc%rec_latt_vecs(:,:), &
                                  zero_of_kpts_scale=zero_of_kpts_scale)

call get_pcbz_dirs_2b_used_for_EBS(all_dirs_used_for_EBS_along_pcbz_dir, crystal_pc, crystal_SC, &
                                   k_starts, k_ends, args)
call print_symm_analysis_for_selected_pcbz_dirs(all_dirs_used_for_EBS_along_pcbz_dir)
call define_pckpts_to_be_checked(pckpts_to_be_checked, all_dirs_used_for_EBS_along_pcbz_dir, n_pckpts_dirs(:))

call get_geom_unfolding_relations(GUR,list_of_SCKPTS,pckpts_to_be_checked,crystal_SC)
call print_message_success_determining_GUR(GUR, stop_if_GUR_fails, is_main_code=.TRUE.) 
if((GUR%n_pckpts /= GUR%n_folding_pckpts) .and. stop_if_GUR_fails) stop
call print_geom_unfolding_relations(GUR, list_of_SCKPTS, crystal_pc, crystal_SC)

call read_energy_info_for_band_search(args%input_file_energies,e_fermi,E_start,E_end,delta_e)
call real_seq(first_term=E_start,last_term=E_end,increment=delta_e, return_list=energy_grid)
call print_last_messages_before_unfolding(file_size,nkpts,crystal_SC%rec_latt_vecs,vbz,E_start,E_end,delta_e,e_fermi, spinor_wf)
call allocate_UnfoldedQuantities(delta_N, pckpts_to_be_checked)

! Main loop
n_folding_pckpts_parsed = 0
do i_SCKPT=1,nkpts
    deallocate(ener_SC_bands, coefficients, iall_G, stat=alloc_stat)
    do i_selec_pcbz_dir=1,n_selec_pcbz_dirs
        do i_needed_dirs=1, size(all_dirs_used_for_EBS_along_pcbz_dir(i_selec_pcbz_dir)%irr_dir(:))
            do ipc_kpt=1, n_pckpts_dirs(i_selec_pcbz_dir)
                pckpt_folds = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds
                if(pckpt_folds)then
                    call update_GUR_indices(GUR, i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt)
                    ! Reading the wavefunction file
                    if(.not. allocated(coefficients))then
                        write(*,"(A,I0,3A)")'Reading plane-wave coefficients for SC-Kpoint K(',i_SCKPT, &
                                            ') from file "', trim(adjustl(args%WF_file)), '"...'
                        call read_from_wavefunc_file(spin_channel=args%spin_channel,i_selected_kpt=i_SCKPT, &
                                                     energies_bands=ener_SC_bands, coeff=coefficients, i_allg_in_file=iall_G, &
                                                     elapsed_time=elapsed,add_elapsed_time_to=times%read_wf)
                        write(*,'(A,f0.1,A)')'    * Done in ',elapsed,'s.'
                    endif
                    call select_coeffs_to_calc_spectral_weights(selected_coeff_indices, iall_G, &
                                                                crystal_pc, crystal_SC, GUR)
                    if(allocated(selected_coeff_indices))then
                        call perform_unfolding(delta_N, times, GUR, coefficients, &
                                               selected_coeff_indices, energy_grid, ener_SC_bands)
                        n_folding_pckpts_parsed = n_folding_pckpts_parsed + 1
                    else
                        call print_message_pckpt_cannot_be_parsed(stop_when_a_pckpt_cannot_be_parsed)
                        if(stop_when_a_pckpt_cannot_be_parsed) stop
                    endif 
                endif 
            enddo
        enddo
    enddo
enddo

!! Preparing output
!! >> delta_N_only_selected_dirs:  Unfolded EBS strictly along the directions requested by the user.
!!                                 Can be good, for instance, to investigate anisotropy effects.
!! >> delta_N_symm_avrgd_for_EBS: Symmetry-averaged unfolded EBS. 
!!                                This is the type of EBS you'll see in my paper [Phys. Rev. B 89, 041407(R) (2014)].
call get_delta_Ns_for_output(delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS, delta_N, &
                             all_dirs_used_for_EBS_along_pcbz_dir, pckpts_to_be_checked)
call say_goodbye_and_save_results(delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS, &
                                  pckpts_to_be_checked,energy_grid, e_fermi, zero_of_kpts_scale, &
                                  GUR%n_pckpts, GUR%n_folding_pckpts, n_folding_pckpts_parsed)
call print_final_times(times)
!!********************************************************************************************************************************
end program BandUP_main
