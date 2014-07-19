!! Copyright (C) 2014 Paulo V. C. Medeiros
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

module options_setup
use kinds
use cla
implicit none
PRIVATE
PUBLIC :: get_commline_args, calc_spec_func_explicitly, &
          stop_when_a_pckpt_cannot_be_parsed, stop_if_GUR_fails, &
          get_all_kpts_needed_for_EBS_averaging, print_GUR_pre_unfolding_utility

!! Hard-coded options I only change for debugging/testing. You probably shouldn't modify this.
logical, parameter :: calc_spec_func_explicitly = .FALSE., &
                      stop_when_a_pckpt_cannot_be_parsed = .TRUE., &
                      stop_if_GUR_fails = .TRUE., &
                      get_all_kpts_needed_for_EBS_averaging = .TRUE., &
                      print_GUR_pre_unfolding_utility = .FALSE.

!! Functions and subroutines
CONTAINS

subroutine get_commline_args(WF_file, input_file_prim_cell, input_file_supercell, &
                             input_file_pc_kpts, input_file_energies, &
                             output_file_symm_averaged_EBS, output_file_only_user_selec_direcs, &
                             spin_channel, stop_if_not_commensurate, &
                             write_attempted_pc_corresp_to_input_pc, write_attempted_pc_corresp_to_SC, &
                             out_file_SC_kpts)
implicit none
character(len=*), intent(out) :: WF_file, input_file_prim_cell, input_file_supercell, &
                                 input_file_pc_kpts, input_file_energies, &
                                 output_file_symm_averaged_EBS, output_file_only_user_selec_direcs
character(len=*), intent(out), optional :: out_file_SC_kpts
integer, intent(out) :: spin_channel
logical, intent(out) :: stop_if_not_commensurate, write_attempted_pc_corresp_to_input_pc, &
                        write_attempted_pc_corresp_to_SC

    WF_file = ''; input_file_prim_cell = ''; input_file_supercell = ''
    input_file_pc_kpts = ''; input_file_energies = ''
    output_file_symm_averaged_EBS = ''; output_file_only_user_selec_direcs = ''
    if(present(out_file_SC_kpts)) out_file_SC_kpts = ''
    ! Optional command line arguments. 
    call cla_init
    call cla_register('-wf_file', 'Path to the wavefunction file.', cla_char, 'WAVECAR')
    call cla_register('-pc_file', 'Path to the primitive cell file.', cla_char, 'prim_cell_lattice.in')
    call cla_register('-sc_file', 'Path to the supercell file.', cla_char, 'supercell_lattice.in')
    call cla_register('-pckpts_file', 'Path to the pc-kpts file.', cla_char, 'KPOINTS_prim_cell.in')
    call cla_register('-out_sckpts_file', 'Path to the SC-kpts file. Used only by the pre-unfolding utility', &
                      cla_char, 'KPOINTS_supercell.out')
    call cla_register('-energy_file', 'Path to the energy configuration file.', cla_char, 'energy_info.in')
    call cla_register('-out_file_nosymm', 'Name of the output file (not symmetry-averaged).', & 
                      cla_char, 'unfolded_EBS_not-symmetry_averaged.dat')
    call cla_register('-out_file_symm', 'Name of the output file (symmetry-averaged).', &
                      cla_char, 'unfolded_EBS_symmetry-averaged.dat')
    call cla_register('-spin_channel', 'Either 1 or 2 (if the wavefunction contains a 2nd spin channel).', &
                      cla_int, '1')
    call cla_register('-stop_if_not_commensurate', &
                      'Stop if the SC and pc are not commensurate.', cla_flag, 'F')
    call cla_register('-skip_propose_pc_for_given_pc', 'Do not attempt to propose a pc &
                                                        for the input reference unit cell.', &
                      cla_flag, 'F')
    call cla_register('-skip_propose_pc_for_given_sc', 'Do not attempt to propose a pc &
                                                        for the input SC', &
                      cla_flag, 'F')
    call cla_validate


   call cla_get('-wf_file', WF_file)
   call cla_get('-pc_file', input_file_prim_cell)
   call cla_get('-sc_file', input_file_supercell)
   call cla_get('-pckpts_file', input_file_pc_kpts)
   call cla_get('-energy_file', input_file_energies)
   call cla_get('-out_file_symm', output_file_symm_averaged_EBS)
   call cla_get('-out_file_nosymm', output_file_only_user_selec_direcs)
   call cla_get('-spin_channel', spin_channel)
   if(present(out_file_SC_kpts)) call cla_get('-out_sckpts_file', out_file_SC_kpts)
   spin_channel = abs(spin_channel)
    if(spin_channel == 0 .or. spin_channel > 2)then 
        if(spin_channel == 0) spin_channel = 1
        if(spin_channel > 2) spin_channel = 2
        write(*,'(A,I0,A)')'WARNING (get_commline_args): Resetting spin channel to ', spin_channel, '.'
    endif

    stop_if_not_commensurate = cla_key_present('-stop_if_not_commensurate')
    write_attempted_pc_corresp_to_input_pc = .not. cla_key_present('-skip_propose_pc_for_given_pc')
    write_attempted_pc_corresp_to_SC = .not. cla_key_present('-skip_propose_pc_for_given_sc')

end subroutine get_commline_args


end module options_setup
