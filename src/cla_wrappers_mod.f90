!! Copyright (C) 2014-2017 Paulo V. C. Medeiros
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
! MODULE: cla_wrappers 
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION:
!> Wrappers to the CLA-F90 module, used for command-line arguments.
!> See http://web.cecs.pdx.edu/~zaron/pub/CLA-F90.html
!==============================================================================

module cla_wrappers
use general_io
use cla
use kinds
use constants_and_types
implicit none
PRIVATE
PUBLIC :: get_commline_args, args

type(comm_line_args), save :: args

!! Functions and subroutines
CONTAINS


function valid_vector_from_cla(cla_key, default, null_allowed) result(rtn)
implicit none
real(kind=dp), dimension(1:3) :: rtn
character(len=*), intent(in) :: cla_key
real(kind=dp), dimension(1:3), intent(in) :: default
logical, intent(in), optional :: null_allowed

character(len=STRLEN) :: str_vec
integer :: ios, i
logical :: is_null, zero_allowed

    rtn = default
    zero_allowed = .TRUE.
    if(present(null_allowed)) zero_allowed = null_allowed
    if(cla_key_present(cla_key))then
        call  cla_get(cla_key, str_vec)
        ios = -1
        read(str_vec, *, iostat=ios) (rtn(i), i=1,3)
        if(ios /= 0) rtn = default
        is_null = sqrt(dot_product(rtn, rtn)) < epsilon(1.0_dp)
        if(is_null .and. .not. zero_allowed) rtn = default
    endif

end function valid_vector_from_cla


subroutine get_abinit_wf_file_from_files_file(args, iostat)
implicit none
type(comm_line_args), intent(inout) :: args
integer, intent(out), optional :: iostat
integer :: io_unit, ios
logical :: abinit_files_file_exists
character(len=256) :: abinit_out_base

    abinit_out_base = ''
    if(present(iostat)) iostat = -1
    inquire(file=args%abinit_files_file, exist=abinit_files_file_exists)
    if(.not. abinit_files_file_exists)then
        if(cla_key_present('-files_file'))then
            write(*,'(3A)')"ERROR: You've chosen ABINIT by passing the option &
                            '-files_file', but the files file " // &
                            new_line('A') // &
                            "       '", &
                            trim(adjustl(args%abinit_files_file)), &
                           "' has not been found."
        else
            write(*,'(A)')&
                "ERROR: You've chosen ABINIT, but you have neither given a &
                        path to a wavefunction file (option '-wf_file')," &
                        // new_line('A') // "       nor a 'files file' name &
                        (option '-files_file')."
        endif
    else
        io_unit = available_io_unit()
        open(unit=io_unit, file=args%abinit_files_file, action='read')
            read(io_unit, *)
            read(io_unit, *)
            read(io_unit, *)
            read(io_unit, *, iostat=ios) abinit_out_base
        close(io_unit) 
        if(ios==0)then
            args%wf_file = trim(adjustl(abinit_out_base)) // '_WFK'
        else
            write(*,'(3A)')"ERROR: Problems reading ABINIT's files file '", &
                            trim(adjustl(args%abinit_files_file)), "'."
        endif
        if(present(iostat)) iostat = ios
    endif

end subroutine get_abinit_wf_file_from_files_file

!==============================================================================
!> \ingroup changes_upon_new_interface
!> Takes care of command line arguments.
!!
!! @param[out] args Variable of the "comm_line_args" type that holds all args
!! @param[in] running_from_main_code Should be .TRUE. if running from main code
!==============================================================================
subroutine get_commline_args(args, running_from_main_code)
use iso_c_binding
implicit none
type(comm_line_args), intent(out) :: args
logical, intent(in), optional :: running_from_main_code
integer :: spin_channel, env_var_stat, ios
logical :: main_code

    main_code = .TRUE.
    if(present(running_from_main_code)) main_code = running_from_main_code 

    args%pw_code = ''
    args%qe_outdir = ''
    args%qe_prefix = ''
    args%abinit_files_file = ''
    args%castep_seed = ''
    args%WF_file = ''
    args%input_file_prim_cell = ''
    args%input_file_supercell = ''
    args%input_file_pc_kpts = ''
    args%input_file_energies = ''
    args%output_file_symm_averaged_EBS = ''
    args%output_file_only_user_selec_direcs = ''
    args%output_file_symm_averaged_unf_dens_op = ''
    args%output_file_only_user_selec_direcs_unf_dens_op = ''
    args%out_file_SC_kpts = ''
    args%origin_for_spin_proj_passed_in_rec = .TRUE.
    ! Optional command line arguments.
    call cla_init
    call cla_register(&
             key='-qe', &
             description='Parse wavefunctions generated by Quantum ESPRESSO.',&
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-abinit', &
             description='Parse wavefunctions generated by ABINIT.', &
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-castep', &
             description='Parse wavefunctions generated by CASTEP.', &
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-outdir', &
             description="Quantum ESPRESSO's outdir. Like in QE,"&
                 //new_line('A')// &
                 "     the default is the value of the environment variable" &
                 //new_line('A')// &
                 "     ESPRESSO_TMPDIR if set, or the current dir otherwise.",&
             kkind=cla_char, default='.' &
         )
    call cla_register(&
             key='-prefix', &
             description="Quantum ESPRESSO's 'prefix' used in your calc.", &
             kkind=cla_char, default='pwscf' &
         )
    call cla_register(&
             key='-files_file', &
             description="Path to ABINIT's 'files file' you used.", &
             kkind=cla_char, default='files.file' &
         )
    call cla_register(&
             key='-wf_file', & 
             description='Path to the wavefunction file.', &
             kkind=cla_char, default='WAVECAR' &
         )
    call cla_register(&
             key='-seed', &
             description='CASTEP seed for the ".orbitals" file.', &
             kkind=cla_char, default='' &
         )
    call cla_register(&
             key='-pc_file', &
             description='Path to the reference primitive cell file.', &
             kkind=cla_char, default='prim_cell_lattice.in' &
         )
    call cla_register(&
             key='-sc_file', &
             description='Path to the supercell file.', &
             kkind=cla_char, default='supercell_lattice.in' &
         )
    call cla_register(&
             key='-pckpts_file', &
             description='Path to the pc-kpts file.', &
             kkind=cla_char, default='KPOINTS_prim_cell.in' &
         )
    call cla_register(&
             key='-out_sckpts_file', &
             description='Path to the SC-kpts file (pre-unfolding utility)', &
             kkind=cla_char, default='KPOINTS_supercell.out' &
          )
    call cla_register(&
             key='-energy_file', &
             description='Path to the energy configuration file.', &
             kkind=cla_char, default='energy_info.in' &
         )
    call cla_register(&
             key='-out_file_nosymm', &
             description='Name of the output file (not symmetry-averaged).', & 
             kkind=cla_char, &
             default='unfolded_EBS_not-symmetry_averaged.dat' &
         )
    call cla_register(&
             key='-out_file_symm', &
             description='Name of the output file (symmetry-averaged).', &
             kkind=cla_char, &
             default='unfolded_EBS_symmetry-averaged.dat' &
         )
    call cla_register(&
             key='-write_unf_dens_op', &
             description='Write the unfolding-density operator to a file.',&
             kkind=cla_flag, default='F' &
          )
    call cla_register(&
             key='-output_file_symm_averaged_unf_dens_op', &
             description='File where the symmetry-averaged '//&
                 'unfolding-density operator will be saved (if requested).', &
             kkind=cla_char, &
             default='unfolding_density_operator_symm_avgd.dat' &
         )
    call cla_register(&
             key='-output_file_only_user_selec_direcs_unf_dens_op', &
             description='File where the non-symm avgd unfolding-density '//&
                 'operator will be saved (if requested).', &
             kkind=cla_char, &
             default='unfolding_density_operator_not-symm_avgd.dat' &
         )
    call cla_register(&
             key='-spin_channel', &
             description='Either 1 or 2 &
                 (if the wavefunction has a 2nd spin channel).', &
             kkind=cla_int, default='1' &
         )
    call cla_register(&
             key='-n_sckpts_to_skip', &
             description='Tells BandUP to skip the first "n_sckpts_to_skip" &
                 KPTS'//new_line('A')//&
                 '     found in the wavefunction file.'&
                 //new_line('A')//'     This can be used, e.g., to do &
                 unfolding with hybrid functional'//new_line('A')//&
                 '     calcs. using VASP.', &
                 kkind=cla_int, &
                 default='0' &
         )
    call cla_register(&
             key='-continue_if_not_commensurate', &
             description='Continue if the SC and pc are not commensurate.', &
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-skip_propose_pc_for_given_pc', &
             description='Do not attempt to propose a pc for the input &
                 reference unit cell.', &
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-skip_propose_pc_for_given_sc', &
             description='Do not attempt to propose a pc for the input SC', &
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-no_symm_sckpts', &
             description="Do not use symmetry to reduce the number of &
                 SC K-points needed." // new_line('A') // &
                 "     Works only if used *both* with BandUP's main &
                 unfolding task *and* when"// new_line('A')//&
                 "     getting the needed SC-Kpts (prior to unfolding).",&
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-no_symm_avg', &
             description="Do not calculate the symmetry-averaged EBS. &
                 If the SC-Kpts needed for unfolding"//new_line('A')//&
                 "     were determined using this option, &
                 then BandUP's main unfolding"//new_line('A')// &
                 "     task *must* be run with this flag as well.", &
                 kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-saxis', &
             description='Quantization axis; applies to noncollinear &
                 calculations. Still under test. ' //new_line('A')//&
                 '     Only the default works at the moment.', &
             kkind=cla_char, default='0, 0, 1' &
          )
    call cla_register(&
             key='-normal_to_proj_plane', &
             description='Applies only to noncollinear calculations. &
                 The projection axis used for'//new_line('A')//&
                 '     "#Spin _|_ k" will be also _|_ to &
                 "normal_to_proj_plane".' //new_line('A')//&
                 '     Usage example: -normal_to_proj_plane "0 0 1"'&
                 //new_line('A'), &
             kkind=cla_char, default='0, 0, 1' &
         )
    call cla_register(&
             key='-origin_for_spin_proj_rec', &
             description='', &
             kkind=cla_char, default='0, 0, 0' &
         )
    call cla_register(&
             key='-origin_for_spin_proj_cart', &
             description='', &
             kkind=cla_char, default='0, 0, 0' &
         )
    call cla_register(&
             key='--continue_if_npw_smaller_than_expected', &
             description='Ignore error if the estimated npws are smaller than &
                 the ones found.', &
             kkind=cla_flag, default='F' &
         )
    call cla_register(&
             key='-dont_unfold', &
             description='Return folded bands (kpts/eigenvals from &
                 wavefunction file). Mostly for tests.', &
             kkind=cla_flag, default='F' &
         )

    call cla_validate

    call cla_get('-wf_file', args%WF_file)
    ! Setting up BandUP for VASP, Quantum Espresso, ABINIT or CASTEP
    ! The default is VASP.
    if(cla_key_present('-qe') .or. cla_key_present('-outdir') .or. &
       cla_key_present('-prefix'))then
#       if defined (__QE_SUPPORT__)
            args%pw_code = 'qe'
            call cla_get('-prefix', args%qe_prefix)
            if(cla_key_present('-outdir'))then
                call cla_get('-outdir', args%qe_outdir)
            else
                call get_environment_variable('ESPRESSO_TMPDIR', &
                                              args%qe_outdir, &
                                              status=env_var_stat)
                if(env_var_stat/=0) args%qe_outdir = '.'
            endif
#       else
            write(*,'(A)')&
                'ERROR: You have compiled BandUP without support for QE!'
            write(*,'(A)')'Cannot continue. Stopping now.'
            stop
#       endif
    else if(cla_key_present('-abinit') .or. cla_key_present('-files_file'))then
        args%pw_code = 'abinit'
        if(main_code)then
            call cla_get('-files_file', args%abinit_files_file)
            if(.not. cla_key_present('-wf_file'))then
                call get_abinit_wf_file_from_files_file(args, ios)
                if(ios/=0)then
                    write(*,'(A)')&
                        'Cannot determine wavefunction file. Stopping now.'
                    stop
                endif
            endif
        endif
    else if(cla_key_present('-castep'))then
#       if defined (__CASTEP_SUPPORT__)
            args%pw_code = 'castep'
            if(cla_key_present('-seed'))then
                call cla_get('-seed', args%castep_seed)
                args%WF_file = trim(adjustl(args%castep_seed)) // '.orbitals'
            else
                write(*,'(A)')&
                    'You must pass a CASTEP seed using "-seed". Stopping now.'
                stop
            endif
#       else
            write(*,'(A)')&
                'ERROR: You have compiled BandUP without support for CASTEP!'
            write(*,'(A)')'Cannot continue. Stopping now.'
            stop
#       endif
    else
        args%pw_code = 'vasp'
    endif

    call cla_get('-pc_file', args%input_file_prim_cell)
    call cla_get('-sc_file', args%input_file_supercell)
    call cla_get('-pckpts_file', args%input_file_pc_kpts)
    call cla_get('-energy_file', args%input_file_energies)
    call cla_get('-out_file_symm', args%output_file_symm_averaged_EBS)
    call cla_get('-out_file_nosymm', args%output_file_only_user_selec_direcs)
    call cla_get('-output_file_symm_averaged_unf_dens_op', &
                  args%output_file_symm_averaged_unf_dens_op)
    call cla_get('-output_file_only_user_selec_direcs_unf_dens_op', &
                  args%output_file_only_user_selec_direcs_unf_dens_op)
    call cla_get('-out_sckpts_file', args%out_file_SC_kpts)
    call cla_get('-n_sckpts_to_skip', args%n_sckpts_to_skip)

    call cla_get('-spin_channel', spin_channel)
    spin_channel = abs(spin_channel)
    if(spin_channel < 1 .or. spin_channel > 2)then 
        if(spin_channel < 1) spin_channel = 1
        if(spin_channel > 2) spin_channel = 2
        write(*,'(A,I0,A)')&
            'WARNING (get_commline_args): Resetting spin channel to ', &
            spin_channel, '.'
    endif
    args%spin_channel = spin_channel

    args%stop_if_not_commensurate = &
        .not. cla_key_present('-continue_if_not_commensurate')
    args%write_attempted_pc_corresp_to_input_pc = &
        .not. cla_key_present('-skip_propose_pc_for_given_pc')
    args%write_attempted_pc_corresp_to_SC = &
        .not. cla_key_present('-skip_propose_pc_for_given_sc')
    args%no_symm_avg = cla_key_present('-no_symm_avg')
    args%no_symm_sckpts = cla_key_present('-no_symm_sckpts')

    args%saxis = valid_vector_from_cla('-saxis', real((/0, 0, 1/), kind=dp), &
                                       null_allowed=.FALSE.)
    args%saxis = real((/0, 0, 1/), kind=dp) ! To be removed in the future

    args%normal_to_proj_plane = valid_vector_from_cla(&
                                    '-normal_to_proj_plane', &
                                    real((/0, 0, 1/), kind=dp), &
                                    null_allowed=.FALSE. &
                                )

    if(cla_key_present('-origin_for_spin_proj_cart'))then
        args%origin_for_spin_proj_cartesian = &
            valid_vector_from_cla(&
                '-origin_for_spin_proj_cart', &
                real((/0, 0, 0/), kind=dp) &
            )
        args%origin_for_spin_proj_passed_in_rec = .FALSE.
    endif
    if(cla_key_present('-origin_for_spin_proj_rec'))then
        args%origin_for_spin_proj_rec = &
            valid_vector_from_cla(&
                '-origin_for_spin_proj_rec', &
                real((/0, 0, 0/), kind=dp) &
            )
        args%origin_for_spin_proj_passed_in_rec = .TRUE.
        if(cla_key_present('-origin_for_spin_proj_cart'))then
            write(*,'(A)')&
                'WARNING: You passed both "-origin_for_spin_proj_rec" and &
                "-origin_for_spin_proj_cart" as parameters.'
            write(*,'(A)')&
                '         Only "-origin_for_spin_proj_rec" will be used.'
        endif
    endif


    args%perform_unfold = .not. cla_key_present('-dont_unfold')
    args%write_unf_dens_op = cla_key_present('-write_unf_dens_op')
    args%continue_if_npw_smaller_than_expected = &
        cla_key_present('--continue_if_npw_smaller_than_expected')

end subroutine get_commline_args



end module cla_wrappers
