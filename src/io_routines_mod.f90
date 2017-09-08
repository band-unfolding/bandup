!! Copyright (C) 2013-2017 Paulo V. C. Medeiros
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
! MODULE: io_routines
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION:
!> Contains most of the specific I/O routines used by BandUP.
!==============================================================================

module io_routines
use constants_and_types
use cla_wrappers
use strings
use general_io
use read_abinit_wavefunctions
use read_vasp_wavecar
use read_vasp_files
use write_vasp_files
use time, only: time_now, formatted_time
use math, only: norm, coords_cart_vec_in_new_basis, cross, same_vector, &
                n_digits_integer
use lists_and_seqs, only: list_index
use units, only: to_angstrom
# if defined (__QE_SUPPORT__)
use qexml_module
use read_qe_wavefunctions
# endif
# if defined (__CASTEP_SUPPORT__)
use read_castep_wavefunctions
# endif
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: print_welcome_messages,print_message_commens_test, &
          read_energy_info_for_band_search, &
          print_message_success_determining_GUR, &
          print_geom_unfolding_relations, &
          print_last_messages_before_unfolding, &
          get_list_of_SCKPTS, read_pckpts_selected_by_user, &
          read_wavefunction, &
          print_symm_analysis_for_selected_pcbz_dirs, &
          say_goodbye_and_save_results, & 
          print_final_times, write_band_struc, &
          print_message_pckpt_cannot_be_parsed, &
          write_attempted_pc_assoc_with_input_unit_cell_and_SC


CONTAINS 


subroutine print_welcome_messages(package_version)
implicit none
character(len=*), intent(in), optional :: package_version
character(len=127) :: git_hash, git_branch, aux_source_info_str

git_hash = get_git_info_compiled_files('hash_latest_commit')
git_branch = get_git_info_compiled_files('branch_name')

write(*,'(A,/,A)') &
'===========================================================================',&
'        BandUP: Band Unfolding code for Plane-wave based calculations      '
if(present(package_version))then
    aux_source_info_str="    "//trim(adjustl(package_version))
    if(len(trim(adjustl(git_branch)))>0)then
        aux_source_info_str = trim(aux_source_info_str) // NEW_LINE('A') // &
                              '    From Git branch: "' // &
                              trim(adjustl(git_branch)) // '"'
    endif 
    write(*,'(A)') trim(aux_source_info_str)
endif
    write(*,'(A)') &
"Compiled using "//trim(adjustl(compiler_version()))//&
" on "//trim(adjustl(compilation_time()))
write(*,'(17(A,/),A)') &
'===========================================================================',&
'Copyright (C) 2013-2017 Paulo V. C. Medeiros                               ',&
'                                                                           ',&
'                        Computational Physics Division                     ',&
'                        Department of Physics, Chemistry and Biology - IFM ',&
'                        Linköping University                               ',&
'                        Sweden                                             ',&
'                                                                           ',&
'                        Current address:                                   ',&
'                        University of Cambridge                            ',&
'                        Theory of Condensed Matter (TCM) Group             ',&
'                        Department of Physics                              ',&
'                        Cavendish Laboratory                               ',&
'                        Cambridge, UK                                      ',&
'                        Email: pvm20@cam.ac.uk                             ',&
'                                                                           ',&
'Please visit www.ifm.liu.se/theomod/compphys/band-unfolding                ',&
'==========================================================================='
write(*,*)

end subroutine print_welcome_messages


subroutine print_message_commens_test(commensurate,M,stop_if_not_commens)
implicit none
logical, intent(in) :: commensurate, stop_if_not_commens
real(kind=dp), dimension(1:3,1:3), intent(in) :: M
character(len=str_len) :: message_header,message_footer,message_footer2, &
                          str_n_char_float_format,str_n_decimals,format_string
character(len=str_len), dimension(1:3) :: float_format
integer :: i, j, n_decimals, max_n_digits_before_dec_point, n_digits

if(stop_if_not_commens)then
    message_header = &
        '                                     ERROR!!                       '
    message_footer = &
        ' Stopping now.'
else
    message_header = &
        '                               >>> WARNING!! <<<                   '
    message_footer = &
        ' >>> The results might not be what you expect. &
        Continue only if you are really sure. '
    message_footer2 = &
        "     My advice: Always use commensurate SC and PC, &
        and never continue if they're not."
endif

write(*,*)''
if(.not. commensurate)then
    write(*,*)''
    write(*,'(10(A,/))') &
'===========================================================================',&
                              message_header                                 ,&
'===========================================================================',&
'  The SC and the reference PC that you have chosen are not commensurate!  ', &
'  The choice of the reference PC vectors is very important and can change ', &
'  a lot the results of unfolding. It is particularly important to always  ', &
'  verify that the SC and the reference PC are commensurate. Otherwise, the', &
'  calculated spectral weigths will most likely be very small, which will  ', &
'  certainly cause the values of the unfolded delta_Ns to be also small.   ', &
'==========================================================================='
    write(*,*)''
endif

if(commensurate)then
    n_decimals = 1
else
    n_decimals = 6
endif
write(str_n_decimals,*) n_decimals

do j=1,3 ! Allowing for a different format for each column
    max_n_digits_before_dec_point = 0
    do i=1,3
        n_digits = n_digits_integer(nint(M(i,j)), add_one_if_negative=.TRUE.)
        if(n_digits > max_n_digits_before_dec_point)then
            max_n_digits_before_dec_point = n_digits
        endif
    enddo
    write(str_n_char_float_format,*) &
        max_n_digits_before_dec_point + n_decimals + 1    
    float_format(j) = 'f' // trim(adjustl(str_n_char_float_format)) // '.' // &
                      trim(adjustl(str_n_decimals))
enddo
format_string = '(2(A,/),3(A,' // adjustl(trim(float_format(1))) // &
                         ',A,' // adjustl(trim(float_format(2))) // &
                         ',A,' // adjustl(trim(float_format(3))) // &
                         ',A,/))'
write(*,format_string) &
' Relation between the latt vecs of the chosen SC and reference unit cell:', &
'                                                                         ', &
'      A[1] = (',M(1,1),')*a[1] + (',M(1,2),')*a[2] + (',M(1,3),')*a[3]   ', & 
'      A[2] = (',M(2,1),')*a[1] + (',M(2,2),')*a[2] + (',M(2,3),')*a[3]   ', & 
'      A[3] = (',M(3,1),')*a[1] + (',M(3,2),')*a[2] + (',M(3,3),')*a[3]   '
if(commensurate)then
    write(*,'(A)') &
    '      * The SC and the reference unit cell are commensurate. Good!'
else
    write(*,*)''
    write(*,'(A)')message_footer
    write(*,'(A)')message_footer2
    write(*,'(A)') &
'==========================================================================='
endif
write(*,*)''

end subroutine print_message_commens_test


subroutine write_attempted_pc_assoc_with_input_unit_cell_and_SC(&
               crystal_pc, crystal_SC &
           )
implicit none
type(crystal_3D), intent(in) :: crystal_pc, crystal_SC

    write(*,*)
    write(*,'(A)')' Checking if you are working with the smallest possible SC:'
    if(crystal_SC%is_prim_cell)then
        write(*,'(A)') '     * OK!'
    else
        write(*,'(3(A,/), A)') &
          '     * Your SC seems to be a perfect repetition of a smaller PC!', &
          "       > That's OK if you really want it, but using the smallest", &
          "         possible SC will certainly save you some time!", &
          '     * We attempted to find a possible PC associated with it.'
        write(*,'(3A)')'       Please take a look at the file "', &
                       trim(adjustl(file_for_SC_reduced_to_prim_cell)),'".'
        call write_crystal_to_file(&
                 crystal_SC%corresp_pc, file_for_SC_reduced_to_prim_cell &
             )
    endif

    write(*,*)
    write(*,'(A)')&
        ' Checking if working with the smallest possible ref. unit cell:'
    if(crystal_pc%is_prim_cell)then
        write(*,'(A)') '     * OK!'
    else
        call write_crystal_to_file(&
                 crystal_pc%corresp_pc, file_for_pc_reduced_to_prim_cell &
             )
        write(*,'(A)')&
            "     * Your reference unit cell doesn't seem to be a prim. cell!"
        write(*,'(A)')&
            "       > That's fine if you really want it."
        write(*,'(A)')&
            '       > Anyway, we attempted to find a PC associated with it.'
        write(*,'(3A)')&
            '         Please take a look at the file "', &
            trim(adjustl(file_for_pc_reduced_to_prim_cell)),'".'
    endif

end subroutine write_attempted_pc_assoc_with_input_unit_cell_and_SC


subroutine print_geom_unfolding_relations(&
               GUR, list_SC_kpts_in_wavecar, crystal_pc, crystal_SC &
           )
implicit none
!! Geometric Unfolding Relations
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR
type(vec3d), dimension(:), intent(in) :: list_SC_kpts_in_wavecar
type(crystal_3D), intent(in) :: crystal_pc, crystal_SC
integer :: nkpts, i_SCKPT, i_selc_pcbz_dir,i_needed_dirs,ipc_kpt
integer, dimension(:), allocatable :: n_pckpts_dirs
real(kind=dp), dimension(1:3) :: current_SCKPT, actual_folding_SCKPT, pc_kpt, &
                                 folding_G
real(kind=dp), dimension(1:3,1:3) :: b_matrix_pc,B_matrix_SC


    b_matrix_pc = crystal_pc%rec_latt_vecs
    B_matrix_SC = crystal_SC%rec_latt_vecs
    nkpts = size(list_SC_kpts_in_wavecar)
    allocate(n_pckpts_dirs(1:size(GUR%SCKPT(1)%selec_pcbz_dir(:))))
    do i_selc_pcbz_dir=1,size(GUR%SCKPT(1)%selec_pcbz_dir(:))
        n_pckpts_dirs(i_selc_pcbz_dir) = &
            size(&
                GUR%SCKPT(1)%selec_pcbz_dir(i_selc_pcbz_dir)%&
                    needed_dir(1)%pckpt(:) &
            )
    enddo

    write(*,"(A)")"Summary of the points checked on the SCBZ and the pcbz:"
    write(*,'(A,I0)')'    * Total # of SC-KPTS found: ',nkpts
    if((GUR%n_pckpts - sum(n_pckpts_dirs)) > 0)then
        write(*,'(A,I0)')&
            '    * Total # of pc-kpts requested: ', sum(n_pckpts_dirs)
        write(*,'(A,I0)')&
            '    * Total # of complementary pc-kpts determined by the &
            symmetry analysis: ', GUR%n_pckpts - sum(n_pckpts_dirs)
    endif
    write(*,'(A,I0)')'    * Total # of pc-kpts to be checked: ',GUR%n_pckpts
    write(*,*)
    write(*,"(A)")"Summary of the geometric folding relations:"
    write(*,'(2(A,I0),A)')&
        '    * A total of ',GUR%n_folding_pckpts,' pc-kpts (out of the ', &
        GUR%n_pckpts,' checked) satisfy the folding condition:'
    do i_SCKPT=1,nkpts
        do i_selc_pcbz_dir=1,size(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(:))
            do i_needed_dirs=1, &
               size(&
                   GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)% &
                       needed_dir(:) &
               )
                do ipc_kpt=1, &
                   size(&
                       GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)% &
                           needed_dir(i_needed_dirs)%pckpt(:) &
                   )
                    if(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)%&
                           needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds)then
                        current_SCKPT = &
                            list_SC_kpts_in_wavecar(i_SCKPT)%coord(:)
                        actual_folding_SCKPT = &
                            GUR%SCKPT(i_SCKPT)%&
                                selec_pcbz_dir(i_selc_pcbz_dir)%&
                                needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%&
                                coords_actual_unfolding_K
                        pc_kpt(:) = GUR%SCKPT(i_SCKPT)%&
                                        selec_pcbz_dir(i_selc_pcbz_dir)%&
                                        needed_dir(i_needed_dirs)%&
                                        pckpt(ipc_kpt)%coords(:)
                        folding_G(:) = pc_kpt(:) - actual_folding_SCKPT(:)
                        call print_message_pckpt_folds(&
                                 pc_kpt, current_SCKPT, i_SCKPT, &
                                 actual_folding_SCKPT, folding_G, &
                                 b_matrix_pc, B_matrix_SC &
                             )
                    endif
                enddo 
            enddo
        enddo
    enddo
    write(*,'(A)')''
    write(*,'(A)')''

end subroutine print_geom_unfolding_relations


subroutine read_energy_info_for_band_search(&
               input_file, e_fermi, E_start, E_end, delta_e &
           )
implicit none
character(len=*), intent(in) :: input_file
real(kind=dp), intent(out) :: e_fermi,E_start,E_end,delta_e
integer :: unt, ios
real(kind=dp) :: E_start_minus_ef, E_end_minus_ef, aux_ener, dE_factor
character(len=str_len) :: str_delta_e, aux_str
character(len=1) :: first_char
logical :: verbose

    verbose = .TRUE.
    dE_factor = 0.4E-2_dp
    unt = available_io_unit()
    open(unit=unt, file=input_file)
        first_char = '#'
        do while (first_char=='#')
            read(unt,*) aux_str
            read(aux_str, *) first_char
        enddo
        read(aux_str,*) e_fermi
        read(unt,*) E_start_minus_ef
        read(unt,*) E_end_minus_ef
        read(unt,*) str_delta_e
    close(unt)

    E_start = E_start_minus_ef + e_fermi
    E_end = E_end_minus_ef + e_fermi
    if(E_end < E_start)then
        aux_ener = E_end
        E_end = E_start
        E_start = aux_ener
    endif

    read(str_delta_e,*,iostat=ios) delta_e
    if(ios==0)then
        delta_e = dabs(delta_e)
    else
        delta_e = dabs(dE_factor*(E_end - E_start))
        if(verbose)then
            if(upper_case(str_delta_e(1:1)) /= 'A')then
                write(*,'(A,E9.2,A)')&
                    'WARNING: Could not read the size of the energy &
                    intervals. Using the default dE = ', dE_factor, &
                    '*(Emax - Emin).'
            else
                write(*,'(A,E9.2,A)')&
                    'Automatically setting the size of the energy intervals &
                    to dE = ', dE_factor, '*(Emax - Emin).' 
            endif
        endif
    endif

end subroutine read_energy_info_for_band_search


subroutine read_unit_cell(input_file,latt,latt_vol)
implicit none
character(len=*), intent(in) :: input_file
real(kind=dp), dimension(1:3,1:3), intent(out) :: latt
real(kind=dp), intent(out), optional :: latt_vol
real(kind=dp) :: scaling_factor
integer :: ivec, icomp

    open(unit=01, file=input_file)
        read(01,*)
        read(01,*)scaling_factor
        do ivec=1,3
            read(01,*) (latt(ivec,icomp), icomp=1,3)
            latt(ivec,:) = scaling_factor*latt(ivec,:)
        enddo
    close(01)
    if(present(latt_vol))then
        latt_vol = dabs(dot_product(latt(1,:),cross(latt(2,:),latt(3,:))))
    endif

end subroutine read_unit_cell


!==============================================================================
!> \ingroup changes_upon_new_interface
!> Prints the last messages before the unfolding process starts.
!> If you introduce an interface to a new ab initio code, and don't modify this
!! routine, no specific information about the ab initio code will be shown, but
!! BandUP will work just fine.
!==============================================================================
subroutine print_last_messages_before_unfolding(&
               args, list_of_SCKPTS, B_matrix_SC, vbz, &
               E_start, E_end, delta_e, e_fermi, spinor_wf &
           )
implicit none
type(comm_line_args), intent(in) :: args
type(vec3d), dimension(:), intent(in) :: list_of_SCKPTS
real(kind=dp), dimension(1:3,1:3), intent(in) :: B_matrix_SC
real(kind=dp), intent(in) :: vbz,E_start,E_end,delta_e,e_fermi
logical, intent(in) :: spinor_wf
real(kind=dp) :: file_size_in_MB, file_size_in_GB, file_size, VSBZ, &
                 approx_mem_per_kpt_in_bytes, mem_per_kpt_in_MB, &
                 mem_per_kpt_in_GB, mem_per_kpt
integer :: nkpts
integer(long_int_kind) :: file_size_in_bytes
character(len=2) :: file_size_units, mem_per_kpt_units
logical :: using_omp, warn_lack_omp, print_using_omp_msg

    nkpts = size(list_of_SCKPTS)
    select case(lower_case(trim(args%pw_code)))
        case('vasp', 'abinit', 'castep') 
            call get_file_size_in_bytes(file_size_in_bytes, file=args%WF_file)
            file_size_in_MB = real(file_size_in_bytes,kind=dp)/(2.0**20)
            file_size_in_GB = real(file_size_in_bytes,kind=dp)/(2.0**30)
            file_size = file_size_in_GB
            file_size_units = 'GB'
            if(file_size_in_GB < 1.0_dp)then
                file_size = file_size_in_MB
                file_size_units = 'MB'
            endif

            approx_mem_per_kpt_in_bytes = &
                real(file_size_in_bytes,kind=dp)/real(nkpts,kind=dp)
            mem_per_kpt_in_MB = approx_mem_per_kpt_in_bytes/(2.0**20)
            mem_per_kpt_in_GB = approx_mem_per_kpt_in_bytes/(2.0**30)
            mem_per_kpt = mem_per_kpt_in_GB
            mem_per_kpt_units = 'GB'
            if(mem_per_kpt_in_GB < 1.0_dp)then
                mem_per_kpt = mem_per_kpt_in_MB
                mem_per_kpt_units = 'MB'
            endif
            write(*,*)
            write(*,'(A,1X,A,1X,A)')'Parsing', &
                                    upper_case(trim(args%pw_code)), &
                                    'wavefunctions.'
            write(*,'(3A)')'    * File = "', trim(adjustl(args%WF_file)),'"'
            if(lower_case(trim(args%pw_code))=='abinit')then
                write(*,'(3(A,/), A)')&
'=========================================================================', &
'    The interface to ABINIT is still under test. Use it with caution.    ', &
'                    Feedback is appreciated!!                            ', &
'========================================================================='
    else if(lower_case(trim(args%pw_code))=='castep')then
        write(*,'(3(A,/), A)')&
'=========================================================================', &
'    The interface to CASTEP is still under test. Use it with caution.    ', &
'                    Feedback is appreciated!!                            ', &
'========================================================================='
            endif
            write(*,'(A,f0.2,1X,2A)')&
                'The wavefunction file is ', &
                file_size,file_size_units, &
                ' big. Only the necessary data will be read.'
            write(*,'(A,f0.2,1X,2A)')&
                '    * Max. of approx. ', &
                mem_per_kpt,mem_per_kpt_units, &
                ' at a time.'

        case('qe')
            write(*,*)
            write(*,'(A)')'Parsing Quantum ESPRESSO (pwscf) wavefunctions.'
            write(*,'(5A)')'Using prefix = "', trim(args%qe_prefix), &
                          '" and outdir = "', trim(args%qe_outdir),'".'
    end select
    write(*,*)

    VSBZ = dabs(&
               dot_product(&
                   B_matrix_SC(1,:), cross(B_matrix_SC(2,:), B_matrix_SC(3,:))&
               ) &
           )
    if(VSBZ < vbz)then
        write(*,'(A,f0.2,A)')&
            'The used supercell is ', vbz/VSBZ, ' times bigger than &
            the primitive cell.'
    else
        write(*,'(A)')&
            'WARNING: The informed primitive cell is not smaller than &
            the supercell found in the WAVECAR file.'
        write(*,'(A,f0.8,A)')'         vpc = ',VSBZ/vbz,' VSC'
    endif

    write(*,'(3(A,f0.5),A)')&
        'Working within the energy interval ', E_start-e_fermi, &
        ' < E-EF < ', E_end-e_fermi, ' in increments of ', delta_e, ' eV'
    write(*,'(A,f0.5,A)')&
        'The Fermi energy EF = ', e_fermi, &
        ' will be set to the zero of the energy scale.'

    select case (trim(adjustl(upper_case(args%pw_code))))
        case default
            warn_lack_omp = .FALSE.
            print_using_omp_msg = .FALSE.
        case('VASP') 
            warn_lack_omp = .TRUE.
            print_using_omp_msg = .TRUE.
    end select

    using_omp = .FALSE.
    write(*,*)
    !$ using_omp = .TRUE.
    if(using_omp .and. print_using_omp_msg)then
        write(*,'(A)')'Some parts of BandUP have been parallelized with &
                       OpenMP and will be running using'
        write(*,'(A,I0,A)')'a maximum of ',omp_get_max_threads(),' thread(s).'
        write(*,'(A)')"You can choose the maximum number of threads by &
                       setting the environment variable"
        write(*,'(A)')'"OMP_NUM_THREADS". Ex.: export OMP_NUM_THREADS=4'
        continue
    endif
    if(.not. using_omp .and. warn_lack_omp)then
        write(*,'(5(A,/), A)')&
            "WARNING: BandUP seems to have been compiled without support to", &
            "         OpenMP shared-memory parallelization. While this is an",&
            "         optional feature and won't change your results, using ",&
            "         OpenMP will certainlt save you some time. If you have ",&
            "         not modified the provided Makefile and 'build' script,",&
            "         please contact your system's admin about this.        "
    endif

    write(*,*)
    !! The following message must always be written to the standard output in
    !! any version of BandUP, no matter whether modified or not
    write(*,'(A)') &      
'==========================================================================='
    write(*,'(8(A,/),A)') &
"NOTE: If you use BandUP or any modified/adapted version/part of it, you    ",&
"      SHOULD EXPLICITLY ACKNOWLEDGE THE USE OF BandUP IN YOUR PUBLICATIONS.",&
"      You should also read and cite                                        ",&
"       ------------------------------------------------------              ",&
"  >>> │ Paulo V. C. Medeiros, Sven Stafström and Jonas Björk |             ",&
"  >>> │ Phys. Rev. B 89, 041407(R) (2014)                    |             ",&
"  >>> │ http://dx.doi.org/10.1103/PhysRevB.89.041407         |             ",&
"       ------------------------------------------------------              ",&
"      and the appropriate references therein.                              "
    if(spinor_wf .or. args%write_unf_dens_op)then
        write(*,*)
        if(spinor_wf)then
            write(*,'(3(A,/))') &
"      Additionally, since your calculation involves spinor eigenstates ",&
"      (SOC, noncollinear magnetism), and BandUP uses the unfolding-density", &
"      operator formalism in this case, you should as well read and cite"
        elseif(args%write_unf_dens_op)then
            write(*,'(2(A,/))') &
"      Additionally, as you have requested the unfolding-density operators",&
"      to be calculated and written out, you should as well read and cite"
        endif
        write(*,'(5(A,/),A)') &
"       ------------------------------------------------------------------- ",&
"  >>> | Paulo V. C. Medeiros, S. S. Tsirkin, Sven Stafström, Jonas Björk  |",&
"  >>> | Phys. Rev. B 91, 041116(R) (2015)                                 |",&
"  >>> | http://dx.doi.org/10.1103/PhysRevB.91.041116                      |",&
"       ------------------------------------------------------------------- ",&
"      and the appropriate references therein.                              "
    endif
    write(*,'(A)') &      
'==========================================================================='
    !! End of message

    write(*,*)
    write(*,'(A)')'Unfolding starts now.'
    write(*,*)

end subroutine print_last_messages_before_unfolding


subroutine print_message_pckpt_folds(&
               pc_kpt, SCKPT, i_SCKPT, actual_folding_SCKPT, &
               folding_G,b_matrix_pc,B_matrix_SC &
           )
implicit none
integer, intent(in) :: i_SCKPT
real(kind=dp), dimension(1:3), intent(in) :: pc_kpt,SCKPT, &
                                             actual_folding_SCKPT, folding_G
real(kind=dp), dimension(1:3, 1:3), intent(in) :: b_matrix_pc,B_matrix_SC
real(kind=dp), dimension(1:3) :: pckpt_coords, coords
logical :: using_symm_to_get_coeffs
character(len=10) :: str_SCKPT_number

    using_symm_to_get_coeffs = &
        .not. same_vector(actual_folding_SCKPT(:), SCKPT)    
    write(*,*)
    write(*,*)
    str_SCKPT_number = ''
    if(.not. using_symm_to_get_coeffs)then
        write(str_SCKPT_number,"(A,I0,A)")'(',i_SCKPT,')'
    endif
    coords(:) = coords_cart_vec_in_new_basis(&
                    cart_vec=actual_folding_SCKPT(:), &
                    new_basis=B_matrix_SC &
                )
    write(*,'(A,3(f9.5,A))')&
        '    SCBZ wave-vector K'//trim(adjustl(str_SCKPT_number))//' = ',  &
        coords(1),'*B1 + ',coords(2),'*B2 +',coords(3),'*B3'
        
    pckpt_coords(:) = coords_cart_vec_in_new_basis(&
                          cart_vec=pc_kpt(:), &
                          new_basis=b_matrix_pc &
                      )
    write(*,'(A,3(f9.5,A))')&
        '    unfolds onto pcbz wave-vector k = ',  &
        pckpt_coords(1), '*b1 + ', pckpt_coords(2),'*b2 + ', &
        pckpt_coords(3), '*b3 = ' 
    pckpt_coords(:) = coords_cart_vec_in_new_basis(&
                          cart_vec=pc_kpt(:), &
                          new_basis=B_matrix_SC &
                      )
    write(*,'(A,3(f9.5,A))')&
        '                                    = ',  &
        pckpt_coords(1), '*B1 + ', pckpt_coords(2), '*B2+ ', &
        pckpt_coords(3), '*B3'

    coords(:) = coords_cart_vec_in_new_basis(&
                    cart_vec=folding_G, &
                    new_basis=B_matrix_SC &
                )
    write(*,'(A,A,3(I0,A))')&
        '    with the unfolding vector ',  &
        'G = ', nint(coords(1)), '*B1 + ', nint(coords(2)), '*B2 + ', &
        nint(coords(3)),'*B3.'
    if(using_symm_to_get_coeffs)then
        coords(:) = coords_cart_vec_in_new_basis(&
                        cart_vec=SCKPT(:), &
                        new_basis=B_matrix_SC &
                    )
        write(*,'(A,/,A,I0,A,3(f9.5,A),/,A,I0,A)') &
            '        *This SCBZ wave-vector belongs to the star of &
            the SCBZ wave-vector',  &
            '                 K(',i_SCKPT,') = ', coords(1),'*B1 + ', &
            coords(2),'*B2 +', coords(3),'*B3.', &
            '         The plane-wave coefficients for K(', i_SCKPT, &
            ') will be used.'
    endif

end subroutine print_message_pckpt_folds


!==============================================================================
!> \ingroup changes_upon_new_interface
!> Gets a list of the SC-KPTs present in the wavefunction file.
!> Similar routines are needed to implement interfaces to other PW codes.
!==============================================================================
subroutine get_SCKPTS_contained_in_wavecar(&
               list_SC_kpts_in_wavecar, args, crystal_SC &
           )
implicit none
type(vec3d), dimension(:), allocatable, intent(out) :: list_SC_kpts_in_wavecar
type(comm_line_args), intent(in) :: args
type(crystal_3D), intent(in) :: crystal_SC
integer(kind=selected_int_kind(18)) ::  nrecl
integer :: i_SCKPT, alloc_stat, i, input_file_unit, iost, nband, irec, nkpts
real(kind=dp) :: xnwk, xnband, xnplane
real(kind=dp), dimension(1:3) :: SCKPT_coords
real(kind=dp), dimension(1:3,1:3) :: B_matrix_SC

    call get_wavecar_record_lenght(nrecl, args%WF_file, iost)
    B_matrix_SC = crystal_SC%rec_latt_vecs
    ! Reading WAVECAR file
    input_file_unit = available_io_unit()
    open(unit=input_file_unit,file=args%WF_file,access='direct', &
         recl=nrecl,iostat=iost,status='old')
        read(unit=input_file_unit,rec=2) xnwk,xnband
        nband=nint(xnband) ! Number of bands
        nkpts = nint(xnwk)

        deallocate(list_SC_kpts_in_wavecar,stat=alloc_stat)
        allocate(list_SC_kpts_in_wavecar(1:nkpts))
        do i_SCKPT=1,nkpts
            ! Positioning the register at the correct k-point
            irec = 2 + (i_SCKPT-1)*(1+nband) + 1 
            read(unit=input_file_unit,rec=irec) xnplane,(SCKPT_coords(i),i=1,3)
            list_SC_kpts_in_wavecar(i_SCKPT)%coord(:) = &
                SCKPT_coords(1)*B_matrix_SC(1,:) + &
                SCKPT_coords(2)*B_matrix_SC(2,:) + &
                SCKPT_coords(3)*B_matrix_SC(3,:)
        enddo
    close(unit=input_file_unit)

end subroutine get_SCKPTS_contained_in_wavecar


!==============================================================================
!> \ingroup changes_upon_new_interface
!> Gets a list of the SC-KPTs present in the wavefunction file.
!> Similar routines are needed to implement interfaces to other PW codes.
!==============================================================================
subroutine get_SCKPTS_QE(list_of_SCKPTS, args)
implicit none
type(vec3d), dimension(:), allocatable, intent(out) :: list_of_SCKPTS
type(comm_line_args), intent(in) :: args
integer :: input_file_unit, ios, ik, nkpts, j, alloc_stat
character(len=256) :: prefix, outdir, xml_data_file_folder, xml_data_file, &
                      length_units, length_units_kpts
real(kind=dp), dimension(:,:), allocatable :: cart_coords_all_kpts
real(kind=dp) :: alat

#   if defined (__QE_SUPPORT__)
    input_file_unit = available_io_unit()
    outdir = trim(adjustl(args%qe_outdir))
    prefix = trim(adjustl(args%qe_prefix)) 
    xml_data_file_folder = trim(outdir) // '/' // trim(prefix) // '.save/'
    xml_data_file = trim(xml_data_file_folder) // "data-file.xml"

    ! Initializing the qexml library and opening the xml data-file
    call qexml_init(input_file_unit, dir=xml_data_file_folder) 
    call qexml_openfile(xml_data_file, "read", ierr=ios)
    if(ios/=0)then
        write(*,'(A)')"ERROR (get_SCKPTS_QE): Problems opening XML data-file!"
    endif  

    call qexml_read_bz(num_k_points=nkpts, ierr=ios)
    deallocate(list_of_SCKPTS, cart_coords_all_kpts, stat=alloc_stat)
    allocate(list_of_SCKPTS(1:nkpts))
    allocate(cart_coords_all_kpts(1:3,1:nkpts))

    ! Reading lattice scaling factor
    call qexml_read_cell(alat=alat, a_units=length_units, ierr=ios)
    if(ios/=0)then
        write(*,'(A)')"ERROR (get_SCKPTS_QE): Could not read lattice param!"
    endif
    alat = to_angstrom(alat, length_units)

    call qexml_read_bz(&
             xk=cart_coords_all_kpts, k_units=length_units_kpts, ierr=ios &
         )
    do ik=1, nkpts
        do j=1,3
            list_of_SCKPTS(ik)%coord(j) = (twopi/alat) * &
                                          cart_coords_all_kpts(j,ik)
        enddo
    enddo

    ! Done. Closing file now.
    call qexml_closefile("read", ierr=ios)
    if(ios/=0)then
        write(*,'(A)')"WARNING (get_SCKPTS_QE): Error closing xml data-file."
    endif
#   else
        write(*,'(A)')'ERROR (get_SCKPTS_QE): &
                      QE interface not compiled!'
        write(*,'(A)')'Cannot continue. Stopping now.'
        stop
#   endif

end subroutine get_SCKPTS_QE


!==============================================================================
!> \ingroup changes_upon_new_interface
!> Gets a list of the SC-KPTs present in the wavefunction file.
!> Similar routines are needed to implement interfaces to other PW codes.
!==============================================================================
subroutine get_SCKPTS_ABINIT(list_of_SCKPTS, args)
implicit none
type(vec3d), dimension(:), allocatable, intent(out) :: list_of_SCKPTS
type(comm_line_args), intent(in) :: args
integer :: io_unit, ikpt, nkpts, alloc_stat
type(vec3d), dimension(:), allocatable :: skpts_rec
real(kind=dp) :: Vcell
real(kind=dp), dimension(1:3) :: a1, a2, a3, b1, b2, b3
real(kind=dp), dimension(1:3, 1:3) :: A_matrix
logical :: wf_file_exists

    inquire(file=args%wf_file, exist=wf_file_exists)
    if(.not. wf_file_exists)then
        write(*,'(3A)')&
            "ERROR: Wavefunction file '", trim(adjustl(args%wf_file)), &
            "' not found."
        write(*,'(A)')'       Cannot continue. Stopping now.'
        stop
    endif

    io_unit = available_io_unit() 
    open(unit=io_unit,file=args%wf_file,form='unformatted',recl=1,status='old')
        call read_abinit_wfk_header(&
                 unit=io_unit, A_matrix=A_matrix, kpts=skpts_rec &
             )
    close(io_unit)
    ! Latt. vecs.
    A_matrix = bohr * A_matrix
    a1 = A_matrix(1,:)
    a2 = A_matrix(2,:)
    a3 = A_matrix(3,:)
    Vcell= dabs(dot_product(a1, cross(a2,a3)))
    ! Vectors of the reciprocal lattice.
    b1 = twopi*cross(a2,a3)/Vcell
    b2 = twopi*cross(a3,a1)/Vcell
    b3 = twopi*cross(a1,a2)/Vcell

    nkpts = size(skpts_rec)
    deallocate(list_of_SCKPTS, stat=alloc_stat)
    allocate(list_of_SCKPTS(1:nkpts))
    do ikpt=1, nkpts
        list_of_SCKPTS(ikpt)%coord(:) = skpts_rec(ikpt)%coord(1) * b1 + &
                                        skpts_rec(ikpt)%coord(2) * b2 + &
                                        skpts_rec(ikpt)%coord(3) * b3
    enddo

end subroutine get_SCKPTS_ABINIT


!==============================================================================
!> \ingroup changes_upon_new_interface
!> Gets a list of the SC-KPTs present in the wavefunction file.
!> Similar routines are needed to implement interfaces to other PW codes.
!==============================================================================
subroutine get_SCKPTS_CASTEP(list_of_SCKPTS, args)
implicit none
type(vec3d), dimension(:), allocatable, intent(out) :: list_of_SCKPTS
type(comm_line_args), intent(in) :: args
logical :: wf_file_exists

#   if defined (__CASTEP_SUPPORT__)
        inquire(file=args%wf_file, exist=wf_file_exists)
        if(.not. wf_file_exists)then
            write(*,'(3A)')&
                "ERROR: Wavefunction file '", trim(adjustl(args%wf_file)), &
                "' not found."
            write(*,'(A)')'       Cannot continue. Stopping now.'
            stop
        endif
        call get_list_of_kpts_in_orbitals_file(list_of_SCKPTS, args%wf_file)
#   else
        write(*,'(A)')'ERROR (get_SCKPTS_CASTEP): &
                      CASTEP interface not compiled!'
        write(*,'(A)')'Cannot continue. Stopping now.'
        stop
#   endif

end subroutine get_SCKPTS_CASTEP

!==============================================================================
!> \ingroup changes_upon_new_interface
!> Gets a list of all SC-kpoints contained in the wavefunction file(s).
!> The k-points should be in cartesian coordinates.
!==============================================================================
subroutine get_list_of_SCKPTS(list_of_SCKPTS, args, crystal_SC)
implicit none
type(vec3d), dimension(:), allocatable, intent(out) :: list_of_SCKPTS
type(comm_line_args), intent(in) :: args
type(crystal_3D), intent(in) :: crystal_SC

    select case(trim(adjustl(args%pw_code)))
        case('qe')
            call get_SCKPTS_QE(list_of_SCKPTS, args)
        case('abinit')
            call get_SCKPTS_ABINIT(list_of_SCKPTS, args)
        case('castep')
            call get_SCKPTS_CASTEP(list_of_SCKPTS, args)
        case default
            call get_SCKPTS_contained_in_wavecar(&
                     list_of_SCKPTS, args, crystal_SC &
                 )
    end select

end subroutine get_list_of_SCKPTS


subroutine write_band_struc(out_file,pckpts_to_be_checked,energy_grid, &
                            delta_N,EF,zero_of_kpts_scale,add_elapsed_time_to)
implicit none
character(len=*), intent(in) :: out_file
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
real(kind=dp), dimension(:), intent(in) :: energy_grid
type(UnfoldedQuantitiesForOutput), intent(in) :: delta_N
real(kind=dp), dimension(1:3) :: pckpt, first_pckpt_dir
real(kind=dp), intent(in), optional :: EF, zero_of_kpts_scale
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer :: nener, iener, idir, ndirs, nkpts, ikpt
real(kind=dp) :: e_fermi, origin_of_kpts_line, coord_first_k_in_dir, coord_k, &
                 stime, ftime
logical :: write_spin_info

    stime = time_now()

    e_fermi = 0.0_dp
    if(present(EF))then
         e_fermi = EF
    endif
    origin_of_kpts_line = 0.0_dp
    if(present(zero_of_kpts_scale))then
        origin_of_kpts_line = zero_of_kpts_scale
    endif


    write_spin_info = allocated(delta_N%pcbz_dir(1)%pckpt(1)%spin_proj_perp)
    nener = size(energy_grid)
    ndirs = size(pckpts_to_be_checked%selec_pcbz_dir(:))
    coord_first_k_in_dir = origin_of_kpts_line
    open(unit=12,file=out_file)
        write(12,'(A)')file_header_BandUP()
        if(write_spin_info)then
            write(12, '(A)')&
                '# NB.: Unfolding of Pauli matrices eigenvs is a test feature.'
            write(12, '(A, 2(ES10.3,", "), ES10.3, A)')&
                '# The quantization axis assumed was saxis = [',args%saxis,']'
            write(12, '(A, 2(ES10.3,", "), ES10.3, A)')&
                '# The projection axis used for "#Spin _|_ k" is also _|_ to &
                [',args%normal_to_proj_plane,'].'
            write(12,'(3(A,1X), 2X, 5(A,1X))')&
                "#KptCoord", "#E-E_Fermi","#delta_N", &
                "#sigma_x   ", "#sigma_y   ", "#sigma_z   ", &
                "#Spin _|_ k", "#Spin // k "
        else
            write(12,'(3(A,1X))')"#KptCoord", "#E-E_Fermi","#delta_N"
        endif
        do idir=1,ndirs
            first_pckpt_dir(:) = pckpts_to_be_checked%selec_pcbz_dir(idir)%&
                                     needed_dir(1)%pckpt(1)%coords(:)
            nkpts = size(&
                        pckpts_to_be_checked%selec_pcbz_dir(idir)%&
                            needed_dir(1)%pckpt(:) &
                    )
            do ikpt=1,nkpts
                pckpt(:) = pckpts_to_be_checked%selec_pcbz_dir(idir)%&
                               needed_dir(1)%pckpt(ikpt)%coords(:)
                do iener=1,nener
                    coord_k = coord_first_k_in_dir + &
                              norm(pckpt(:) - first_pckpt_dir(:))
                    if(write_spin_info)then
                        write(12,'(2(f8.4,2X),6(ES10.3, 2X))') &
                            coord_k, energy_grid(iener) - e_fermi, &
                            delta_N%pcbz_dir(idir)%pckpt(ikpt)%&
                                dN(iener), &
                            delta_N%pcbz_dir(idir)%pckpt(ikpt)%&
                                sigma(iener,:), &
                            delta_N%pcbz_dir(idir)%pckpt(ikpt)%&
                                spin_proj_perp(iener), &
                            delta_N%pcbz_dir(idir)%pckpt(ikpt)%&
                                spin_proj_para(iener)
                    else
                        write(12,'(2(f8.4,2X),ES10.3)')&
                            coord_k, energy_grid(iener) - e_fermi, &
                            delta_N%pcbz_dir(idir)%pckpt(ikpt)%dN(iener)
                    endif
                enddo
            enddo
            coord_first_k_in_dir = coord_k
        enddo
    close(12)

    ftime = time_now()
    if(present(add_elapsed_time_to))then
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine write_band_struc



subroutine read_pckpts_selected_by_user(k_starts,k_ends,ndirs,n_kpts_dirs, &
                                        input_file,b_matrix_pc, & 
                                        zero_of_kpts_scale,verbose)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
real(kind=dp), dimension(:,:), allocatable, intent(out) :: k_starts, k_ends
integer, intent(out) :: ndirs
integer, dimension(:), allocatable, intent(out) :: n_kpts_dirs
character(len=*), intent(in) :: input_file
real(kind=dp), dimension(1:3,1:3), intent(in) :: b_matrix_pc
real(kind=dp), intent(out), optional :: zero_of_kpts_scale
logical, intent(in), optional :: verbose
integer :: ios,idir,i
character(len=1) :: coords_type_1st_letter
real(kind=dp), dimension(:,:), allocatable :: read_k_start, read_k_end
real(kind=dp) :: a0, origin_of_kpts_line
character(len=str_len) :: coords_type_line, coords_type_flag, &
                          char_n_kpts_dirs, char_line_mode, aux_char, &
                          char_kstart, char_kend
logical :: opt_for_auto_pkpt_search, print_stuff, a0_informed_in_new_format, &
           a0_informed_in_old_format, give_tip_a0_for_reciprocal, &
           warn_old_format_cartesian

print_stuff = .TRUE.
if(present(verbose)) print_stuff = verbose
!! Getting the number of directions
ndirs = 0
open(unit=03,file=input_file)
    read(03,"(A)") aux_char
    read(03,"(A)") aux_char
    read(03,"(A)") aux_char
    read(03,"(A)") aux_char
    ios=0
    aux_char = ''
    do while (ios==0)
        read(03,"(A)",iostat=ios)char_kstart
        if(ios == 0 .and. trim(adjustl(char_kstart)) /= '')then
            read(03,"(A)",iostat=ios)char_kend
        endif
        if(ios == 0 .and. trim(adjustl(char_kstart)) /= ''  .and. &
           trim(adjustl(char_kend)) /= '')then        
            read(03,"(A)",iostat=ios)aux_char
            if((ios == 0 .and. trim(adjustl(aux_char)) == '') .or. ios < 0)then
                ndirs = ndirs + 1
            else
                if(ios==0) ios = 58
            endif
        endif
        if(ios > 0)then
            ios = 58
            write(*,'(A)')&
                'ERROR reading input pc-kpts file. Please check the format.'
            write(*,'(A)')'Stopping now.'
            stop
        endif
    enddo
close(03)
!! Actually reading the file now
allocate(read_k_start(1:ndirs, 1:3), read_k_end(1:ndirs, 1:3))
open(unit=03,file=input_file)
    read(03,*) aux_char
    ! Old style of passing the scaling factor "a0", needed when the k-points 
    ! are given in cartesian coordinates
    ! a0 is normally not present in the k-points file, but I chose to require
    ! it if the k-points are given in cartesian coordinates. 
    ! By doing so, one doesn't need a separate file to specify a0. 
    read(aux_char,*,iostat=ios) a0
    a0_informed_in_old_format = (ios==0)
    read(03,'(A)')char_n_kpts_dirs
    read(03,'(A)')char_line_mode
    read(03,'(A)')coords_type_line
    coords_type_line = trim(adjustl(coords_type_line))
    coords_type_1st_letter = coords_type_line(1:1)
    ! Passing a0 after "cartesian", recommended format now.
    read(coords_type_line,*,iostat=ios)coords_type_flag,a0
    a0_informed_in_new_format = (ios==0)

    do idir=1,ndirs
        read(03,*,iostat=ios)(read_k_start(idir,i), i=1,3)
        if(ios == 0) read(03,*,iostat=ios)(read_k_end(idir,i), i=1,3)
        if(ios/=0)then
            write(*,'(A)')'ERROR reading input pc-kpts file.'
            write(*,'(A)')&
                '      * Please check if there are numbers missing or if &
                there are extra non-number characters.'
            write(*,'(A)')'Stopping now.'
            stop
        endif
        read(03,*,iostat=ios)
    enddo
close(03)
!! Reading next the position where the user wants the kpts string to begin at.
!! This might be useful if you make separate calculations for each pcbz
!! direction and want to combine them in a single file afterwards
char_line_mode = trim(adjustl(char_line_mode))
read(char_line_mode,*,iostat=ios)aux_char, origin_of_kpts_line
if(ios /= 0)then
    origin_of_kpts_line = 0.0_dp
endif
if(present(zero_of_kpts_scale))then
    zero_of_kpts_scale = origin_of_kpts_line
endif

!! Making sure the k-points are returned in cartesian coordinates
allocate(k_starts(1:ndirs, 1:3), k_ends(1:ndirs, 1:3))
give_tip_a0_for_reciprocal = .FALSE.
warn_old_format_cartesian = .FALSE.
do idir=1,ndirs
    k_starts(idir,:) = 0.0_dp
    k_ends(idir,:) = 0.0_dp
    if(upper_case(coords_type_1st_letter) == 'C')then
        warn_old_format_cartesian = a0_informed_in_old_format
        if(.not. a0_informed_in_old_format .and. &
           .not. a0_informed_in_new_format)then
            write(*,'(A)')&
                'ERROR: You have selected cartesian coordinates in your &
                        input k-points file, but you have not passed a &
                        scaling parameter "a0".'
            write(*,'(A)')&
                '       The actuall coordiates of the k-points are given by: &
                ki[actual] = two_pi*ki[passed in file]/a0.'
            write(*,'(A)')&
                '       Please write the value of a0 after your tag "' // &
                trim(adjustl(coords_type_flag)) // '", and run BandUP again.'
            write(*,'(A)')'Stopping now.'
            stop
        endif
        k_starts(idir,:) = twopi*read_k_start(idir,:)/a0
        k_ends(idir,:) = twopi*read_k_end(idir,:)/a0
    else
        if((upper_case(coords_type_1st_letter) /= 'R') .and. &
           (idir==1) .and. print_stuff)then
            write(*,*)''
            write(*,'(A)')&
                'WARNING: Assuming that the pc-kpts have been informed in &
                fractional (reciprocal) coordinates.'
            write(*,*)''
        endif
        if(a0_informed_in_old_format .or. a0_informed_in_new_format)then
            give_tip_a0_for_reciprocal = .TRUE.
        endif
        do i=1,3
            k_starts(idir,:) = k_starts(idir,:) + &
                               read_k_start(idir,i) * b_matrix_pc(i,:)
            k_ends(idir,:) = k_ends(idir,:) + &
                             read_k_end(idir,i) * b_matrix_pc(i,:)
        enddo
    endif
enddo

!! Defining the number of k-points requested along each direction
allocate(n_kpts_dirs(1:ndirs))
char_n_kpts_dirs = trim(adjustl(char_n_kpts_dirs))
read(char_n_kpts_dirs,*,iostat=ios) n_kpts_dirs(1)
if(ios == 0)then
    opt_for_auto_pkpt_search = .FALSE.
    read(char_n_kpts_dirs,*,iostat=ios)(n_kpts_dirs(idir), idir=1,ndirs)
    if(ios /= 0)then
        n_kpts_dirs(:) = n_kpts_dirs(1)
    endif
else
    !! If BandUP cannot find a (sequence of) number(s) on "char_n_kpts_dirs", 
    !! then it assumes you want to perform an automatic scan for folding pckpts
    !! This will probably not be used if you use the BandUP's 
    !! pre-processing tool to get exactly the pckpts you want
    opt_for_auto_pkpt_search = .TRUE.
endif
!! If you want to parse only 1 pckpt, then use it as both the start and the end
!! of the kpts line.
do idir=1,ndirs
    if(norm(k_starts(idir,:)-k_ends(idir,:))<epsilon(1.0_dp))then
        n_kpts_dirs(idir) = 1
    endif
enddo
!! Setting now the number of pckpts along each direction,
!! according to whether the user has actually defined them or not
if(opt_for_auto_pkpt_search)then
    if(print_stuff)then
        write(*,'(A)')'Checking automatically for folding pc-kpts along &
                       the directions specified in the input file.'
        write(*,'(A,E10.3,A)')'The automatic scan for pc-kpts will be &
                               performed in intervals of ',min_dk,' A^-1.'
    endif
    do idir=1, ndirs
        n_kpts_dirs(idir) = &
            ceiling(&
                1.0_dp + norm(k_ends(idir,:) - k_starts(idir,:))/dabs(min_dk) &
            )
    enddo
endif
if(print_stuff)then
    do idir=1, ndirs
        if(opt_for_auto_pkpt_search)then
            write(*,'(2(A,I0))')&
                '      # of pc-kpts requested along pcbz direction ', &
                idir, ': ', n_kpts_dirs(idir)
        else
            write(*,'(2(A,I0))')&
                '# of pc-kpts requested along pcbz direction ', &
                idir, ': ', n_kpts_dirs(idir)
        endif
    enddo
endif

if(warn_old_format_cartesian)then
    write(*,*)
    write(*,'(A)')&
        '>>> NOTICE: Passing a0 in the first line of the k-points file is &
        deprecated and no longer recommended.'
    if(a0_informed_in_new_format)then
        write(*,'(A)')&
            '            * You have also passed a0 after "' // &
            trim(adjustl(coords_type_flag)) // '".'
        write(*,'(A,f0.5,A)')&
            '            * Only this value will be used (a0 =',a0,').'
    endif
    write(*,'(A)')&
        '            * It will work, but we now recommend that you specify a0 &
        after the tag "cartesian" (separated by a space).'
endif

if(give_tip_a0_for_reciprocal)then
    write(*,*)
    write(*,'(A)')&
        ">>> Tip: You don't need to pass the scaling parameter a0 when &
        the k-points are informed in fractional (reciprocal) coordinates."
endif

write(*,*)''

end subroutine read_pckpts_selected_by_user


!==============================================================================
!> \ingroup changes_upon_new_interface
!> General interface to the routines that read wavefunction file(s) from the
!! various codes supported by BandUP. It checks for the existence of the files,
!! decides on whether read the coefficients or not, checks the consistence of
!! the ikpt and spin choices and (optionally) renormalizes the wavefunctions.
!==============================================================================
subroutine read_wavefunction(wf, args, ikpt, read_coeffs, &
                             elapsed_time, add_elapsed_time_to, &
                             iostat, verbose, stop_if_not_found)
implicit none
type(pw_wavefunction), intent(inout) :: wf
type(comm_line_args), intent(in) :: args
integer, intent(in) :: ikpt
logical, intent(in), optional :: read_coeffs, verbose, stop_if_not_found
real(kind=dp), intent(out), optional :: elapsed_time
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
integer, intent(out), optional :: iostat
integer :: ios, i_kpt, iband, i
real(kind=dp) :: stime, ftime, inner_prod
logical :: file_exists, spin_reset, read_coefficients, print_stuff, &
           stop_if_wf_not_found

    stime = time_now()

    print_stuff = .FALSE.
    if(present(verbose)) print_stuff = verbose
    read_coefficients = .TRUE. ! Reading the coeffs by default
    if(present(read_coeffs)) read_coefficients = read_coeffs
    stop_if_wf_not_found = .FALSE.
    if(present(stop_if_not_found)) stop_if_wf_not_found = stop_if_not_found


    ios = 0
    ! Check if file exists
    select case(trim(adjustl(args%pw_code)))
        case('qe')
            continue ! Still to implement this test
        case default ! If VASP, ABINIT or CASTEP
            inquire(file=args%WF_file, exist=file_exists)
            if(.not. file_exists)then
                write(*,'(3A)')"ERROR (read_wavefunction): File ", &
                               trim(adjustl(args%WF_file)), " doesn't exist."
                ios = -1
                if(stop_if_wf_not_found)then
                    write(*,'(A)')'Stopping now.'
                    stop
                endif
                if(present(iostat)) iostat = ios
                return
            endif
    end select

    i_kpt = ikpt
    ! If the user doesn't specify a valid kpt number, then ikpt=1 will be read
    if(i_kpt < 1) i_kpt = 1

    if(read_coefficients .and. print_stuff)then
        write(*,"(A,I0,A)")&
            'Reading plane-wave coefficients for SC-Kpoint K(',i_kpt,')...'
    endif

    wf%i_spin = args%spin_channel
    spin_reset = (wf%i_spin < 1 .or. wf%i_spin > 2)
    if(spin_reset)then
        if(wf%i_spin < 1) wf%i_spin = 1
        if(wf%i_spin > 2) wf%i_spin = 2
        write(*,'(A,I1)')&
            'WARNING (read_wavefunction): Spin channel reset to ', wf%i_spin
    endif

    select case(trim(adjustl(args%pw_code)))
        case default ! Using VASP as default
            call read_wavecar(wf, file=args%WF_file, &
                              continue_if_npw_smaller_than_expected=&
                              args%continue_if_npw_smaller_than_expected, &
                              ikpt=i_kpt, &
                              read_coeffs=read_coefficients, iostat=ios)
        case('abinit')
            call read_abinit_wfk_file(&
                     wf, file=args%WF_file, ikpt=i_kpt, &
                     read_coeffs=read_coefficients, iostat=ios &
                 )
        case('qe')
#           if defined (__QE_SUPPORT__)
                call read_qe_evc_file(wf, args, i_kpt, read_coefficients, ios)
#           else
                write(*,'(A)')&
                    'ERROR (read_wavefunction): QE interface not compiled!'
                write(*,'(A)')'Cannot continue. Stopping now.'
                stop
#           endif
        case('castep')
#           if defined (__CASTEP_SUPPORT__)
                call read_castep_orbitals_file(wf, args%WF_file, i_kpt, &
                                               read_coefficients, ios)
#           else
                write(*,'(A)')&
                    'ERROR (read_wavefunction): CASTEP interface not compiled!'
                write(*,'(A)')'Cannot continue. Stopping now.'
                stop
#           endif
    end select
    if(read_coefficients .and. renormalize_wf)then
        !$omp parallel do default(none) schedule(guided) &
        !$omp private(iband, inner_prod) &
        !$omp shared(wf)
        do iband=1, wf%n_bands
            inner_prod = sum((/(dot_product(wf%pw_coeffs(i,:,iband), &
                                            wf%pw_coeffs(i,:,iband)), &
                                            i=1, wf%n_spinor)/))
            if(inner_prod > epsilon(1.0_dp))then
                wf%pw_coeffs(:,:,iband) = (1.0_dp/sqrt(abs(inner_prod))) * &
                                           wf%pw_coeffs(:,:,iband)
            endif
        enddo
    endif

    ftime = time_now()
    if(present(elapsed_time)) elapsed_time = ftime - stime
    if(present(add_elapsed_time_to))then
        add_elapsed_time_to = add_elapsed_time_to  + (ftime - stime)
    endif
    if(present(iostat)) iostat = ios

    if(read_coefficients .and. print_stuff)then
        write(*,'(A,f0.1,A)')'    * Done in ', (ftime - stime),'s.'
    endif

    return

end subroutine read_wavefunction


subroutine print_symm_analysis_for_selected_pcbz_dirs(&
               all_dirs_used_for_EBS_along_pcbz_dir &
           ) 
implicit none
type(irr_bz_directions), dimension(:), intent(in) :: &
    all_dirs_used_for_EBS_along_pcbz_dir
integer :: ndirs, idir, n_needed_dirs, neqv_dirs_pcbz, neqv_dirs_SCBZ, &
           ncompl_dirs,n_irr_compl_dirs

if(args%no_symm_avg)then
    write(*,"(A)")'Running with the flag "-no_symm_avg".'
    write(*,"(A)")&
        '    * No avg over complementary pcbz directions will be performed.'
else
    ndirs = size(all_dirs_used_for_EBS_along_pcbz_dir(:))
    write(*,"(A)")"Symmetry analysis for the selected pcbz directions:"
    do idir=1,ndirs
        neqv_dirs_pcbz = all_dirs_used_for_EBS_along_pcbz_dir(idir)%neqv
        neqv_dirs_SCBZ = all_dirs_used_for_EBS_along_pcbz_dir(idir)%neqv_SCBZ
        write(*,"(A,I0,A)")"    >>> Direction #",idir,":"
        write(*,"(A,I0,A)")&
            "        * Found ", neqv_dirs_pcbz, &
            " equivalent directions w.r.t. symmetry operations of the pc"
        write(*,"(A,I0,A)")&
            "        * Found ", neqv_dirs_SCBZ, &
            " equivalent directions w.r.t. symmetry operations of the SC"

        n_needed_dirs = size(&
                            all_dirs_used_for_EBS_along_pcbz_dir(idir)%&
                                irr_dir(:) &
                        )
        if(n_needed_dirs > 1)then
            ncompl_dirs = all_dirs_used_for_EBS_along_pcbz_dir(idir)%&
                              ncompl_dirs 
            n_irr_compl_dirs = all_dirs_used_for_EBS_along_pcbz_dir(idir)%&
                                   n_irr_compl_dirs
            write(*,"(A,I0,A)")&
                "        * ",ncompl_dirs, &
                " complementary pcbz directions will be considered in order to"
            write(*,"(A)")"          get a symmetry-averaged EBS."
            if(ncompl_dirs /= n_irr_compl_dirs)then
                write(*,"(A,I0,A)")&
                    "        * The number of irreducible complementary &
                    directions is ",n_irr_compl_dirs,"." 
            endif
        else
            write(*,"(A)")&
                "        * No complementary pcbz directions are needed."
        endif
    enddo
endif
write(*,*)

end subroutine print_symm_analysis_for_selected_pcbz_dirs


subroutine print_simm_ops(isym,trans,rot)
implicit none
integer, intent(in) :: isym
integer, intent(in), dimension(1:3) :: trans
integer, intent(in), dimension(1:3,1:3) :: rot

    write(*,'(8X,A,I0,A)')'Symmetry op. #',isym,':'
    write(*,'(12X,A,4(I0,A))')&
        'T[',isym,'] = (',trans(1),', ',trans(2),', ',trans(3),')'
    write(*,'(12X,A,4(I0,A))')&
        'R[',isym,'](1,:) = (',rot(1,1),', ',rot(1,2),', ',rot(1,3),')'
    write(*,'(12X,A,4(I0,A))')&
        'R[',isym,'](2,:) = (',rot(2,1),', ',rot(2,2),', ',rot(2,3),')'
    write(*,'(12X,A,4(I0,A))')&
        'R[',isym,'](3,:) = (',rot(3,1),', ',rot(3,2),', ',rot(3,3),')'
    write(*,*)

end subroutine print_simm_ops


subroutine print_message_pckpt_cannot_be_parsed(&
               stop_when_a_pckpt_cannot_be_parsed &
           )
implicit none
logical, intent(in) :: stop_when_a_pckpt_cannot_be_parsed

    if(stop_when_a_pckpt_cannot_be_parsed)then
        write(*,'(A)')&
            "    ERROR: Could not calculate the spectral weight for the pair &
            (k, K) being parsed now."
        write(*,'(A)')"           The cause is probably that either:"
        write(*,'(A)')&
            "               * The folding vector G = k - K has been estimated &
            incorrectly by the code (symmetry issue), or"
        write(*,'(A)')"               * K doesn't unfold onto k after all."
        write(*,'(A)')&
            "           Have you used the pre-procesing tool to &
            find the SC-Kpts you needed?"
        write(*,'(A)')"               * If not, please try that."
        write(*,'(A)')"               * If yes, please contact me (Paulo)."
        write(*,'(A)')"           Stopping now."
    else
        write(*,'(A)')&
            "    WARNING: Not enough coefficients in the wavefunction file &
            to unfold this pc wave-vector."
        write(*,'(A)')&
            "             Be careful with your results: They might be &
            incomplete or even wrong."
    endif

end subroutine print_message_pckpt_cannot_be_parsed


subroutine print_message_success_determining_GUR(&
               GUR, stop_if_GUR_fails, is_main_code &
           )
implicit none
!! Geometric Unfolding Relations
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR
logical, intent(in) :: stop_if_GUR_fails, is_main_code

    if(GUR%n_pckpts /= GUR%n_folding_pckpts)then ! if GUR failed
        if(stop_if_GUR_fails)then
            write(*,'(A)')&
                'ERROR: Could not determine all gemetric unfolding relations &
                between the chosen SC-Kpts and pc-kpts.'
        else
            write(*,'(A)')&
                'WARNING: Could not determine all gem. unfolding relations &
                between the chosen SC-Kpts and pc-kpts.'
            write(*,'(A)')'         The code may fail!'
        endif
        write(*,'(A,I0,A)')&
            '       * ',GUR%n_pckpts,' pc-kpts have been checked.'
        write(*,'(A,I0,A)')&
            '       * ',GUR%n_folding_pckpts, &
            ' pc-kpts satisfied the geometric unfolding relations.'
        if(is_main_code)then
            write(*,'(A)')&
                '       Have you used the pre-procesing tool &
                to find the SC-Kpts you needed?'
            write(*,'(A)')'       * If not, please try that.'
            write(*,'(A)')'       * If yes, please contact me (Paulo).'
        endif
        if(stop_if_GUR_fails)then
            write(*,'(A)')'       Stopping now.'
        endif
    else
        write(*,'(A)')&
            'The geometric unfolding relations have successfully been &
            determined. Good!'
        write(*,'(A,I0,A)')&
            '    * ',GUR%n_folding_pckpts, &
            ' pc-kpts satisfied the geometric unfolding relations.'
    endif

end subroutine print_message_success_determining_GUR


subroutine say_goodbye_and_save_results(&
               delta_N_only_selected_dirs, &
               delta_N_symm_avrgd_for_EBS, &
               GUR, &
               pckpts_to_be_checked, energy_grid, &
               e_fermi, zero_of_kpts_scale, &
               n_input_pc_kpts,n_folding_pckpts,&
               n_folding_pckpts_parsed, times &
           )
implicit none
type(UnfoldedQuantitiesForOutput),  intent(in) :: delta_N_only_selected_dirs, &
                                                  delta_N_symm_avrgd_for_EBS
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
real(kind=dp), dimension(:), intent(in) :: energy_grid
real(kind=dp), intent(in) :: e_fermi, zero_of_kpts_scale
integer, intent(in) :: n_input_pc_kpts,n_folding_pckpts,n_folding_pckpts_parsed
type(timekeeping), intent(inout) :: times


    write(*,*)
    write(*,'(A)')'Band unfolding process finished. &
                   The output files will now be written.'
    write(*,'(2(A,I0),A)')'A total of ',n_folding_pckpts, &
                          ' pc-kpts (out of the ',n_input_pc_kpts, &
                          ' checked) satisfied the folding condition.'
    if(n_folding_pckpts_parsed > 0)then
        if(n_folding_pckpts_parsed /= n_folding_pckpts)then
            write(*,'(A,I0,A)')&
                'From these points, ', n_folding_pckpts_parsed, &
                ' could be used.'
        endif

        if(args%no_symm_avg)then
            write(*,"(A)")'>>> No symmetry-averaged EBS has been &
                               calculated ("-no_symm_avg" flag used).'
        else
            write(*,'(A)')&
                '>>> Writing the symm-averaged unfolded N(k,E) for the EBS...'
            call write_band_struc(&
                     args%output_file_symm_averaged_EBS, &
                     pckpts_to_be_checked, energy_grid, &
                     delta_N_symm_avrgd_for_EBS, &
                     EF=e_fermi, zero_of_kpts_scale=zero_of_kpts_scale, &
                     add_elapsed_time_to=times%write_dN_files &
                 )
            write(*,'(A)')&
                '    Done. The symmetry-averaged unfolded N(k,E) were stored &
                in the file listed below:'
            write(*,'(2A)')&
                '          * ', &
                trim(adjustl(args%output_file_symm_averaged_EBS))
        endif

        !! Writing delta_N_only_selected_dirs
        write(*,'(A)')&
            '>>> Writing the not-symm-averaged unfolded N(k,E) for the EBS...'
        call write_band_struc(&
                 args%output_file_only_user_selec_direcs, &
                 pckpts_to_be_checked, energy_grid, &
                 delta_N_only_selected_dirs, &
                 EF=e_fermi, zero_of_kpts_scale=zero_of_kpts_scale, &
                 add_elapsed_time_to=times%write_dN_files &
             )
        write(*,'(A)')&
            '    Done. The unfolded N(k,E) calculated strictly along the &
            direction(s) you'
        write(*,'(A)')&
            '          requested were stored in the file listed below:'
        write(*,'(2A)')&
            '          * ', &
            trim(adjustl(args%output_file_only_user_selec_direcs))
         
        if(zero_of_kpts_scale > 1E-4_dp)then
            write(*,'(A,f8.4,A)')&
                'The zero of the k-points line has been set to ', &
                zero_of_kpts_scale, '.'
        endif
        
        if(args%write_unf_dens_op)then
            if(.not. args%no_symm_avg)then
                write(*,'(A)')&
                    '>>> Writing the symm-averaged unfolding-density &
                    operators...'
                call write_unf_dens_op(&
                         args%output_file_symm_averaged_unf_dens_op, &
                         pckpts_to_be_checked, energy_grid, &
                         GUR, delta_N_symm_avrgd_for_EBS, &
                         EF=e_fermi, &
                         add_elapsed_time_to=times%write_unf_dens_op_files &
                     )
                write(*,'(A)')&
                    '    Done. The symmetry-averaged unfolding-density &
                    operators were stored in the file'
                write(*,'(A)')'          listed below:'
                write(*,'(2A)')&
                    '          * ', &
                    trim(adjustl(args%output_file_symm_averaged_unf_dens_op))
            endif
            write(*,'(A)')&
                '>>> Writing the not-symm-averaged unfolding-density &
                operators...'
            call write_unf_dens_op(&
                     args%output_file_only_user_selec_direcs_unf_dens_op, &
                     pckpts_to_be_checked, energy_grid, &
                     GUR, delta_N_only_selected_dirs, &
                     EF=e_fermi, &
                     add_elapsed_time_to=times%write_unf_dens_op_files &
                 )
            write(*,'(A)')&
                '    Done. The unfolding-density operators calculated &
                strictly along the directions'
            write(*,'(A)')&
                '          you requested have been saved to the file &
                listed below:'
            write(*,'(2A)')&
                '          * ', &
                trim(&
                    adjustl(&
                        args%output_file_only_user_selec_direcs_unf_dens_op &
                    ) &
                )
        endif
        write(*,'(A)')'Done writing output files.'
    else
        write(*,'(A)')&
            'No pc-kpts parsed. Nothing to be written to output files.'
    endif
    write(*,*)
    write(*,'(A)')'BandUP has now finished running. Good bye.'

end subroutine say_goodbye_and_save_results


subroutine write_unf_dens_op(out_file, pckpts, energy_grid, GUR, delta_N, &
                             EF, add_elapsed_time_to)
implicit none
character(len=*), intent(in) :: out_file
type(selected_pcbz_directions), intent(in) :: pckpts
real(kind=dp), dimension(:), intent(in) :: energy_grid
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR
type(UnfoldedQuantitiesForOutput),  intent(in), target :: delta_N
real(kind=dp), intent(inout), optional :: add_elapsed_time_to
real(kind=dp), intent(in), optional :: EF
! Internal variables
integer :: n_SCKPTS, n_pc_direcs, unf_dens_file_unit, idir, ipc_kpt, &
           n_pc_kpts, m1, m2, i_SCKPT, i_rho, iener,  &
           ipc_kpt_general, linearized_band_index
logical :: folds
character(len=127) :: fmt_str, str_trace
real(kind=dp) :: stime, ftime, coord_first_k_in_dir, coord_k, trace, e_fermi
real(kind=dp), dimension(1:3) :: pckpt, actual_folding_pckpt, sckpt, &
                                 aux_coords, first_pckpt_dir
type(UnfoldDensityOpContainer), dimension(:), pointer :: rhos
real(kind=dp), dimension(:), pointer :: dN
real(kind=dp), parameter :: min_allowed_dN = 1E-03_dp

    stime = time_now()

    e_fermi = 0.0_dp
    if(present(EF))then
         e_fermi = EF
    endif

    n_SCKPTS = size(GUR%SCKPT)
    n_pc_direcs = size(GUR%SCKPT(1)%selec_pcbz_dir)
    unf_dens_file_unit = available_io_unit()
    open(unit=unf_dens_file_unit, file=out_file, action='write')
        write(unf_dens_file_unit, '(A)')&
            '############################################&
             ############################################'
        write(unf_dens_file_unit, '(A)')file_header_BandUP()
        write(unf_dens_file_unit, '(A)')&
            '############################################&
             ############################################'
        write(unf_dens_file_unit, '(A)')&
            '# This file contains all information needed to &
            unfold the expectation values of any'
        write(unf_dens_file_unit, '(A)')&
            '# operator defined in the (k,E) space. &
            In particular, it reports, for each point in'
        write(unf_dens_file_unit, '(A,ES10.3,A)')&
            '# the (k,E) grid satisfying N(k,E) >= ',min_allowed_dN, &
            ', the following information:'
        write(unf_dens_file_unit, '(A)')&
            '#     * Unfolding-density operators (UnfDensOp)'
        write(unf_dens_file_unit, '(A)')&
            '#     * Generalized spectral weight matrices &
            (GenSpecWeight)'
        write(unf_dens_file_unit, '(A)')&
            '#     * k, E and N(k,E) -- used for "regular" &
            band unfolding (unfolding eigenvalues of H)'
        write(unf_dens_file_unit, '(A)')&
            '# Both UnfDensOp and GenSpecWeight are &
            hermitian; only the upper-triangular matrix'
        write(unf_dens_file_unit, '(A)')&
            '# elements are present. They are represented &
            as (RealPart, ImagPart), and were'
        write(unf_dens_file_unit, '(A)')&
            '# evaluated between SC states with band &
            indices indicated by m1 and m2.'
        write(unf_dens_file_unit,'(A)')'#'
        write(unf_dens_file_unit, '(A)')&
            '# All quantities reported in this file are &
            introduced and discussed in detail in'
        write(unf_dens_file_unit, '(A,25X,A)')&
            '#','Phys. Rev. B 91, 041116(R) (2015)'
        write(unf_dens_file_unit, '(A,25X,A)')&
            '#','Phys. Rev. B 89, 041407(R) (2014)'
        write(unf_dens_file_unit,'(A)')&
            '# You should read and cite these papers, as &
            well as the appropriate references therein.'
        write(unf_dens_file_unit,'(A)')'#'
        write(unf_dens_file_unit, '(A)')&
            '# N.B.: If symmetry has been used for the &
            unfolding, then the folding/unfolding'
        write(unf_dens_file_unit, '(A)')&
            '# relations between the ScKpts and PcKpts &
            shown here might not seem obvious.'
        write(unf_dens_file_unit, '(A)')&
            '# For full details about the geometric &
            unfolding relations, please see the log'
        write(unf_dens_file_unit, '(A)')&
            '# produced by BandUP when you run &
            the code (normally printed to stdout).'
        write(unf_dens_file_unit, '(A)')'#'
        write(unf_dens_file_unit, '(A,1X,I0)')&
            '# SpinChannel =', args%spin_channel
        write(unf_dens_file_unit, '(A,1X,I0)')&
            '# nScBands =', delta_N%n_SC_bands
        write(unf_dens_file_unit, '(A,1X,f0.5,1X,A)')&
            '# emin =', minval(energy_grid)-e_fermi, 'eV'
        write(unf_dens_file_unit, '(A,1X,f0.5,1X,A)')&
            '# emax =', maxval(energy_grid)-e_fermi, 'eV'
        write(unf_dens_file_unit, '(A,1X,I0)')'# nEner =', size(energy_grid)
        write(unf_dens_file_unit, '(A,1X,f0.5,1X,A)')&
            '# OriginalEfermi =',e_fermi,'eV'
        write(unf_dens_file_unit, '(A)')&
            '# The energies were shifted so that the Fermi level is at 0 eV.'
        write(unf_dens_file_unit, '(A)')&
            '############################################&
             ############################################'

        ipc_kpt_general = 0
        coord_first_k_in_dir = 0.0_dp
        do idir=1, n_pc_direcs
            first_pckpt_dir(:) = pckpts%selec_pcbz_dir(idir)%needed_dir(1)%&
                                        pckpt(1)%coords(:)
            n_pc_kpts=size(pckpts%selec_pcbz_dir(idir)%needed_dir(1)%pckpt(:))
            do ipc_kpt=1,n_pc_kpts
                ipc_kpt_general = ipc_kpt_general + 1
                ! Finding the SC-KPT this PC-KPT folds into
                folds = .FALSE.
                do i_SCKPT=1, n_SCKPTS
                    if(.not. GUR%SCKPT_used_for_unfolding(i_SCKPT)) cycle
                    folds = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(idir)% &
                                needed_dir(1)%pckpt(ipc_kpt)%folds
                    if(folds) exit
                enddo
                if(.not. folds) cycle
                sckpt(:) = GUR%list_of_SCKPTS(i_SCKPT)%coord(:)
                pckpt(:) = pckpts%selec_pcbz_dir(idir)% &
                                  needed_dir(1)%pckpt(ipc_kpt)%coords(:)
                actual_folding_pckpt = GUR%SCKPT(i_SCKPT)%&
                                           selec_pcbz_dir(idir)%&
                                           needed_dir(1)% &
                                           pckpt(ipc_kpt)% &
                                           Scoords(:)
                coord_k = coord_first_k_in_dir + &
                          norm(pckpt(:) - first_pckpt_dir(:))
                write(unf_dens_file_unit,'(A)')
                write(unf_dens_file_unit, '(A,1X,I0)')&
                    '# PcKptNumber =',ipc_kpt_general
                fmt_str = '(A,3(2X,f0.8),2X,A)'
                write(unf_dens_file_unit, fmt_str)&
                    '#     PcKptCartesianCoords =', pckpt,'(A^{-1})'
                fmt_str = '(A,3(2X,f0.8),2X,A)'
                aux_coords = coords_cart_vec_in_new_basis(&
                                 cart_vec=pckpt, &
                                 new_basis=GUR%B_matrix_SC &
                             )
                write(unf_dens_file_unit, fmt_str)&
                    '#     PcKptFractionalCoords =', &
                    aux_coords, '(w.r.t. SCRL basis)'
                fmt_str = '(A,3(2X,f0.8),2X,A)'
                aux_coords = coords_cart_vec_in_new_basis(&
                                 cart_vec=pckpt, &
                                 new_basis=GUR%b_matrix_pc &
                             )
                write(unf_dens_file_unit, fmt_str)&
                    '#     PcKptFractionalCoords =', aux_coords, &
                    '(w.r.t. PCRL basis)'
                fmt_str = '(A,1X,f0.4,2X,A)'
                write(unf_dens_file_unit, fmt_str)&
                    '#     PcKptLinearCoordsInBandPlot =', coord_k, '(A^{-1})'
                write(unf_dens_file_unit, '(A,1X,I0)')&
                    '#     Folds into ScKptNumber =', i_SCKPT
                fmt_str = '(A,3(2X,f0.8),2X,A)'
                write(unf_dens_file_unit, fmt_str)&
                    '#         ScKptCartesianCoords =', sckpt, '(A^{-1})'
                fmt_str = '(A,3(2X,f0.8),2X,A)'
                aux_coords = coords_cart_vec_in_new_basis(&
                                 cart_vec=sckpt, &
                                 new_basis=GUR%B_matrix_SC &
                             )
                write(unf_dens_file_unit, fmt_str)&
                    '#         ScKptFractionalCoords =', aux_coords, &
                    '(w.r.t. SCRL basis)'
                write(unf_dens_file_unit, '(A)')'#'

                rhos => delta_N%pcbz_dir(idir)%pckpt(ipc_kpt)%rhos
                dN => delta_N%pcbz_dir(idir)%pckpt(ipc_kpt)%dN
                do i_rho=1, size(rhos)
                    if(.not. allocated(rhos(i_rho)%band_indices)) cycle
                    iener = rhos(i_rho)%iener_in_full_pc_egrid
                    if(dN(iener) < min_allowed_dN) cycle
                    trace = 0.0_dp
                    do m1=1,rhos(i_rho)%nbands
                        linearized_band_index = &
                            list_index(&
                                item=[m1,m1], &
                                list=rhos(i_rho)%band_indices &
                            )
                        if(linearized_band_index<1) cycle
                        trace = trace + &
                                real(rhos(i_rho)%rho(linearized_band_index), &
                                    kind=dp &
                                )
                    enddo
                    if(abs(trace-1.0_dp)<1E-1)then
                        ! The trace of the unfolding-density operator should
                        ! equal 1, as discussed in 
                        ! Phys. Rev. B 91, 041116(R) (2015)
                        write(str_trace, '(ES8.1)') trace
                    else
                        write(str_trace, '(ES8.1,1X,A)') trace, '(!=1)'
                    endif
                    fmt_str = &
                    '(A,5X,A,1X,I0,2X,A,1X,f0.4,1X,A,2X,A,1X,ES10.3,2X,A,1X,A)'
                    write(unf_dens_file_unit, fmt_str)&
                        '#','iEnerPC =', &
                        iener,'E =', &
                        energy_grid(iener), 'eV', &
                        'N(k,E) =', dN(iener), &
                        'Tr{UnfDensOp} =', trim(str_trace)
                    fmt_str='(A,8X,A,6X,A,8X,A,9X,A)'
                    write(unf_dens_file_unit, fmt_str)&
                        '#','m1','m2', &
                        'GenSpecWeight[m1,m2]', &
                        'UnfDensOp[m1,m2]'
                    fmt_str = &
                        '(3X,2(3X,I5),3X,2(3X,"(",ES10.3,",",1X,ES10.3,")"))'
                    do m1=1,rhos(i_rho)%nbands
                        do m2=m1,rhos(i_rho)%nbands
                            linearized_band_index = &
                                list_index(item=[m1,m2], &
                                           list=rhos(i_rho)%band_indices)
                            if(linearized_band_index<1) cycle
                            if(abs(rhos(i_rho)%rho(linearized_band_index)) < &
                               1E-4) cycle
                            write(unf_dens_file_unit, fmt_str) m1, m2, &
                                rhos(i_rho)%rho(linearized_band_index) * &
                                    dN(iener), &
                                rhos(i_rho)%rho(linearized_band_index)
                        enddo
                    enddo
                enddo
                nullify(rhos, dN)
            enddo
            coord_first_k_in_dir = coord_k
        enddo

    close(unf_dens_file_unit)

    ftime = time_now()
    if(present(add_elapsed_time_to))then
        add_elapsed_time_to = add_elapsed_time_to + (ftime - stime)
    endif

end subroutine write_unf_dens_op


subroutine print_time(t_in_sec, message, ref_time)
implicit none
real(kind=dp), intent(in) :: t_in_sec
character(len=*), intent(in) :: message
real(kind=dp), intent(in), optional :: ref_time
real(kind=dp) :: time_percentual

    time_percentual = 100.0_dp
    if(present(ref_time))then
        if(ref_time > epsilon(1.0_dp))then
            time_percentual = 100.0_dp*t_in_sec/ref_time
        endif
    endif

    if(abs(time_percentual - 100.0_dp) <= epsilon(1.0_dp))then
        if(t_in_sec < 60.0_dp)then
            write(*,'(A,f7.2,A)') message, t_in_sec, 's.'
        else
            write(*,'(A,f7.2,3A)') message, t_in_sec, 's (', &
                                   trim(adjustl(formatted_time(t_in_sec))),').'
        endif
    else
        if(t_in_sec >= 0.1_dp .or. (time_percentual >= 1.0_dp))then
            if(t_in_sec < 60.0_dp)then
                write(*,'(A,f7.2,A,f0.1,A)')&
                    message, t_in_sec, 's (', time_percentual, '%).'
            else
                    write(*,'(A,f7.2,A,f0.1,3A)') &
                        message, t_in_sec, 's (', &
                        time_percentual, '%, ', &
                        trim(adjustl(formatted_time(t_in_sec))), ').'
            endif
        endif
    endif

end subroutine print_time


subroutine print_final_times(times)
implicit none
type(timekeeping), intent(inout) :: times
real(kind=dp) :: elapsed_time

    times%end = time_now()
    elapsed_time = times%end - times%start
    call print_time(elapsed_time, &
         'Total elapsed time:                                               ',&
         elapsed_time)
    call print_time(times%read_wf, &
         'Time spent reading the wavefunctions file:                        ',&
         elapsed_time)
    call print_time(times%calc_spec_weights, &
         'Time spent calculating spectral weights:                          ',&
         elapsed_time)
    call print_time(times%calc_SF, &
         'Time spent calculating spectral functions:                        ',&
         elapsed_time)
    call print_time(times%calc_dN, &
         'Time spent calculating the unfolded N(k,E):                       ',&
         elapsed_time)
    call print_time(times%calc_rho, &
         'Time spent calculating unfolding-density operators:               ',&
         elapsed_time)
    call print_time(times%calc_pauli_vec, &
         'Time spent calculating SC matrix elements of Pauli spin matrices: ',&
         elapsed_time)
    call print_time(times%calc_pauli_vec_projs, &
         'Time spent calculating projections of Pauli spin matrices:        ',&
         elapsed_time)
    call print_time(times%write_dN_files, &
         'Time spent writing N(k,E) to output file(s):                      ',&
         elapsed_time)
    call print_time(times%write_unf_dens_op_files, &
         'Time spent writing unfolding-density operators to output file(s): ',&
         elapsed_time)


end subroutine print_final_times



end module io_routines
