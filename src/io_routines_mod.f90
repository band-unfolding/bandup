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
module io_routines
use strings
use general_io
use read_wavecar
use read_vasp_files
use write_vasp_files
use math
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: print_welcome_messages,print_message_commens_test, &
          read_energy_info_for_band_search, print_message_success_determining_GUR, &
          print_geom_unfolding_relations, print_last_messages_before_unfolding, &
          get_SCKPTS_contained_in_wavecar, read_pckpts_selected_by_user, &
          print_symm_analysis_for_selected_pcbz_dirs, say_goodbye_and_save_results, & 
          print_final_times, write_band_struc, print_message_pckpt_cannot_be_parsed, &
          write_attempted_pc_assoc_with_input_unit_cell_and_SC


CONTAINS 


subroutine print_welcome_messages(package_version)
implicit none
character(len=*), intent(in), optional :: package_version
write(*,'(A,/,A)')   '=====================================================================================', &
                     '             BandUP: Band Unfolding code for Plane-wave based calculations           '
if(present(package_version))then
    write(*,'(A)')   "                            V. "//trim(adjustl(package_version))
endif
write(*,'(8(A,/),A)')'=====================================================================================', &
                     'Copyright (C) 2013, 2014 Paulo V. C. Medeiros                                        ', &
                     '                         paume@ifm.liu.se                                            ', &
                     '                         Computational Physics Division                              ', &
                     '                         Department of Physics, Chemistry and Biology - IFM          ', &
                     '                         Linköping University                                        ', &
                     '                         Sweden                                                      ', & 
                     'Please visit www.ifm.liu.se/theomod/compphys/band-unfolding                          ', &
                     '====================================================================================='
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
    message_header = '                                     ERROR!!                                         '
    message_footer = ' Stopping now.'
else
    message_header = '                               >>> WARNING!! <<<                                     '
    message_footer = ' >>> The results might not be what you expect. Continue only if you are really sure. '
    message_footer2= "     My advice: Always use commensurate SC and PC, and never continue if they're not."
endif

write(*,*)''
if(.not. commensurate)then
    write(*,*)''
    write(*,'(10(A,/))')'=====================================================================================', &
                                                          message_header                                       , &
                        '=====================================================================================', &
                        '  The SC and the reference PC that you have chosen are not commensurate!             ', &
                        '  The choice of the reference PC vectors is very important and can change a lot the  ', &
                        '  results of the unfolding. It is particularly important to always verify that the SC', & 
                        '  and the reference PC are commensurate. If this condition is not fulfilled, then the', &
                        '  calculated spectral weigths will most likely be very small, which might therefore  ', & 
                        '  cause the values of the unfolded delta_Ns to be also quite small.                  ', &
                        '====================================================================================='
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
        if(n_digits > max_n_digits_before_dec_point) max_n_digits_before_dec_point = n_digits
    enddo
    write(str_n_char_float_format,*) max_n_digits_before_dec_point + n_decimals + 1    
    float_format(j) = 'f' // trim(adjustl(str_n_char_float_format)) // '.' // trim(adjustl(str_n_decimals))
enddo
format_string = '(2(A,/),3(A,' // adjustl(trim(float_format(1))) // &
                         ',A,' // adjustl(trim(float_format(2))) // &
                         ',A,' // adjustl(trim(float_format(3))) // &
                         ',A,/))'
write(*,format_string)' The following relation holds between the real space vectors of the chosen SC and PC:', &
                      '                                                                                     ', &
                      '      A[1] = (',M(1,1),')*a[1] + (',M(1,2),')*a[2] + (',M(1,3),')*a[3]               ', & 
                      '      A[2] = (',M(2,1),')*a[1] + (',M(2,2),')*a[2] + (',M(2,3),')*a[3]               ', & 
                      '      A[3] = (',M(3,1),')*a[1] + (',M(3,2),')*a[2] + (',M(3,3),')*a[3]               '
if(commensurate)then
    write(*,'(A)')    '      * The SC and the reference PC are commensurate. Good!'
else
    write(*,*)''
    write(*,'(A)')message_footer
    write(*,'(A)')message_footer2
    write(*,'(A)')'====================================================================================='
endif
write(*,*)''
write(*,*)''

end subroutine print_message_commens_test


subroutine write_attempted_pc_assoc_with_input_unit_cell_and_SC(crystal_pc_reduced_to_prim_cell, &
                                                                crystal_SC_reduced_to_prim_cell, &
                                                                pc_is_prim_cell,SC_is_prim_cell, &
                                                                write_attempted_pc_corresp_to_input_pc, &
                                                                write_attempted_pc_corresp_to_SC)
implicit none
type(crystal_3D), intent(in) :: crystal_pc_reduced_to_prim_cell, crystal_SC_reduced_to_prim_cell
logical, intent(in) :: pc_is_prim_cell, SC_is_prim_cell,  &
                       write_attempted_pc_corresp_to_input_pc, write_attempted_pc_corresp_to_SC

    write(*,*)
    if((.not. SC_is_prim_cell) .and. write_attempted_pc_corresp_to_SC)then
        call write_crystal_to_file(crystal_SC_reduced_to_prim_cell,file_for_SC_reduced_to_prim_cell)
        write(*,'(A)')'>>> We attempted to find a possible primitive cell associated with your SC.'
        write(*,'(3A)')'    Please take a look at the file "',trim(adjustl(file_for_SC_reduced_to_prim_cell)),'".'
    endif
    if((.not. pc_is_prim_cell) .and. write_attempted_pc_corresp_to_input_pc)then
        call write_crystal_to_file(crystal_pc_reduced_to_prim_cell,file_for_pc_reduced_to_prim_cell)
        write(*,'(A)')'>>> Just in case:'
        write(*,'(A)')'    We attempted to find a possible primitive cell associated with your input reference unit cell.'
        write(*,'(3A)')'    Please take a look at the file "',trim(adjustl(file_for_pc_reduced_to_prim_cell)),'".'
    endif

end subroutine write_attempted_pc_assoc_with_input_unit_cell_and_SC


subroutine print_geom_unfolding_relations(GUR,list_SC_kpts_in_wavecar,b_matrix_pc,B_matrix_SC)
implicit none
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR !! Geometric Unfolding Relations
type(vec3d), dimension(:), intent(in) :: list_SC_kpts_in_wavecar
real(kind=dp), dimension(1:3,1:3), intent(in) :: b_matrix_pc,B_matrix_SC
integer :: nkpts, i_SCKPT, i_selc_pcbz_dir,i_needed_dirs,ipc_kpt
integer, dimension(:), allocatable :: n_pckpts_dirs
real(kind=dp), dimension(1:3) :: current_SCKPT, actual_folding_SCKPT, pc_kpt, folding_G


    nkpts = size(list_SC_kpts_in_wavecar)
    allocate(n_pckpts_dirs(1:size(GUR%SCKPT(1)%selec_pcbz_dir(:))))
    do i_selc_pcbz_dir=1,size(GUR%SCKPT(1)%selec_pcbz_dir(:))
        n_pckpts_dirs(i_selc_pcbz_dir) = size(GUR%SCKPT(1)%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(1)%pckpt(:))
    enddo

    write(*,"(A)")"Summary of the points checked on the SCBZ and the pcbz:"
    write(*,'(A,I0)')'    * Total # of SC-KPTS found: ',nkpts
    if((GUR%n_pckpts - sum(n_pckpts_dirs)) > 0)then
        write(*,'(A,I0)')'    * Total # of pc-kpts requested: ', sum(n_pckpts_dirs)
        write(*,'(A,I0)')'    * Total # of complementary pc-kpts determined by the symmetry analysis: ', GUR%n_pckpts - sum(n_pckpts_dirs)
    endif
    write(*,'(A,I0)')'    * Total # of pc-kpts to be checked: ',GUR%n_pckpts
    write(*,*)
    write(*,"(A)")"Summary of the geometric folding relations:"
    write(*,'(2(A,I0),A)')'    * A total of ',GUR%n_folding_pckpts,' pc-kpts (out of the ',GUR%n_pckpts,' checked) satisfy the folding condition:'
    do i_SCKPT=1,nkpts
        do i_selc_pcbz_dir=1,size(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(:))
            do i_needed_dirs=1,size(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(:))
                do ipc_kpt=1, size(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(:))
                    if(GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%folds)then
                        current_SCKPT = list_SC_kpts_in_wavecar(i_SCKPT)%coord(:)
                        actual_folding_SCKPT = &
                            GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords_actual_unfolding_K
                        pc_kpt(:) = GUR%SCKPT(i_SCKPT)%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%coords(:)
                        folding_G(:) = pc_kpt(:) - actual_folding_SCKPT(:)
                        call print_message_pckpt_folds(pc_kpt,current_SCKPT,i_SCKPT,actual_folding_SCKPT,folding_G,b_matrix_pc,B_matrix_SC)
                    endif
                enddo 
            enddo
        enddo
    enddo
    write(*,'(A)')''
    write(*,'(A)')''

end subroutine print_geom_unfolding_relations


subroutine read_energy_info_for_band_search(input_file,e_fermi,E_start,E_end,delta_e)
implicit none
character(len=*), intent(in) :: input_file
real(kind=dp), intent(out) :: e_fermi,E_start,E_end,delta_e
integer :: unt, ios
real(kind=dp) :: E_start_minus_ef, E_end_minus_ef, aux_ener, dE_factor
character(len=str_len) :: str_delta_e
logical :: verbose

    verbose = .TRUE.
    dE_factor = 0.4E-2_dp
    unt = available_io_unit()
    open(unit=unt, file=input_file)
        read(unt,*) e_fermi
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
                write(*,'(A,E9.2,A)')'WARNING: Could not read the size of the energy intervals. Using the default dE = ',dE_factor,'*(Emax - Emin).'
            else
                write(*,'(A,E9.2,A)')'Automatically setting the size of the energy intervals to dE = ',dE_factor,'*(Emax - Emin).' 
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


subroutine print_last_messages_before_unfolding(file_size_in_bytes,nkpts,B_matrix_SC,vbz,E_start,E_end,delta_e,e_fermi)
implicit none
integer(8), intent(in) :: file_size_in_bytes
integer, intent(in) :: nkpts
real(kind=dp), dimension(1:3,1:3), intent(in) :: B_matrix_SC
real(kind=dp), intent(in) :: vbz,E_start,E_end,delta_e,e_fermi
real(kind=dp) :: file_size_in_MB, file_size_in_GB, file_size, &
                 approx_mem_per_kpt_in_bytes, mem_per_kpt_in_MB, mem_per_kpt_in_GB, mem_per_kpt
character(len=2) :: file_size_units, mem_per_kpt_units
real(kind=dp) :: VSBZ
logical :: using_omp

    file_size_in_MB = real(file_size_in_bytes,kind=dp)/(2.0**20)
    file_size_in_GB = real(file_size_in_bytes,kind=dp)/(2.0**30)
    file_size = file_size_in_GB
    file_size_units = 'GB'
    if(file_size_in_GB < 1.0_dp)then
        file_size = file_size_in_MB
        file_size_units = 'MB'
    endif

    approx_mem_per_kpt_in_bytes = real(file_size_in_bytes,kind=dp)/real(nkpts,kind=dp)
    mem_per_kpt_in_MB = approx_mem_per_kpt_in_bytes/(2.0**20)
    mem_per_kpt_in_GB = approx_mem_per_kpt_in_bytes/(2.0**30)
    mem_per_kpt = mem_per_kpt_in_GB
    mem_per_kpt_units = 'GB'
    if(mem_per_kpt_in_GB < 1.0_dp)then
        mem_per_kpt = mem_per_kpt_in_MB
        mem_per_kpt_units = 'MB'
    endif
    write(*,*)
    write(*,'(A,f0.2,X,2A)')'The wavefunction file is ',file_size,file_size_units,' big. Only the necessary data will be read.'
    write(*,'(A,f0.2,X,2A)')'    * Max. of approx. ',mem_per_kpt,mem_per_kpt_units,' at a time.'
    write(*,*)

    VSBZ = dabs(dot_product(B_matrix_SC(1,:),cross(B_matrix_SC(2,:),B_matrix_SC(3,:))))
    if(VSBZ < vbz)then
        write(*,'(A,f0.2,A)') 'The used supercell is ', vbz/VSBZ, ' times bigger than the primitive cell.'
    else
        write(*,'(A)')       'WARNING: The informed primitive cell is not smaller than the supercell found in the WAVECAR file.'
        write(*,'(A,f0.8,A)')'         vpc = ',VSBZ/vbz,' VSC'
    endif

    write(*,'(3(A,f0.5),A)')'Working within the energy interval ',E_start-e_fermi,' < E-EF < ',E_end-e_fermi,' in increments of ',delta_e,' eV'
    write(*,'(A,f0.5,A)')'The Fermi energy EF = ',e_fermi,' will be set to the zero of the energy scale in the output file.'

    using_omp = .FALSE.
    write(*,*)
    write(*,*)
    !$ using_omp = .TRUE.
    !$ write(*,'(A,A,I0,A)')'Some parts of BandUP have been parallelized with OpenMP and will be running using ', &
    !$                      'a maximum of ',omp_get_max_threads(),' thread(s).'
    !$ write(*,'(A)')"You can choose the maximum number of threads by setting the environment variable 'OMP_NUM_THREADS'"
    if(.not. using_omp)then
        write(*,'(A)')"OpenMP shared-memory parallelization has not been enabled (it is optional, don't worry).", &
                      "To enable it, uncomment the line with the flag '-openmp' on the Makefile at BandUP/srci and run build.sh again.", & 
                      "It won't change the results, but it might save you some time."
    endif

    write(*,*)
    !! The following message should always be written to the standad output (also in any modified version of BandUP)
    write(*,'(A)')       '=========================================================================================================================='
    write(*,'(5(A,/),A)')"NOTICE: If you use BandUP or any modified/adapted version/part of it, please don't forget to mention this in your paper!",& 
                         '        You should also read and cite', &
                         '                                                                                                   ', &
                         '  >>>   Paulo V. C. Medeiros, Sven Stafström and Jonas Björk, Phys. Rev. B 89, 041407(R) (2014)', &
                         '                                                                                               ', &
                         '        (http://dx.doi.org/10.1103/PhysRevB.89.041407) and the appropriate references therein.' 
    write(*,'(A)')       '=========================================================================================================================='
    !! End of message
    write(*,*)
    write(*,'(A)')'Band unfolding starts now.'
    write(*,*)

end subroutine print_last_messages_before_unfolding


subroutine print_message_pckpt_folds(pc_kpt,SCKPT,i_SCKPT,actual_folding_SCKPT,folding_G,b_matrix_pc,B_matrix_SC)
implicit none
integer, intent(in) :: i_SCKPT
real(kind=dp), dimension(1:3), intent(in) :: pc_kpt,SCKPT,actual_folding_SCKPT,folding_G
real(kind=dp), dimension(1:3, 1:3), intent(in) :: b_matrix_pc,B_matrix_SC
real(kind=dp), dimension(1:3) :: pckpt_coords, coords
logical :: using_symm_to_get_coeffs
character(len=10) :: str_SCKPT_number

    using_symm_to_get_coeffs = .not. same_vector(actual_folding_SCKPT(:), SCKPT)    
    write(*,*)
    write(*,*)
    str_SCKPT_number = ''
    if(.not. using_symm_to_get_coeffs)then
        write(str_SCKPT_number,"(A,I0,A)")'(',i_SCKPT,')'
    endif
    coords(:) = coords_cart_vec_in_new_basis(cart_vec=actual_folding_SCKPT(:),new_basis=B_matrix_SC)
    write(*,'(A,3(f9.5,A))')'    SCBZ wave-vector K'//trim(adjustl(str_SCKPT_number))//' = ',  &
                                  coords(1),'*B1 + ',coords(2),'*B2 +',coords(3),'*B3'
        
    pckpt_coords(:) = coords_cart_vec_in_new_basis(cart_vec=pc_kpt(:),new_basis=b_matrix_pc)
    write(*,'(A,3(f9.5,A))')'    unfolds onto pcbz wave-vector k = ',  &
                                  pckpt_coords(1),'*b1 + ',pckpt_coords(2),'*b2 + ',pckpt_coords(3),'*b3 = ' 
    pckpt_coords(:) =  coords_cart_vec_in_new_basis(cart_vec=pc_kpt(:),new_basis=B_matrix_SC)
    write(*,'(A,3(f9.5,A))')'                                    = ',  &
                                  pckpt_coords(1),'*B1 + ',pckpt_coords(2),'*B2+ ',pckpt_coords(3),'*B3'

    coords(:) = coords_cart_vec_in_new_basis(cart_vec=folding_G,new_basis=B_matrix_SC)
    write(*,'(A,A,3(I0,A))')'    with the unfolding vector ',  &
                                  'G = ',nint(coords(1)),'*B1 + ',nint(coords(2)),'*B2 + ',nint(coords(3)),'*B3.'
    if(using_symm_to_get_coeffs)then
        coords(:) = coords_cart_vec_in_new_basis(cart_vec=SCKPT(:),new_basis=B_matrix_SC)
        write(*,'(A,/,A,I0,A,3(f9.5,A),/,A,I0,A)')'        *This SCBZ wave-vector belongs to the star of the SCBZ wave-vector',  &
                                                  '                 K(',i_SCKPT,') = ', coords(1),'*B1 + ',coords(2),'*B2 +',coords(3),'*B3.', &
                                                  '         The plane-wave coefficients for K(',i_SCKPT,') will be used.'
    endif

end subroutine print_message_pckpt_folds


subroutine get_SCKPTS_contained_in_wavecar(nkpts,B_matrix_SC,list_SC_kpts_in_wavecar,nearest_neigh_dist_SC_kpts)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
implicit none
integer, intent(in) :: nkpts
real(kind=dp), dimension(1:3,1:3), intent(in) :: B_matrix_SC
type(vec3d), dimension(:), allocatable, intent(out) :: list_SC_kpts_in_wavecar
real(kind=dp), intent(out), optional :: nearest_neigh_dist_SC_kpts
integer :: i_SCKPT, alloc_stat, i
real(kind=dp), dimension(1:3) :: SCKPT_coords
real(kind=dp) :: dist

    ! Reading WAVECAR file
    ! nplane is the number of wave-vectors of the SCBZ that differ from SCBZ_kpt(ikpt) by vectors of the SC rec. latt.
    deallocate(list_SC_kpts_in_wavecar,stat=alloc_stat)
    allocate(list_SC_kpts_in_wavecar(1:nkpts))
    do i_SCKPT=1,nkpts
        call read_from_wavefunc_file(spin_channel=1,i_selected_kpt=i_SCKPT, kpt_frac_coords=SCKPT_coords)
        list_SC_kpts_in_wavecar(i_SCKPT)%coord(:) = SCKPT_coords(1)*B_matrix_SC(1,:) + &
                                                    SCKPT_coords(2)*B_matrix_SC(2,:) + &
                                                    SCKPT_coords(3)*B_matrix_SC(3,:)
    enddo

    if(present(nearest_neigh_dist_SC_kpts).and.(nkpts>1))then
        nearest_neigh_dist_SC_kpts = huge(real(1.0, kind=dp))
        do i_SCKPT=1,nkpts-1
            do i=i_SCKPT+1,nkpts
                dist = norm(list_SC_kpts_in_wavecar(i_SCKPT)%coord(:) - list_SC_kpts_in_wavecar(i)%coord(:))
               if((dist < nearest_neigh_dist_SC_kpts).and.(dist > 5d-4))then
                   nearest_neigh_dist_SC_kpts = dist
               endif 
            enddo
        enddo
    endif

end subroutine get_SCKPTS_contained_in_wavecar


subroutine write_band_struc(out_file,pckpts_to_be_checked,energy_grid,delta_N,EF,zero_of_kpts_scale)
implicit none
character(len=*), intent(in) :: out_file
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
real(kind=dp), dimension(:), intent(in) :: energy_grid
type(UnfoldedQuantitiesForOutput), intent(in) :: delta_N
real(kind=dp), dimension(1:3) :: pckpt, first_pckpt_dir
real(kind=dp), intent(in), optional :: EF, zero_of_kpts_scale
integer :: nener, iener, idir, ndirs, nkpts, ikpt
real(kind=dp) :: e_fermi, origin_of_kpts_line, coord_first_k_in_dir, coord_k
logical :: write_spin_info


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
        write(12,'(A)')trim(adjustl(file_header_BandUP))
        if(write_spin_info)then
            write(12, '(A)')'# Please mind that the support to two-component spinor-like wavefunctions is still under test.' 
            write(12, '(A)')'# Remember this when checking the spin-related unfolded quantities.'
            write(12, '(A, 2(ES10.3,", "), ES10.3, A)')'# The quantization axis used was saxis = [', saxis, ']'
            write(12,'(3(A,X), 2X, 2(A,X))')"#KptCoord", "#E-E_Fermi","#delta_N", "#Spin " // '_|_' // " k", "#Spin " // '//' // " k"
        else
            write(12,'(3(A,X))')"#KptCoord", "#E-E_Fermi","#delta_N"
        endif
        do idir=1,ndirs
            first_pckpt_dir(:) = pckpts_to_be_checked%selec_pcbz_dir(idir)%needed_dir(1)%pckpt(1)%coords(:)
            nkpts = size(pckpts_to_be_checked%selec_pcbz_dir(idir)%needed_dir(1)%pckpt(:))
            do ikpt=1,nkpts
                pckpt(:) = pckpts_to_be_checked%selec_pcbz_dir(idir)%needed_dir(1)%pckpt(ikpt)%coords(:)
                do iener=1,nener
                    coord_k = coord_first_k_in_dir + norm(pckpt(:) - first_pckpt_dir(:))
                    if(write_spin_info)then
                        write(12,'(2(f8.4,2X),3(ES10.3, 2X))')coord_k, energy_grid(iener) - e_fermi, &
                                                              delta_N%pcbz_dir(idir)%pckpt(ikpt)%dN(iener), &
                                                              delta_N%pcbz_dir(idir)%pckpt(ikpt)%spin_proj_perp(iener), &
                                                              delta_N%pcbz_dir(idir)%pckpt(ikpt)%spin_proj_para(iener)
                    else
                        write(12,'(2(f8.4,2X),ES10.3)')coord_k, energy_grid(iener) - e_fermi, &
                                                       delta_N%pcbz_dir(idir)%pckpt(ikpt)%dN(iener)
                    endif
                enddo
            enddo
            coord_first_k_in_dir = coord_k
        enddo
    close(12)

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
character(len=str_len) :: coords_type_line, coords_type_flag, char_n_kpts_dirs, char_line_mode, aux_char, &
                          char_kstart, char_kend
logical :: opt_for_auto_pkpt_search, print_stuff, a0_informed_in_new_format, a0_informed_in_old_format, & 
           give_tip_a0_for_reciprocal, warn_old_format_cartesian

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
        if(ios == 0 .and. trim(adjustl(char_kstart)) /= '') read(03,"(A)",iostat=ios)char_kend
        if(ios == 0 .and. trim(adjustl(char_kstart)) /= ''  .and. trim(adjustl(char_kend)) /= '')then        
            read(03,"(A)",iostat=ios)aux_char
            if((ios == 0 .and. trim(adjustl(aux_char)) == '') .or. ios < 0)then
                ndirs = ndirs + 1
            else
                if(ios==0) ios = 58
            endif
        endif
        if(ios > 0)then
            ios = 58
            write(*,'(A)')'ERROR reading input pc-kpts file. Please check the format.'
            write(*,'(A)')'Stopping now.'
            stop
        endif
    enddo
close(03)
!! Actually reading the file now
allocate(read_k_start(1:ndirs, 1:3), read_k_end(1:ndirs, 1:3))
open(unit=03,file=input_file)
    read(03,*) aux_char
    ! Old style of passing the scaling factor "a0", needed when the k-points are given in cartesian coordinates
    ! a0 is normally not present in the k-points file, but I chose to require it if the k-points are given in 
    ! cartesian coordinates. By doing so, one doesn't need a separate file to specify a0. 
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
            write(*,'(A)')'      * Please check if there are numbers missing or if there are extra non-number characters.'
            write(*,'(A)')'Stopping now.'
            stop
        endif
        read(03,*,iostat=ios)
    enddo
close(03)
!! Reading next the position where the user wants the kpoints string to begin at.
!! This might be useful if you make separate calculations for each pcbz direction
!! and want to combine them in a single file afterwards
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
        if(.not. a0_informed_in_old_format .and. .not. a0_informed_in_new_format)then
            write(*,'(A)')'ERROR: You have selected cartesian coordinates in your input k-points file, but you have not passed a scaling parameter "a0".'
            write(*,'(A)')'       The actuall coordiates of the k-points are given by: ki[actual] = two_pi*ki[passed in file]/a0.'
            write(*,'(A)')'       Please write the value of a0 after your tag "' // trim(adjustl(coords_type_flag))  // '", and run the code again.'
            write(*,'(A)')'Stopping now.'
            stop
        endif
        k_starts(idir,:) = twopi*read_k_start(idir,:)/a0
        k_ends(idir,:) = twopi*read_k_end(idir,:)/a0
    else
        if((upper_case(coords_type_1st_letter) /= 'R').and.(idir==1).and.print_stuff)then
            write(*,*)''
            write(*,'(A)')'WARNING: Assuming that the pc-kpts have been informed in fractional (reciprocal) coordinates.'
            write(*,*)''
        endif
        if(a0_informed_in_old_format .or. a0_informed_in_new_format) give_tip_a0_for_reciprocal = .TRUE.
        do i=1,3
            k_starts(idir,:) = k_starts(idir,:) + read_k_start(idir,i)*b_matrix_pc(i,:)
            k_ends(idir,:) = k_ends(idir,:) + read_k_end(idir,i)*b_matrix_pc(i,:)
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
    !! If the BandUP cannot find a (sequence of) number(s) on "char_n_kpts_dirs", then it assumes you want to 
    !! perform an automatic scan for folding pckpts
    !! This will probably not be used if you use the BandUP's pre-processing tool to get exactly the pckpts you want
    opt_for_auto_pkpt_search = .TRUE.
endif
!! If you want to parse only 1 pckpt, then use it as both the start and the end of the kpts line.
do idir=1,ndirs
    if(norm(k_starts(idir,:)-k_ends(idir,:))<epsilon(1.0_dp))then
        n_kpts_dirs(idir) = 1
    endif
enddo
!! Setting now the number of pckpts along each direction,
!! according to whether the user has actually defined them or not
if(opt_for_auto_pkpt_search)then
    if(print_stuff)then
        write(*,'(A)')'Checking automatically for folding pc-kpts along the directions specified in the input file.'
        write(*,'(A,E10.3,A)')'The automatic scan for pc-kpts will be performed in intervals of ',min_dk,' A^-1.'
    endif
    do idir=1, ndirs
        n_kpts_dirs(idir) = ceiling(1.0_dp + norm(k_ends(idir,:) - k_starts(idir,:))/dabs(min_dk))
    enddo
endif
if(print_stuff)then
    do idir=1, ndirs
        if(opt_for_auto_pkpt_search)then
            write(*,'(2(A,I0))')'      # of pc-kpts requested along pcbz direction ',idir,': ',n_kpts_dirs(idir)
        else
            write(*,'(2(A,I0))')'# of pc-kpts requested along pcbz direction ',idir,': ',n_kpts_dirs(idir)
        endif
    enddo
endif

if(warn_old_format_cartesian)then
    write(*,*)
    write(*,'(A)')'>>> NOTICE: Passing a0 in the first line of the k-points file is deprecated and no longer recommended.'
    if(a0_informed_in_new_format)then
        write(*,'(A)')'            * You have also passed a0 after "' // trim(adjustl(coords_type_flag)) // '".'
        write(*,'(A,f0.5,A)')'            * Only this value will be used (a0 =',a0,').'
    endif
    write(*,'(A)')'            * It will work, but we now recommend that you specify a0 after the tag "cartesian" (separated by a space).'
endif

if(give_tip_a0_for_reciprocal)then
    write(*,*)
    write(*,'(A)')">>> Tip: You don't need to pass the scaling parameter a0 when the k-points are informed in fractional (reciprocal) coordinates."
endif

write(*,*)''

end subroutine read_pckpts_selected_by_user


subroutine print_symm_analysis_for_selected_pcbz_dirs(dirs_req_for_symmavgd_EBS_along_pcbz_dir, & 
                                                      neqv_dirs_pcbz, neqv_dirs_SCBZ, ncompl_dirs, n_irr_compl_dirs)
implicit none
type(irr_bz_directions), dimension(:), intent(in) :: dirs_req_for_symmavgd_EBS_along_pcbz_dir
integer, dimension(:), intent(in) :: neqv_dirs_pcbz,neqv_dirs_SCBZ,ncompl_dirs,n_irr_compl_dirs
integer :: ndirs,idir,n_needed_dirs

ndirs = size(dirs_req_for_symmavgd_EBS_along_pcbz_dir(:))
write(*,*)
write(*,"(A)")"Symmetry analysis for the selected pcbz directions:"
do idir=1,ndirs
    n_needed_dirs = size(dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(:))
    write(*,"(A,I0,A)")"    >>> Direction #",idir,":"
    write(*,"(A,I0,A)")"        * Found ",neqv_dirs_pcbz(idir)," equivalent directions w.r.t. symmetry operations of the pc"
    write(*,"(A,I0,A)")"        * Found ",neqv_dirs_SCBZ(idir)," equivalent directions w.r.t. symmetry operations of the SC"
    if(n_needed_dirs > 1)then
        write(*,"(A,I0,A)")"        * ",ncompl_dirs(idir)," complementary pcbz directions will need to be considered in order to get a proper &
                                                            symmetry-averaged EBS."
        if(ncompl_dirs(idir) /= n_irr_compl_dirs(idir))then
            write(*,"(A,I0,A)")"        * The number of irreducible complementary directions is ",n_irr_compl_dirs(idir),"." 
        endif
    else
        write(*,"(A)")"        * No complementary pcbz directions are needed."
    endif
enddo
write(*,*)

end subroutine print_symm_analysis_for_selected_pcbz_dirs


subroutine print_simm_ops(isym,trans,rot)
implicit none
integer, intent(in) :: isym
integer, intent(in), dimension(1:3) :: trans
integer, intent(in), dimension(1:3,1:3) :: rot

    write(*,'(8X,A,I0,A)')'Symmetry op. #',isym,':'
    write(*,'(12X,A,4(I0,A))')'T[',isym,'] = (',trans(1),', ',trans(2),', ',trans(3),')'
    write(*,'(12X,A,4(I0,A))')'R[',isym,'](1,:) = (',rot(1,1),', ',rot(1,2),', ',rot(1,3),')'
    write(*,'(12X,A,4(I0,A))')'R[',isym,'](2,:) = (',rot(2,1),', ',rot(2,2),', ',rot(2,3),')'
    write(*,'(12X,A,4(I0,A))')'R[',isym,'](3,:) = (',rot(3,1),', ',rot(3,2),', ',rot(3,3),')'
    write(*,*)

end subroutine print_simm_ops


subroutine print_message_pckpt_cannot_be_parsed(stop_when_a_pckpt_cannot_be_parsed)
implicit none
logical, intent(in) :: stop_when_a_pckpt_cannot_be_parsed

    if(stop_when_a_pckpt_cannot_be_parsed)then
        write(*,'(A)')"    ERROR: Could not calculate the spectral weight for the pair (k, K) being parsed right now."
        write(*,'(A)')"           The cause is probably that either:"
        write(*,'(A)')"               * The folding vector G = k - K has been estimated incorrectly by the code (symmetry issue), or"
        write(*,'(A)')"               * K doesn't unfold onto k after all."
        write(*,'(A)')"           Have you used the pre-procesing tool to find the SC-Kpts you needed?"
        write(*,'(A)')"               * If not, please try that. This might solve the problem."
        write(*,'(A)')"               * If yes, please contact me (Paulo)."
        write(*,'(A)')"           Stopping now."
    else
        write(*,'(A)')"    WARNING: Not enough coefficients in the WAVECAR file to unfold this pc wave-vector."
        write(*,'(A)')"             Be careful with your results: they might be incomplete or even wrong."
    endif

end subroutine print_message_pckpt_cannot_be_parsed


subroutine print_message_success_determining_GUR(GUR, stop_if_GUR_fails, is_main_code)
implicit none
type(geom_unfolding_relations_for_each_SCKPT), intent(in) :: GUR !! Geometric Unfolding Relations
logical, intent(in) :: stop_if_GUR_fails, is_main_code

    if(GUR%n_pckpts /= GUR%n_folding_pckpts)then ! if GUR failed
        if(stop_if_GUR_fails)then
            write(*,'(A)')'ERROR: Could not determine all gemetric unfolding relations between the chosen SC-Kpts and pc-kpts.'
        else
            write(*,'(A)')'WARNING: Could not determine all gemetric unfolding relations between the chosen SC-Kpts and pc-kpts.'
            write(*,'(A)')'         The code may fail!'
        endif
        write(*,'(A,I0,A)')'       * ',GUR%n_pckpts,' pc-kpts have been checked.'
        write(*,'(A,I0,A)')'       * ',GUR%n_folding_pckpts,' pc-kpts satisfied the geometric unfolding relations.'
        if(is_main_code)then
            write(*,'(A)')'       Have you used the pre-procesing tool to find the SC-Kpts you needed?'
            write(*,'(A)')'       * If not, please try that. This might solve the problem.'
            write(*,'(A)')'       * If yes, please contact me (Paulo).'
        endif
        if(stop_if_GUR_fails)then
            write(*,'(A)')'       Stopping now.'
        endif
    else
        write(*,'(A)')'The geometric unfolding relations have been successfully determined. Good!'
        write(*,'(A,I0,A)')'    * ',GUR%n_folding_pckpts,' pc-kpts satisfied the geometric unfolding relations.'
    endif

end subroutine print_message_success_determining_GUR


subroutine say_goodbye_and_save_results(delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS, &
                                        pckpts_to_be_checked,energy_grid,e_fermi,zero_of_kpts_scale, &
                                        n_input_pc_kpts,n_folding_pckpts,n_folding_pckpts_parsed)
implicit none
type(UnfoldedQuantitiesForOutput),  intent(in) :: delta_N_only_selected_dirs, delta_N_symm_avrgd_for_EBS
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked
real(kind=dp), dimension(:), intent(in) :: energy_grid
real(kind=dp), intent(in) :: e_fermi, zero_of_kpts_scale
integer, intent(in) :: n_input_pc_kpts,n_folding_pckpts,n_folding_pckpts_parsed


    write(*,*)
    write(*,*)
    write(*,'(A)')'Band unfolding process finished.'
    write(*,'(2(A,I0),A)')'A total of ',n_folding_pckpts,' pc-kpts (out of the ',n_input_pc_kpts,' checked) satisfied the folding condition.'
    if(n_folding_pckpts_parsed > 0)then
        if(n_folding_pckpts_parsed /= n_folding_pckpts)then
            write(*,'(A,I0,A)')'From these points, ', n_folding_pckpts_parsed, ' could be used.'
        endif
        !! Writing delta_N_only_selected_dirs
        call write_band_struc(output_file_only_user_selec_direcs, &
                              pckpts_to_be_checked, energy_grid, &
                              delta_N_only_selected_dirs, &
                              EF=e_fermi, zero_of_kpts_scale=zero_of_kpts_scale)
        !! Writing the delta_N_symm_avrgd_for_EBS
        call write_band_struc(output_file_symm_averaged_EBS, &
                              pckpts_to_be_checked, energy_grid, &
                              delta_N_symm_avrgd_for_EBS, &
                              EF=e_fermi, zero_of_kpts_scale=zero_of_kpts_scale)

        write(*,'(A)')'>>> The unfolded delta_Ns for the EBS strictly along the direction you requested have been saved to the file listed below:'
        write(*,'(2A)')'    * ', trim(adjustl(output_file_only_user_selec_direcs))

        write(*,'(A)')'>>> The symmetry-averaged unfolded delta_Ns for the EBS have been saved to the file listed below:'
        write(*,'(2A)')'    * ', trim(adjustl(output_file_symm_averaged_EBS))
         
        if(zero_of_kpts_scale > 1E-4_dp)then
            write(*,'(A,f8.4,A)')'The zero of the k-points line has been set to ', zero_of_kpts_scale,'.'
        endif
    else
        write(*,'(A)')'No pc-kpts could be parsed. Nothing to be written to the output file.'
    endif

end subroutine say_goodbye_and_save_results


subroutine print_final_times(stime,ftime, &
                       total_time_reading_wavecar,time_calc_spectral_weights,&
                       calc_spec_func,time_calc_spectral_function,&
                       time_spent_calculating_delta_Ns)
implicit none
real(kind=dp), intent(in) :: stime,ftime,total_time_reading_wavecar,time_calc_spectral_weights, &
                             time_calc_spectral_function,time_spent_calculating_delta_Ns
logical, intent(in) :: calc_spec_func
real(kind=dp) :: elapsed_time, time_percentual

    elapsed_time = ftime-stime
    write(*,*)
    write(*,'(A,f0.1,A)')           'Total elapsed time:                                                ',elapsed_time,' s.'

    time_percentual = 100.0_dp*total_time_reading_wavecar/elapsed_time
    if(time_percentual>=1.0_dp)then
        write(*,'(2(A,f0.1),A)')    'Time spent reading the WAVECAR file:                               ',&
                                     total_time_reading_wavecar,'s (',time_percentual,'%).'
    endif
    time_percentual = 100.0_dp*time_calc_spectral_weights/elapsed_time
    if(time_percentual>=1.0_dp)then
        write(*,'(2(A,f0.1),A)')    'Time spent calculating spectral weights:                           ',&
                                     time_calc_spectral_weights,'s (',time_percentual,'%).'
    endif
    if(calc_spec_func)then
        time_percentual = 100.0_dp*time_calc_spectral_function/elapsed_time
        if(time_percentual>=1.0_dp)then
            write(*,'(2(A,f0.1),A)')'Time spent calculating spectral functions:                         ',&
                                     time_calc_spectral_function,'s (',time_percentual,'%).'
        endif
    endif
    time_percentual = 100.0_dp*time_spent_calculating_delta_Ns/elapsed_time
    write(*,'(2(A,f0.1),A)')        'Time spent calculating the delta_Ns:                               ',&
                                     time_spent_calculating_delta_Ns,'s (',time_percentual,'%).'

end subroutine print_final_times


end module io_routines
