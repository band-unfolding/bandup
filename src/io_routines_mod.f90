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

module io_routines
use strings
use general_io
use read_wavecar
use math
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: print_welcome_messages,read_energy_info_for_band_search, read_unit_cell, &
          print_folding_table_results, print_last_messages_before_unfolding, &
          get_SCKPTS_contained_in_wavecar, define_pckpts_to_be_checked, &
          say_goodbye_and_save_results, print_final_times
CONTAINS 


subroutine print_welcome_messages(package_version, version_date)
implicit none
character(len=*), intent(in), optional :: package_version,version_date
write(*,'(A)')'=============================================================================='
write(*,'(A)')'BandUP: Band Unfolding code for Plane-wave based calculations'
if(present(package_version))then
    write(*,'(A)') "Version "//trim(adjustl(package_version))
endif
if(present(version_date))then
    write(*,'(A)') "Date: "//trim(adjustl(version_date))
endif
write(*,'(A)')'=============================================================================='
write(*,'(6(A,/),A)')'Copyright (C) 2013 Paulo V. C. Medeiros',                                        &
                   '                   paume@ifm.liu.se',                                            &
                   '                   Computational Physics Division',                              &
                   '                   Department of Physics, Chemistry and Biology - IFM',          &
                   '                   Linköping University',                                        &
                   '                   Sweden',                                                      & 
                   '                   Visit http://www.ifm.liu.se/theomod/compphys/band-unfolding'
write(*,'(A)')'=============================================================================='
write(*,*)

end subroutine print_welcome_messages


subroutine print_folding_table_results(fold_into_each_other,list_SC_kpts_in_wavecar,pckpts_to_be_checked,b_matrix_pc,B_matrix_SC)
implicit none
logical, dimension(:,:), intent(in) :: fold_into_each_other
type(vec3d), dimension(:), intent(in) :: pckpts_to_be_checked, list_SC_kpts_in_wavecar
real(kind=dp), dimension(1:3,1:3), intent(in) :: b_matrix_pc,B_matrix_SC
integer :: n_input_pc_kpts, nkpts, n_folding_pckpts, &
           i_SCKPT, ipc_kpt
real(kind=dp), dimension(1:3) :: current_SCKPT, pc_kpt, folding_G


n_folding_pckpts = count(fold_into_each_other==.TRUE.)
n_input_pc_kpts = size(pckpts_to_be_checked)
nkpts = size(list_SC_kpts_in_wavecar)
write(*,'(A,I0)')'Total # of pc-kpts to be checked: ',n_input_pc_kpts
write(*,'(A,I0)')'Total # of SC-KPTS in  WAVECAR  : ',nkpts
write(*,'(2(A,I0),A)')'A total of ',n_folding_pckpts,' pc-kpts (out of the ',n_input_pc_kpts,' checked) satisfy the folding condition:'
do i_SCKPT=1,nkpts
    do ipc_kpt=1,n_input_pc_kpts
        if(fold_into_each_other(i_SCKPT,ipc_kpt))then
            current_SCKPT = list_SC_kpts_in_wavecar(i_SCKPT)%coord(:)
            pc_kpt(:) = pckpts_to_be_checked(ipc_kpt)%coord(:)
            folding_G(:) = pc_kpt(:) - current_SCKPT(:)
            call print_message_pckpt_folds(pc_kpt,current_SCKPT,folding_G,b_matrix_pc,B_matrix_SC)
        endif
    enddo
enddo
write(*,'(A)')''
write(*,'(A)')''

end subroutine print_folding_table_results


subroutine read_energy_info_for_band_search(input_file,e_fermi,E_start,E_end,delta_e)
implicit none
character(len=*), intent(in) :: input_file
real(kind=dp), intent(out) :: e_fermi,E_start,E_end,delta_e
integer :: unt, ios
real(kind=dp) :: E_start_minus_ef, E_end_minus_ef, aux_ener, dE_factor
character(len=127) :: str_delta_e
logical :: verbose

    verbose = .TRUE.
    dE_factor = 0.4D-2
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


subroutine print_last_messages_before_unfolding(file_size,B_matrix_SC,vbz,E_start,E_end,delta_e,e_fermi)
implicit none
integer(8), intent(in) :: file_size
real(kind=dp), dimension(1:3,1:3), intent(in) :: B_matrix_SC
real(kind=dp), intent(in) :: vbz,E_start,E_end,delta_e,e_fermi
real(kind=dp) :: VSBZ

    write(*,*)
    write(*,'(A,f0.2,A)')'The wavefunction file is ',file_size/(2.0**30),' GB big. Only the necessary data will be read.'
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

    !$ write(*,*)
    !$ write(*,*)
    !$ write(*,'(A,A,I0,A)')'Some parts of the code have been parallelized with openmp and will be running using ', &
    !$                      'a maximum of ',omp_get_max_threads(),' thread(s).'
    !$ write(*,'(A)')"You can choose the maximum number of threads by setting the environment variable 'OMP_NUM_THREADS'"

    write(*,*)
    !! The following message should always be written to the standad output (in any modified version of BandUP)
    !write(*,'(A)')         '=========================================================================================================================='
    !write(*,'(5(A,/))')'NOTICE: If you use BandUP or any modified version of it, you should read and cite                       ', &
    !                   '                                                                                                   ', &
    !                   '  >>>   Paulo V. C. Medeiros, Sven Stafström and Jonas Björk, Physical Review B: Rapid Communications, ... (2014)', &
    !                   '                                                                                               ', &
    !                   '        and the appropriate references therein.'
    !write(*,'(A)')         '=========================================================================================================================='
    !! End of message
    write(*,*)
    write(*,'(A)')'Band unfolding starts now.'
    write(*,*)

end subroutine print_last_messages_before_unfolding



subroutine print_message_pckpt_folds(pc_kpt,SCKPT,folding_G,b_matrix_pc,B_matrix_SC)
implicit none
real(kind=dp), dimension(1:3), intent(in) :: pc_kpt,SCKPT,folding_G
real(kind=dp), dimension(1:3, 1:3), intent(in) :: b_matrix_pc,B_matrix_SC
real(kind=dp), dimension(1:3) :: pckpt_coords, coords

    write(*,*)
    write(*,*)

    coords(:) = coords_cart_vec_in_new_basis(cart_vec=SCKPT(:),new_basis=B_matrix_SC)
    write(*,'(A,3(f0.5,A))')'SCBZ wave-vector K = ',  &
                                  coords(1),'*B1 + ',coords(2),'*B2 +',coords(3),'*B3'

    pckpt_coords(:) = coords_cart_vec_in_new_basis(cart_vec=pc_kpt(:),new_basis=b_matrix_pc)
    write(*,'(A,3(f0.5,A))')'unfolds onto pcbz wave-vector k = ',  &
                                  pckpt_coords(1),'*b1 + ',pckpt_coords(2),'*b2 + ',pckpt_coords(3),'*b3 = ' 
    pckpt_coords(:) =  coords_cart_vec_in_new_basis(cart_vec=pc_kpt(:),new_basis=B_matrix_SC)
    write(*,'(A,3(f0.5,A))')'                                = ',  &
                                  pckpt_coords(1),'*B1 + ',pckpt_coords(2),'*B2+ ',pckpt_coords(3),'*B3'

    coords(:) = coords_cart_vec_in_new_basis(cart_vec=folding_G,new_basis=B_matrix_SC)
    write(*,'(A,3(I0,A))')'k = K + G, unfolding vector G = ',  &
                                  nint(coords(1)),'*B1 + ',nint(coords(2)),'*B2 + ',nint(coords(3)),'*B3.'

end subroutine print_message_pckpt_folds


subroutine get_SCKPTS_contained_in_wavecar(nkpts,B_matrix_SC,list_SC_kpts_in_wavecar,nearest_neigh_dist_SC_kpts)
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


subroutine write_band_struc(out_file,pc_kpts, E_vs_k, n_kpts_dirs, delta_Ns,EF,zero_of_kpts_scale)
implicit none
character(len=*) :: out_file
type(vec3d), dimension(:), intent(in) :: pc_kpts
real*8, dimension(:,:), intent(in) :: E_vs_k
integer, dimension(:), intent(in) :: n_kpts_dirs
real*8, dimension(:,:), intent(in), optional :: delta_Ns
real*8, intent(in), optional :: EF, zero_of_kpts_scale
integer :: nbands, ndirections, iband, idir, i_1st_kpt_of_dir, ikpt_current_dir, ikpt
real*8 :: coord_k, last_end_coord,e_fermi, origin_of_kpts_line

e_fermi = 0.0d0
if(present(EF))then
     e_fermi = EF
endif
origin_of_kpts_line = 0.0d0
if(present(zero_of_kpts_scale))then
    origin_of_kpts_line = zero_of_kpts_scale
endif

nbands = size(E_vs_k(:,:), dim=1)
ndirections = size(n_kpts_dirs, dim=1)

open(unit=12,file=out_file)
    write(12,'(3(A,X))')"#KptCoord", "#E-E_Fermi","#delta_N"
    do iband=1, nbands
        last_end_coord = origin_of_kpts_line
        i_1st_kpt_of_dir = 1
        do idir=1, ndirections
            if(idir > 1)then
                i_1st_kpt_of_dir = 1 + sum(n_kpts_dirs(1:idir-1))
            endif
            do ikpt_current_dir=1, n_kpts_dirs(idir)
                ikpt = i_1st_kpt_of_dir + ikpt_current_dir - 1
                coord_k = last_end_coord + norm(pc_kpts(ikpt)%coord(:) - pc_kpts(i_1st_kpt_of_dir)%coord(:))
                    write(12,'(2(f8.4,2X),ES10.3)')coord_k, E_vs_k(iband,ikpt)-e_fermi, delta_Ns(iband,ikpt)
            enddo
            last_end_coord = coord_k
        enddo
    enddo
close(12)

end subroutine write_band_struc



subroutine define_pckpts_to_be_checked(input_file, b_matrix_pc, list, n_kpts_dirs, min_dk,verbose,zero_of_kpts_scale)
implicit none
character(len=*), intent(in) :: input_file
real(kind=dp), dimension(1:3,1:3), intent(in) :: b_matrix_pc
type(vec3d), dimension(:), allocatable, intent(out) :: list
integer, intent(out), dimension(:), allocatable :: n_kpts_dirs
real(kind=dp), intent(in), optional :: min_dk
logical, intent(in), optional :: verbose
real(kind=dp), intent(out), optional :: zero_of_kpts_scale
integer :: ndirections, ios, last_ikpt, idir, ikpt, i
character(len=1) :: coords_type_1st_letter
real(kind=dp), dimension(1:3) :: line_vec
real(kind=dp), dimension(:,:), allocatable :: read_k_start, k_start, read_k_end, k_end
real(kind=dp) :: a0, origin_of_kpts_line
character(len=127) :: char_n_kpts_dirs, char_line_mode, aux_char
logical :: opt_for_auto_pkpt_search, print_stuff

print_stuff = .TRUE.
if(present(verbose)) print_stuff = verbose

open(unit=03,file=input_file)
    read(03,*)
    read(03,*)
    read(03,*)
    read(03,*)
    ndirections = 0
    ios=0
    do while (ios==0)
        read(03,*,iostat=ios)
        read(03,*,iostat=ios)
        read(03,*,iostat=ios)
        ndirections = ndirections + 1
    enddo
close(03)

allocate(n_kpts_dirs(1:ndirections))
allocate(read_k_start(1:ndirections, 1:3), read_k_end(1:ndirections, 1:3))
allocate(k_start(1:ndirections, 1:3), k_end(1:ndirections, 1:3))
open(unit=03,file=input_file)
    read(03,*)a0 ! Not normally present in the KPOINTS file, but it is needed here.
    read(03,'(A)')char_n_kpts_dirs
    read(03,'(A)')char_line_mode
    read(03,'(A1)')coords_type_1st_letter
    do idir=1,ndirections
        read(03,*)(read_k_start(idir,i), i=1,3)
        read(03,*)(read_k_end(idir,i), i=1,3)
        read(03,*,iostat=ios)
    enddo
close(03)

char_line_mode = trim(adjustl(char_line_mode))
read(char_line_mode,*,iostat=ios)aux_char, origin_of_kpts_line
if(ios /= 0)then
    origin_of_kpts_line = 0.0d0
endif
if(present(zero_of_kpts_scale))then
    zero_of_kpts_scale = origin_of_kpts_line
endif

do idir=1,ndirections
    k_start(idir,:) = 0.0d0
    k_end(idir,:) = 0.0d0
    if(upper_case(coords_type_1st_letter) == 'R')then
        do i=1,3
            k_start(idir,:) = k_start(idir,:) + read_k_start(idir,i)*b_matrix_pc(i,:)
            k_end(idir,:) = k_end(idir,:) + read_k_end(idir,i)*b_matrix_pc(i,:)
        enddo
    else
        if((upper_case(coords_type_1st_letter) /= 'C').and.(idir==1))then
            if(print_stuff) write(*,'(A)')'Assuming that the pc-kpts have been informed in cartesian coordinates.'
        endif
        k_start(idir,:) = twopi*read_k_start(idir,:)/a0
        k_end(idir,:) = twopi*read_k_end(idir,:)/a0
    endif
enddo

char_n_kpts_dirs = trim(adjustl(char_n_kpts_dirs))
read(char_n_kpts_dirs,*,iostat=ios) n_kpts_dirs(1)
if(ios == 0)then
    opt_for_auto_pkpt_search = .FALSE.
    read(char_n_kpts_dirs,*,iostat=ios)(n_kpts_dirs(idir), idir=1,ndirections)
    if(ios /= 0)then
        n_kpts_dirs(:) = n_kpts_dirs(1)
    endif
else
    opt_for_auto_pkpt_search = .TRUE.
endif

do idir=1,ndirections
    if(norm(k_start(idir,:)-k_end(idir,:))<epsilon(1.0d0))then
        n_kpts_dirs(idir) = 1
    endif
enddo

if(present(min_dk).and.opt_for_auto_pkpt_search)then
    if(print_stuff)then
        write(*,'(A)')'Performing automatic search of folding pc-kpts in the directions specified in the input file.'
        write(*,'(A,E10.3,A)')'The search of pc-kpts will be performed in intervals of ',min_dk,' A^-1.'
    endif
    do idir=1, ndirections
        n_kpts_dirs(idir) = ceiling(1.0d0 + norm(k_end(idir,:) - k_start(idir,:))/dabs(min_dk))
    enddo
else
    if(print_stuff) write(*,'(A)')'Searching folding pc-kpts only among the points specified in the input file.'
endif
if(print_stuff)then
    do idir=1, ndirections
        write(*,'(2(A,I0))')'               # of pc-kpts searched in direction ',idir,': ',n_kpts_dirs(idir)
    enddo
endif

if(allocated(list))then
    deallocate(list)
endif
allocate(list(1:sum(n_kpts_dirs(:))))

last_ikpt = 0
do idir=1, ndirections
    if(n_kpts_dirs(idir)>1)then
        line_vec = (k_end(idir,:) - k_start(idir,:))/real(n_kpts_dirs(idir) - 1, kind=dp)
    else
        line_vec = 0.0d0
    endif
    do ikpt=1,n_kpts_dirs(idir)
        list(last_ikpt + ikpt)%coord(:) = k_start(idir,:) + real(ikpt - 1,kind=dp)*line_vec(:)
    enddo
    last_ikpt = last_ikpt + n_kpts_dirs(idir)
enddo

end subroutine define_pckpts_to_be_checked


subroutine say_goodbye_and_save_results(n_input_pc_kpts,n_folding_pckpts,n_folding_pckpts_parsed, &
                                  pc_E_vs_k, pc_delta_Ns, pckpts_to_be_checked, nkpts_dirs, zero_of_kpts_scale, &
                                  e_fermi, output_file)
implicit none
integer, intent(in) :: n_input_pc_kpts,n_folding_pckpts,n_folding_pckpts_parsed
real(kind=dp), dimension(:,:), allocatable, intent(in) :: pc_E_vs_k, pc_delta_Ns
type(vec3d), dimension(:), intent(in) :: pckpts_to_be_checked
integer, dimension(:) :: nkpts_dirs
real(kind=dp), intent(in) :: e_fermi, zero_of_kpts_scale
character(len=*), intent(in) :: output_file

    write(*,*)
    write(*,*)
    write(*,'(A)')'Band unfolding process finished.'
    write(*,'(2(A,I0),A)')'A total of ',n_folding_pckpts,' pc-kpts (out of the ',n_input_pc_kpts,' checked) satisfied the folding condition.'
    if(n_folding_pckpts_parsed /= n_folding_pckpts)then
        write(*,'(A,I0,A)')'From these points, ',n_folding_pckpts_parsed,' could be used.'
    endif
    if(allocated(pc_E_vs_k).and.allocated(pc_delta_Ns))then
        write(*,'(2A)')'Writing results to file ',trim(adjustl(output_file))
        write(*,'(A,f8.4,A)')'The zero of the k-points line has been set to ',zero_of_kpts_scale,'.'
        call write_band_struc(output_file,pckpts_to_be_checked(:), pc_E_vs_k(:,:), nkpts_dirs, pc_delta_Ns(:,:), EF=e_fermi, &
                              zero_of_kpts_scale=zero_of_kpts_scale)
    else
        write(*,'(A)')'Nothing to be written to the output file.'
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

    time_percentual = 100.0d0*total_time_reading_wavecar/elapsed_time
    if(time_percentual>=1.0d0)then
        write(*,'(2(A,f0.1),A)')    'Time spent reading the WAVECAR file:                               ',&
                                     total_time_reading_wavecar,'s (',time_percentual,'%).'
    endif
    time_percentual = 100.0d0*time_calc_spectral_weights/elapsed_time
    if(time_percentual>=1.0d0)then
        write(*,'(2(A,f0.1),A)')    'Time spent calculating spectral weights:                           ',&
                                     time_calc_spectral_weights,'s (',time_percentual,'%).'
    endif
    if(calc_spec_func)then
        time_percentual = 100.0d0*time_calc_spectral_function/elapsed_time
        if(time_percentual>=1.0d0)then
            write(*,'(2(A,f0.1),A)')'Time spent calculating spectral functions:                         ',&
                                     time_calc_spectral_function,'s (',time_percentual,'%).'
        endif
    endif
    time_percentual = 100.0d0*time_spent_calculating_delta_Ns/elapsed_time
    write(*,'(2(A,f0.1),A)')        'Time spent calculating the delta_Ns:                               ',&
                                     time_spent_calculating_delta_Ns,'s (',time_percentual,'%).'

end subroutine print_final_times

end module io_routines
