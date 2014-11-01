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

!===============================================================================
! MODULE: read_qe_wavefunctions
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Provides a routine to read Quantum Espresso's wavefunctions and g-vectors.
!===============================================================================

module read_qe_wavefunctions
use qexml_module
use cla_wrappers
use general_io
use constants_and_types
use math
implicit none
PRIVATE
PUBLIC :: read_qe_evc_file

CONTAINS

subroutine read_qe_evc_file(wf, args, ikpt, read_coeffs, iostat)
implicit none
! Mandatory input and output variables
type(pw_wavefunction), intent(inout) :: wf
type(comm_line_args), intent(in) :: args
integer, intent(in) :: ikpt
! Optional input and output variables
logical, optional, intent(in) :: read_coeffs
integer, optional, intent(out) :: iostat
! Local variables
integer, parameter :: qe_dp = kind(1.0_dp)
integer, dimension(:,:), allocatable :: G_frac, G_frac_full_grid
complex(kind=qe_dp), allocatable, dimension(:,:) :: pw_coeffs
integer :: nkpts, n_pw, n_bands, n_bands_up, n_bands_down, n_spin, n_spinor
integer :: input_file_unit, ios, alloc_stat
integer :: ik, i, j, ipw, iband
integer, dimension(:), allocatable :: index_G
logical :: is_spinor, two_efs, read_coefficients
real(kind=dp) :: e_fermi, ef_up, ef_dw, alat, encut
real(kind=dp), dimension(1:3,1:3) :: A_matrix, B_matrix
real(kind=dp), dimension(:,:), allocatable :: cart_coords_all_kpts, eigenvales, occupations
character(len=256) :: prefix, outdir, xml_data_file_folder, xml_data_file
character(len=256) :: length_units_kpts, length_units, e_units, encut_units

    if(present(iostat)) iostat = -1
    read_coefficients = .TRUE.
    if(present(read_coeffs)) read_coefficients = read_coeffs

    input_file_unit = available_io_unit()
    outdir = trim(adjustl(args%qe_outdir))
    prefix = trim(adjustl(args%qe_prefix))
    xml_data_file_folder = trim(outdir) // '/' // trim(prefix) // '.save/'
    xml_data_file = trim(xml_data_file_folder) // "data-file.xml"

    ! Initializing the qexml library and opening the xml data-file
    call qexml_init(input_file_unit, dir=xml_data_file_folder) 
    call qexml_openfile(xml_data_file, "read", ierr=ios)
    if(ios/=0) write(*,'(3A)')'ERROR (read_qe_evc_file): Problems opening the XML data-file "',trim(adjustl(xml_data_file)),'"!'
  
    ! Getting some generic info
    call qexml_read_bz(num_k_points=nkpts, ierr=ios)
    ! Reading lattice data
    call qexml_read_cell(alat=alat, a_units=length_units, &
                         a1=A_matrix(1,:), a2=A_matrix(2,:), a3=A_matrix(3,:), &
                         ierr=ios)
    if(ios/=0) write(*,'(A)')"ERROR (read_qe_evc_file): Could not read lattice info!"
    alat = to_angstrom(alat, length_units)
    do j=1,3
        do i=1,3
            A_matrix(i,j) = to_angstrom(A_matrix(i,j), length_units)
        enddo
    enddo
    length_units = 'Angstrom'
    call get_rec_latt(latt=A_matrix, rec_latt=B_matrix)

    allocate(cart_coords_all_kpts(1:3, 1:nkpts))
    call qexml_read_bz(xk=cart_coords_all_kpts, k_units=length_units_kpts, ierr=ios)
    cart_coords_all_kpts(:,:) = (twopi/alat) * cart_coords_all_kpts(:,:)

    call qexml_read_bands_info(nbnd=n_bands, nbnd_up=n_bands_up, nbnd_down=n_bands_down, &
                               nspin=n_spin, noncolin=is_spinor, &
                               ef=e_fermi, two_fermi_energies=two_efs, &
                               ef_up=ef_up, ef_dw=ef_dw, &
                               energy_units=e_units, ierr=ios)

    allocate(eigenvales(1:n_bands, 1:nkpts), occupations(1:n_bands, 1:nkpts))
    call qexml_read_bands_pw(num_k_points=nkpts, nkstot=nkpts, lsda=(n_spin==2), nbnd=n_bands, &
                             lkpoint_dir=.TRUE., filename='default', &
                             et=eigenvales, wg=occupations , ierr=ios)

    call qexml_read_planewaves(ecutwfc=encut, cutoff_units=encut_units, ierr=ios)
    if(ios/=0) write(*,'(A)')"ERROR (read_qe_evc_file): Could not read general info regarding wavefunctions!"

    n_spinor = 1
    if(is_spinor) n_spinor = 2

    if(read_coefficients)then
        ! Getting number of pw for current kpt
        call qexml_read_gk(ik=ikpt, npwk=n_pw, ierr=ios)
        if(ios/=0) write(*,'(A)')"ERROR (read_qe_evc_file): Could not determine the number of plane-waves for the current Kpt."

        ! Allocating arrays
        allocate(G_frac(1:3,1:n_pw), stat=alloc_stat)
        if(alloc_stat/=0) write(*,'(A)')"ERROR (read_qe_evc_file): Could not allocate G-vectors for the current Kpt."
        allocate(pw_coeffs(1:n_spinor*n_pw, 1:n_bands), stat=alloc_stat)
        if(alloc_stat/=0) write(*,'(A)')"ERROR (read_qe_evc_file): Could not allocate plane-wave coefficients for the current Kpt."

        ! Reading the G vectors for the current kpt
        call qexml_read_gk(ik=ikpt, igk=G_frac, ierr=ios)
        if(ios/=0) write(*,'(A)')"ERROR (read_qe_evc_file): Could not read G-vectors for the current Kpt."

        ! Reading pw coeffs
        call qexml_read_wfc(ibnds=1, ibnde=n_bands, ik=ikpt, ispin=wf%i_spin, wf=pw_coeffs, ierr=ios)
        if(ios/=0) call qexml_read_wfc(ibnds=1, ibnde=n_bands, ik=ikpt, wf=pw_coeffs, ierr=ios)
        if(ios/=0) write(*,'(A)')"ERROR (read_qe_evc_file): Could not read plane-wave coefficients for the current Kpt."
    endif

    ! Done. Closing file now.
    call qexml_closefile("read", ierr=ios)
    if(ios/=0) write(*,'(A)')"WARNING (read_qe_evc_file): Error closing xml data-file."


    ! Transferring data to wf
    wf%is_spinor = is_spinor
    wf%n_spin = n_spin
    wf%n_bands = n_bands
    wf%n_bands_up = n_bands_up
    wf%n_bands_down = n_bands_down
    wf%n_pw = n_pw

    wf%encut = to_ev(encut, encut_units)
    wf%e_fermi = to_ev(e_fermi, e_units)
    wf%e_fermi_up = to_ev(ef_up, e_units)
    wf%e_fermi_down = to_ev(ef_dw, e_units)
    wf%two_efs = two_efs
    if(.not. two_efs)then
        wf%e_fermi_up = wf%e_fermi
        wf%e_fermi_down = wf%e_fermi
    endif

    deallocate(wf%band_occupations, wf%band_energies, stat=alloc_stat)
    allocate(wf%band_occupations(1:n_bands), wf%band_energies(1:n_bands), stat=alloc_stat)
    wf%band_occupations(:) = occupations(:,ikpt)
    do j=1, size(eigenvales, dim=1)
        wf%band_energies(j) = to_ev(eigenvales(j,ikpt), e_units)
    enddo

    wf%A_matrix = A_matrix
    wf%B_matrix = B_matrix
    wf%Vcell = abs(triple_product(A_matrix(1,:), A_matrix(2,:), A_matrix(3,:)))
    wf%kpt_cart_coords = cart_coords_all_kpts(:,ikpt)
    wf%kpt_frac_coords = coords_cart_vec_in_new_basis(wf%kpt_cart_coords, new_basis=B_matrix)

    if(read_coefficients)then
        ! Transferring G-vecs
        deallocate(wf%G_frac, wf%G_cart, stat=alloc_stat)
        allocate(wf%G_frac(1:n_pw), wf%G_cart(1:n_pw))
        do ipw=1, n_pw
            wf%G_frac(ipw)%coord(:) = G_frac(:,ipw)
            wf%G_cart(ipw)%coord(:) = 0.0_dp
            do j=1,3
                wf%G_cart(ipw)%coord(:) = wf%G_cart(ipw)%coord(:) + &
                                          real(G_frac(j,ipw), kind=dp)*B_matrix(j,:)
            enddo
        enddo
        deallocate(G_frac, stat=alloc_stat)

        ! Transferring pw coeffs
        deallocate(wf%pw_coeffs, stat=alloc_stat)
        allocate(wf%pw_coeffs(1:n_spinor, 1:n_pw, 1:n_bands), stat=alloc_stat)
        if(wf%is_spinor)then
            wf%pw_coeffs(1,:,:) = pw_coeffs(1:n_pw/2,:)
            wf%pw_coeffs(2,:,:) = pw_coeffs(n_pw/2+1:,:)
        else
            wf%pw_coeffs(1,:,:) = pw_coeffs(:,:)
        endif
        deallocate(pw_coeffs, stat=alloc_stat)
    endif

end subroutine read_qe_evc_file


end module read_qe_wavefunctions
