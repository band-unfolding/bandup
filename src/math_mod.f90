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
module math
use spglib_f08
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: sp, dp, pi, twopi, min_dk, default_tol_for_vec_equality, default_tol_for_int_commens_test, &
          default_symprec, max_tol_for_vec_equality, vec3d, &
          star, time, n_digits_integer, cross, norm, angle, real_seq, delta, integral_delta_x_minus_x0, &
          symmetry_operation, get_rec_latt, coords_cart_vec_in_new_basis, vec_in_latt, &
          reduce_point_to_bz, same_vector, get_symm, get_star, get_irr_kpts, &
          pt_eqv_by_point_group_symop, bz_direction, eqv_bz_directions, irr_bz_directions, &
          get_all_irr_dirs_req_for_symmavgd_EBS, kpts_line, selected_pcbz_directions, &
          geom_unfolding_relations_for_each_SCKPT, delta_Ns, allocate_delta_Ns_type, & 
          delta_Ns_for_output, allocate_delta_Ns_for_output_type, check_if_pc_and_SC_are_commensurate, &
          crystal_3D, create_crystal, append, operator(+)

interface append
  module procedure append_integer_list, append_character_list
end interface append

interface operator(+)
    module procedure sum_delta_Ns_for_output_type
end interface operator(+)

!! Constants used throughout the code
integer, parameter :: sp = selected_real_kind(6, 37)    ! Single precision
integer, parameter :: dp = selected_real_kind(15, 307)  ! Double precision
real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp), twopi = 2.0_dp*pi
real(kind=dp), parameter :: default_tol_for_vec_equality=1E-5_dp, &
                            max_tol_for_vec_equality=1E-3_dp, &
                            min_dk=2.0_dp*default_tol_for_vec_equality, &
                            default_tol_for_int_commens_test=1E-5_dp, & 
                            default_symprec=1E-6_dp
real(kind=dp), parameter :: two_m_over_hbar_sqrd = 0.262465831 ! c = 2m/hbar**2 in units of 1/eV Ang^2 (from WaveTrans)

!! Derived type definitions
type :: vec3d
  real(kind=dp), dimension(1:3) :: coord
end type vec3d

type :: crystal_3D
    character(len=256) :: description=''
    real(kind=dp), dimension(1:3,1:3) :: latt_vecs, rec_latt_vecs
    real(kind=dp) :: vol, rec_latt_vol
    real(kind=dp), dimension(:,:), allocatable :: coords_basis_atoms, & ! Shall be dimension(:,1:3)
                                                  fractional_coords_basis_atoms
    logical, dimension(:,:), allocatable :: unconstrained_dof_basis_atoms
    character(len=3), dimension(:), allocatable :: atomic_symbols_basis_atoms
    integer, dimension(:), allocatable :: integer_types_basis_atoms
end type crystal_3D

type :: symmetry_operation
    integer, dimension(1:3) :: translation_fractional_coords
    integer, dimension(1:3,1:3) :: rotation_fractional_coords
    real(kind=dp), dimension(1:3) :: translation_cartesian_coords
    real, dimension(1:3,1:3) :: rotation_cartesian_coords
    real(kind=dp), dimension(1:3,1:3) :: basis
end type symmetry_operation

type :: star_point_properties
    real(kind=dp), dimension(1:3) :: coord
    integer :: symop  ! The index of the symmetry operation
end type star_point_properties

type :: star
   integer :: neqv
   type(star_point_properties), dimension(:), allocatable :: eqv_pt
end type star

type :: bz_direction
    real(kind=dp), dimension(1:3) :: kstart, kend
    integer :: neqv=1
    real(kind=dp) :: weight=0.0_dp
end type bz_direction

type :: eqv_bz_directions
    integer :: neqv=1 ! Number of eqv. directions in the set.
    type(bz_direction), dimension(:), allocatable :: eqv_dir
end type eqv_bz_directions

type :: irr_bz_directions
    type(bz_direction), dimension(:), allocatable :: irr_dir
end type irr_bz_directions

type :: trial_folding_pckpt
    ! Message from Paulo:
    ! The prefix "S" means "symmetrized". This is a little trick I came up with
    ! that allows me to use the coefficients of a SC wavefunction psi(K',n) to
    ! calculate the spectral weights associated with a SC wavefunction psi(K,n),
    ! where K' = SK and S is a symmetry operation of the crystal's point group.
    real(kind=dp), dimension(1:3) :: coords_actual_unfolding_K,coords,Scoords,Sfolding_vec
    logical :: folds
end type trial_folding_pckpt

type :: list_of_trial_folding_pckpts
    type(trial_folding_pckpt), dimension(:), allocatable :: pckpt
end type list_of_trial_folding_pckpts

type :: needed_pcbz_dirs_for_EBS
    type(list_of_trial_folding_pckpts), dimension(:), allocatable :: needed_dir
end type needed_pcbz_dirs_for_EBS

type :: selected_pcbz_directions
    type(needed_pcbz_dirs_for_EBS), dimension(:), allocatable :: selec_pcbz_dir
end type selected_pcbz_directions

type :: geom_unfolding_relations_for_each_SCKPT
    type(selected_pcbz_directions), dimension(:), allocatable :: SCKPT
    logical, dimension(:), allocatable :: SCKPT_used_for_unfolding
    integer :: n_pckpts,n_folding_pckpts
end type geom_unfolding_relations_for_each_SCKPT

!! Defining derived types to support a "delta_Ns" type.
type :: delta_Ns_and_SFs_at_every_energy_grid_point_for_a_pckpt
    real(kind=dp), dimension(:), allocatable :: dN,SF  !! delta_Ns and spectral functions at each point of the energy grid
end type delta_Ns_and_SFs_at_every_energy_grid_point_for_a_pckpt

type :: list_of_pckpts_for_SFs_and_delta_Ns
    type(delta_Ns_and_SFs_at_every_energy_grid_point_for_a_pckpt), dimension(:), allocatable :: pckpt
end type list_of_pckpts_for_SFs_and_delta_Ns

type :: delta_N_and_SF_info_for_needed_pcbz_dirs_for_EBS
    type(list_of_pckpts_for_SFs_and_delta_Ns), dimension(:), allocatable :: needed_dir
end type delta_N_and_SF_info_for_needed_pcbz_dirs_for_EBS
!! Defining now the "delta_Ns" type, which can also hold info about spectral functions SFs
!! The structure of a variable delta_N of type(delta_Ns) is:
!! delta_N%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(ipc_kpt)%[ dN(iener) or SF(iener)]
type :: delta_Ns
    type(delta_N_and_SF_info_for_needed_pcbz_dirs_for_EBS), dimension(:), allocatable :: selec_pcbz_dir
end type delta_Ns

type :: delta_Ns_for_output
    type(list_of_pckpts_for_SFs_and_delta_Ns), dimension(:), allocatable :: pcbz_dir
end type delta_Ns_for_output


!! Functions and subroutines
CONTAINS

subroutine append_integer_list(item, list)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
integer, intent(in) :: item
integer, dimension(:), allocatable, intent(inout) :: list
integer, dimension(:), allocatable :: new_list

if(.not.allocated(list))then
    allocate(list(1:1))
    list(1) = item
else
    allocate(new_list(1:size(list)+1))
    new_list(1:size(list)) = list
    new_list(size(list)+1) = item
    deallocate(list)
    allocate(list(1:size(new_list)))
    list = new_list
    deallocate(new_list)
endif

end subroutine append_integer_list

subroutine append_character_list(item, list)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
character(len=*), intent(in) :: item
character(len=*), dimension(:), allocatable, intent(inout) :: list
character(len=len(list)), dimension(:), allocatable :: new_list

if(.not.allocated(list))then
    allocate(list(1:1))
    list(1) = item
else
    allocate(new_list(1:size(list)+1))
    new_list(1:size(list)) = list
    new_list(size(list)+1) = item
    deallocate(list)
    allocate(list(1:size(new_list)))
    list = new_list
    deallocate(new_list)
endif

end subroutine append_character_list



subroutine create_crystal(crystal, description, latt_vecs, coords_basis_atoms, &
                          atomic_symbols_basis_atoms, unconstrained_dof_basis_atoms)
!! Copyright (C) 2014 Paulo V. C. Medeiros
type(crystal_3D), intent(out), allocatable :: crystal
character(len=*), intent(in), optional :: description
real(kind=dp), dimension(1:3,1:3), intent(in) :: latt_vecs
real(kind=dp), dimension(:,:), intent(in), optional :: coords_basis_atoms ! dimension(:,1:3)
character(len=3), dimension(:), optional :: atomic_symbols_basis_atoms
logical, dimension(:,:), optional :: unconstrained_dof_basis_atoms ! dimension(:,1:3)
integer :: alloc_stat, iatom, isymb
real(kind=dp), dimension(1:3) :: default_positions
logical, dimension(1:1,1:3) :: default_unconstrained_dof
real(kind=dp), dimension(:,:), allocatable :: atom_pos ! dimension(:,1:3)
character(len=3), dimension(:), allocatable :: symbols, unique_symbols
logical, dimension(:,:), allocatable :: dof_basis_atoms ! dimension(:,1:3)
   
    ! Defaults of optional arguments 
    default_positions(:) = 0.0_dp * (latt_vecs(1,:) + latt_vecs(2,:) + latt_vecs(3,:))
    default_unconstrained_dof(1,:) = (/.TRUE.,.TRUE.,.TRUE./) 

    ! 'Cleaning' crystal
    deallocate(crystal,stat=alloc_stat)
    allocate(crystal)

    ! Handling optional arguments
    if(present(coords_basis_atoms))then
        allocate(atom_pos(1:size(coords_basis_atoms,dim=1), &
                          1:size(coords_basis_atoms,dim=2)))
        atom_pos = coords_basis_atoms
    else
        allocate(atom_pos(1:1,1:3))
        atom_pos(1,:) = default_positions
    endif
    if(present(unconstrained_dof_basis_atoms))then
        allocate(dof_basis_atoms(1:size(unconstrained_dof_basis_atoms,dim=1), &
                                 1:size(unconstrained_dof_basis_atoms,dim=2)))
        dof_basis_atoms = unconstrained_dof_basis_atoms
    else
        allocate(dof_basis_atoms(1:1,1:3))
        dof_basis_atoms = default_unconstrained_dof
    endif
    if(present(atomic_symbols_basis_atoms))then
        allocate(symbols(1:size(atomic_symbols_basis_atoms)))
        symbols = atomic_symbols_basis_atoms
    else
        allocate(symbols(1:size(atom_pos,dim=1)))
        symbols(1) = ''
    endif
    if(present(description))then
        crystal%description = trim(adjustl(description))
    else
        crystal%description = ''
    endif
    ! End of handling optional arguments
    
    ! Lattice
    crystal%latt_vecs = latt_vecs
    crystal%vol = abs(triple_product(latt_vecs(1,:), latt_vecs(2,:), latt_vecs(3,:)))
    call get_rec_latt(latt_vecs,crystal%rec_latt_vecs,crystal%rec_latt_vol)
    ! Basis
    allocate(crystal%coords_basis_atoms(1:size(atom_pos,dim=1), &
                                        1:size(atom_pos,dim=2)))
    crystal%coords_basis_atoms = atom_pos
    allocate(crystal%fractional_coords_basis_atoms(1:size(atom_pos,dim=1), &
                                                   1:size(atom_pos,dim=2)))
    do iatom=1, size(atom_pos,dim=1)
        crystal%fractional_coords_basis_atoms(iatom,:) = &
            coords_cart_vec_in_new_basis(atom_pos(iatom,:), latt_vecs(:,:))
    enddo
    ! Atomic symbols
    allocate(crystal%atomic_symbols_basis_atoms(1:size(symbols)))
    crystal%atomic_symbols_basis_atoms = symbols
    call append(item=symbols(1), list=unique_symbols)
    if(size(symbols)>1)then
        do iatom=2, size(symbols)
            if(all(unique_symbols /= symbols(iatom)))then
                call append(item=symbols(iatom), list=unique_symbols)
            endif
        enddo
    endif
    allocate(crystal%integer_types_basis_atoms(1:size(symbols)))
    do iatom=1, size(symbols)
        do isymb=1, size(unique_symbols)
            if(symbols(iatom) == unique_symbols(isymb))then
                crystal%integer_types_basis_atoms(iatom) = isymb
                exit
            endif
        enddo
    enddo

    ! Degrees of freedom
    allocate(crystal%unconstrained_dof_basis_atoms(1:size(dof_basis_atoms,dim=1), &
                                                   1:size(dof_basis_atoms,dim=2)))
    crystal%unconstrained_dof_basis_atoms = dof_basis_atoms


end subroutine create_crystal



subroutine check_if_pc_and_SC_are_commensurate(commensurate,M,b_matrix_pc,B_matrix_SC,tol)
!! Copyright (C) 2014 Paulo V. C. Medeiros
! Calculates the matrix M so that  
! A[i] = sum(M_{ij}*a[j]; j=1,2,3),
! where 'A' and 'a' are the real space lattice vectors of the SC and pc, respectively. 
! The subroutine then uses M to check if the PC and SC are commensurate.
! I decided to work with the reciprocal lattice vectors instead just because it makes it easier to 
! incorporate the subroutine into the main code of BandUP.
! The only additional op. is a 3x3 matrix transposition, which is really not expensive.
implicit none
logical, intent(out) :: commensurate
real(kind=dp), dimension(1:3,1:3), intent(out) :: M
real(kind=dp), dimension(1:3,1:3), intent(in) :: b_matrix_pc,B_matrix_SC
real(kind=dp), intent(in), optional :: tol
integer, dimension(1:3,1:3) :: int_M
real(kind=dp), dimension(1:3,1:3) :: residue_M
real(kind=dp) :: max_residue

    commensurate = .TRUE.
    max_residue = default_tol_for_int_commens_test
    if(present(tol))then
        max_residue = dabs(tol)
    endif
    M = transpose(matmul(b_matrix_pc,inverse_of_3x3_matrix(B_matrix_SC)))
    int_M = nint(M)
    residue_M = dabs(M - int_M)
    if(any(residue_M > max_residue))then
        commensurate = .FALSE.
    endif

end subroutine check_if_pc_and_SC_are_commensurate



subroutine allocate_delta_Ns_type(delta_N,pckpts_to_be_checked)
implicit none
type(delta_Ns), intent(out) :: delta_N
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked !! Geometric Unfolding Relations
integer :: n_selec_pcbz_dir,i_selc_pcbz_dir,n_needed_dirs,i_needed_dirs,nkpts

    n_selec_pcbz_dir = size(pckpts_to_be_checked%selec_pcbz_dir(:))
    allocate(delta_N%selec_pcbz_dir(1:n_selec_pcbz_dir))
    do i_selc_pcbz_dir=1,n_selec_pcbz_dir
        n_needed_dirs = size(pckpts_to_be_checked%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(:))
        allocate(delta_N%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(1:n_needed_dirs))
        do i_needed_dirs=1,n_needed_dirs
            nkpts = size(pckpts_to_be_checked%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(:))
            allocate(delta_N%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(i_needed_dirs)%pckpt(1:nkpts))
        enddo
    enddo

end subroutine allocate_delta_Ns_type

subroutine allocate_delta_Ns_for_output_type(delta_N,pckpts_to_be_checked)
implicit none
type(delta_Ns_for_output), intent(out) :: delta_N
type(selected_pcbz_directions), intent(in) :: pckpts_to_be_checked !! Geometric Unfolding Relations
integer :: n_selec_pcbz_dir,i_selc_pcbz_dir,nkpts

    n_selec_pcbz_dir = size(pckpts_to_be_checked%selec_pcbz_dir(:))
    allocate(delta_N%pcbz_dir(1:n_selec_pcbz_dir))
    do i_selc_pcbz_dir=1,n_selec_pcbz_dir
        nkpts = size(pckpts_to_be_checked%selec_pcbz_dir(i_selc_pcbz_dir)%needed_dir(1)%pckpt(:))
        allocate(delta_N%pcbz_dir(i_selc_pcbz_dir)%pckpt(1:nkpts))
    enddo

end subroutine allocate_delta_Ns_for_output_type


function sum_delta_Ns_for_output_type(delta_N1, delta_N2) result(delta_N_sum)
implicit none
type(delta_Ns_for_output) :: delta_N_sum
type(delta_Ns_for_output), intent(in) :: delta_N1, delta_N2
integer :: idir, ikpt, n_ener
    
    allocate(delta_N_sum%pcbz_dir(1:size(delta_N1%pcbz_dir)))
    do idir=1, size(delta_N1%pcbz_dir)
        allocate(delta_N_sum%pcbz_dir(idir)%pckpt(1:size(delta_N1%pcbz_dir(idir)%pckpt)))
        do ikpt=1, size(delta_N1%pcbz_dir(idir)%pckpt)
            n_ener = size(delta_N1%pcbz_dir(idir)%pckpt(ikpt)%dN(:))
            allocate(delta_N_sum%pcbz_dir(idir)%pckpt(ikpt)%dN(1:n_ener))
            delta_N_sum%pcbz_dir(idir)%pckpt(ikpt)%dN(:) = delta_N1%pcbz_dir(idir)%pckpt(ikpt)%dN(:) + &
                                                           delta_N2%pcbz_dir(idir)%pckpt(ikpt)%dN(:)
            if(allocated(delta_N1%pcbz_dir(idir)%pckpt(ikpt)%SF))then
                n_ener = size(delta_N1%pcbz_dir(idir)%pckpt(ikpt)%SF(:))
                allocate(delta_N_sum%pcbz_dir(idir)%pckpt(ikpt)%SF(1:n_ener))
                delta_N_sum%pcbz_dir(idir)%pckpt(ikpt)%SF(:) = delta_N1%pcbz_dir(idir)%pckpt(ikpt)%SF(:) + &
                                                               delta_N2%pcbz_dir(idir)%pckpt(ikpt)%SF(:)
            endif
        enddo
    enddo
    return

end function sum_delta_Ns_for_output_type


function time() result(rtn)
implicit none
real(kind=dp) :: rtn
integer :: count, count_rate

!$ rtn = omp_get_wtime()
!$ return
   CALL system_clock(count, count_rate)
   rtn = count/count_rate

end function time


function n_digits_integer(int_number, add_one_if_negative) result(n_digits)
implicit none
integer :: n_digits
integer, intent(in) :: int_number
logical, intent(in), optional :: add_one_if_negative
logical :: count_minus_sign

    count_minus_sign = .FALSE.
    if(present(add_one_if_negative)) count_minus_sign = add_one_if_negative

    n_digits = 0
    do 
        n_digits = n_digits + 1
        if (mod(int_number, 10**n_digits) == int_number) exit
    enddo
    if(count_minus_sign .and. (int_number<0))then
        n_digits = n_digits + 1
    endif

end function n_digits_integer


function cross(b,c) result(a)
implicit none
  real(kind=dp), dimension(1:3) :: a, b, c

  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)

end function cross

function triple_product(a,b,c) result(triple_prod)
implicit none
real(kind=dp) :: triple_prod
real(kind=dp), dimension(1:3), intent(in) :: a, b, c

    triple_prod = dot_product(a,cross(b,c))

end function triple_product

function norm(v) result(rtn)
implicit none
  real(kind=dp), dimension(1:3) :: v
  real(kind=dp) :: rtn

  rtn = dsqrt(dot_product(v,v))

end function norm

function angle(v1, v2) result(rtn)
implicit none
  real(kind=dp), dimension(1:3) :: v1, v2
  real(kind=dp) :: rtn

  rtn = acos(dot_product(v1, v2)/(norm(v1)*norm(v2)))

end function angle

subroutine real_seq(first_term,last_term,increment,return_list)
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: return_list
real(kind=dp), intent(in) :: first_term,last_term,increment
integer :: i, n_terms

    n_terms = nint((last_term - first_term)/increment + 1.0_dp)
    allocate(return_list(1:n_terms))
    do i = 1, n_terms
         return_list(i) = first_term + real(i - 1, kind=dp)*increment
    enddo

end subroutine real_seq

function kpts_line(kstart,kend,nk) result(line)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
type(vec3d), dimension(:), allocatable :: line
real(kind=dp), dimension(1:3), intent(in) :: kstart, kend
integer, intent(in) :: nk
integer :: alloc_stat,ikpt
real(kind=dp), dimension(1:3) :: line_vec

    deallocate(line, stat=alloc_stat)
    allocate(line(1:nk))
    if(nk>1)then
        line_vec = (kend - kstart)/real(nk - 1, kind=dp)
    else
        line_vec = 0.0_dp
    endif
    do ikpt=1,nk
        line(ikpt)%coord(:) = kstart + real(ikpt - 1,kind=dp)*line_vec(:)
    enddo
    return

end function kpts_line

function delta(x,FWHM,std_dev) result(rtn)
implicit none
! It is actually a Gaussian function with a default FWHM = def_FWHM 
real(kind=dp) :: rtn, stddev
real(kind=dp), parameter :: def_FWHM = 0.05
real(kind=dp), intent(in) :: x
real(kind=dp), intent(in), optional :: std_dev, FWHM

    stddev = def_FWHM/2.35482
    if(present(FWHM))then
        stddev = dabs(FWHM)/2.35482
    else
        if(present(std_dev))then
            stddev = dabs(std_dev)
        endif
    endif

    rtn = exp(-((x/stddev)**2.0_dp)/2.0_dp)/(stddev*dsqrt(twopi))

end function delta


function integral_delta_x_minus_x0(x0,lower_lim,upper_lim,std_dev) result(rtn)
! Consistently with the employed approx. by a Gaussian distribution
implicit none
real(kind=dp) :: rtn
real(kind=dp), intent(in) :: x0,lower_lim,upper_lim,std_dev
real(kind=dp) :: x1, x2


    x1 = (lower_lim - x0)/(sqrt(2.0_dp)*std_dev)
    x2 = (upper_lim - x0)/(sqrt(2.0_dp)*std_dev)
    rtn = 0.5_dp*(erf(x2) - erf(x1))

end function integral_delta_x_minus_x0


function inverse_of_3x3_matrix(matrix, success) result(inverse_matrix)
implicit none
real(kind=dp), dimension(3,3) :: inverse_matrix
real(kind=dp), dimension(3,3), intent(in) :: matrix
logical, intent(out), optional :: success
real(kind=dp), dimension(3,3) :: cofactor
real(kind=dp) :: det
real(kind=dp), parameter :: eps = 1.0E-10_dp

if(present(success))then
  success = .FALSE.
endif
cofactor(1,1) =  (matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))
cofactor(1,2) = -(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))
cofactor(1,3) =  (matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
! Cofactor expansion along the 1st row
det = dot_product(matrix(1,:), cofactor(1,:))
if(abs(det) <= eps)then
   inverse_matrix = 0.0D0
else
    cofactor(2,1) = -(matrix(1,2)*matrix(3,3)-matrix(1,3)*matrix(3,2))
    cofactor(2,2) =  (matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))
    cofactor(2,3) = -(matrix(1,1)*matrix(3,2)-matrix(1,2)*matrix(3,1))
    cofactor(3,1) =  (matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))
    cofactor(3,2) = -(matrix(1,1)*matrix(2,3)-matrix(1,3)*matrix(2,1))
    cofactor(3,3) =  (matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))

    inverse_matrix = transpose(cofactor)/det
    if(present(success))then
       success = .TRUE.
    endif
endif

return

end function inverse_of_3x3_matrix

subroutine get_rec_latt(latt,rec_latt,rec_latt_vol)
implicit none
real(kind=dp), dimension(1:3,1:3), intent(in) :: latt
real(kind=dp), dimension(1:3,1:3), intent(out) :: rec_latt
real(kind=dp), intent(out), optional :: rec_latt_vol
real(kind=dp) :: vpc

    vpc = dabs(dot_product(latt(1,:), cross(latt(2,:),latt(3,:))))
    rec_latt(1,:) = twopi*cross(latt(2,:),latt(3,:))/vpc;
    rec_latt(2,:) = twopi*cross(latt(3,:),latt(1,:))/vpc;
    rec_latt(3,:) = twopi*cross(latt(1,:),latt(2,:))/vpc;

    if(present(rec_latt_vol))then
        rec_latt_vol = (8.0*(Pi**3))/vpc
    endif

end subroutine get_rec_latt

function coords_cart_vec_in_new_basis(cart_vec, new_basis) result(coords)
! Returns the coordinates of the vector cart_vec = (x,y,z) in the basis "new_basis"
real*8, dimension(1:3) :: coords
real*8, dimension(1:3), intent(in) :: cart_vec
real*8, dimension(1:3,1:3), intent(in) :: new_basis
real*8, dimension(1:3,1:3) :: aux_matrix
integer :: i

    aux_matrix = inverse_of_3x3_matrix(transpose(new_basis))
    do i=1,3
        coords(i) = dot_product(aux_matrix(i,:), cart_vec)
    enddo

end function coords_cart_vec_in_new_basis

function same_vector(v1,v2,tolerance) result(rtn)
!! Copyright (C) 2013 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
! Checks if v1 and v2 are the same vector (within a given tolerance)
implicit none
logical :: rtn
real(kind=dp), dimension(1:3), intent(in) :: v1, v2
real(kind=dp), intent(in), optional :: tolerance
real(kind=dp) :: tol,dist 
integer :: i

    tol = default_tol_for_vec_equality
    if(present(tolerance))then
        tol = dabs(tolerance)
    endif

    ! Test 1: a sufficient condition for d_vec > tolerance
    rtn = .TRUE.
    do i=1,3
        if(dabs(v2(i)-v1(i)) > tol)then
            rtn = .FALSE.
            return
        endif
    enddo
    ! Test 2: if it passes by test 1
    rtn = .TRUE.
    dist = norm(v2(:)-v1(:))
    if(dist > tol)then
        rtn = .FALSE.
        return
    endif

end function same_vector


function vec_in_latt(vec, latt, tolerance) result(rtn)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
!! Checks if the vector "vec" belongs to the Bravais lattice spanned by the
!! vectors latt(i,:), i=1,2,3.
implicit none
logical :: rtn
real*8, dimension(1:3), intent(in) :: vec
real*8, dimension(1:3,1:3), intent(in) :: latt
real*8, intent(in), optional :: tolerance
real*8, dimension(1:3) :: reduced_vec,frac_coords,reduced_frac_coords,g
real*8 :: tol
integer :: i,ig1,ig2,ig3,ig1p,ig2p,ig3p,nb1max,nb2max,nb3max

tol = default_tol_for_vec_equality
if(present(tolerance))then
     tol = dabs(tolerance)
endif

rtn = .FALSE.
frac_coords(:) = coords_cart_vec_in_new_basis(cart_vec=vec, new_basis=latt)
reduced_frac_coords(:) = mod(frac_coords(:),1.0_dp)
reduced_vec(:) = 0.0_dp
do i=1,3
    reduced_vec(:) = reduced_vec(:) + reduced_frac_coords(i)*latt(i,:)
enddo

nb1max = 1
nb2max = 1
nb3max = 1
do ig3=0,2*nb3max
    ig3p=ig3
    if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1  ! Trick to make the indices ig3p vary from 0 to nb3max and then from -nb3max to -1
    do ig2=0,2*nb2max
        ig2p=ig2
        if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
        do ig1=0,2*nb1max
            ig1p=ig1
            if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
            g(:) = ig1p*latt(1,:) + ig2p*latt(2,:) + ig3p*latt(3,:)
            if(same_vector(reduced_vec,g,tol))then
                rtn = .TRUE.
                return
            endif
        enddo
    enddo
enddo
return

end function vec_in_latt


function point_is_in_bz(point,rec_latt,origin_point) result(rtn)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Checks if the point "point" belongs to the 1st Brillouin zone of a Bravais
!! lattice for which the reciprocal lettice vectors are "rec_latt(i,:)", i=1,2,3
!! The optional variable "origin_point" is to be used if you want the mesh to be
!! centered in other point than (0,0,0).
!! Not fully tested, but works for BandUP
implicit none
logical :: rtn
real(kind=dp), dimension(1:3), intent(in) :: point
real(kind=dp), dimension(1:3,1:3), intent(in) :: rec_latt
real(kind=dp), dimension(1:3), optional, intent(in) :: origin_point
real(kind=dp), dimension(1:3) :: origin, g
integer :: ig1, ig2, ig3

rtn = .TRUE.
if(present(origin_point))then
    origin(:) = origin_point(:)
else
    origin(:) = 0.0_dp
endif

do ig3=-1,1
    do ig2=-1,1
        do ig1=-1,1
            if(ig1==0 .and. ig2==0 .and. ig3==0) cycle
            g(:) = origin(:) + real(ig1,kind=dp)*rec_latt(1,:) + &
                               real(ig2,kind=dp)*rec_latt(2,:) + & 
                               real(ig3,kind=dp)*rec_latt(3,:)
            if(norm(point(:) - origin(:)) > norm(point(:) - g(:)))then
                rtn = .FALSE.
                return
            endif

        enddo
    enddo
enddo

end function point_is_in_bz


recursive subroutine reduce_point_to_bz(point,rec_latt,point_reduced_to_bz,max_index_lin_comb_RL_vecs)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Takes a point "point" and returns another point "point_reduced_to_bz" which
!! is equivalent to "point", but lies inside the 1st BZ of the Bravais lattice
!! with reciprocal lattice vectors "rec_latt(i,:)", i=1,2,3.
!! Not fully tested, but works for BandUP
implicit none
real(kind=dp), dimension(1:3), intent(in) :: point
real(kind=dp), dimension(1:3,1:3), intent(in) :: rec_latt
real(kind=dp), dimension(1:3), intent(out) :: point_reduced_to_bz ! Returned in cartesian coordinates!
integer, intent(in), optional :: max_index_lin_comb_RL_vecs
logical :: point_has_been_reduced
real(kind=dp), dimension(1:3) :: frac_coords_point, nearest_G, ref_G, trial_nearest_G
integer :: ig1, ig2, ig3, igmax

    igmax = 0
    if(present(max_index_lin_comb_RL_vecs))then
        igmax = max_index_lin_comb_RL_vecs
    endif

    frac_coords_point = coords_cart_vec_in_new_basis(cart_vec=point, new_basis=rec_latt)
    ! Finding the rec. latt. point "nearest_G" that is closest to the input point "point"
    ref_G = matmul(real(nint(frac_coords_point), kind=dp),rec_latt) ! In cartesian coords
    nearest_G = ref_G
    do ig3=-igmax,igmax
        do ig2=-igmax,igmax
            do ig1=-igmax,igmax
                trial_nearest_G = ref_G + ig1*rec_latt(1,:) + ig2*rec_latt(2,:) + ig3*rec_latt(3,:)
                if(norm(point - nearest_G) > norm(point - trial_nearest_G))then
                    nearest_G = trial_nearest_G
                endif
            enddo 
        enddo 
    enddo
    ! Reducing the point, now that we have the nearest rec. latt. point
    point_reduced_to_bz(:) = point(:) - nearest_G(:)

    ! Double-checking if the point has really been reduced 
    ! The routine will raise an error and stop otherwise
    point_has_been_reduced = point_is_in_bz(point_reduced_to_bz,rec_latt)
    if(.not. point_has_been_reduced)then
        if(igmax>=10)then !  Recursion limit
            write(*,'(A)')'ERROR (reduce_point_to_bz): Failed to reduce k-point to the 1st Brillouin zone.'
            write(*,'(A,I0,A)')'                            A total of ',igmax+1,' layer(s) of neighboring rec. latt. points were checked for nearest neighbors.'
            write(*,'(A)')'                            The nearest rec. latt. point was not found.'
            write(*,'(3(A,f0.5),A)')'                            * Cartesian coordinates of the k-point:  (', &
                                    point(1),', ',point(2),', ',point(2),').'
            write(*,'(3(A,f0.5),A)')'                            * Fractional coordinates of the k-point: (', &
                                    frac_coords_point(1),', ',frac_coords_point(2),', ',frac_coords_point(2),').'
            write(*,'(A)')'Stopping now.'
            stop
        else
            ! If the previous scan fails, go for an extra layer
            call reduce_point_to_bz(point,rec_latt,point_reduced_to_bz,igmax+1) 
        endif
    endif

end subroutine reduce_point_to_bz


subroutine get_prim_cell(crystal_pc, crystal, &
                         cystal_is_prim_cell, symmetryze_and_refine_cell, symprec)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! If the crystal is already a primitive cell, then a copy of it is
!! returned in the crystal_pc variable
!! Not fully tested, but works for BandUP
implicit none
type(crystal_3D), intent(out), allocatable :: crystal_pc
type(crystal_3D), intent(in) :: crystal
logical, intent(out), optional :: cystal_is_prim_cell
logical, intent(in), optional :: symmetryze_and_refine_cell
real(kind=dp), intent(in), optional :: symprec
type(crystal_3D), allocatable :: aux_crystal
real(kind=dp), dimension(:,:), allocatable :: transpose_frac_coords_SC, transpose_coords_pc, &
                                              aux_transpose_frac_coords
integer :: num_atoms_pc, num_atoms_SC,iatom, iatom2, alloc_stat
integer, dimension(:), allocatable :: aux_integer_types_basis_atoms
real(kind=dp) :: symm_prec
character(len=3), dimension(:), allocatable :: atom_symbols
logical :: symmetryze_and_ref_cell

    symmetryze_and_ref_cell = .FALSE.
    if(present(symmetryze_and_refine_cell))then
        symmetryze_and_ref_cell = symmetryze_and_refine_cell
    endif
 
    symm_prec=default_symprec
    if(present(symprec))then 
        symm_prec=dabs(symprec)
    endif

    allocate(aux_crystal)
    aux_crystal = crystal

    if(symmetryze_and_ref_cell)then
        deallocate(aux_integer_types_basis_atoms, stat=alloc_stat)
        aux_integer_types_basis_atoms = 0 
        deallocate(aux_transpose_frac_coords, stat=alloc_stat)
        ! spglib requires space for 4*n_atoms to be allocated here, because
        ! 4 times more atoms are generated for face-centered crystals.
        allocate(aux_integer_types_basis_atoms(1:4*size(crystal%integer_types_basis_atoms)))
        allocate(aux_transpose_frac_coords(1:size(crystal%fractional_coords_basis_atoms,dim=2), &
                                          1:4*size(crystal%fractional_coords_basis_atoms,dim=1)))
        aux_transpose_frac_coords(:,1:size(crystal%fractional_coords_basis_atoms,dim=1)) = &
            transpose(crystal%fractional_coords_basis_atoms(:,:))
        aux_integer_types_basis_atoms(1:size(crystal%integer_types_basis_atoms)) = &
            crystal%integer_types_basis_atoms(:)
        num_atoms_SC = size(crystal%fractional_coords_basis_atoms,dim=1)
        num_atoms_SC = spg_refine_cell(aux_crystal%latt_vecs, &
                                       aux_transpose_frac_coords, &
                                       aux_integer_types_basis_atoms, &
                                       num_atoms_SC, symm_prec)
        deallocate(aux_crystal%fractional_coords_basis_atoms, stat=alloc_stat)
        deallocate(aux_crystal%integer_types_basis_atoms, stat=alloc_stat)
        allocate(aux_crystal%fractional_coords_basis_atoms(1:num_atoms_SC,1:3))
        allocate(aux_crystal%integer_types_basis_atoms(1:num_atoms_SC))
        aux_crystal%fractional_coords_basis_atoms(:,:) = &
            transpose(aux_transpose_frac_coords(:,1:num_atoms_SC))
        aux_crystal%integer_types_basis_atoms = aux_integer_types_basis_atoms(1:num_atoms_SC) 
        deallocate(aux_transpose_frac_coords, stat=alloc_stat)
        deallocate(aux_integer_types_basis_atoms, stat=alloc_stat)
    endif

    allocate(transpose_frac_coords_SC(1:size(aux_crystal%fractional_coords_basis_atoms,dim=2), &
                                      1:size(aux_crystal%fractional_coords_basis_atoms,dim=1)))
    transpose_frac_coords_SC = transpose(aux_crystal%fractional_coords_basis_atoms)
    num_atoms_pc = spg_find_primitive(aux_crystal%latt_vecs, &
                                      transpose_frac_coords_SC, &
                                      aux_crystal%integer_types_basis_atoms, &
                                      size(aux_crystal%fractional_coords_basis_atoms, dim=1), &
                                      symm_prec)
    if(num_atoms_pc>0)then 
       ! This means that the input crystal celll is a perfect supercell
       ! corresponding to a primitive unit cell with num_atoms_pc atoms
        if(present(cystal_is_prim_cell)) cystal_is_prim_cell = .FALSE.
        allocate(atom_symbols(1:num_atoms_pc))
        allocate(transpose_coords_pc(1:3,1:num_atoms_pc))
        do iatom=1,num_atoms_pc
            do iatom2=1,size(crystal%atomic_symbols_basis_atoms)
                if(aux_crystal%integer_types_basis_atoms(iatom) == &
                   crystal%integer_types_basis_atoms(iatom2))then
                    atom_symbols(iatom) = crystal%atomic_symbols_basis_atoms(iatom2)
                    exit
                endif
            enddo
            ! Making sure the atoms are inside the simulation box
            ! modulo(a,p) returns a - floor(a / p) * p
            transpose_frac_coords_SC(:,iatom) = modulo(transpose_frac_coords_SC(:,iatom), &
                                                       real(1,kind=dp))
            ! Changing to cartesian coords
            transpose_coords_pc(:,iatom) = matmul(transpose_frac_coords_SC(:,iatom), & 
                                                  aux_crystal%latt_vecs)
        enddo
        call create_crystal(crystal=crystal_pc, description=crystal%description, &
                            latt_vecs=aux_crystal%latt_vecs, & 
                            coords_basis_atoms=transpose(transpose_coords_pc), &
                            atomic_symbols_basis_atoms=atom_symbols)
    else
        ! If the crystal is already a primitive cell, then a copy of it is
        ! returned in the crystal_pc variable
        if(present(cystal_is_prim_cell)) cystal_is_prim_cell = .TRUE.
        allocate(crystal_pc)
        crystal_pc = crystal
    endif

end subroutine get_prim_cell


subroutine get_symm(nsym, symops, schoenflies, international_symb, symprec, crystal, &
                    cystal_is_prim_cell, try_to_reduce_to_a_pc, crystal_reduced_prim_cell)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
integer, intent(out), optional :: nsym
type(symmetry_operation), dimension(:), allocatable, optional, intent(out) :: symops
character(len=10), intent(out), optional :: schoenflies
character(len=11), intent(out), optional :: international_symb
real(kind=dp), intent(in), optional :: symprec
type(crystal_3D), intent(in) :: crystal
logical, intent(out), optional :: cystal_is_prim_cell
logical, intent(in), optional :: try_to_reduce_to_a_pc
type(crystal_3D), intent(out), allocatable, optional :: crystal_reduced_prim_cell
! Local variables
type(crystal_3D), allocatable :: work_crystal, aux_crystal
type(SpglibDataset) :: symm_dataset
real(kind=dp) :: symm_prec
real(kind=dp), dimension(1:3,1:3) :: M, M_inv
integer :: isym, alloc_stat, space_group_num, n_operations
character(len=10) :: schoenflies_notation
logical :: work_with_pc, already_pc

    symm_prec = default_symprec
    if(present(symprec))then 
        symm_prec=dabs(symprec)
    endif
    work_with_pc = .FALSE.
    if(present(try_to_reduce_to_a_pc))then
        work_with_pc = try_to_reduce_to_a_pc
    endif

    already_pc = .TRUE.
    if(work_with_pc)then
        call get_prim_cell(crystal_pc=work_crystal,crystal=crystal, &
                           cystal_is_prim_cell=already_pc, symprec=symm_prec)
        if(present(crystal_reduced_prim_cell))then
            deallocate(crystal_reduced_prim_cell, stat=alloc_stat)
            allocate(crystal_reduced_prim_cell)
            crystal_reduced_prim_cell = work_crystal
        endif
    else
        ! The use of aux_crystal is just to know whether the crystal cell is a
        ! pc or not. I'll be working with the original input crysral anyway.
        call get_prim_cell(crystal_pc=aux_crystal,crystal=crystal, &
                           cystal_is_prim_cell=already_pc, symprec=symm_prec)
        if(present(crystal_reduced_prim_cell))then
            deallocate(crystal_reduced_prim_cell, stat=alloc_stat)
            allocate(crystal_reduced_prim_cell)
            crystal_reduced_prim_cell = aux_crystal
        endif
        deallocate(aux_crystal, stat=alloc_stat)

        allocate(work_crystal)
        work_crystal = crystal
    endif
    if(present(cystal_is_prim_cell)) cystal_is_prim_cell = already_pc

    symm_dataset = spg_get_dataset(work_crystal%latt_vecs, &
                                   transpose(work_crystal%fractional_coords_basis_atoms), &
                                   work_crystal%integer_types_basis_atoms, &
                                   size(work_crystal%fractional_coords_basis_atoms, dim=1), &
                                   symm_prec)

    n_operations = symm_dataset % n_operations
    if(present(nsym))then
        nsym = n_operations
    endif

    if(present(symops))then
        M = work_crystal%latt_vecs ! Matrix that changes from the lattice basis to the canonical basis in 3D
        M_inv = inverse_of_3x3_matrix(M) 
        deallocate(symops,stat=alloc_stat)
        allocate(symops(1:n_operations))
        do isym=1, n_operations
            symops(isym) % translation_fractional_coords(:) = symm_dataset % translations(:,isym)
            symops(isym) % rotation_fractional_coords(:,:) = symm_dataset % rotations(:,:,isym)
            symops(isym) % translation_cartesian_coords(:) = matmul(symops(isym) % translation_fractional_coords(:), M)
            symops(isym) % rotation_cartesian_coords(:,:) = matmul(M_inv,matmul(symops(isym) % rotation_fractional_coords(:,:),M))
        enddo
    endif

    if(present(schoenflies))then
        schoenflies_notation = '     '
        space_group_num = spg_get_schoenflies(schoenflies_notation, work_crystal%latt_vecs, &
                                              transpose(work_crystal%fractional_coords_basis_atoms), &
                                              work_crystal%integer_types_basis_atoms, &
                                              size(work_crystal%fractional_coords_basis_atoms, dim=1), &
                                              symm_prec)
        schoenflies = schoenflies_notation(:len(schoenflies_notation)-1)
    endif

    if(present(international_symb))then
        international_symb = symm_dataset % international_symbol(:len(symm_dataset % international_symbol)-1)
    endif

end subroutine get_symm


function pt_eqv_by_point_group_symop(point,symops,isym,fractional_coords,invert_symop) result(eqv_point)
!! Copyright (C) 2013 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
real(kind=dp), dimension(1:3) :: eqv_point
real(kind=dp), dimension(1:3), intent(in) :: point
type(symmetry_operation), dimension(:), intent(in) :: symops
logical, intent(in), optional :: fractional_coords,invert_symop
integer, intent(in) :: isym
logical :: frac_coords, invert

        frac_coords = .FALSE.
        if(present(fractional_coords))then
            frac_coords = fractional_coords
        endif
        invert=.FALSE.
        if(present(invert_symop))then
            invert = invert_symop
        endif

        if(frac_coords)then
            if(.not. invert)then
                eqv_point(:) = matmul(point(:),symops(isym) % rotation_fractional_coords(:,:))
            else
                eqv_point(:) = matmul(point(:),inverse_of_3x3_matrix(1.0_dp*symops(isym) % rotation_fractional_coords(:,:)))
            endif
        else
            if(.not. invert)then
                eqv_point(:) = matmul(point(:),symops(isym) % rotation_cartesian_coords(:,:))
            else
                eqv_point(:) = matmul(point(:),inverse_of_3x3_matrix(1.0_dp*symops(isym) % rotation_cartesian_coords(:,:)))
            endif
        endif

end function pt_eqv_by_point_group_symop

subroutine get_star(star_of_pt, list_of_all_generated_points, points, crystal, symm_ops, &
                    tol_for_vec_equality, symprec, reduce_to_bz, try_to_reduce_to_a_pc)
!! Copyright (C) 2013 Paulo V. C. Medeiros
!! The input list of points shall be in cartesian coordinates.
!! The output list of points will be also given in cartesian coordinates.
!! Not fully tested, but works for BandUP
implicit none
type :: symop_properties
    real(kind=dp), dimension(1:3) :: coord
    logical :: is_eqv_to_a_previous_symop
    logical, dimension(:), allocatable :: eqv_to_symop
end type symop_properties
type :: degenerate_star
    type(symop_properties), dimension(:), allocatable :: symop
end type degenerate_star


type(star), dimension(:), allocatable, intent(out), optional :: star_of_pt 
type(vec3d), dimension(:), allocatable, intent(out), optional :: list_of_all_generated_points
type(vec3d), dimension(:), intent(in) :: points
type(crystal_3D), intent(in) :: crystal
type(symmetry_operation), dimension(:), allocatable, optional, intent(inout) :: symm_ops
real(kind=dp), intent(in), optional :: tol_for_vec_equality,symprec
logical, intent(in), optional :: reduce_to_bz, try_to_reduce_to_a_pc

type(degenerate_star), dimension(:), allocatable :: possibly_degenerate_pts_eqv_to_pt
type(star), dimension(:), allocatable :: reduced_points_eqv_to_pt 
real(kind=dp), dimension(1:3) :: point, point_reduced_to_bz
real(kind=dp), dimension(1:3,1:3) :: rec_lattice
integer :: size_points,ipt,ipt2,nsym,isym,isym2,iupt,neqv,alloc_stat,n_generated_points,ieqv_pt
real(kind=dp) :: tol,sym_prec
logical :: reduce2bz, reduce2pc
type(symmetry_operation), dimension(:), allocatable :: symops

    tol = default_tol_for_vec_equality
    if(present(tol_for_vec_equality))then
        tol = dabs(tol_for_vec_equality)
    endif
    sym_prec = default_symprec
    if(present(symprec))then
        sym_prec = dabs(symprec)
    endif
    reduce2bz = .FALSE.
    if(present(reduce_to_bz))then
        reduce2bz = reduce_to_bz
    endif
    reduce2pc = .FALSE.
    if(present(try_to_reduce_to_a_pc))then
        reduce2pc = try_to_reduce_to_a_pc
    endif

    deallocate(symops, stat=alloc_stat)
    if(present(symm_ops))then
        if(size(symm_ops)>0)then
            allocate(symops(1:size(symm_ops)), stat=alloc_stat)
            symops = symm_ops
        else
            call get_symm(symops=symops, symprec=sym_prec, crystal=crystal,try_to_reduce_to_a_pc=reduce2pc)
            deallocate(symm_ops, stat=alloc_stat)
            allocate(symm_ops(1:size(symops)), stat=alloc_stat)
            symm_ops = symops
        endif
    else
        call get_symm(symops=symops, symprec=sym_prec,crystal=crystal,try_to_reduce_to_a_pc=reduce2pc)
    endif
    nsym = size(symops)

    size_points = size(points)
    rec_lattice = crystal%rec_latt_vecs

    deallocate(possibly_degenerate_pts_eqv_to_pt,stat=alloc_stat)
    allocate(possibly_degenerate_pts_eqv_to_pt(1:size_points))
    do ipt=1, size_points
        deallocate(possibly_degenerate_pts_eqv_to_pt(ipt)%symop,stat=alloc_stat)
        allocate(possibly_degenerate_pts_eqv_to_pt(ipt)%symop(1:nsym))
        point(:) = points(ipt)%coord(:) ! In cartesian coords.
        if(reduce2bz)then
            call reduce_point_to_bz(point, rec_lattice, point_reduced_to_bz) ! Mind that the subroutine returns the vector in cartesian coords.
            point = point_reduced_to_bz
        endif
        do isym=1, nsym
            deallocate(possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%eqv_to_symop,stat=alloc_stat)
            allocate(possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%eqv_to_symop(1:nsym))
            ! The symmetry ops have to be performed using cartesian coordinates,
            ! because they've been obtained from the real-space lattice, but
            ! they'll be applied to the points in the reciprocal space
            possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%coord(:) = &
                pt_eqv_by_point_group_symop(point=point, symops=symops, isym=isym, &
                                            fractional_coords=.FALSE.)
        enddo

        do isym=1,nsym
            possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%is_eqv_to_a_previous_symop = .FALSE.
            do isym2=1,nsym
                possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%eqv_to_symop(isym2) = (isym==isym2)
            enddo
        enddo

        do isym=1,nsym-1
            do isym2=isym+1,nsym
                if(same_vector(possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym2)%coord(:), &
                               possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%coord(:), &
                               tol))then
                    possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%eqv_to_symop(isym2) = .TRUE.
                    possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym2)%eqv_to_symop(isym) = .TRUE.
                    possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym2)%is_eqv_to_a_previous_symop = .TRUE.
                endif
            enddo
        enddo
    enddo
 
    deallocate(reduced_points_eqv_to_pt,stat=alloc_stat)
    allocate(reduced_points_eqv_to_pt(1:size_points))
    do ipt=1,size_points
        neqv = count(.not. possibly_degenerate_pts_eqv_to_pt(ipt)%symop(:)%is_eqv_to_a_previous_symop)
        deallocate(reduced_points_eqv_to_pt(ipt)%eqv_pt,stat=alloc_stat)
        allocate(reduced_points_eqv_to_pt(ipt)%eqv_pt(1:neqv))
        reduced_points_eqv_to_pt(ipt)%neqv = neqv
        iupt = 0
        do isym=1, nsym
            if(.not. possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%is_eqv_to_a_previous_symop)then
                iupt = iupt + 1 !! iupt = index of "unique" pt
                reduced_points_eqv_to_pt(ipt)%eqv_pt(iupt)%coord(:) = possibly_degenerate_pts_eqv_to_pt(ipt)%symop(isym)%coord(:)
                reduced_points_eqv_to_pt(ipt)%eqv_pt(iupt)%symop = isym
            endif
        enddo
    enddo

    if(present(star_of_pt))then
        deallocate(star_of_pt,stat=alloc_stat)
        allocate(star_of_pt(1:size(reduced_points_eqv_to_pt)))
        star_of_pt = reduced_points_eqv_to_pt
    endif

    if(present(list_of_all_generated_points))then
        n_generated_points = sum(reduced_points_eqv_to_pt(:)%neqv)
        allocate(list_of_all_generated_points(1:n_generated_points))
        ipt2=0
        do ipt=1,size(reduced_points_eqv_to_pt(:))
            do ieqv_pt=1,reduced_points_eqv_to_pt(ipt)%neqv
                ipt2 = ipt2 + 1
                list_of_all_generated_points(ipt2)%coord(:) = reduced_points_eqv_to_pt(ipt)%eqv_pt(ieqv_pt)%coord(:)
            enddo
        enddo
    endif

end subroutine get_star

subroutine get_irr_kpts(n_irr_kpts,irr_kpts_list,irr_kpts_list_frac_coords, &
                        kpts_list,crystal,reduce_to_bz,symprec)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
integer, intent(out) :: n_irr_kpts
type(crystal_3D), intent(in) :: crystal
type(vec3d), dimension(:), intent(in) :: kpts_list
type(vec3d), dimension(:), allocatable, intent(out) :: irr_kpts_list, irr_kpts_list_frac_coords
logical, intent(in), optional :: reduce_to_bz
real(kind=dp), optional :: symprec

real(kind=dp), dimension(1:3) :: point, previous_point, point_reduced_to_SCBZ
real(kind=dp), dimension(1:3,1:3) :: rec_latt
type(symmetry_operation), dimension(:), allocatable :: symops
integer :: ikpt,ikpt2,ieqv_kpt,i_irr_kpt
logical :: reduce2bz
logical, dimension(:), allocatable :: kpt_is_eqv_to_previous_kpt 
real(kind=dp) :: sprec
type(star), dimension(:), allocatable :: SKPTS_eqv_to_SKPT
type(vec3d), dimension(:), allocatable :: work_kpts_list

    reduce2bz = .FALSE.
    if(present(reduce_to_bz))then
        reduce2bz = reduce_to_bz
    endif

    sprec = default_symprec
    if(present(symprec))then
        sprec = symprec
    endif

    rec_latt = crystal%rec_latt_vecs

    allocate(work_kpts_list(1:size(kpts_list)))
    if(reduce2bz)then
        do ikpt=1, size(kpts_list)
            call reduce_point_to_bz(kpts_list(ikpt)%coord(:), rec_latt, point_reduced_to_SCBZ)  ! In cartesian coords.
            work_kpts_list(ikpt)%coord(:) = point_reduced_to_SCBZ(:)
        enddo
    else
        work_kpts_list = kpts_list
    endif

    call get_star(star_of_pt=SKPTS_eqv_to_SKPT, points=work_kpts_list, &
                  crystal=crystal, symm_ops=symops, symprec=sprec, reduce_to_bz=reduce2bz)
   

    allocate(kpt_is_eqv_to_previous_kpt(1:size(work_kpts_list))) 
    kpt_is_eqv_to_previous_kpt(:) = .FALSE.
    do ikpt=size(work_kpts_list),1,-1
        point(:) = work_kpts_list(ikpt)%coord(:)
        if(ikpt>1)then
            do ikpt2=ikpt-1,1,-1
                do ieqv_kpt=1,size(SKPTS_eqv_to_SKPT(ikpt2)%eqv_pt(:))
                    previous_point(:) = SKPTS_eqv_to_SKPT(ikpt2)%eqv_pt(ieqv_kpt)%coord(:)
                    if(same_vector(point,previous_point))then
                        kpt_is_eqv_to_previous_kpt(ikpt) = .TRUE.
                        exit
                    endif
                enddo
            enddo 
        endif
    enddo
    n_irr_kpts = count(.not. kpt_is_eqv_to_previous_kpt(:))

   ! Returning the lists
    allocate(irr_kpts_list(1:n_irr_kpts), irr_kpts_list_frac_coords(1:n_irr_kpts))
    i_irr_kpt = 0
    do ikpt=1, size(work_kpts_list)
        if(.not. kpt_is_eqv_to_previous_kpt(ikpt))then
            i_irr_kpt = i_irr_kpt + 1
            irr_kpts_list(i_irr_kpt)%coord(:) = work_kpts_list(ikpt)%coord(:)
            irr_kpts_list_frac_coords(i_irr_kpt)%coord(:) = &
                coords_cart_vec_in_new_basis(cart_vec=work_kpts_list(ikpt)%coord(:), &
                                             new_basis=rec_latt)
        endif
    enddo

end subroutine get_irr_kpts


function trace(square_matrix) result(rtn)
implicit none
real(kind=dp) :: rtn
real(kind=dp), dimension(:,:) :: square_matrix
integer :: i, dimens

    dimens=size(square_matrix,dim=1)
    rtn = 0.0_dp
    do i=1,dimens
        rtn = rtn + square_matrix(i,i)
    enddo

end function trace


function frobenius_inner_product(A,B) result(rtn)
implicit none
real(kind=dp), intent(in), dimension(:,:) :: A, B
real(kind=dp) :: rtn

   rtn = trace(matmul(A,transpose(B))) 

end function frobenius_inner_product


function frobenius_norm(A) result(rtn)
implicit none
real(kind=dp), intent(in), dimension(:,:) :: A
real(kind=dp) :: rtn

    rtn = dsqrt(frobenius_inner_product(A,A))

end function frobenius_norm


subroutine get_compl_pcbz_direcs(compl_pcdirs,ncompl,dirs_eqv_in_pcbz,dirs_eqv_in_SCBZ,symops_SC,symprec)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
type(eqv_bz_directions), intent(out) :: compl_pcdirs
integer, intent(out) :: ncompl
type(eqv_bz_directions), intent(in) :: dirs_eqv_in_pcbz,dirs_eqv_in_SCBZ
type(symmetry_operation), dimension(:), intent(in) :: symops_SC
real(kind=dp), intent(in), optional :: symprec
real(kind=dp) :: sprec
integer :: neqv_pc,neqv_SC,ieqv_pc,ieqv_SC,icompl
real(kind=dp), dimension(1:3) :: kstart_pc,kend_pc,kstart_SC,kend_SC
logical, dimension(:), allocatable :: direc_also_eqv_in_SC


    ncompl = 0
    sprec = default_symprec
    if(present(symprec))then
        sprec=symprec
    endif
    neqv_pc = size(dirs_eqv_in_pcbz%eqv_dir(:))
    neqv_SC = size(dirs_eqv_in_SCBZ%eqv_dir(:))

    allocate(direc_also_eqv_in_SC(1:neqv_pc))
    direc_also_eqv_in_SC(:) = .FALSE.
    do ieqv_pc=1,neqv_pc
        kstart_pc(:) = dirs_eqv_in_pcbz%eqv_dir(ieqv_pc)%kstart(:)
        kend_pc(:) = dirs_eqv_in_pcbz%eqv_dir(ieqv_pc)%kend(:)
        do ieqv_SC=1,neqv_SC
            kstart_SC = dirs_eqv_in_SCBZ%eqv_dir(ieqv_SC)%kstart(:)
            kend_SC = dirs_eqv_in_SCBZ%eqv_dir(ieqv_SC)%kend(:)
            if(direcs_are_eqv(kstart_pc,kend_pc,kstart_SC,kend_SC,symops_SC,sprec))then
                direc_also_eqv_in_SC(ieqv_pc) = .TRUE.
                exit
            endif
        enddo
    enddo
    ncompl = count(.not. direc_also_eqv_in_SC)
    if(ncompl /= 0)then
        allocate(compl_pcdirs%eqv_dir(1:ncompl))
        icompl = 0
        do ieqv_pc=1,neqv_pc
            if(.not. direc_also_eqv_in_SC(ieqv_pc))then
                icompl = icompl + 1
                compl_pcdirs%eqv_dir(icompl) = dirs_eqv_in_pcbz%eqv_dir(ieqv_pc)
            endif
        enddo
    endif

end subroutine get_compl_pcbz_direcs

subroutine get_irr_bz_directions(irr_dirs,nirr,list_of_dirs,symm_ops,crystal,symprec)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
type(irr_bz_directions), intent(out) :: irr_dirs
integer, intent(out) :: nirr
type(eqv_bz_directions), intent(in) :: list_of_dirs
type(symmetry_operation), dimension(:), optional, intent(in) :: symm_ops
type(crystal_3D), intent(in) :: crystal
real(kind=dp), intent(in), optional :: symprec
real(kind=dp) :: sprec
logical, dimension(:,:), allocatable ::  dirs_are_eqv
logical, dimension(:), allocatable :: dir_eqv_to_a_previous_one
integer :: ndirs,idir,idir2
real(kind=dp), dimension(1:3) :: current_kstart,current_kend,next_kstart,next_kend
type(symmetry_operation), dimension(:), allocatable :: symops

    nirr = 0
    sprec = default_symprec
    if(present(symprec))then
        sprec=symprec
    endif
    if(present(symm_ops))then
        allocate(symops(1:size(symm_ops)))
        symops=symm_ops
    else
        call get_symm(symops=symops, symprec=sprec, crystal=crystal)
    endif

    ndirs = size(list_of_dirs%eqv_dir(:))
    allocate(dirs_are_eqv(1:ndirs,1:ndirs))
    dirs_are_eqv(:,:) = .FALSE.
    allocate(dir_eqv_to_a_previous_one(1:ndirs))
    dir_eqv_to_a_previous_one(:) = .FALSE.
    do idir=1,ndirs-1
        dirs_are_eqv(idir,idir) = .TRUE.
        if(dir_eqv_to_a_previous_one(idir)) cycle
        current_kstart = list_of_dirs%eqv_dir(idir)%kstart
        current_kend = list_of_dirs%eqv_dir(idir)%kend
        do idir2=idir+1,ndirs
            if(dir_eqv_to_a_previous_one(idir2)) cycle
            next_kstart = list_of_dirs%eqv_dir(idir2)%kstart
            next_kend = list_of_dirs%eqv_dir(idir2)%kend
            if(direcs_are_eqv(current_kstart,current_kend,next_kstart,next_kend,symops))then
                dir_eqv_to_a_previous_one(idir2) = .TRUE.
                dirs_are_eqv(idir,idir2) = .TRUE.
                dirs_are_eqv(idir2,idir) = dirs_are_eqv(idir,idir2)
            endif
        enddo
    enddo
    dirs_are_eqv(ndirs,ndirs) = .TRUE.
    nirr = count(.not. dir_eqv_to_a_previous_one)
    if(nirr > 0)then
        allocate(irr_dirs%irr_dir(1:nirr))
        idir2 = 0
        do idir=1,ndirs
            if(.not. dir_eqv_to_a_previous_one(idir))then
                idir2 = idir2 + 1
                irr_dirs%irr_dir(idir2) = list_of_dirs%eqv_dir(idir)
                irr_dirs%irr_dir(idir2)%neqv = count(dirs_are_eqv(idir,:))
            endif
        enddo
    endif

end subroutine get_irr_bz_directions


subroutine get_eqv_bz_dirs(eqv_dirs,kstart,kend,symm_ops,crystal,symprec) 
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
type(eqv_bz_directions), intent(out) :: eqv_dirs
real(kind=dp), dimension(1:3), intent(in) :: kstart,kend
type(symmetry_operation), dimension(:), optional, intent(in) :: symm_ops
type(crystal_3D), intent(in) :: crystal
real(kind=dp), intent(in), optional :: symprec
real(kind=dp) :: sprec
logical, dimension(:,:), allocatable :: same_dirs,pairs_start_end_lead_to_eqv_dir
logical, dimension(:), allocatable :: repeated_dir
integer :: isym,nsym,n_eqv_dirs_full_star,ieqv_dir,ieqv_dir2,n_eqv_dirs,i_eqv_dirs_full_star,istart,iend
real(kind=dp), dimension(1:3) :: current_kstart,current_kend,next_kstart,next_kend
type(eqv_bz_directions) :: full_star_of_eqv_dirs
type(symmetry_operation), dimension(:), allocatable :: symops
type(vec3d), dimension(:), allocatable :: full_star_of_eqv_starts, full_star_of_eqv_ends

    sprec = default_symprec
    if(present(symprec))then
        sprec=symprec
    endif
    if(present(symm_ops))then
        allocate(symops(1:size(symm_ops)))
        symops=symm_ops
        nsym = size(symops)
    else
        call get_symm(nsym=nsym, symops=symops, symprec=sprec, crystal=crystal) 
    endif

    allocate(full_star_of_eqv_starts(1:nsym))
    allocate(full_star_of_eqv_ends(1:nsym))
    do isym=1,nsym
        full_star_of_eqv_starts(isym)%coord(:) = pt_eqv_by_point_group_symop(point=kstart,symops=symops,isym=isym)
        full_star_of_eqv_ends(isym)%coord(:) = pt_eqv_by_point_group_symop(point=kend,symops=symops,isym=isym)
    enddo
    
    allocate(pairs_start_end_lead_to_eqv_dir(1:nsym,1:nsym))
    pairs_start_end_lead_to_eqv_dir(:,:) = .FALSE.
    do istart=1,nsym
        current_kstart = full_star_of_eqv_starts(istart)%coord(:)
        do iend=1,nsym
            current_kend = full_star_of_eqv_ends(iend)%coord(:)
            pairs_start_end_lead_to_eqv_dir(istart,iend) = direcs_are_eqv(kstart,kend,current_kstart,current_kend,symops,sprec) 
        enddo
    enddo
    n_eqv_dirs_full_star = count(pairs_start_end_lead_to_eqv_dir == .TRUE.)

    allocate(full_star_of_eqv_dirs%eqv_dir(1:n_eqv_dirs_full_star))
    ieqv_dir=0
    do istart=1,nsym
        do iend=1,nsym
            if(pairs_start_end_lead_to_eqv_dir(istart,iend))then
                ieqv_dir = ieqv_dir + 1
                full_star_of_eqv_dirs%eqv_dir(ieqv_dir)%kstart(:) = full_star_of_eqv_starts(istart)%coord(:)
                full_star_of_eqv_dirs%eqv_dir(ieqv_dir)%kend(:) = full_star_of_eqv_ends(iend)%coord(:)
            endif
        enddo
    enddo

    allocate(repeated_dir(1:n_eqv_dirs_full_star))
    allocate(same_dirs(1:n_eqv_dirs_full_star,1:n_eqv_dirs_full_star))
    repeated_dir(:) = .FALSE.
    same_dirs(:,:) = .FALSE.
    do ieqv_dir=1,n_eqv_dirs_full_star-1
        if(repeated_dir(ieqv_dir)) cycle
        current_kstart = full_star_of_eqv_dirs%eqv_dir(ieqv_dir)%kstart(:)
        current_kend = full_star_of_eqv_dirs%eqv_dir(ieqv_dir)%kend(:)
        do ieqv_dir2=ieqv_dir+1,n_eqv_dirs_full_star
            if(repeated_dir(ieqv_dir2)) cycle
            next_kstart = full_star_of_eqv_dirs%eqv_dir(ieqv_dir2)%kstart(:)
            next_kend = full_star_of_eqv_dirs%eqv_dir(ieqv_dir2)%kend(:)
            if(same_vector(current_kstart,next_kstart,sprec) .and. same_vector(current_kend,next_kend,sprec))then
                repeated_dir(ieqv_dir2) = .TRUE.
                same_dirs(ieqv_dir,ieqv_dir2) = .TRUE.
                same_dirs(ieqv_dir2,ieqv_dir) = same_dirs(ieqv_dir,ieqv_dir2)
            endif
        enddo        
    enddo
    n_eqv_dirs = count(.not. repeated_dir)
    eqv_dirs%neqv = n_eqv_dirs
    allocate(eqv_dirs%eqv_dir(1:n_eqv_dirs))
    ieqv_dir = 0
    do i_eqv_dirs_full_star=1,n_eqv_dirs_full_star
        if(.not. repeated_dir(i_eqv_dirs_full_star))then
            ieqv_dir = ieqv_dir + 1
            eqv_dirs%eqv_dir(ieqv_dir) = full_star_of_eqv_dirs%eqv_dir(i_eqv_dirs_full_star)
        endif
    enddo
    
end subroutine get_eqv_bz_dirs

function points_are_eqv(v1,v2,symops,symprec) result(rtn)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
logical :: rtn
real(kind=dp), dimension(1:3), intent(in) :: v1, v2
type(symmetry_operation), dimension(:), intent(in) :: symops
real(kind=dp), intent(in), optional :: symprec
integer :: nsym,isym
real(kind=dp), dimension(1:3) :: eqv_v1
real(kind=dp) :: sprec

    rtn = .FALSE.

    sprec = default_symprec
    if(present(symprec))then
        sprec = symprec
    endif

    nsym = size(symops)
    do isym=1,nsym
        eqv_v1 = pt_eqv_by_point_group_symop(point=v1,symops=symops,isym=isym)
        if(same_vector(eqv_v1,v2,sprec))then
            rtn = .TRUE.
            return
        endif
    enddo

end function points_are_eqv

function direcs_are_eqv(start1,end1,start2,end2,symops,symprec) result(rtn)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
logical :: rtn
real(kind=dp), dimension(1:3), intent(in) :: start1,end1,start2,end2
type(symmetry_operation), dimension(:), intent(in) :: symops
real(kind=dp), intent(in), optional :: symprec
real(kind=dp) :: sprec

    rtn = .FALSE.

    sprec = default_symprec
    if(present(symprec))then
        sprec = symprec
    endif

    if(dabs(norm(end1 - start1) - norm(end2 - start2)) <= sprec)then
        rtn = (points_are_eqv(start1,start2,symops) .and. points_are_eqv(end1,end2,symops))
    endif


end function direcs_are_eqv



subroutine get_all_irr_dirs_req_for_symmavgd_EBS(dirs_req_for_symmavgd_EBS_along_pcbz_dir, &
                                                 n_dirs_for_EBS_along_pcbz_dir, &
                                                 neqv_dirs_pcbz, neqv_dirs_SCBZ, &
                                                 ncompl_dirs, n_irr_compl_dirs, &
                                                 crystal_pc_reduced_to_prim_cell, &
                                                 pc_is_prim_cell, &
                                                 crystal_SC_reduced_to_prim_cell, &
                                                 SC_is_prim_cell, &
                                                 crystal_pc,crystal_SC, &
                                                 k_starts,k_ends,verbose, &
                                                 symprec)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
type(irr_bz_directions), dimension(:), allocatable, intent(out) :: dirs_req_for_symmavgd_EBS_along_pcbz_dir
integer, dimension(:), allocatable, intent(out) :: n_dirs_for_EBS_along_pcbz_dir,neqv_dirs_pcbz, neqv_dirs_SCBZ, ncompl_dirs, n_irr_compl_dirs
type(crystal_3D), allocatable, intent(out) :: crystal_SC_reduced_to_prim_cell, &
                                              crystal_pc_reduced_to_prim_cell
logical, intent(out) :: pc_is_prim_cell, SC_is_prim_cell
type(crystal_3D), intent(in) :: crystal_pc, crystal_SC
real(kind=dp), dimension(:,:), intent(in) :: k_starts,k_ends
logical, intent(in), optional :: verbose
real(kind=dp), intent(in), optional :: symprec
type(symmetry_operation), dimension(:), allocatable :: symops_SC, symops_pc
type(eqv_bz_directions), dimension(:), allocatable :: dirs_eqv_in_SC_to_selec_pcbz_dir,dirs_eqv_in_pc_to_selec_pcbz_dir,complementary_pcdirs_for_pcdir
type(irr_bz_directions), dimension(:), allocatable :: irr_compl_pcdirs_for_pcdir
integer :: idir,ndirs,n_eqv_dirs_in_pcbz,n_eqv_dirs_in_SCBZ,n_req_dirs,i_req_dir,nsym_SC,nsym_pc
character(len=10) :: schoenflies_SC, schoenflies_pc
character(len=11) :: int_symb_pcbz, int_symb_SCBZ
logical :: print_stuff, try_to_reduce_SC_to_a_pc_for_symmetry, &
           try_to_reduce_pc_to_a_pc_for_symmetry
real(kind=dp) :: sym_prec

print_stuff = .TRUE.
if(present(verbose))then
    print_stuff = verbose
endif
sym_prec = default_symprec
if(present(symprec))then
    sym_prec = dabs(symprec)
endif
! Change the two variables below to .FALSE. if there are problems with symmetry
try_to_reduce_SC_to_a_pc_for_symmetry = .TRUE.
try_to_reduce_pc_to_a_pc_for_symmetry = .TRUE.


! Symm analysis, SC
write(*,*)
if(print_stuff)then
    write(*,'(A)')'Checking the symmetry of the SC...'
endif
call get_symm(symops=symops_SC,crystal=crystal_SC,nsym=nsym_SC,&
              schoenflies=schoenflies_SC,international_symb=int_symb_SCBZ, &
              cystal_is_prim_cell=SC_is_prim_cell, &
              try_to_reduce_to_a_pc=try_to_reduce_SC_to_a_pc_for_symmetry, &
              crystal_reduced_prim_cell = crystal_SC_reduced_to_prim_cell, &
              symprec=sym_prec)
if(print_stuff)then
    write(*,'(A)')         '    * Done.'
    write(*,'(A,I0,A)')    '    * Found ',nsym_SC,' symmetry operations for the SC.'
    write(*,'(5A)')        '      * Symmetry group of the SC: ',trim(adjustl(schoenflies_SC)),'(',trim(adjustl(int_symb_SCBZ)),').'
endif
! This always gets printed out, as it serves as a warning.
if(.not. SC_is_prim_cell)then
    write(*,'(A)')     '    * Your SC seems to be a perfect repetition of a smaller unit cell!'
    if(try_to_reduce_SC_to_a_pc_for_symmetry)then
        write(*,'(A)') '      * One of the possible corresponding primitive unit cells has been used to get the symmetry.'
    endif
endif

! Symm analysis, pc
write(*,*)
if(print_stuff)then
    write(*,'(A)')'Checking the symmetry of the reference unit cell...'
endif
call get_symm(symops=symops_pc, crystal=crystal_pc,nsym=nsym_pc,&
              schoenflies=schoenflies_pc,international_symb=int_symb_pcbz, &
              cystal_is_prim_cell=pc_is_prim_cell, &
              try_to_reduce_to_a_pc=try_to_reduce_pc_to_a_pc_for_symmetry, &
              crystal_reduced_prim_cell = crystal_pc_reduced_to_prim_cell, &
              symprec=sym_prec)
if(print_stuff)then
    write(*,'(A)')     '    * Done.'
    write(*,'(A,I0,A)')'    * Found ',nsym_pc,' symmetry operations for the pc.'
    write(*,'(5A)')    '      * Symmetry group of the pc: ',trim(adjustl(schoenflies_pc)),'(',trim(adjustl(int_symb_pcbz)),').'
endif
! This always gets printed out, as it is a warning.
if(.not. pc_is_prim_cell)then
    write(*,'(A)') "      WARNING: Your reference unit cell doesn't seem to be a primitive cell!"
    if(try_to_reduce_pc_to_a_pc_for_symmetry)then
        write(*,'(A)') '               * One of the possible corresponding primitive unit cells has been used to get the symmetry.'
    else
        write(*,'(A)') "               * It should be fine if that's really what you want."
        write(*,'(A)') "                 > It will certainly require more virtual memory though..."
    endif
endif

ndirs = size(k_starts(:,:),dim=1)
allocate(dirs_eqv_in_SC_to_selec_pcbz_dir(1:ndirs), dirs_eqv_in_pc_to_selec_pcbz_dir(1:ndirs))
!! Getting all the pcbz directions that are equivalent to each of the selected pcbz directions:
do idir=1,ndirs
    !! (1) By symmetry operations of the SCBZ
    call get_eqv_bz_dirs(eqv_dirs=dirs_eqv_in_SC_to_selec_pcbz_dir(idir),&
                         kstart=k_starts(idir,:),kend=k_ends(idir,:),symm_ops=symops_SC,crystal=crystal_SC)
    !! (2) By symmetry operations of the pcbz
    call get_eqv_bz_dirs(eqv_dirs=dirs_eqv_in_pc_to_selec_pcbz_dir(idir),&
                         kstart=k_starts(idir,:),kend=k_ends(idir,:),symm_ops=symops_pc,crystal=crystal_pc)
enddo

!! Getting all the pcbz directions that are:
!!     >>> Eqv. to the selected ones by sym. ops. of the pcbz 
!!                         and 
!!     >>> NOT eqv. to the selected ones by sym. ops. of the SCBZ
!! I call them >>> complementary pcbz directions <<<
allocate(complementary_pcdirs_for_pcdir(1:ndirs),ncompl_dirs(1:ndirs))
do idir=1,ndirs
    call get_compl_pcbz_direcs(compl_pcdirs=complementary_pcdirs_for_pcdir(idir), ncompl=ncompl_dirs(idir), &
                               dirs_eqv_in_pcbz=dirs_eqv_in_pc_to_selec_pcbz_dir(idir), &
                               dirs_eqv_in_SCBZ=dirs_eqv_in_SC_to_selec_pcbz_dir(idir),symops_SC=symops_SC)
enddo
!! Getting now the irreducible complementary pcbz directions
allocate(irr_compl_pcdirs_for_pcdir(1:ndirs),n_irr_compl_dirs(1:ndirs))
do idir=1,ndirs
    if(ncompl_dirs(idir) > 0)then
        call get_irr_bz_directions(irr_dirs=irr_compl_pcdirs_for_pcdir(idir), nirr=n_irr_compl_dirs(idir), &
                                   list_of_dirs=complementary_pcdirs_for_pcdir(idir), &
                                   symm_ops=symops_SC, crystal=crystal_SC)
    else
        n_irr_compl_dirs(idir) = 0
    endif
enddo
!! Finally getting all directions required for a symmetry-averaged EBS along the selected pcbz directions
allocate(neqv_dirs_pcbz(1:ndirs), neqv_dirs_SCBZ(1:ndirs))
allocate(dirs_req_for_symmavgd_EBS_along_pcbz_dir(1:ndirs))
allocate(n_dirs_for_EBS_along_pcbz_dir(1:ndirs))
do idir=1,ndirs
    n_eqv_dirs_in_pcbz = size(dirs_eqv_in_pc_to_selec_pcbz_dir(idir)%eqv_dir(:))
    n_req_dirs = 1 + n_irr_compl_dirs(idir)
    n_dirs_for_EBS_along_pcbz_dir(idir) = n_req_dirs

    allocate(dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(1:n_req_dirs))
    
    dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(1)%kstart = k_starts(idir,:) 
    dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(1)%kend = k_ends(idir,:) 
    n_eqv_dirs_in_SCBZ = size(dirs_eqv_in_SC_to_selec_pcbz_dir(idir)%eqv_dir(:))

    dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(1)%neqv = n_eqv_dirs_in_SCBZ
    dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(1)%weight = &
        real(n_eqv_dirs_in_pcbz - ncompl_dirs(idir),kind=dp)/real(n_eqv_dirs_in_pcbz,kind=dp)
    if(n_req_dirs > 1)then
        do i_req_dir=2,n_req_dirs
            dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(i_req_dir) = irr_compl_pcdirs_for_pcdir(idir)%irr_dir(i_req_dir-1)
            n_eqv_dirs_in_SCBZ = dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(i_req_dir)%neqv
            dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(i_req_dir)%weight = &
                real(irr_compl_pcdirs_for_pcdir(idir)%irr_dir(i_req_dir-1)%neqv,kind=dp)/real(n_eqv_dirs_in_pcbz,kind=dp)
        enddo
    endif

    if(abs(sum(dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(:)%weight) - 1.0_dp) > 1E-3_dp)then
        write(*,"(A,I0,A)")  "ERROR: The sum of the weights of the pcbz directions that are equivalent to the selected pcbz direction #",idir," doesn't equal 1."
        write(*,"(A)")       "       This is necessary for the symmetry-averaged EBS."
        write(*,"(A,f0.3,A)")"       The sum was ",sum(dirs_req_for_symmavgd_EBS_along_pcbz_dir(idir)%irr_dir(:)%weight)," instead."
        write(*,"(A)")       "Stopping now."
        stop
    endif
    neqv_dirs_pcbz(idir) = n_eqv_dirs_in_pcbz
    neqv_dirs_SCBZ(idir) = size(dirs_eqv_in_SC_to_selec_pcbz_dir(idir)%eqv_dir)
enddo


end subroutine get_all_irr_dirs_req_for_symmavgd_EBS


end module math
