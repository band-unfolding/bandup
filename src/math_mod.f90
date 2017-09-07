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
! MODULE: math 
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Defines some mathematical functions and subroutines used by BandUP.
!==============================================================================

module math
use constants_and_types, only : sp, dp, kind_cplx_coeffs, &
                                default_tol_for_vec_equality, twopi, &
                                UnfoldDensityOpContainer
use lists_and_seqs, only : list_index
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: n_digits_integer, cross, norm, angle, delta, &
          integral_delta_x_minus_x0, triple_product, &
          coords_cart_vec_in_new_basis, same_vector, &
          versor, trace_AB, trace, inverse_of_3x3_matrix

interface trace_AB
    module procedure trace_AB_real, trace_AB_cplx, trace_AB_A_is_rho_B_is_cplx
end interface trace_AB

interface trace
    module procedure trace_cplx
end interface trace

interface norm
    module procedure norm_sp, norm_dp, cplx_norm_dp, cplx_norm_sp
end interface norm

CONTAINS

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


function norm_sp(v) result(rtn)
implicit none
  real(kind=sp), dimension(:) :: v
  real(kind=sp) :: rtn

  rtn = sqrt(dot_product(v,v))

end function norm_sp


function cplx_norm_sp(v) result(rtn)
implicit none
  complex(kind=sp), dimension(:) :: v
  real(kind=sp) :: rtn

  rtn = sqrt(abs(dot_product(v,v)))

end function cplx_norm_sp

function norm_dp(v) result(rtn)
implicit none
  real(kind=dp), dimension(:) :: v
  real(kind=dp) :: rtn

  rtn = sqrt(dot_product(v,v))

end function norm_dp


function cplx_norm_dp(v) result(rtn)
implicit none
  complex(kind=dp), dimension(:) :: v
  real(kind=dp) :: rtn

  rtn = sqrt(abs(dot_product(v,v)))

end function cplx_norm_dp


function versor(v) result(rtn)
implicit none
real(kind=dp), dimension(1:3), intent(in) :: v
real(kind=dp), dimension(1:3) :: rtn

    rtn = 0.0_dp
    if(norm(v) > 0.0_dp)then
        rtn = v/norm(v)
    endif
    return

end function versor

function angle(v1, v2) result(rtn)
implicit none
  real(kind=dp), dimension(1:3) :: v1, v2
  real(kind=dp) :: rtn

  rtn = acos(dot_product(v1, v2)/(norm(v1)*norm(v2)))

end function angle

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


function coords_cart_vec_in_new_basis(cart_vec, new_basis) result(coords)
! Returns the coordinates of the vector cart_vec = (x,y,z)
! in the basis "new_basis"
real(kind=dp), dimension(1:3) :: coords
real(kind=dp), dimension(1:3), intent(in) :: cart_vec
real(kind=dp), dimension(1:3,1:3), intent(in) :: new_basis
real(kind=dp), dimension(1:3,1:3) :: aux_matrix
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


function trace_AB_real(A,B) result(rtn)
real(kind=dp) :: rtn
real(kind=dp), dimension(:,:), intent(in) :: A, B
integer :: msize, i

    rtn = 0.0_dp
    msize = size(A, dim=1)    
    !$omp parallel do default(none) schedule(guided) &
    !$omp private(i) &
    !$omp shared(msize, A, B) &
    !$omp reduction(+:rtn)
    do i=1, msize
        rtn = rtn + sum(A(i,:) * B(:,i))
    enddo

end function trace_AB_real


function trace_AB_cplx(A,B) result(rtn)
complex(kind=kind_cplx_coeffs) :: rtn
complex(kind=kind_cplx_coeffs), dimension(:,:), intent(in) :: A, B
integer :: msize, i

    rtn = 0.0_kind_cplx_coeffs
    msize = size(A, dim=1)    
    !$omp parallel do default(none) schedule(guided) &
    !$omp private(i) &
    !$omp shared(msize, A, B) &
    !$omp reduction(+:rtn)
    do i=1, msize
        rtn = rtn + sum(A(i,:) * B(:,i))
    enddo

end function trace_AB_cplx


function trace_AB_A_is_rho_B_is_cplx(A,B) result(rtn)
complex(kind=kind_cplx_coeffs) :: rtn
type(UnfoldDensityOpContainer), intent(in) :: A
complex(kind=kind_cplx_coeffs), dimension(:,:), intent(in) :: B
integer :: msize, i, j, linear_index

    rtn = 0.0_kind_cplx_coeffs
    msize = size(B, dim=1)    
    !$omp parallel do default(none) schedule(guided) &
    !$omp private(i,j,linear_index) &
    !$omp shared(msize, A, B) &
    !$omp reduction(+:rtn)
    do i=1, msize
        do j=1, msize
            linear_index = list_index([i,j], A%band_indices)
            if(linear_index<1) cycle
            rtn = rtn + A%rho(linear_index) * B(j,i)
        enddo
    enddo

end function trace_AB_A_is_rho_B_is_cplx


function trace_AB_hermit(A,B) result(rtn)
real(kind=dp) :: rtn
complex(kind=kind_cplx_coeffs), dimension(:,:), intent(in) :: A, B
integer :: msize, i

    rtn = 0.0_dp
    msize = size(A, dim=1)    
    !$omp parallel do default(none) schedule(guided) &
    !$omp private(i) &
    !$omp shared(msize, A, B) &
    !$omp reduction(+:rtn)
    do i=1, msize
        rtn = rtn + sum(real(A(i,:)) * real(B(:,i))) - &
                    sum(aimag(A(i,:)) * aimag(B(:,i)))
    enddo

end function trace_AB_hermit


function trace_cplx(A) result(rtn)
complex(kind=kind_cplx_coeffs) :: rtn
complex(kind=kind_cplx_coeffs), dimension(:,:), intent(in) :: A
integer :: i

    rtn = 0.0_kind_cplx_coeffs
    !$omp parallel do default(none) schedule(guided) &
    !$omp private(i) &
    !$omp shared(A) &
    !$omp reduction(+:rtn)
    do i=1, size(A, dim=1)
        rtn = rtn + A(i,i)
    enddo

end function trace_cplx

end module math
