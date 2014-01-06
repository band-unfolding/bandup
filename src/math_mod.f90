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

module math
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: sp, dp, qp, pi, twopi, vec3d, &
          time,cross, norm, angle, real_seq, delta,  &
          get_rec_latt, coords_cart_vec_in_new_basis, vec_in_latt, &
          reduce_point_to_bz, same_vector, divide_SCKPTS_evenly_among_MPI_tasks 
integer, parameter :: sp = selected_real_kind(6, 37)    ! Single precision
integer, parameter :: dp = selected_real_kind(15, 307)  ! Double precision
integer, parameter :: qp = selected_real_kind(33, 4931) ! Quadruple precision

type :: vec3d
  real(kind=dp), dimension(1:3) :: coord
end type vec3d
real(kind=dp), parameter :: pi = 4.*atan(1.), twopi = 2.*pi
real(kind=dp), parameter :: two_m_over_hbar_sqrd = 0.262465831 ! c = 2m/hbar**2 in units of 1/eV Ang^2 (from WaveTrans)

CONTAINS

function time() result(rtn)
implicit none
real(kind=dp) :: rtn
integer :: count, count_rate

!$ rtn = omp_get_wtime()
!$ return
   CALL system_clock(count, count_rate)
   rtn = count/count_rate

end function time

function cross(b,c) result(a)
implicit none
  real(kind=dp), dimension(1:3) :: a, b, c

  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)

end function cross

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

    n_terms = nint((last_term - first_term)/increment + 1.0d0)
    allocate(return_list(1:n_terms))
    do i = 1, n_terms
         return_list(i) = first_term + real(i - 1, kind=dp)*increment
    enddo

end subroutine real_seq


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

    rtn = exp(-((x/stddev)**2.0d0)/2.0d0)/(stddev*dsqrt(twopi))

end function delta

function integral_delta_of_x_minus_x0(x0,low_x,upp_x) result(rtn)
! Consistently with the employed approx. by a Gaussian distribution
implicit none
real(kind=dp) :: rtn
real(kind=dp), intent(in) :: x0,low_x,upp_x

   rtn = 0.5d0*(erf(upp_x-x0) - erf(low_x-x0))

end function integral_delta_of_x_minus_x0


function inverse_of_3x3_matrix(matrix, success) result(inverse_matrix)
implicit none
real(kind=dp), dimension(3,3) :: inverse_matrix
real(kind=dp), dimension(3,3), intent(in) :: matrix
logical, intent(out), optional :: success
real(kind=dp), dimension(3,3) :: cofactor
real(kind=dp) :: det
real(kind=dp), parameter :: eps = 1.0D-10

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
! Checks if v1 and v2 are the same vector (within a given tolerance)
implicit none
logical :: rtn
real(kind=dp), dimension(1:3), intent(in) :: v1, v2
real(kind=dp), intent(in), optional :: tolerance
real(kind=dp) :: tol,dist 
integer :: i

    tol = real(5E-3, kind=dp)
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
implicit none
logical :: rtn
real*8, dimension(1:3), intent(in) :: vec
real*8, dimension(1:3,1:3), intent(in) :: latt
real*8, intent(in), optional :: tolerance
real*8, dimension(1:3) :: reduced_vec,frac_coords,reduced_frac_coords,g
real*8 :: tol
integer :: i,ig1,ig2,ig3,ig1p,ig2p,ig3p,nb1max,nb2max,nb3max

tol = 1E-5
if(present(tolerance))then
     tol = dabs(tolerance)
endif

rtn = .FALSE.
frac_coords(:) = coords_cart_vec_in_new_basis(cart_vec=vec, new_basis=latt)
reduced_frac_coords(:) = mod(frac_coords(:),1.0d0)
reduced_vec(:) = 0.0d0
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
        if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1 ! Trick to make the indices ig2p vary from 0 to nb2max and then from -nb2max to -1
        do ig1=0,2*nb1max
            ig1p=ig1
            if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1 ! Trick to make the indices ig1p vary from 0 to nb1max and then from -nb1max to -1
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
    origin(:) = 0.0d0
endif

do ig3=-1,1
    do ig2=-1,1
        do ig1=-1,1

            if(ig1==0 .and. ig2==0 .and. ig3==0) cycle
            g(:) = origin(:) + ig1*rec_latt(1,:) + ig2*rec_latt(2,:) + ig3*rec_latt(3,:)
            if(norm(point(:) - g(:)) < norm(point(:) - origin(:)))then
                rtn = .FALSE.
                return
            endif

        enddo
    enddo
enddo

end function point_is_in_bz


subroutine reduce_point_to_bz(point,rec_latt,point_reduced_to_bz,frac_coords_reduc_vec)
implicit none
real(kind=dp), dimension(1:3), intent(in) :: point
real(kind=dp), dimension(1:3,1:3), intent(in) :: rec_latt
real(kind=dp), dimension(1:3), intent(out) :: point_reduced_to_bz
real(kind=dp), dimension(1:3), optional, intent(out) :: frac_coords_reduc_vec
logical :: point_has_been_reduced
real(kind=dp), dimension(1:3) :: frac_coords_point, g, reduc_vec
integer :: icoord, ig1, ig2, ig3

    point_has_been_reduced = .FALSE.
    if(.not.point_is_in_bz(point,rec_latt))then
        frac_coords_point = coords_cart_vec_in_new_basis(cart_vec=point, new_basis=rec_latt)
        point_reduced_to_bz(:) = 0.0d0
        do icoord=1,3
            point_reduced_to_bz(:) = point_reduced_to_bz(:) + modulo(frac_coords_point(icoord),1.0d0)*rec_latt(icoord,:)
        enddo
        do ig3=-1,1
            do ig2=-1,1
                do ig1=-1,1
                    g(:) = ig1*rec_latt(1,:) + ig2*rec_latt(2,:) + ig3*rec_latt(3,:)
                    if(point_is_in_bz(point_reduced_to_bz - g(:),rec_latt))then
                        point_reduced_to_bz(:) = point_reduced_to_bz(:) - g(:)
                        point_has_been_reduced = .TRUE.
                        if(present(frac_coords_reduc_vec))then
                            reduc_vec(:) = point(:) - point_reduced_to_bz(:)
                            frac_coords_reduc_vec(:) = coords_cart_vec_in_new_basis(cart_vec=reduc_vec,new_basis=rec_latt)
                        endif
                        exit
                    endif
                enddo
            enddo
        enddo
    else
        point_reduced_to_bz(:) = point(:)
        reduc_vec(:) = 0.0d0
        point_has_been_reduced = .TRUE.
    endif
    if(.not.point_has_been_reduced)then
        write(*,*)'ERROR (reduce_point_to_bz): Point has not been reduced. Stopping now.'
        stop
    endif

end subroutine reduce_point_to_bz


subroutine divide_SCKPTS_evenly_among_MPI_tasks(SCKPTS_for_task,nkpts,ntasks)
implicit none
integer, intent(in) :: nkpts,ntasks
integer, dimension(:,:), allocatable, intent(out) :: SCKPTS_for_task
integer :: max_n_SCKPTS_per_MPI_task, task_kpt, i_SCKPT, itask

max_n_SCKPTS_per_MPI_task = ceiling(real(nkpts)/real(ntasks))
allocate(SCKPTS_for_task(0:ntasks-1,1:max_n_SCKPTS_per_MPI_task))
SCKPTS_for_task(:,:) = 0

task_kpt = 0
do i_SCKPT=1,nkpts
    itask = modulo(i_SCKPT-1,ntasks)
    if(itask==0) task_kpt = task_kpt + 1
    SCKPTS_for_task(itask,task_kpt) = i_SCKPT
enddo

end subroutine divide_SCKPTS_evenly_among_MPI_tasks

end module math
