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
! MODULE: constants_and_types 
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION:
!> Provides most constants and derived types used in BandUP.
!==============================================================================

module constants_and_types
implicit none
PRIVATE
! Public parameters
! Character
PUBLIC :: package_version
! Integer
PUBLIC :: sp, dp, long_int_kind, kind_cplx_coeffs, qe_dp, str_len
! Real
PUBLIC :: pi, twopi, min_dk, default_tol_for_vec_equality, &
          max_tol_for_vec_equality, default_tol_for_int_commens_test, &
          default_symprec, two_m_over_hbar_sqrd, identity_3D, Hartree, Ry, bohr
! Logical
PUBLIC :: calc_spec_func_explicitly, stop_if_pckpt_cannot_be_parsed, &
          stop_if_GUR_fails, get_all_kpts_needed_for_EBS_averaging, &
          print_GUR_pre_unfolding_utility, renormalize_wf 

! Derived types
PUBLIC :: timekeeping, vec3d, vec3d_int, symmetry_operation, crystal_3D, &
          star_point_properties, star, bz_direction, eqv_bz_directions, &
          irr_bz_directions, selected_pcbz_directions, &
          geom_unfolding_relations_for_each_SCKPT, UnfoldedQuantities, & 
          UnfoldedQuantitiesForOutput, comm_line_args, pw_wavefunction, &
          UnfoldDensityOpContainer, MatrixIndices

! I only change the following options for debugging/testing. 
! You probably shouldn't modify these either.
integer, parameter :: sp = selected_real_kind(6, 37), &    ! Single precision
                      dp = selected_real_kind(15, 307), &  ! Double precision
                      long_int_kind = selected_int_kind(10), &
                      qe_dp = kind(1.0_dp), & ! QE's double precision
                      kind_cplx_coeffs = sp, & 
                      str_len=256
real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp), twopi = 2.0_dp*pi, &
                            default_tol_for_vec_equality=1E-5_dp, &
                            max_tol_for_vec_equality=1E-3_dp, &
                            min_dk=2.0_dp*default_tol_for_vec_equality, &
                            default_tol_for_int_commens_test=1E-5_dp, & 
                            default_symprec=1E-5_dp, &
                            ! c = 2m/hbar**2 in 1/eV Ang^2 (from WaveTrans)
                            two_m_over_hbar_sqrd = 0.262465831, & 
                            ! Energy convertion to eV
                            Ry =  13.60569172, Hartree = 27.2113834, & 
                            bohr = 0.52917721092 ! Length conv. to Angstrom
real(kind=dp), dimension(1:3,1:3), parameter :: identity_3D = &
    reshape(real((/ 1, 0, 0, 0, 1, 0, 0, 0, 1 /), kind=dp), shape(identity_3D))
logical, parameter :: calc_spec_func_explicitly = .FALSE., &
                      stop_if_pckpt_cannot_be_parsed = .TRUE., &
                      stop_if_GUR_fails = .TRUE., &
                      get_all_kpts_needed_for_EBS_averaging = .TRUE., &
                      print_GUR_pre_unfolding_utility = .FALSE., &
                      renormalize_wf = .TRUE.
! Code version. If you modify this with no good reason and without following
! what is described in the CONTRIBUTING.md file, then it will be
! *very* difficult for me to offer you any support.
character(len=str_len), parameter :: tag_for_push="3.0.0-beta.6"
#if defined (__COMMIT_TAG__)
    character(len=str_len), parameter :: &
        version_tag=trim(adjustl(__COMMIT_TAG__))
#else
    character(len=str_len), parameter :: version_tag='v'//tag_for_push
#endif
#if defined (__COMMIT_DATE__)
    character(len=str_len), parameter :: &
        commit_date = trim(adjustl(__COMMIT_DATE__))
#else
    character(len=30), parameter :: commit_date=" (unknown commit date)"
#endif
character(len=str_len), parameter :: &
    package_version = trim(adjustl(version_tag)) // &
                      ", " // trim(adjustl(commit_date))
!! Derived type definitions
type :: comm_line_args
    character(len=256) :: WF_file, input_file_prim_cell, input_file_supercell,&
                          input_file_pc_kpts, input_file_energies, &
                          out_file_SC_kpts, &
                          output_file_symm_averaged_EBS, &
                          output_file_only_user_selec_direcs, &
                          output_file_symm_averaged_unf_dens_op, &
                          output_file_only_user_selec_direcs_unf_dens_op, &
                          pw_code, qe_outdir, qe_prefix, abinit_files_file, &
                          castep_seed 
    integer :: spin_channel, n_sckpts_to_skip
    real(kind=dp), dimension(1:3) :: saxis, normal_to_proj_plane, &
                                     origin_for_spin_proj_cartesian, &
                                     origin_for_spin_proj_rec
    logical :: stop_if_not_commensurate, &
               write_attempted_pc_corresp_to_input_pc, &
               write_attempted_pc_corresp_to_SC, no_symm_avg, no_symm_sckpts, &
               perform_unfold, origin_for_spin_proj_passed_in_rec, &
               continue_if_npw_smaller_than_expected, write_unf_dens_op
end type comm_line_args

type :: timekeeping 
    real(kind=dp) :: start, end, read_wf, calc_spec_weights, calc_SF, calc_dN,&
                     calc_rho, calc_pauli_vec, calc_pauli_vec_projs, &
                     write_dN_files, write_unf_dens_op_files
end type timekeeping 

type :: vec3d
  real(kind=dp), dimension(1:3) :: coord
end type vec3d

type :: vec3d_int
  integer, dimension(1:3) :: coord
end type vec3d_int

type :: MatrixIndices
    integer, dimension(1:2) :: indices
end type MatrixIndices

type :: pw_wavefunction
    integer :: i_spin, n_pw, n_spin, n_bands, n_spinor, &
               n_bands_up, n_bands_down
    real(kind=dp) :: encut, Vcell, e_fermi, e_fermi_up, e_fermi_down
    real(kind=dp), dimension(1:3) :: kpt_frac_coords, kpt_cart_coords
    ! Direct and reciprocal lattice vectors
    real(kind=dp), dimension(1:3,1:3) :: A_matrix, B_matrix
    real(kind=dp), dimension(:), allocatable :: band_energies, band_occupations
    ! G(ipw)%coord(1:3) := coords of RL vec associated with pw_coeff(ipw)
    type(vec3d), dimension(:), allocatable :: G_cart
    type(vec3d_int), dimension(:), allocatable :: G_frac
    complex(kind=kind_cplx_coeffs), dimension(:,:,:), allocatable :: pw_coeffs
    logical :: is_spinor, two_efs
end type pw_wavefunction

type :: symmetry_operation
    integer, dimension(1:3) :: translation_fractional_coords
    integer, dimension(1:3,1:3) :: rotation_fractional_coords
    real(kind=dp), dimension(1:3) :: translation_cartesian_coords
    real, dimension(1:3,1:3) :: rotation_cartesian_coords
    real(kind=dp), dimension(1:3,1:3) :: basis
end type symmetry_operation

type :: crystal_3D
    character(len=256) :: description=''
    real(kind=dp), dimension(1:3,1:3) :: latt_vecs, rec_latt_vecs
    real(kind=dp) :: vol, rec_latt_vol
    real(kind=dp), dimension(:,:), allocatable :: &
        ! Shall be dimension(:,1:3)
        coords_basis_atoms, fractional_coords_basis_atoms
    logical, dimension(:,:), allocatable :: unconstrained_dof_basis_atoms
    character(len=3), dimension(:), allocatable :: atomic_symbols_basis_atoms
    integer, dimension(:), allocatable :: integer_types_basis_atoms
    ! Symmetry
    character(len=10) :: schoenflies
    character(len=11) :: international_symb
    integer :: space_group_num, nsym
    type(symmetry_operation), dimension(:), allocatable :: symops
    logical :: is_prim_cell
    type(crystal_3D), pointer :: corresp_pc => null()
end type crystal_3D

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
    real(kind=dp) :: weight=0.0_dp
    integer :: neqv
end type bz_direction

type :: eqv_bz_directions
    integer :: neqv=1 ! Number of eqv. directions in the set.
    type(bz_direction), dimension(:), allocatable :: eqv_dir
end type eqv_bz_directions

type :: irr_bz_directions
    type(bz_direction), dimension(:), allocatable :: irr_dir
    integer :: neqv=1, neqv_SCBZ=1, ncompl_dirs=0, n_irr_compl_dirs=0
end type irr_bz_directions

type :: trial_folding_pckpt
    ! Message from Paulo:
    ! The prefix "S" means "symmetrized". This is a trick that allows the use
    ! of the the coefficients of a SC wavefunction psi(K',n) to calculate the
    ! spectral weights associated with a SC wavefunction psi(K,n), where
    ! K' = SK and S is a symmetry operation of the crystal's point group.
    real(kind=dp), dimension(1:3) :: coords_actual_unfolding_K, coords, &
                                     coords_SCKPT_used_for_coeffs, &
                                     Scoords, Sfolding_vec, &
                                     Sorigin_for_spin_proj
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

type :: GUR_indices
    integer :: i_SCKPT, i_selec_pcbz_dir, i_needed_dirs, ipc_kpt
end type  GUR_indices

type :: geom_unfolding_relations_for_each_SCKPT
    type(selected_pcbz_directions), dimension(:), allocatable :: SCKPT
    logical, dimension(:), allocatable :: SCKPT_used_for_unfolding
    integer :: n_pckpts,n_folding_pckpts
    type(vec3d), dimension(:), allocatable :: list_of_SCKPTS
    type(GUR_indices) :: current_index
    real(kind=dp), dimension(1:3,1:3) :: B_matrix_SC, b_matrix_pc
end type geom_unfolding_relations_for_each_SCKPT


!! Defining derived types to support a "UnfoldedQuantities" type.
type :: UnfoldDensityOpContainer
    complex(kind=kind_cplx_coeffs), dimension(:), allocatable :: rho
    type(MatrixIndices), dimension(:), allocatable :: band_indices
    integer :: iener_in_full_pc_egrid, nbands
end type UnfoldDensityOpContainer

type :: unfolded_quantities_for_given_pckpt
    !! This type holds the unfolded quantities, for a given pc-kpt, 
    !! at ever point of the energy grid.
    ! delta_Ns and spectral functions at each point of the energy grid
    real(kind=dp), dimension(:), allocatable :: dN, SF, &
                                                spin_proj_perp, spin_proj_para
    ! sigma(1:nener,1:3) holds the expected values of the Pauli matrices 
    ! sigma_i; i=x, y, z (used if spinor WF)
    real(kind=dp), dimension(:,:), allocatable :: sigma
    ! rho will hold the unfolding density operator, if needed
    ! The dimension(:) statement is needed because there is one such operator
    ! for each energy E in the PC energy grid such that delta_N(k,E)>0.
    type(UnfoldDensityOpContainer), dimension(:), allocatable :: rhos 
end type unfolded_quantities_for_given_pckpt

type :: list_of_pckpts_for_unfolded_quantities
    type(unfolded_quantities_for_given_pckpt), dimension(:), allocatable :: &
        pckpt
end type list_of_pckpts_for_unfolded_quantities

type :: UnfoldedQuantities_info_for_needed_pcbz_dirs
    type(list_of_pckpts_for_unfolded_quantities), dimension(:), allocatable ::&
        needed_dir
end type UnfoldedQuantities_info_for_needed_pcbz_dirs
!! Defining now the "UnfoldedQuantities" type
!! It holds info about unfodled properties 
!! [delta_Ns (dN), spectral functions (SF), and all other components of the
!! previously defined type "unfolded_quantities_for_given_pckpt"]
!! The structure of a variable unf of type(UnfoldedQuantities) is:
!! unf%selec_pcbz_dir(i_selec_pcbz_dir)%needed_dir(i_needed_dirs)%&
!!     pckpt(ipc_kpt)%PROPERTY(iener)
type :: UnfoldedQuantities
    integer :: n_SC_bands
    type(UnfoldedQuantities_info_for_needed_pcbz_dirs), dimension(:), &
        allocatable :: &
        selec_pcbz_dir
end type UnfoldedQuantities

type :: UnfoldedQuantitiesForOutput
    integer :: n_SC_bands
    type(list_of_pckpts_for_unfolded_quantities), dimension(:), allocatable ::&
        pcbz_dir
end type UnfoldedQuantitiesForOutput

end module constants_and_types
