#define CPASSERT(cond) IF (.NOT. (cond)) ERROR STOP 'CPASSERT'
MODULE parallel_kernel

 USE kinds, ONLY: dp
 USE qs_grid_atom, ONLY: grid_atom_type
 USE qs_harmonics_atom, ONLY: harmonics_atom_type
 USE orbital_pointers, ONLY: indco, indso, nco, ncoset, nsoset
 USE orbital_transformation_matrices, ONLY: orbtramat
 USE spherical_harmonics, ONLY: y_lm
 USE cell_types, ONLY: cell_type
 USE particle_types, ONLY: particle_type
 IMPLICIT NONE
 CONTAINS
   SUBROUTINE radial_node_derivative_coefficients( &
      grid_atom, descending, node, logical_start, slope, curvature)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      LOGICAL, INTENT(IN)                                :: descending
      INTEGER, INTENT(IN)                                :: node, logical_start
      REAL(dp), DIMENSION(4), INTENT(OUT)                :: slope, curvature

      INTEGER                                            :: center, hi, lo, n
      REAL(dp)                                           :: h_hi, h_lo, x_center, x_hi, x_lo

      n = grid_atom%nr
      slope = 0.0_dp
      curvature = 0.0_dp
      IF (node == 1) THEN
         lo = 1
         hi = 2
      ELSE IF (node == n) THEN
         lo = n - 1
         hi = n
      ELSE
         lo = node - 1
         hi = node + 1
      END IF
      x_lo = grid_atom%rad(MERGE(n + 1 - lo, lo, descending))
      x_hi = grid_atom%rad(MERGE(n + 1 - hi, hi, descending))
      slope(lo - logical_start + 1) = -1.0_dp/(x_hi - x_lo)
      slope(hi - logical_start + 1) = 1.0_dp/(x_hi - x_lo)

      IF (n < 3) RETURN
      center = MIN(MAX(node, 2), n - 1)
      lo = center - 1
      hi = center + 1
      x_lo = grid_atom%rad(MERGE(n + 1 - lo, lo, descending))
      x_center = grid_atom%rad(MERGE(n + 1 - center, center, descending))
      x_hi = grid_atom%rad(MERGE(n + 1 - hi, hi, descending))
      h_lo = x_center - x_lo
      h_hi = x_hi - x_center
      curvature(lo - logical_start + 1) = 2.0_dp/(h_lo*(h_lo + h_hi))
      curvature(center - logical_start + 1) = &
         -2.0_dp*(1.0_dp/h_lo + 1.0_dp/h_hi)/(h_lo + h_hi)
      curvature(hi - logical_start + 1) = 2.0_dp/(h_hi*(h_lo + h_hi))

   END SUBROUTINE radial_node_derivative_coefficients
   SUBROUTINE atom_grid_interpolation_weights( &
      grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
      radial_derivative_weights, nradial, angular_weights, angular_derivative_weights, active)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      REAL(dp), INTENT(IN)                               :: cutoff
      INTEGER, DIMENSION(4), INTENT(OUT)                 :: radial_indices
      REAL(dp), DIMENSION(4), INTENT(OUT)                :: radial_weights
      REAL(dp), DIMENSION(4), INTENT(OUT), OPTIONAL      :: radial_derivative_weights
      INTEGER, INTENT(OUT)                               :: nradial
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: angular_weights
      REAL(dp), DIMENSION(:, :), INTENT(OUT), OPTIONAL   :: angular_derivative_weights
      LOGICAL, INTENT(OUT)                               :: active

      INTEGER                                            :: ia, ic, inode, iso, l, left, left_pos, &
                                                            logical_end, logical_start, lx, ly, &
                                                            lz, n, right_pos, shell_index
      LOGICAL                                            :: descending
      REAL(dp)                                           :: dh00, dh01, dh10, dh11, dh20, dh21, h00, &
                                                            h01, h10, h11, h20, h21, h_interval, &
                                                            monomial, radius, solid_derivative, t, &
                                                            x1, x2
      REAL(dp), DIMENSION(3)                             :: direction
      REAL(dp), DIMENSION(harmonics%max_s_harm)          :: angular_values
      REAL(dp), DIMENSION(4)                             :: curvature_left, curvature_right, &
                                                            slope_left, slope_right
      REAL(dp), DIMENSION(3, harmonics%max_s_harm)       :: angular_value_derivatives

      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(harmonics))
      CPASSERT(SIZE(angular_weights) == grid_atom%ng_sphere)
      IF (PRESENT(angular_derivative_weights)) THEN
         CPASSERT(SIZE(angular_derivative_weights, 1) == 3)
         CPASSERT(SIZE(angular_derivative_weights, 2) == grid_atom%ng_sphere)
      END IF

      radial_indices = 0
      radial_weights = 0.0_dp
      IF (PRESENT(radial_derivative_weights)) radial_derivative_weights = 0.0_dp
      angular_weights = 0.0_dp
      IF (PRESENT(angular_derivative_weights)) angular_derivative_weights = 0.0_dp
      nradial = 0
      active = .FALSE.
      n = grid_atom%nr
      descending = grid_atom%rad(1) > grid_atom%rad(n)
      radius = SQRT(SUM(displacement**2))
      IF (radius > cutoff .OR. &
          radius > grid_atom%rad(MERGE(1, n, descending))) RETURN
      IF (radius <= 1.0E-12_dp) RETURN
      direction = displacement/radius

      CPASSERT(n >= 2)
      left = n - 1
      IF (radius <= grid_atom%rad(MERGE(n, 1, descending))) THEN
         left = 1
      ELSE
         DO inode = 1, n - 1
            IF (radius <= grid_atom%rad(MERGE(n - inode, inode + 1, descending))) THEN
               left = inode
               EXIT
            END IF
         END DO
      END IF

      logical_start = MAX(1, left - 1)
      logical_end = MIN(n, left + 2)
      nradial = logical_end - logical_start + 1
      DO inode = 1, nradial
         radial_indices(inode) = MERGE(n + 2 - logical_start - inode, &
                                       logical_start + inode - 1, descending)
      END DO
      left_pos = left - logical_start + 1
      right_pos = left_pos + 1
      x1 = grid_atom%rad(radial_indices(left_pos))
      x2 = grid_atom%rad(radial_indices(right_pos))
      h_interval = x2 - x1
      t = (radius - x1)/h_interval

      h00 = 1.0_dp - 10.0_dp*t**3 + 15.0_dp*t**4 - 6.0_dp*t**5
      h10 = t - 6.0_dp*t**3 + 8.0_dp*t**4 - 3.0_dp*t**5
      h20 = 0.5_dp*(t**2 - 3.0_dp*t**3 + 3.0_dp*t**4 - t**5)
      h01 = 10.0_dp*t**3 - 15.0_dp*t**4 + 6.0_dp*t**5
      h11 = -4.0_dp*t**3 + 7.0_dp*t**4 - 3.0_dp*t**5
      h21 = 0.5_dp*(t**3 - 2.0_dp*t**4 + t**5)

      CALL radial_node_derivative_coefficients( &
         grid_atom, descending, left, logical_start, slope_left, curvature_left)
      CALL radial_node_derivative_coefficients( &
         grid_atom, descending, left + 1, logical_start, slope_right, curvature_right)
      radial_weights(left_pos) = radial_weights(left_pos) + h00
      radial_weights(right_pos) = radial_weights(right_pos) + h01
      radial_weights(1:nradial) = radial_weights(1:nradial) + &
                                  h_interval*(h10*slope_left(1:nradial) + &
                                              h11*slope_right(1:nradial)) + &
                                  h_interval**2*(h20*curvature_left(1:nradial) + &
                                                 h21*curvature_right(1:nradial))
      IF (PRESENT(radial_derivative_weights)) THEN
         dh00 = -30.0_dp*t**2 + 60.0_dp*t**3 - 30.0_dp*t**4
         dh10 = 1.0_dp - 18.0_dp*t**2 + 32.0_dp*t**3 - 15.0_dp*t**4
         dh20 = 0.5_dp*(2.0_dp*t - 9.0_dp*t**2 + 12.0_dp*t**3 - 5.0_dp*t**4)
         dh01 = 30.0_dp*t**2 - 60.0_dp*t**3 + 30.0_dp*t**4
         dh11 = -12.0_dp*t**2 + 28.0_dp*t**3 - 15.0_dp*t**4
         dh21 = 0.5_dp*(3.0_dp*t**2 - 8.0_dp*t**3 + 5.0_dp*t**4)
         radial_derivative_weights(left_pos) = &
            radial_derivative_weights(left_pos) + dh00/h_interval
         radial_derivative_weights(right_pos) = &
            radial_derivative_weights(right_pos) + dh01/h_interval
         radial_derivative_weights(1:nradial) = radial_derivative_weights(1:nradial) + &
                                                dh10*slope_left(1:nradial) + &
                                                dh11*slope_right(1:nradial) + &
                                                h_interval*(dh20*curvature_left(1:nradial) + &
                                                            dh21*curvature_right(1:nradial))
      END IF

      angular_values = 0.0_dp
      IF (PRESENT(angular_derivative_weights)) angular_value_derivatives = 0.0_dp
      DO iso = 1, harmonics%max_s_harm
         l = indso(1, iso)
         CALL y_lm(direction, angular_values(iso), l, indso(2, iso))
         IF (l == 0 .OR. .NOT. PRESENT(angular_derivative_weights)) CYCLE
         shell_index = iso - nsoset(l - 1)
         DO ic = 1, nco(l)
            lx = indco(1, ic + ncoset(l - 1))
            ly = indco(2, ic + ncoset(l - 1))
            lz = indco(3, ic + ncoset(l - 1))
            IF (lx > 0) THEN
               monomial = REAL(lx, dp)*direction(1)**(lx - 1)* &
                          direction(2)**ly*direction(3)**lz
               angular_value_derivatives(1, iso) = angular_value_derivatives(1, iso) + &
                                                   orbtramat(l)%slm(shell_index, ic)*monomial
            END IF
            IF (ly > 0) THEN
               monomial = direction(1)**lx*REAL(ly, dp)*direction(2)**(ly - 1)* &
                          direction(3)**lz
               angular_value_derivatives(2, iso) = angular_value_derivatives(2, iso) + &
                                                   orbtramat(l)%slm(shell_index, ic)*monomial
            END IF
            IF (lz > 0) THEN
               monomial = direction(1)**lx*direction(2)**ly* &
                          REAL(lz, dp)*direction(3)**(lz - 1)
               angular_value_derivatives(3, iso) = angular_value_derivatives(3, iso) + &
                                                   orbtramat(l)%slm(shell_index, ic)*monomial
            END IF
         END DO
         DO inode = 1, 3
            solid_derivative = angular_value_derivatives(inode, iso)
            angular_value_derivatives(inode, iso) = (solid_derivative - &
                                                     REAL(l, dp)*angular_values(iso)* &
                                                     direction(inode))/radius
         END DO
      END DO

      DO ia = 1, grid_atom%ng_sphere
         angular_weights(ia) = grid_atom%wa(ia)* &
                               DOT_PRODUCT(harmonics%slm(ia, 1:harmonics%max_s_harm), &
                                           angular_values)
         IF (PRESENT(angular_derivative_weights)) THEN
            DO inode = 1, 3
               angular_derivative_weights(inode, ia) = grid_atom%wa(ia)* &
                                                       DOT_PRODUCT(harmonics%slm(ia, 1:harmonics%max_s_harm), &
                                                                   angular_value_derivatives(inode, :))
            END DO
         END IF
      END DO
      active = .TRUE.

   END SUBROUTINE atom_grid_interpolation_weights
   SUBROUTINE interpolate_gapw_atom_grid_fields( &
      grid_atom, harmonics, displacement, cutoff, nspins, rho_h, rho_s, drho_h, drho_s, &
      tau_h, tau_s, density, gradient, kin, density_spatial, gradient_spatial, kin_spatial, calculate_spatial)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      REAL(dp), INTENT(IN)                               :: cutoff
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: rho_h, rho_s
      REAL(dp), DIMENSION(:, :, :, :), INTENT(IN)        :: drho_h, drho_s
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: tau_h, tau_s
      REAL(dp), DIMENSION(2), INTENT(OUT)                :: density
      REAL(dp), DIMENSION(3, 2), INTENT(OUT)             :: gradient
      REAL(dp), DIMENSION(2), INTENT(OUT)                :: kin
      REAL(dp), DIMENSION(3, 2), INTENT(OUT)             :: density_spatial
      REAL(dp), DIMENSION(3, 3, 2), INTENT(OUT)          :: gradient_spatial
      REAL(dp), DIMENSION(3, 2), INTENT(OUT)             :: kin_spatial
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_spatial

      INTEGER                                            :: ia, idir, inode, ir, ispin, jdir, nradial
      INTEGER, DIMENSION(4)                              :: radial_indices
      LOGICAL                                            :: active, need_spatial
      REAL(dp)                                           :: derivative_weight, weight
      REAL(dp), DIMENSION(grid_atom%ng_sphere)           :: angular_weights
      REAL(dp), DIMENSION(4)                             :: radial_derivative_weights, radial_weights
      REAL(dp), DIMENSION(3, grid_atom%ng_sphere)        :: angular_derivative_weights

      density = 0.0_dp
      gradient = 0.0_dp
      kin = 0.0_dp
      density_spatial = 0.0_dp
      gradient_spatial = 0.0_dp
      kin_spatial = 0.0_dp
      need_spatial = .TRUE.
      IF (PRESENT(calculate_spatial)) need_spatial = calculate_spatial
      IF (need_spatial) THEN
         CALL atom_grid_interpolation_weights( &
            grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
            radial_derivative_weights, nradial, angular_weights, angular_derivative_weights, active)
      ELSE
         CALL atom_grid_interpolation_weights( &
            grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
            nradial=nradial, angular_weights=angular_weights, active=active)
      END IF
      IF (active) THEN
         DO inode = 1, nradial
            ir = radial_indices(inode)
            DO ia = 1, grid_atom%ng_sphere
               weight = radial_weights(inode)*angular_weights(ia)
               DO ispin = 1, nspins
                  density(ispin) = density(ispin) + &
                                   weight*(rho_h(ia, ir, ispin) - rho_s(ia, ir, ispin))
                  kin(ispin) = kin(ispin) + &
                               weight*(tau_h(ia, ir, ispin) - tau_s(ia, ir, ispin))
                  DO idir = 1, 3
                     gradient(idir, ispin) = gradient(idir, ispin) + &
                                             weight*(drho_h(idir, ia, ir, ispin) - &
                                                     drho_s(idir, ia, ir, ispin))
                     IF (.NOT. need_spatial) CYCLE
                     derivative_weight = radial_derivative_weights(inode)* &
                                         displacement(idir)/SQRT(SUM(displacement**2))* &
                                         angular_weights(ia) + radial_weights(inode)* &
                                         angular_derivative_weights(idir, ia)
                     density_spatial(idir, ispin) = density_spatial(idir, ispin) + &
                                                    derivative_weight*(rho_h(ia, ir, ispin) - rho_s(ia, ir, ispin))
                     kin_spatial(idir, ispin) = kin_spatial(idir, ispin) + &
                                                derivative_weight*(tau_h(ia, ir, ispin) - tau_s(ia, ir, ispin))
                     DO jdir = 1, 3
                        gradient_spatial(jdir, idir, ispin) = &
                           gradient_spatial(jdir, idir, ispin) + derivative_weight*( &
                           drho_h(jdir, ia, ir, ispin) - drho_s(jdir, ia, ir, ispin))
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
   END SUBROUTINE interpolate_gapw_atom_grid_fields
   SUBROUTINE add_gapw_atom_grid_interpolation_adjoint( &
      grid_atom, harmonics, displacement, cutoff, nspins, density_adjoint, gradient_adjoint, &
      kin_adjoint, vxc_h, vxc_s, vxg_h, vxg_s, vtau_h, vtau_s)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      REAL(dp), INTENT(IN)                               :: cutoff
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: density_adjoint
      REAL(dp), DIMENSION(3, 2), INTENT(IN)              :: gradient_adjoint
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: kin_adjoint
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: vxc_h, vxc_s
      REAL(dp), DIMENSION(:, :, :, :), INTENT(INOUT)     :: vxg_h, vxg_s
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: vtau_h, vtau_s

      INTEGER                                            :: ia, idir, inode, ir, ispin, nradial
      INTEGER, DIMENSION(4)                              :: radial_indices
      LOGICAL                                            :: active
      REAL(dp)                                           :: value
      REAL(dp), DIMENSION(4)                             :: radial_weights
      REAL(dp), DIMENSION(grid_atom%ng_sphere)           :: angular_weights

      CALL atom_grid_interpolation_weights( &
         grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
         nradial=nradial, angular_weights=angular_weights, active=active)
      IF (active) THEN
         DO inode = 1, nradial
            ir = radial_indices(inode)
            DO ia = 1, grid_atom%ng_sphere
               value = radial_weights(inode)*angular_weights(ia)
               DO ispin = 1, nspins
                  ! CP2K applies the hard-minus-soft sign when the two one-center
                  ! matrices are assembled, so both stored potentials carry the
                  ! same transpose-interpolation coefficient.
                  vxc_h(ia, ir, ispin) = vxc_h(ia, ir, ispin) + value*density_adjoint(ispin)
                  vxc_s(ia, ir, ispin) = vxc_s(ia, ir, ispin) + value*density_adjoint(ispin)
                  vtau_h(ia, ir, ispin) = vtau_h(ia, ir, ispin) + value*kin_adjoint(ispin)
                  vtau_s(ia, ir, ispin) = vtau_s(ia, ir, ispin) + value*kin_adjoint(ispin)
                  DO idir = 1, 3
                     vxg_h(idir, ia, ir, ispin) = vxg_h(idir, ia, ir, ispin) + &
                                                  value*gradient_adjoint(idir, ispin)
                     vxg_s(idir, ia, ir, ispin) = vxg_s(idir, ia, ir, ispin) + &
                                                  value*gradient_adjoint(idir, ispin)
                  END DO
               END DO
            END DO
         END DO
      END IF
   END SUBROUTINE add_gapw_atom_grid_interpolation_adjoint
 SUBROUTINE project(cell, particle_set, grid_atom, harmonics, nspins, flags, &
   composite_grid_coords, composite_grid_atom, composite_local_atoms, &
   composite_atom_start, composite_atom_end, composite_density_grad, &
   composite_grad_grad, composite_kin_grad, cross_cutoff, &
   vxc_h, vxc_s, vxg_h, vxg_s, vtau_h, vtau_s)
 TYPE(cell_type), POINTER :: cell
 TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
 TYPE(grid_atom_type), POINTER :: grid_atom
 TYPE(harmonics_atom_type), POINTER :: harmonics
 INTEGER, INTENT(IN) :: nspins
 LOGICAL, INTENT(IN) :: flags(3)
 REAL(dp), INTENT(IN) :: composite_grid_coords(:, :), composite_density_grad(:, :)
 REAL(dp), INTENT(IN) :: composite_grad_grad(:, :, :), composite_kin_grad(:, :)
 INTEGER, INTENT(IN) :: composite_grid_atom(:), composite_local_atoms(:)
 INTEGER, INTENT(IN) :: composite_atom_start(:), composite_atom_end(:)
 REAL(dp), INTENT(IN) :: cross_cutoff
 REAL(dp), POINTER :: vxc_h(:, :, :), vxc_s(:, :, :), vtau_h(:, :, :), vtau_s(:, :, :)
 REAL(dp), POINTER :: vxg_h(:, :, :, :), vxg_s(:, :, :, :)
 REAL(dp), ALLOCATABLE :: vxc_h_local(:, :, :), vxc_s_local(:, :, :), vtau_h_local(:, :, :), vtau_s_local(:, :, :)
 REAL(dp), ALLOCATABLE :: vxg_h_local(:, :, :, :), vxg_s_local(:, :, :, :)
 INTEGER :: base_shift(3), composite_row, idir, image_i1, image_i2, image_i3
 INTEGER :: image_shift(3), image_shell(3), jdir, target_atom, iatom
 INTEGER :: composite_local_atom, composite_local_natom, composite_nflat
 REAL(dp) :: cross_density_adjoint(2), cross_grad_adjoint(3, 2), cross_kin_adjoint(2)
 REAL(dp) :: cross_displacement(3), fractional(3), image_translation(3)
 LOGICAL :: lsd, use_atom_composite_density, use_atom_composite_gradient, use_atom_composite_tau
 iatom = 1
 composite_nflat = SIZE(composite_grid_atom)
 composite_local_natom = SIZE(composite_local_atoms)
 lsd = nspins == 2
 use_atom_composite_density = flags(1)
 use_atom_composite_gradient = flags(2)
 use_atom_composite_tau = flags(3)
 image_shell = 0
 DO idir = 1, 3
   IF (cell%perd(idir) == 1) image_shell(idir) = &
     CEILING(cross_cutoff*SQRT(SUM(cell%h_inv(idir, :)**2))) + 1
 END DO
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(base_shift, composite_row, cross_density_adjoint, cross_displacement, &
!$OMP         cross_grad_adjoint, cross_kin_adjoint, fractional, idir, image_i1, &
!$OMP         image_i2, image_i3, image_shift, image_translation, jdir, target_atom, &
!$OMP         vxc_h_local, vxc_s_local, vxg_h_local, vxg_s_local, vtau_h_local, vtau_s_local) &
!$OMP SHARED(cell, composite_density_grad, composite_grad_grad, composite_grid_atom, &
!$OMP        composite_grid_coords, composite_kin_grad, composite_nflat, cross_cutoff, &
!$OMP        grid_atom, harmonics, iatom, image_shell, lsd, nspins, particle_set, &
!$OMP        use_atom_composite_density, use_atom_composite_gradient, use_atom_composite_tau, &
!$OMP        vxc_h, vxc_s, vxg_h, vxg_s, vtau_h, vtau_s)
                        ALLOCATE (vxc_h_local, MOLD=vxc_h)
                        ALLOCATE (vxc_s_local, MOLD=vxc_s)
                        ALLOCATE (vxg_h_local, MOLD=vxg_h)
                        ALLOCATE (vxg_s_local, MOLD=vxg_s)
                        ALLOCATE (vtau_h_local, MOLD=vtau_h)
                        ALLOCATE (vtau_s_local, MOLD=vtau_s)
                        vxc_h_local = 0.0_dp
                        vxc_s_local = 0.0_dp
                        vxg_h_local = 0.0_dp
                        vxg_s_local = 0.0_dp
                        vtau_h_local = 0.0_dp
                        vtau_s_local = 0.0_dp
!$OMP DO SCHEDULE(STATIC)
                        DO composite_row = 1, composite_nflat
                           target_atom = composite_grid_atom(composite_row)
                           cross_density_adjoint = 0.0_dp
                           cross_grad_adjoint = 0.0_dp
                           cross_kin_adjoint = 0.0_dp
                           IF (lsd) THEN
                              IF (use_atom_composite_density) THEN
                                 cross_density_adjoint(1:2) = &
                                    composite_density_grad(composite_row, 1:2)
                              END IF
                              IF (use_atom_composite_gradient) THEN
                                 cross_grad_adjoint(:, 1:2) = &
                                    composite_grad_grad(composite_row, :, 1:2)
                              END IF
                              IF (use_atom_composite_tau) THEN
                                 cross_kin_adjoint(1:2) = &
                                    composite_kin_grad(composite_row, 1:2)
                              END IF
                           ELSE
                              IF (use_atom_composite_density) THEN
                                 cross_density_adjoint(1) = 0.5_dp* &
                                                            SUM(composite_density_grad(composite_row, :))
                              END IF
                              IF (use_atom_composite_gradient) THEN
                                 DO idir = 1, 3
                                    cross_grad_adjoint(idir, 1) = 0.5_dp* &
                                                                  SUM(composite_grad_grad(composite_row, idir, :))
                                 END DO
                              END IF
                              IF (use_atom_composite_tau) THEN
                                 cross_kin_adjoint(1) = 0.5_dp* &
                                                        SUM(composite_kin_grad(composite_row, :))
                              END IF
                           END IF
                           fractional = 0.0_dp
                           DO idir = 1, 3
                              DO jdir = 1, 3
                                 fractional(idir) = fractional(idir) + cell%h_inv(idir, jdir)* &
                                                    (composite_grid_coords(jdir, composite_row) - &
                                                     particle_set(iatom)%r(jdir))
                              END DO
                           END DO
                           DO idir = 1, 3
                              base_shift(idir) = cell%perd(idir)*NINT(fractional(idir))
                           END DO
                           DO image_i3 = base_shift(3) - image_shell(3), &
                              base_shift(3) + image_shell(3)
                              DO image_i2 = base_shift(2) - image_shell(2), &
                                 base_shift(2) + image_shell(2)
                                 DO image_i1 = base_shift(1) - image_shell(1), &
                                    base_shift(1) + image_shell(1)
                                    image_shift = [image_i1, image_i2, image_i3]
                                    IF (target_atom == iatom .AND. ALL(image_shift == 0)) CYCLE
                                    image_translation = MATMUL( &
                                                        cell%hmat, REAL(image_shift, dp))
                                    cross_displacement = &
                                       composite_grid_coords(:, composite_row) - &
                                       particle_set(iatom)%r - image_translation
                                    ! Reject inactive images before allocating interpolation scratch.
                                    IF (SQRT(SUM(cross_displacement**2)) > cross_cutoff) CYCLE
                                    CALL add_gapw_atom_grid_interpolation_adjoint( &
                                       grid_atom, harmonics, cross_displacement, cross_cutoff, &
                                       nspins, cross_density_adjoint, cross_grad_adjoint, &
                                       cross_kin_adjoint, vxc_h_local, vxc_s_local, vxg_h_local, vxg_s_local, &
                                       vtau_h_local, vtau_s_local)
                                 END DO
                              END DO
                           END DO
                        END DO
!$OMP END DO
!$OMP CRITICAL(skala_atom_composite_adjoint_reduction)
                        vxc_h = vxc_h + vxc_h_local
                        vxc_s = vxc_s + vxc_s_local
                        vxg_h = vxg_h + vxg_h_local
                        vxg_s = vxg_s + vxg_s_local
                        vtau_h = vtau_h + vtau_h_local
                        vtau_s = vtau_s + vtau_s_local
!$OMP END CRITICAL(skala_atom_composite_adjoint_reduction)
                        DEALLOCATE (vxc_h_local, vxc_s_local, vxg_h_local, vxg_s_local, vtau_h_local, vtau_s_local)
!$OMP END PARALLEL
 END SUBROUTINE project

END MODULE parallel_kernel
