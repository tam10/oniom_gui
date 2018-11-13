module FMathMod
    implicit none
    integer, parameter:: dp=kind(0.d0)
    !integer*4 :: num_atoms, num_stretches, num_bends, &
    !num_torsions, num_non_bonding

    !integer*4 :: stretch_atom_nums(num_stretches * 2), &
    !bend_atom_nums(num_bends * 3), torsion_atom_nums(num_torsions * 4), &
    !non_bonding_atom_nums(num_non_bonding * 2)

    !real(dp) :: stretch_es_tot(2), bend_es_tot(2), &
    !torsions_es_tot(2), vdw_es_tot(2), coul_es_tot(2)

    !real(dp) :: stretch_reqs(num_stretches), stretch_keqs(num_stretches), &
    !bend_aeqs(num_bends), bend_keqs(num_bends), &
    !torsions_vs(num_torsions * 4), torsions_gammas(num_torsions * 4), &
    !torsion_paths(num_torsions), &
    !vdw_vs(num_non_bonding), vdw_reqs(num_non_bonding), q0q1s(num_non_bonding)

    

    contains

    !Move parameters to here and initialise

    subroutine e_amber(Coords, Forces, num_atoms, &
        num_stretches, stretch_atom_nums, &
        stretch_es_tot, stretch_reqs, stretch_keqs, &
        num_bends, bend_atom_nums, &
        bend_es_tot, bend_aeqs, bend_keqs, &
        num_torsions, torsion_atom_nums, &
        torsions_es_tot, torsions_vs, torsions_gammas, torsion_paths, &
        num_non_bonding, non_bonding_atom_nums, &
        vdw_es_tot, vdw_vs, vdw_reqs, vdw_cutoff, &
        coul_es_tot, q0q1s, diel, coul_cutoff, &
        abs_gradient, crash)
        real(dp), intent(in) :: Coords(3, num_atoms)
        real(dp), intent(out) :: Forces(3, num_atoms)
        integer*4, intent(in) :: num_atoms, num_stretches, &
        num_bends, num_torsions, num_non_bonding
        integer*4, intent(in) :: stretch_atom_nums(num_stretches * 2), &
        bend_atom_nums(num_bends * 3), torsion_atom_nums(num_torsions * 4), &
        non_bonding_atom_nums(num_non_bonding * 2)
        real(dp), intent(out) :: stretch_es_tot(2), bend_es_tot(2), &
        torsions_es_tot(2), vdw_es_tot(2), coul_es_tot(2)
        real(dp), intent(in) :: stretch_reqs(num_stretches), &
        stretch_keqs(num_stretches), bend_aeqs(num_bends), &
        bend_keqs(num_bends), torsions_vs(num_torsions * 4), & 
        torsions_gammas(num_torsions * 4), torsion_paths(num_torsions), &
        vdw_vs(num_non_bonding), vdw_reqs(num_non_bonding), q0q1s(num_non_bonding)
        real(dp), intent(in) :: vdw_cutoff, diel, coul_cutoff
        real(dp), intent(out) :: abs_gradient
        integer*4, intent(inout) :: crash

        Forces = 0.0d0

        call e_stretches(Coords, Forces, num_atoms, &
        num_stretches, stretch_atom_nums, &
        stretch_es_tot, stretch_reqs, stretch_keqs, &
        abs_gradient)

        call e_bends(Coords, Forces, num_atoms, &
        num_bends, bend_atom_nums, &
        bend_es_tot, bend_aeqs, bend_keqs, &
        abs_gradient, crash)

        call e_torsions(Coords, Forces, num_atoms, &
        num_torsions, torsion_atom_nums, &
        torsions_es_tot, torsions_vs, torsions_gammas, &
        torsion_paths, abs_gradient, crash)

        call e_nonbonding(Coords, Forces, num_atoms, &
        num_non_bonding, non_bonding_atom_nums, &
        vdw_es_tot, vdw_vs, vdw_reqs, vdw_cutoff, &
        coul_es_tot, q0q1s, diel, coul_cutoff, &
        abs_gradient, crash)

    end subroutine

    subroutine e_stretches(Coords, Forces, num_atoms, &
        num_stretches, stretch_atom_nums, &
        stretch_es_tot, stretch_reqs, stretch_keqs, &
        abs_gradient)
        real(dp), intent(in) :: Coords(3, num_atoms)
        real(dp), intent(inout) :: Forces(3, num_atoms)
        integer*4, intent(in) :: num_atoms, num_stretches
        integer*4, intent(in) :: stretch_atom_nums(num_stretches * 2)
        real(dp), intent(out) :: stretch_es_tot(2)
        real(dp), intent(in) :: stretch_reqs(num_stretches), &
        stretch_keqs(num_stretches)
        real(dp), intent(out) :: abs_gradient

        integer*4 :: stretch_num, j, a0, a1
        real(dp) :: v10(3), force(3)
        real(dp) :: stretch_es(2)
        real(dp) :: r10, dr, v

        integer*4, parameter :: int_one=1, int_two=2, int_three=3
        real(dp), parameter :: two=2.d0

        do stretch_num = 1, num_stretches
            
            a0 = stretch_atom_nums(stretch_num * int_two - int_one) + int_one
            a1 = stretch_atom_nums(stretch_num * int_two) + int_one

            v = stretch_keqs(stretch_num)
            if (v.eq.0) then
                cycle
            end if
            
            !Get positions, vector and distance
            call get_r(Coords, num_atoms, a1, a0, v10, r10)
            
            dr = r10 - stretch_reqs(stretch_num)
            
            stretch_es(1) = dr  * dr * v
            stretch_es(2) = two * dr * v
            
            stretch_es_tot(1) = stretch_es_tot(1) + stretch_es(1)
            stretch_es_tot(2) = stretch_es_tot(2) + stretch_es(2)
            abs_gradient = abs_gradient + abs(stretch_es(2))
            
            do j = 1, 3
                force(j) = v10(j) * stretch_es(2) / r10
                Forces(j,a0) = Forces(j,a0) - force(j)
                Forces(j,a1) = Forces(j,a1) + force(j)
            end do
            
        end do

    end subroutine e_stretches

    subroutine e_bends(Coords, Forces, num_atoms, &
        num_bends, bend_atom_nums, &
        bend_es_tot, bend_aeqs, bend_keqs, &
        abs_gradient, crash)
        real(dp), intent(in) :: Coords(3, num_atoms)
        real(dp), intent(inout) :: Forces(3, num_atoms)
        integer*4, intent(in) :: num_atoms, num_bends
        integer*4, intent(in) :: bend_atom_nums(num_bends * 2)
        real(dp), intent(out) :: bend_es_tot(2)
        real(dp), intent(in) :: bend_aeqs(num_bends), &
        bend_keqs(num_bends)
        real(dp), intent(out) :: abs_gradient
        integer*4, intent(out) :: crash

        integer*4 :: bend_num, j, a0, a1, a2, bend_index
        real(dp) :: v01n(3), v21n(3), force0(3), force2(3), w012(3)
        real(dp) :: bend_es(2)
        real(dp) :: r01, r21, v, a012, da

        integer*4, parameter :: int_one=1, int_two=2, int_three=3
        real(dp), parameter :: two=2.d0, almost_pi=3.95d0*atan(1.d0)

        crash = 0

        do bend_num = 1, num_bends
            bend_index = bend_num * int_three

            v = bend_keqs(bend_num)
            if (v.eq.0) then
                cycle
            end if

            a0 = bend_atom_nums(bend_index - int_two) + int_one
            a1 = bend_atom_nums(bend_index - int_one) + int_one
            a2 = bend_atom_nums(bend_index) + int_one

            call get_a(Coords, num_atoms, a0, a1, a2, &
            v01n, v21n, r01, r21, a012, w012)

            if (a012.ge.almost_pi) then
                crash = 1
                cycle
            end if
            
            da = a012 - bend_aeqs(bend_num)

            bend_es(1) =  da * da * v
            bend_es(2) = two * da * v

            bend_es_tot(1) = bend_es_tot(1) + bend_es(1)
            bend_es_tot(2) = bend_es_tot(2) + bend_es(2)
            abs_gradient = abs_gradient + abs(bend_es(2))

            force0 = cross_product(v01n, w012)
            force2 = cross_product(v21n, w012)
            
            call multiply_vec(force0, bend_es(2) / r01)
            call multiply_vec(force2, - bend_es(2) / r21)

            do j = 1, 3
                Forces(j, a0) = Forces(j, a0) + force0(j)
                Forces(j, a1) = Forces(j, a1) - force0(j) - force2(j)
                Forces(j, a2) = Forces(j, a2) + force2(j)
            end do

        end do

    end subroutine e_bends

    subroutine e_torsions(Coords, Forces, num_atoms, &
        num_torsions, torsion_atom_nums, &
        torsions_es_tot, torsions_vs, torsions_gammas, torsion_paths, &
        abs_gradient, crash)
        real(dp), intent(in) :: Coords(3, num_atoms)
        real(dp), intent(inout) :: Forces(3, num_atoms)
        integer*4, intent(in) :: num_atoms, num_torsions
        integer*4, intent(in) :: torsion_atom_nums(num_torsions * 2)
        real(dp), intent(out) :: torsions_es_tot(2)
        real(dp), intent(out) :: abs_gradient
        real(dp), intent(in) :: torsions_vs(num_torsions * 4), &
        torsions_gammas(num_torsions * 4), torsion_paths(num_torsions)
        integer*4, intent(out) :: crash

        integer*4 :: torsion_num, torsion_index, periodicity, j, &
        a0, a1, a2, a3
        real(dp) :: v01n(3), v21n(3), v23n(3), w012n(3), w123n(3), &
        force0(3), force1(3), force2(3), force3(3), torque3(3), vo3(3)
        real(dp) :: torsion_es(2)
        real(dp) :: r01, r21, r23, a012, a123, d0123, dd, v
        real(dp) :: e_da_dr10, e_da_dr23, ro3ro3

        integer*4, parameter :: int_one=1, int_two=2, int_three=3, int_four=4
        real(dp), parameter :: zero=0.d0, half=0.5d0, one=1.0d0, two=2.0d0

        !Reference:
		!http://xray.bmc.uu.se/~aqwww/q_legacy/documents/qman5.pdf

        do torsion_num = 1, num_torsions

            torsion_index = (torsion_num - int_one) * int_four
            
            a0 = torsion_atom_nums(torsion_index + int_one) + int_one
            a1 = torsion_atom_nums(torsion_index + int_two) + int_one
            a2 = torsion_atom_nums(torsion_index + int_three) + int_one
            a3 = torsion_atom_nums(torsion_index + int_four) + int_one

            !Get dihedral, angles, vectors, norms and distances
            call get_d(Coords, num_atoms, a0, a1, a2, a3, &
            v01n, v21n, v23n, r01, r21, r23, a012, a123, w012n, w123n, d0123)

            !Compute energy and gradient for dihedral
            torsion_es(1) = zero
            torsion_es(2) = zero
            do periodicity = 1, 4
                v = half * torsions_vs(torsion_index + periodicity) / &
                (torsion_paths(torsion_index + int_one))

                dd = periodicity * d0123 - &
                torsions_gammas(torsion_index + periodicity)

                torsion_es(1) = torsion_es(1) + v * (one + cos(dd))
                torsion_es(2) = torsion_es(2) - v * sin(dd)
            end do

            torsions_es_tot(1) = torsions_es_tot(1) + torsion_es(1)
            torsions_es_tot(2) = torsions_es_tot(2) + torsion_es(2)
            abs_gradient = abs_gradient + abs(torsion_es(2))

            if (a012.eq.zero.or.a123.eq.zero) then
                crash = 1
                cycle
            end if

            !Energy times da/dr to convert to cartesian
            e_da_dr10 = torsion_es(1) / (r01 * sin(a012))
            e_da_dr23 = torsion_es(1) / (r23 * sin(a123))

            !Compute forces on outside atoms and 
            ! distance from centre of 1-2 to 3 squared
            ro3ro3 = 0
            do j = 1, 3
                force0(j) = - w012n(j) * e_da_dr10
                force3(j) = - w123n(j) * e_da_dr23

                vo3(j) = Coords(j, a3) - &
                    (Coords(j, a1) + Coords(j, a2)) * half

                ro3ro3 = ro3ro3 + vo3(j) ** 2
            end do

            torque3 = (cross_product(vo3, force3) * two + &
                cross_product(v23n, force3) * r23 + &
                cross_product(v21n, force0) * r21)

            force2 = cross_product(torque3, vo3) / &
                ( - ro3ro3 * two)

            do j = 1, 3
                force1(j) = - force0(j) - force2(j) - force3(j)

                Forces(j,a0) = Forces(j,a0) + force0(j)
                Forces(j,a1) = Forces(j,a1) + force1(j)
                Forces(j,a2) = Forces(j,a2) + force2(j)
                Forces(j,a3) = Forces(j,a3) + force3(j)

            end do 

        end do

    end subroutine e_torsions

    subroutine e_nonbonding(Coords, Forces, num_atoms, &
        num_non_bonding, non_bonding_atom_nums, &
        vdw_es_tot, vdw_vs, vdw_reqs, vdw_cutoff, &
        coul_es_tot, q0q1s, diel, coul_cutoff, &
        abs_gradient, crash)
        real(dp), intent(in) :: Coords(3, num_atoms)
        real(dp), intent(inout) :: Forces(3, num_atoms)
        integer*4, intent(in) :: num_atoms, num_non_bonding
        integer*4, intent(in) :: non_bonding_atom_nums(num_non_bonding * 2)
        real(dp), intent(out) :: vdw_es_tot(2), coul_es_tot(2)
        real(dp), intent(out) :: abs_gradient
        real(dp), intent(in) :: vdw_vs(num_non_bonding), &
        vdw_reqs(num_non_bonding), q0q1s(num_non_bonding)
        real(dp), intent(in) :: vdw_cutoff, coul_cutoff
        real(dp), intent(in) :: diel
        integer*4, intent(out) :: crash

        integer*4 :: non_bonding_num, j, a0, a1
        real(dp) :: v10(3), force(3)
        real(dp) :: vdw_es(2), coul_es(2)
        real(dp) :: r10, r1, r2, r6, r12, v

        integer*4, parameter :: int_one=1, int_two=2
        real(dp), parameter :: two=2.d0, twelve=12.0d0, one=1.0d0

        vdw_es(1) = 0
        vdw_es(2) = 0
        coul_es(1) = 0
        coul_es(2) = 0
        crash = 0
        abs_gradient = 0

        do non_bonding_num = 1, num_non_bonding
            a0 = non_bonding_atom_nums(non_bonding_num * int_two - int_one) + int_one
            a1 = non_bonding_atom_nums(non_bonding_num * int_two) + int_one

            !Get positions, vector and distance
            call get_r(Coords, num_atoms, a1, a0, v10, r10)

            if (r10.eq.0) then
                crash = 1
                cycle
            end if

            !Compute VdW
            v = vdw_vs(non_bonding_num)
            if (r10.le.vdw_cutoff.and.v.gt.0) then
                r1 = vdw_reqs(non_bonding_num) / r10
                r2 = r1 * r1
                r6 = r2 * r2 * r2
                r12 = v * r6 * r6
                r6 = v * r6

                vdw_es(1) = r12 - two * r6
                r12 = twelve * r12 * r1 / r10
                r6 = twelve * r6 * r1 / r10
                vdw_es(2) = r6 - r12
                !r6 = - 7 * r6 * r1 / r10
                !r12 = 13 * r12 / r10
                !vdw_es(3) = vdw_es(3) + r12 + r6

                vdw_es_tot(1) = vdw_es_tot(1) + vdw_es(1)
                vdw_es_tot(2) = vdw_es_tot(2) + vdw_es(2)
                abs_gradient = abs_gradient + abs(vdw_es(2))


            else if (r10.ge.coul_cutoff) then
                !skip if outside of both vdw_cutoff and coul_cutoff
                cycle
            end if

            !Compute electrostatic
            v = q0q1s(non_bonding_num)
            if (v.ne.0.and.diel.gt.0) then
                r1 = 1 / r10
                coul_es(1) = v * r1 / diel
                coul_es(2) = - coul_es(1) * r1
                !coul_es(3) = - 2 * coul_es(2) * r1

                coul_es_tot(1) = coul_es_tot(1) + coul_es(1)
                coul_es_tot(2) = coul_es_tot(2) + coul_es(2)
                abs_gradient = abs_gradient + abs(coul_es(2))
            end if

            do j = 1, 3
                force(j) = v10(j) * (vdw_es(2) + coul_es(2)) / r10
                Forces(j,a0) = Forces(j,a0) + force(j)
                Forces(j,a1) = Forces(j,a1) - force(j)
            end do

        end do

    end subroutine

    subroutine get_r(Coords, num_atoms, a0, a1, v01, r01)
        real(dp), intent(in) :: Coords(3, num_atoms)
        integer*4, intent(in) :: num_atoms, a0, a1
        real(dp), intent(out) :: v01(3)
        real(dp), intent(out) :: r01

        real(dp) :: r01s
        integer*4 :: j

        r01s = 0
        do j = 1, 3
            v01(j) = Coords(j, a1) - Coords(j, a0)
            r01s = r01s + v01(j) * v01(j)
        end do
        r01 = sqrt(r01s)
        
    end subroutine get_r

    subroutine get_a(Coords, num_atoms, a0, a1, a2, &
        v01n, v21n, r01, r21, a012, w012n)
        real(dp), intent(in) :: Coords(3, num_atoms)
        integer*4, intent(in) :: num_atoms, a0, a1, a2
        real(dp), intent(out) :: v01n(3), v21n(3)
        real(dp), intent(out) :: r01, r21
        real(dp), intent(out) :: a012
        real(dp), intent(out) :: w012n(3)

        call get_r(Coords, num_atoms, a0, a1, v01n, r01)
        call get_r(Coords, num_atoms, a2, a1, v21n, r21)

        call divide_vec(v01n, r01)
        call divide_vec(v21n, r21)
        w012n = cross_product(v01n, v21n)
        a012 = unsigned_angle(v01n, v21n)
        
    end subroutine get_a

    subroutine get_d(Coords, num_atoms, a0, a1, a2, a3, &
        v01n, v21n, v23n, r01, r21, r23, a012, a123, w012n, w123n, d0123)
        real(dp), intent(in) :: Coords(3, num_atoms)
        integer*4, intent(in) :: num_atoms, a0, a1, a2, a3
        real(dp), intent(out) :: v01n(3), v21n(3), v23n(3)
        real(dp), intent(out) :: r01, r21, r23
        real(dp), intent(out) :: a012, a123
        real(dp), intent(out) :: w012n(3), w123n(3)
        real(dp), intent(out) :: d0123

        real(dp) :: c0123(3)

        call get_a(Coords, num_atoms, a0, a1, a2, &
        v01n, v21n, r01, r21, a012, w012n)
        call get_a(Coords, num_atoms, a1, a2, a3, &
        v21n, v23n, r21, r23, a123, w123n)

        d0123 = unsigned_angle(w012n, w123n)
        c0123 = cross_product(w012n, w123n)

        if (dot_product(v21n, c0123) < 0) then
            d0123 = - d0123
        end if

    end subroutine get_d

    subroutine divide_vec(vec, r) 
        real(dp), intent(inout) :: vec(3)
        real(dp), intent(in) :: r

        real(dp) :: r_inv

        r_inv = 1 / r
        vec(1) = vec(1) * r_inv
        vec(2) = vec(2) * r_inv
        vec(3) = vec(3) * r_inv

    end subroutine divide_vec

    subroutine multiply_vec(vec, r) 
        real(dp), intent(inout) :: vec(3)
        real(dp), intent(in) :: r

        vec(1) = vec(1) * r
        vec(2) = vec(2) * r
        vec(3) = vec(3) * r

    end subroutine multiply_vec

    function cross_product(v0, v1)
        real(dp), intent(in) :: v0(3), v1(3)
        real(dp) :: cross_product(3)
        
        cross_product(1) = v0(2) * v1(3) - v0(3) * v1(2)
        cross_product(2) = v0(3) * v1(1) - v0(1) * v1(3)
        cross_product(3) = v0(1) * v1(2) - v0(2) * v1(1)
    end function cross_product

    function unsigned_angle(v0, v1)
        real(dp), intent(in) :: v0(3), v1(3)
        real(dp) :: unsigned_angle

        unsigned_angle = acos(dot_product(v0, v1))
    end function unsigned_angle


    subroutine get_best_camera_view(Coords, num_atoms, &
        camera_position, field_of_view_rad, fill_amount, &
        min_distance)
        real(dp), intent(in) :: Coords(3, num_atoms)
        real(dp), intent(out) :: camera_position(3)
        integer*4, intent(in) :: num_atoms
        real(dp), intent(in) :: field_of_view_rad, fill_amount, min_distance

        real(dp) :: centre(3)
        real(dp) :: squared_radius, largest_squared_radius, radius
        integer*4 :: atom_num
        real(dp) :: ratio
        
        call get_centre(Coords, num_atoms, centre)
        !Align along x and y axes
        camera_position(1) = centre(1)
        camera_position(2) = centre(2)

        !Get largest distance from aligned axis
        largest_squared_radius = 0
        do atom_num = 1, num_atoms
            squared_radius = &
            (Coords(1, atom_num) - centre(1)) ** 2 + &
            (Coords(2, atom_num) - centre(2)) ** 2

            if (squared_radius.gt.largest_squared_radius) then
                largest_squared_radius = squared_radius
            end if
        end do

        radius = sqrt(largest_squared_radius) / fill_amount
        ratio = tan(field_of_view_rad / 2)
        camera_position(3) = min(-radius / ratio, min_distance)

    end subroutine get_best_camera_view

    subroutine get_centre(Coords, num_atoms, centre)
        real(dp), intent(in) :: Coords(3, num_atoms)
        real(dp), intent(out) :: centre(3)
        integer*4, intent(in) :: num_atoms

        integer*4 :: atom_num, j

        do j = 1, 3
            centre(j) = 0
        end do

        do atom_num = 1, num_atoms
            do j = 1, 3
                centre(j) = centre(j) + Coords(j, atom_num)
            end do
        end do

        do j = 1, 3
            centre(j) = centre(j) / num_atoms
        end do

    end subroutine get_centre

    subroutine e_vdw(r, v, req, energies)
      real(dp), intent(in) :: r, v, req
      real(dp), intent(out) :: energies(3)
      real(dp) :: r1, r2, r6, r12

    r1 = req / r
    r2 = r1 * r1
    r6 = r2 * r2 * r2
    r12 = r6 * r6

    !energies(1) = v * (r12 - 2 * r6)
    !energies(2) = 12 * v * (- r12 * r1 + r6 * r1) / r
    !energies(3) = 12 * v * (13 * r12 * r1 - 7 * r6 * r2) / (r * r)

    r12 = r12 * v
    r6 = r6 * v
    energies(1) = r12 - 2 * r6
    r12 = 12 * r12 * r1 / r
    r6 = 12 * r6 * r1 / r
    energies(2) = r6 - r12
    r6 = - 7 * r6 * r1 / r
    r12 = 13 * r12 / r
    energies(3) = r12 + r6

    end subroutine e_vdw

end module
