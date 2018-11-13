module AmberFortMod
    implicit none
    integer, parameter:: dp=kind(0.d0)
    integer*4 :: m_num_atoms, m_num_stretches=0, m_num_bends=0, &
    m_num_torsions=0, m_num_non_bonding=0

    integer*4, allocatable :: m_stretch_atom_nums(:), &
    m_bend_atom_nums(:), m_torsion_atom_nums(:), &
    m_non_bonding_atom_nums(:)

    real(dp), allocatable :: m_stretch_reqs(:), m_stretch_keqs(:), &
    m_bend_aeqs(:), m_bend_keqs(:), &
    m_torsions_vs(:), m_torsions_gammas(:), m_torsion_paths(:), &
    m_vdw_vs(:), m_vdw_reqs(:), m_q0q1s(:), m_masses(:)

    real(dp), allocatable :: m_Coords(:,:), m_Forces(:,:), & 
    m_Acceleration(:,:), m_Velocity(:,:)

    real(dp) :: m_vdw_cutoff, m_diel, m_coul_cutoff

    integer*4 :: bad_a0, bad_a1, bad_a2, bad_a3

    enum, bind(c)
        enumerator :: OK=0, ANGLE_EXC=1, STRETCH_EXC=2, BEND_EXC=3, &
        DIHEDRAL_EXC=4, NONBON_EXC1=5, NONBON_EXC2=6, NOATOMS_EXC=7, &
        NOMASS_EXC=8
    endenum

    !Bools to check if things are allocated or not.
    !Supposedly safer than allocated()
    logical :: m_Coords_allocated=.false.
    logical :: m_Forces_allocated=.false.
    logical :: m_Acceleration_allocated=.false.
    logical :: m_Velocity_allocated=.false.
    logical :: m_masses_allocated=.false.

    logical :: m_stretch_atom_nums_allocated=.false.
    logical :: m_stretch_reqs_allocated=.false.
    logical :: m_stretch_keqs_allocated=.false.

    logical :: m_bend_atom_nums_allocated=.false.
    logical :: m_bend_aeqs_allocated=.false.
    logical :: m_bend_keqs_allocated=.false.

    logical :: m_torsion_atom_nums_allocated=.false.
    logical :: m_torsions_vs_allocated=.false.
    logical :: m_torsions_gammas_allocated=.false.
    logical :: m_torsion_paths_allocated=.false.

    logical :: m_non_bonding_atom_nums_allocated=.false.
    logical :: m_vdw_vs_allocated=.false.
    logical :: m_vdw_reqs_allocated=.false.
    logical :: m_q0q1s_allocated=.false.


    integer*4, parameter :: int_one=1, int_two=2, int_three=3, int_four=4
    real(dp), parameter :: zero=0.d0, half=0.5d0, one=1.d0, &
    two=2.d0, almost_pi=3.95d0*atan(1.d0), twelve=12.d0

    contains

    subroutine reallocate_i1(A, is_allocated, dim1)
        integer*4, intent(in) :: dim1
        integer*4, allocatable, intent(inout) :: A(:)
        logical, intent(inout) :: is_allocated

        if (is_allocated) then
            deallocate(A)
        end if
        allocate(A(dim1))
        is_allocated = .true.

    end subroutine reallocate_i1

    subroutine reallocate_r1(A, is_allocated, dim1)
        integer*4, intent(in) :: dim1
        real(dp), allocatable, intent(inout) :: A(:)
        logical, intent(inout) :: is_allocated

        if (is_allocated) then
            deallocate(A)
        end if
        allocate(A(dim1))
        is_allocated = .true.

    end subroutine reallocate_r1

    subroutine reallocate_r2(A, is_allocated, dim1, dim2)
        integer*4, intent(in) :: dim1, dim2
        real(dp), allocatable, intent(inout) :: A(:, :)
        logical, intent(inout) :: is_allocated

        if (is_allocated) then
            deallocate(A)
        end if
        allocate(A(dim1, dim2))
        is_allocated = .true.

    end subroutine reallocate_r2

    subroutine get_bad_atoms(a0, a1, a2, a3)
        integer*4, intent(out) :: a0, a1, a2, a3
        a0 = bad_a0
        a1 = bad_a1
        a2 = bad_a2
        a3 = bad_a3

    end subroutine get_bad_atoms

    subroutine set_geometry(Coords, num_atoms)
        integer*4, intent(in) :: num_atoms
        real(dp), intent(in) :: Coords(3, num_atoms)

        m_num_atoms = num_atoms
        bad_a0 = -1
        bad_a1 = -1
        bad_a2 = -1
        bad_a3 = -1

        call reallocate_r2(m_Coords, m_Coords_allocated, 3, num_atoms)
        m_Coords(:,:) = Coords(:,:)

        call reallocate_r2(m_Forces, m_Forces_allocated, 3, num_atoms)
        m_Forces = 0.d0

        call reallocate_r2(m_Acceleration, m_Acceleration_allocated, 3, num_atoms)
        m_Acceleration = 0.d0

        call reallocate_r2(m_Velocity, m_Velocity_allocated, 3, num_atoms)
        m_Velocity = 0.d0

    end subroutine set_geometry
        
    subroutine set_masses(masses, status)
        real(dp), intent(in) :: masses(m_num_atoms)
        integer*4, intent(out) :: status

        integer*4 :: atom_num
        real(dp) :: mass

        if (m_num_atoms.eq.0) then
            status = NOATOMS_EXC
            return
        end if

        call reallocate_r1(m_masses, &
            m_masses_allocated, m_num_atoms)

        do atom_num = 1, m_num_atoms
            mass = masses(atom_num)
            if(mass.le.zero) then
                status = NOMASS_EXC
                bad_a0 = atom_num - 1
                return
            end if
            m_masses(atom_num) = mass
        end do

    end subroutine set_masses

    subroutine set_stretches(num_stretches, stretch_atom_nums, &
        stretch_reqs, stretch_keqs)
        integer*4, intent(in) :: num_stretches
        integer*4, intent(in) :: stretch_atom_nums(num_stretches * 2)
        real(dp), intent(in) :: stretch_reqs(num_stretches), &
        stretch_keqs(num_stretches)

        m_num_stretches = num_stretches

        call reallocate_i1(m_stretch_atom_nums, &
            m_stretch_atom_nums_allocated, num_stretches * 2)
        m_stretch_atom_nums(:) = stretch_atom_nums(:)

        call reallocate_r1(m_stretch_reqs, &
            m_stretch_reqs_allocated, num_stretches)
        m_stretch_reqs(:) = stretch_reqs(:)

        call reallocate_r1(m_stretch_keqs, &
            m_stretch_keqs_allocated, num_stretches)
        m_stretch_keqs(:) = stretch_keqs(:)
    end subroutine set_stretches

    subroutine set_bends(num_bends, bend_atom_nums, &
        bend_aeqs, bend_keqs)
        integer*4, intent(in) :: num_bends
        integer*4, intent(in) :: bend_atom_nums(num_bends * 3)
        real(dp), intent(in) :: bend_aeqs(num_bends), &
        bend_keqs(num_bends)

        m_num_bends = num_bends

        call reallocate_i1(m_bend_atom_nums, &
            m_bend_atom_nums_allocated, num_bends * 3)
        m_bend_atom_nums(:) = bend_atom_nums(:)

        call reallocate_r1(m_bend_aeqs, &
            m_bend_aeqs_allocated, num_bends)
        m_bend_aeqs(:) = bend_aeqs(:)

        call reallocate_r1(m_bend_keqs, & 
            m_bend_keqs_allocated, num_bends)
        m_bend_keqs(:) = bend_keqs(:)
    end subroutine set_bends

    subroutine set_torsions(num_torsions, torsion_atom_nums, &
        torsions_vs, torsions_gammas, torsion_paths)
        integer*4, intent(in) :: num_torsions
        integer*4, intent(in) :: torsion_atom_nums(num_torsions * 4)
        real(dp), intent(in) :: torsions_vs(num_torsions * 4), & 
        torsions_gammas(num_torsions * 4), torsion_paths(num_torsions)
        
        m_num_torsions = num_torsions

        call reallocate_i1(m_torsion_atom_nums, &
            m_torsion_atom_nums_allocated, num_torsions * 4)
        m_torsion_atom_nums(:) = torsion_atom_nums(:)

        call reallocate_r1(m_torsions_vs, &
            m_torsions_vs_allocated, num_torsions * 4)
        m_torsions_vs(:) = torsions_vs(:)

        call reallocate_r1(m_torsions_gammas, &
            m_torsions_gammas_allocated, num_torsions * 4)
        m_torsions_gammas(:) = torsions_gammas(:)

        call reallocate_r1(m_torsion_paths, &
            m_torsion_paths_allocated, num_torsions)
        m_torsion_paths(:) = torsion_paths(:)
    end subroutine set_torsions

    subroutine set_non_bonding(num_non_bonding, non_bonding_atom_nums, &
        vdw_vs, vdw_reqs, vdw_cutoff, q0q1s, diel, coul_cutoff)
        integer*4, intent(in) :: num_non_bonding
        integer*4, intent(in) :: non_bonding_atom_nums(num_non_bonding * 2)
        real(dp), intent(in) :: vdw_vs(num_non_bonding), &
        vdw_reqs(num_non_bonding), q0q1s(num_non_bonding)
        real(dp), intent(in) :: vdw_cutoff, diel, coul_cutoff

        m_num_non_bonding = num_non_bonding

        call reallocate_i1(m_non_bonding_atom_nums, &
            m_non_bonding_atom_nums_allocated, num_non_bonding * 2)
        m_non_bonding_atom_nums(:) = non_bonding_atom_nums(:)

        call reallocate_r1(m_vdw_vs, &
            m_vdw_vs_allocated, num_non_bonding)
        m_vdw_vs(:) = vdw_vs(:)

        call reallocate_r1(m_vdw_reqs, &
            m_vdw_reqs_allocated, num_non_bonding)
        m_vdw_reqs(:) = vdw_reqs(:)

        call reallocate_r1(m_q0q1s, &
            m_q0q1s_allocated, num_non_bonding)
        m_q0q1s(:) = q0q1s(:)

        m_vdw_cutoff = vdw_cutoff
        m_diel = diel
        m_coul_cutoff = coul_cutoff

    end subroutine set_non_bonding

    subroutine e_amber(Forces, &
        stretch_es_tot, bend_es_tot, &
        torsions_es_tot, vdw_es_tot, coul_es_tot, &
        abs_gradient, status)
        real(dp), intent(inout) :: Forces(3, m_num_atoms)
        real(dp), intent(out) :: stretch_es_tot(2), bend_es_tot(2), &
        torsions_es_tot(2), vdw_es_tot(2), coul_es_tot(2)
        real(dp), intent(out) :: abs_gradient
        integer*4, intent(inout) :: status

        Forces = 0.0d0

        call e_stretches(Forces, stretch_es_tot, abs_gradient, status)

        call e_bends(Forces, bend_es_tot, abs_gradient, status)
        if (status.ne.OK) then
            return
        end if
        call e_torsions(Forces, torsions_es_tot, abs_gradient, status)
        if (status.ne.OK) then
            return
        end if
        call e_nonbonding(Forces, vdw_es_tot, coul_es_tot, abs_gradient, status)
        if (status.ne.OK) then
            return
        end if

    end subroutine

    subroutine md_verlet(Coords, steps, dt, damping_factor, &
        stretch_es_tot, bend_es_tot, &
        torsions_es_tot, vdw_es_tot, coul_es_tot, &
        kinetic_energy, abs_gradient, status)
        real(dp), intent(inout) :: Coords(3, m_num_atoms)
        real(dp), intent(out) :: stretch_es_tot(2), bend_es_tot(2), &
        torsions_es_tot(2), vdw_es_tot(2), coul_es_tot(2)
        real(dp), intent(out) :: kinetic_energy, abs_gradient, &
        dt, damping_factor
        integer*4, intent(inout) :: steps, status

        integer*4 :: step_num, atom_num, j

        m_Coords(:,:) = Coords(:,:)
        stretch_es_tot = 0
        bend_es_tot = 0
        torsions_es_tot = 0
        vdw_es_tot = 0
        coul_es_tot = 0
        kinetic_energy = 0
        abs_gradient = 0

        do step_num = 1, steps
            call e_amber(m_Forces, stretch_es_tot, bend_es_tot, &
                torsions_es_tot, vdw_es_tot, coul_es_tot, &
                abs_gradient, status)

            do j = 1, 3
                m_Acceleration(:,j) = m_Forces(:,j) / m_masses(:)
            end do

            !Velocity Verlet update
            m_Coords = m_Coords + m_Velocity * dt + &
                half * (dt ** 2) * m_Acceleration
            
            m_Velocity = m_Velocity + half * dt * &
                m_Acceleration * damping_factor

        end do

        do atom_num = 1, m_num_atoms
            kinetic_energy = kinetic_energy + half * m_masses(atom_num) &
                * ( m_Velocity(1, atom_num) ** 2 + &
                    m_Velocity(2, atom_num) ** 2 + &
                    m_Velocity(3, atom_num) ** 2 )
        end do

        Coords(:,:) = m_Coords(:,:)

    end subroutine

    subroutine e_stretches(Forces, stretch_es_tot, abs_gradient, status)
        real(dp), intent(inout) :: Forces(3, m_num_atoms)
        real(dp), intent(out) :: stretch_es_tot(2)
        real(dp), intent(out) :: abs_gradient
        integer*4, intent(out) :: status

        integer*4 :: stretch_num, j, a0, a1
        real(dp) :: v10(3), force(3)
        real(dp) :: stretch_es(2)
        real(dp) :: r10, dr, v

        status=OK

        do stretch_num = 1, m_num_stretches
            
            a0 = m_stretch_atom_nums(stretch_num * int_two - int_one) + int_one
            a1 = m_stretch_atom_nums(stretch_num * int_two) + int_one

            v = m_stretch_keqs(stretch_num)
            if (v.eq.0) then
                cycle
            end if
            
            !Get positions, vector and distance
            call get_r(a1, a0, v10, r10)

            if (r10.eq.0) then
                bad_a0 = a0 - 1
                bad_a1 = a1 - 1
                status=STRETCH_EXC
            end if
            
            dr = r10 - m_stretch_reqs(stretch_num)
            
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

    subroutine e_bends(Forces, bend_es_tot, abs_gradient, status)
        real(dp), intent(inout) :: Forces(3, m_num_atoms)
        real(dp), intent(out) :: bend_es_tot(2)
        real(dp), intent(out) :: abs_gradient
        integer*4, intent(out) :: status

        integer*4 :: bend_num, j, a0, a1, a2, bend_index
        real(dp) :: v01n(3), v21n(3), force0(3), force2(3), w012(3)
        real(dp) :: bend_es(2)
        real(dp) :: r01, r21, v, a012, da

        status=OK

        do bend_num = 1, m_num_bends
            bend_index = bend_num * int_three

            v = m_bend_keqs(bend_num)
            if (v.eq.0) then
                cycle
            end if

            a0 = m_bend_atom_nums(bend_index - int_two) + int_one
            a1 = m_bend_atom_nums(bend_index - int_one) + int_one
            a2 = m_bend_atom_nums(bend_index) + int_one

            call get_a(a0, a1, a2, v01n, v21n, r01, r21, a012, w012, status)

            if (a012.ge.almost_pi) then
                bad_a0 = a0
                bad_a1 = a1
                bad_a2 = a2
                status = BEND_EXC
            end if

            if (status.ne.OK) then
                return
            end if
            
            da = a012 - m_bend_aeqs(bend_num)

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

    subroutine e_torsions(Forces, torsions_es_tot, abs_gradient, status)
        real(dp), intent(inout) :: Forces(3, m_num_atoms)
        real(dp), intent(out) :: torsions_es_tot(2)
        real(dp), intent(out) :: abs_gradient
        integer*4, intent(out) :: status

        integer*4 :: torsion_num, torsion_index, periodicity, j, &
        a0, a1, a2, a3
        real(dp) :: v01n(3), v21n(3), v23n(3), w012n(3), w123n(3), &
        force0(3), force1(3), force2(3), force3(3), torque3(3), vo3(3)
        real(dp) :: torsion_es(2)
        real(dp) :: r01, r21, r23, a012, a123, d0123, dd, v
        real(dp) :: e_da_dr10, e_da_dr23, ro3ro3

        !Reference:
		!http://xray.bmc.uu.se/~aqwww/q_legacy/documents/qman5.pdf

        do torsion_num = 1, m_num_torsions

            torsion_index = (torsion_num - int_one) * int_four
            
            a0 = m_torsion_atom_nums(torsion_index + int_one) + int_one
            a1 = m_torsion_atom_nums(torsion_index + int_two) + int_one
            a2 = m_torsion_atom_nums(torsion_index + int_three) + int_one
            a3 = m_torsion_atom_nums(torsion_index + int_four) + int_one

            !Get dihedral, angles, vectors, norms and distances
            call get_d(a0, a1, a2, a3, &
                v01n, v21n, v23n, r01, r21, r23, &
                a012, a123, w012n, w123n, d0123, status)


            if (a012.eq.zero.or.a123.eq.zero) then
                bad_a0 = a0
                bad_a1 = a1
                bad_a2 = a2
                bad_a3 = a3
                status = DIHEDRAL_EXC
            end if

            if (status.ne.OK) then
                return
            end if

            !Compute energy and gradient for dihedral
            torsion_es(1) = zero
            torsion_es(2) = zero
            do periodicity = 1, 4
                v = half * m_torsions_vs(torsion_index + periodicity) / &
                (m_torsion_paths(torsion_index + int_one))

                dd = periodicity * d0123 - &
                m_torsions_gammas(torsion_index + periodicity)

                torsion_es(1) = torsion_es(1) + v * (one + cos(dd))
                torsion_es(2) = torsion_es(2) - v * sin(dd)
            end do

            torsions_es_tot(1) = torsions_es_tot(1) + torsion_es(1)
            torsions_es_tot(2) = torsions_es_tot(2) + torsion_es(2)
            abs_gradient = abs_gradient + abs(torsion_es(2))

            !Energy times da/dr to convert to cartesian
            e_da_dr10 = torsion_es(1) / (r01 * sin(a012))
            e_da_dr23 = torsion_es(1) / (r23 * sin(a123))

            !Compute forces on outside atoms and 
            ! distance from centre of 1-2 to 3 squared
            ro3ro3 = 0
            do j = 1, 3
                force0(j) = - w012n(j) * e_da_dr10
                force3(j) = - w123n(j) * e_da_dr23

                vo3(j) = m_Coords(j, a3) - &
                    (m_Coords(j, a1) + m_Coords(j, a2)) * half

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

    subroutine e_nonbonding(Forces, vdw_es_tot, coul_es_tot, &
        abs_gradient, status)
        real(dp), intent(inout) :: Forces(3, m_num_atoms)
        real(dp), intent(out) :: vdw_es_tot(2), coul_es_tot(2)
        real(dp), intent(out) :: abs_gradient
        integer*4, intent(out) :: status

        integer*4 :: non_bonding_num, j, a0, a1
        real(dp) :: v10(3), force(3)
        real(dp) :: vdw_es(2), coul_es(2)
        real(dp) :: r10, r_1, r_2, r_6, r_12, v

        vdw_es(1) = 0
        vdw_es(2) = 0
        coul_es(1) = 0
        coul_es(2) = 0
        status = OK
        abs_gradient = 0

        if (m_diel.eq.0) then
            status = NONBON_EXC2
        end if

        do non_bonding_num = 1, m_num_non_bonding
            a0 = m_non_bonding_atom_nums(non_bonding_num * int_two - int_one) + int_one
            a1 = m_non_bonding_atom_nums(non_bonding_num * int_two) + int_one

            !Get positions, vector and distance
            call get_r(a1, a0, v10, r10)

            if (r10.eq.0) then
                bad_a0 = a0
                bad_a1 = a1
                status = NONBON_EXC1
            end if

            if (status.ne.OK) then
                return
            end if

            !Compute VdW
            v = m_vdw_vs(non_bonding_num)
            if (r10.le.m_vdw_cutoff.and.v.gt.0) then
                r_1 = m_vdw_reqs(non_bonding_num) / r10
                r_2 = r_1 * r_1
                r_6 = r_2 * r_2 * r_2
                r_12 = v * r_6 * r_6
                r_6 = v * r_6

                vdw_es(1) = r_12 - two * r_6
                vdw_es(2) = r_6 - r_12
                vdw_es(2) = vdw_es(2) * twelve * r_1 / r10

                vdw_es_tot(1) = vdw_es_tot(1) + vdw_es(1)
                vdw_es_tot(2) = vdw_es_tot(2) + vdw_es(2)
                abs_gradient = abs_gradient + abs(vdw_es(2))


            else if (r10.ge.m_coul_cutoff) then
                !skip if outside of both vdw_cutoff and coul_cutoff
                cycle
            end if

            !Compute electrostatic
            v = m_q0q1s(non_bonding_num)
            if (v.ne.0.and.m_diel.gt.0) then
                r_1 = 1 / r10
                coul_es(1) = v * r_1 / m_diel
                coul_es(2) = - coul_es(1) * r_1
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

    subroutine get_r(a0, a1, v01, r01)
        integer*4, intent(in) :: a0, a1
        real(dp), intent(out) :: v01(3)
        real(dp), intent(out) :: r01

        real(dp) :: r01s
        integer*4 :: j

        r01s = 0
        do j = 1, 3
            v01(j) = m_Coords(j, a1) - m_Coords(j, a0)
            r01s = r01s + v01(j) * v01(j)
        end do
        r01 = sqrt(r01s)
        
    end subroutine get_r

    subroutine get_a(a0, a1, a2, v01n, v21n, r01, r21, a012, w012n, status)
        integer*4, intent(in) :: a0, a1, a2
        real(dp), intent(out) :: v01n(3), v21n(3)
        real(dp), intent(out) :: r01, r21
        real(dp), intent(out) :: a012
        real(dp), intent(out) :: w012n(3)
        integer*4, intent(out) :: status

        call get_r(a0, a1, v01n, r01)
        call get_r(a2, a1, v21n, r21)

        if (r01.eq.0.or.r21.eq.0) then
            bad_a0 = a0
            bad_a1 = a1
            bad_a2 = a2
            status=ANGLE_EXC
        end if

        call divide_vec(v01n, r01)
        call divide_vec(v21n, r21)
        w012n = cross_product(v01n, v21n)
        a012 = unsigned_angle(v01n, v21n)
        
    end subroutine get_a

    subroutine get_d(a0, a1, a2, a3, &
        v01n, v21n, v23n, r01, r21, r23, &
        a012, a123, w012n, w123n, d0123, status)
        integer*4, intent(in) :: a0, a1, a2, a3
        real(dp), intent(out) :: v01n(3), v21n(3), v23n(3)
        real(dp), intent(out) :: r01, r21, r23
        real(dp), intent(out) :: a012, a123
        real(dp), intent(out) :: w012n(3), w123n(3)
        real(dp), intent(out) :: d0123
        integer*4, intent(out) :: status

        real(dp) :: c0123(3)

        call get_a(a0, a1, a2, v01n, v21n, r01, r21, a012, w012n, status)
        call get_a(a1, a2, a3, v21n, v23n, r21, r23, a123, w123n, status)

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

    function unsigned_angle(v0n, v1n)
        real(dp), intent(in) :: v0n(3), v1n(3)
        real(dp) :: unsigned_angle

        unsigned_angle = acos(dot_product(v0n, v1n))
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

end module AmberFortMod
