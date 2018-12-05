module AmberFortMod
    implicit none
    integer, parameter:: dp=kind(0.d0)
    integer*4 :: m_num_atoms=0, m_num_stretches=0, m_num_bends=0, &
    m_num_torsions=0, m_num_non_bonding=0
    real(dp) :: m_energy, m_prev_energy
    real(dp) :: m_step_size = 0.001d0
    real(dp) :: m_lambda = 1d0
    real(dp) :: m_mag0
    real(dp) :: m_step_reduce = 0.5d0
    real(dp) :: m_step_increase = 1.25d0
    integer*4 :: m_step_num = 1
    integer*4 :: m_num_steps

    integer*4, allocatable :: m_stretch_atom_nums(:), &
    m_bend_atom_nums(:), m_torsion_atom_nums(:), &
    m_non_bonding_atom_nums(:)

    real(dp), allocatable :: m_stretch_reqs(:), m_stretch_keqs(:), &
    m_bend_aeqs(:), m_bend_keqs(:), &
    m_torsions_vs(:), m_torsions_gammas(:), m_torsion_paths(:), &
    m_vdw_vs(:), m_vdw_reqs(:), m_q0q1s(:), m_masses(:)

    real(dp), allocatable :: m_Forces(:), & 
    m_Acceleration(:), m_Velocity(:)

    real(dp) :: m_vdw_cutoff, m_diel, m_coul_cutoff

    integer*4 :: bad_a0, bad_a1, bad_a2, bad_a3

    enum, bind(c)
        enumerator :: OK=0, ANGLE_EXC=1, STRETCH_EXC=2, BEND_EXC=3, &
        DIHEDRAL_EXC=4, NONBON_EXC1=5, NONBON_EXC2=6, NOATOMS_EXC=7, &
        NOMASS_EXC=8, NPATHS_EXC=9, OTHER_EXC = 10
    endenum

    !Bools to check if things are allocated or not.
    !Supposedly safer than allocated()
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
    real(dp), parameter :: almost_zero=0.00000001d0, zero=0.d0, half=0.5d0, &
    one=1.d0, two=2.d0, almost_pi=(4-almost_zero)*atan(1.d0), twelve=12.d0

    contains

    subroutine safe_deallocate_i1(A, is_allocated)
        integer*4, allocatable, intent(inout) :: A(:)
        logical, intent(inout) :: is_allocated
        if (is_allocated) then
            deallocate(A)
        end if
        is_allocated = .false.
    end subroutine

    subroutine reallocate_i1(A, is_allocated, dim1)
        integer*4, intent(in) :: dim1
        integer*4, allocatable, intent(inout) :: A(:)
        logical, intent(inout) :: is_allocated

        call safe_deallocate_i1(A, is_allocated)
        allocate(A(dim1))
        is_allocated = .true.

    end subroutine reallocate_i1

    subroutine safe_deallocate_r1(A, is_allocated)
        real(dp), allocatable, intent(inout) :: A(:)
        logical, intent(inout) :: is_allocated
        if (is_allocated) then
            deallocate(A)
        end if
        is_allocated = .false.
    end subroutine

    subroutine reallocate_r1(A, is_allocated, dim1)
        integer*4, intent(in) :: dim1
        real(dp), allocatable, intent(inout) :: A(:)
        logical, intent(inout) :: is_allocated

        call safe_deallocate_r1(A, is_allocated)
        allocate(A(dim1))
        is_allocated = .true.

    end subroutine reallocate_r1

    subroutine deallocate_all()
        call safe_deallocate_i1( &
            m_stretch_atom_nums, m_stretch_atom_nums_allocated)
        call safe_deallocate_i1( &
            m_bend_atom_nums, m_bend_atom_nums_allocated)
        call safe_deallocate_i1( &
            m_torsion_atom_nums, m_torsion_atom_nums_allocated)
        call safe_deallocate_i1( &
            m_non_bonding_atom_nums, m_non_bonding_atom_nums_allocated)
    
        call safe_deallocate_r1( &
            m_stretch_reqs, m_stretch_reqs_allocated)
        call safe_deallocate_r1( &
            m_stretch_keqs, m_stretch_keqs_allocated)
        call safe_deallocate_r1( &
            m_bend_aeqs, m_bend_aeqs_allocated)
        call safe_deallocate_r1( &
            m_bend_keqs, m_bend_keqs_allocated)
        call safe_deallocate_r1( &
            m_torsions_vs, m_torsions_vs_allocated)
        call safe_deallocate_r1( &
            m_torsions_gammas, m_torsions_gammas_allocated)
        call safe_deallocate_r1( &
            m_torsion_paths, m_torsion_paths_allocated)
        call safe_deallocate_r1( &
            m_vdw_vs, m_vdw_vs_allocated)
        call safe_deallocate_r1( &
            m_vdw_reqs, m_vdw_reqs_allocated)
        call safe_deallocate_r1( &
            m_q0q1s, m_q0q1s_allocated)
        call safe_deallocate_r1( &
            m_masses, m_masses_allocated)
    
        call safe_deallocate_r1( &
            m_Forces, m_Forces_allocated)
        call safe_deallocate_r1( &
            m_Acceleration, m_Acceleration_allocated)
        call safe_deallocate_r1( &
            m_Velocity, m_Velocity_allocated)

        m_num_atoms=0
        m_num_stretches=0
        m_num_bends=0
        m_num_torsions=0
        m_num_non_bonding=0
    end subroutine

    subroutine get_bad_atoms(a0, a1, a2, a3)
        integer*4, intent(out) :: a0, a1, a2, a3
        a0 = bad_a0
        a1 = bad_a1
        a2 = bad_a2
        a3 = bad_a3

    end subroutine get_bad_atoms

    subroutine allocate_atom_arrays(num_atoms)
        integer*4, intent(in) :: num_atoms

        m_step_num = 1
        m_num_atoms = num_atoms
        bad_a0 = -1
        bad_a1 = -1
        bad_a2 = -1
        bad_a3 = -1

        call reallocate_r1(m_Forces, m_Forces_allocated, 3 * m_num_atoms)
        m_Forces = 0.d0

        call reallocate_r1(m_Acceleration, m_Acceleration_allocated, &
            3 * m_num_atoms)
        m_Acceleration = 0.d0

        call reallocate_r1(m_Velocity, m_Velocity_allocated, 3 * m_num_atoms)
        m_Velocity = 0.d0

    end subroutine allocate_atom_arrays
        
    subroutine set_masses(masses, status)
        real(dp), intent(in) :: masses(m_num_atoms)
        integer*4, intent(inout) :: status

        integer*4 :: atom_num
        real(dp) :: mass

        if (m_num_atoms.eq.0) then
            status = NOATOMS_EXC
            return
        end if

        call reallocate_r1(m_masses, &
            m_masses_allocated, 3 * m_num_atoms)

        do atom_num = 1, m_num_atoms
            mass = masses(atom_num)
            if(mass.le.almost_zero) then
                status = NOMASS_EXC
                bad_a0 = atom_num - 1
                return
            end if
            m_masses(3 * atom_num - 2) = mass
            m_masses(3 * atom_num - 1) = mass
            m_masses(3 * atom_num) = mass
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

    subroutine e_amber(Coords, Forces, &
        stretch_es_tot, bend_es_tot, &
        torsions_es_tot, vdw_es_tot, coul_es_tot, &
        norm_forces, status)
        real(dp), intent(inout) :: Coords(3 * m_num_atoms)
        real(dp), intent(inout) :: Forces(3 * m_num_atoms)
        real(dp), intent(out) :: stretch_es_tot(2), bend_es_tot(2), &
        torsions_es_tot(2), vdw_es_tot(2), coul_es_tot(2), norm_forces
        integer*4, intent(inout) :: status

        Forces = 0.0d0

        call e_stretches(Coords, Forces, stretch_es_tot, status)
        if (status.ne.OK) then
            return
        end if

        call e_bends(Coords, Forces, bend_es_tot, status)
        if (status.ne.OK) then
            return
        end if

        call e_torsions(Coords, Forces, torsions_es_tot, status)
        if (status.ne.OK) then
            return
        end if

        call e_nonbonding(Coords, Forces, vdw_es_tot, coul_es_tot, status)
        if (status.ne.OK) then
            return
        end if
        
        norm_forces = sqrt(scalar_product1(Forces, 3 * m_num_atoms))

    end subroutine

    subroutine md_verlet(Coords, steps, dt, damping_factor, &
        stretch_es_tot, bend_es_tot, &
        torsions_es_tot, vdw_es_tot, coul_es_tot, &
        kinetic_energy, norm_forces, status)
        real(dp), intent(inout) :: Coords(3 * m_num_atoms)
        real(dp), intent(out) :: stretch_es_tot(2), bend_es_tot(2), &
        torsions_es_tot(2), vdw_es_tot(2), coul_es_tot(2)
        real(dp), intent(out) :: kinetic_energy, dt, damping_factor, norm_forces
        integer*4, intent(in) :: steps
        integer*4, intent(inout) :: status

        integer*4 :: step_num, atom_num, j

        kinetic_energy = 0

        do step_num = 1, steps

            !Velocity Verlet update
            Coords = Coords + m_Velocity * dt + &
                half * (dt ** 2) * m_Acceleration
            
            m_Velocity = m_Velocity * damping_factor + half * dt * &
                m_Acceleration

            stretch_es_tot = 0
            bend_es_tot = 0
            torsions_es_tot = 0
            vdw_es_tot = 0
            coul_es_tot = 0
            call e_amber(Coords, m_Forces, stretch_es_tot, bend_es_tot, &
                torsions_es_tot, vdw_es_tot, coul_es_tot, norm_forces, status)

            if (status.ne.OK) then
                return
            end if

            m_Acceleration = m_Forces / m_masses
            
            m_Velocity = m_Velocity * damping_factor + half * dt * &
                m_Acceleration

        end do

        kinetic_energy = sum(half * m_masses * m_Velocity ** 2)

    end subroutine

    subroutine optimise(Coords, steps, &
        stretch_es_tot, bend_es_tot, &
        torsions_es_tot, vdw_es_tot, coul_es_tot, norm_forces, status)

        ! 
        ! OPTIMISE ROUTINE
        ! 
        ! Minimise energy wrt coordinates using AMBER forcefield
        !
        ! step_size (t): Maximum distance to step
        ! lambda (l): step_size scaler 
            
        ! 1:  Compute energy (e(n)) + forces (F(n))
        ! 2:  Set e(n-1)=e(n), F(n-1)=F(n), and g(n-1)=g(n)
        ! 3:  Set Prev_Coords = Coords
        ! 4:  Take step: Coords = Coords + F(n) * l * t
        ! 5:  Compute (e(n)), (F(n)), (g(n))
        ! 6:  If e(n) < e(n-1) goto 7, else goto 9
        ! 7:  Set l=1
        ! 8:  Goto 2
        ! 9:  Compute l as the minimum of a cubic function (see below)
        ! 10: Revert e(n)=e(n-1), F(n)=F(n-1), and g(n)=g(n-1)
        ! 11: Revert Coords = Prev_Coords
        ! 12: Goto 4
        ! 
        ! Cubic function minimisation justification:
        ! 
        ! Energy must go down in the forward direction (F(n))
        ! If the energy after a step has increased, there must be some lambda
        !   that yields the lowest energy in the forward direction.
        ! Forces are calculated "for free" in each step, so we can set
        !   a function of energy as a cubic function using the gradient
        !   information at the starting and finishing point. The minimum
        !   of this is the solution of a quadratic equation.
        !
        ! Let x be a value between 0 and 1, corresponding to the points
        !   between the current position and the stepped position where
        !   the energy went up. 
        ! 
        ! e = a*x**3 + b*x**2 + c*x + d
        ! g = 3*a*x**2 + 2*b*x + c
        !
        ! -g(0) = c
        !  g(1) = 3*a + 2*b + c
        !  de = e(1) - e(0) = a + b + c
        ! 
        ! This yields:
        !
        ! a = g(1) - g(0) - 2de
        ! b = 3de + 2g(0) - g(1)
        ! c = g(0)
        ! d = 0
        !
        ! We want to solve g(l) = 0 for l
        !
        ! 3 * (g(1) - g(0) - 2de) * l**2 + 
        !   2 * (3de + 2g(0) - g(1)) * l - g(0) = 0
        !
        ! or simplifying:
        !
        ! v0 = 3 * (g(1) - g(0) - 2de) = 3a
        ! v1 = 2 * (3de + 2g(0) - g(1)) = 2b
        ! v0 * l**2 + v1 * l + g(0) = 0
        !
        ! which can be solved with the quadratic equation:
        ! 
        ! l = (- v1 - sqrt(v1**2 - 4 * v0 * g(0))) / (2 * v0)
        ! l = (- v1 + sqrt(v1**2 - 4 * v0 * g(0))) / (2 * v0)
        ! 
        ! The discriminant can be computed as:
        !
        !
        !  36*de**2   - 24*de*g(0) - 24*de*g(1)
        !   4g(0)**2 +  4g(0)g(1) +  4g(1)**2
        ! 
        ! This has a minimum at:
        ! 
        ! g(1) = - (4g(0)-24d(0)) / 8
        !
        ! which in turn has a value of:
        !
        ! min(discriminant) = 3g(0)**2 - 12*de*g(0)
        !
        ! As long as de is positive (energy increased, which is a criterion)
        !   and g(0) is negative, which it is always set as, this quadratic
        !   will always have a solution

        real(dp), intent(inout) :: Coords(3 * m_num_atoms)
        real(dp), intent(out) :: stretch_es_tot(2), bend_es_tot(2), &
        torsions_es_tot(2), vdw_es_tot(2), coul_es_tot(2), norm_forces
        integer*4, intent(in) :: steps
        integer*4, intent(inout) :: status

        real(dp) :: Test_Coords(3 * m_num_atoms), Test_Forces(3 * m_num_atoms)
        real(dp) :: test_energy, test_norm_forces, &
        mag1, de, a, b, c, s0, s1, ep0, ep1, discriminant
        
        m_num_steps = m_step_num + steps - 1

        !TEMP
        !m_step_size = 0.0001

        if (m_step_num.eq.1) then

            m_lambda = 1

            stretch_es_tot = 0
            bend_es_tot = 0
            torsions_es_tot = 0
            vdw_es_tot = 0
            coul_es_tot = 0
            call e_amber(Coords, m_Forces, stretch_es_tot, bend_es_tot, & 
                torsions_es_tot, vdw_es_tot, coul_es_tot, norm_forces, status)
            if (status.ne.OK) then
                return
            end if
            m_energy = stretch_es_tot(1) + bend_es_tot(1) + &
            torsions_es_tot(1) + vdw_es_tot(1) + coul_es_tot(1)
        end if

        do m_step_num = m_step_num, m_num_steps


            write(*,*) "Prev Coords"
            write(*,*) Coords
            write(*,*) "Prev Energy"
            write(*,*) m_energy
            write(*,*) "Current Forces"
            write(*,*) m_Forces
            write(*,*) "Current Norm Forces"
            write(*,*) norm_forces
            write(*,*) "Step Vector"
            write(*,*) m_Forces * m_step_size * m_lambda
            write(*,*) "Step Size", m_step_size, "lambda", m_lambda

            Test_Coords = Coords + m_Forces * m_step_size * m_lambda

            stretch_es_tot = 0
            bend_es_tot = 0
            torsions_es_tot = 0
            vdw_es_tot = 0
            coul_es_tot = 0
            call e_amber(Test_Coords, Test_Forces, stretch_es_tot, &
                bend_es_tot, torsions_es_tot, vdw_es_tot, coul_es_tot, &
                test_norm_forces, status)
            if (status.ne.OK) then
                return
            end if

            test_energy = stretch_es_tot(1) + bend_es_tot(1) + &
            torsions_es_tot(1) + vdw_es_tot(1) + coul_es_tot(1)

            de = test_energy - m_energy

            if (de.le.0) then
                write(*,*) "Energy dropped. de:", de
                call flush()
                m_lambda = 1

                m_prev_energy = m_energy
                m_energy = test_energy
                norm_forces = test_norm_forces
                Coords = Test_Coords
                m_Forces = Test_Forces

                m_mag0 = - norm_forces

                write(*,*) "a"
                write(*,*) 0
                write(*,*) "b"
                write(*,*) 0
                write(*,*) "c"
                write(*,*) 0
                
            else
                write(*,*) "Energy rose. de:", de

                !Energy went up. Use a cubic interpolation to choose
                !a ratio to multiply the step by
!
                m_mag0 = - norm_forces
                mag1 = scalar_product2(Test_Forces, m_Forces, 3 * m_num_atoms) &
                    / m_mag0
!
                a = mag1 + m_mag0 - 2 * de
                b = 3 * de - 2 * m_mag0 - mag1
                c = m_mag0

                write(*,*) "a"
                write(*,*) a
                write(*,*) "b"
                write(*,*) b
                write(*,*) "c"
                write(*,*) c

                if (abs(a).lt.almost_zero) then
                !Rare case when interpolated function is quadratic.
                !Use the solution of a quadratic instead
                    m_lambda = m_lambda * ( - c / 2 * b)
                    write(*,*) "a approached 0. Used quadratic. lambda:", m_lambda
                else
                    discriminant = 4 * b ** 2 - 12 * a * c
                    if (discriminant.lt.almost_zero) then
                        write(*,*) "discriminant approached 0. Used quadratic. lambda:", m_lambda
                        m_lambda = m_lambda * ( - c / 2 * b)
                    else
                        !Two solutions corresponding to the predicted stationary
                        !points of the energy along the vector
                        s0 = (- 2 * b + sqrt(discriminant)) / (6 * a)
                        s1 = (- 2 * b - sqrt(discriminant)) / (6 * a)

                        !Predicted energies of the geometries at these points
                        ep0 = a * s0 ** 3 + b * s0 ** 2 + c * s0
                        ep1 = a * s1 ** 3 + b * s1 ** 2 + c * s1
                        if (ep0.lt.ep1) then
                            m_lambda = m_lambda * s0
                        else 
                            m_lambda = m_lambda * s1
                        end if
                        write(*,*) "g0", m_mag0, "g1", mag1, "lambda", m_lambda
                    end if 
                end if
            end if

            write(*,*) "New Coords"
            write(*,*) Test_Coords
            write(*,*) "New Energy"
            write(*,*) test_energy
            write(*,*) "Test Forces"
            write(*,*) Test_Forces
            write(*,*) "New Norm Forces"
            write(*,*) sqrt(scalar_product1(Test_Forces, 3 * m_num_atoms))
            write(*,*) ""

            call flush()
        end do
    end subroutine

    subroutine e_stretches(Coords, Forces, stretch_es_tot, status)
        real(dp), intent(inout) :: Coords(3 * m_num_atoms)
        real(dp), intent(inout) :: Forces(3 * m_num_atoms)
        real(dp), intent(out) :: stretch_es_tot(2)
        integer*4, intent(inout) :: status

        integer*4 :: stretch_num, j, a0, a1
        real(dp) :: v10(3), force(3)
        real(dp) :: stretch_es(2)
        real(dp) :: r10, dr, v

        status=OK

        if (m_num_stretches.eq.0) then
            return
        end if

        do stretch_num = 1, m_num_stretches
            
            a0 = m_stretch_atom_nums(stretch_num * int_two - int_one) + int_one
            a1 = m_stretch_atom_nums(stretch_num * int_two) + int_one

            v = m_stretch_keqs(stretch_num)
            if (v.eq.0) then
                cycle
            end if
            
            !Get positions, vector and distance
            call get_r(Coords, a1, a0, v10, r10)

            if (r10.le.almost_zero) then
                bad_a0 = a0 - 1
                bad_a1 = a1 - 1
                status=STRETCH_EXC
                return
            end if
            
            dr = r10 - m_stretch_reqs(stretch_num)
            
            stretch_es(1) =  v * dr ** 2
            stretch_es(2) = two * v * dr
        
            stretch_es_tot(1) = stretch_es_tot(1) + stretch_es(1)
            stretch_es_tot(2) = stretch_es_tot(2) + stretch_es(2)
            
            do j = 1, 3
                !Use normalised vector
                force(j) = v10(j) * stretch_es(2) / r10
                Forces(a0*3 - 3 + j) = Forces(a0*3 - 3 + j) - force(j)
                Forces(a1*3 - 3 + j) = Forces(a1*3 - 3 + j) + force(j)
            end do
            
        end do

    end subroutine e_stretches

    subroutine e_bends(Coords, Forces, bend_es_tot, status)
        real(dp), intent(inout) :: Coords(3 * m_num_atoms)
        real(dp), intent(inout) :: Forces(3 * m_num_atoms)
        real(dp), intent(out) :: bend_es_tot(2)
        integer*4, intent(inout) :: status

        integer*4 :: bend_num, j, a0, a1, a2, bend_index
        real(dp) :: v01n(3), v21n(3), force0(3), force2(3), w012(3)
        real(dp) :: bend_es(2)
        real(dp) :: r01, r21, v, a012, da

        status=OK

        if (m_num_bends.eq.0) then
            return
        end if

        do bend_num = 1, m_num_bends
            bend_index = bend_num * int_three

            v = m_bend_keqs(bend_num)
            if (v.eq.0) then
                cycle
            end if

            a0 = m_bend_atom_nums(bend_index - int_two) + int_one
            a1 = m_bend_atom_nums(bend_index - int_one) + int_one
            a2 = m_bend_atom_nums(bend_index) + int_one

            call get_a(Coords, a0, a1, a2, v01n, v21n, &
            r01, r21, a012, w012, status)

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

            bend_es(1) = v * da ** 2
            bend_es(2) = two * v * da

            bend_es_tot(1) = bend_es_tot(1) + bend_es(1)
            bend_es_tot(2) = bend_es_tot(2) + bend_es(2)

            force0 = cross_product(v01n, w012)
            force2 = cross_product(v21n, w012)
            
            force0 = force0 * bend_es(2)
            force2 = - force2 * bend_es(2)

            do j = 1, 3
                Forces(a0*3 - 3 + j) = Forces(a0*3 - 3 + j) + force0(j)
                Forces(a1*3 - 3 + j) = Forces(a1*3 - 3 + j) - force0(j) - &
                    force2(j)
                Forces(a2*3 - 3 + j) = Forces(a2*3 - 3 + j) + force2(j)
            end do

        end do

    end subroutine e_bends

    subroutine e_torsions(Coords, Forces, torsions_es_tot, status)
        real(dp), intent(inout) :: Coords(3 * m_num_atoms)
        real(dp), intent(inout) :: Forces(3 * m_num_atoms)
        real(dp), intent(out) :: torsions_es_tot(2)
        integer*4, intent(inout) :: status

        integer*4 :: torsion_num, torsion_index, periodicity, j, &
        a0, a1, a2, a3
        real(dp) :: v01n(3), v21n(3), v23n(3), w012n(3), w123n(3), &
        force0(3), force1(3), force2(3), force3(3), tc(3)
        real(dp) :: torsion_es(2)
        real(dp) :: r01, r21, r23, a012, a123, d0123, dd, v, npaths
        real(dp) :: e_da_dr10, e_da_dr23

        !References:
        ! (1) http://xray.bmc.uu.se/~aqwww/q_legacy/documents/qman5.pdf
        ! (2) https://arxiv.org/pdf/1401.1181.pdf

        ! Dihedral forces on abcd:
        ! Force on a = fa etc
        ! o is the centre of bond bc
        ! 
        ! From now on, everything is a vector (eg ab is the vector from a to b)
        !
        ! CONSTRAINTS:
        ! (1) Forces cancel:
        ! fa + fb + fc + fd = 0
        ! 
        ! (2) Torques cancel:
        ! oa x fa + ob x fb + oc x fc + od x fd = 0
        !
        ! Determine torque on c (tc):
        ! from 2:
        ! (ob - ab) x fa + ob x fb + oc x fc + (oc + cd) x fd = 0
        ! 
        ! ob = -oc
        ! (- oc - ab) x fa - oc x fb + oc x fc + (oc + cd) x fd = 0
        ! oc x (- fa - fb + fc + fd) - ab x fa + cd x fd = 0
        ! 
        ! From 1:
        ! oc x (2*fc + 2*fd) - ab x fa + cd x fd = 0
        ! 2 * oc x fc + 2 * oc x fd - ab x fa + cd x fd = 0
        ! oc x fc = - oc x fd + 0.5 * ab x fa - 0.5 * cd x fd = tc
        !
        ! Determine fc and thus fb using tc:
        ! Using the perpendicular solution from reference 2:
        ! fc = tc x oc / (|oc||oc|)
        !
        ! fc = (- oc x fd + 0.5 * ab x fa - 0.5 * cd x fd) x oc / (|oc||oc|)
        !
        ! Substitutions:
        ! cb is already computed as v21n * r21 = - 2 * oc:
        ! oc = - v21n * r21 / 2
        ! ab = v01n * r01
        ! cd = v23n * r23
        ! |oc| = r21 / 2
        !
        ! fc = - (2 / r21) *
        !        ( (r21 / 2) * v21n x fd + 
        !          (r01 / 2) * v01n x fa -
        !          (r23 / 2) * v23n x fd   ) x 
        !        v21n
        !    = ( - r21 * v21n x fd
        !        - r01 * v01n x fa
        !        + r23 * v23n x fd ) x 
        !        v21n
        ! 
        ! fb, from 1 is simply:
        ! fb = - fa - fc - fd


        status=OK

        if (m_num_torsions.eq.0) then
            return
        end if

        do torsion_num = 1, m_num_torsions

            torsion_index = (torsion_num - int_one) * int_four
            
            a0 = m_torsion_atom_nums(torsion_index + int_one) + int_one
            a1 = m_torsion_atom_nums(torsion_index + int_two) + int_one
            a2 = m_torsion_atom_nums(torsion_index + int_three) + int_one
            a3 = m_torsion_atom_nums(torsion_index + int_four) + int_one

            npaths = m_torsion_paths(torsion_num)
            if (npaths.le.0) then
                bad_a0 = a0
                bad_a1 = a1
                bad_a2 = a2
                bad_a3 = a3
                status = NPATHS_EXC
            end if

            !Get dihedral, angles, vectors, norms and distances
            call get_d(Coords, a0, a1, a2, a3, &
                v01n, v21n, v23n, r01, r21, r23, &
                a012, a123, w012n, w123n, d0123, status)

            write(*,*) "dihedral", d0123
            write(*,*) "angles", a012, a123
            write(*,*) "lengths", r01, r21, r23

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
                npaths

                dd = periodicity * d0123 - &
                m_torsions_gammas(torsion_index + periodicity)

                torsion_es(1) = torsion_es(1) + v * (one + cos(dd))
                torsion_es(2) = torsion_es(2) - v * sin(dd)
            end do

            torsions_es_tot(1) = torsions_es_tot(1) + torsion_es(1)
            torsions_es_tot(2) = torsions_es_tot(2) + torsion_es(2)

            !Energy times da/dr to convert to cartesian
            e_da_dr10 = torsion_es(1) / (r01 * sin(a012))
            e_da_dr23 = torsion_es(1) / (r23 * sin(a123))

            !Compute forces on outside atoms and 
            ! distance from centre of 1-2 to 3 squared
            force0 = - w012n * e_da_dr10
            force3 = + w123n * e_da_dr23

            tc = (- cross_product(v21n, force3) - &
                    cross_product(v01n, force0) * (r01 / r21) + &
                    cross_product(v23n, force3) * (r23 / r21))

            force2 = cross_product(tc, v21n)

            do j = 1, 3
                force1(j) = - force0(j) - force2(j) - force3(j)

                Forces(a0*3 - 3 + j) = Forces(a0*3 - 3 + j) - force0(j)
                Forces(a1*3 - 3 + j) = Forces(a1*3 - 3 + j) - force1(j)
                Forces(a2*3 - 3 + j) = Forces(a2*3 - 3 + j) - force2(j)
                Forces(a3*3 - 3 + j) = Forces(a3*3 - 3 + j) - force3(j)

            end do 

        end do

    end subroutine e_torsions

    subroutine e_nonbonding(Coords, Forces, vdw_es_tot, coul_es_tot, status)
        real(dp), intent(inout) :: Coords(3 * m_num_atoms)
        real(dp), intent(inout) :: Forces(3 * m_num_atoms)
        real(dp), intent(out) :: vdw_es_tot(2), coul_es_tot(2)
        integer*4, intent(inout) :: status

        integer*4 :: non_bonding_num, j, a0, a1
        real(dp) :: v10(3), force(3)
        real(dp) :: vdw_es(2), coul_es(2)
        real(dp) :: r10, r_1, r_2, r_6, r_12, v

        vdw_es(1) = 0
        vdw_es(2) = 0
        coul_es(1) = 0
        coul_es(2) = 0

        status = OK

        if (m_num_non_bonding.eq.0) then
            return
        end if

        if (m_diel.le.almost_zero) then
            status = NONBON_EXC2
        end if

        do non_bonding_num = 1, m_num_non_bonding
            a0 = m_non_bonding_atom_nums(non_bonding_num * int_two - int_one) &
                + int_one
            a1 = m_non_bonding_atom_nums(non_bonding_num * int_two) + int_one

            !Get positions, vector and distance
            call get_r(Coords, a1, a0, v10, r10)

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

            else if (r10.ge.m_coul_cutoff) then
                !skip if outside of both vdw_cutoff and coul_cutoff
                cycle
            end if

            !Compute electrostatic
            v = m_q0q1s(non_bonding_num)
            if (v.ne.0) then
                r_1 = 1 / r10
                coul_es(1) = v * r_1 / m_diel
                coul_es(2) = - coul_es(1) * r_1

                coul_es_tot(1) = coul_es_tot(1) + coul_es(1)
                coul_es_tot(2) = coul_es_tot(2) + coul_es(2)
            end if

            do j = 1, 3
                force(j) = v10(j) * (vdw_es(2) + coul_es(2)) / r10
                Forces(a0*3 - 3 + j) = Forces(a0*3 - 3 + j) - force(j)
                Forces(a1*3 - 3 + j) = Forces(a1*3 - 3 + j) + force(j)
            end do

            call flush()

        end do

    end subroutine

    subroutine get_r(Coords, a1, a2, v12, r12)
        real(dp), intent(in) :: Coords(3 * m_num_atoms)
        integer*4, intent(in) :: a1, a2
        real(dp), intent(out) :: v12(3)
        real(dp), intent(out) :: r12

        real(dp) :: r12s
        integer*4 :: j, c_index1, c_index2

        r12s = zero
        c_index1 = a1 * 3 - 3
        c_index2 = a2 * 3 - 3
        do j = 1, 3
            v12(j) = Coords(c_index1 + j) - Coords(c_index2 + j)
            r12s = r12s + v12(j) * v12(j)
        end do
        r12 = sqrt(r12s)
        
    end subroutine get_r

    subroutine get_a(Coords, a1, a2, a3, v12n, v32n, &
        r12, r32, a123, w123n, status)
        real(dp), intent(in) :: Coords(3 * m_num_atoms)
        integer*4, intent(in) :: a1, a2, a3
        real(dp), intent(out) :: v12n(3), v32n(3)
        real(dp), intent(out) :: r12, r32
        real(dp), intent(out) :: a123
        real(dp), intent(out) :: w123n(3)
        integer*4, intent(inout) :: status

        call get_r(Coords, a1, a2, v12n, r12)
        call get_r(Coords, a3, a2, v32n, r32)

        if (r12.le.almost_zero.or.r32.le.almost_zero) then
            bad_a1 = a1
            bad_a2 = a2
            bad_a3 = a3
            status=ANGLE_EXC
        end if

        call divide_vec(v12n, r12)
        call divide_vec(v32n, r32)
        w123n = cross_product(v12n, v32n)
        a123 = unsigned_angle(v12n, v32n)
        
    end subroutine get_a

    subroutine get_d(Coords, a1, a2, a3, a4, &
        v12n, v32n, v34n, r12, r32, r34, &
        a123, a234, w123n, w234n, d1234, status)
        real(dp), intent(in) :: Coords(3 * m_num_atoms)
        integer*4, intent(in) :: a1, a2, a3, a4
        real(dp), intent(out) :: v12n(3), v32n(3), v34n(3)
        real(dp), intent(out) :: r12, r32, r34
        real(dp), intent(out) :: a123, a234
        real(dp), intent(out) :: w123n(3), w234n(3)
        real(dp), intent(out) :: d1234
        integer*4, intent(inout) :: status

        real(dp) :: c1234(3)

        call get_a(Coords, a1, a2, a3, v12n, v32n, &
        r12, r32, a123, w123n, status)
        call get_a(Coords, a2, a3, a4, v32n, v34n, &
        r32, r34, a234, w234n, status)

        d1234 = unsigned_angle(w123n, w234n)
        c1234 = cross_product(w123n, w234n)

        if (dot_product(v32n, c1234) < 0) then
            d1234 = - d1234
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

    function cross_product(v1, v2)
        real(dp), intent(in) :: v1(3), v2(3)
        real(dp) :: cross_product(3)
        
        cross_product(1) = v1(2) * v2(3) - v1(3) * v2(2)
        cross_product(2) = v1(3) * v2(1) - v1(1) * v2(3)
        cross_product(3) = v1(1) * v2(2) - v1(2) * v2(1)
    end function cross_product

    function unsigned_angle(v1n, v2n)
        real(dp), intent(in) :: v1n(3), v2n(3)
        real(dp) :: unsigned_angle

        unsigned_angle = acos(dot_product(v1n, v2n))
    end function unsigned_angle

    function scalar_product1(A1, dim1)
        real(dp), intent(in) :: A1(dim1)
        integer*4, intent(in) :: dim1
        real(dp) :: scalar_product1

        integer*4 :: i, j

        scalar_product1 = zero
        do j = 1, dim1
            scalar_product1 = scalar_product1 + A1(j) ** 2
        end do
        
    end function scalar_product1

    function scalar_product2(A1, A2, dim1)
        real(dp), intent(in) :: A1(dim1), A2(dim1)
        integer*4, intent(in) :: dim1
        real(dp) :: scalar_product2

        integer*4 :: i, j

        scalar_product2 = zero
        do i = 1, dim1
            scalar_product2 = scalar_product2 + A1(j) * A2(j)
        end do

    end function scalar_product2


    subroutine get_best_camera_view(Coords, num_atoms, &
        camera_position, field_of_view_rad, fill_amount, &
        min_distance)
        real(dp), intent(in) :: Coords(3 * num_atoms)
        real(dp), intent(out) :: camera_position(3)
        integer*4, intent(in) :: num_atoms
        real(dp), intent(in) :: field_of_view_rad, fill_amount, min_distance

        real(dp) :: centre(3)
        real(dp) :: squared_radius, largest_squared_radius, radius
        integer*4 :: c_index
        real(dp) :: ratio
        
        call get_centre(Coords, num_atoms, centre)
        !Align along x and y axes
        camera_position(1) = centre(1)
        camera_position(2) = centre(2)

        !Get largest distance from aligned axis
        largest_squared_radius = 0
        do c_index = 1, num_atoms * 3, 3
            squared_radius = &
            (Coords(c_index) - centre(1)) ** two + &
            (Coords(c_index + int_one) - centre(2)) ** two

            if (squared_radius.gt.largest_squared_radius) then
                largest_squared_radius = squared_radius
            end if
        end do

        radius = sqrt(largest_squared_radius) / fill_amount
        ratio = tan(field_of_view_rad / 2)
        camera_position(3) = min(-radius / ratio, min_distance)

    end subroutine get_best_camera_view

    subroutine get_centre(Coords, num_atoms, centre)
        real(dp), intent(in) :: Coords(3 * num_atoms)
        real(dp), intent(out) :: centre(3)
        integer*4, intent(in) :: num_atoms

        integer*4 :: c_index, j

        centre = zero
        do c_index = 1, num_atoms * 3, 3
            centre(j) = centre(j) + Coords(c_index)
            centre(j) = centre(j) + Coords(c_index + int_one)
            centre(j) = centre(j) + Coords(c_index + int_two)
        end do

        centre = centre / num_atoms

    end subroutine get_centre

    subroutine get_sizes(num_atoms, num_stretches, num_bends, &
        num_torsions, num_non_bonding)
        integer*4, intent(out) :: num_atoms, num_stretches, num_bends, &
            num_torsions, num_non_bonding

        num_atoms = m_num_atoms
        num_stretches = m_num_stretches
        num_bends = m_num_bends
        num_torsions = m_num_torsions
        num_non_bonding = m_num_non_bonding

    end subroutine get_sizes

    subroutine get_atoms_allocated(forces_allocated, acceleration_allocated, &
        velocity_allocated, masses_allocated)
        logical, intent(out) :: forces_allocated, acceleration_allocated, &
            velocity_allocated, masses_allocated

        forces_allocated = m_forces_allocated
        acceleration_allocated = m_acceleration_allocated
        velocity_allocated = m_velocity_allocated
        masses_allocated = m_masses_allocated

    end subroutine get_atoms_allocated

    subroutine get_calc_state(energy, prev_energy, step_size, step_reduce, &
        step_increase, step_num, num_steps)
        real(dp), intent(out) :: energy, prev_energy, step_size, &
            step_reduce, step_increase
        integer*4, intent(out) :: step_num, num_steps

        energy = m_energy
        prev_energy = m_prev_energy
        step_size = m_step_size
        step_reduce = m_step_reduce
        step_increase = m_step_increase
        step_num = m_step_num
        num_steps = m_num_steps

    end subroutine get_calc_state

    subroutine get_stretches_allocated(stretch_atom_nums_allocated, &
        stretch_reqs_allocated, stretch_keqs_allocated)
        logical, intent(out) :: stretch_atom_nums_allocated, &
            stretch_reqs_allocated, stretch_keqs_allocated

        stretch_atom_nums_allocated = m_stretch_atom_nums_allocated
        stretch_reqs_allocated = m_stretch_reqs_allocated
        stretch_keqs_allocated = m_stretch_keqs_allocated

    end subroutine get_stretches_allocated

    subroutine get_bends_allocated(bend_atom_nums_allocated, &
        bend_aeqs_allocated, bend_keqs_allocated)
        logical, intent(out) :: bend_atom_nums_allocated, &
            bend_aeqs_allocated, bend_keqs_allocated

        bend_atom_nums_allocated = m_bend_atom_nums_allocated
        bend_aeqs_allocated = m_bend_aeqs_allocated
        bend_keqs_allocated = m_bend_keqs_allocated

    end subroutine get_bends_allocated

    subroutine get_torsions_allocated(torsion_atom_nums_allocated, &
        torsions_vs_allocated, torsions_gammas_allocated, &
        torsion_paths_allocated)
        logical, intent(out) :: torsion_atom_nums_allocated, &
        torsions_vs_allocated, torsions_gammas_allocated, &
        torsion_paths_allocated

        torsion_atom_nums_allocated = m_torsion_atom_nums_allocated
        torsions_vs_allocated = m_torsions_vs_allocated
        torsions_gammas_allocated = m_torsions_gammas_allocated
        torsion_paths_allocated = m_torsion_paths_allocated

    end subroutine get_torsions_allocated

    subroutine get_non_bondings_allocated(non_bonding_atom_nums_allocated, &
        vdw_vs_allocated, vdw_reqs_allocated, q0q1s_allocated)
        logical, intent(out) :: non_bonding_atom_nums_allocated, &
        vdw_vs_allocated, vdw_reqs_allocated, q0q1s_allocated

        non_bonding_atom_nums_allocated = m_non_bonding_atom_nums_allocated
        vdw_vs_allocated = m_vdw_vs_allocated
        vdw_reqs_allocated = m_vdw_reqs_allocated
        q0q1s_allocated = m_q0q1s_allocated

    end subroutine get_non_bondings_allocated

    subroutine get_derivative_arrays(Forces, Acceleration, Velocity, Masses)
        real(dp), intent(out) :: Forces(3 * m_num_atoms), &
        Acceleration(3 * m_num_atoms), Velocity(3 * m_num_atoms), & 
        Masses(m_num_atoms)
        integer*4 :: atom_num

        Forces = m_Forces
        Acceleration = m_Acceleration
        Velocity = m_Velocity

        do atom_num = 1, m_num_atoms
            Masses(atom_num) = m_Masses(atom_num * 3)
        end do

    end subroutine get_derivative_arrays

    subroutine get_stretch_params(stretch_atom_nums, stretch_reqs, &
        stretch_keqs)
        integer*4, intent(out) :: stretch_atom_nums(m_num_stretches * 2)
        real(dp), intent(out) :: stretch_reqs(m_num_stretches), &
        stretch_keqs(m_num_stretches)

        stretch_atom_nums(:) = m_stretch_atom_nums(:)
        stretch_reqs(:) = m_stretch_reqs(:)
        stretch_keqs(:) = m_stretch_keqs(:)

    end subroutine get_stretch_params

    subroutine get_bend_params(bend_atom_nums, bend_aeqs, bend_keqs)
        integer*4, intent(out) :: bend_atom_nums(m_num_bends * 3)
        real(dp), intent(out) :: bend_aeqs(m_num_bends), &
        bend_keqs(m_num_bends)

        bend_atom_nums(:) = m_bend_atom_nums(:)
        bend_aeqs(:) = m_bend_aeqs(:)
        bend_keqs(:) = m_bend_keqs(:)

    end subroutine get_bend_params

    subroutine get_torsion_params(torsion_atom_nums, torsions_vs, &
        torsions_gammas, torsion_paths)
        integer*4, intent(out) :: torsion_atom_nums(m_num_torsions * 4)
        real(dp), intent(out) :: torsions_vs(m_num_torsions * 4), &
        torsions_gammas(m_num_torsions * 4), torsion_paths(m_num_torsions)

        torsion_atom_nums(:) = m_torsion_atom_nums(:)
        torsions_vs(:) = m_torsions_vs(:)
        torsions_gammas(:) = m_bend_keqs(:)
        torsion_paths(:) = m_torsion_paths(:)

    end subroutine get_torsion_params

    subroutine get_non_bonding_params(non_bonding_atom_nums, vdw_vs, &
        vdw_reqs, q0q1s, vdw_cutoff, diel, coul_cutoff)
        integer*4, intent(out) :: non_bonding_atom_nums(m_num_non_bonding * 2)
        real(dp), intent(out) :: vdw_vs(m_num_non_bonding), &
        vdw_reqs(m_num_non_bonding), q0q1s(m_num_non_bonding)
        real(dp), intent(out) :: vdw_cutoff, diel, coul_cutoff

        non_bonding_atom_nums(:) = m_non_bonding_atom_nums(:)
        vdw_vs(:) = m_vdw_vs(:)
        vdw_reqs(:) = m_vdw_reqs(:)
        q0q1s(:) = q0q1s(:)

        vdw_cutoff = m_vdw_cutoff
        diel = m_diel
        coul_cutoff = m_coul_cutoff

    end subroutine get_non_bonding_params
    
! This part of the code is a modified version of L-BFGS-B included in full in
! the LBFGSB_Opt directory

!*==DRIVER.spg  processed by SPAG 6.72Dc at 11:17 on  3 Dec 2018
    subroutine l_bfgs_opt(Coords, max_iterations, stretch_es_tot, &
        bend_es_tot, torsions_es_tot, vdw_es_tot, coul_es_tot, &
        convergence_threshold, precision, iprint, lim_mat_corrections, &
        status)
        
        ! precision: 1: Low
        !            2: Medium
        !            3: High
        !            4: Very High

        ! lim_mat_corrections: Limited Matrix Corrections
        !     set between 3 and 20
 
        integer*4, intent(in) :: max_iterations, precision, iprint, lim_mat_corrections
        real(dp), intent(in) :: convergence_threshold
        real(dp), intent(inout) :: Coords(m_num_atoms * 3)
        real(dp) :: stretch_es_tot(2), bend_es_tot(2), torsions_es_tot(2), &
            vdw_es_tot(2), coul_es_tot(2), norm_forces
        integer*4, intent(out) :: status

        integer*4 :: MMAX
        parameter (MMAX=17)

        character*60 :: task, csave
        logical :: lsave(4)
        integer*4 :: n, nbd(m_num_atoms * 3), iwa(m_num_atoms * 9), &
            isave(44), step_num
        real(dp) :: factr, l(m_num_atoms * 3), g(m_num_atoms * 3), &
            u(m_num_atoms * 3), dsave(29), &
            wa(6*MMAX*m_num_atoms + 12*m_num_atoms + 11*MMAX*MMAX+8*MMAX)   
   
        if (precision.eq.1) then
            factr = 1.0D+12
        else if (precision.eq.2) then
            factr = 1.0D+7
        else if (precision.eq.3) then
            factr = 1.0D+4
        else if (precision.eq.4) then
            factr = 1.0D+1
        end if
   
        n = m_num_atoms * 3

        ! Not touching l and u, so setting ndb array to 0

        nbd = 0
   
  !     We start the iteration by initializing task.
        task = 'START'
   
  !     ------- The beginning of the loop ----------
        do step_num = 1, max_iterations

            ! This is the call to the L-BFGS-B code.
            call setulb(n, lim_mat_corrections, Coords, l, u, nbd, &
                m_energy, m_Forces, factr, convergence_threshold, &
                wa, iwa, task, iprint, csave, lsave, isave, dsave)

            write(*,*) task

            if ( task(1:2).eq.'FG' ) then
    
            ! The minimization routine has returned to request the
            ! function f and gradient g values at the current x.
    
                stretch_es_tot = 0
                bend_es_tot = 0
                torsions_es_tot = 0
                vdw_es_tot = 0
                coul_es_tot = 0
                call e_amber(Coords, m_Forces, stretch_es_tot, &
                    bend_es_tot, torsions_es_tot, vdw_es_tot, coul_es_tot, &
                    norm_forces, status)

                m_energy = stretch_es_tot(1) + bend_es_tot(1) + &
                torsions_es_tot(1) + vdw_es_tot(1) + coul_es_tot(1)

            else if ( task(1:5).EQ.'NEW_X' ) then
                cycle
            else
    
    !        We terminate execution when task is neither FG nor NEW_X.
    !        We print the information contained in the string task
    !        if the default output is not used and the execution is
    !        not stopped intentionally by the user.
    
                exit
            end if
        end do
    end subroutine

end module AmberFortMod
