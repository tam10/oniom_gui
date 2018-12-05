using System;
using System.Runtime.InteropServices;
using System.Diagnostics;

public static class Test {

    const string dllName = "AmberFort.bundle";
    const CallingConvention cc = CallingConvention.Cdecl;
    const CharSet cs = CharSet.Ansi;

    enum _statusEnum : int {
        OK=0, ANGLE_EXC=1, STRETCH_EXC=2, BEND_EXC=3,
        DIHEDRAL_EXC=4, NONBON_EXC1=5, NONBON_EXC2=6, NOATOMS_EXC=7,
        NOMASS_EXC=8
    }
    static void Main() {
        //eVdW_Test();
        //nonBondingTest();
        //stretchTest();
        //bendTest();
        //torsionTest();
        //getRTest();
        //hooh_test();
        //h2_test();
        hooh_lbfgs_test();
    }

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_allocate_atom_arrays", CharSet = cs)]
	public static extern void allocate_atom_arrays(ref int num_atoms);

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_set_masses", CharSet = cs)]
	public static extern void set_masses(
        [In] double[] masses, ref int status
    );

    [DllImport(dllName, CallingConvention = cc,EntryPoint = "__amberfortmod_MOD_set_stretches", CharSet = cs)]
	public static extern void set_stretches(
        ref int num_stretchs, [In] int[] stretch_atom_nums,
		[In] double[] stretch_reqs, [In] double[] stretch_keqs  
    );

    [DllImport(dllName, CallingConvention = cc,EntryPoint = "__amberfortmod_MOD_set_bends", CharSet = cs)]
	public static extern void set_bends(
        ref int num_bends, [In] int[] bend_atom_nums,
		[In] double[] bend_aeqs, [In] double[] bend_keqs 
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_set_torsions", CharSet = cs)]
	public static extern void set_torsions(
		ref int num_torsions, [In] int[] torsion_atom_nums,
		[In] double[] torsion_vs, [In] double[] torsion_gammas, 
        [In] double[] torsionPaths
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_set_non_bonding", CharSet = cs)]
	public static extern void set_non_bonding(
        ref int num_non_bonding, [In] int[] non_bonding_atom_nums,
        [In] double[] vdw_vs, [In] double[] vdw_reqs, ref double vdw_cutoff,
        [In] double[] q0q1s, ref double diel, ref double coul_cutoff
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_e_amber", CharSet = cs)]
	public static extern void e_amber([In] double[] Coords,
		[In, Out] double[] forces, [Out] double[] stretch_es_tot,
		[Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
		[Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
		ref double rms_forces, ref int status
    );

	[DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_md_verlet", CharSet = cs)]
    public static extern void md_verlet(
        [In, Out] double[] Coords, ref int steps, ref double dt,
        ref double damping_factor, [Out] double[] stretch_es_tot,
        [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
        [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
        ref double kinetic_energy, ref double rms_forces, ref int status
    );

	[DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_optimise", CharSet = cs)]
    public static extern void optimise(
        [In, Out] double[] Coords, ref int steps, [Out] double[] stretch_es_tot,
        [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
        [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
        ref double rms_forces, ref int status
    );

	[DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_l_bfgs_opt", CharSet = cs)]
    public static extern void lbfgs_optimise(
        [In, Out] double[] Coords, ref int steps, [Out] double[] stretch_es_tot,
        [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
        [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
        ref double convergence_threshold, ref int precision, ref int iprint, 
        ref int matrix_corrections, ref int status
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_sizes", CharSet = cs)]
    public static extern void get_sizes(ref int num_atoms, ref int num_stretches, 
        ref int num_bends, ref int num_torsions, ref int num_non_bonding
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_atoms_allocated", CharSet = cs)]
    public static extern void get_atoms_allocated(ref bool forces_allocated, 
        ref bool acceleration_allocated, ref bool velocity_allocated, 
        ref bool masses_allocated
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_calc_state", CharSet = cs)]
    public static extern void get_calc_state(ref double energy, ref double prev_energy, 
        ref double step_size, ref double step_reduce,
        ref double step_increase, ref int step_num, ref int num_steps
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_stretches_allocated", CharSet = cs)]
    public static extern void get_stretches_allocated(ref bool stretch_atom_nums_allocated, 
        ref bool stretch_reqs_allocated, ref bool stretch_keqs_allocated
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_bends_allocated", CharSet = cs)]
    public static extern void get_bends_allocated(ref bool bend_atom_nums_allocated, 
        ref bool bend_aeqs_allocated, ref bool bend_keqs_allocated
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_torsions_allocated", CharSet = cs)]
    public static extern void get_torsions_allocated(ref bool torsion_atom_nums_allocated, 
        ref bool torsions_vs_allocated, ref bool torsions_gammas_allocated, ref bool torsion_paths_allocated
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_non_bondings_allocated", CharSet = cs)]
    public static extern void get_non_bondings_allocated(ref bool non_bonding_atom_nums_allocated, 
        ref bool vdw_vs_allocated, ref bool vdw_reqs_allocated, ref bool q0q1s_allocated
    );

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_derivative_arrays", CharSet = cs)]
    private static extern void _get_derivative_arrays([Out] double[] Forces, 
        [Out] double[] Acceleration, [Out] double[] Velocity, [Out] double[] Masses
    );
    public static void get_derivative_arrays(double[] Forces, double[] Acceleration, 
        double[] Velocity, double[] Masses) {

        bool forces_allocated = false;
        bool acceleration_allocated = false;
        bool velocity_allocated = false;
        bool masses_allocated = false;

        int num_atoms = 0;
        int num_stretches = 0;
        int num_bends = 0;
        int num_torsions = 0;
        int num_non_bonding = 0;


        get_atoms_allocated(ref forces_allocated, ref acceleration_allocated,
            ref velocity_allocated, ref masses_allocated
        );

        get_sizes(ref num_atoms, ref num_stretches, 
            ref num_bends, ref num_torsions, ref num_non_bonding
        );

        _check_1D_real_array(Forces, num_atoms * 3, "Forces");
        _check_1D_real_array(Acceleration, num_atoms * 3, "Acceleration");
        _check_1D_real_array(Velocity, num_atoms * 3, "Velocity");
        _check_1D_real_array(Masses, num_atoms, "Masses");

        _check_allocated(forces_allocated, "Forces");
        _check_allocated(acceleration_allocated, "Acceleration");
        _check_allocated(velocity_allocated, "Velocity");
        _check_allocated(masses_allocated, "Masses");

        _get_derivative_arrays(Forces, Acceleration, Velocity, Masses);
    }

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_stretch_params", CharSet = cs)]
    private static extern void _get_stretch_params([Out] int[] stretch_atom_nums, 
        [Out] double[] stretch_reqs, [Out] double[] stretch_keqs
    );
    public static void get_stretch_params(int[] stretch_atom_nums, 
        double[] stretch_reqs, double[] stretch_keqs) {

        bool stretch_atom_nums_allocated = false;
        bool stretch_reqs_allocated = false;
        bool stretch_keqs_allocated = false;

        int num_atoms = 0;
        int num_stretches = 0;
        int num_bends = 0;
        int num_torsions = 0;
        int num_non_bonding = 0;

        get_stretches_allocated(ref stretch_atom_nums_allocated, 
            ref stretch_reqs_allocated, ref stretch_keqs_allocated
        );

        get_sizes(ref num_atoms, ref num_stretches, 
            ref num_bends, ref num_torsions, ref num_non_bonding
        );

        _check_1D_int_array(stretch_atom_nums, num_stretches * 2, "stretch_atom_nums");
        _check_1D_real_array(stretch_reqs, num_stretches, "stretch_reqs");
        _check_1D_real_array(stretch_keqs, num_stretches, "stretch_keqs");

        _check_allocated(stretch_atom_nums_allocated, "stretch_atom_nums");
        _check_allocated(stretch_reqs_allocated, "stretch_reqs");
        _check_allocated(stretch_keqs_allocated, "stretch_keqs");

        _get_stretch_params(stretch_atom_nums, stretch_reqs, stretch_keqs);
        
    }

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_bend_params", CharSet = cs)]
    private static extern void _get_bend_params([Out] int[] bend_atom_nums, 
        [Out] double[] bend_aeqs, [Out] double[] bend_keqs
    );
    public static void get_bend_params(int[] bend_atom_nums, 
        double[] bend_aeqs, double[] bend_keqs) {

        bool bend_atom_nums_allocated = false;
        bool bend_aeqs_allocated = false;
        bool bend_keqs_allocated = false;

        int num_atoms = 0;
        int num_stretches = 0;
        int num_bends = 0;
        int num_torsions = 0;
        int num_non_bonding = 0;

        get_bends_allocated(ref bend_atom_nums_allocated, 
            ref bend_aeqs_allocated, ref bend_keqs_allocated
        );

        get_sizes(ref num_atoms, ref num_stretches, 
            ref num_bends, ref num_torsions, ref num_non_bonding
        );

        _check_1D_int_array(bend_atom_nums, num_bends * 3, "bend_atom_nums");
        _check_1D_real_array(bend_aeqs, num_bends, "bend_aeqs");
        _check_1D_real_array(bend_keqs, num_bends, "bend_keqs");

        _check_allocated(bend_atom_nums_allocated, "bend_atom_nums");
        _check_allocated(bend_aeqs_allocated, "bend_aeqs");
        _check_allocated(bend_keqs_allocated, "bend_keqs");

        _get_bend_params(bend_atom_nums, bend_aeqs, bend_keqs);
    }

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_torsion_params", CharSet = cs)]
    private static extern void _get_torsion_params([Out] int[] torsion_atom_nums, 
        [Out] double[] torsions_vs, [Out] double[] torsions_gammas, [Out] double[] torsion_paths
    );
    public static void get_torsion_params(int[] torsion_atom_nums, 
        double[] torsions_vs, double[] torsions_gammas, double[] torsion_paths) {

        bool torsion_atom_nums_allocated = false;
        bool torsions_vs_allocated = false;
        bool torsions_gammas_allocated = false;
        bool torsion_paths_allocated = false;

        int num_atoms = 0;
        int num_stretches = 0;
        int num_bends = 0;
        int num_torsions = 0;
        int num_non_bonding = 0;

        get_torsions_allocated(ref torsion_atom_nums_allocated, ref torsions_vs_allocated, 
        ref torsions_gammas_allocated, ref torsion_paths_allocated
        );

        get_sizes(ref num_atoms, ref num_stretches, 
            ref num_bends, ref num_torsions, ref num_non_bonding
        );

        _check_1D_int_array(torsion_atom_nums, num_torsions * 4, "torsion_atom_nums");
        _check_1D_real_array(torsions_vs, num_torsions, "torsions_vs");
        _check_1D_real_array(torsions_gammas, num_torsions, "torsions_gammas");
        _check_1D_real_array(torsion_paths, num_torsions, "torsion_paths");

        _check_allocated(torsion_atom_nums_allocated, "torsion_atom_nums");
        _check_allocated(torsions_vs_allocated, "torsions_vs");
        _check_allocated(torsions_gammas_allocated, "torsions_gammas");
        _check_allocated(torsion_paths_allocated, "torsion_paths");

        _get_torsion_params(torsion_atom_nums, torsions_vs, torsions_gammas, torsion_paths);
    }

    [DllImport(dllName, CallingConvention = cc, EntryPoint = "__amberfortmod_MOD_get_non_bonding_params", CharSet = cs)]
    private static extern void _get_non_bonding_params([Out] int[] non_bonding_atom_nums, 
        [Out] double[] vdw_vs, [Out] double[] vdw_reqs, [Out] double[] q0q1s, 
        ref double vdw_cutoff, ref double diel, ref double coul_cutoff
    );
    public static void get_non_bonding_params(int[] non_bonding_atom_nums, 
        double[] vdw_vs, double[] vdw_reqs, double[] q0q1s, double vdw_cutoff, 
        double diel, double coul_cutoff) {

        bool non_bonding_atom_nums_allocated = false;
        bool vdw_vs_allocated = false;
        bool vdw_reqs_allocated = false;
        bool q0q1s_allocated = false;

        int num_atoms = 0;
        int num_stretches = 0;
        int num_bends = 0;
        int num_torsions = 0;
        int num_non_bonding = 0;

        get_non_bondings_allocated(ref non_bonding_atom_nums_allocated, 
            ref vdw_vs_allocated, ref vdw_reqs_allocated, ref q0q1s_allocated
        );

        get_sizes(ref num_atoms, ref num_stretches, 
            ref num_bends, ref num_torsions, ref num_non_bonding
        );

        _check_1D_int_array(non_bonding_atom_nums, num_non_bonding * 2, "non_bonding_atom_nums");
        _check_1D_real_array(vdw_vs, num_non_bonding, "vdw_vs");
        _check_1D_real_array(vdw_reqs, num_non_bonding, "vdw_reqs");
        _check_1D_real_array(q0q1s, num_non_bonding, "q0q1s");

        _check_allocated(non_bonding_atom_nums_allocated, "non_bonding_atom_nums");
        _check_allocated(vdw_vs_allocated, "vdw_vs");
        _check_allocated(vdw_reqs_allocated, "vdw_reqs");
        _check_allocated(q0q1s_allocated, "q0q1s");

        _get_non_bonding_params(non_bonding_atom_nums, vdw_vs, vdw_reqs, q0q1s,
            ref vdw_cutoff, ref diel, ref coul_cutoff
        );
    }





    private static void _check_2D_real_array(double[,] A, int dim0, int dim1, string name) {

        int aDim0 = A.GetLength(0);
        int aDim1 = A.GetLength(1);

        if (aDim0 != dim0 || aDim1 != dim1) {
            throw new System.Exception(string.Format(
                "Wrong size for {0}: [{1},{2}]. Should be [{3},{4}]",
                name, aDim0, aDim1, dim0, dim1
            ));
        }
    }

    private static void _check_1D_real_array(double[] A, int dim0, string name) {

        int aDim0 = A.GetLength(0);

        if (aDim0 != dim0) {
            throw new System.Exception(string.Format(
                "Wrong size for {0}: [{1}]. Should be [{2}]",
                name, aDim0, dim0
            ));
        }
    }
    
    private static void _check_1D_int_array(int[] A, int dim0, string name) {

        int aDim0 = A.GetLength(0);

        if (aDim0 != dim0) {
            throw new System.Exception(string.Format(
                "Wrong size for {0}: [{1}]. Should be [{2}]",
                name, aDim0, dim0
            ));
        }
    }

    private static void _check_allocated(bool allocated, string name) {
        if (!allocated) {
            throw new System.Exception(string.Format("Array '{}' not allocated", name));
        }
    }


    [DllImport(
        dllName, CallingConvention = cc,
        EntryPoint = "__amberfortmod_MOD_get_bad_atoms", CharSet = cs)
    ]
    public static extern void get_bad_atoms(ref int a0, ref int a1, ref int a2, ref int a3);

    public static void CheckStatus(int status) {
        if(status == (int)_statusEnum.OK) {
            return;
        }

        int a0 = 0, a1 = 0, a2 = 0, a3 = 0;
        get_bad_atoms(ref a0, ref a1, ref a2, ref a3);

        switch (status) {
            case (int)_statusEnum.ANGLE_EXC:
                throw new System.Exception( string.Format(
                    "0 distance in atoms during angle calculation: {0} {1} {2}",
                    a0, a1, a2
                ));
            case (int)_statusEnum.STRETCH_EXC:
                throw new System.Exception( string.Format(
                    "0 distance in atoms during stretch calculation: {0} {1}",
                    a0, a1
                ));
            case (int)_statusEnum.BEND_EXC:
                throw new System.Exception( string.Format(
                    "Linear angle in bend: {0} {1} {2}",
                    a0, a1, a2
                ));
            case (int)_statusEnum.DIHEDRAL_EXC:
                throw new System.Exception( string.Format(
                    "Linear angle in dihedral: {0} {1} {2} {3}",
                    a0, a1, a2, a3
                ));
            case (int)_statusEnum.NONBON_EXC1:
                throw new System.Exception( string.Format(
                    "0 distance in atoms during non-bonding calculation: {0} {1}",
                    a0, a1
                ));
            case (int)_statusEnum.NONBON_EXC2:
                throw new System.Exception( string.Format(
                    "Cannot use 0 for dielectric constant."
                ));
            case (int)_statusEnum.NOATOMS_EXC:
                throw new System.Exception( string.Format(
                    "Geometry not set; cannot assign masses."
                ));
            case (int)_statusEnum.NOMASS_EXC:
                throw new System.Exception( string.Format(
                    "Cannot use 0 for mass of atom {0}",
                    a0
                ));
            default:
                throw new System.Exception( string.Format(
                    "Fortran Error: {0}",
                    status
                ));
        }
    }



    public static void hooh_test_old() {

        double[] stretchEnergies = new double[2];
        double[] bendEnergies = new double[2];
        double[] torsionEnergies = new double[2];
        double[] vdwEnergies = new double[2];
        double[] coulombicEnergies = new double[2];

		const int numAtoms = 4;
		int _numAtoms = numAtoms;

		double abs_gradient = 0f;
        int status = 0;

        double[] positions = new double[numAtoms * 3] {
            -1.00000,     1.00000,   0.000000,
            -1.00000,     0.00000,   0.000000,
             1.00000,     0.00000,   0.000000,
             1.00000,     0.00000,   1.000000
        };
        double[] forces = new double[numAtoms * 3];

        allocate_atom_arrays(ref _numAtoms);

        //Stretches
        const int numStretches = 3;
        int _numStretches = numStretches;
        int[] stretchAtomNums = new int[numStretches * 2] {0, 1, 1, 2, 2, 3};
        double[] stretchKeqs = new double[numStretches] {553.0, 353.0, 553.0};
        double[] stretchReqs = new double[numStretches] {0.950, 1.475, 0.950};
        set_stretches(ref _numStretches, stretchAtomNums, stretchReqs, stretchKeqs);


        //Bends
        const int numBends = 2;
        int _numBends = numBends;
        int[] bendAtomNums = new int[numBends * 3] {0, 1, 2, 1, 2, 3};
        double[] bendAEqs = new double[numBends] {94.8, 94.8};
        double[] bendKEqs = new double[numBends] {60.0, 60.0};
        for (int bendIndex = 0; bendIndex < numBends; bendIndex++) {
            bendAEqs[bendIndex] = (Math.PI / 180) * bendAEqs[bendIndex];
        }
        set_bends(ref _numBends, bendAtomNums, bendAEqs, bendKEqs);

        //Torsions
        const int numTorsions = 1;
        int _numTorsions = numTorsions;
        int[] torsionAtomNums = new int[numTorsions * 4] {0, 1, 2, 3};
        double[] torsionVs = new double[numTorsions * 4] {0.0, 0.0, 1.4, 0.0};
        double[] torsionGammas = new double[numTorsions * 4] {0.0, 0.0, 0.0, 0.0};
        double[] torsionPaths = new double[numTorsions] {2.0};
        for (int torsionIndex = 0; torsionIndex < numTorsions; torsionIndex++) {
            torsionGammas[torsionIndex] = (Math.PI / 180) * torsionGammas[torsionIndex];
        }
        set_torsions(ref _numTorsions, torsionAtomNums, torsionVs, torsionGammas, torsionPaths);

        //Non bonding
        const int numNonBonding = 1;
        int _numNonBonding = numNonBonding;
        int[] nonBondingAtomNums = new int[numNonBonding * 2] {0, 3};
        double[] vdwVs = new double[numNonBonding] {0.0};
        double[] vdwReqs = new double[numNonBonding] {0.0};
        double vCutoff = 15.0;
        double cCutoff = 15.0;
        double dielectricConstant = 1.0;
        double[] q0q1s = new double[numNonBonding] {0.6000 * 0.6000};
        set_non_bonding(ref _numNonBonding, nonBondingAtomNums, 
            vdwVs, vdwReqs, ref vCutoff, 
            q0q1s, ref dielectricConstant, ref cCutoff
        );

        for (int i=0;i<2;i++ ){
            e_amber(positions, forces, 
                stretchEnergies, bendEnergies, torsionEnergies,
                vdwEnergies, coulombicEnergies,
                ref abs_gradient, ref status
            );
            CheckStatus(status);
        }

        Console.WriteLine(forces.Length);
        for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
            Console.WriteLine(
                "Atom: {0}. Force: {1} {2} {3}",
                atomNum,
                forces[atomNum * 3],
                forces[atomNum * 3 + 1],
                forces[atomNum * 3 + 2]
            );
        }

        Console.WriteLine("Crash        : {0}", status);
        Console.WriteLine("Stretches ({0}): {1} {2}", _numStretches, stretchEnergies[0], stretchEnergies[1]);
        Console.WriteLine("Bends     ({0}): {1} {2}", _numBends, bendEnergies[0], bendEnergies[1]);
        Console.WriteLine("Torsions  ({0}): {1} {2}", _numTorsions, torsionEnergies[0], torsionEnergies[1]);
        Console.WriteLine("VdWs      ({0}): {1} {2}", _numNonBonding, vdwEnergies[0], vdwEnergies[1]);
        Console.WriteLine("Coulombic ({0}): {1} {2}", _numNonBonding, coulombicEnergies[0], coulombicEnergies[1]);

    }

    public static void hooh_opt_test() {

        double[] stretchEnergies = new double[2];
        double[] bendEnergies = new double[2];
        double[] torsionEnergies = new double[2];
        double[] vdwEnergies = new double[2];
        double[] coulombicEnergies = new double[2];

		const int numAtoms = 4;
		int _numAtoms = numAtoms;

		double rms_forces = 0f;
        int status = 0;

        double[] positions = new double[numAtoms * 3] {
            -1.00000,     Math.Sqrt(2.0) * 0.5,   Math.Sqrt(2.0) * 0.5,
            -1.00000,     0.00000,   0.000000,
             1.00000,     0.00000,   0.000000,
             1.00000,     0.00000,   1.000000
        };

        allocate_atom_arrays(ref _numAtoms);

        double[] masses = new double[numAtoms] {1.0, 16.0, 16.0, 1.0};
        set_masses(masses, ref status);
        CheckStatus(status);

        //Stretches
        const int numStretches = 3;
        int _numStretches = numStretches;
        int[] stretchAtomNums = new int[numStretches * 2] {0, 1, 1, 2, 2, 3};double[] stretchKeqs = new double[numStretches] {553.0, 353.0, 553.0};
        double[] stretchReqs = new double[numStretches] {0.950, 1.475, 0.950};


        //Bends
        const int numBends = 2;
        int _numBends = numBends;
        int[] bendAtomNums = new int[numBends * 3] {0, 1, 2, 1, 2, 3};
        double[] bendAEqs = new double[numBends] {94.8, 94.8};
        double[] bendKEqs = new double[numBends] {60.0, 60.0};
        for (int bendIndex = 0; bendIndex < numBends; bendIndex++) {
            bendAEqs[bendIndex] = (Math.PI / 180) * bendAEqs[bendIndex];
        }

        //Torsions
        const int numTorsions = 1;
        int _numTorsions = numTorsions;
        int[] torsionAtomNums = new int[numTorsions * 4] {0, 1, 2, 3};
        double[] torsionVs = new double[numTorsions * 4] {0.0, 0.0, 1.4, 0.0};
        double[] torsionGammas = new double[numTorsions * 4] {0.0, 0.0, 0.0, 0.0};
        double[] torsionPaths = new double[numTorsions] {2.0};
        for (int torsionIndex = 0; torsionIndex < numTorsions; torsionIndex++) {
            torsionGammas[torsionIndex] = (Math.PI / 180) * torsionGammas[torsionIndex];
        }

        //Calc state
        double energy = 0;
        double prev_energy = 0;
        double step_size = 0;
        double step_reduce = 0;
        double step_increase = 0;
        int step_num = 0;
        int num_steps = 0;

        //Non bonding
        const int numNonBonding = 1;
        int _numNonBonding = numNonBonding;
        int[] nonBondingAtomNums = new int[numNonBonding * 2] {0, 3};
        double[] vdwVs = new double[numNonBonding] {0.0};
        double[] vdwReqs = new double[numNonBonding] {0.0};
        double vCutoff = 15.0;
        double cCutoff = 15.0;
        double dielectricConstant = 1.0;
        double[] q0q1s = new double[numNonBonding] {0.6 * 0.6};

        //_numStretches = 0;
        //_numBends = 0;
        //_numTorsions = 0;
        //_numNonBonding = 0;

        set_stretches(ref _numStretches, stretchAtomNums, stretchReqs, stretchKeqs);
        set_bends(ref _numBends, bendAtomNums, bendAEqs, bendKEqs);
        set_torsions(ref _numTorsions, torsionAtomNums, torsionVs, torsionGammas, torsionPaths);
        set_non_bonding(ref _numNonBonding, nonBondingAtomNums, 
            vdwVs, vdwReqs, ref vCutoff, 
            q0q1s, ref dielectricConstant, ref cCutoff
        );

        double[] Forces = new double[numAtoms * 3];
        double[] Acceleration = new double[numAtoms * 3];
        double[] Velocity = new double[numAtoms * 3];
        double[] Masses = new double[numAtoms];

        int stepNum = 0;
        int updateAfterSteps = 5;
        int steps = 50 * updateAfterSteps;
        while (stepNum < steps) {
            optimise(positions, ref updateAfterSteps,
                stretchEnergies, bendEnergies, torsionEnergies,
                vdwEnergies, coulombicEnergies,
                ref rms_forces, ref status
            );
            CheckStatus(status);

            get_calc_state(ref energy, ref prev_energy, ref step_size, ref step_reduce,
                ref step_increase, ref step_num, ref num_steps);

            stepNum += updateAfterSteps;

            Console.WriteLine("{0}\tTotal: {1}\tGradient: {2}",
                stepNum, 
                stretchEnergies[0] + bendEnergies[0] +
                torsionEnergies[0] + vdwEnergies[0] + coulombicEnergies[0],
                rms_forces
            );

            Console.WriteLine("m_energy: {0}\tm_prev_energy: {1}\tm_step_size: {2}",
                energy, prev_energy, step_size
            );

            get_derivative_arrays(Forces, Acceleration, Velocity, Masses);
            
            System.Text.StringBuilder sb = new System.Text.StringBuilder();
            for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
                sb.Append( string.Format("{0} ", Math.Sqrt(
                    Forces[atomNum * 3] * Forces[atomNum * 3 ] +
                    Forces[atomNum * 3 + 1] * Forces[atomNum * 3 + 1] +
                    Forces[atomNum * 3 + 2] * Forces[atomNum * 3 + 2]
                )));
            }

            Console.WriteLine("Forces: {0}\n",
                sb.ToString()
            );

        }

    }

    public static void hooh_lbfgs_test() {

        double[] stretchEnergies = new double[2];
        double[] bendEnergies = new double[2];
        double[] torsionEnergies = new double[2];
        double[] vdwEnergies = new double[2];
        double[] coulombicEnergies = new double[2];

		const int numAtoms = 4;
		int _numAtoms = numAtoms;

		double rms_forces = 0f;
        int status = 0;

        double[] positions = new double[numAtoms * 3] {
            -1.00000,     Math.Sqrt(2.0) * 0.5,   Math.Sqrt(2.0) * 0.5,
            -1.00000,     0.00000,   0.000000,
             1.00000,     0.00000,   0.000000,
             1.00000,     0.00000,   1.000000
        };

        allocate_atom_arrays(ref _numAtoms);

        double[] masses = new double[numAtoms] {1.0, 16.0, 16.0, 1.0};
        set_masses(masses, ref status);
        CheckStatus(status);

        //Stretches
        const int numStretches = 3;
        int _numStretches = numStretches;
        int[] stretchAtomNums = new int[numStretches * 2] {0, 1, 1, 2, 2, 3};double[] stretchKeqs = new double[numStretches] {553.0, 353.0, 553.0};
        double[] stretchReqs = new double[numStretches] {0.950, 1.475, 0.950};


        //Bends
        const int numBends = 2;
        int _numBends = numBends;
        int[] bendAtomNums = new int[numBends * 3] {0, 1, 2, 1, 2, 3};
        double[] bendAEqs = new double[numBends] {94.8, 94.8};
        double[] bendKEqs = new double[numBends] {60.0, 60.0};
        for (int bendIndex = 0; bendIndex < numBends; bendIndex++) {
            bendAEqs[bendIndex] = (Math.PI / 180) * bendAEqs[bendIndex];
        }

        //Torsions
        const int numTorsions = 1;
        int _numTorsions = numTorsions;
        int[] torsionAtomNums = new int[numTorsions * 4] {0, 1, 2, 3};
        double[] torsionVs = new double[numTorsions * 4] {0.0, 0.0, 1.4, 0.0};
        double[] torsionGammas = new double[numTorsions * 4] {0.0, 0.0, 0.0, 0.0};
        double[] torsionPaths = new double[numTorsions] {2.0};
        for (int torsionIndex = 0; torsionIndex < numTorsions; torsionIndex++) {
            torsionGammas[torsionIndex] = (Math.PI / 180) * torsionGammas[torsionIndex];
        }

        //Calc state
        double energy = 0;
        double prev_energy = 0;
        double step_size = 0;
        double step_reduce = 0;
        double step_increase = 0;
        int step_num = 0;
        int num_steps = 0;

        //Non bonding
        const int numNonBonding = 1;
        int _numNonBonding = numNonBonding;
        int[] nonBondingAtomNums = new int[numNonBonding * 2] {0, 3};
        double[] vdwVs = new double[numNonBonding] {0.0};
        double[] vdwReqs = new double[numNonBonding] {0.0};
        double vCutoff = 15.0;
        double cCutoff = 15.0;
        double dielectricConstant = 1.0;
        double[] q0q1s = new double[numNonBonding] {0.6 * 0.6};

        //_numStretches = 0;
        //_numBends = 0;
        //_numTorsions = 0;
        //_numNonBonding = 0;

        set_stretches(ref _numStretches, stretchAtomNums, stretchReqs, stretchKeqs);
        set_bends(ref _numBends, bendAtomNums, bendAEqs, bendKEqs);
        set_torsions(ref _numTorsions, torsionAtomNums, torsionVs, torsionGammas, torsionPaths);
        set_non_bonding(ref _numNonBonding, nonBondingAtomNums, 
            vdwVs, vdwReqs, ref vCutoff, 
            q0q1s, ref dielectricConstant, ref cCutoff
        );

        double[] Forces = new double[numAtoms * 3];
        double[] Acceleration = new double[numAtoms * 3];
        double[] Velocity = new double[numAtoms * 3];
        double[] Masses = new double[numAtoms];

        int stepNum = 0;
        int updateAfterSteps = 5;
        int steps = 50 * updateAfterSteps;

        double convergence_threshold = 0.00001;
        int precision = 3;
        int iprint = 0;
        int matrix_corrections = 17;

        while (stepNum < steps) {

            lbfgs_optimise(positions, ref updateAfterSteps,
                stretchEnergies, bendEnergies, torsionEnergies,
                vdwEnergies, coulombicEnergies, ref convergence_threshold,
                ref precision, ref iprint, ref matrix_corrections, ref status
            );
            CheckStatus(status);

            get_calc_state(ref energy, ref prev_energy, ref step_size, ref step_reduce,
                ref step_increase, ref step_num, ref num_steps);

            stepNum += updateAfterSteps;

            Console.WriteLine("{0}\tTotal: {1}\tGradient: {2}",
                stepNum, 
                stretchEnergies[0] + bendEnergies[0] +
                torsionEnergies[0] + vdwEnergies[0] + coulombicEnergies[0],
                rms_forces
            );

            Console.WriteLine("m_energy: {0}\tm_prev_energy: {1}\tm_step_size: {2}",
                energy, prev_energy, step_size
            );

            get_derivative_arrays(Forces, Acceleration, Velocity, Masses);
            
            System.Text.StringBuilder sb = new System.Text.StringBuilder();
            for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
                sb.Append( string.Format("{0} ", Math.Sqrt(
                    Forces[atomNum * 3] * Forces[atomNum * 3 ] +
                    Forces[atomNum * 3 + 1] * Forces[atomNum * 3 + 1] +
                    Forces[atomNum * 3 + 2] * Forces[atomNum * 3 + 2]
                )));
            }

            Console.WriteLine("Forces: {0}\n",
                sb.ToString()
            );

        }

    }

    public static void h2_test() {

        double[] stretchEnergies = new double[2];
        double[] bendEnergies = new double[2];
        double[] torsionEnergies = new double[2];
        double[] vdwEnergies = new double[2];
        double[] coulombicEnergies = new double[2];

		const int numAtoms = 2;
		int _numAtoms = numAtoms;

		double abs_gradient = 0f;
        int status = 0;

        double[] positions = new double[numAtoms * 3] {
            -0.50000,     0.00000,   0.000000,
             0.50000,     0.00000,   0.000000
        };

        allocate_atom_arrays(ref _numAtoms);

        double[] masses = new double[numAtoms] {1.0, 1.0};
        set_masses(masses, ref status);
        CheckStatus(status);

        //Stretches
        const int numStretches = 1;
        int _numStretches = numStretches;
        int[] stretchAtomNums = new int[numStretches * 2] {0, 1,};
        double[] stretchKeqs = new double[numStretches] {600.0};
        double[] stretchReqs = new double[numStretches] {0.740};
        set_stretches(ref _numStretches, stretchAtomNums, stretchReqs, stretchKeqs);


        //Bends
        const int numBends = 0;
        int _numBends = numBends;
        int[] bendAtomNums = new int[numBends * 3] {};
        double[] bendAEqs = new double[numBends] {};
        double[] bendKEqs = new double[numBends] {};
		set_bends(ref _numBends, bendAtomNums, bendAEqs, bendKEqs);

        //Torsions
        const int numTorsions = 0;
        int _numTorsions = numTorsions;
        int[] torsionAtomNums = new int[numTorsions * 4] {};
        double[] torsionVs = new double[numTorsions * 4] {};
        double[] torsionGammas = new double[numTorsions * 4] {};
        double[] torsionPaths = new double[numTorsions] {};
        set_torsions (
			ref _numTorsions, 
			torsionAtomNums,
			torsionVs,
			torsionGammas,
			torsionPaths
		);

        //Non bonding
        const int numNonBonding = 0;
        int _numNonBonding = numNonBonding;
        int[] nonBondingAtomNums = new int[numNonBonding * 2] {};
        double[] vdwVs = new double[numNonBonding] {};
        double[] vdwReqs = new double[numNonBonding] {};

        double vCutoff = 15.0;
        double cCutoff = 15.0;
        double dielectricConstant = 1.0;

        double[] q0q1s = new double[numNonBonding] {};
        set_non_bonding(ref _numNonBonding, nonBondingAtomNums, 
            vdwVs, vdwReqs, ref vCutoff, 
            q0q1s, ref dielectricConstant, ref cCutoff
        );

        double distance;
        int stepNum = 0;
        int steps = 200;
        int updateAfterSteps = 1;
        double dt = 0.001;
        double damping_factor = 1.0;
        double kineticEnergy = 0.0;
        while (stepNum < steps) {
            md_verlet(positions, ref updateAfterSteps, ref dt, ref damping_factor,
                stretchEnergies, bendEnergies, torsionEnergies,
                vdwEnergies, coulombicEnergies,
                ref kineticEnergy, ref abs_gradient, ref status
            );
            CheckStatus(status);

            distance = Math.Sqrt(
                Math.Pow(positions[0] - positions[3], 2) +
                Math.Pow(positions[1] - positions[4], 2) +
                Math.Pow(positions[2] - positions[5], 2)
            );

            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}",stepNum, distance, stretchEnergies[0], kineticEnergy, stretchEnergies[0] + kineticEnergy, stretchEnergies[1]);
            //Console.WriteLine("Status       : {0}", status);
            //Console.WriteLine("Stretches ({0}): {1} {2}", _numStretches, stretchEnergies[0], stretchEnergies[1]);
            //Console.WriteLine("Bends     ({0}): {1} {2}", _numBends, bendEnergies[0], bendEnergies[1]);
            //Console.WriteLine("Torsions  ({0}): {1} {2}", _numTorsions, torsionEnergies[0], torsionEnergies[1]);
            //Console.WriteLine("VdWs      ({0}): {1} {2}", _numNonBonding, vdwEnergies[0], vdwEnergies[1]);
            //Console.WriteLine("Coulombic ({0}): {1} {2}", _numNonBonding, coulombicEnergies[0], coulombicEnergies[1]);
            //Console.WriteLine("Kinetic      : {0}", kineticEnergy);
            //Console.WriteLine("Total        : {0}\n", kineticEnergy + stretchEnergies[0] + bendEnergies[0] +
            //torsionEnergies[0] + vdwEnergies[0] + coulombicEnergies[0]);

            stepNum += updateAfterSteps;
        }

    }

    public static void hooh_test() {
        double[] stretchEnergies = new double[2];
        double[] bendEnergies = new double[2];
        double[] torsionEnergies = new double[2];
        double[] vdwEnergies = new double[2];
        double[] coulombicEnergies = new double[2];

		const int numAtoms = 2;
		int _numAtoms = numAtoms;

		double abs_gradient = 0f;
        int status = 0;

        double[] positions = new double[numAtoms * 3] {
            -0.50000,     0.00000,   0.000000,
             0.50000,     0.00000,   0.000000
        };

        allocate_atom_arrays(ref _numAtoms);

        double[] masses = new double[numAtoms] {1.0, 1.0};
        set_masses(masses, ref status);
        CheckStatus(status);

        //Stretches
        const int numStretches = 1;
        int _numStretches = numStretches;
        int[] stretchAtomNums = new int[numStretches * 2] {0, 1};
        double[] stretchKeqs = new double[numStretches] {600.0};
        double[] stretchReqs = new double[numStretches] {0.740};
        set_stretches(ref _numStretches, stretchAtomNums, stretchReqs, stretchKeqs);


        //Bends
        const int numBends = 0;
        int _numBends = numBends;
        int[] bendAtomNums = new int[numBends * 3] {};
        double[] bendAEqs = new double[numBends] {};
        double[] bendKEqs = new double[numBends] {};
		set_bends(ref _numBends, bendAtomNums, bendAEqs, bendKEqs);

        //Torsions
        const int numTorsions = 0;
        int _numTorsions = numTorsions;
        int[] torsionAtomNums = new int[numTorsions * 4] {};
        double[] torsionVs = new double[numTorsions * 4] {};
        double[] torsionGammas = new double[numTorsions * 4] {};
        double[] torsionPaths = new double[numTorsions] {};
        set_torsions (
			ref _numTorsions, 
			torsionAtomNums,
			torsionVs,
			torsionGammas,
			torsionPaths
		);

        //Non bonding
        const int numNonBonding = 0;
        int _numNonBonding = numNonBonding;
        int[] nonBondingAtomNums = new int[numNonBonding * 2] {};
        double[] vdwVs = new double[numNonBonding] {};
        double[] vdwReqs = new double[numNonBonding] {};

        double vCutoff = 15.0;
        double cCutoff = 15.0;
        double dielectricConstant = 1.0;

        double[] q0q1s = new double[numNonBonding] {};
        set_non_bonding(ref _numNonBonding, nonBondingAtomNums, 
            vdwVs, vdwReqs, ref vCutoff, 
            q0q1s, ref dielectricConstant, ref cCutoff
        );

        double distance;
        int stepNum = 0;
        int steps = 200;
        int updateAfterSteps = 1;
        double dt = 0.001;
        double damping_factor = 1.0;
        double kineticEnergy = 0.0;
        while (stepNum < steps) {
            md_verlet(positions, ref updateAfterSteps, ref dt, ref damping_factor,
                stretchEnergies, bendEnergies, torsionEnergies,
                vdwEnergies, coulombicEnergies,
                ref kineticEnergy, ref abs_gradient, ref status
            );
            CheckStatus(status);

            distance = Math.Sqrt(
                Math.Pow(positions[0] - positions[3], 2) +
                Math.Pow(positions[1] - positions[4], 2) +
                Math.Pow(positions[2] - positions[5], 2)
            );

            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}",stepNum, distance, stretchEnergies[0], kineticEnergy, stretchEnergies[0] + kineticEnergy, stretchEnergies[1]);
            //Console.WriteLine("Status       : {0}", status);
            //Console.WriteLine("Stretches ({0}): {1} {2}", _numStretches, stretchEnergies[0], stretchEnergies[1]);
            //Console.WriteLine("Bends     ({0}): {1} {2}", _numBends, bendEnergies[0], bendEnergies[1]);
            //Console.WriteLine("Torsions  ({0}): {1} {2}", _numTorsions, torsionEnergies[0], torsionEnergies[1]);
            //Console.WriteLine("VdWs      ({0}): {1} {2}", _numNonBonding, vdwEnergies[0], vdwEnergies[1]);
            //Console.WriteLine("Coulombic ({0}): {1} {2}", _numNonBonding, coulombicEnergies[0], coulombicEnergies[1]);
            //Console.WriteLine("Kinetic      : {0}", kineticEnergy);
            //Console.WriteLine("Total        : {0}\n", kineticEnergy + stretchEnergies[0] + bendEnergies[0] +
            //torsionEnergies[0] + vdwEnergies[0] + coulombicEnergies[0]);

            stepNum += updateAfterSteps;
        }
    }

}