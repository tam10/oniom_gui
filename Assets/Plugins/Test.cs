using System;
using System.Runtime.InteropServices;
using System.Diagnostics;

public static class Test {

    static void Main() {
        //eVdW_Test();
        //nonBondingTest();
        //stretchTest();
        //bendTest();
        //torsionTest();
        //getRTest();
        hooh_test();
    }

    [DllImport(
		"AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
		EntryPoint = "__amberfortmod_MOD_set_geometry", CharSet = CharSet.Ansi)
	]
	public static extern void set_geometry(
        [In] double[,] coords, ref int num_atoms
    );

    [DllImport(
		"AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
		EntryPoint = "__amberfortmod_MOD_set_masses", CharSet = CharSet.Ansi)
	]
	public static extern void set_masses(
        [In] double[] masses, ref int crash
    );

    [DllImport(
		"AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
		EntryPoint = "__amberfortmod_MOD_set_stretches", CharSet = CharSet.Ansi)
	]
	public static extern void set_stretches(
        ref int num_stretchs, [In] int[] stretch_atom_nums,
		[In] double[] stretch_reqs, [In] double[] stretch_keqs  
    );

    [DllImport(
		"AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
		EntryPoint = "__amberfortmod_MOD_set_bends", CharSet = CharSet.Ansi)
	]
	public static extern void set_bends(
        ref int num_bends, [In] int[] bend_atom_nums,
		[In] double[] bend_aeqs, [In] double[] bend_keqs 
    );

    [DllImport(
		"AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
		EntryPoint = "__amberfortmod_MOD_set_torsions", CharSet = CharSet.Ansi)
	]
	public static extern void set_torsions(
		ref int num_torsions, [In] int[] torsion_atom_nums,
		[In] double[] torsion_vs, [In] double[] torsion_gammas, 
        [In] double[] torsionPaths
    );

    [DllImport(
		"AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
		EntryPoint = "__amberfortmod_MOD_set_non_bonding", CharSet = CharSet.Ansi)
	]
	public static extern void set_non_bonding(
        ref int num_non_bonding, [In] int[] non_bonding_atom_nums,
        [In] double[] vdw_vs, [In] double[] vdw_reqs, ref double vdw_cutoff,
        [In] double[] q0q1s, ref double diel, ref double coul_cutoff
    );

    [DllImport(
		"AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
		EntryPoint = "__amberfortmod_MOD_e_amber", CharSet = CharSet.Ansi)
	]
	public static extern void e_amber(
		[Out] double[,] forces, [Out] double[] stretch_es_tot,
		[Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
		[Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
		ref double abs_gradient, ref int crash
    );

	[DllImport(
        "AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "__amberfortmod_MOD_md_verlet", CharSet = CharSet.Ansi)
    ]
    public static extern void md_verlet(
        [In, Out] double[,] Coords, ref int steps, ref double dt,
        ref double damping_factor, [Out] double[] stretch_es_tot,
        [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
        [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
        ref double kinetic_energy, ref double abs_gradient, ref int crash
    );

    public static void hooh_test_old() {

        double[] stretchEnergies = new double[2];
        double[] bendEnergies = new double[2];
        double[] torsionEnergies = new double[2];
        double[] vdwEnergies = new double[2];
        double[] coulombicEnergies = new double[2];

		const int numAtoms = 4;
		int _numAtoms = numAtoms;

		double abs_gradient = 0f;
        int crash = 0;

        double[,] positions = new double[numAtoms, 3] {
            {-1.00000,     1.00000,   0.000000},
            {-1.00000,     0.00000,   0.000000},
            { 1.00000,     0.00000,   0.000000},
            { 1.00000,     0.00000,   1.000000}
        };
        double[,] forces = new double[numAtoms, 3];

        set_geometry(positions, ref _numAtoms);

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
            e_amber(
                forces, stretchEnergies, bendEnergies, torsionEnergies,
                vdwEnergies, coulombicEnergies,
                ref abs_gradient, ref crash
            );
        }

        Console.WriteLine(forces.Length);
        for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
            Console.WriteLine(
                "Atom: {0}. Force: {1} {2} {3}",
                atomNum,
                forces[atomNum, 0],
                forces[atomNum, 1],
                forces[atomNum, 2]
            );
        }

        Console.WriteLine("Crash        : {0}", crash);
        Console.WriteLine("Stretches ({0}): {1} {2}", _numStretches, stretchEnergies[0], stretchEnergies[1]);
        Console.WriteLine("Bends     ({0}): {1} {2}", _numBends, bendEnergies[0], bendEnergies[1]);
        Console.WriteLine("Torsions  ({0}): {1} {2}", _numTorsions, torsionEnergies[0], torsionEnergies[1]);
        Console.WriteLine("VdWs      ({0}): {1} {2}", _numNonBonding, vdwEnergies[0], vdwEnergies[1]);
        Console.WriteLine("Coulombic ({0}): {1} {2}", _numNonBonding, coulombicEnergies[0], coulombicEnergies[1]);

    }

    public static void hooh_test() {

        double[] stretchEnergies = new double[2];
        double[] bendEnergies = new double[2];
        double[] torsionEnergies = new double[2];
        double[] vdwEnergies = new double[2];
        double[] coulombicEnergies = new double[2];

		const int numAtoms = 4;
		int _numAtoms = numAtoms;

		double abs_gradient = 0f;
        int crash = 0;

        double[,] positions = new double[numAtoms, 3] {
            {-1.00000,     1.00000,   0.000000},
            {-1.00000,     0.00000,   0.000000},
            { 1.00000,     0.00000,   0.000000},
            { 1.00000,     0.00000,   1.000000}
        };

        set_geometry(positions, ref _numAtoms);

        int zero_atom = -1;
        double[] masses = new double[numAtoms] {1.0, 16.0, 16.0, 1.0};
        set_masses(masses, ref crash);
        Console.WriteLine("{0} {1}", zero_atom, crash);
		if (crash > 0) {
			throw new System.Exception(string.Format("Geometry not set. Stopping calculation"));
		}
		if (zero_atom > -1) {
			throw new System.Exception(string.Format("Atom {0} has no mass. Stopping calculation", zero_atom));
		}

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

        int stepNum = 0;
        int steps = 250;
        int updateAfterSteps = 50;
        double dt = 0.05;
        double damping_factor = 0.2;
        double kineticEnergy = 0.0;
        while (stepNum < steps) {
            md_verlet(positions, ref updateAfterSteps, ref dt, ref damping_factor,
                stretchEnergies, bendEnergies, torsionEnergies,
                vdwEnergies, coulombicEnergies,
                ref kineticEnergy, ref abs_gradient, ref crash
            );

            Console.WriteLine("Crash        : {0}", crash);
            Console.WriteLine("Stretches ({0}): {1} {2}", _numStretches, stretchEnergies[0], stretchEnergies[1]);
            Console.WriteLine("Bends     ({0}): {1} {2}", _numBends, bendEnergies[0], bendEnergies[1]);
            Console.WriteLine("Torsions  ({0}): {1} {2}", _numTorsions, torsionEnergies[0], torsionEnergies[1]);
            Console.WriteLine("VdWs      ({0}): {1} {2}", _numNonBonding, vdwEnergies[0], vdwEnergies[1]);
            Console.WriteLine("Coulombic ({0}): {1} {2}", _numNonBonding, coulombicEnergies[0], coulombicEnergies[1]);
            Console.WriteLine("Kinetic      : {0}", kineticEnergy);
            Console.WriteLine("Total        : {0}", kineticEnergy + stretchEnergies[0], bendEnergies[0],
            torsionEnergies[0], vdwEnergies[0], coulombicEnergies[0]);

            stepNum += updateAfterSteps;
        }

    }

}