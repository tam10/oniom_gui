using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;

public static class Fortran {
    enum _statusEnum : int {
        OK=0, ANGLE_EXC=1, STRETCH_EXC=2, BEND_EXC=3,
        DIHEDRAL_EXC=4, NONBON_EXC1=5, NONBON_EXC2=6, NOATOMS_EXC=7,
        NOMASS_EXC=8, NPATHS_EXC=9
    }

    private const string dllName = "AmberFort";
    private const CallingConvention cc = CallingConvention.Cdecl;
    private const CharSet cs = CharSet.Ansi;

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_allocate_atom_arrays", CharSet = CharSet.Ansi)
    ]
    public static extern void allocate_atom_arrays(ref int num_atoms);

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_set_masses", CharSet = CharSet.Ansi)
    ]
    public static extern void set_masses([In] double[] masses, ref int crash);

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_get_bad_atoms", CharSet = CharSet.Ansi)
    ]
    public static extern void get_bad_atoms(ref int a0, ref int a1, ref int a2, ref int a3);

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_set_stretches", CharSet = CharSet.Ansi)
    ]
    public static extern void set_stretches(
            ref int num_stretchs, [In] int[] stretch_atom_nums,
            [In] double[] stretch_reqs, [In] double[] stretch_keqs  
    );

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_set_bends", CharSet = CharSet.Ansi)
    ]
    public static extern void set_bends(
            ref int num_bends, [In] int[] bend_atom_nums,
            [In] double[] bend_aeqs, [In] double[] bend_keqs 
    );

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_set_torsions", CharSet = CharSet.Ansi)
    ]
    public static extern void set_torsions(
            ref int num_torsions, [In] int[] torsion_atom_nums,
            [In] double[] torsion_vs, [In] double[] torsion_gammas, 
            [In] double[] torsionPaths
    );

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_set_non_bonding", CharSet = CharSet.Ansi)
    ]
    public static extern void set_non_bonding(
            ref int num_non_bonding, [In] int[] non_bonding_atom_nums,
            [In] double[] vdw_vs, [In] double[] vdw_reqs, ref double vdw_cutoff,
            [In] double[] q0q1s, ref double diel, ref double coul_cutoff
    );

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_e_amber", CharSet = CharSet.Ansi)
    ]
    public static extern void e_amber([In] double[] Coords,
		    [In, Out] double[] forces, [Out] double[] stretch_es_tot,
            [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
            [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
            ref double rms_forces, ref int status
    );

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_md_verlet", CharSet = CharSet.Ansi)
    ]
    public static extern void md_verlet(
            [In, Out] double[] Coords, ref int steps, ref double dt,
            ref double damping_factor, [Out] double[] stretch_es_tot,
            [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
            [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
            ref double kinetic_energy, ref double rms_forces, ref int status
    );

	[DllImport(
        "AmberFort.bundle", CallingConvention = CallingConvention.Cdecl,
        EntryPoint = "__amberfortmod_MOD_optimise", CharSet = CharSet.Ansi)
    ]
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

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_best_camera_view", CharSet = cs)
    ]
    public static extern void getBestCameraView(
            [In] double[] coords, ref int num_Atoms, [Out] double[] camera_position, 
            ref double field_of_view_rad, ref double fill_amount, ref double min_distance
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_deallocate_all", CharSet = cs)
    ]
    public static extern void deallocate_all();

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_sizes", CharSet = cs)
    ]
    public static extern void get_sizes(ref int num_atoms, ref int num_stretches, 
        ref int num_bends, ref int num_torsions, ref int num_non_bonding
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_atoms_allocated", CharSet = cs)
    ]
    public static extern void get_atoms_allocated(ref bool forces_allocated, 
        ref bool acceleration_allocated, ref bool velocity_allocated, 
        ref bool masses_allocated
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_calc_state", CharSet = cs)
    ]
    public static extern void get_calc_state(ref double energy, ref double prev_energy, 
        ref double step_size, ref double step_reduce,
        ref double step_increase, ref int step_num, ref int num_steps
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_stretches_allocated", CharSet = cs)
    ]
    public static extern void get_stretches_allocated(ref bool stretch_atom_nums_allocated, 
        ref bool stretch_reqs_allocated, ref bool stretch_keqs_allocated
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_bends_allocated", CharSet = cs)
    ]
    public static extern void get_bends_allocated(ref bool bend_atom_nums_allocated, 
        ref bool bend_aeqs_allocated, ref bool bend_keqs_allocated
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_torsions_allocated", CharSet = cs)
    ]
    public static extern void get_torsions_allocated(ref bool torsion_atom_nums_allocated, 
        ref bool torsions_vs_allocated, ref bool torsions_gammas_allocated, ref bool torsion_paths_allocated
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_non_bondings_allocated", CharSet = cs)
    ]
    public static extern void get_non_bondings_allocated(ref bool non_bonding_atom_nums_allocated, 
        ref bool vdw_vs_allocated, ref bool vdw_reqs_allocated, ref bool q0q1s_allocated
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_derivative_arrays", CharSet = cs)
    ]
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

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_stretch_params", CharSet = cs)
    ]
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

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_bend_params", CharSet = cs)
    ]
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

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_torsion_params", CharSet = cs)
    ]
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

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_non_bonding_params", CharSet = cs)
    ]
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
            case (int)_statusEnum.NPATHS_EXC:
                throw new System.Exception( string.Format(
                    "NPaths <= 0 for atoms: {0} {1} {2} {3}",
                    a0, a1, a2, a3
                ));
            default:
                throw new System.Exception(string.Format(
                    "Unknown error: {0}",
                    status
                ));
        }
    }

}
