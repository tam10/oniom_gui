using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;

public static class Fortran {
    enum _statusEnum : int {
        OK=0, ANGLE_EXC=1, STRETCH_EXC=2, BEND_EXC=3,
        DIHEDRAL_EXC=4, NONBON_EXC1=5, NONBON_EXC2=6, NOATOMS_EXC=7,
        NOMASS_EXC=8
    }

    private const string dllName = "AmberFort";
    private const CallingConvention cc = CallingConvention.Cdecl;
    private const CharSet cs = CharSet.Ansi;

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_set_geometry", CharSet = CharSet.Ansi)
    ]
    public static extern void set_geometry([In] double[,] coords, ref int num_atoms);

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
    public static extern void e_amber(
            [Out] double[,] forces, [Out] double[] stretch_es_tot,
            [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
            [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
            ref double abs_gradient, ref int crash
    );

    [DllImport(
            dllName, CallingConvention = CallingConvention.Cdecl,
            EntryPoint = "__amberfortmod_MOD_md_verlet", CharSet = CharSet.Ansi)
    ]
    public static extern void md_verlet(
            [In, Out] double[,] Coords, ref int steps, ref double dt,
            ref double damping_factor, [Out] double[] stretch_es_tot,
            [Out] double[] bend_es_tot, [Out] double[] torsion_es_tot,
            [Out] double[] vdw_es_tot, [Out] double[] coul_es_tot,
            ref double kinetic_energy, ref double abs_gradient, ref int crash
    );

    [DllImport(
            dllName, CallingConvention = cc, 
            EntryPoint = "__amberfortmod_MOD_get_best_camera_view", CharSet = cs)
    ]
    public static extern void getBestCameraView(
            [In] double[,] coords, ref int num_Atoms, [Out] double[] camera_position, 
            ref double field_of_view_rad, ref double fill_amount, ref double min_distance
    );

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
        }
    }

}
