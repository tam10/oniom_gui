using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Graph : MonoBehaviour {

	public List<Connection> connections;

	//Matrix of connectivity.
	//0: Not bonded
	//1: Single
	//2: Double
	//3: Triple
	//4: Quadruple
	//5: Aromatic
	//6: Non-MM Connection
	public int[,] connectivityMatrix;
	public int[,] connectivityDistanceMatrix;
	public Atoms atoms;
	public Connection connectionPrefab;

	public int size {
		get {
			return atoms != null ? atoms.size : 0;
		}
	}

	// Update is called once per frame
	void Update () {

	}

	public void SetAtoms(Atoms atoms) {
		this.atoms = atoms;
		connections = new List<Connection> ();
		connectivityMatrix = new int[size,size];
		connectivityDistanceMatrix = new int[size,size];
	}

	public void Connect(int a0, int a1, int bondOrder=1) {

		connectivityMatrix [a0, a1] = bondOrder;
		connectivityMatrix [a1, a0] = bondOrder;

		Connection connection = Instantiate<Connection> (connectionPrefab);
		connection.Initialise (atoms, a0, a1, bondOrder);
		connections.Add (connection);
	}

	public void Disconnect(int a0, int a1) {
		connectivityMatrix [a0, a1] = 0;
		connectivityMatrix [a1, a0] = 0;
	}

	List<int> ExpandSelection(List<int> selection, int numberOfBonds, bool excludeInitialSelection=false) {

		List<int> initialSelection = new List<int>(selection);

		for (int expansionLevel = 0; expansionLevel < numberOfBonds; expansionLevel++) {
			List<int> newSelection = new List<int> ();
			foreach (int selectedAtom in selection) {
				for (int atomNum = 0; atomNum < size; atomNum++) {
					if (!selection.Contains (atomNum)) {
						if (connectivityMatrix [selectedAtom, atomNum] > 0) {
							newSelection.Add (atomNum);
						}
					}
					
				}
			}

			foreach (int newAtomNum in newSelection) {
				selection.Add (newAtomNum);
			}
		}

		if (excludeInitialSelection) {
			foreach (int oldAtomNum in initialSelection) {
				selection.Remove (oldAtomNum);
			}
		}

		return selection;

	}

	void GenerateConnectivityDistanceMatrix(int maxDistance=4) {
		for (int a0 = 0; a0 < size; a0++) {
			List<int> selection = new List<int> () { a0 };
			for (int iteration = 0; iteration < maxDistance; iteration++) {
				List<int> newSelection = ExpandSelection (selection, 1);
				foreach (int a1 in newSelection) {
					if (!selection.Contains (a1)) {
						connectivityDistanceMatrix [a0, a1] = iteration + 1;
					}
				}
			}
		}
	}
}



public class AngleConnection {

	public Atom atom0;
	public Atom atom1;
	public Atom atom2;
	public int a0;
	public int a1;
	public int a2;
	public Atoms atoms;

	public AngleConnection (Atoms atoms, int a0, int a1, int a2) {
		this.atoms = atoms;
		this.atom0 = atoms[a0];
		this.atom1 = atoms[a1];
		this.atom2 = atoms[a2];
		this.a0 = a0;
		this.a1 = a1;
		this.a2 = a2;
	}

	public Bend GetParameter() {
		foreach (Bend param in atoms.parameters.bends) {
			if (this.atom1.amberName != param.t1)
				continue;
			if (this.atom0.amberName == param.t0 && this.atom2.amberName == param.t2)
				return param;
			if (this.atom0.amberName == param.t2 && this.atom2.amberName == param.t0)
				return param;
		}
		return null;
	}

	public float EAngle(int order, bool suppress=false) {
		Bend param = GetParameter ();
		if (param == null) {
			if (!suppress)
				Debug.LogError (string.Format ("No parameter for Angle: {0}-{1}-{2}", this.atom0.amberName, this.atom1.amberName, this.atom2.amberName));
			return 0f;
		}

		float angle = this.atoms.getAngle (this.a0, this.a1, this.a2);

		if (order == 0) {
			return param.keq * Mathf.Pow (angle - param.req, 2f);
		} else if (order == 1) {
			return 2f * param.keq * (angle - param.req);
		} else if (order == 2) {
			return 2f * param.keq;
		} else {
			return 0f;
		}

	}
}

public class DihedralConnection {

	public Atom atom0;
	public Atom atom1;
	public Atom atom2;
	public Atom atom3;
	public int a0;
	public int a1;
	public int a2;
	public int a3;
	public Atoms atoms;

	public DihedralConnection (Atoms atoms, int a0, int a1, int a2, int a3) {
		this.atoms = atoms;
		this.atom0 = atoms[a0];
		this.atom1 = atoms[a1];
		this.atom2 = atoms[a2];
		this.atom3 = atoms[a3];
		this.a0 = a0;
		this.a1 = a1;
		this.a2 = a2;
		this.a3 = a3;
	}

	public Torsion GetTorsionParameter() {
		foreach (Torsion param in atoms.parameters.torsions) {
			if (
				this.atom0.amberName == param.t0 &&
				this.atom1.amberName == param.t1 &&
				this.atom2.amberName == param.t2 &&
				this.atom3.amberName == param.t3
			)
				return param;
			if (
				this.atom0.amberName == param.t3 &&
				this.atom1.amberName == param.t2 &&
				this.atom2.amberName == param.t1 &&
				this.atom3.amberName == param.t0
			)
				return param;
		}
		return null;
	}

	public ImproperTorsion GetImproperTorsionParameter() {
		foreach (ImproperTorsion param in atoms.parameters.improperTorsions) {
			if (
				this.atom0.amberName == param.t0 &&
				this.atom1.amberName == param.t1 &&
				this.atom2.amberName == param.t2 &&
				this.atom3.amberName == param.t3
			)
				return param;
			if (
				this.atom0.amberName == param.t3 &&
				this.atom1.amberName == param.t2 &&
				this.atom2.amberName == param.t1 &&
				this.atom3.amberName == param.t0
			)
				return param;
		}
		return null;
	}

	public float ETorsion(int order, bool suppress=false) {
		Torsion param = GetTorsionParameter ();
		if (param == null) {
			if (!suppress)
				Debug.LogError (string.Format ("No parameter for Torsion: {0}-{1}-{2}-{3}", this.atom0.amberName, this.atom1.amberName, this.atom2.amberName, this.atom3.amberName));
			return 0f;
		}

		float dihedral = this.atoms.getDihedral (this.a0, this.a1, this.a2, this.a3);

		float e = 0f;
		if (order == 0) {
			if (param.v0 != 0)
				e += 0.5f * param.v0 * (1f + Mathf.Cos (dihedral - param.gamma0));
			if (param.v1 != 0)
				e += 0.5f * param.v1 * (1f + Mathf.Cos (2f * dihedral - param.gamma1));
			if (param.v2 != 0)
				e += 0.5f * param.v2 * (1f + Mathf.Cos (3f * dihedral - param.gamma2));
			if (param.v3 != 0)
				e += 0.5f * param.v3 * (1f + Mathf.Cos (4f * dihedral - param.gamma3));
		} else if (order == 1) {
			if (param.v0 != 0)
				e -= 0.5f * param.v0 * (Mathf.Sin (dihedral - param.gamma0));
			if (param.v1 != 0)
				e -= 0.5f * param.v1 * 2f * (Mathf.Sin (2f * dihedral - param.gamma1));
			if (param.v2 != 0)
				e -= 0.5f * param.v2 * 3f * (Mathf.Sin (3f * dihedral - param.gamma2));
			if (param.v3 != 0)
				e -= 0.5f * param.v3 * 4f * (Mathf.Sin (4f * dihedral - param.gamma3));
		} else if (order == 2) {
			if (param.v0 != 0)
				e -= 0.5f * param.v0 * (Mathf.Cos (dihedral - param.gamma0));
			if (param.v1 != 0)
				e -= 0.5f * param.v1 * 4f * (Mathf.Cos (2f * dihedral - param.gamma1));
			if (param.v2 != 0)
				e -= 0.5f * param.v2 * 9f * (Mathf.Cos (3f * dihedral - param.gamma2));
			if (param.v3 != 0)
				e -= 0.5f * param.v3 * 16f * (Mathf.Cos (4f * dihedral - param.gamma3));
		}
		return e / param.npaths;
	}

	public float EImproperTorsion(int order, bool suppress=false) {
		ImproperTorsion param = GetImproperTorsionParameter ();
		if (param == null) {
			if (!suppress)
				Debug.LogError (string.Format ("No parameter for Improper Torsion: {0}-{1}-{2}-{3}", this.atom0.amberName, this.atom1.amberName, this.atom2.amberName, this.atom3.amberName));
			return 0f;
		}

		float dihedral = this.atoms.getDihedral (this.a0, this.a1, this.a2, this.a3);

		float e = 0f;
		if (order == 0) {
			e += 0.5f * param.v * (1f + Mathf.Cos (param.periodicity * dihedral - param.gamma));
		} else if (order == 1) {
			e -= 0.5f * param.v * param.periodicity * (Mathf.Sin (param.periodicity * dihedral - param.gamma));
		} else if (order == 2) {
			e -= 0.5f * param.v * param.periodicity * param.periodicity * (Mathf.Cos (param.periodicity * dihedral - param.gamma));
		}
		return e;

	}
}