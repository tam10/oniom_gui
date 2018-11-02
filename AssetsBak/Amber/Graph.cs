using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.LinearAlgebra.Double;


public class Graph : MonoBehaviour {

	public List<Connection> connections;

	//Matrix of connectivity.
	//0.0: Not bonded
	//1.0: Single
	//1.5: Aromatic
	//2.0: Double
	//3.0: Triple
	//4.0: Quadruple
	public int[] neighbourCount;
	private Matrix<double> adjacencyMatrix;
	public Matrix<double> connectivityMatrix;
	private Matrix<double> connectivityDistanceMatrix;

	private Vector<double> selection;

	public Vector<double> fiveMemberedRings;
	public Vector<double> sixMemberedRings;

	public Atoms atoms;
	public Connection connectionPrefab;
	public GameObject connectionsHolder;

	public AmberEnvironment amberEnvironment;
	public BondTypes bondTypes;

	private bool _showLines;
	public bool showLines;

	public int size {
		get {
			return atoms != null ? atoms.size : 0;
		}
	}

	//_busy is a counter, counting the depth of demanding scripts currently running in this object
	private int _busy;
	public bool busy {
		get { return (_busy > 0); }
	}

	// Update is called once per frame
	void Awake () {

		connectionsHolder = new GameObject ("Connections");
		connectionsHolder.transform.parent = transform;

		amberEnvironment.GetEnvironmentFromXML ();
		bondTypes.GetBondTypesFromXML ();
	}

	void Update() {
		if (atoms == null)
			return;
		if (_showLines != showLines) {
			if (showLines) 
				atoms.postRenderer.AddAtoms (atoms);
			else 
				atoms.postRenderer.RemoveAtoms (atoms);
			_showLines = showLines;
		}
	}

	public void SetAtoms(Atoms atoms) {
		_busy++;
		this.atoms = atoms;
		connections = new List<Connection> ();
		neighbourCount = new int[size];
		adjacencyMatrix = Matrix<double>.Build.Sparse (size, size);
		connectivityMatrix = Matrix<double>.Build.Sparse (size, size);
		connectivityDistanceMatrix = Matrix<double>.Build.Sparse (size, size);

		foreach (Atom atom in atoms.atomList) {
			atom.GetDefaults ();
		}

		selection = Vector<double>.Build.Sparse (size);

	}

	public double GetBondOrder(int a0, int a1) {
		return connectivityMatrix [a0, a1];
	}

	public void Connect(int a0, int a1, double bondOrder=1.0) {
		_busy++;

		if (bondOrder == 0.0) {
			Disconnect (a0, a1);
			return;
		}

		if (adjacencyMatrix [a0, a1] == 0.0 && adjacencyMatrix [a1, a0] == 0.0) {
			Connection connection = Instantiate<Connection> (connectionPrefab);
			connection.Initialise (atoms, a0, a1, bondOrder, resolution: atoms.globalSettings.stickResolution);
			connection.transform.SetParent (connectionsHolder.transform);
			connections.Add (connection);
			atoms [a0].valency += bondOrder;
			atoms [a1].valency += bondOrder;

			if (atoms.gaussianCalculator.layers.GetAtomLayerFromAtomNum(a0) == 0 || atoms.gaussianCalculator.layers.GetAtomLayerFromAtomNum(a1) == 0)
				connection.Render ();
		}

		connectivityMatrix [a0, a1] = connectivityMatrix [a1, a0] = bondOrder;
		adjacencyMatrix [a0, a1] = adjacencyMatrix [a1, a0] = 1.0;


		_busy--;
	}

	public void Disconnect(int a0, int a1) {
		_busy++;
		connectivityMatrix [a0, a1] = 0.0;
		connectivityMatrix [a1, a0] = 0.0;
		adjacencyMatrix [a0, a1] = 0.0;
		adjacencyMatrix [a1, a0] = 0.0;
		_busy--;
	}

	public List<int> SelectIf(
		List<int> selection=null,
		List<string> elements=null,
		List<string> ambers=null,
		List<string> pdbs=null,
		List<string> residues=null,
		List<int> resnums=null,
		List<int> neighbourCount=null,
		bool heavy=false
	) {
		_busy++;

		if (selection == null) {
			selection = Enumerable.Range (0, size).ToList();
		}

		IEnumerable<int> newSelection = new List<int> ();

		newSelection = selection.Where (atomNum => !heavy || atoms [atomNum].element != "H")
			.Where (atomNum => elements == null || elements.Contains (atoms [atomNum].element))
			.Where (atomNum => ambers == null || ambers.Contains (atoms [atomNum].amberName))
			.Where (atomNum => pdbs == null || pdbs.Contains (atoms [atomNum].pdbName))
			.Where (atomNum => residues == null || residues.Contains (atoms [atomNum].residueName))
			.Where (atomNum => resnums == null || resnums.Contains (atoms [atomNum].residueNumber))
			.Where (atomNum => neighbourCount == null || neighbourCount.Contains (neighbourCount [atomNum]));

		_busy--;
		return newSelection.ToList();
	}
	//instruction.AttachmentURLS = ToList().Where(Attachment => !String.IsNullOrEmpty(Attachment)).ToList();

	public void ExpandSelection(int numberOfBonds=1, bool excludeInitialSelection=false) {
		_busy++;
		
		Vector<double> initialSelection = selection.Clone();
		Vector<double> newSelection = Vector<double>.Build.Sparse(size);
		int i;

		for (int expansionLevel = 0; expansionLevel < numberOfBonds; expansionLevel++) {
			newSelection = selection * adjacencyMatrix;

			for ( i = 0; i < size; i++) {
				if (newSelection[i] != 0 && selection[i] == 0)
					selection[i] = 1;
			}
		}

		if (excludeInitialSelection) {
			for (i = 0; i < size; i++) {
				if (initialSelection [i] == 1) {
					selection [i] = 0;
				}
			}
		}
		_busy--;
	}

	public List<int> GetNeighbours(int a0, List<int> excludeList=null) {
		_busy++;


		List<int> neighbours = new List<int>();
		int index = 0;

		if (excludeList == null) {
			foreach (double n in adjacencyMatrix.Row(a0)) {
				if (n > 0.0)
					neighbours.Add (index);
				index++;
			}
		} else {
			foreach (double n in adjacencyMatrix.Row(a0)) {
				if (n > 0.0 && !excludeList.Contains(index))
					neighbours.Add (index);
				index++;
			}
		}

		_busy--;
		return neighbours;
	}

	public IEnumerator SolveGraph(int maxDistance=3) {
		_busy++;

		int i;
		int j;

		// An element of the adjacency matrix A, A[i,j] is 1 if the associated element of the connectivity matrix is non-zero, otherwise 0.
		// Each successive power n of A gives the number of paths of length n between nodes i and j.
		adjacencyMatrix = Matrix<double>.Build.Sparse (size, size);
		for (i = 0; i < size; i++) {
			for (j = i; j < size; j++) {
				if (connectivityMatrix [i, j] == 0.0) {
					adjacencyMatrix [i, j] = adjacencyMatrix [j, i] = 0.0;
				} else {
					adjacencyMatrix [i, j] = adjacencyMatrix [j, i] = 1.0;
				}
					
			}
		}

		yield return new WaitForSeconds (0.1f);

		fiveMemberedRings = Vector<double>.Build.Sparse (size);
		sixMemberedRings = Vector<double>.Build.Sparse (size);

		connectivityDistanceMatrix = adjacencyMatrix.Clone();
		Matrix<double> pows =  adjacencyMatrix.Clone ();

		// Update the connectivityDistanceMatrix iteratively, while solving for rings and neighbour count
		for (int distance = 2; distance < maxDistance + 1; distance++) {
			pows = adjacencyMatrix.Multiply (pows);

			yield return new WaitForSeconds (0.1f);

			for (i = 0; i < size - 1; i++) {
				for (j = i + 1; j < size; j++) {

					// 5-membered rings have pairs of paths of length 2 and 3 (a->b->c->d, a->e->d)
					if (distance == 3 && connectivityDistanceMatrix [i, j] == 2.0) {
						if (pows [i, j] >= 1.0) {
							fiveMemberedRings[i] = fiveMemberedRings[j] = 1.0;
						}
					}

					// pows is empty at this index, nothing to update
					if (pows [i, j] == 0.0)
						continue;
					// connectivityDistanceMatrix already populated at this index, no need to update
					if (connectivityDistanceMatrix [i, j] != 0.0)
						continue;

					//Set the connectivityDistanceMatrix at this index
					connectivityDistanceMatrix [i, j] = connectivityDistanceMatrix [j, i] = (double)distance;

					// 6-membered rings have pairs of paths of length 3
					if (distance == 3) {
						if (pows[i,j] >= 2.0) {
							sixMemberedRings[i] = sixMemberedRings[j] = 1.0;
						}
					}

				}
			}

			// Number of neighbours is the 2nd power of the Adjacency Matrix ( i -> j -> i )
			if (distance == 2) {
				for (i = 0; i < size; i++)
					neighbourCount [i] = (int)pows [i, i];
			}

		}


		_busy--;
	}
		
	//Environment
	public IEnumerator CalculateEnvironment(bool overwrite=true, int maxIterations=4) {
		_busy++;

		int atomNum;
		int environmentCalculatedCount = 0;

		for (atomNum = 0; atomNum < size; atomNum++) {
			selection [atomNum] = 1.0;
		}

		if (!overwrite) {
			for (atomNum = 0; atomNum < size; atomNum++) {
				if (atoms[atomNum].amberName != "")
					selection [atomNum] = 0.0;
			}
		}

		yield return new WaitForSeconds (0.1f);

		for (int iteration = 0; iteration < maxIterations; iteration++) {
			if (environmentCalculatedCount == size)
				break;
			for (atomNum = 0; atomNum < size; atomNum++) {
				if (selection [atomNum] == 0.0)
					continue;

				string element = atoms [atomNum].element;

				if (!amberEnvironment.elementEnvironments.ContainsKey(element)) {
					Debug.LogErrorFormat ("Element {0} not included in data file", element);
					selection [atomNum] = 0.0;
					continue;
				}
				ElementEnvironment elementEnvironment = amberEnvironment.elementEnvironments [element];

				foreach (string key in elementEnvironment.environmentDict.Keys) {
					AtomEnvironment atomEnvironment = elementEnvironment.environmentDict [key];

					foreach (EnvironmentCondition environmentCondition in atomEnvironment.environmentConditions) {

						//Meets condition
						if (environmentCondition.EvaluateConditions (atoms, atomNum) == 0) {
							atoms [atomNum].amberName = atomEnvironment.amber;
							atoms [atomNum].atomState = atomEnvironment.atomState;
							atoms [atomNum].vdwRadius = atomEnvironment.radius;
							atoms [atomNum].color = atomEnvironment.color;
							selection [atomNum] = 0.0;
							goto CONDITIONMET;

						}
					}
				}
				CONDITIONMET:;

			}
			yield return new WaitForSeconds (0.1f);
		}
		_busy--;
	}

	public IEnumerator CalculateConnectivity(int maxIterations) {
		_busy++;

		int atomNum;
		int valencyUsedCount = 0;
		bondTypes.leeway = 0f;

		for (atomNum = 0; atomNum < size; atomNum++) {
			selection [atomNum] = 1.0;
			atoms [atomNum].connectivityCollider.enabled = true;
		}

		for (int iteration = 0; iteration < maxIterations; iteration++) {

			if (valencyUsedCount == size)
				break;

			//Expand spheres to new radius
			for (atomNum = 0; atomNum < size; atomNum++) {
				if (selection [atomNum] == 0.0)
					continue;

				if (atoms [atomNum].valency >= atoms [atomNum].maxValency) {
					atoms [atomNum].connectivityCollider.enabled = false;
					selection [atomNum] = 0.0;
				}
				atoms [atomNum].connectivityCollider.radius = atoms [atomNum].vdwRadius * 0.6f;
				
			}
			//Debug.Break ();

			yield return new WaitForSeconds (0.1f);

			//Set existing spheres back to 0.1f to restart collisions
			for (atomNum = 0; atomNum < size; atomNum++) {
				if (selection [atomNum] == 0.0)
					continue;

				atoms [atomNum].connectivityCollider.radius = 0.1f;


			}

			//Debug.Break ();

			yield return new WaitForSeconds (0.1f);

			//Last iteration. Destroy sphere colliders
			if (iteration == maxIterations - 1) {
				for (atomNum = 0; atomNum < size; atomNum++) {
					atoms [atomNum].connectivityCollider.enabled = false;
					selection [atomNum] = 0.0;
				}
			}

			//Debug.Break ();
			bondTypes.leeway += 0.2f;

			yield return new WaitForSeconds (0.1f);

		}
		_busy--;
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