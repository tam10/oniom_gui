using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class Graph : MonoBehaviour {

	private int maxAtoms = 1000000;
	public Dictionary<long, Connection> connectionsDict;
	private Dictionary<long, AngleConnection> anglesDict;
	private Dictionary<long, DihedralConnection> dihedralsDict;

	//POSITIONS + FORCES
	private double[] positions;
	private double[] forces;
/* 
	//STRETCHES
	private int numStretches;
	private double[] stretchREqs;
	private double[] stretchKEqs;
	private int[] stretchAtomNums;

	//BENDS
	private int numBends;
	private double[] bendAEqs;
	private double[] bendKEqs;
	private int[] bendAtomNums;

	//DIHEDRALS

	private int numTorsions;
	private double[] torsionVs;
	private double[] torsionGammas;
	private int[] torsionAtomNums;
	private double[] torsionPaths;

	//NONBONDING

	private int numNonBonding;
	private double[] vdwREqs;
	private double[] vdwVs;
	private double[] q0q1s;
	private int[] nonBondingAtomNums;

	private double amberNonbondingCutoff;
	private double amberVCutoff;
	private double amberCCutoff;*/

	private int[,] distanceMatrix;

	public int[] neighbourCount;
	private List<List<int>> adjacencyList;

	private List<int> selection;

	public List<Atom> fiveMemberedRings;
	public List<Atom> sixMemberedRings;

	public Atoms atoms;
	public GameObject connectionsHolder;
	public Connection connectionPrefab;

	//Parameters
	public Parameters parameters;
	public AmberEnvironment amberEnvironment;
	public BondTypes bondTypes;


	public double amberEnergy;
	public double allStretchesEnergy;
	public double allBendsEnergy;
	public double allTorsionsEnergy;
	public double allVdwsEnergy;
	public double allCoulombicEnergy;
	public double rmsForces;
	public double kineticEnergy;
	public double totalEnergy;


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

		amberEnvironment = Instantiate<AmberEnvironment> (GameObject.FindObjectOfType<PrefabManager>().amberEnvironmentPrefab, transform);

		connectionPrefab = Instantiate<Connection> (GameObject.FindObjectOfType<PrefabManager>().connectionPrefab, transform);

		amberEnvironment.GetEnvironmentFromXML ();
		bondTypes.GetBondTypesFromXML ();

		anglesDict = new Dictionary<long, AngleConnection>();
		dihedralsDict = new Dictionary<long, DihedralConnection>();

		selection = new List<int>();

		parameters = Instantiate<Parameters>(GameObject.FindObjectOfType<PrefabManager>().parametersPrefab, transform);
		
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

	public Connection GetConnection(int a0, int a1) {
		Connection connection;
		return connectionsDict.TryGetValue( GetConnectionKey(a0, a1), out connection) ? connection : null;
	}

	long GetConnectionKey(int a0, int a1) {
		long key = (a0 > a1) ? a0 * maxAtoms + a1 : a1 * maxAtoms + a0;
		return key;
	}

	public AngleConnection GetAngleConnection(int a0, int a1, int a2) {
		AngleConnection angleConnection;
		return anglesDict.TryGetValue( GetAngleConnectionKey(a0, a1, a2), out angleConnection) ? angleConnection : null;
	}

	long GetAngleConnectionKey(int a0, int a1, int a2) {
		long key = (a0 > a2) ? (a0 * maxAtoms + a1) * maxAtoms + a2 : (a2 * maxAtoms + a1) * maxAtoms + a0;
		return key;
	}

	long GetDihedralConnectionKey(int a0, int a1, int a2, int a3) {
		long key = (a0 > a3) ? ((a0 * maxAtoms + a1) * maxAtoms + a2) * maxAtoms + a3 : ((a3 * maxAtoms + a2) * maxAtoms + a1) * maxAtoms + a0;
		return key;
	}

	public void SetAtoms(Atoms atoms) {
		_busy++;
		this.atoms = atoms;
		connectionsDict = new Dictionary<long, Connection>();
		neighbourCount = new int[size];
		distanceMatrix = new int[size,size];

		adjacencyList = new List<List<int>>();

		for (int i=0;i<size;i++){
			adjacencyList.Add(new List<int>());
		}

		foreach (Atom atom in atoms.atomList) {
			atom.GetDefaults ();
		}

		selection = new List<int>();

	}

	public double GetBondOrder(int a0, int a1) {
		return GetConnection(a0, a1).bondOrder;
	}

	public void Connect(int a0, int a1, double bondOrder=1.0) {
		_busy++;

		if (bondOrder == 0.0) {
			Disconnect (a0, a1);
			return;
		}

		long key = GetConnectionKey(a0, a1);
		if (!connectionsDict.ContainsKey(key)) {

			adjacencyList[a0].Add(a1);
			adjacencyList[a1].Add(a0);

			Connection connection = Instantiate<Connection> (connectionPrefab);
			connection.Initialise (atoms, a0, a1, bondOrder, resolution: atoms.globalSettings.stickResolution);
			connection.transform.SetParent (connectionsHolder.transform);
			
			connectionsDict[key] = connection;

			atoms [a0].connections.Add(connection);
			atoms [a1].connections.Add(connection);
			atoms [a0].valency += bondOrder;
			atoms [a1].valency += bondOrder;
		}

		_busy--;
	}

	public void Disconnect(int a0, int a1) {
		_busy++;

		long key = GetConnectionKey(a0, a1);
		if (connectionsDict.ContainsKey(key)) {

			adjacencyList[a0].Remove(a1);
			adjacencyList[a1].Remove(a0);

			Connection connection = connectionsDict[key];
			connectionsDict.Remove(key);
			
			atoms [a0].connections.Remove(connection);
			atoms [a1].connections.Remove(connection);
			atoms [a0].valency += connection.bondOrder;
			atoms [a1].valency += connection.bondOrder;

			GameObject.Destroy(connection.gameObject);

		}
		_busy--;
	}

	public void AngleConnect(int a0, int a1, int a2) {
		_busy++;

		long key = GetAngleConnectionKey(a0, a1, a2);
		if (!anglesDict.ContainsKey(key)) {

			AngleConnection angleConnection = new AngleConnection(atoms, a0, a1, a2);
			anglesDict[key] = angleConnection;
		}

		_busy--;
	}

	public void DihedralConnect(int a0, int a1, int a2, int a3) {
		_busy++;

		long key = GetDihedralConnectionKey(a0, a1, a2, a3);
		if (!dihedralsDict.ContainsKey(key)) {

			DihedralConnection dihedralConnection = new DihedralConnection(atoms, a0, a1, a2, a3);
			dihedralsDict[key] = dihedralConnection;
		}

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

	public void ExpandSelection(int numberOfBonds=1, bool excludeInitialSelection=false) {
		_busy++;
		
		List<int> initialSelection = new List<int> (selection);
		IEnumerable<int> newSelection = new List<int>(selection);

		for (int expansionLevel = 0; expansionLevel < numberOfBonds; expansionLevel++) {
			foreach (int atomNum in selection) {
				newSelection = newSelection.Union(adjacencyList[atomNum]);
			}
			selection = newSelection.Distinct().ToList();
		}

		if (excludeInitialSelection) {
			selection = selection.Except(initialSelection).ToList();
		}
		_busy--;
	}

	public List<int> GetNeighbours(int atomNum, List<int> excludeList=null) {
		_busy++;
		List<int> neighbours = adjacencyList[atomNum];

		if (excludeList != null) {
			neighbours = neighbours.Except(excludeList).ToList();
		}

		_busy--;
		return neighbours;
	}

	public IEnumerator InitialiseParameters(bool suppressStretches=false, bool suppressBends=false, bool suppressTorsions=true, bool suppressNonBonding=false) {

		InitialiseMasses();
		InitialiseStretches(suppressStretches);
		InitialiseBends(suppressBends);
		InitialiseTorsions(suppressTorsions);
		InitialiseNonBonding(suppressNonBonding);

		yield return null;
	}

	public void InitialiseMasses() {

		int status = 0;
		Fortran.set_masses(atoms.masses, ref status);
		Fortran.CheckStatus(status);

	}

	public void InitialiseStretches(bool suppress=false) {

		Connection connection;
		Stretch stretch;
		Connection[] connections = connectionsDict.Values.ToArray();
		int numStretches = connections.Length;

		int[] stretchAtomNums = new int[numStretches * 2];
		double[] stretchREqs = new double[numStretches];
		double[] stretchKEqs = new double[numStretches];

		for (int i = 0; i < numStretches; i++) {
			connection = connections[i];
			stretch = GetStretchParameter(connection);
			if (stretch == null) {
				if (!suppress) {
				throw new NoParameterException(typeof(Stretch), connection.atom0, connection.atom1);
				}
				continue;
			}

			stretchREqs[i] = stretch.req;
			stretchKEqs[i] = stretch.keq;
			stretchAtomNums[i * 2] = connection.atom0.index;
			stretchAtomNums[i * 2 + 1] = connection.atom1.index;
		}

		Fortran.set_stretches(ref numStretches, stretchAtomNums, stretchREqs, stretchKEqs);
	}

	public void InitialiseBends(bool suppress=false) {

		AngleConnection angle;
		Bend bend;
		AngleConnection[] connections = anglesDict.Values.ToArray();
		int numBends = connections.Length;

		int[] bendAtomNums = new int[numBends * 3];
		double[] bendAEqs = new double[numBends];
		double[] bendKEqs = new double[numBends];

		for (int i = 0; i < numBends; i++) {
			angle = connections[i];
			bend = GetBendParameter(angle);
			if (bend == null) {
				if (!suppress) {
				throw new NoParameterException(typeof(Bend), angle.atom0, angle.atom1, angle.atom2);
				}
				continue;
			}

			bendAEqs[i] = Mathf.Deg2Rad * bend.req;
			bendKEqs[i] = bend.keq;
			bendAtomNums[i * 3] = angle.atom0.index;
			bendAtomNums[i * 3 + 1] = angle.atom1.index;
			bendAtomNums[i * 3 + 2] = angle.atom2.index;
		}

		Fortran.set_bends(ref numBends,	bendAtomNums, bendAEqs, bendKEqs);
	}

	public void InitialiseTorsions(bool suppress=true) {

		DihedralConnection dihedral;
		Torsion torsion;
		DihedralConnection[] dihedrals = dihedralsDict.Values.ToArray();

		List<int> torsionAtomNumList = new List<int>();
		List<double> torsionVList = new List<double>();
		List<double> torsionGammaList = new List<double>();
		List<double> torsionPathList = new List<double>();
		int numTorsions = 0;
		double npaths;

		for (int i = 0; i < dihedrals.Length; i++) {
			dihedral = dihedrals[i];
			torsion = GetTorsionParameter(dihedral);
			if (torsion == null) {
				if (!suppress) {
				throw new NoParameterException(typeof(Torsion), dihedral.atom0, dihedral.atom1, dihedral.atom2, dihedral.atom3);
				}
				continue;
			}

			torsionVList.Add(torsion.v0);
			torsionVList.Add(torsion.v1);
			torsionVList.Add(torsion.v2);
			torsionVList.Add(torsion.v3);
			torsionGammaList.Add(Mathf.Deg2Rad * torsion.gamma0);
			torsionGammaList.Add(Mathf.Deg2Rad * torsion.gamma1);
			torsionGammaList.Add(Mathf.Deg2Rad * torsion.gamma2);
			torsionGammaList.Add(Mathf.Deg2Rad * torsion.gamma3);

			if (torsion.npaths > 0.0){
				npaths = torsion.npaths;
			} else {
				//Determine number of paths from connectivity
				npaths = (double)
					(neighbourCount[dihedral.atom1.index] - 1) * 
					(neighbourCount[dihedral.atom2.index] - 1)
				;
				if (npaths <= 0.0) {
					throw new System.Exception(string.Format(
						"Error defining torsion: Atom {0} has {1} neighbours. Atom {2} has {3} neighbours",
						dihedral.atom1.index,
						neighbourCount[dihedral.atom1.index],
						dihedral.atom2.index,
						neighbourCount[dihedral.atom2.index]
						)
					);
				}
			}

			torsionPathList.Add(npaths);

			torsionAtomNumList.Add(dihedral.atom0.index);
			torsionAtomNumList.Add(dihedral.atom1.index);
			torsionAtomNumList.Add(dihedral.atom2.index);
			torsionAtomNumList.Add(dihedral.atom3.index);


			numTorsions++;
		}

		Fortran.set_torsions (
			ref numTorsions, 
			torsionAtomNumList.ToArray(),
			torsionVList.ToArray(),
			torsionGammaList.ToArray(),
			torsionPathList.ToArray()
		);
	}

	public void InitialiseNonBonding(bool suppress=false) {

		int a0;
		int a1;
		int bondDistance;

		double amberVCutoff = parameters.nonbonding.vCutoff == 0f ? atoms.globalSettings.maxNonBondingCutoff : parameters.nonbonding.vCutoff;
		double amberCCutoff = parameters.nonbonding.cCutoff == 0f ? atoms.globalSettings.maxNonBondingCutoff : parameters.nonbonding.cCutoff;

		float vScale;
		float cScale;
		float[] vScales = parameters.nonbonding.vScales;
		float[] cScales = parameters.nonbonding.cScales;

		Dictionary<string, VdW> vdwDict = new Dictionary<string, VdW>();

		foreach (VdW vdw in parameters.vdws) {
			vdwDict.Add(vdw.t, vdw);
		}

		VdW vdw0;
		VdW vdw1;
		Atom atom0;
		Atom atom1;

		int numNonBonding = size * (size - 1) / 2;
		int nonBondingIndex = 0;

		double[] q0q1s = new double[numNonBonding];
		double[] vdwREqs = new double[numNonBonding];
		double[] vdwVs = new double[numNonBonding];
		int[] nonBondingAtomNums = new int[numNonBonding * 2];

		for (a0 = 0; a0 < size - 1; a0++) {
			atom0 = atoms[a0];

			if (!vdwDict.TryGetValue(atom0.amberName, out vdw0)) {
				if (!suppress) {
					throw new NoParameterException(typeof(VdW), atom0);
				}
				for (a1 = a0 + 1; a1 < size; a1++) {
					nonBondingAtomNums[nonBondingIndex * 2] = a0;
					nonBondingAtomNums[nonBondingIndex * 2 + 1] = a1;
					nonBondingIndex++;
				}
				continue;
			}

			for (a1 = a0 + 1; a1 < size; a1++) {
				atom1 = atoms[a1];
				bondDistance = distanceMatrix[a0, a1];

				vScale = vScales[bondDistance];
				cScale = cScales[bondDistance];

				nonBondingAtomNums[nonBondingIndex * 2] = a0;
				nonBondingAtomNums[nonBondingIndex * 2 + 1] = a1;

				if (vScale != 0f) {
					if (vdwDict.TryGetValue(atom1.amberName, out vdw1)) {
						
						q0q1s[nonBondingIndex] = cScale * atom0.partialCharge * atom0.partialCharge;
						vdwREqs[nonBondingIndex] = (vdw0.r + vdw1.r) * 0.5f;
						vdwVs[nonBondingIndex] = Mathf.Sqrt (vdw0.v * vdw1.v) * vScale;
					};
				} else {
					
				}

				nonBondingIndex++;
			}

		}

		Fortran.set_non_bonding(
			ref numNonBonding,
			nonBondingAtomNums,
			vdwVs,
			vdwREqs,
			ref amberVCutoff,
			q0q1s,
			ref parameters.dielectricConstant,
			ref amberCCutoff
		);
	}

	public IEnumerator GetAmberEGH(bool suppressStretches=false, bool suppressBends=false, bool suppressTorsions=true, bool suppressNonBonding=false) {
		_busy++;

		rmsForces = 0f;
		int c;
		for (int atomNum = 0; atomNum < atoms.size; atomNum++) {
			for (c=0; c<3; c++){
				positions[atomNum * 3 + c] = atoms[atomNum].p[c];
			}
		}
		int numAtoms = size;
		forces = new double[numAtoms * 3];
		Fortran.allocate_atom_arrays(ref numAtoms);

		int status = 0;

		double[] stretchesEnergy = new double[2];
		double[] bendsEnergy = new double[2];
		double[] torsionsEnergy = new double[2];
		double[] vdwsEnergy = new double[2];
		double[] coulombicEnergy = new double[2];
		
		Fortran.e_amber(positions, forces, 
			stretchesEnergy, bendsEnergy, torsionsEnergy, vdwsEnergy,
			coulombicEnergy, ref rmsForces, ref status
		);
		Fortran.CheckStatus(status);
		
		//GetStretchesEGH();
		//yield return null;
		//GetBendsEGH();
		//yield return null;
		//GetTorsionEGH();
		///yield return null;
		//GetNonBondingEGH();
		//yield return null;

		allStretchesEnergy = stretchesEnergy[0];
		allBendsEnergy = bendsEnergy[0];
		allTorsionsEnergy = torsionsEnergy[0];
		allVdwsEnergy = vdwsEnergy[0];
		allCoulombicEnergy = coulombicEnergy[0];

		amberEnergy = allStretchesEnergy + allBendsEnergy + allTorsionsEnergy + allVdwsEnergy + allCoulombicEnergy;
		totalEnergy = amberEnergy + kineticEnergy;
		for (int atomNum = 0; atomNum < atoms.size; atomNum++) {
			for (c=0; c<3; c++){
				atoms[atomNum].force[c] = (float)forces[3 * atomNum + c];
			}
		}
		_busy--;
		yield return null;
	}

	public IEnumerator MD2(int steps) {
		_busy++;

		IEnumerator iEnumerator;
		iEnumerator = InitialiseParameters();
			while (iEnumerator.MoveNext()) {}

		foreach (Atom atom in atoms) {
			atom.mdTimeStep = atoms.globalSettings.mdTimeStep;
			atom.mdDampingFactor = atoms.globalSettings.mdDampingFactor;
			atom.force = new Vector3();
			atom.mobile = true;
		}

		yield return null;

		//Run MD
		for (int stepNum = 0; stepNum < steps; stepNum++) {

			iEnumerator = GetAmberEGH();
			while (iEnumerator.MoveNext()) {yield return new WaitForSeconds(0.1f);}

			
			yield return null;
		}
		yield return null;

		foreach (Atom atom in atoms)
			atom.mobile = false;

		_busy--;
		yield return null;
	}

	public IEnumerator MDVerlet(int steps, int updateAfterSteps=5) {
		_busy++;

		int status = 0;
		int c = 0;
		int numAtoms = size;
		positions = new double[numAtoms * 3];

		for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
			for (c=0; c<3; c++){
				positions[3 * atomNum + c] = atoms[atomNum].p[c];
			}
		}

		

		Fortran.allocate_atom_arrays(ref numAtoms);

		IEnumerator iEnumerator;
		iEnumerator = InitialiseParameters();
			while (iEnumerator.MoveNext()) {}
		yield return null;

		double[] stretchesEnergy = new double[2];
		double[] bendsEnergy = new double[2];
		double[] torsionsEnergy = new double[2];
		double[] vdwsEnergy = new double[2];
		double[] coulombicEnergy = new double[2];

		double mdTimeStep = atoms.globalSettings.mdTimeStep;
		double mdDampingFactor = atoms.globalSettings.mdDampingFactor;

		string logPath = "/Users/tristanmackenzie/Unity/ONIOM/Assets/Plugins/out.log";
		System.IO.File.Delete(logPath);

		//Run MD
		int stepNum = 0;
		Vector3 vPosition = new Vector3();
		while (stepNum < steps) {

			writeArray(positions, logPath);
			Fortran.md_verlet(positions, ref updateAfterSteps, 
				ref mdTimeStep, ref mdDampingFactor,
				stretchesEnergy, bendsEnergy, torsionsEnergy, vdwsEnergy, coulombicEnergy,
				ref kineticEnergy, ref rmsForces, ref status);
			writeArray(positions, logPath);
			Fortran.CheckStatus(status);

			allStretchesEnergy = stretchesEnergy[0];
			allBendsEnergy = bendsEnergy[0];
			allTorsionsEnergy = torsionsEnergy[0];
			allVdwsEnergy = vdwsEnergy[0];
			allCoulombicEnergy = coulombicEnergy[0];

			amberEnergy = allStretchesEnergy + allBendsEnergy + allTorsionsEnergy + allVdwsEnergy + allCoulombicEnergy;
			totalEnergy = amberEnergy + kineticEnergy;
			
			stepNum += updateAfterSteps;
			for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
				for (c = 0; c < 3; c++) {
					vPosition[c] = (float)positions[3 * atomNum + c];
				}
				atoms[atomNum].transform.localPosition = vPosition;
			}
			
			yield return null;
		}

		_busy--;
		yield return null;
	}

	public IEnumerator AMBEROpt(int steps, int updateAfterSteps=5) {
		_busy++;

		int status = 0;
		int c = 0;
		int numAtoms = size;
		positions = new double[numAtoms * 3];

		for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
			for (c=0; c<3; c++){
				positions[3 * atomNum + c] = atoms[atomNum].p[c];
			}
		}

		Fortran.allocate_atom_arrays(ref numAtoms);

		IEnumerator iEnumerator;
		iEnumerator = InitialiseParameters();
			while (iEnumerator.MoveNext()) {}
		yield return null;

		double[] stretchesEnergy = new double[2];
		double[] bendsEnergy = new double[2];
		double[] torsionsEnergy = new double[2];
		double[] vdwsEnergy = new double[2];
		double[] coulombicEnergy = new double[2];

		string logPath = "/Users/tristanmackenzie/Unity/ONIOM/Assets/Plugins/out.log";
		System.IO.File.Delete(logPath);

		double convergence_threshold = 0.00001;
		int precision = 3;
		int iprint = 0;
		int matrix_corrections = 17;

		int stepNum = 0;
		Vector3 vPosition = new Vector3();
		while (stepNum < steps) {

			writeArray(positions, logPath);
			Fortran.lbfgs_optimise(positions, ref updateAfterSteps,
                stretchesEnergy, bendsEnergy, torsionsEnergy,
                vdwsEnergy, coulombicEnergy, ref convergence_threshold,
                ref precision, ref iprint, ref matrix_corrections, ref status
            );
			writeArray(positions, logPath);
			Fortran.CheckStatus(status);

			allStretchesEnergy = stretchesEnergy[0];
			allBendsEnergy = bendsEnergy[0];
			allTorsionsEnergy = torsionsEnergy[0];
			allVdwsEnergy = vdwsEnergy[0];
			allCoulombicEnergy = coulombicEnergy[0];

			amberEnergy = allStretchesEnergy + allBendsEnergy + allTorsionsEnergy + allVdwsEnergy + allCoulombicEnergy;
			totalEnergy = amberEnergy + kineticEnergy;
			
			stepNum += updateAfterSteps;
			for (int atomNum = 0; atomNum < numAtoms; atomNum++) {
				for (c = 0; c < 3; c++) {
					vPosition[c] = (float)positions[3 * atomNum + c];
				}
				atoms[atomNum].transform.localPosition = vPosition;
			}
			
			yield return null;
		}

		_busy--;
		yield return null;
	}

	void writeArray(double[] A, string fn) {
		System.Text.StringBuilder sb = new System.Text.StringBuilder();

		for (int ln=0;ln<A.GetLength(0) / 3;ln++) {
			sb.Append(string.Format(
				"{0:00000} {1,12:0.000000} {2,12:0.000000} {3,12:0.000000}\n",
				ln,
				A[ln],
				A[ln + 1],
				A[ln + 2]
			));
		}
		sb.Append("\n");

		System.IO.File.AppendAllText(fn, sb.ToString());
	}

	public IEnumerator SolveGraph() {
		_busy++;

		int a0;
		int a1;
		int a2;
		int a3;

		int a1Num;
		int a2Num;
		int a3Num;

		List<int> a1s;
		List<int> a2s;
		List<int> a3s;

		yield return null;
		//yield return new WaitForSeconds (0.1f);

		fiveMemberedRings = new List<Atom>();
		sixMemberedRings = new List<Atom>();

		int i;
		int j;
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				distanceMatrix[i,j] = 0;
			}
		}

		for (a0 = 0; a0 < adjacencyList.Count; a0++) {

			a1s = adjacencyList[a0];
			neighbourCount[a0] = a1s.Count;
			for (a1Num = 0; a1Num < a1s.Count; a1Num++) {
				a1 = a1s[a1Num];
				distanceMatrix[a0,a1] = 1;
				distanceMatrix[a1,a0] = 1;
				
				a2s = adjacencyList[a1];
				for (a2Num = 0; a2Num < a2s.Count; a2Num++) {
					a2 = a2s[a2Num];
					if (a2 == a0)
						continue;

					if (distanceMatrix[a0,a2] == 0){
						distanceMatrix[a0,a2] = 2;
						distanceMatrix[a2,a0] = 2;
					}

					AngleConnect(a0, a1, a2);

					a3s = adjacencyList[a2];
					for (a3Num = 0; a3Num < a3s.Count; a3Num++) {
						a3 = a3s[a3Num];
						if (a3 == a1)
							continue;

						if (distanceMatrix[a0, a3] == 2) {
							fiveMemberedRings.Add(atoms[a0]);
						} else if (distanceMatrix[a0, a3] == 3) {
							sixMemberedRings.Add(atoms[a0]);
						} else if (distanceMatrix[a0,a3] == 0) {
							distanceMatrix[a0,a3] = 3;
							distanceMatrix[a3,a0] = 3;
						}
						
						DihedralConnect(a0, a1, a2, a3);

					}
				}
			}

		}

		_busy--;
	}
		
	//Environment
	public IEnumerator CalculateEnvironment(bool overwrite=true, int maxIterations=4) {
		_busy++;

		int selectionNum;
		int atomNum;
		int environmentCalculatedCount = 0;
		Atom atom;

		selection = SelectIf();

		if (!overwrite) {
			for (atomNum = 0; atomNum < size; atomNum++) {
				if (atoms[atomNum].amberName != "")
					selection.Remove(atomNum);
			}
		}

		yield return new WaitForSeconds (0.1f);

		for (int iteration = 0; iteration < maxIterations; iteration++) {
			if (environmentCalculatedCount == size)
				break;
			for (selectionNum = selection.Count - 1; selectionNum >= 0; selectionNum--){
				atomNum = selection[selectionNum];
				atom = atoms[atomNum];

				string element = atom.element;

				if (!amberEnvironment.elementEnvironments.ContainsKey(element)) {
					Debug.LogErrorFormat ("Element {0} not included in data file", element);
					selection.RemoveAt (atomNum);
					continue;
				}
				ElementEnvironment elementEnvironment = amberEnvironment.elementEnvironments [element];
				atom.mass = elementEnvironment.mass;

				bool conditionMet = false;

				foreach (string key in elementEnvironment.environmentDict.Keys) {
					AtomEnvironment atomEnvironment = elementEnvironment.environmentDict [key];

					foreach (EnvironmentCondition environmentCondition in atomEnvironment.environmentConditions) {

						//Meets condition
						if (environmentCondition.EvaluateConditions (atoms, atomNum) == 0) {
							atom.amberName = atomEnvironment.amber;
							atom.atomState = atomEnvironment.atomState;
							atom.vdwRadius = atomEnvironment.radius;
							atom.color = atomEnvironment.color;
							selection.RemoveAt (atomNum);

							conditionMet = true;
							break;
						}
					}

					if (conditionMet)
						break;
				}

			}
			yield return new WaitForSeconds (0.1f);
		}
		_busy--;
	}

	public IEnumerator CalculateConnectivity(int maxIterations) {
		_busy++;

		int selectionNum;
		int atomNum;
		int valencyUsedCount = 0;
		bondTypes.leeway = 0f;

		selection = SelectIf();

		for (selectionNum = selection.Count - 1; selectionNum >= 0; selectionNum--){
			atomNum = selection[selectionNum];
			atoms [atomNum].connectivityCollider.enabled = true;
		}

		for (int iteration = 0; iteration < maxIterations; iteration++) {

			if (valencyUsedCount == size)
				break;

			//Expand spheres to new radius
			for (selectionNum = selection.Count - 1; selectionNum >= 0; selectionNum--){
				atomNum = selection[selectionNum];

				if (atoms [atomNum].valency >= atoms [atomNum].maxValency) {
					atoms [atomNum].connectivityCollider.radius = 0.1f;
					atoms [atomNum].connectivityCollider.enabled = false;
					selection.RemoveAt(selectionNum);
				}
				atoms [atomNum].connectivityCollider.radius = atoms [atomNum].vdwRadius * 0.6f;
				
			}

			yield return new WaitForSeconds (0.1f);

			//Set existing spheres back to 0.1f to restart collisions
			for (selectionNum = selection.Count - 1; selectionNum >= 0; selectionNum--){
				atomNum = selection[selectionNum];

				atoms [atomNum].connectivityCollider.radius = 0.1f;

			}

			yield return new WaitForSeconds (0.1f);

			//Last iteration. Destroy sphere colliders
			if (iteration == maxIterations - 1) {
				for (atomNum = 0; atomNum < size; atomNum++) {
					atoms [atomNum].connectivityCollider.enabled = false;
					atoms [atomNum].connectivityCollider.radius = 0.1f;
					selection = new List<int>();
				}
			}

			bondTypes.leeway += 0.2f;

			yield return new WaitForSeconds (0.1f);

		}
		_busy--;
	}

	//Params
	public Bend GetBendParameter(AngleConnection angleConnection) {
		foreach (Bend param in parameters.bends) {
			if (angleConnection.atom1.amberName != param.t1)
				continue;
			if (angleConnection.atom0.amberName == param.t0 && angleConnection.atom2.amberName == param.t2)
				return param;
			if (angleConnection.atom0.amberName == param.t2 && angleConnection.atom2.amberName == param.t0)
				return param;
		}
		return null;
	}

	public Stretch GetStretchParameter(Connection connection) {
		foreach (Stretch param in parameters.stretches) {
			if (connection.atom0.amberName == param.t0 && connection.atom1.amberName == param.t1)
				return param;
			if (connection.atom0.amberName == param.t1 && connection.atom1.amberName == param.t0)
				return param;
		}
		return null;
	}

	public Torsion GetTorsionParameter(DihedralConnection dihedralConnection) {
		List<string> types = new List<string>();

		types.Add(dihedralConnection.atom0.amberName);
		types.Add(dihedralConnection.atom1.amberName);
		types.Add(dihedralConnection.atom2.amberName);
		types.Add(dihedralConnection.atom3.amberName);

		return parameters.torsions.Where(p => (types.SequenceEqual(p.types, new WildCardEqualityComparer()) || types.SequenceEqual(p.reverseTypes, new WildCardEqualityComparer()))).OrderBy(p => p.wildcardCount).FirstOrDefault();
	}

	public ImproperTorsion GetImproperTorsionParameter(DihedralConnection dihedralConnection) {
		List<string> types = new List<string>();

		types.Add(dihedralConnection.atom0.amberName);
		types.Add(dihedralConnection.atom1.amberName);
		types.Add(dihedralConnection.atom2.amberName);
		types.Add(dihedralConnection.atom3.amberName);
		
		return parameters.improperTorsions.Where(p => (types.SequenceEqual(p.types, new WildCardEqualityComparer()) || types.SequenceEqual(p.reverseTypes, new WildCardEqualityComparer()))).OrderBy(p => p.wildcardCount).FirstOrDefault();
	}

}


public class AngleConnection {

	public Atom atom0;
	public Atom atom1;
	public Atom atom2;

	public float angle {
		get {return Mathematics.GetAngleDeg(atom0.p, atom1.p, atom2.p);}
	}

	public AngleConnection (Atoms atoms, int a0, int a1, int a2) {
		this.atom0 = atoms[a0];
		this.atom1 = atoms[a1];
		this.atom2 = atoms[a2];
	}
}

public class DihedralConnection {

	public Atom atom0;
	public Atom atom1;
	public Atom atom2;
	public Atom atom3;

	public float dihedral {
		get {return Mathematics.GetDihedralDeg(atom0.p, atom1.p, atom2.p, atom3.p);}
	}

	public DihedralConnection (Atoms atoms, int a0, int a1, int a2, int a3) {
		this.atom0 = atoms[a0];
		this.atom1 = atoms[a1];
		this.atom2 = atoms[a2];
		this.atom3 = atoms[a3];
	}

	
}


public class NoParameterException : System.Exception {
	public NoParameterException(System.Type type, params Atom[] atoms) : base (CustomMessage(type, atoms)) {}
	private static string CustomMessage(System.Type type, params Atom[] atoms) {
		
		string[] ambers = atoms.Select(atom => atom.amberName).ToArray();
		string[] numbers = atoms.Select(atom => atom.index.ToString()).ToArray();
		string plural = ambers.Length == 1 ? "" : "s";
		return string.Format("No {0} Parameter for AMBER Type{3}: ({1}), Atom Number{3}: ({2})", type.Name, string.Join("-", ambers), string.Join("-", numbers), plural);
	}
	public NoParameterException(System.Type type, params string[] amberNames) : base (CustomMessage(type, amberNames)) {}
	private static string CustomMessage(System.Type type, params string[] amberNames) {
		
		string plural = amberNames.Length == 1 ? "" : "s";
		return string.Format("No {0} Parameter for AMBER Type{2}: ({1})", type.Name, string.Join("-", amberNames), plural);
	}
	

}
