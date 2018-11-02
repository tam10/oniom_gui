using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;

public class Graph : MonoBehaviour {

	private int maxAtoms = 1000000;
	public Dictionary<long, Connection> connectionsDict;
	private Dictionary<long, AngleConnection> anglesDict;
	private Dictionary<long, DihedralConnection> dihedralsDict;


	//STRETCHES
	private int numStretches;
	private float[] stretchREqs;
	private float[] stretchKEqs;
	private int[] stretchAtomNums;

	//BENDS
	private int numBends;
	private float[] bendAEqs;
	private float[] bendKEqs;
	private int[] bendAtomNums;

	//DIHEDRALS

	private int numTorsions;
	private float[] torsionVs;
	private float[] torsionGammas;
	private int[] torsionAtomNums;

	//NonBonding

	private int numNonBonding;
	private float[] vdwREqs;
	private float[] vdwVs;
	private float[] q0q1s;
	private int[] nonBondingAtomNums;


	private int numPrecomputedBondings;
	private PrecomputedNonBonding[] precomputedNonBondings;
	private float amberNonbondingCutoff;
	private float amberVCutoff;
	private float amberCCutoff;

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


	public float amberEnergy;
	public float allStretchesEnergy;
	public float allBendsEnergy;
	public float allTorsionsEnergy;
	public float allVdwsEnergy;
	public float allCoulombicEnergy;
	public float amberGradient;


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
		return connectionsDict[GetConnectionKey(a0, a1)];
	}

	long GetConnectionKey(int a0, int a1) {
		long key = (a0 > a1) ? a0 * maxAtoms + a1 : a1 * maxAtoms + a0;
		return key;
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

		InitialiseStretches(suppressStretches);
		InitialiseBends(suppressBends);
		InitialiseTorsions(suppressTorsions);
		InitialiseNonBonding(suppressNonBonding);

		yield return null;
	}

	public void InitialiseStretches(bool suppress=false) {

		Connection connection;
		Stretch stretch;
		Connection[] connections = connectionsDict.Values.ToArray();
		numStretches = connections.Length;

		stretchAtomNums = new int[numStretches * 2];
		stretchREqs = new float[numStretches];
		stretchKEqs = new float[numStretches];

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
	}

	public void InitialiseBends(bool suppress=false) {

		AngleConnection angle;
		Bend bend;
		AngleConnection[] connections = anglesDict.Values.ToArray();
		numBends = connections.Length;

		bendAtomNums = new int[numBends * 3];
		bendAEqs = new float[numBends];
		bendKEqs = new float[numBends];

		for (int i = 0; i < numBends; i++) {
			angle = connections[i];
			bend = GetBendParameter(angle);
			if (bend == null) {
				if (!suppress) {
				throw new NoParameterException(typeof(Bend), angle.atom0, angle.atom1, angle.atom2);
				}
				continue;
			}

			bendAEqs[i] = bend.req;
			bendKEqs[i] = bend.keq;
			bendAtomNums[i * 3] = angle.atom0.index;
			bendAtomNums[i * 3 + 1] = angle.atom1.index;
			bendAtomNums[i * 3 + 2] = angle.atom2.index;
		}
	}

	public void InitialiseTorsions(bool suppress=true) {

		DihedralConnection dihedral;
		Torsion torsion;
		DihedralConnection[] dihedrals = dihedralsDict.Values.ToArray();

		List<int> torsionAtomNumList = new List<int>();
		List<float> torsionVList = new List<float>();
		List<float> torsionGammaList = new List<float>();
		numTorsions = 0;

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
			torsionGammaList.Add(torsion.gamma0);
			torsionGammaList.Add(torsion.gamma1);
			torsionGammaList.Add(torsion.gamma2);
			torsionGammaList.Add(torsion.gamma3);
			
			torsionAtomNumList.Add(dihedral.atom0.index);
			torsionAtomNumList.Add(dihedral.atom1.index);
			torsionAtomNumList.Add(dihedral.atom2.index);
			torsionAtomNumList.Add(dihedral.atom3.index);

			numTorsions++;
		}

		torsionAtomNums = torsionAtomNumList.ToArray();
		torsionVs = torsionVList.ToArray();
		torsionGammas = torsionGammaList.ToArray();
	}

	public void InitialiseNonBonding(bool suppress=false) {

		int a0;
		int a1;
		int bondDistance;
		numPrecomputedBondings = 0;

		amberVCutoff = parameters.nonbonding.vCutoff == 0f ? atoms.globalSettings.maxNonBondingCutoff : parameters.nonbonding.vCutoff;
		amberCCutoff = parameters.nonbonding.cCutoff == 0f ? atoms.globalSettings.maxNonBondingCutoff : parameters.nonbonding.cCutoff;

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

		numNonBonding = size * (size - 1) / 2;
		int nonBondingIndex = 0;

		q0q1s = new float[numNonBonding];
		vdwREqs = new float[numNonBonding];
		vdwVs = new float[numNonBonding];
		nonBondingAtomNums = new int[numNonBonding * 2];

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
	}


	public void GetStretchesEGH() {
		float[] energies = new float[3];
		Atom atom0;
		Atom atom1;
		Vector3 v10;
		float r10;
		Vector3 force;
		allStretchesEnergy = 0f;
		for (int i = 0; i < numStretches; i++) {

			atom0 = atoms[stretchAtomNums[i * 2]];
			atom1 = atoms[stretchAtomNums[i * 2 + 1]];
			v10 = atom1.p - atom0.p;
			r10 = v10.magnitude;

			Mathematics.EStretch(
				r10,
				stretchKEqs[i],
				stretchREqs[i],
				energies
			);

			allStretchesEnergy += energies[0];
			amberGradient += Mathf.Abs(energies[1]);

			force = v10 * energies[1] / r10;
			atom0.force += force;
			atom1.force -= force;
		}
	}

	public void GetBendsEGH() {
		float[] energies = new float[3];
		Atom atom0;
		Atom atom1;
		Atom atom2;

		Vector3 perp;
		Vector3 v10;
		Vector3 v12;

		float r10;
		float r12;

		Vector3 force10;
		Vector3 force12;

		allBendsEnergy = 0f;
		for (int i = 0; i < numBends; i++) {

			atom0 = atoms[bendAtomNums[i * 3]];
			atom1 = atoms[bendAtomNums[i * 3 + 1]];
			atom2 = atoms[bendAtomNums[i * 3 + 2]];
			v10 = (atom1.p - atom0.p);
			v12 = (atom1.p - atom2.p);
			r10 = v10.magnitude;
			r12 = v12.magnitude;

			Mathematics.EAngle(
				Mathf.Deg2Rad * (Vector3.Angle(v10, v12) - bendAEqs[i]),
				bendKEqs[i],
				energies
			);

			allBendsEnergy += energies[0];
			amberGradient += Mathf.Abs(energies[1]);

			perp = Vector3.Cross(v10, v12);

			force10 = Vector3.Cross(v10, perp) * energies[1] / r10;
			force12 = - Vector3.Cross(v12, perp) * energies[1] / r12;

			atom0.force += force10;
			atom1.force -= force10 + force12;
			atom2.force += force12;

		}
	}

	public void GetTorsionEGH() {
		float[] energies = new float[3];
		Atom atom0;
		Atom atom1;
		Atom atom2;
		Atom atom3;

		Vector3 v10;
		Vector3 v12;
		Vector3 v23;
		Vector3 w1;
		Vector3 w2;
		Vector3 c;
		Vector3 vc2;

		float dihedral;
		float r10;
		float r23;
		float rc2;
		float da_dr10;
		float da_dr23;
		float s1;
		float s2;

		Vector3 force0;
		Vector3 force1;
		Vector3 force2;
		Vector3 force3;

		allTorsionsEnergy = 0f;
		int torsionIndex;
		for (int i = 0; i < numTorsions; i++) {
			torsionIndex = i * 4;
			//if (i!=55)
			//	continue;
			atom0 = atoms[torsionAtomNums[i * 4]];
			atom1 = atoms[torsionAtomNums[i * 4 + 1]];
			atom2 = atoms[torsionAtomNums[i * 4 + 2]];
			atom3 = atoms[torsionAtomNums[i * 4 + 3]];

			v10 = (atom0.p - atom1.p);
			v12 = (atom2.p - atom1.p).normalized;
			v23 = (atom3.p - atom2.p);

			r10 = v10.magnitude;
			r23 = v23.magnitude;

			w1 = Vector3.Cross(v10 / r10, v12);
			w2 = Vector3.Cross(v23 / r23, v12);

			//dihedral = Vector3.Angle(w1, w2);
			//if (Vector3.Cross(w1, w2).y > 0) dihedral = -dihedral;

			dihedral = Mathf.Atan2(Vector3.Dot(v12, Vector3.Cross(w1, w2)), Vector3.Dot(w1, w2));

			s1 = Mathf.Sin(Mathf.Deg2Rad * Vector3.Angle(v10, v12));
			s2 = Mathf.Sin(Mathf.Deg2Rad * Vector3.Angle(- v12, v23));
			da_dr10 = 1f / (r10 * s1);
			da_dr23 = 1f / (r23 * s2);

			Mathematics.ETorsion(
				dihedral,
				torsionVs[torsionIndex],
				torsionVs[torsionIndex + 1],
				torsionVs[torsionIndex + 2],
				torsionVs[torsionIndex + 3],
				torsionGammas[torsionIndex],
				torsionGammas[torsionIndex + 1],
				torsionGammas[torsionIndex + 2],
				torsionGammas[torsionIndex + 3],
				energies
			);

			allTorsionsEnergy += energies[0];
			amberGradient += Mathf.Abs(energies[1]);

			force0 = - w1 * da_dr10 * energies[1];
			force3 = w2 * da_dr23 * energies[1];

			c = (atom1.p + atom2.p) * 0.5f;
			vc2 = c - atom3.p;
			rc2 = vc2.magnitude;

			force2 = Vector3.Cross( - ((Vector3.Cross(vc2, force3) + 0.5f * Vector3.Cross(v23, force3) - 0.5f * Vector3.Cross(v12, force0)) / (rc2 * rc2)), vc2);
			force1 = - force0 - force2 - force3;

			atom0.force += force0;
			atom1.force += force1;
			atom2.force += force2;
			atom3.force += force3;

			//Debug.LogFormat("{0} {1} {2} {3} {4}", i, atom0.index, atom1.index, atom2.index,atom3.index);
			//if (i==55){
			/* if (false){
				Debug.LogFormat("F: {0}, E: {1}, G: {2}, D: {3}", force0, energies[0], energies[1], dihedral);

				Debug.LogFormat(
					"{0} {1} {2} {3}   {4} {5} {6} {7}   {8} {9} {10} {11}",
					atom0.amberName,
					atom1.amberName,
					atom2.amberName,
					atom3.amberName,
					torsionVs[torsionIndex],
					torsionVs[torsionIndex + 1],
					torsionVs[torsionIndex + 2],
					torsionVs[torsionIndex + 3],
					torsionGammas[torsionIndex],
					torsionGammas[torsionIndex + 1],
					torsionGammas[torsionIndex + 2],
					torsionGammas[torsionIndex + 3]
				); 

				Debug.DrawLine(atom0.p, atom0.p + force0, Color.red);
				Debug.DrawLine(atom3.p, atom3.p + force3, Color.red);

				Debug.DrawLine(atom1.p, atom1.p + w1, Color.white);
				Debug.DrawLine(atom2.p, atom2.p + w2, Color.white);

				Debug.DrawLine(atom1.p, atom1.p + v10, Color.blue);
				Debug.DrawLine(atom1.p, atom1.p + v12, Color.blue);
				Debug.DrawLine(atom2.p, atom2.p + v23, Color.blue);
				//Debug.DrawLine(atom1.p, atom1.p + force1, Color.green);
				//Debug.DrawLine(atom2.p, atom2.p + force2, Color.blue);
			}
			*/
		}
	}

	public void GetNonBondingEGH() {
		float[] vdwEnergies = new float[3];
		float[] coulombicEnergies = new float[3];
		Atom atom0;
		Atom atom1;
		Vector3 v10;
		float r10;
		Vector3 force;
		

		float dielectricConstant = parameters.dielectricConstant;

		allCoulombicEnergy = 0f;
		allVdwsEnergy = 0f;
		for (int i = 0; i < numNonBonding; i++) {

			atom0 = atoms[nonBondingAtomNums[i * 2]];
			atom1 = atoms[nonBondingAtomNums[i * 2 + 1]];
			v10 = atom1.p - atom0.p;
			r10 = v10.magnitude;

			Mathematics.EElectrostaticR1(
				r10,
				q0q1s[i],
				dielectricConstant,
				coulombicEnergies
			);

			Mathematics.EVdWAmber(
				r10,
				vdwVs[i],
				vdwREqs[i],
				vdwEnergies
			);


			allCoulombicEnergy += coulombicEnergies[0];
			allVdwsEnergy += vdwEnergies[0];

			amberGradient += Mathf.Abs(coulombicEnergies[1]) + Mathf.Abs(vdwEnergies[1]);

			force = v10 * (coulombicEnergies[1] + vdwEnergies[1]) / r10;
			atom0.force += force;
			atom1.force -= force;
		}
	}

	public IEnumerator GetAmberEGH(bool suppressStretches=false, bool suppressBends=false, bool suppressTorsions=true, bool suppressNonBonding=false) {
		_busy++;

		amberGradient = 0f;
		GetStretchesEGH();
		GetBendsEGH();
		GetTorsionEGH();
		GetNonBondingEGH();

		amberEnergy = allStretchesEnergy + allBendsEnergy + allTorsionsEnergy + allVdwsEnergy + allCoulombicEnergy;

		_busy--;
		yield return null;
	}

	public IEnumerator MD2(int steps) {
		_busy++;

		IEnumerator iEnumerator;
		iEnumerator = InitialiseParameters();
			while (iEnumerator.MoveNext()) {}

		foreach (Atom atom in atoms)
			atom.mobile = true;

		yield return null;

		//Run MD
		float mdDampingFactor;
		for (int stepNum = 0; stepNum < steps; stepNum++) {

			mdDampingFactor = atoms.globalSettings.mdDampingFactor;
			foreach (Atom atom in atoms) {
				atom.force = Vector3.zero;
				atom.velocity *= mdDampingFactor;
			}

			iEnumerator = GetAmberEGH();
			while (iEnumerator.MoveNext()) {}

			yield return null;
		}
		yield return null;

		foreach (Atom atom in atoms)
			atom.mobile = false;

		_busy--;
		yield return null;
	}

	public IEnumerator MD(int steps, bool suppress=false) {
		_busy++;

		//Initialise Non-bonding
		List<PrecomputedNonBonding> precomputedNonBondingsList = new List<PrecomputedNonBonding>();

		int a0;
		int a1;
		int bondDistance;
		numPrecomputedBondings = 0;

		amberNonbondingCutoff = parameters.nonbonding.vCutoff == 0f ? atoms.globalSettings.maxNonBondingCutoff : parameters.nonbonding.vCutoff;
		float vScale1 = parameters.nonbonding.vScales[1];
		float vScale2 = parameters.nonbonding.vScales[2];
		float vScale3 = parameters.nonbonding.vScales[3];
		float vScale;
		float cScale1 = parameters.nonbonding.cScales[1];
		float cScale2 = parameters.nonbonding.cScales[2];
		float cScale3 = parameters.nonbonding.cScales[3];
		float cScale;
		float dielectricConstant = parameters.dielectricConstant;

		Dictionary<string, VdW> vdwDict = new Dictionary<string, VdW>();

		foreach (VdW vdw in parameters.vdws) {
			vdwDict.Add(vdw.t, vdw);
		}

		VdW vdw0;
		VdW vdw1;
		Atom atom0;
		Atom atom1;

		for (a0 = 0; a0 < size - 1; a0++) {
			atom0 = atoms[a0];

			if (!vdwDict.TryGetValue(atom0.amberName, out vdw0)) {
				if (!suppress)
					throw new NoParameterException(typeof(VdW), atom0);
				continue;
			}

			for (a1 = a0 + 1; a1 < size; a1++) {
				atom1 = atoms[a1];

				bondDistance = distanceMatrix[a0, a1];

				if (bondDistance == 1) {
					vScale = vScale1;
					cScale = cScale1;
				} else if (bondDistance == 2) {
					vScale = vScale2;
					cScale = cScale2;
				} else if (bondDistance == 3) {
					vScale = vScale3;
					cScale = cScale3;
				} else {
					vScale = 1f;
					cScale = 1f;
				}

				if (!vdwDict.TryGetValue(atoms[a1].amberName, out vdw1)) {
					if (!suppress)
					throw new NoParameterException(typeof(VdW), atom1);
					continue;
				}

				if (vScale == 0f)
					continue;

				PrecomputedNonBonding precomputedNonBonding = new PrecomputedNonBonding(vdw0, vdw1, atom0, atom1, vScale, cScale, dielectricConstant);
				//If we're ignoring errors, don't add this param

				precomputedNonBondingsList.Add(precomputedNonBonding);
				numPrecomputedBondings++;
			}
		}

		precomputedNonBondings = precomputedNonBondingsList.ToArray();

		yield return null;

		//Mobilise atoms
		foreach (Atom atom in atoms)
			atom.mobile = true;

		yield return null;

		//Run MD
		for (int stepNum = 0; stepNum < steps; stepNum++) {
			//IEnumerator iEnumerator = GetAmberEnergy(0,true);
			//while (iEnumerator.MoveNext()) {}

			Coroutine coroutine = StartCoroutine( MDStep(suppress));

			Debug.Log(amberEnergy);
			yield return null;
			yield return coroutine;
		}
		yield return null;

		foreach (Atom atom in atoms)
			atom.mobile = false;
		_busy--;
	}

	public IEnumerator MDStep(bool suppress=false) {
		_busy++;

		foreach (Atom atom in atoms)
			atom.force = Vector3.zero;
		
		IEnumerator iEnumerator;

		//Stretches
		iEnumerator = GetStretchesForces(suppress);
		while (iEnumerator.MoveNext()) {}

		//Bends
		iEnumerator = GetBendsForces(suppress);
		while (iEnumerator.MoveNext()) {}

		//Non bonding
		iEnumerator = GetNonBondingForces(suppress);
		while (iEnumerator.MoveNext()) {}

		yield return null;

		_busy--;
	}

	public IEnumerator GetStretchesForces(bool suppress=false) {
		_busy++;

		float dE_dr;
		Vector3 connectionVector;

		foreach (Connection connection in connectionsDict.Values) {
			dE_dr = connection.EStretch(1, suppress);
			connectionVector = (connection.atom1.p - connection.atom0.p).normalized;
			connection.atom0.force += connectionVector * dE_dr;
			connection.atom1.force -= connectionVector * dE_dr;

		}

		_busy--;
		yield return null;
	}

	public IEnumerator GetBendsForces(bool suppress=false) {
		_busy++;

		float dE_da;
		Vector3 connectionVector10;
		Vector3 connectionVector12;
		Vector3 perp;

		foreach (AngleConnection angle in anglesDict.Values) {
			dE_da = EAngle(angle, 1, suppress);

			connectionVector10 = (angle.atom1.p - angle.atom0.p).normalized;
			connectionVector12 = (angle.atom1.p - angle.atom2.p).normalized;
			perp = Vector3.Cross(connectionVector10, connectionVector12);

			angle.atom0.force += Vector3.Cross(connectionVector10, perp) * dE_da;
			angle.atom2.force -= Vector3.Cross(connectionVector12, perp) * dE_da;

		}

		_busy--;
		yield return null;
	}

	public IEnumerator GetNonBondingForces(bool suppress=false) {
		_busy++;
		Vector3 connectionVector;

		int nonBondingIndex;
		PrecomputedNonBonding precomputedNonBonding;
		float[] vdwEnergies = new float[3];
		float[] coulombEnergies  = new float[3];
		float r;

		allVdwsEnergy = 0f;
		allCoulombicEnergy = 0f;

		for (nonBondingIndex = 0; nonBondingIndex < numPrecomputedBondings; nonBondingIndex++) {
			precomputedNonBonding = precomputedNonBondings[nonBondingIndex];
			r = Vector3.Distance(precomputedNonBonding.atom0.p, precomputedNonBonding.atom1.p);
			if (r > amberNonbondingCutoff) {
				continue;
			}

			connectionVector = (precomputedNonBonding.atom1.p - precomputedNonBonding.atom0.p).normalized;
			
			precomputedNonBonding.GetEnergy(r, atoms.globalSettings.angstromToBohr, 1, vdwEnergies, coulombEnergies);

			allVdwsEnergy += vdwEnergies[0];
			allCoulombicEnergy += coulombEnergies[0];

			precomputedNonBonding.atom0.force += connectionVector * (vdwEnergies[1] + coulombEnergies[1]);
			precomputedNonBonding.atom1.force -= connectionVector * (vdwEnergies[1] + coulombEnergies[1]);
			
		}

		_busy--;
		yield return null;
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

	public float EAngle(AngleConnection angleConnection, int order, bool suppress=false) {
		Bend param = GetBendParameter (angleConnection);
		if (param == null) {
			if (!suppress)
				Debug.LogError (string.Format ("No parameter for Angle: {0}-{1}-{2}", angleConnection.atom0.amberName, angleConnection.atom1.amberName, angleConnection.atom2.amberName));
			return 0f;
		}

		float angle = Mathematics.GetAngleDeg(angleConnection.atom0.p, angleConnection.atom1.p, angleConnection.atom2.p);

		if (order == 0) {
			return param.keq * Mathf.Pow (Mathf.Deg2Rad * (angle - param.req), 2f);
		} else if (order == 1) {
			return 2f * param.keq * (Mathf.Deg2Rad * (angle - param.req));
		} else if (order == 2) {
			return 2f * param.keq;
		} else {
			return 0f;
		}
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

	public float ETorsion(DihedralConnection dihedralConnection, int order, bool suppress=false) {
		Torsion param = GetTorsionParameter (dihedralConnection);
		if (param == null) {
			if (!suppress)
				Debug.LogError (string.Format ("No parameter for Torsion: {0}-{1}-{2}-{3}", dihedralConnection.atom0.amberName, dihedralConnection.atom1.amberName, dihedralConnection.atom2.amberName, dihedralConnection.atom3.amberName));
			return 0f;
		}

		float dihedral = dihedralConnection.dihedral;

		float e = 0f;
		if (order == 0) {
			if (param.v0 != 0)
				e += 0.5f * param.v0 * (1f + Mathf.Cos (Mathf.Deg2Rad *  (dihedral - param.gamma0)));
			if (param.v1 != 0)
				e += 0.5f * param.v1 * (1f + Mathf.Cos (2f * Mathf.Deg2Rad *  (dihedral - param.gamma1)));
			if (param.v2 != 0)
				e += 0.5f * param.v2 * (1f + Mathf.Cos (3f * Mathf.Deg2Rad *  (dihedral - param.gamma2)));
			if (param.v3 != 0)
				e += 0.5f * param.v3 * (1f + Mathf.Cos (4f * Mathf.Deg2Rad *  (dihedral - param.gamma3)));
		} else if (order == 1) {
			if (param.v0 != 0)
				e -= 0.5f * param.v0 * (Mathf.Sin (Mathf.Deg2Rad *  (dihedral - param.gamma0)));
			if (param.v1 != 0)
				e -= 0.5f * param.v1 * 2f * (Mathf.Sin (2f * Mathf.Deg2Rad *  (dihedral - param.gamma1)));
			if (param.v2 != 0)
				e -= 0.5f * param.v2 * 3f * (Mathf.Sin (3f * Mathf.Deg2Rad *  (dihedral - param.gamma2)));
			if (param.v3 != 0)
				e -= 0.5f * param.v3 * 4f * (Mathf.Sin (4f * Mathf.Deg2Rad *  (dihedral - param.gamma3)));
		} else if (order == 2) {
			if (param.v0 != 0)
				e -= 0.5f * param.v0 * (Mathf.Cos (Mathf.Deg2Rad *  (dihedral - param.gamma0)));
			if (param.v1 != 0)
				e -= 0.5f * param.v1 * 4f * (Mathf.Cos (2f * Mathf.Deg2Rad *  (dihedral - param.gamma1)));
			if (param.v2 != 0)
				e -= 0.5f * param.v2 * 9f * (Mathf.Cos (3f * Mathf.Deg2Rad *  (dihedral - param.gamma2)));
			if (param.v3 != 0)
				e -= 0.5f * param.v3 * 16f * (Mathf.Cos (4f * Mathf.Deg2Rad *  (dihedral - param.gamma3)));
		}
		return e / param.npaths;
	}

	public float EImproperTorsion(DihedralConnection dihedralConnection, int order, bool suppress=false) {
		ImproperTorsion param = GetImproperTorsionParameter (dihedralConnection);
		if (param == null) {
			if (!suppress)
				Debug.LogError (string.Format ("No parameter for Improper Torsion: {0}-{1}-{2}-{3}", dihedralConnection.atom0.amberName, dihedralConnection.atom1.amberName, dihedralConnection.atom2.amberName, dihedralConnection.atom3.amberName));
			return 0f;
		}

		float dihedral = dihedralConnection.dihedral;

		float e = 0f;
		if (order == 0) {
			e += 0.5f * param.v * (1f + Mathf.Cos (Mathf.Deg2Rad *  (param.periodicity * dihedral - param.gamma)));
		} else if (order == 1) {
			e -= 0.5f * param.v * param.periodicity * (Mathf.Sin (Mathf.Deg2Rad *  (param.periodicity * dihedral - param.gamma)));
		} else if (order == 2) {
			e -= 0.5f * param.v * param.periodicity * param.periodicity * (Mathf.Cos (Mathf.Deg2Rad *  (param.periodicity * dihedral - param.gamma)));
		}
		return e;

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

public class PrecomputedNonBonding {
	public float vdwR;
	public float vdwV;

	public float q0q1;
	public float dielectricConstant;
	public int cType;
	public Atom atom0;
	public Atom atom1;

	public PrecomputedNonBonding(VdW vdw0, VdW vdw1, Atom atom0, Atom atom1, float vScale, float cScale, float dielectricConstant) {

		this.atom0 = atom0;
		this.atom1 = atom1;
		this.dielectricConstant = dielectricConstant * cScale;
		this.q0q1 = atom0.partialCharge * atom0.partialCharge;
		this.vdwR = (vdw0.r + vdw1.r) * 0.5f;
		this.vdwV = Mathf.Sqrt (vdw0.v * vdw1.v) * vScale;
	}

	public void GetEnergy(float r, float angstromToBohr, int order, float[] vdwEnergies, float[] coulombEnergies) {
		Mathematics.EVdWAmber(r, vdwV, vdwR, vdwEnergies);
		Mathematics.EElectrostaticR1(r * angstromToBohr, q0q1, dielectricConstant, coulombEnergies);
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
