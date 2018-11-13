using System.Collections;
using System.Collections.Generic;
using UnityEngine;

#if UNITY_EDITOR
using UnityEditor;
#endif

public class Atoms : MonoBehaviour {

	public Camera activeCamera;

	//Atoms
	public List<Atom> atomList;
	public GameObject atomsHolder;
	public List<int> selection;
	public List<int> userSelection;

	//Graph
	public Graph graph;

	//Gaussian
	public GaussianCalculator gaussianCalculator;

	//Settings
	public Settings globalSettings;
	public string sourceFilename;

	//Mesh
	public Mesh combinedMesh;

	//Amber Reps
	public MeshFilter amberRepMeshFilter;
	public Mesh amberRepMesh;

	private Connection stretchConnection;
	private bool showStretchRep;
	public Vector3 stretchRepForward;
	public Vector3 stretchRepUp;



	public PostRender postRenderer;

	public GameObject haloHolder;
	private SelectionHalo hoverHalo;
	public float haloZSyncTime;


	private Dictionary<int, SelectionHalo> selectionDict;
	public Color hoverColor = new Color(1.2f, 0.7f, 0.2f, 0.5f);
	public Color selectionColor = new Color(0.3f, 0.3f, 1.2f, 0.5f);

	private int selectionCount;

	//_busy is a counter, counting the depth of demanding scripts currently running in this object
	private int _busy;
	public bool busy {
		get { return (_busy > 0); }
	}

	private bool _active;
	public bool active {
		get { return _active;}
		set { 
			if (value) {
				if (!_active) {

				}
			} else {
				if (active) {
					ClearSelection();
				}
			}
			_active = value;
		}
	}

	private int _hoveredAtom;
	public int hoveredAtom {
		get {return _hoveredAtom;}
		set {
			_hoveredAtom = value;
			if (value == -1) {
				hoverHalo.ClearAtom();
			} else {
				hoverHalo.SetAtom(this[value]);
			}
		}
	}

	void Awake() {
		_busy++;

		active = false;

		atomList = new List<Atom> ();
		userSelection = new List<int>();

		atomsHolder = new GameObject ("Atoms");
		atomsHolder.transform.parent = transform;

		graph = Instantiate<Graph> (GameObject.FindObjectOfType<PrefabManager>().graphPrefab, transform);

		combinedMesh = GetComponent<MeshFilter> ().mesh;
		amberRepMesh = amberRepMeshFilter.mesh;

		selectionDict = new Dictionary<int, SelectionHalo>();

		_hoveredAtom = -1;

		activeCamera = Camera.main;

		hoverHalo = Instantiate<SelectionHalo>(GameObject.FindObjectOfType<PrefabManager>().selectionHaloPrefab, haloHolder.transform);
		hoverHalo.parent = this;
		hoverHalo.SetColor(hoverColor);

		selectionCount = 0;

		_busy--;
	}

	public void setSourceFile(string filename) {
		sourceFilename = filename;
	}

/* 
	//Render
	public void RenderAll(bool optimise=true) {
		_busy++;
		
		CombineInstance[] combine = new CombineInstance[atomList.Count + graph.connections.Count];

		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			atomList [atomNum].showSphere = true;

			if (optimise) {
				combine [atomNum].mesh = atomList [atomNum].meshFilter.sharedMesh;
				combine [atomNum].transform = atomList [atomNum].transform.localToWorldMatrix;
				Destroy (atomList [atomNum].meshFilter);
				Destroy (atomList [atomNum].mesh);
			}
		}

		for (int connectionNum = 0; connectionNum < graph.connections.Count; connectionNum++) {
			graph.connections [connectionNum].Render ();

			if (optimise) {
				combine [atomList.Count + connectionNum].mesh = graph.connections [connectionNum].meshFilter.sharedMesh;
				combine [atomList.Count + connectionNum].transform = graph.connections [connectionNum].transform.localToWorldMatrix;
				Destroy (graph.connections [connectionNum].meshFilter);
				Destroy (graph.connections [connectionNum].mesh);
			}
		}

		if (optimise) {
			transform.GetComponent<MeshFilter> ().mesh.CombineMeshes (combine);
		}

		_busy--;
	}*/

	public void GetDefaults() {
		_busy++;

		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			atomList [atomNum].GetDefaults ();
		}
	}

	//Vectors
	public Vector3 GetPosition(int a0) {
		return atomList[a0].p;
	}
	public Vector3 GetVector(int a0, int a1) {
		return GetPosition(a1) - GetPosition(a0);
	}

	public float GetDistance(int a0, int a1) {
		return Vector3.Distance (GetPosition(a0), GetPosition(a1));
	}

	public float GetAngleDeg(int a0, int a1, int a2) {
		return Mathematics.GetAngleDeg(GetPosition(a0), GetPosition(a1), GetPosition(a2));
	}

	public float GetDihedralDeg(int a0, int a1, int a2, int a3) {
		return Mathematics.GetDihedralDeg(GetPosition(a0), GetPosition(a1), GetPosition(a2), GetPosition(a3));
	}


	public Vector3 centre {
		get {
			return GetCentre ();
		}
	}

	public double[] masses {
		get {
			int _size = size;
			double[] masses = new double[_size];
			for (int atomNum=0; atomNum < _size; atomNum++) {
				masses[atomNum] = this[atomNum].mass;
			}
			return masses;
		}
	}

	public Vector3 GetCentre() {
		if (size == 0)
			return Vector3.zero;

		Vector3 totalV = Vector3.zero;
		
		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			totalV += atomList [atomNum].transform.position;
		}
		return totalV / size;
	}

	//Atoms
	[SerializeField]
	public int size {
		get {
			return atomList.Count;
		}
	}

	public void AddAtom(Atom atom) {
		atom.resolution = globalSettings.ballResolution;
		atom.index = size;
		atom.transform.SetParent(atomsHolder.transform);
		atom.parent = this;
		atomList.Add (atom);

		//Resize connectivity matrices
		/*
		int[,] connectivityMatrix = new int[size,size];
		int[,] connectivityDistanceMatrix = new int[size,size];
		for (int i = 0; i < size - 1; i++) {
			for (int j = 0; j < size - 1; j++) {
				connectivityMatrix [i, j] = graph.connectivityMatrix [i, j];
				connectivityDistanceMatrix [i, j] = graph.connectivityDistanceMatrix [i, j];
			}
		}

		graph.connectivityMatrix = connectivityMatrix;
		graph.connectivityDistanceMatrix = connectivityDistanceMatrix;
		*/
	}

	public void InsertAtom(Atom atom, int index) {
		atomList.Insert (index, atom);

		//Resize connectivity matrices
		/*
		int[,] connectivityMatrix = new int[size,size];
		int[,] connectivityDistanceMatrix = new int[size,size];
		for (int i = 0; i < size - 1; i++) {
			for (int j = 0; j < size - 1; j++) {
				graph.connectivityMatrix [i, j] = connectivityMatrix [i + (i >= index ? 1 : 0), j + (j >= index ? 1 : 0)];
				graph.connectivityDistanceMatrix [i, j] = connectivityDistanceMatrix [i + (i >= index ? 1 : 0), j + (j >= index ? 1 : 0)];
			}
		}
		*/

	}

	public void DeleteAtom(int index) {
		atomList.RemoveAt (index);

		//Resize connectivity matrices
		/*
		int[,] connectivityMatrix = new int[size,size];
		int[,] connectivityDistanceMatrix = new int[size,size];
		for (int i = 0; i < size + 1; i++) {
			for (int j = 0; j < size + 1; j++) {
				graph.connectivityMatrix [i, j] = connectivityMatrix [i - (i > index ? 1 : 0), j - (j > index ? 1 : 0)];
				graph.connectivityDistanceMatrix [i, j] = connectivityDistanceMatrix [i - (i > index ? 1 : 0), j - (j > index ? 1 : 0)];
			}
		}
		*/

	}

	public void Reindex() {
		for (int index = 0; index < atomList.Count; index++) {
			atomList [index].index = index;
		}
	}

	//Indexing
	public Atom this[int index] {
		get {
			return atomList [index];
		}
		set {
			atomList [index] = value;
		}
	}

	public Atoms this[List<int> indices] {
		get {
			Atoms newAtoms = Instantiate<Atoms>(GameObject.FindObjectOfType<PrefabManager>().atomsPrefab);
			newAtoms.gaussianCalculator = Instantiate<GaussianCalculator> (GameObject.FindObjectOfType<PrefabManager>().gaussianCalculatorPrefab, newAtoms.transform);
			newAtoms.globalSettings = this.globalSettings;


			foreach (int index in indices) {
				newAtoms.AddAtom(this[index]);
			}
			return newAtoms;
		}
	}

	public IEnumerator GetEnumerator() {
		foreach (Atom atom in atomList) {
			yield return atom;
		}
	}

	//Selection
	public void Select(int a0, bool updateAmberReps=true) {
		if (!selection.Contains(a0)) {
			SelectionHalo newHalo = Instantiate<SelectionHalo>(hoverHalo, haloHolder.transform);
			newHalo.SetColor(selectionColor);
			newHalo.SetAtom(this[a0]);
			newHalo.sizeRatio = 2.2f;
			newHalo.parent = this;
			
			selectionDict.Add(a0, newHalo);
			selection.Add(a0);
			selectionCount++;

			if (updateAmberReps)
				UpdateAmberReps();
		}
	}

	public void Deselect(int a0, bool updateAmberReps=true) {
		if (selection.Contains(a0)) { 
			GameObject.Destroy(selectionDict[a0].gameObject);
			selectionDict.Remove(a0);
			selection.Remove(a0);
			selectionCount--;

			if (updateAmberReps)
				UpdateAmberReps();
		}
	}

	public void ToggleSelect(int a0) {
		if (selection.Contains(a0)) {
			Deselect(a0);
		} else {
			Select(a0);
		}
	}

	public void ClearSelection() {
		List<int> tSelection = new List<int>(selection);
		foreach(int i in tSelection) {
			Deselect(i, false);
		}
		UpdateAmberReps();
	}

	public void UpdateAmberReps() {
		showStretchRep = false;
		if (selectionCount == 2) {
			DrawStretchRep(selection[0], selection[1]);
		} else {
			amberRepMesh.Clear();
		}
	}

	void Update() {
		if (active) {
			haloZSyncTime = Time.deltaTime + haloZSyncTime % 1f;

			if (showStretchRep) {
				amberRepMeshFilter.transform.LookAt(stretchConnection.atom1.transform);
				stretchRepForward = amberRepMeshFilter.transform.forward;
				stretchRepUp = Vector3.Cross(activeCamera.transform.rotation.eulerAngles, stretchRepForward);
				//amberRepMeshFilter.transform.localRotation = Quaternion.LookRotation(stretchRepForward, stretchRepUp);
			}
		}
	}

	public void SetCentre(Vector3 newCentre) {
		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			atomList [atomNum].transform.position -= newCentre;
		}
	}

	void DrawStretchRep(int a0, int a1) {

		List<float> lengths = new List<float>();
		List<float> energies = new List<float>();

		stretchConnection = graph.GetConnection(a0, a1);
		if (stretchConnection == null) {
			return;
		}

		Stretch stretch = graph.GetStretchParameter(stretchConnection);
		if (stretch == null) {
			return;
		}

		showStretchRep = true;
		
		float length = stretch.req - (globalSettings.amberStretchSteps + 1) * globalSettings.amberStretchInterval;

		for (int segmentNum = 0; segmentNum < (globalSettings.amberStretchSteps * 2) + 1; segmentNum++) {
			lengths.Add(length);
			energies.Add(Mathematics.EStretch(length, stretch.keq, stretch.req, 0));

			length += globalSettings.amberStretchInterval;
		}

		StretchRep.GenerateStretchRep(globalSettings.amberRepResolution, globalSettings.amberRepThickness, globalSettings.amberStretchRepOffset, lengths, energies, amberRepMesh);

		amberRepMeshFilter.transform.localPosition = stretchConnection.atom0.p;
	}

}
