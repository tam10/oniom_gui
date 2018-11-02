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

	//Graph
	public Graph graph;
	public Graph graphPrefab;

	//Gaussian
	public GaussianCalculator gaussianCalculator;

	//Parameters
	public Parameters parameters;

	//Settings
	public Settings globalSettings;
	public string sourceFilename;

	//Mesh
	public Mesh combinedMesh;

	public PostRender postRenderer;

	//_busy is a counter, counting the depth of demanding scripts currently running in this object
	private int _busy;
	public bool busy {
		get { return (_busy > 0); }
	}

	public int hoveredAtom;

	void Awake() {
		_busy++;

		atomList = new List<Atom> ();
		atomsHolder = new GameObject ("Atoms");
		atomsHolder.transform.parent = transform;

		graph = Instantiate<Graph> (graphPrefab, transform);

		combinedMesh = GetComponent<MeshFilter> ().mesh;

		hoveredAtom = -1;

		_busy--;
	}

	public void setSourceFile(string filename) {
		sourceFilename = filename;
	}

	//Render
	public void RenderAll(bool optimise=true) {
		_busy++;
		
		CombineInstance[] combine = new CombineInstance[atomList.Count + graph.connections.Count];

		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			atomList [atomNum].Render ();

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
	}

	public void GetDefaults() {
		_busy++;

		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			atomList [atomNum].GetDefaults ();
		}
	}

	//Vectors
	public Vector3 getVector(int a0, int a1) {
		return atomList [a1].transform.localPosition - atomList [a0].transform.localPosition;
	}

	public float getDistance(int a0, int a1) {
		return Vector3.Distance (atomList [a0].transform.localPosition, atomList [a1].transform.localPosition);
	}

	public float getAngle(int a0, int a1, int a2) {
		return Vector3.Angle (getVector (a1, a0), getVector (a1, a2));
	}

	public float getDihedral(int a0, int a1, int a2, int a3) {
		Vector3 v10 = getVector (a1, a0);
		Vector3 v12 = getVector (a2, a1);
		Vector3 v23 = getVector (a3, a2);

		v10.Normalize ();

		Vector3 w10 = v10 - v12 * Vector3.Dot (v10, v12);
		Vector3 w23 = v23 - v12 * Vector3.Dot (v23, v12);

		float x = Vector3.Dot (w10, w23);
		float y = Vector3.Dot (Vector3.Cross (v12, w10), w23);

		return Mathf.Atan2 (y, x);
	}


	public Vector3 centre {
		get {
			return GetCentre ();
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

	//Connectivity
	public void Connect(int a0, int a1, int bondOrder=1) {
		if (bondOrder == 0) {

		}
	}

	public void Disconnect(int a0, int a1) {
		
	}


	//Selection
	public void Select(int a0) {

	}

	public void SetCentre(Vector3 newCentre) {
		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			atomList [atomNum].transform.position -= newCentre;
		}
	}


}
