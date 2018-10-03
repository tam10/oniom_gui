using System.Collections;
using System.Collections.Generic;
using UnityEngine;

#if UNITY_EDITOR
using UnityEditor;
#endif

public class Atoms : MonoBehaviour {


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
	private string sourceFiletype;

	//Mesh
	public Mesh combinedMesh;

	void Awake() {
		atomsHolder = new GameObject ("Atoms");
		atomsHolder.transform.parent = transform;

		graph = Instantiate<Graph> (graphPrefab, transform);
		graph.SetAtoms(this);

		combinedMesh = GetComponent<MeshFilter> ().mesh;
	}

	public void setSourceFile(string filename, string filetype) {
		sourceFilename = filename;
		sourceFiletype = filetype;
	}

	//Render
	public void RenderAll() {
		
		CombineInstance[] combine = new CombineInstance[atomList.Count + graph.connections.Count];

		for (int atomNum = 0; atomNum < atomList.Count; atomNum++) {
			atomList [atomNum].Render ();

			combine [atomNum].mesh = atomList [atomNum].meshFilter.sharedMesh;
			combine [atomNum].transform = atomList [atomNum].transform.localToWorldMatrix;
			Destroy(atomList [atomNum].meshFilter);
			Destroy(atomList [atomNum].mesh);
		}

		for (int connectionNum = 0; connectionNum < graph.connections.Count; connectionNum++) {
			graph.connections [connectionNum].Render ();

			combine [atomList.Count + connectionNum].mesh = graph.connections [connectionNum].meshFilter.sharedMesh;
			combine [atomList.Count + connectionNum].transform = graph.connections [connectionNum].transform.localToWorldMatrix;
			Destroy(graph.connections [connectionNum].meshFilter);
			Destroy(graph.connections [connectionNum].mesh);
		}

		transform.GetComponent<MeshFilter>().mesh.CombineMeshes(combine);


	}

	//Vectors
	public Vector3 getVector(int a0, int a1) {
		return atomList [a1].transform.position - atomList [a0].transform.position;
	}

	public float getDistance(int a0, int a1) {
		return Vector3.Distance (atomList [a0].transform.position, atomList [a1].transform.position);
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

	//Atoms
	[SerializeField]
	public int size {
		get {
			return atomList.Count;
		}
	}

	public void AddAtom(Atom atom) {
		atomList.Add (atom);

		//Resize connectivity matrices
		int[,] connectivityMatrix = new int[size,size];
		int[,] connectivityDistanceMatrix = new int[size,size];
		for (int i = 0; i < size - 1; i++) {
			for (int j = 0; i < size - 1; i++) {
				graph.connectivityMatrix [i, j] = connectivityMatrix [i, j];
				graph.connectivityDistanceMatrix [i, j] = connectivityDistanceMatrix [i, j];
			}
		}

	}

	public void InsertAtom(Atom atom, int index) {
		atomList.Insert (index, atom);

		//Resize connectivity matrices
		int[,] connectivityMatrix = new int[size,size];
		int[,] connectivityDistanceMatrix = new int[size,size];
		for (int i = 0; i < size - 1; i++) {
			for (int j = 0; i < size - 1; i++) {
				graph.connectivityMatrix [i, j] = connectivityMatrix [i + (i >= index ? 1 : 0), j + (j >= index ? 1 : 0)];
				graph.connectivityDistanceMatrix [i, j] = connectivityDistanceMatrix [i + (i >= index ? 1 : 0), j + (j >= index ? 1 : 0)];
			}
		}

	}

	public void DeleteAtom(int index) {
		atomList.RemoveAt (index);

		//Resize connectivity matrices
		int[,] connectivityMatrix = new int[size,size];
		int[,] connectivityDistanceMatrix = new int[size,size];
		for (int i = 0; i < size + 1; i++) {
			for (int j = 0; i < size + 1; i++) {
				graph.connectivityMatrix [i, j] = connectivityMatrix [i - (i > index ? 1 : 0), j - (j > index ? 1 : 0)];
				graph.connectivityDistanceMatrix [i, j] = connectivityDistanceMatrix [i - (i > index ? 1 : 0), j - (j > index ? 1 : 0)];
			}
		}

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


}
