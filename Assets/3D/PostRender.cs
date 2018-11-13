using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PostRender : MonoBehaviour {

	public List<Atoms> atomsList;
	public Material lineMaterial;
	public InputHandler inputHandler;

	public Camera activeCamera;

	private int atomsNum;
	private Atoms currentAtoms;
	private int connectionNum;
	private Connection connection;

	//Fog
	private float fogRatio;
	private float fogStartDistance;
	private float fogEndDistance;
	private float distance;
	private float fogAmount;

	Vector3 cameraVector;
	Vector3 connectionVector;
	Vector3 offset;

	// Use this for initialization
	void Start () {
		atomsList = new List<Atoms> ();
	}

	public void AddAtoms(Atoms atoms) {
		if (!atomsList.Contains(atoms))
			atomsList.Add (atoms);
	}

	public void RemoveAtoms(Atoms atoms) {
		if (atomsList.Contains(atoms))
			atomsList.Remove (atoms);
	}

	void OnPostRender() {
		try {
			DrawGLConnections();
		} catch {
			return;
		}
	}

	void OnDrawGizmos() {
 
		try {
			DrawGLConnections();
		} catch {
			return;
		}

	}

	void DrawGLConnections() {
		if (atomsList == null)
			return;

		lineMaterial.SetPass (0);
		GL.Begin (GL.LINES);

		cameraVector = activeCamera.transform.forward;
		
		for (atomsNum=0;atomsNum<atomsList.Count;atomsNum++) {
			currentAtoms = atomsList [atomsNum];

			fogStartDistance = currentAtoms.globalSettings.fogStartDistance;
			fogEndDistance = currentAtoms.globalSettings.fogEndDistance;
			fogRatio = currentAtoms.globalSettings.fogRatio;

			foreach (Connection connection in currentAtoms.graph.connectionsDict.Values) {

				distance = Vector3.Distance(activeCamera.transform.position, connection.mid);
				fogAmount = Mathf.Lerp(0f, fogRatio, (distance - fogStartDistance) / fogEndDistance);


				if (connection.bondOrder == 1.0) {
					DrawSingle (connection);
				} else {
					DrawDouble (connection, cameraVector);
				}
			}
		}
		GL.End ();
	}

	void DrawSingle(Connection connection) {
		GL.Color (Color.Lerp(connection.color0, Color.black, fogAmount));
		GL.Vertex (connection.start);
		GL.Vertex (connection.mid);

		GL.Color (Color.Lerp(connection.color1, Color.black, fogAmount));
		GL.Vertex (connection.mid);
		GL.Vertex (connection.end);

	}

	void DrawDouble(Connection connection, Vector3 cameraVector) {

		connectionVector = connection.end - connection.start;
		offset = Vector3.Cross (connectionVector, cameraVector).normalized * 0.05f;

		GL.Color (Color.Lerp(connection.color0, Color.black, fogAmount));
		GL.Vertex (connection.start);
		GL.Vertex (connection.mid);
		GL.Vertex (connection.PointAlong(0.1f) + offset);
		GL.Vertex (connection.mid + offset);

		GL.Color (Color.Lerp(connection.color1, Color.black, fogAmount));
		GL.Vertex (connection.mid);
		GL.Vertex (connection.end);
		GL.Vertex (connection.mid + offset);
		GL.Vertex (connection.PointAlong(0.9f) + offset);
	}
}
