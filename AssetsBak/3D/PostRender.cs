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
		DrawGLConnections ();
	}

	void OnDrawGizmos() {
		DrawGLConnections();
	}

	void DrawGLConnections() {
		if (atomsList == null)
			return;

		lineMaterial.SetPass (0);
		GL.Begin (GL.LINES);

		cameraVector = activeCamera.transform.forward;
		
		for (atomsNum=0;atomsNum<atomsList.Count;atomsNum++) {
			currentAtoms = atomsList [atomsNum];

			for (connectionNum=0;connectionNum<currentAtoms.graph.connections.Count;connectionNum++) {
				connection = currentAtoms.graph.connections [connectionNum];

				if (connection.bondOrder == 1.0) {
					DrawSingle (connection);
				} else {
					DrawDouble (connection, cameraVector);
				}
			}

			//Draw a magenta line from cursor to hovered atom
			//if (atoms.hoveredAtom > -1) {
			//	GL.Color (Color.magenta);
			//	GL.Vertex (inputHandler.worldMousePosition);
			//	GL.Vertex (atoms [atoms.hoveredAtom].transform.position);
			//}
		}
		GL.End ();
	}

	void DrawSingle(Connection connection) {
		GL.Color (connection.color0);
		GL.Vertex (connection.start);
		GL.Vertex (connection.mid);

		GL.Color (connection.color1);
		GL.Vertex (connection.mid);
		GL.Vertex (connection.end);

	}

	void DrawDouble(Connection connection, Vector3 cameraVector) {

		connectionVector = connection.end - connection.start;
		offset = Vector3.Cross (connectionVector, connectionVector).normalized * 0.1f;

		GL.Color (connection.color0);
		GL.Vertex (connection.start);
		GL.Vertex (connection.mid);
		GL.Vertex (connection.PointAlong(0.1f) + offset);
		GL.Vertex (connection.mid + offset);

		GL.Color (connection.color1);
		GL.Vertex (connection.start);
		GL.Vertex (connection.end);
		GL.Vertex (connection.mid + offset);
		GL.Vertex (connection.PointAlong(0.9f) + offset);
	}
}
