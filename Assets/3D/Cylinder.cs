using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class CylinderGenerator {

	private static Dictionary<int, CylinderGeometry> cachedGeometries = new Dictionary<int, CylinderGeometry>();

	public static void GenerateCylinder(int resolution, Mesh mesh, Color32 color0, Color32 color1) {
		if (cachedGeometries.ContainsKey (resolution)) {
			cachedGeometries [resolution].SetMesh (mesh, color0, color1);
		} else {

			List<Vector3> vertices = new List<Vector3> ();
			List<int> triangleList = new List<int> ();

			int divisions = 4 * (int)Mathf.Pow (2, resolution);
			for (int division = 0; division < divisions; division++) {
				
				float angle = 2 * Mathf.PI * division / divisions;

				float x = Mathf.Cos (angle);
				float y = Mathf.Sin (angle);

				//Add vertices
				vertices.Add(new Vector3(x, y, -0.5f));
				vertices.Add(new Vector3(x, y, 0.5f));

				int v0 = division * 2;
				int v1 = division * 2 + 1;
				int v2 = (division == divisions - 1 ? 0 : division * 2 + 2);
				int v3 = (division == divisions - 1 ? 1 : division * 2 + 3);

				triangleList.Add (v0);
				triangleList.Add (v2);
				triangleList.Add (v1);
				triangleList.Add (v1);
				triangleList.Add (v2);
				triangleList.Add (v3);

			}

			Vector3[] verts = vertices.ToArray ();

			CylinderGeometry cylinder = new CylinderGeometry (triangleList.ToArray(), verts, verts);
			cylinder.SetMesh (mesh, color0, color1);

			cachedGeometries [resolution] = cylinder;
		}
	}

	private static int[] FacesToTriangles(List<Face> faces) {
		int[] triangles = new int[faces.Count * 3];
		int j;
		for (int i=0;i<faces.Count;i++) {
			for (j = 0; j < 3; j++) {
				triangles [i * 3 + j] = faces [i] [j];
			}
		}
		return triangles;
	}
}


class CylinderGeometry {
	public int[] triangles;
	public Vector3[] vertices;
	public Vector3[] normals;

	public CylinderGeometry(int[] triangles, Vector3[] vertices, Vector3[] normals) {
		this.triangles = triangles;
		this.vertices = vertices;
		this.normals = normals;
	}

	public void SetColors(Mesh mesh, Color color0, Color color1) {

		Color32[] colors = new Color32[vertices.Length];
		for (int i = 0; i < vertices.Length / 2; i++) {
			colors [2 * i] = color0;
			colors [2 * i + 1] = color1;
		}

		mesh.colors32 = colors;
	}

	public void SetMesh(Mesh mesh, Color color0, Color color1) {
		mesh.vertices = vertices;
		mesh.normals = normals;
		mesh.triangles = triangles;

		SetColors(mesh, color0, color1);
	}
}