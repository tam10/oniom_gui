using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class IcosphereGenerator {

	private static Dictionary<int, SphereGeometry> cachedGeometries = new Dictionary<int, SphereGeometry>();

	public static void GenerateSphere(int resolution, Mesh mesh, Color color) {
		if (cachedGeometries.ContainsKey (resolution)) {
			cachedGeometries [resolution].SetMesh (mesh, color);
		} else {
			
			//Generate new sphere
			Dictionary<long, int> cachedVertIndices = new Dictionary<long, int> ();
			List<Vector3> vertices = GetIcosahedronVertices();
			List<Face> faces = GetIcosahedronFaces ();

			for (int subdivision = 0; subdivision < resolution; subdivision++) {
				int faceCount = faces.Count;
				for (int faceIndex=0; faceIndex < faceCount; faceIndex++) {
					RefineFace (faceIndex, faces, vertices, cachedVertIndices);

				}
			}

			int[] triangles = FacesToTriangles (faces);

			Vector3[] verts = vertices.ToArray ();

			SphereGeometry sphere = new SphereGeometry (triangles, verts, verts);
			sphere.SetMesh (mesh, color);

			cachedGeometries [resolution] = sphere;
		}
	}

	private static List<Vector3> GetIcosahedronVertices() {

		float a = 1f;
		float b = (1f + Mathf.Sqrt (5f)) / 2f;
		float c = 0f;

		List<Vector3> vertices = new List<Vector3> ();

		vertices.Add(new Vector3(-a, b, c).normalized);
		vertices.Add(new Vector3(a, b, c).normalized);
		vertices.Add(new Vector3(-a, -b, c).normalized);
		vertices.Add(new Vector3(a, -b, c).normalized);
						 
		vertices.Add(new Vector3(c, -a, b).normalized);
		vertices.Add(new Vector3(c, a, b).normalized);
		vertices.Add(new Vector3(c, -a, -b).normalized);
		vertices.Add(new Vector3(c, a, -b).normalized);
						 
		vertices.Add(new Vector3(b, c, -a).normalized);
		vertices.Add(new Vector3(b, c, a).normalized);
		vertices.Add(new Vector3(-b, c, -a).normalized);
		vertices.Add(new Vector3(-b, c, a).normalized);

		return vertices;
	}

	private static List<Face> GetIcosahedronFaces() {

		List<Face> faces = new List<Face>();

		faces.Add(new Face(0, 1, 7));
		faces.Add(new Face(7, 1, 8));
		faces.Add(new Face(0, 7, 10));
		faces.Add(new Face(10, 7, 6));
		faces.Add(new Face(0, 10, 11));
		faces.Add(new Face(11, 10, 2));
		faces.Add(new Face(0, 11, 5));
		faces.Add(new Face(5, 11, 4));
		faces.Add(new Face(0, 5, 1));
		faces.Add(new Face(1, 5, 9));

		faces.Add(new Face(3, 6, 8));
		faces.Add(new Face(8, 6, 7));
		faces.Add(new Face(3, 2, 6));
		faces.Add(new Face(6, 2, 10));
		faces.Add(new Face(3, 4, 2));
		faces.Add(new Face(2, 4, 11));
		faces.Add(new Face(3, 9, 4));
		faces.Add(new Face(4, 9, 5)); 
		faces.Add(new Face(3, 8, 9));
		faces.Add(new Face(9, 8, 1));

		return faces;
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

	private static void RefineFace(int faceIndex, List<Face> faces, List<Vector3> verts, Dictionary<long, int> cachedVertIndices) {

		int i0;
		int i1;
		int vertIndex0;
		int vertIndex1;
		Vector3 vert0;
		Vector3 vert1;
		Vector3 midVert;
		int index;
		Face outerFace = faces [faceIndex];
		Face innerFace = new Face();

		//Add inner vertices if needed
		for (i0 = 0; i0 < 3; i0++) {
			
			i1 = i0 == 2 ? 0 : i0 + 1;

			vertIndex0 = outerFace [i0] > outerFace [i1] ? outerFace [i1] : outerFace [i0];
			vertIndex1 = outerFace [i0] > outerFace [i1] ? outerFace [i0] : outerFace [i1];
			
			long key = ((long)vertIndex0 << 32) + vertIndex1;

			//Get vertex if it's in the cache, or calculate and add it to cache
			int cachedInt;
			if (cachedVertIndices.TryGetValue (key, out cachedInt)) {
				index = cachedInt;
			} else {
				vert0 = verts [vertIndex0];
				vert1 = verts [vertIndex1];

				midVert = (vert0 + vert1).normalized;
				index = verts.Count;
				verts.Add (midVert);
				cachedVertIndices.Add (key, index);
			}
			innerFace [i0] = index;
		}

		//Change the original face to inner face and add faces
		faces.Add(new Face(innerFace[0], outerFace[1], innerFace[1]));
		faces.Add(new Face(innerFace[1], outerFace[2], innerFace[2]));
		faces.Add(new Face(innerFace[2], outerFace[0], innerFace[0])); 
		faces[faceIndex] = new Face(innerFace);
	}
}

class SphereGeometry {
	public int[] triangles;
	public Vector3[] vertices;
	public Vector3[] normals;

	public SphereGeometry(int[] triangles, Vector3[] vertices, Vector3[] normals) {
		this.triangles = triangles;
		this.vertices = vertices;
		this.normals = normals;
	}

	public void SetMesh(Mesh mesh, Color color) {
		mesh.vertices = vertices;
		mesh.normals = normals;
		mesh.triangles = triangles;

		Color32[] colors = new Color32[vertices.Length];
		for (int c = 0; c < colors.Length; c++)
			colors [c] = color;
		mesh.colors32 = colors;

	}
}
