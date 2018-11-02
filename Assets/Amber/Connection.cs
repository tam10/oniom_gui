using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Connection : MonoBehaviour {

	public Mesh mesh;	
	public MeshFilter meshFilter;
	public Transform cylinderTransform;
	public Atoms parent;

	public Atom atom0;
	public Atom atom1;
	public double bondOrder;

	public int resolution;


	//Geometry
	public float radius {
		get {return _localScale.x;}
		set {
			_localScale.x = value;
			_localScale.y = value;
			UpdateLocalScale();
		}
	}

	public float length {
		get {return _localScale.z;}
		set {
			_localScale.z = value;
			UpdateLocalScale();
		}
	}

	private Vector3 _localScale;
	public Vector3 localScale {
		get {return _localScale;}
		set {
			_localScale = value;
			UpdateLocalScale();
		}
	}

	private void UpdateLocalScale() {
		transform.localScale = _localScale;
	}

	private void UpdateLocalPosition() {
		transform.localPosition = mid;
	}

	private void UpdateLocalRotation() {
		transform.localRotation = Quaternion.LookRotation (end - start);
	}

	public void UpdateLocalTransform() {
		_start = atom0.p;
		_end = atom1.p;

		length = Vector3.Distance(start, end);

		UpdateLocalPosition();
		UpdateLocalScale();
		UpdateLocalRotation();
	}

	public float distance {
		get {return Vector3.Distance(atom1.p, atom0.p);}
	}

	public Vector3 vector {
		get {return atom1.p - atom0.p;}
	}

	private Vector3 _start;
	public Vector3 start {
		get {return _start;}
		set {
			_start = value;
			length = Vector3.Distance(start, end);
		}
	}

	public Vector3 mid {
		get {return (start + end) * 0.5f;}
	}

	private Vector3 _end;
	public Vector3 end {
		get {return _end;}
		set {
			_end = value;
			length = Vector3.Distance(start, end);
		}
	}

	public Vector3 PointAlong(float ratio) {
		return start * (1f - ratio) + end * ratio;
	}


	public Color color0;
	public Color color1;

	public void UpdateColors() {
		color0 = atom0.color;
		color1 = atom1.color;
		UpdateCylinderColors();
	}

	public void UpdateCylinderColors() {
		if (_cylinderRendered) {
			Color32[] colors = new Color32[mesh.colors.Length];
			for (int i = 0; i < mesh.colors.Length / 2; i++) {
				colors [2 * i] = color0;
				colors [2 * i + 1] = color1;
			}

			mesh.colors32 = colors;
		}
	}

	public void UpdateRenderType() {
		showCylinder = (atom0.showSphere && atom1.showSphere);
	}

	private bool _cylinderRendered;
	private bool _showCylinder;
	public bool showCylinder {
		get {return _showCylinder;}
		set {
			if (value) {
				if (!_cylinderRendered) {
					Render();
				} 
				cylinderTransform.gameObject.SetActive(true);
			} else {
				cylinderTransform.gameObject.SetActive(true);
			}
			_showCylinder = value;
		}
	}

	public void Initialise (Atoms atoms, int a0, int a1, double bondOrder, int resolution=1, float radius=0f) {
		this.parent = atoms;
		this.atom0 = atoms [a0];
		this.atom1 = atoms [a1];
		this.bondOrder = bondOrder;
		this.resolution = resolution;

		this.mesh = this.meshFilter.mesh;

		UpdateRenderType();
		UpdateColors();
		this.radius = radius == 0f ? ((this.atom0.vdwRadius * this.atom0.vdwToSphereRatio) + (this.atom1.vdwRadius * this.atom1.vdwToSphereRatio)) / 8f : radius;

	}

	public Stretch GetStretchParameter() {
		foreach (Stretch param in parent.graph.parameters.stretches) {
			if (this.atom0.amberName == param.t0 && this.atom1.amberName == param.t1)
				return param;
			if (this.atom0.amberName == param.t1 && this.atom1.amberName == param.t0)
				return param;
		}
		return null;
	}

	public float EStretch(int order, bool suppress=false) {

		float e = 0f;

		//Not connected
		if (bondOrder == 0.0)
			e = 0f;

		Stretch param = GetStretchParameter ();
		if (param == null) {
			if (!suppress) {
				throw new NoParameterException(typeof(Stretch), atom0, atom1);
			}
			e = 0f;
		} else {
			e = Mathematics.EStretch(distance, param.keq, param.req, order);
		}

		return e;
	}

	public void Render () {

		CylinderGenerator.GenerateCylinder (resolution, mesh, atom0.color, atom1.color);
		_cylinderRendered = true;
		UpdateLocalTransform();

	}

}
