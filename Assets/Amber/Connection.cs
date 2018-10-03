using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Connection : MonoBehaviour {

	public Atom atom0;
	public Atom atom1;
	public int a0;
	public int a1;
	public int bondOrder;
	public Atoms atoms;

	public int resolution;
	public MeshFilter meshFilter;
	public Mesh mesh;

	public float radius;

	public void Initialise (Atoms atoms, int a0, int a1, int bondOrder, int resolution=1, float radius=0f) {
		this.atoms = atoms;
		this.atom0 = atoms [a0];
		this.atom1 = atoms [a1];
		this.a0 = a0;
		this.a1 = a1;
		this.bondOrder = bondOrder;
		this.resolution = resolution;

		this.meshFilter = GetComponent<MeshFilter> ();
		this.mesh = this.meshFilter.mesh;

		this.radius = radius == 0f ? ((this.atom0.vdwRadius * this.atom0.vdwToSphereRatio) + (this.atom1.vdwRadius * this.atom1.vdwToSphereRatio)) / 8f : radius;

	}

	public VdW GetVdWParameter(int index) {
		foreach (VdW param in atoms.parameters.vdws) {
			if ((index == 0 && this.atom0.amberName == param.t) || (index == 1 && this.atom1.amberName == param.t))
				return param;
		}
		return null;
	}

	public Stretch GetStretchParameter() {
		foreach (Stretch param in atoms.parameters.stretches) {
			if (this.atom0.amberName == param.t0 && this.atom1.amberName == param.t1)
				return param;
			if (this.atom0.amberName == param.t1 && this.atom1.amberName == param.t0)
				return param;
		}
		return null;
	}

	public float EStretch(int order, bool suppress=false) {
		//Not connected
		if (bondOrder == 0 || bondOrder > 5)
			return 0f;

		Stretch param = GetStretchParameter ();
		if (param == null) {
			if (!suppress)
				Debug.LogError (string.Format ("No parameter for Stretch: {0}-{1}", this.atom0.amberName, this.atom1.amberName));
			return 0f;
		}

		float length = this.atoms.getDistance (this.a0, this.a1);

		if (order == 0) {
			return param.keq * Mathf.Pow (length - param.req, 2f);
		} else if (order == 1) {
			return 2f * param.keq * (length - param.req);
		} else if (order == 2) {
			return 2f * param.keq;
		} else {
			return 0f;
		}
	}

	public float EVdW(int order, bool suppress) {
		VdW vdw0 = GetVdWParameter (0);
		VdW vdw1 = GetVdWParameter (1);

		float r = atoms.getDistance (this.a0, this.a1);
		float r0 = (vdw0.r + vdw1.r) * 0.5f;
		float v = Mathf.Sqrt (vdw0.v * vdw1.v);

		//Amber type VdW
		if (atoms.parameters.nonbonding.vType == 3) {
			if (order == 0) {
				return v * (Mathf.Pow (r0 / r, 12f) - 2f * Mathf.Pow (r0 / r, 6f));
			} else if (order == 1) {
				return 12f * v * (Mathf.Pow (r0 / r, 7f) * -Mathf.Pow (r0 / r, 13f)) / r0;
			} else if (order == 2) {
				return 12f * v * (13f * Mathf.Pow (r0 / r, 14f) - 7f * Mathf.Pow (r0 / r, 8f)) / (r0 * r0);
			}
		} else {
			Debug.LogError (string.Format ("VdW Nonbonding type {0} not available", atoms.parameters.nonbonding.vType));
		}

		return 0f;
	}

	public float EElectrostatic(int order, bool suppress) {
		float eps = atoms.parameters.dielectricConstant * Mathf.PI * 4f;

		float r = atoms.getDistance (this.a0, this.a1);

		//No coulomb term
		if (atoms.parameters.nonbonding.cType == 0) {
			return 0f;
			// Amber style 1/R
		} else if (atoms.parameters.nonbonding.cType == 1) {
			if (order == 0) {
				return (atom0.partialCharge * atom1.partialCharge) / (eps * r);
			} else if (order == 1) {
				return -(atom0.partialCharge * atom1.partialCharge) / (eps * r * r);
			} else if (order == 2) {
				return 2f * (atom0.partialCharge * atom1.partialCharge) / (eps * r * r * r);
			}
			// 1/R^2
		} else if (atoms.parameters.nonbonding.cType == 2) {
			if (order == 0) {
				return (atom0.partialCharge * atom1.partialCharge) / (eps * r * r);
			} else if (order == 1) {
				return - 2f * (atom0.partialCharge * atom1.partialCharge) / (eps * r * r * r);
			} else if (order == 2) {
				return 6f * (atom0.partialCharge * atom1.partialCharge) / (eps * r * r * r * r);
			}
		} else {
			Debug.LogError (string.Format ("VdW Coulomb type {0} not available", atoms.parameters.nonbonding.cType));
		}
		return 0f;
	}

	public void Render () {

		CylinderGenerator.GenerateCylinder (resolution, mesh, atom0.color, atom1.color);
		float l = atoms.getDistance (a0, a1);

		transform.localPosition = (atom0.transform.position + atom1.transform.position) / 2f;
		transform.localScale = new Vector3 (radius, radius, l);
		transform.localRotation = Quaternion.Euler (atoms.getVector (a0, a1));

	}

}
