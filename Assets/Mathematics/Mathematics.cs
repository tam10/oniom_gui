using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class Mathematics {

	public static float GetAngleDeg(Vector3 p0, Vector3 p1, Vector3 p2) {
		return Vector3.Angle(p1 - p0, p1 - p2);
	}

	public static float GetDihedralDeg(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3) {

		Vector3 v10 = (p1 - p0).normalized;
		Vector3 v12 = (p1 - p2).normalized;
		Vector3 v23 = (p2 - p3).normalized;

		Vector3 w1 = v10 - v12 * Vector3.Dot (v10, v12);
		Vector3 w2 = v23 - v12 * Vector3.Dot (v23, v12);

		return Mathf.Atan2(Vector3.Dot(v12, Vector3.Cross(w1, w2)), Vector3.Dot(w1, w2));
	}

	public static float EStretch(float length, float keq, float req, int order) {

		float e;

		if (order == 0) {
			e = keq * Mathf.Pow (length - req, 2f);
		} else if (order == 1) {
			e = 2f * keq * (length - req);
		} else if (order == 2) {
			e = 2f * keq;
		} else {
			e = 0f;
		}
		
		return e;
	}

	public static void EStretch(float r, float keq, float req, float[] energies) {

		float _r = r - req;
		energies[2] = 2f * keq;
		energies[1] = energies[2] * _r;
		energies[0] = energies[1] * _r * 0.5f;

	}
	public static float EVdWAmber(float r, float v, float req, int order){

		float e;
		if (v == 0f) {
			e = 0f;
		} else if (order == 0) {
			e = v * (Mathf.Pow (req / r, 12f) - 2f * Mathf.Pow (req / r, 6f));
		} else if (order == 1) {
			e = 12f * v * (Mathf.Pow (req / r, 7f) - Mathf.Pow (req / r, 13f)) / r;
		} else if (order == 2) {
			e = 12f * v * (13f * Mathf.Pow (req / r, 14f) - 7f * Mathf.Pow (req / r, 8f)) / (r * r);
		} else {
			e = 0f;
		}

		return e;
	}

	public static void EVdWAmber(float r, float v, float req, float[] energies){

		float _r = v * req / r;
		float _r2 = _r * _r;
		float _r6 = _r2 * _r2 * _r2;
		float _r12 = _r6 * _r6;

		energies[0] = _r12 - 2 * _r6;
		energies[1] = 12f * (_r6 * _r - _r12 * _r) / r;
		energies[2] = 12f * (13f * _r12 * _r - 7f * _r6 * _r2) / (r * r);
	}

	public static float EElectrostatic(float r, float q0q1, float eps, int cType, int order) {
		
		float e = 0f;

		//No coulomb term
		if (cType == 0 || q0q1 == 0f) {
			e = 0f;
			// Amber style 1/R
		} else if (cType == 1) {
			if (order == 0) {
				e = q0q1 / (eps * r);
			} else if (order == 1) {
				e = -q0q1 / (eps * r * r);
			} else if (order == 2) {
				e = 2f * q0q1 / (eps * r * r * r);
			}
			// 1/R^2
		} else if (cType == 2) {
			if (order == 0) {
				e = q0q1 / (eps * r * r);
			} else if (order == 1) {
				e = - 2f * q0q1 / (eps * r * r * r);
			} else if (order == 2) {
				e = 6f * q0q1 / (eps * r * r * r * r);
			}
		} else {
			Debug.LogError (string.Format ("VdW Coulomb type {0} not available", cType));
		}
		return e;
	}

	public static void EElectrostaticR1(float r, float q0q1, float eps, float[] energies) {
		
		if (q0q1 == 0) {return;}
		float _r = 1f / r;
		energies[0] = q0q1 * _r / eps;
		energies[1] = - energies[0] * _r;
		energies[2] = - 2f * energies[1] * _r;
	}

	public static void EElectrostaticR2(float r, float q0q1, float eps, float[] energies) {
		
		if (q0q1 == 0) {return;}
		float _r = 1f / r;
		energies[0] = q0q1 * _r * _r / eps;
		energies[1] = -2f * energies[0] * _r;
		energies[2] = - 3f * energies[1] * _r;
	}

	public static void EAngle(float da, float keq, float[] energies) {
		if (keq == 0f) {return;}
		energies[2] = 2f * keq;
		energies[1] = energies[2] * da;
		energies[0] = energies[1] * da * 0.5f;
	}

	public static void EImproperTorsion(float dihedral, float v, float gamma, float periodicity, float[] energies) {
		if (v == 0f) {return;}
		float t = Mathf.Deg2Rad * (periodicity * dihedral - gamma);
		v *= 0.5f;
		energies[0] = v * (1f + Mathf.Cos(t));
		energies[1] = - v * periodicity * (Mathf.Sin (t));
		energies[2] = - v * periodicity * periodicity * (Mathf.Cos (t));
	}

	public static void EImproperTorsionAdditive(float dihedral, float v, float gamma, float periodicity, float[] energies) {
		if (v == 0f) {return;}
		float t = Mathf.Deg2Rad * (periodicity * dihedral - gamma);
		v *= 0.5f;
		energies[0] += v * (1f + Mathf.Cos(t));
		energies[1] -= v * periodicity * (Mathf.Sin (t));
		energies[2] -= v * periodicity * periodicity * (Mathf.Cos (t));

	}

	public static void ETorsion(float dihedral, float v0, float v1, float v2, float v3, float gamma0, float gamma1, float gamma2, float gamma3, float[] energies) {

		EImproperTorsion(dihedral, v0, gamma0, 1f, energies);
		EImproperTorsionAdditive(dihedral, v1, gamma1, 2f, energies);
		EImproperTorsionAdditive(dihedral, v2, gamma2, 3f, energies);
		EImproperTorsionAdditive(dihedral, v3, gamma3, 4f, energies);
	}

}
