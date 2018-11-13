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

	public static void EVdWAmberOld(float r, float v, float req, float[] energies){

		float _r = v * req / r;
		float _r2 = _r * _r;
		float _r6 = _r2 * _r2 * _r2;
		float _r12 = _r6 * _r6;

		energies[0] = _r12 - 2 * _r6;
		energies[1] = 12f * (_r6 * _r - _r12 * _r) / r;
		energies[2] = 12f * (13f * _r12 * _r - 7f * _r6 * _r2) / (r * r);
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
		energies[1] = - 2f * energies[0] * _r;
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
		float t = periodicity * dihedral - gamma;
		v *= 0.5f;
		energies[0] = v * (1f + Mathf.Cos(t));
		energies[1] = - v * periodicity * (Mathf.Sin (t));
		energies[2] = - v * periodicity * periodicity * (Mathf.Cos (t));
	}

	public static void EImproperTorsionAdditive(float dihedral, float v, float gamma, float periodicity, float[] energies) {
		if (v == 0f) {return;}
		float t = periodicity * dihedral - gamma;
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

	public static void ItemFromArray(float[,] A, int index, float[] result) {
		result[0] = A[index, 0];
		result[1] = A[index, 1];
		result[2] = A[index, 2];
	}

	public static void AverageFromArray(float[] v, float[,] A, int[] indices, int index, int numIndices, int itemNum) {
		
		v[0] = 0f;
		v[1] = 0f;
		v[2] = 0f;
		for (itemNum = 0; itemNum < numIndices; itemNum++) {
			index = indices[index];
			v[0] += A[index, 0];
			v[1] += A[index, 1];
			v[2] += A[index, 2];
		}
		v[0] /= numIndices;
		v[1] /= numIndices;
		v[2] /= numIndices;
	}

	public static void AverageFromArray(float[,] A, int index0, int index1, float[] result) {
		
		result[0] = (A[index0, 0] + A[index1, 0]) / 2f;
		result[1] = (A[index0, 1] + A[index1, 1]) / 2f;
		result[2] = (A[index0, 2] + A[index1, 2]) / 2f;
	}

	public static void VectorFromArray(float[,] A, int from, int to, float[] result) {
		result[0] = A[to, 0] - A[from, 0];
		result[1] = A[to, 1] - A[from, 1];
		result[2] = A[to, 2] - A[from, 2];
	}

	public static void VectorFromVector3(Vector3 v, float[] result) {
		result[0] = v[0];
		result[1] = v[1];
		result[2] = v[2];
	}

	public static void Copy(float[] from, float[] to) {
		to[0] = from[0];
		to[1] = from[1];
		to[2] = from[2];
	}

	public static bool IsEqual(float[] v0, float[] v1) {
		if (v0[0] != v1[0]) return false;
		if (v0[1] != v1[1]) return false;
		if (v0[2] != v1[2]) return false;
		return true;
	}

	public static float Magnitude3(float[] v) {
		return Mathf.Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	}

	public static void Add3(float[] v0, float[] v1, float[] result) {
		result[0] = v0[0] + v1[0];
		result[1] = v0[1] + v1[1];
		result[2] = v0[2] + v1[2];
	}

	public static void Add3(float[] v0, float[] v1) {
		v0[0] += v1[0];
		v0[1] += v1[1];
		v0[2] += v1[2];
	}

	public static void Subtract3(float[] v0, float[] v1, float[] result) {
		result[0] = v0[0] - v1[0];
		result[1] = v0[1] - v1[1];
		result[2] = v0[2] - v1[2];
	}

	public static void Subtract3(float[] v0, float[] v1) {
		v0[0] -= v1[0];
		v0[1] -= v1[1];
		v0[2] -= v1[2];
	}

	public static void Divide3(float[] v, float d, float[] result) {
		result[0] = v[0] / d;
		result[1] = v[1] / d;
		result[2] = v[2] / d;
	}

	public static void Divide3(float[] v, float d) {
		v[0] /= d;
		v[1] /= d;
		v[2] /= d;
	}

	public static void Multiply3(float[] v, float m, float[] result) {
		result[0] = v[0] * m;
		result[1] = v[1] * m;
		result[2] = v[2] * m;
	}

	public static void Multiply3(float[] v, float m) {
		v[0] *= m;
		v[1] *= m;
		v[2] *= m;
	}

	public static void Normalise3(float[] v) {
		float mag = Magnitude3(v);
		Divide3(v, mag);
	}

	public static void Normalise3(float[] v, float[] result) {
		float mag = Magnitude3(v);
		Divide3(v, mag, result);
	}

	public static void Cross3(float[] v0, float[] v1, float[] result) {
		result[0] = v0[1] * v1[2] - v0[2] * v1[1];
		result[1] = v0[2] * v1[0] - v0[0] * v1[2];
		result[2] = v0[0] * v1[1] - v0[1] * v1[0];
	}

	public static float Dot3(float[] v0, float[] v1) {
		return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
	}

	public static float SignedAngleRad3(float[] v0, float[] v1, float[] v01) {
		float[] cross = new float[3];
		float x = 0f;
		float y = 0f;

		return SignedAngleRad3(v0, v1, v01, cross, x, y);
	}

	public static float SignedAngleRad3(float[] v0, float[] v1, float[] v01, float[] cross, float x, float y) {
		Cross3(v0, v1, cross);
		x = Dot3(v0, v1);
		y = Dot3(v01, cross);
		return Mathf.Atan2(y, x);
	}

	public static float UnsignedAngleRad3(float[] v0_, float[] v1_) {
		float dot = 0f;
		return UnsignedAngleRad3(v0_, v1_, dot);
	}

	public static float UnsignedAngleRad3(float[] v0_, float[] v1_, float dot) {
		dot = Dot3(v0_, v1_);
		return Mathf.Acos(dot);
	}


	//BENDS
	public static void GetBendForce(int a0, int a1, int a2, float[,] positions, float[] energies, float[,] forces, float bendAEq, float bendKEq) {
		int c;

		float[] perp = new float[3];
		float[] v10 = new float[3];
		float[] v12 = new float[3];
		float dot = 0f;

		float r10;
		float r12;

		float[] force0 = new float[3];
		float[] force2 = new float[3];

		Mathematics.VectorFromArray(positions, a1, a0, v10);
		Mathematics.VectorFromArray(positions, a1, a2, v12);
		r10 = Mathematics.Magnitude3(v10);
		r12 = Mathematics.Magnitude3(v12);

		Mathematics.Divide3(v10, r10);
		Mathematics.Divide3(v12, r12);

		Mathematics.EAngle(
			Mathematics.UnsignedAngleRad3(v10, v12, dot) - bendAEq,
			bendKEq,
			energies
		);

		Mathematics.Cross3(v10, v12, perp);
		Mathematics.Cross3(v10, perp, force0);
		Mathematics.Cross3(v12, perp, force2);

		Mathematics.Multiply3(force0, - energies[1] / r10);
		Mathematics.Multiply3(force2, energies[1] / r12);

		for (c = 0; c < 3; c++) {
			forces[a0, c] += force0[c];
			forces[a1, c] -= force0[c] + force2[c];
			forces[a2, c] += force2[c];
		}
	}

	public static void GetBendForce(int a0, int a1, int a2, float[,] positions, float[] energies, float[] force0, float[] force1, float[] force2, float bendAEq, float bendKEq) {
		int c;

		float[] perp = new float[3];
		float[] v10 = new float[3];
		float[] v12 = new float[3];
		float dot = 0f;

		float r10;
		float r12;

		VectorFromArray(positions, a1, a0, v10);
		VectorFromArray(positions, a1, a2, v12);
		r10 = Magnitude3(v10);
		r12 = Magnitude3(v12);

		Divide3(v10, r10);
		Divide3(v12, r12);

		EAngle(
			UnsignedAngleRad3(v10, v12, dot) - bendAEq,
			bendKEq,
			energies
		);

		Cross3(v10, v12, perp);
		Cross3(v10, perp, force0);
		Cross3(v12, perp, force2);

		Multiply3(force0, - energies[1] / r10);
		Multiply3(force2, energies[1] / r12);

		for (c = 0; c < 3; c++) {
			force1[c] = - force0[c] - force2[c];
		}
	}

	public static void GetBendForce(float[] p0, float[] p1, float[] p2, float[] energies, float[] force0, float[] force1, float[] force2, float bendAEq, float bendKEq) {
		int c;

		float[] perp = new float[3];
		float[] v10 = new float[3];
		float[] v12 = new float[3];
		float dot = 0f;

		float r10;
		float r12;

		Subtract3(p0, p1, v10);
		Subtract3(p2, p1, v12);

		r10 = Magnitude3(v10);
		r12 = Magnitude3(v12);

		Divide3(v10, r10);
		Divide3(v12, r12);

		EAngle(
			UnsignedAngleRad3(v10, v12, dot) - bendAEq,
			bendKEq,
			energies
		);

		Cross3(v10, v12, perp);
		Cross3(v10, perp, force0);
		Cross3(v12, perp, force2);

		Multiply3(force0, - energies[1] / r10);
		Multiply3(force2, energies[1] / r12);

		for (c = 0; c < 3; c++) {
			force1[c] = - force0[c] - force2[c];
		}
	}

}
