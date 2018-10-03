using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Parameters : MonoBehaviour {

	public float dielectricConstant;

	public List<VdW> vdws;
	public List<Stretch> stretches;
	public List<Bend> bends;
	public List<Torsion> torsions;
	public List<ImproperTorsion> improperTorsions;
	public NonBonding nonbonding;

	// Use this for initialization
	void Awake () {
		nonbonding = new NonBonding ();
		vdws = new List<VdW> ();
		stretches = new List<Stretch> ();
		bends = new List<Bend> ();
		torsions = new List<Torsion> ();
		improperTorsions = new List<ImproperTorsion> ();
		dielectricConstant = 1f;
	}

	public bool ContainsVdW(VdW otherVdW) {
		foreach (VdW thisVdW in vdws) {
			if (thisVdW.TypeEquivalent(otherVdW))
				return true;
		}
		return false;
	}

	//Find index of type equivalent VdW
	//Return -1 if no type equivalent VdWs are in this parameters set
	public int IndexVdW(VdW otherVdW) {
		for (int index=0;index<vdws.Count;index++) {
			if (vdws[index].TypeEquivalent(otherVdW))
				return index;
		}
		return -1;
	}
		
	public bool ContainsStretch(Stretch otherStretch) {
		foreach (Stretch thisStretch in stretches) {
			if (thisStretch.TypeEquivalent(otherStretch))
				return true;
		}
		return false;
	}

	//Find index of type equivalent stretch
	//Return -1 if no type equivalent stretches are in this parameters set
	public int IndexStretch(Stretch otherStretch) {
		for (int index=0;index<stretches.Count;index++) {
			if (stretches[index].TypeEquivalent(otherStretch))
				return index;
		}
		return -1;
	}

	public bool ContainsBend(Bend otherBend) {
		foreach (Bend thisBend in bends) {
			if (thisBend.TypeEquivalent(otherBend))
				return true;
		}
		return false;
	}

	//Find index of type equivalent bend
	//Return -1 if no type equivalent bends are in this parameters set
	public int IndexBend(Bend otherBend) {
		for (int index=0;index<bends.Count;index++) {
			if (bends[index].TypeEquivalent(otherBend))
				return index;
		}
		return -1;
	}

	public bool ContainsTorsion(Torsion otherTorsion) {
		foreach (Torsion thisTorsion in torsions) {
			if (thisTorsion.TypeEquivalent(otherTorsion))
				return true;
		}
		return false;
	}

	//Find index of type equivalent torsion
	//Return -1 if no type equivalent torsions are in this parameters set
	public int IndexTorsion(Torsion otherTorsion) {
		for (int index=0;index<torsions.Count;index++) {
			if (torsions[index].TypeEquivalent(otherTorsion))
				return index;
		}
		return -1;
	}

	public bool ContainsImproperTorsion(ImproperTorsion otherImproperTorsion) {
		foreach (ImproperTorsion thisImproperTorsion in improperTorsions) {
			if (thisImproperTorsion.TypeEquivalent(otherImproperTorsion))
				return true;
		}
		return false;
	}

	//Find index of type equivalent improper torsion
	//Return -1 if no type equivalent improper torsions are in this parameters set
	public int IndexImproperTorsion(ImproperTorsion otherImproperTorsion) {
		for (int index=0;index<improperTorsions.Count;index++) {
			if (improperTorsions[index].TypeEquivalent(otherImproperTorsion))
				return index;
		}
		return -1;
	}

	public string GetGaussianParamsStr() {
		System.Text.StringBuilder paramsSB = new System.Text.StringBuilder();

		paramsSB.Append(nonbonding.GetGaussianParamStr());
		foreach (Stretch stretch in stretches) 
			paramsSB.Append (stretch.GetGaussianParamStr ());
		foreach (Bend bend in bends) 
			paramsSB.Append (bend.GetGaussianParamStr ());
		foreach (Torsion torsion in torsions)
			paramsSB.Append (torsion.GetGaussianParamStr ());
		foreach (ImproperTorsion improperTorsion in improperTorsions) 
			paramsSB.Append (improperTorsion.GetGaussianParamStr ());
		foreach (VdW vdw in vdws)
			paramsSB.Append (vdw.GetGaussianParamStr ());
		return paramsSB.ToString();
	}

	public void UpdateParameters(Parameters other, bool replace=false) {
		if (replace)
			this.nonbonding = other.nonbonding;
		
		for (int otherIndex = 0; otherIndex < other.vdws.Count; otherIndex++) {
			int thisIndex = IndexVdW (other.vdws [otherIndex]);
			if (thisIndex == -1) {
				this.vdws.Add (other.vdws [otherIndex].Copy ());
			} else if (replace) {
				this.vdws [thisIndex] = other.vdws [otherIndex].Copy ();
			}
		}

		for (int otherIndex = 0; otherIndex < other.stretches.Count; otherIndex++) {
			int thisIndex = IndexStretch (other.stretches [otherIndex]);
			if (thisIndex == -1) {
				this.stretches.Add (other.stretches [otherIndex].Copy ());
			} else if (replace) {
				this.stretches [thisIndex] = other.stretches [otherIndex].Copy ();
			}
		}

		for (int otherIndex = 0; otherIndex < other.bends.Count; otherIndex++) {
			int thisIndex = IndexBend (other.bends [otherIndex]);
			if (thisIndex == -1) {
				this.bends.Add (other.bends [otherIndex].Copy ());
			} else if (replace) {
				this.bends [thisIndex] = other.bends [otherIndex].Copy ();
			}
		}

		for (int otherIndex = 0; otherIndex < other.torsions.Count; otherIndex++) {
			int thisIndex = IndexTorsion (other.torsions [otherIndex]);
			if (thisIndex == -1) {
				Debug.Log ("New Torsion: " + other.torsions [otherIndex].ToString ());
				this.torsions.Add (other.torsions [otherIndex].Copy ());
			} else if (replace) {
				Debug.Log ("Replaced Torsion: " + other.torsions [otherIndex].ToString ());
				this.torsions [thisIndex] = other.torsions [otherIndex].Copy ();
			}
		}

		for (int otherIndex = 0; otherIndex < other.improperTorsions.Count; otherIndex++) {
			int thisIndex = IndexImproperTorsion (other.improperTorsions [otherIndex]);
			if (thisIndex == -1) {
				this.improperTorsions.Add (other.improperTorsions [otherIndex].Copy ());
			} else if (replace) {
				this.improperTorsions [thisIndex] = other.improperTorsions [otherIndex].Copy ();
			}
		}

	}

}

public class VdW {
	public string t;
	public float r;
	public float v;

	public VdW(string t, float r, float v) {
		this.t = t;
		this.r = r;
		this.v = v;
	}

	public bool TypeEquivalent(VdW other) {
		return (this.t == other.t);
	}

	public string GetGaussianParamStr() {
		return string.Format ("VDW {0,-2} {1,7:F4} {2,7:F4}\n", t, r, v);
	}

	public override string ToString () {
		return string.Format ("VdW(t = {0}, r = {1}, v = {2}", t, r, v);
	}

	public VdW Copy() {
		return new VdW (t, r, v);
	}
}

public class Stretch {

	public string t0;
	public string t1;
	public float req;
	public float keq;

	public Stretch(string t0, string t1, float req, float keq) {
		this.t0 = t0;
		this.t1 = t1;
		this.req = req;
		this.keq = keq;
	}

	public bool TypeEquivalent(Stretch other) {
		return (this.t0 == other.t0 && this.t1 == other.t1 || this.t0 == other.t1 && this.t1 == other.t0);
	}

	public string GetGaussianParamStr() {
		return string.Format ("HrmStr1 {0,-2} {1,-2} {2,7:F4} {3,7:F4}\n", t0, t1, req, keq);
	}

	public override string ToString () {
		return string.Format ("Stretch(t0 = {0}, t1 = {1}), req = {2}, keq = {3}", t0, t1, req, keq);
	}

	public Stretch Copy() {
		return new Stretch (t0, t1, req, keq);
	}
}

public class Bend {

	public string t0;
	public string t1;
	public string t2;
	public float req;
	public float keq;

	public Bend(string t0, string t1, string t2, float req, float keq) {
		this.t0 = t0;
		this.t1 = t1;
		this.t2 = t2;
		this.req = req;
		this.keq = keq;
	}

	public bool TypeEquivalent(Bend other) {
		if (this.t1 != other.t1)
			return false;
		if (this.t0 == other.t0 && this.t2 == other.t2)
			return true;
		if (this.t0 == other.t2 && this.t2 == other.t0)
			return true;
		return false;
	}

	public string GetGaussianParamStr() {
		return string.Format ("HrmBnd1 {0,-2} {1,-2} {2,-2} {3,7:F4} {4,7:F4}\n", t0, t1, t2, req, keq);
	}

	public override string ToString () {
		return string.Format ("Bend(t0 = {0}, t1 = {1}), t2 = {2}, req = {3}, keq = {4}", t0, t1, t2, req, keq);
	}

	public Bend Copy() {
		return new Bend (t0, t1, t2, req, keq);
	}
}

public class Torsion {
	public string t0;
	public string t1;
	public string t2;
	public string t3;
	public float v0;
	public float v1;
	public float v2;
	public float v3;
	public float gamma0;
	public float gamma1;
	public float gamma2;
	public float gamma3;
	public int npaths;

	public Torsion(string t0, string t1, string t2, string t3, float v0=0f, float v1=0f, float v2=0f, float v3=0f, float gamma0=0f, float gamma1=0f, float gamma2=0f, float gamma3=0f, int npaths=0) {
		this.t0 = t0;
		this.t1 = t1;
		this.t2 = t2;
		this.t3 = t3;
		this.v0 = v0;
		this.v1 = v1;
		this.v2 = v2;
		this.v3 = v3;
		this.gamma0 = gamma0;
		this.gamma1 = gamma1;
		this.gamma2 = gamma2;
		this.gamma3 = gamma3;
		this.npaths = npaths;
	}

	public void Modify(int periodicity, float v, float gamma) {
		if (periodicity == 1) {
			this.v0 = v;
			this.gamma0 = gamma;
		} else if (periodicity == 2) {
			this.v1 = v;
			this.gamma1 = gamma;
		} else if (periodicity == 3) {
			this.v2 = v;
			this.gamma2 = gamma;
		} else if (periodicity == 4) {
			this.v3 = v;
			this.gamma3 = gamma;
		}
	}

	public  bool TypeEquivalent(Torsion other) {
		if (this.t0 == other.t0 && this.t1 == other.t1  && this.t2 == other.t2  && this.t3 == other.t3)
			return true;
		if (this.t0 == other.t3 && this.t1 == other.t2  && this.t2 == other.t1  && this.t3 == other.t0)
			return true;
		return false;
	}

	public string GetGaussianParamStr() {
		return string.Format ("AmbTrs {0,-2} {1,-2} {2,-2} {3,-2} {4,6:F1} {5,6:F1} {6,6:F1} {7,6:F1} {8,7:F4} {9,7:F4} {10,7:F4} {11,7:F4} {12,4:F1}\n", t0, t1, t2, t3, v0, v1, v2, v3, gamma0, gamma1, gamma2, gamma3, npaths);
	}

	public override string ToString () {
		return string.Format ("Torsion(t0 = {0}, t1 = {1}, t2 = {2}, t3 = {3})", t0, t1, t2, t3);
	}

	public Torsion Copy() {
		return new Torsion (t0, t1, t2, t3, v0, v1, v2, v3, gamma0, gamma1, gamma2, gamma3, npaths);
	}
}

public class ImproperTorsion {
	public string t0;
	public string t1;
	public string t2;
	public string t3;
	public float v;
	public float gamma;
	public int periodicity;

	public ImproperTorsion(string t0, string t1, string t2, string t3, float v, float gamma, int periodicity) {
		this.t0 = t0;
		this.t1 = t1;
		this.t2 = t2;
		this.t3 = t3;
		this.v = v;
		this.gamma = gamma;
		this.periodicity = periodicity;
	}

	public  bool TypeEquivalent(ImproperTorsion other) {
		if (this.t0 == other.t0 && this.t1 == other.t1  && this.t2 == other.t2  && this.t3 == other.t3)
			return true;
		if (this.t0 == other.t3 && this.t1 == other.t2  && this.t2 == other.t1  && this.t3 == other.t0)
			return true;
		return false;
	}

	public string GetGaussianParamStr() {
		return string.Format ("ImpTrs {0,-2} {1,-2} {2,-2} {3,-2} {4,7:F4} {5,6:F1} {6,4:F1}\n", t0, t1, t2, t3, v, gamma, periodicity);
	}

	public override string ToString () {
		return string.Format ("ImproperTorsion(t0 = {0}, t1 = {1}, t2 = {2}, t3 = {3})", t0, t1, t2, t3);
	}

	public ImproperTorsion Copy() {
		return new ImproperTorsion (t0, t1, t2, t3, v, gamma, periodicity);
	}
}

public class NonBonding {
	public int vType;
	public int cType;
	public int vCutoff;
	public int cCutoff;
	public float vScale1;
	public float vScale2;
	public float vScale3;
	public float cScale1;
	public float cScale2;
	public float cScale3;

	public NonBonding(
		int vType=3, int cType=1, int vCutoff=0, int cCutoff=0, float vScale1=0f, float vScale2=0f, float vScale3=0.5f, float cScale1=0f, float cScale2=0f, float cScale3=-1.2f) {
		this.vType = vType;
		this.cType = cType;
		this.vCutoff = vCutoff;
		this.cCutoff = cCutoff;
		this.vScale1 = vScale1;
		this.vScale2 = vScale2;
		this.vScale3 = vScale3;
		this.cScale1 = cScale1;
		this.cScale2 = cScale2;
		this.cScale3 = cScale3;
	}

	public string GetGaussianParamStr() {
		return string.Format ("NonBon {0:D} {1:D} {2:D} {3:D} {4,7:F4} {5,7:F4} {6,7:F4} {7,7:F4} {8,7:F4} {9,7:F4}\n", vType, cType, vCutoff, cCutoff, vScale1, vScale2, vScale3, cScale1, cScale2, cScale3);
	}

	public override string ToString () {
		return string.Format ("NonBonding(vType = {0}, cType = {1}, vCutoff = {2}, cCutoff = {3}, vScale1 = {4}, vScale2 = {5}, vScale3 = {6}, cScale1 = {7}, cScale2 = {8}, cScale3 = {9})", vType, cType, vCutoff, cCutoff, vScale1, vScale2, vScale3, cScale1, cScale2, cScale3);
	}
}