﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PRMReader : MonoBehaviour {

	public IEnumerator ParametersFromPRMFile(string filename, Parameters parameters) {
		string[] lines = FileIO.Readlines (filename);
		Coroutine coroutine = StartCoroutine( ParametersFromPRMLines (lines, parameters));
		
		yield return coroutine;
	}

	//This can be used to read parameters from gaussian input files as well
	public IEnumerator ParametersFromPRMLines(string[] lines, Parameters parameters) {

		int vType;
		int cType;
		int vCutoff;
		int cCutoff;
		float vScale1;
		float vScale2;
		float vScale3;
		float cScale1;
		float cScale2;
		float cScale3;
		string t0;
		string t1;
		string t2;
		string t3;
		float req;
		float keq;
		float v;
		float v0;
		float v1;
		float v2;
		float v3;
		float gamma;
		float gamma0;
		float gamma1;
		float gamma2;
		float gamma3;
		int npaths;
		int periodicity;

		foreach (string line in lines) {
			string[] splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
			if (splitLine.Length == 0) {
				continue;
			}
			
			if (splitLine [0].ToUpper () == "NONBON") {
				vType = int.Parse (splitLine [1]);
				cType = int.Parse (splitLine [2]);
				vCutoff = int.Parse (splitLine [3]);
				cCutoff = int.Parse (splitLine [4]);
				vScale1 = float.Parse (splitLine [5]);
				vScale2 = float.Parse (splitLine [6]);
				vScale3 = float.Parse (splitLine [7]);
				cScale1 = float.Parse (splitLine [8]);
				cScale2 = float.Parse (splitLine [9]);
				cScale3 = float.Parse (splitLine [10]);

				parameters.SetNonBonding (vType, cType, vCutoff, cCutoff, vScale1, vScale2, vScale3, cScale1, cScale2, cScale3);

			}  else if (splitLine [0].ToUpper ().StartsWith ("VDW")) {
				t0 = splitLine [1];
				req = float.Parse (splitLine [2]);
				v = float.Parse (splitLine [3]);

				VdW stretch = new VdW (t0, req, v);
				parameters.AddVdW (stretch);

			} else if (splitLine [0].ToUpper ().StartsWith ("HRMSTR")) {
				t0 = splitLine [1];
				t1 = splitLine [2];
				keq = float.Parse (splitLine [3]);
				req = float.Parse (splitLine [4]);

				Stretch stretch = new Stretch (t0, t1, req, keq);
				parameters.AddStretch (stretch);

			} else if (splitLine [0].ToUpper ().StartsWith ("HRMBND")) {
				 t0 = splitLine [1];
				 t1 = splitLine [2];
				 t2 = splitLine [3];
				keq = float.Parse (splitLine [4]);
				req = float.Parse (splitLine [5]);

				Bend bend = new Bend (t0, t1, t2, req, keq);
				parameters.AddBend (bend);

			} else if (splitLine [0].ToUpper ().StartsWith ("AMBTRS")) {
				t0 = splitLine [1];
				t1 = splitLine [2];
				t2 = splitLine [3];
				t3 = splitLine [4];
				gamma0 = float.Parse (splitLine [5]);
				gamma1 = float.Parse (splitLine [6]);
				gamma2 = float.Parse (splitLine [7]);
				gamma3 = float.Parse (splitLine [8]);
				v0 = float.Parse (splitLine [9]);
				v1 = float.Parse (splitLine [10]);
				v2 = float.Parse (splitLine [11]);
				v3 = float.Parse (splitLine [12]);
				npaths = (int)float.Parse (splitLine [13]);

				Torsion torsion = new Torsion (t0, t1, t2, t3, v0, v1, v2, v3, gamma0, gamma1, gamma2, gamma3, npaths);
				parameters.AddTorsion (torsion);

			} else if (splitLine [0].ToUpper ().StartsWith ("IMPTRS")) {
				t0 = splitLine [1];
				t1 = splitLine [2];
				t2 = splitLine [3];
				t3 = splitLine [4];
				v = float.Parse (splitLine [5]);;
				gamma = float.Parse (splitLine [6]);
				periodicity = (int)float.Parse (splitLine [7]);

				ImproperTorsion improperTorsion = new ImproperTorsion (t0, t1, t2, t3, v, gamma, periodicity);
				parameters.AddImproperTorsion (improperTorsion);

			}
		}

		yield return null;

	}
}
