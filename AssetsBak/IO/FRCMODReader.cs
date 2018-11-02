using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FRCMODReader : MonoBehaviour {
	
	public Parameters parametersPrefab;
	public char[] trimChars = new char[] { ' ', '-' };

	private int lineNumber;
	private string line;

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
	string t;
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
	float radius;

	public Parameters ParametersFromFRCMODFile(string filename) {

		Parameters parameters = Instantiate<Parameters>(parametersPrefab);

		string[] lines = FileIO.Readlines (filename);
		int maxLines = lines.Length;

		List<string> sections = new List<string> {
			"MASS",
			"BOND",
			"ANGL",
			"DIHE",
			"IMPR",
			"NONB",
			"IPOL",
			"CMAP",
			"%"
		};
		string sectionName;
		int sectionNum = 0;
		int numSections = sections.Count;

		string section = "";
		lineNumber = 0;

		string restOfLine;
		string[] splitLine;


		Stretch newStretchParam;
		Bend newBendParam;
		Torsion newTorsionParam;
		ImproperTorsion newImproperTorsionParam;
		VdW newVdWParam;

		Torsion existingTorsionParam;
		int existingTorsionNum;
		bool torsionModified;

		int offset;

		while (lineNumber < maxLines) {

			line = lines [lineNumber];

			if (line == "")
				continue;

			for (sectionNum = 0; sectionNum < numSections; sectionNum++) {
				sectionName = sections [sectionNum];
				if (line.StartsWith (sectionName)) {
					section = sectionName;
					goto SKIP;
				}
			}

			//Ignore MASS section

			//Parse Bond
			if (section == "BOND") {

				t0 = GetAmberFromString (line, 0);
				t1 = GetAmberFromString (line, 3);

				restOfLine = line.Substring (5, line.Length - 6);
				splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				keq = float.Parse (splitLine [0]);
				req = float.Parse (splitLine [1]);

				newStretchParam = new Stretch (t0, t1, req, keq);
				parameters.stretches.Add (newStretchParam);

			} else if (section == "ANGL") {

				t0 = GetAmberFromString (line, 0);
				t1 = GetAmberFromString (line, 3);
				t2 = GetAmberFromString (line, 6);

				restOfLine = line.Substring (8, line.Length - 9);
				splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				keq = float.Parse (splitLine [0]);
				req = float.Parse (splitLine [1]);

				newBendParam = new Bend (t0, t1, t2, req, keq);
				parameters.bends.Add (newBendParam);

			} else if (section == "DIHE") {

				t0 = GetAmberFromString (line, 0);
				t1 = GetAmberFromString (line, 3);
				t2 = GetAmberFromString (line, 6);
				t3 = GetAmberFromString (line, 9);

				newTorsionParam = new Torsion (t0, t1, t2, t3);

				restOfLine = line.Substring (11, line.Length - 12);
				splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				npaths = int.Parse (splitLine [0]);
				v = float.Parse (splitLine [1]);
				gamma = float.Parse (splitLine [2]);
				periodicity = (int)Mathf.Abs (float.Parse (splitLine [3]));

				//The same param can be defined on multiple lines
				torsionModified = false;

				for (existingTorsionNum = 0; existingTorsionNum < parameters.torsions.Count; existingTorsionNum++) {
					existingTorsionParam = parameters.torsions [existingTorsionNum];
					if (existingTorsionParam.TypeEquivalent (newTorsionParam)) {
						existingTorsionParam.Modify (periodicity, v, gamma);
						torsionModified = true;
					}
				}


				if (!torsionModified) {
					newTorsionParam.npaths = npaths;
					if (periodicity == 1) {
						newTorsionParam.v0 = v;
						newTorsionParam.gamma0 = gamma;
					} else if (periodicity == 2) {
						newTorsionParam.v1 = v;
						newTorsionParam.gamma1 = gamma;
					} else if (periodicity == 3) {
						newTorsionParam.v2 = v;
						newTorsionParam.gamma2 = gamma;
					} else if (periodicity == 4) {
						newTorsionParam.v3 = v;
						newTorsionParam.gamma3 = gamma;
					}

					parameters.torsions.Add (newTorsionParam);
				}

			} else if (section == "IMPR") {

				t0 = GetAmberFromString (line, 0);
				t1 = GetAmberFromString (line, 3);
				t2 = GetAmberFromString (line, 6);
				t3 = GetAmberFromString (line, 9);

				restOfLine = line.Substring (11, line.Length - 12);
				splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				offset = 1;

				//Sometimes npaths is included here, even though ignored. Sometimes it's not.
				if (splitLine [0].Contains ("."))
					offset = 0;

				v = float.Parse (splitLine [offset]);
				gamma = float.Parse (splitLine [offset + 1]);
				periodicity = (int)Mathf.Abs (float.Parse (splitLine [offset + 2]));

				newImproperTorsionParam = new ImproperTorsion (t0, t1, t2, t3, v, gamma, periodicity);

				parameters.improperTorsions.Add (newImproperTorsionParam);

			} else if (section == "NONB") {
				splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				t = splitLine [0];
				radius = float.Parse (splitLine [1]);
				v = float.Parse (splitLine [2]);

				newVdWParam = new VdW (t, radius, v);

				parameters.vdws.Add (newVdWParam);

			}

			SKIP:
			;
		}

		return parameters;

	}

	private string GetAmberFromString(string str, int start) {
		string amber = str.Substring (start, 2).Trim (trimChars);
		if (amber.ToUpper () == "X") {
			return "*";
		}
		return amber;
	}
}
