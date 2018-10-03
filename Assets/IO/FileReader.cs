using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class FileReader: MonoBehaviour {

	//Prefabs
	public Atom atomPrefab;
	public Atoms atomsPrefab;
	public Parameters parametersPrefab;
	public GaussianCalculator gaussianCalculatorPrefab;

	public Data data;
	public Settings globalSettings;

	public char[] trimChars = new char[] { ' ', '-' };

	public Atoms AtomsFromGaussianInput(string path) {

		Atoms atoms = Instantiate<Atoms> (atomsPrefab);
		atoms.globalSettings = globalSettings;
		atoms.setSourceFile(path, "com");

		atoms.name = Path.GetFileName (path);
		GaussianCalculator gc = Instantiate<GaussianCalculator> (gaussianCalculatorPrefab, atoms.transform);
		atoms.gaussianCalculator = gc;

		string[] lines = FileIO.Readlines (path, "!");
		int lineNumber = 0;
		int maxLines = lines.Length;

		//Link 0 commands
		while (lineNumber < maxLines) {
			string line = lines [lineNumber].Trim();
			if (line.StartsWith ("#")) {
				break;

			} else if (line.ToUpper().StartsWith ("%CHK")) {
				gc.checkpointPath = GetValueFromPair (line);
			} else if (line.ToUpper().StartsWith ("%OLDCHK")) {
				gc.oldCheckpointPath = GetValueFromPair (line);
			} else if (line.ToUpper().StartsWith ("%MEM")) {
				gc.jobMemoryMB = GetMemoryMB (GetValueFromPair (line));
			} else if (line.ToUpper().StartsWith ("%NPROC")) {
				gc.nProcessors = int.Parse (GetValueFromPair (line));
			} else if (line.ToUpper().StartsWith ("%KJOB")) {
				string[] splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
				gc.killJobLink = splitLine [1];

				if (splitLine.Length > 2) {
					gc.killJobAfter = int.Parse (splitLine [2]);
				} else {
					gc.killJobAfter = 1;
				}
			}

			lineNumber++;
		}


		List<string> keywordsList = new List<string> ();

		//Keywords
		while (lineNumber < maxLines) {
			string line = lines [lineNumber].Trim();
			if (line == "") {
				break;
			} 

			string[] splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
			foreach (string keywordItem in splitLine) {
				keywordsList.Add (keywordItem.ToLower());
			}
			lineNumber++;
		}

		//Parse keywords
		foreach (string keywordItem in keywordsList) {

			//Print level
			if (keywordItem.StartsWith ("#")) {
				gc.gaussianPrintLevel = keywordItem.Replace ("#", "");

				//Method
			} else if (keywordItem.StartsWith ("oniom")) {
				string oniomMethodsStr = GetStringInParentheses (keywordItem);
				string oniomOptionsStr = GetStringInParentheses (GetValueFromPair (keywordItem, checkEnclosed: true));

				string[] oniomOptions = oniomOptionsStr.Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries);

				foreach (string oniomOption in oniomOptions)
					gc.oniomOptions.Add (oniomOption);

				string[] methods = oniomMethodsStr.Split (new []{ ":" }, System.StringSplitOptions.RemoveEmptyEntries);

				string[] highMBO = GetMethodFromString (methods [0]);
				gc.layers.Add( new Layer(highMBO [0], highMBO [1], highMBO [2], 'H'));


				if (methods.Length == 2) {
					string[] lowMBO = GetMethodFromString (methods [1]);
					gc.layers.Add( new Layer(lowMBO [0], lowMBO [1], lowMBO [2], 'L'));
				} else if (methods.Length == 3) {
					string[] mediumMBO = GetMethodFromString (methods [1]);
					gc.layers.Add( new Layer(mediumMBO [0], mediumMBO [1], mediumMBO [2], 'M'));

					string[] lowMBO = GetMethodFromString (methods [2]);
					gc.layers.Add( new Layer(lowMBO [0], lowMBO [1], lowMBO [2], 'L'));
				}
			} else if (keywordItem.StartsWith ("guess")) {
				string guessOptionsStr = GetStringInParentheses (GetValueFromPair (keywordItem, checkEnclosed: true));
				string[] guessOptions = guessOptionsStr.Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries);
				foreach (string guessOption in guessOptions)
					gc.guessOptions.Add (guessOption);
			}  else if (keywordItem.StartsWith ("geom")) {
				string geomOptionsStr = GetStringInParentheses (GetValueFromPair (keywordItem, checkEnclosed: true));
				string[] geomOptions = geomOptionsStr.Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries);
				foreach (string geomOption in geomOptions)
					gc.geomOptions.Add (geomOption);
			}
		}

		//Title
		List<string> titleLines = new List<string>();

		if (!gc.geomOptions.Contains ("allcheck")) {
			while (lineNumber < maxLines) {
				string line = lines [lineNumber].Trim ();
				if (line != "") {
					titleLines.Add (line);
				} else if (titleLines.Count != 0) {
					break;
				}

				lineNumber++;
			}
		}
		gc.title = string.Join ("\n", titleLines.ToArray ());

		//Layer Charge/Multiplicity
		if (!gc.geomOptions.Contains ("allcheck")) {
			while (lineNumber < maxLines) {
				string line = lines [lineNumber].Trim ();
				if (line != "") {
					string[] cmStr = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

					int numLayersSpecified = cmStr.Length / 2;

					for (int layerNum = 0; layerNum < gc.layers.Count; layerNum++) {

						//Not all layers need to have charge and multiplicity specified
						if (layerNum < numLayersSpecified) {
							gc.layers [layerNum].charge = int.Parse (cmStr [2 * layerNum]);
							gc.layers [layerNum].multiplicity = int.Parse (cmStr [2 * layerNum + 1]);
						} else {
							gc.layers [layerNum].charge = gc.layers [layerNum - 1].charge;
							gc.layers [layerNum].multiplicity = gc.layers[layerNum].multiplicity;
						}
					}

					lineNumber++;
					break;
				} 
				lineNumber++;
			}
		}

		//Atoms
		if (!(gc.geomOptions.Contains ("allcheck") && gc.geomOptions.Contains ("check") && gc.geomOptions.Contains ("checkpoint"))) {

			int index = 0;
			while (lineNumber < maxLines) {
				string line = lines [lineNumber].Trim ();
				if (line == "") {
					lineNumber++;
					break;
				}
					
				Atom atom = AtomFromGaussianInputLine (line);
				atom.index = index;
				atoms.atomList.Add (atom);
				atom.transform.parent = atoms.atomsHolder.transform;

				index++;
				lineNumber++;
			}
		}

		//Connectivity
		atoms.graph.SetAtoms (atoms);
		if (gc.geomOptions.Contains ("connectivity")) {
			while (lineNumber < maxLines) {
				string line = lines [lineNumber].Trim ();
				if (line == "")
					break;
				ConnectionsFromGaussianInputLine (line, atoms.graph);

				lineNumber++;
			}
		}
			
		//Parameters
		List<string> parameterLines = new List<string>();

		while (lineNumber < maxLines) {
			string line = lines [lineNumber].Trim ();
			parameterLines.Add (line);
			lineNumber++;
		}

		atoms.parameters = ParametersFromPRMLines (parameterLines.ToArray ());
		atoms.parameters.transform.parent = atoms.transform;

		return atoms;

	}

	public Atoms AtomsFromPDBFile(string path, string alternateLocationID="A") {

		Atoms atoms = Instantiate<Atoms> (atomsPrefab);
		atoms.globalSettings = globalSettings;
		atoms.setSourceFile(path, "pdb");

		atoms.name = Path.GetFileName (path);

		int index = 0;

		string[] lines = FileIO.Readlines (path, "#");
		foreach (string line in lines) {
			if (
				(
				    line.StartsWith ("ATOM") ||
				    line.StartsWith ("HETATM")
				) && (
				    line [16].ToString () == " " ||
				    line [16].ToString () == alternateLocationID
				)) {
				Atom atom = Instantiate (atomPrefab, atoms.atomsHolder.transform);

				//See if we can use the last columns to get Atom
				bool amberUsed = false;
				string amber = line.Substring (76, 2).Trim ();
				if (amber != string.Empty) {
					atom.amberName = amber;
					foreach (Data.AmberToElement amberToElement in data.amberToElementData) {
						if (amberToElement.amberName == amber) {
							atom.element = amberToElement.element;
							atom.formalCharge = amberToElement.formalCharge;
							amberUsed = true;
						}
					}
				}

				//PDB
				atom.pdbName = line.Substring (12, 4).Trim ();

				//Use PDB to get Atom
				if (!amberUsed) {
					//If 12th column is a letter, it is a metal
					if (char.IsLetter (line [12])) {
						string element = line.Substring (12, 2).Trim ();
						atom.element = ToAlpha(element);
					} else {
						string element = line.Substring (12, 4).Trim ();
						atom.element = ToAlpha(element) [0].ToString ();
					}
				}

				//Residue
				atom.residueName = line.Substring (17, 3).Trim ();
				atom.residueNumber = int.Parse (line.Substring (22, 4).Trim ());
				atom.chainID = line [21].ToString ();

				//Position
				Vector3 position = new Vector3 ();
				position.x = float.Parse (line.Substring (30, 8));
				position.y = float.Parse (line.Substring (38, 8));
				position.z = float.Parse (line.Substring (46, 8));
				atom.transform.position = position;

				string chargeStr = line [79].ToString ();
				if (chargeStr != " ") {
					if (line [78].ToString () == "-") {
						atom.formalCharge = -int.Parse (chargeStr);
					} else {
						atom.formalCharge = int.Parse (chargeStr);
					}
				}

				atom.index = index;
				index++;
				atoms.atomList.Add (atom);

			}
		}
		atoms.graph.SetAtoms (atoms);
		return atoms;
	}

	public Parameters ParametersFromFRCMODFile(string filename) {

		Parameters parameters = Instantiate<Parameters>(parametersPrefab);

		string[] lines = FileIO.Readlines (filename);

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

		string section = "";

		foreach (string line in lines) {

			if (line == "")
				continue;
			
			foreach (string sectionName in sections) {
				if (line.StartsWith (sectionName)) {
					section = sectionName;
					goto SKIP;
				}
			}

			//Ignore MASS section

			//Parse Bond
			if (section == "BOND") {

				string t0 = GetAmberFromString (line, 0);
				string t1 = GetAmberFromString (line, 3);

				string restOfLine = line.Substring (5, line.Length - 6);
				string[] splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				float keq = float.Parse (splitLine [0]);
				float req = float.Parse (splitLine [1]);

				Stretch stretchParam = new Stretch (t0, t1, req, keq);
				parameters.stretches.Add (stretchParam);

			} else if (section == "ANGL") {

				string t0 = GetAmberFromString (line, 0);
				string t1 = GetAmberFromString (line, 3);
				string t2 = GetAmberFromString (line, 6);

				string restOfLine = line.Substring (8, line.Length - 9);
				string[] splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				float keq = float.Parse (splitLine [0]);
				float req = float.Parse (splitLine [1]);

				Bend bendParam = new Bend (t0, t1, t2, req, keq);
				parameters.bends.Add (bendParam);
				
			} else if (section == "DIHE") {

				string t0 = GetAmberFromString (line, 0);
				string t1 = GetAmberFromString (line, 3);
				string t2 = GetAmberFromString (line, 6);
				string t3 = GetAmberFromString (line, 9);

				Torsion torsionParam = new Torsion (t0, t1, t2, t3);

				string restOfLine = line.Substring (11, line.Length - 12);
				string[] splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				int npaths = int.Parse (splitLine [0]);
				float v = float.Parse (splitLine [1]);
				float gamma = float.Parse (splitLine [2]);
				int nt = (int)Mathf.Abs (float.Parse (splitLine [3]));

				//The same param can be defined on multiple lines
				bool modified = false;
				foreach (Torsion existingTorsion in parameters.torsions) {
					if (existingTorsion.TypeEquivalent (torsionParam)) {
						existingTorsion.Modify (nt, v, gamma);
						modified = true;
					}
				}
					

				if (!modified) {
					if (nt == 1) {
						torsionParam.v0 = v;
						torsionParam.gamma0 = gamma;
					} else if (nt == 2) {
						torsionParam.v1 = v;
						torsionParam.gamma1 = gamma;
					} else if (nt == 3) {
						torsionParam.v2 = v;
						torsionParam.gamma2 = gamma;
					} else if (nt == 4) {
						torsionParam.v3 = v;
						torsionParam.gamma3 = gamma;
					}

					parameters.torsions.Add (torsionParam);
				}

			} else if (section == "IMPR") {

				string t0 = GetAmberFromString (line, 0);
				string t1 = GetAmberFromString (line, 3);
				string t2 = GetAmberFromString (line, 6);
				string t3 = GetAmberFromString (line, 9);

				Torsion torsionParam = new Torsion (t0, t1, t2, t3);

				string restOfLine = line.Substring (11, line.Length - 12);
				string[] splitLine = restOfLine.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				int offset = 1;

				//Sometimes npaths is included here, even though ignored. Sometimes it's not.
				if (splitLine [0].Contains ("."))
					offset = 0;

				float v = float.Parse (splitLine [offset]);
				float gamma = float.Parse (splitLine [offset + 1]);
				int nt = (int)Mathf.Abs (float.Parse (splitLine [offset + 2]));

				ImproperTorsion improperTorsionParam = new ImproperTorsion (t0, t1, t2, t3, v, gamma, nt);

				parameters.improperTorsions.Add (improperTorsionParam);

			} else if (section == "NONB") {
				string[] splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				string t = splitLine [0];
				float r = float.Parse (splitLine [1]);
				float v = float.Parse (splitLine [2]);

				VdW vdwParam = new VdW (t, r, v);

				parameters.vdws.Add (vdwParam);

			}

			SKIP:
			;
		}

		return parameters;

	}

	public Parameters ParametersFromPRMFile(string filename) {
		string[] lines = FileIO.Readlines (filename);
		return ParametersFromPRMLines (lines);
	}

	//This can be used to read parameters from gaussian input files as well
	public Parameters ParametersFromPRMLines(string[] lines) {

		Parameters parameters = Instantiate<Parameters>(parametersPrefab);

		foreach (string line in lines) {
			string[] splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
			if (splitLine.Length == 0) {
				continue;
			} if (splitLine [0].ToUpper () == "NONBON") {
				int vType = int.Parse (splitLine [1]);
				int cType = int.Parse (splitLine [2]);
				int vCutoff = int.Parse (splitLine [3]);
				int cCutoff = int.Parse (splitLine [4]);
				float vScale1 = float.Parse (splitLine [5]);
				float vScale2 = float.Parse (splitLine [6]);
				float vScale3 = float.Parse (splitLine [7]);
				float cScale1 = float.Parse (splitLine [8]);
				float cScale2 = float.Parse (splitLine [9]);
				float cScale3 = float.Parse (splitLine [10]);

				parameters.nonbonding = new NonBonding (vType, cType, vCutoff, cCutoff, vScale1, vScale2, vScale3, cScale1, cScale2, cScale3);

			} else if (splitLine [0].ToUpper ().StartsWith ("HRMSTR")) {
				string t0 = splitLine [1];
				string t1 = splitLine [2];
				float req = float.Parse (splitLine [3]);
				float keq = float.Parse (splitLine [4]);

				Stretch stretch = new Stretch (t0, t1, req, keq);
				parameters.stretches.Add (stretch);

			} else if (splitLine [0].ToUpper ().StartsWith ("HRMBND")) {
				string t0 = splitLine [1];
				string t1 = splitLine [2];
				string t2 = splitLine [3];
				float req = float.Parse (splitLine [4]);
				float keq = float.Parse (splitLine [5]);

				Bend bend = new Bend (t0, t1, t2, req, keq);
				parameters.bends.Add (bend);

			} else if (splitLine [0].ToUpper ().StartsWith ("AMBTRS")) {
				string t0 = splitLine [1];
				string t1 = splitLine [2];
				string t2 = splitLine [3];
				string t3 = splitLine [4];
				float v0 = float.Parse (splitLine [5]);
				float v1 = float.Parse (splitLine [6]);
				float v2 = float.Parse (splitLine [7]);
				float v3 = float.Parse (splitLine [8]);
				float gamma0 = float.Parse (splitLine [9]);
				float gamma1 = float.Parse (splitLine [10]);
				float gamma2 = float.Parse (splitLine [11]);
				float gamma3 = float.Parse (splitLine [12]);
				int npaths = (int)float.Parse (splitLine [13]);

				Torsion torsion = new Torsion (t0, t1, t2, t3, v0, v1, v2, v3, gamma0, gamma1, gamma2, gamma3, npaths);
				parameters.torsions.Add (torsion);

			} else if (splitLine [0].ToUpper ().StartsWith ("IMPTRS")) {
				string t0 = splitLine [1];
				string t1 = splitLine [2];
				string t2 = splitLine [3];
				string t3 = splitLine [4];
				float v = float.Parse (splitLine [5]);;
				float gamma = float.Parse (splitLine [6]);
				int periodicity = (int)float.Parse (splitLine [7]);

				ImproperTorsion improperTorsion = new ImproperTorsion (t0, t1, t2, t3, v, gamma, periodicity);
				parameters.improperTorsions.Add (improperTorsion);

			}
		}

		return parameters;

	}

	//Tools

	private string GetAmberFromString(string str, int start) {
		string amber = str.Substring (start, 2).Trim (trimChars);
		if (amber.ToUpper () == "X") {
			return "*";
		}
		return amber;
	}

	private string GetValueFromPair (string inputStr, string delimiter = "=", bool checkEnclosed=false, char openChar='(', char closeChar=')') {

		//There's the possibility of the delimiter being inside parentheses
		//e.g. key(option1=value1)=value
		if (checkEnclosed) {

			for (int i = 0; i < inputStr.Length; i++) {
				if (inputStr [i].ToString () == delimiter) {
					if (!StringIndexEnclosed (inputStr, i, openChar, closeChar)) {
						return inputStr.Substring (i + 1, inputStr.Length - i - 1);
					}
				}
			}
			return "";
		} 

		string[] keyValue = inputStr.Split (new []{ delimiter }, System.StringSplitOptions.RemoveEmptyEntries);
		return keyValue [1];
	}

	private int GetMemoryMB(string memoryString) {
		int oldMem = int.Parse(ToNumber (memoryString));
		string oldUnits = ToAlpha (memoryString).ToUpper();

		if (oldUnits == "KB") {
			return oldMem / 1024;
		} else if (oldUnits == "MB") {
			return oldMem;
		} else if (oldUnits == "GB") {
			return oldMem * 1024;
		} else if (oldUnits == "TB") {
			return oldMem * 1048576;
		} else {
			Debug.LogError(string.Format("Memory units '{0}' not recognised.", oldUnits));
		}

		return 4000;
	}

	private string GetStringInParentheses(string inputStr) {

		//Get the indices of the pair of parentheses
		int startIndex = -1;
		int endIndex = -1;

		//How many parenthesis pairs are we in?
		int depth = 0;

		for (int i = 0; i < inputStr.Length; i++) {
			switch (inputStr [i]) {
			case '(':
				if (depth == 0 && startIndex == -1)
					startIndex = i + 1;
				depth += 1;
				break;
			case ')':
				if (depth == 1 && endIndex == -1)
					endIndex = i - 1;
				depth -= 1;
				break;
			}
		}

		//No parentheses
		if (startIndex == -1)
			return inputStr;

		return inputStr.Substring (startIndex, endIndex - startIndex + 1);
	}

	private string[] GetMethodFromString(string inputStr) {
		string[] methodBasisSplit = inputStr.Split (new []{ "/" }, System.StringSplitOptions.RemoveEmptyEntries);

		string method = methodBasisSplit [0];
		string basis = "";
		string options = "";

		if (methodBasisSplit.Length == 1) {
			string[] methodOptionsSplit = method.Split (new []{ "=" }, System.StringSplitOptions.RemoveEmptyEntries);
			method = methodOptionsSplit [0];
			if (methodOptionsSplit.Length == 2)
				options = methodOptionsSplit [1];
		} else {
			string[] basisOptionsSplit = methodBasisSplit [1].Split (new []{ "/" }, System.StringSplitOptions.RemoveEmptyEntries);
			basis = basisOptionsSplit [0];
			if (basisOptionsSplit.Length == 2)
				options = basisOptionsSplit [1];
		}

		return new string[]{ method, basis, options };
	}

	private bool StringIndexEnclosed(string inputStr, int index, char openChar='(', char closeChar=')') {

		int depth = 0;

		for (int i = 0; i < inputStr.Length; i++) {
			char chr = inputStr [i];

			if (index == i) {
				return (depth != 0);
			} else if (chr == openChar) {
				depth += 1;
			} else if (chr == closeChar) {
				depth -= 1;
			}
		}

		return false;
	}

	private string ToAlpha(string inputStr) {
		return System.Text.RegularExpressions.Regex.Replace (inputStr, @"[^a-zA-Z -]", string.Empty);
	}

	private string ToNumber(string inputStr) {
		return System.Text.RegularExpressions.Regex.Replace (inputStr, @"[^0-9 -]", string.Empty);
	}

	private Atom AtomFromGaussianInputLine(string inputStr) {

		Atom atom = Instantiate<Atom> (atomPrefab);
		int splitIndex = 0;

		string[] splitLine = inputStr.ToUpper().Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

		string atomSpec = splitLine [splitIndex];

		//Get Element, Amber and partial charge
		string[] splitSpec = atomSpec.ToUpper().Split (new []{ "-" }, 3, System.StringSplitOptions.RemoveEmptyEntries);

		if (splitSpec.Length > 1) {
			atom.element = splitSpec [0];
			float charge;

			bool chargeSpecified = float.TryParse (splitSpec [1], out charge);

			if (splitSpec.Length > 2) {

				if (chargeSpecified) {
					atom.partialCharge = charge;
				} else {
					atom.amberName = splitSpec [1];
				}
			}

			if (splitSpec.Length == 3 && chargeSpecified) {
				if (float.TryParse (splitSpec [1], out charge)) {
					atom.partialCharge = charge;
				}
			}
		}

		//Read the info in parentheses
		if (atomSpec.Contains ("(")) {
			if (splitSpec.Length == 0) {
				atom.element = atomSpec.Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries) [0];
			}

			string pdbStr = GetStringInParentheses (atomSpec);
			string[] splitPdb = pdbStr.Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries);

			foreach (string pdbOption in splitPdb) {
				if (pdbOption.StartsWith ("PDBNAME")) {
					atom.pdbName = GetValueFromPair (pdbOption);
				} else if (pdbOption.StartsWith ("RESNAME")) {
					atom.residueName = GetValueFromPair (pdbOption);
				} else if (pdbOption.StartsWith ("RESNUM")) {
					atom.residueNumber = int.Parse (ToNumber( GetValueFromPair (pdbOption)));
				}
			}
		} else if (splitSpec.Length == 0) {
			atom.element = atomSpec;
		}

		splitIndex += 1;
		
		//Check if frozen, which is optional. Gaussian checks if there is an int here (instead of float)
		int frozen = 0;
		if (int.TryParse (splitLine [splitIndex], out frozen)) {
			splitIndex += 1;
		}
		atom.frozen = frozen;

		Vector3 position = new Vector3 ();
		//Positions
		for (int c = 0; c < 3; c++) {
			position [c] = float.Parse(splitLine [splitIndex]);
			splitIndex += 1;
		}
		atom.transform.position = position;

		return atom;
	}

	private void ConnectionsFromGaussianInputLine(string inputStr, Graph graph) {
		string[] splitConn = inputStr.Split(new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

		int index = int.Parse (splitConn [0]) - 1;
		int connections = (splitConn.Length - 1) / 2;

		for (int connectionNum = 0; connectionNum < connections; connectionNum++) {
			int connectedIndex = int.Parse(splitConn [connectionNum * 2 + 1]) - 1;
			float gaussBondOrder = float.Parse(splitConn [connectionNum * 2 + 2]);

			//Read 1.5 as aromatic, and 0.5 for connectivity purposes
			int bondOrder = gaussBondOrder == 1.5f ? 5 : gaussBondOrder == 0.5f ? 6 : (int)gaussBondOrder;

			graph.Connect (index, connectedIndex, bondOrder);
		}


	}

}
