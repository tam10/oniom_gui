using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class GaussianInputReader : MonoBehaviour {

	private GaussianCalculator gc;
	public PRMReader prmReader;

	public Settings globalSettings;
	public Data data;

	enum atomPhase {
		ELEMENT,
		AMBER_NAME,
		PARTIAL_CHARGE,
		PDB,
		PDB_VALUE,
		FROZEN,
		X,
		Y,
		Z,
		LAYER_NAME,
		LINK,
		LINK_INDEX
	};

	enum pdbOptions {
		PDBNAME,
		RESNAME,
		RESNUM
	}

	public IEnumerator AtomsFromGaussianInput(string path, Atoms atoms, Atom atomPrefab) {

		atoms.globalSettings = globalSettings;
		atoms.setSourceFile(path);

		atoms.name = Path.GetFileName (path);
		gc = Instantiate<GaussianCalculator> (GameObject.FindObjectOfType<PrefabManager>().gaussianCalculatorPrefab, atoms.transform);
		atoms.gaussianCalculator = gc;


		int index = 0;
		int lineNumber = 0;
		string line;

		string[] lines = FileIO.Readlines (path, "!");
		string[] stringArray;
		int stringArrayIndex;
		int maxLines = lines.Length;

		//Link 0 commands
		while (lineNumber < maxLines) {
			line = lines [lineNumber].Trim();
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
				stringArray = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
				gc.killJobLink = stringArray [1];

				if (stringArray.Length > 2) {
					gc.killJobAfter = int.Parse (stringArray [2]);
				} else {
					gc.killJobAfter = 1;
				}
			}

			lineNumber++;
		}


		List<string> keywordsList = new List<string> ();
		string keywordItem;


		//Keywords
		while (lineNumber < maxLines) {
			line = lines [lineNumber].Trim();
			if (line == "") {
				break;
			} 

			stringArray = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
			for (stringArrayIndex = 0; stringArrayIndex < stringArray.Length; stringArrayIndex++) {
				keywordsList.Add (stringArray [stringArrayIndex].ToLower ());
			}
			lineNumber++;
		}

		//Parse keywords - only one line so don't need to worry about optimisation
		for (stringArrayIndex = 0; stringArrayIndex < keywordsList.Count; stringArrayIndex++) {
			keywordItem = keywordsList [stringArrayIndex];
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
				gc.layers.AddLayer (highMBO [0], highMBO [1], highMBO [2], 'H');

				if (methods.Length == 2) {
					string[] lowMBO = GetMethodFromString (methods [1]);
					gc.layers.AddLayer (lowMBO [0], lowMBO [1], lowMBO [2], 'L');
				} else if (methods.Length == 3) {
					string[] mediumMBO = GetMethodFromString (methods [1]);
					gc.layers.AddLayer (mediumMBO [0], mediumMBO [1], mediumMBO [2], 'M');

					string[] lowMBO = GetMethodFromString (methods [2]);
					gc.layers.AddLayer (lowMBO [0], lowMBO [1], lowMBO [2], 'L');
				}
			} else if (keywordItem.StartsWith ("guess")) {
				string guessOptionsStr = GetStringInParentheses (GetValueFromPair (keywordItem, checkEnclosed: true));
				string[] guessOptions = guessOptionsStr.Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries);
				foreach (string guessOption in guessOptions)
					gc.guessOptions.Add (guessOption);
			} else if (keywordItem.StartsWith ("geom")) {
				string geomOptionsStr = GetStringInParentheses (GetValueFromPair (keywordItem, checkEnclosed: true));
				string[] geomOptions = geomOptionsStr.Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries);
				foreach (string geomOption in geomOptions)
					gc.geomOptions.Add (geomOption);
			} else {
				foreach (string methodName in data.gaussianMethods) {
					if (keywordItem.StartsWith (methodName)) {
						string[] highMBO = GetMethodFromString (keywordItem);
						gc.layers.AddLayer (highMBO [0], highMBO [1], highMBO [2], 'H');
						break;
					}
				}
			}
		}

		//Title
		List<string> titleLines = new List<string>();

		if (!gc.geomOptions.Contains ("allcheck")) {
			while (lineNumber < maxLines) {
				line = lines [lineNumber].Trim ();
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
				line = lines [lineNumber].Trim ();
				if (line != "") {
					string[] cmStr = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

					int numLayersSpecified = cmStr.Length / 2;

					for (int layerIndex = 0; layerIndex < numLayersSpecified; layerIndex++) {

						gc.layers.SetLayerCharge (gc.layers.layerNames[layerIndex], int.Parse (cmStr [2 * layerIndex]));
						gc.layers.SetLayerMultiplicity (gc.layers.layerNames[layerIndex], int.Parse (cmStr [2 * layerIndex + 1]));
					}

					lineNumber++;
					break;
				} 
				lineNumber++;
			}
		}

		System.Text.StringBuilder sb = new System.Text.StringBuilder ();
		int charNum = 0;
		int phase = 0;
		char lineChar;
		int pdbOption = 0;
		string linkType = "";
		int linkIndex = 0;
		Vector3 position = new Vector3();
		List<int> highAtoms = new List<int>();
		List<int> mediumAtoms = new List<int>();
		List<int> lowAtoms = new List<int>();

		//Atoms - AtomFromGaussianInputLine should be optimised
		if (!(gc.geomOptions.Contains ("allcheck") && gc.geomOptions.Contains ("check") && gc.geomOptions.Contains ("checkpoint"))) {

			while (lineNumber < maxLines) {
				line = lines [lineNumber].Trim ();
				if (line == "") {
					lineNumber++;
					break;
				}

				Atom atom = Instantiate<Atom> (atomPrefab);
				atom.index = index;
				charNum = 0;
				phase = (int)atomPhase.ELEMENT;
				sb.Length = 0;

				while (charNum <= line.Length) {

					if (charNum == line.Length) {
						lineChar = ' ';
					} else {
						lineChar = line [charNum];
					}


					//Read element
					if (phase == (int)atomPhase.ELEMENT) {
						if (lineChar == '-') {
							phase = (int)atomPhase.AMBER_NAME;
						} else if (lineChar == '(') {
							phase = (int)atomPhase.PDB; 
						} else if (lineChar == ' ') {
							phase = (int)atomPhase.FROZEN;
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						atom.element = sb.ToString ();

						sb.Length = 0;

					} else if (phase == (int)atomPhase.AMBER_NAME) {
						if (lineChar == '-') {
							phase = (int)atomPhase.PARTIAL_CHARGE;
						} else if (lineChar == '(') {
							phase = (int)atomPhase.PDB; 
						} else if (lineChar == ' ') {
							phase = (int)atomPhase.FROZEN;
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						atom.amberName = sb.ToString ();
						sb.Length = 0;

					} else if (phase == (int)atomPhase.PARTIAL_CHARGE) {
						if (lineChar == '(') {
							phase = (int)atomPhase.PDB; 
						} else if (lineChar == ' ') {
							phase = (int)atomPhase.FROZEN;
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						if (sb.Length > 0) 
							atom.partialCharge = float.Parse(sb.ToString ());
						
						sb.Length = 0;

					} else if (phase == (int)atomPhase.PDB) {
						if (lineChar == '=') {
							phase = (int)atomPhase.PDB_VALUE;
						} else if (lineChar == ')') {
							phase = (int)atomPhase.FROZEN;
						} else if (sb.Length == 5) {
							if (sb.ToString ().ToUpper () == "PDBNA")
								pdbOption = (int)pdbOptions.PDBNAME;
							if (sb.ToString ().ToUpper () == "RESNA")
								pdbOption = (int)pdbOptions.RESNAME;
							if (sb.ToString ().ToUpper () == "RESNU")
								pdbOption = (int)pdbOptions.RESNUM;
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						} 

						sb.Length = 0;

					} else if (phase == (int)atomPhase.PDB_VALUE) {
						
						if (lineChar == ',') {
							if (pdbOption == (int)pdbOptions.PDBNAME) {
								atom.pdbName = sb.ToString ();
							} else if (pdbOption == (int)pdbOptions.RESNAME) {
								atom.residueName = sb.ToString ();
							} else if (pdbOption == (int)pdbOptions.RESNUM) {
								atom.residueNumber = int.Parse(ToNumber( sb.ToString ()));
							}
							phase = (int)atomPhase.PDB;
						} else if (lineChar == ')') {
							if (pdbOption == (int)pdbOptions.PDBNAME) {
								atom.pdbName = sb.ToString ();
							} else if (pdbOption == (int)pdbOptions.RESNAME) {
								atom.residueName = sb.ToString ();
							} else if (pdbOption == (int)pdbOptions.RESNUM) {
								atom.residueNumber = int.Parse(ToNumber( sb.ToString ()));
							}
							phase = (int)atomPhase.FROZEN;
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						} 

						sb.Length = 0;
					
					//Frozen is optional, so this could be the X coordinate
					} else if (phase == (int)atomPhase.FROZEN) {
						if (lineChar == ' ') {
							if (sb.Length == 0) {
								;
							} else if (sb.Length == 1) {
								atom.frozen = int.Parse (sb.ToString ());
								phase = (int)atomPhase.X;
							} else {
								position.x = float.Parse (sb.ToString ());
								phase = (int)atomPhase.Y;
							}
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						sb.Length = 0;
						charNum++;
						continue;

					} else if (phase == (int)atomPhase.X) {
						if (lineChar == ' ') {
							if (sb.Length == 0) {
								;
							} else {
								position.x = float.Parse (sb.ToString ());
								phase = (int)atomPhase.Y;
							}
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						sb.Length = 0;
						charNum++;
						continue;

					} else if (phase == (int)atomPhase.Y) {
						if (lineChar == ' ') {
							if (sb.Length == 0) {
								;
							} else {
								position.y = float.Parse (sb.ToString ());
								phase = (int)atomPhase.Z;
							}
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						sb.Length = 0;
						charNum++;
						continue;

					} else if (phase == (int)atomPhase.Z) {
						if (lineChar == ' ') {
							if (sb.Length == 0) {
								;
							} else {
								position.z = float.Parse (sb.ToString ());
								atom.transform.position = position;
								phase = (int)atomPhase.LAYER_NAME;
							}
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						sb.Length = 0;
						charNum++;
						continue;

					} else if (phase == (int)atomPhase.LAYER_NAME) {
						if (lineChar == ' ') {
							if (sb.Length == 0) {
								;
							} else {
								if (sb[0] == 'H') {
									highAtoms.Add (index);
								} else if (sb[0] == 'M') {
									mediumAtoms.Add (index);
								} else if (sb[0] == 'L') {
									lowAtoms.Add (index);
								}

								phase = (int)atomPhase.LINK;
							}
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						sb.Length = 0;
						charNum++;
						continue;

					} else if (phase == (int)atomPhase.LINK) {
						if (lineChar == ' ') {
							if (sb.Length == 0) {
								;
							} else {
								linkType = sb.ToString();
								phase = (int)atomPhase.LINK_INDEX;
							}
						} else if (lineChar == '-') {
							sb.Length = 0;
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						sb.Length = 0;
						charNum++;
						continue;

					}  else if (phase == (int)atomPhase.LINK_INDEX) {
						if (lineChar == ' ') {
							if (sb.Length == 0) {
								;
							} else {
								linkIndex = int.Parse(sb.ToString());
								gc.layers.AddLink (index, linkIndex, linkType);
								charNum++;
								continue;
							}
						} else {
							sb.Append (lineChar);
							charNum++;
							continue;
						}

						sb.Length = 0;
						charNum++;
						continue;

					}


					charNum++;
				}
				atoms.AddAtom (atom);

				index++;
				lineNumber++;

				if (lineNumber % globalSettings.maxLinesToYield == 0)
					yield return null;
			}
		}

		gc.SetAtoms (atoms);
		gc.layers.SetHighAtoms (highAtoms);
		gc.layers.SetMediumAtoms (mediumAtoms);
		gc.layers.SetLowAtoms (lowAtoms);


		//Connectivity
		atoms.graph.SetAtoms (atoms);

		if (gc.geomOptions.Contains ("connectivity")) {
			
			string[] splitConn;

			int connectionNum;
			int numConnections;

			int connectionIndex0;
			int connectionIndex1;
			double bondOrder;

			while (lineNumber < maxLines) {
				
				line = lines [lineNumber].Trim ();
				if (line == "")
					break;

				splitConn = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				connectionIndex0 = int.Parse (splitConn [0]) - 1;
				numConnections = (splitConn.Length - 1) / 2;

				for (connectionNum = 0; connectionNum < numConnections; connectionNum++) {
					connectionIndex1 = int.Parse (splitConn [connectionNum * 2 + 1]) - 1;
					bondOrder = double.Parse (splitConn [connectionNum * 2 + 2]);

					atoms.graph.Connect (connectionIndex0, connectionIndex1, bondOrder);
				}

				lineNumber++;

				if (lineNumber % globalSettings.maxLinesToYield == 0)
					yield return new WaitForSeconds (0.1f);
			}
		}

		//Parameters
		List<string> parameterLines = new List<string>();

		while (lineNumber < maxLines) {
			line = lines [lineNumber].Trim ();
			parameterLines.Add (line);
			lineNumber++;
		}

		IEnumerator iEnumerator = prmReader.ParametersFromPRMLines (parameterLines.ToArray (), atoms.graph.parameters);
		while (iEnumerator.MoveNext ()) {}
		
		atoms.graph.parameters.transform.parent = atoms.transform;

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

	private bool StringIndexEnclosed(string inputStr, int index, char openChar='(', char closeChar=')') {

		int depth = 0;
		char chr;

		for (int i = 0; i < inputStr.Length; i++) {
			chr = inputStr [i];

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


	private string ToAlpha(string inputStr) {
		return System.Text.RegularExpressions.Regex.Replace (inputStr, @"[^a-zA-Z -]", string.Empty);
	}

	private string ToNumber(string inputStr) {
		return System.Text.RegularExpressions.Regex.Replace (inputStr, @"[^0-9 -]", string.Empty);
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
}
