using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class PQRReader : MonoBehaviour {

	public Settings globalSettings;
	public Data data;

	private int lineNumber;
	private string line;
	private GaussianCalculator gc;

	public IEnumerator AtomsFromPQRFile(string path, Atoms atoms, Atom atomPrefab) {
		atoms.globalSettings = globalSettings;
		atoms.setSourceFile(path);

		atoms.name = Path.GetFileName (path);
		gc = Instantiate<GaussianCalculator> (GameObject.FindObjectOfType<PrefabManager>().gaussianCalculatorPrefab, atoms.transform);
		atoms.gaussianCalculator = gc;

		int index = 0;
		int residueNumber;
		int offset;
		Vector3 position;
		string[] splitLine;
		List<int> highAtoms = new List<int>();
		List<int> lowAtoms = new List<int>();

		string[] lines = FileIO.Readlines (path, "#");
		for ( lineNumber=0; lineNumber < lines.Length; lineNumber++) {
			line = lines [lineNumber];
			if ((line.StartsWith ("ATOM") || line.StartsWith ("HETATM") )) {
				Atom atom = Instantiate (atomPrefab, atoms.atomsHolder.transform);

				splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);

				//PDB
				atom.pdbName = line.Substring (12, 4).Trim();

				//Element
				if (data.pdbToElementDict.ContainsKey (atom.pdbName.ToUpper ())) {
					atom.element = data.pdbToElementDict [atom.pdbName];
				} else {
					atom.element = atom.pdbName.Substring (0, 1);
				}

				//Residue
				atom.residueName = splitLine[3];

				if (line.StartsWith ("HETATM") && ! (atom.residueName == "WAT" || atom.residueName == "HOH")) {
					highAtoms.Add (index);
				} else {
					lowAtoms.Add(index);
				}

				//Chain ID is optional in PQR, but can also merge with residue number
				residueNumber = 0;
				offset = 0;
				if (!int.TryParse (splitLine [4], out residueNumber)) {
					if (splitLine [4].Length == 1) {
						offset = 1;
						atom.chainID = splitLine [4];
						residueNumber = int.Parse (splitLine [5]);
					} else {
						atom.chainID = splitLine [4].Substring(0,1);
						residueNumber = int.Parse (splitLine [4].Substring(1));
					}
				}
				atom.residueNumber = residueNumber;

				//Position
				position = new Vector3 ();
				position.x = float.Parse (splitLine[5 + offset]);
				position.y = float.Parse (splitLine[6 + offset]);
				position.z = float.Parse (splitLine[7 + offset]);
				atom.transform.position = position;

				atom.partialCharge = float.Parse (splitLine[8 + offset]);
				atom.vdwRadius = float.Parse (splitLine[9 + offset]);

				atom.index = index;
				atoms.AddAtom (atom);
				//atom.resolution = globalSettings.ballResolution;
				//atom.transform.SetParent(atoms.atomsHolder.transform);
				//atoms.atomList.Add (atom);

				index++;
			}

			if (lineNumber % globalSettings.maxLinesToYield == 0)
				yield return new WaitForSeconds (0.1f);
		}

		atoms.graph.SetAtoms (atoms);
		gc.SetAtoms(atoms);
		gc.layers.SetLowAtoms (lowAtoms);
		gc.layers.SetHighAtoms (highAtoms);

		yield return null;
	}

	public IEnumerator NonStandardResiduesFromPQRFile(string path, List<int> nonStandardResidueList) {

		bool readResidues = false;
		int residueNumber;

		string[] lines = FileIO.Readlines (path);
		string[] splitLine;
		for ( lineNumber=0; lineNumber < lines.Length; lineNumber++) {

			line = lines [lineNumber];
			if (line.StartsWith ("REMARK")) {
				if (line.Contains("(omitted below):")) {
					readResidues = true;
					continue;
				} else if (line.Contains("This is usually due")) {
					break;
				}
			} 
			if (readResidues) {
				splitLine = line.Split(new char[] {' '}, System.StringSplitOptions.RemoveEmptyEntries);
				residueNumber = int.Parse(splitLine[splitLine.Length - 1]);

				if (!nonStandardResidueList.Contains(residueNumber)) {
					nonStandardResidueList.Add(residueNumber);
				}
			}
		}

		yield return null;
	}

}
