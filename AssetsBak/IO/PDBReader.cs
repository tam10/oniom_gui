using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class PDBReader : MonoBehaviour {

	public GaussianCalculator gaussianCalculatorPrefab;
	public Settings globalSettings;
	public Data data;

	private int lineNumber;
	private string line;
	private GaussianCalculator gc;

	public IEnumerator AtomsFromPDBFile(string path, Atoms atoms, Atom atomPrefab, string alternateLocationID="A") {

		atoms.globalSettings = globalSettings;
		atoms.setSourceFile(path);

		atoms.name = Path.GetFileName (path);
		gc = Instantiate<GaussianCalculator> (gaussianCalculatorPrefab, atoms.transform);
		atoms.gaussianCalculator = gc;

		string amber;
		string element;
		Vector3 position;
		string chargeStr;
		List<int> highAtoms = new List<int>();
		List<int> lowAtoms = new List<int>();

		List<Data.AmberToElement> amberToElementList = data.amberToElementData;
		Data.AmberToElement amberToElement;
		int amberToElementNum;

		int index = 0;
		bool amberUsed;

		string[] lines = FileIO.Readlines (path, "#");
		for ( lineNumber=0; lineNumber < lines.Length; lineNumber++) {
			line = lines [lineNumber];
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
				amberUsed = false;
				amber = line.Substring (76, 2).Trim ();
				if (amber != string.Empty) {
					
					atom.amberName = amber;

					for (amberToElementNum = 0; amberToElementNum < amberToElementList.Count; amberToElementNum++) {
						amberToElement = amberToElementList [amberToElementNum];
						if (amberToElement.amberName == amber) {
							atom.element = amberToElement.element;
							atom.formalCharge = amberToElement.formalCharge;
							amberUsed = true;
						}
					}
				}

				//PDB
				atom.pdbName = line.Substring (12, 4).Trim();

				//Use PDB to get Atom
				if (!amberUsed) {
					//If 12th column is a letter, it is a metal
					if (char.IsLetter (line [12])) {
						element = line.Substring (12, 2).Trim ();
						atom.element = ToAlpha(element);
					} else {
						element = line.Substring (12, 4).Trim ();
						atom.element = ToAlpha(element) [0].ToString ();
					}
				}


				//Residue
				atom.residueName = line.Substring (17, 3).Trim ();
				atom.residueNumber = int.Parse (line.Substring (22, 4).Trim ());
				atom.chainID = line [21].ToString ();

				if (line.StartsWith ("HETATM") && ! (atom.residueName == "WAT" || atom.residueName == "HOH")) {
					highAtoms.Add (index);
				} else {
					lowAtoms.Add(index);
				}

				//Position
				position = new Vector3 ();
				position.x = float.Parse (line.Substring (30, 8));
				position.y = float.Parse (line.Substring (38, 8));
				position.z = float.Parse (line.Substring (46, 8));
				atom.transform.position = position;

				chargeStr = line [79].ToString ();
				if (chargeStr != " ") {
					if (line [78].ToString () == "-") {
						atom.formalCharge = -int.Parse (chargeStr);
					} else {
						atom.formalCharge = int.Parse (chargeStr);
					}
				}

				atom.index = index;
				atoms.AddAtom (atom);

				index++;
			}

			if (lineNumber % globalSettings.maxLinesToYield == 0)
				yield return new WaitForSeconds (0.1f);
		}

		atoms.graph.SetAtoms (atoms);
		gc.layers.SetSize (atoms.size);
		gc.layers.SetLowAtoms (lowAtoms);
		gc.layers.SetHighAtoms (highAtoms);
		yield return null;
	}

	private string ToAlpha(string inputStr) {
		return System.Text.RegularExpressions.Regex.Replace (inputStr, @"[^a-zA-Z -]", string.Empty);
	}

	private string ToNumber(string inputStr) {
		return System.Text.RegularExpressions.Regex.Replace (inputStr, @"[^0-9 -]", string.Empty);
	}
}
