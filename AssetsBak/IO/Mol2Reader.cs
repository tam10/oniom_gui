using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class Mol2Reader : MonoBehaviour {

	private GaussianCalculator gc;
	public GaussianCalculator gaussianCalculatorPrefab;

	public Settings globalSettings;
	public Data data;

	public IEnumerator AtomsFromMol2File(string path, Atoms atoms, Atom atomPrefab) {
		
		atoms.globalSettings = globalSettings;
		atoms.setSourceFile(path);

		atoms.name = Path.GetFileName (path);
		gc = Instantiate<GaussianCalculator> (gaussianCalculatorPrefab, atoms.transform);
		atoms.gaussianCalculator = gc;

		int index = 0;

		Vector3 position = new Vector3 ();

		List<Data.AmberToElement> amberToElementList = data.amberToElementData;
		Data.AmberToElement amberToElement;
		int amberToElementNum;

		List<int[]> connectionList = new List<int[]> ();
		int a0;
		int a1;
		int bondOrder;

		bool readAtoms = false;
		bool readConnectivity = false;
		string[] lines = FileIO.Readlines (path);
		string[] splitLine;

		int lineNumber;
		string line;

		for (lineNumber = 0; lineNumber < lines.Length; lineNumber++) {
			line = lines [lineNumber];
			if (line.StartsWith ("@<TRIPOS>")) {
				readAtoms = line.Contains ("ATOM");
				//readConnectivity = line.Contains("BOND");
			} else if (readAtoms) {
				splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
				Atom atom = Instantiate<Atom> (atomPrefab);

				//Atom types
				atom.pdbName = splitLine [1];

				atom.amberName = splitLine [5];

				for (amberToElementNum = 0; amberToElementNum < amberToElementList.Count; amberToElementNum++) {
					amberToElement = amberToElementList [amberToElementNum];
					if (amberToElement.amberName == atom.amberName) {
						atom.element = amberToElement.element;
						atom.formalCharge = amberToElement.formalCharge;
					}
				}


				//Position
				position.x = float.Parse (splitLine [2]);
				position.y = float.Parse (splitLine [3]);
				position.z = float.Parse (splitLine [4]);
				atom.transform.position = position;

				//Residue
				atom.residueNumber = int.Parse (splitLine [6]);
				atom.residueName = splitLine [7];

				if (splitLine.Length > 8)
					atom.partialCharge = float.Parse (splitLine [8]);

				atom.index = index;
				atoms.AddAtom (atom);
				//atom.resolution = globalSettings.ballResolution;
				//atom.transform.SetParent(atoms.atomsHolder.transform);
				//atoms.atomList.Add (atom);

				index++;

			} else if (readConnectivity) {
				splitLine = line.Split (new []{ " " }, System.StringSplitOptions.RemoveEmptyEntries);
				a0 = int.Parse (splitLine [1]) - 1;
				a1 = int.Parse (splitLine [2]) - 1;
				bondOrder = int.Parse (splitLine [3]);
				connectionList.Add (new int[3] { a0, a1, bondOrder });
			}

			if (lineNumber % globalSettings.maxLinesToYield == 0)
				yield return new WaitForSeconds (0.1f);
		}

		atoms.graph.SetAtoms (atoms);
		gc.layers.SetSize (atoms.size);
		foreach (int[] connectionPair in connectionList) {
			atoms.graph.Connect (connectionPair [0], connectionPair [1], (double)connectionPair [2]);
		}
			
	}
}
