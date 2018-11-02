using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Text;
using System.Diagnostics;
using System.IO;

public class Protonation : MonoBehaviour {

	public Atoms protonatedAtoms;
	public Settings globalSettings;
	public FileReader fileReader;
	public FileWriter fileWriter;

	public List<int> nonStandardResidues;

	public Bash bash;

	public IEnumerator GetProtonatedStandardResidues(Atoms atoms, Atoms protonatedAtoms, List<int> nonStandardResidueList) {

		string pdbName = globalSettings.standardResidueProtonationFilename + ".pdb";
		string pqrName = globalSettings.standardResidueProtonationFilename + ".pqr";
		string logName = globalSettings.standardResidueProtonationFilename + ".log";

		string pdbPath = Path.Combine(globalSettings.tempFolder, pdbName);
		string pqrPath = Path.Combine(globalSettings.tempFolder, pqrName);

		fileWriter.WritePDBFile (atoms, pdbPath, false);

		StringBuilder sb = new StringBuilder ();

		sb.Append (globalSettings.pdb2pqrCommand + " -v ");
		foreach (string option in globalSettings.pdb2pqrOptions) {
			sb.Append (option + " ");
		}

		sb.Append(string.Format("--ph-calc-method=propka --with-ph={0:0.0} {1} {2} > {3}", globalSettings.pH, pdbName, pqrName, logName)); 

		Process proc = bash.GetBashProcess (sb.ToString ());
		proc.Start ();
		proc.WaitForExit ();

		string error = proc.StandardError.ReadToEnd ();

		if (error != "") {
			UnityEngine.Debug.LogErrorFormat ("pdb2pqr failed with error message: {0}", error);
		} else {
			yield return StartCoroutine(fileReader.LoadAtomsFileAsync(pqrPath, protonatedAtoms));
			yield return StartCoroutine(fileReader.NonStandardResiduesFromPQRFile(pqrPath, nonStandardResidueList));
		}

		yield return null;
	}


}
