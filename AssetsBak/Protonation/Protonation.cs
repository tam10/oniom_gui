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

	public Bash bash;

	public Atoms GetProtonatedStandardResidues(Atoms atoms) {

		string tempFolder = Path.Combine(UnityEngine.Application.persistentDataPath, "Temp");

		if (!Directory.Exists(tempFolder)) {
			Directory.CreateDirectory(tempFolder);
		}

		string pdbPath = Path.Combine(tempFolder, "protTest.pdb");
		string pqrPath = Path.Combine(tempFolder, "protTest.pqr");
		string logPath = Path.Combine(tempFolder, "protTest.log");

		fileWriter.WritePDBFile (atoms, pdbPath, false);

		StringBuilder sb = new StringBuilder ();

		sb.Append (string.Format("\"{0}\" -v ", globalSettings.pdb2pqrCommand));
		foreach (string option in globalSettings.pdb2pqrOptions) {
			sb.Append (option + " ");
		}

		sb.Append(string.Format("--ph-calc-method=propka --with-ph={0:0.0} \"{1}\" \"{2}\" > \"{3}\"", globalSettings.pH, pdbPath, pqrPath, logPath)); 

		Process proc = bash.GetBashProcess (sb.ToString ());
		proc.Start ();
		proc.WaitForExit ();

		string error = proc.StandardError.ReadToEnd ();

		if (error != "") {
			UnityEngine.Debug.LogErrorFormat ("pdb2pqr failed with error message: {0}", error);
			return null;
		}

		Atoms protonatedAtoms = fileReader.LoadAtomsFile(pqrPath);
		return protonatedAtoms;

	}


}
