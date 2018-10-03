using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Main : MonoBehaviour {

	public Settings globalSettings;
	public FileReader fileReader;
	public Data data;

	// Use this for initialization
	void Start () {
		fileReader.globalSettings = globalSettings;
		fileReader.data = data;
		//Atoms atoms = fileReader.AtomsFromGaussianInput ("/Users/tam10/Notebooks/GFPMutation/protonation/CRA.pdb");
		Atoms atoms = fileReader.AtomsFromGaussianInput ("/Users/tam10/Notebooks/GFPMutation/model.com");
		atoms.transform.parent = transform;
		atoms.RenderAll ();

	}
	
	// Update is called once per frame
	void Update () {
		
	}


	//Parameters
	public void AddParametersFromFile(string filename, Atoms atoms) {
		string[] split = filename.Split(new []{ "." }, System.StringSplitOptions.RemoveEmptyEntries);
		string filetypeA = split[split.Length - 1].ToLower();

		string[] splitB = split [split.Length - 2].Split(new []{ System.IO.Path.DirectorySeparatorChar }, System.StringSplitOptions.RemoveEmptyEntries);
		string filetypeB = splitB [splitB.Length - 1].ToLower ();

		Parameters parameters;

		if (filetypeA == "frcmod" || filetypeB == "frcmod") {
			parameters = fileReader.ParametersFromFRCMODFile (filename);
		} else {
			Debug.LogError("Filetype '" + filetypeA + "' not available for parameters");
			return;
		}

		parameters.transform.parent = atoms.transform;
		atoms.parameters = parameters;
	}
}
