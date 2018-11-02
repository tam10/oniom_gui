using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class FileReader: MonoBehaviour {

	//Prefabs
	public Atom atomPrefab;
	public Atoms atomsPrefab;

	public Settings globalSettings;
	public Data data;

	public GaussianInputReader gaussianInputReader;
	public PDBReader pdbReader;
	public PQRReader pqrReader;
	public Mol2Reader mol2Reader;

	public PRMReader prmReader;
	public FRCMODReader frcmodReader;

	public string alternateLocationID="A";

	private int _busy;
	public bool busy {
		get { return (_busy > 0); }
	}

	void Awake() {
		atomsPrefab = GameObject.FindObjectOfType<PrefabManager>().atomsPrefab;
		atomPrefab = GameObject.FindObjectOfType<PrefabManager>().atomPrefab;
	}

	public void SetGlobalSettings(Settings globalSettings) {
		this.globalSettings = globalSettings;
		gaussianInputReader.globalSettings = globalSettings;
		mol2Reader.globalSettings = globalSettings;
		pdbReader.globalSettings = globalSettings;
		pqrReader.globalSettings = globalSettings;
	}

	public void SetData(Data data) {
		this.data = data;
		gaussianInputReader.data = data;
		mol2Reader.data = data;
		pdbReader.data = data;
		pqrReader.data = data;

	}

	public Atoms LoadAtomsFile(string path) {
		_busy++;

		Atoms atoms = null;

		string filetype = Path.GetExtension (path).ToLower();
		switch (filetype) {
		case ".com":
		case ".gjf":
			atoms = AtomsFromGaussianInput (path);
			break;
		case ".mol2":
			atoms = AtomsFromMol2File (path);
			break;
		case ".pdb":
			atoms = AtomsFromPDBFile (path);
			break;
		case ".pqr":
			atoms = AtomsFromPQRFile (path);
			break;
		default:
			Debug.LogErrorFormat ("Filetype {0} not recognised", filetype);
			break;
		}

		_busy--;
		return atoms;
	}

	public IEnumerator LoadAtomsFileAsync(string path, Atoms atoms) {
		_busy++;
		string filetype = Path.GetExtension (path).ToLower();

		Coroutine coroutine = null;

		switch (filetype) {
		case ".com":
		case ".gjf":
			coroutine = StartCoroutine(AtomsFromGaussianInputCoroutine (path, atoms));
			break;
		case ".mol2":
			coroutine = StartCoroutine(AtomsFromMol2FileCoroutine (path, atoms));
			break;
		case ".pdb":
			coroutine = StartCoroutine(AtomsFromPDBFileCoroutine (path, atoms));
			break;
		case ".pqr":
			coroutine = StartCoroutine(AtomsFromPQRFileCoroutine (path, atoms));
			break;
		default:
			Debug.LogErrorFormat ("Filetype {0} not recognised", filetype);
			break;
		}
		
		_busy--;
		yield return coroutine;
	}

	public IEnumerator LoadParametersFromNameAsync(string paramsName, Parameters parameters) {
		_busy++;
		Coroutine coroutine = null;

		string path = globalSettings.GetPath(globalSettings.baseDataDirectories, Path.Combine("Amber", paramsName));

		if (path == string.Empty) {
			throw new System.Exception(string.Format("Path for {0} could not be found", paramsName));
		}

		if (paramsName.ToLower().StartsWith("frcmod")) {
			coroutine = StartCoroutine(ParametersFromFRCMODFile(path, parameters));
		} else if (paramsName.ToLower().EndsWith("prm")) {
			coroutine = StartCoroutine(ParametersFromPRMFile(path, parameters));
		}

		_busy--;
		yield return coroutine;
	}

	public Atoms AtomsFromGaussianInput(string path) {
		_busy++;

		Atoms atoms = Instantiate<Atoms> (atomsPrefab);
		IEnumerator iEnumerator = gaussianInputReader.AtomsFromGaussianInput (path, atoms, atomPrefab);
		while (iEnumerator.MoveNext ()) {}

		_busy--;
		return atoms;
	}
	public IEnumerator AtomsFromGaussianInputCoroutine(string path, Atoms atoms) {
		_busy++;

		Coroutine coroutine = StartCoroutine(gaussianInputReader.AtomsFromGaussianInput (path, atoms, atomPrefab));
		yield return coroutine;

		_busy--;
		yield return null;

	}

	public Atoms AtomsFromMol2File(string path) {
		_busy++;

		Atoms atoms = Instantiate<Atoms> (atomsPrefab);
		IEnumerator iEnumerator = mol2Reader.AtomsFromMol2File (path, atoms, atomPrefab);
		while (iEnumerator.MoveNext ()) {}

		_busy--;
		return atoms;
	}
	public IEnumerator AtomsFromMol2FileCoroutine(string path, Atoms atoms) {
		_busy++;

		Coroutine coroutine = StartCoroutine(mol2Reader.AtomsFromMol2File (path, atoms, atomPrefab));
		yield return coroutine;

		_busy--;
		yield return null;

	}

	public Atoms AtomsFromPDBFile(string path) {
		_busy++;

		Atoms atoms = Instantiate<Atoms> (atomsPrefab);
		IEnumerator iEnumerator = pdbReader.AtomsFromPDBFile (path, atoms, atomPrefab, alternateLocationID);
		while (iEnumerator.MoveNext ()) {}

		_busy--;
		return atoms;
	}
	public IEnumerator AtomsFromPDBFileCoroutine(string path, Atoms atoms) {
		_busy++;

		Coroutine coroutine = StartCoroutine(pdbReader.AtomsFromPDBFile (path, atoms, atomPrefab, alternateLocationID));
		yield return coroutine;

		_busy--;
		yield return null;
	}

	public Atoms AtomsFromPQRFile(string path) {
		_busy++;

		Atoms atoms = Instantiate<Atoms> (atomsPrefab);
		IEnumerator iEnumerator = pqrReader.AtomsFromPQRFile (path, atoms, atomPrefab);
		while (iEnumerator.MoveNext ()) {}

		_busy--;
		return atoms;
	}
	public IEnumerator AtomsFromPQRFileCoroutine(string path, Atoms atoms) {
		_busy++;
		Coroutine coroutine = StartCoroutine(pqrReader.AtomsFromPQRFile (path, atoms, atomPrefab));
		yield return coroutine;
		_busy--;
		yield return null;
	}

	public IEnumerator NonStandardResiduesFromPQRFile(string path, List<int> nonStandardResidueList) {
		_busy++;
		Coroutine coroutine = StartCoroutine(pqrReader.NonStandardResiduesFromPQRFile (path, nonStandardResidueList));
		yield return coroutine;
		_busy--;
		yield return null;
	}

	public IEnumerator ParametersFromFRCMODFile(string filename, Parameters parameters) {
		_busy++;
		Coroutine coroutine = StartCoroutine(frcmodReader.ParametersFromFRCMODFile (filename, parameters));
		_busy--;
		yield return coroutine;
	}

	public IEnumerator ParametersFromPRMFile(string filename, Parameters parameters) {
		_busy++;
		Coroutine coroutine = StartCoroutine(prmReader.ParametersFromPRMFile (filename, parameters));
		_busy--;
		yield return coroutine;
	}
		
}
