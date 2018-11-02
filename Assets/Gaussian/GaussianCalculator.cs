using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Text;

public class GaussianCalculator : MonoBehaviour {

	private Atoms parent;

	public int nProcessors;
	public int jobMemoryMB;

	public string checkpointPath;
	public string oldCheckpointPath;

	public string killJobLink;
	public int killJobAfter;

	//LAYERS
	public Layers layers;
	
	//KEYWORDS
	//Character after # in keywords
	private string _gaussianPrintLevel;
	public string gaussianPrintLevel {
		set {
			List<string> allowed = new List<string> { "P", "T", "N", "" };
			if (allowed.Contains (value.ToUpper ())) {
				this._gaussianPrintLevel = value.ToUpper ();
			} else {
				this._gaussianPrintLevel = "";
			}
		}

		get {
			return this._gaussianPrintLevel;
		}
	}

	//Methods, basis, options
	public List<string> oniomOptions;

	//Guess
	public List<string> guessOptions;

	//Geometry
	public List<string> geomOptions;

	//Frequency
	public bool doFreq;
	public List<string> freqOptions;

	public List<string> additionalKeywords;

	//Title
	public string title;

	// Use this for initialization
	void Awake () {
		
		nProcessors = 1;
		jobMemoryMB = 4000;

		checkpointPath = "";
		oldCheckpointPath = "";

		killJobLink = "";
		killJobAfter = 1;

		gaussianPrintLevel = "";

		doFreq = false;

		layers = new Layers ();

		oniomOptions = new List<string> ();
		guessOptions = new List<string> ();
		geomOptions = new List<string> ();
		freqOptions = new List<string> ();
		additionalKeywords = new List<string> ();

		title = "Title";

		parent = GetComponentInParent<Atoms>();

		if (parent == null) {
			Debug.Log("Gaussian Calculator parent is not Atoms object");
		}

	}

	public void SetAtoms(Atoms atoms) {
		parent = atoms;
		layers.SetAtoms(atoms);
	}
}

public class Layers {
	public Dictionary<int, Dictionary<int, string>> linkAtomAmberType = new Dictionary<int, Dictionary<int, string>>();
	private int[] _layerList;

	private Atoms atoms;

	private int[] layerList {

		get {
			if (_layerList == null) {
				return new int[] {};
			}
			return _layerList;
		}
		set {_layerList = value;}
	}

	//private int layerList[int index] {
	//	get {return layerList[index];}
	//}

	private Dictionary<char, Layer> layerDict = new Dictionary<char, Layer>();
	public Dictionary<char, int> nameToNum = new Dictionary<char, int> { { 'H', 0 }, { 'M', 1 }, { 'L', 2 } };
	public Dictionary<int, char> numToName = new Dictionary<int, char> { { 0, 'H' }, { 1, 'M' }, { 2, 'L' } };
	public List<char> layerNames = new List<char> {};

	public void SetAtoms(Atoms atoms) {
		this.atoms = atoms;
		int[] tempArray = new int[atoms.size];
		for (int atomNum =0; atomNum<Mathf.Min( layerList.Length, atoms.size);atomNum++) {
			tempArray[atomNum] = layerList[atomNum];
		}
		layerList = tempArray;
	}

	public void SetHighAtoms(List<int> highAtomNums) {
		foreach (int atomNum in highAtomNums)
			AddAtomToLayerByNumber(0, atomNum);
	}

	public void SetMediumAtoms(List<int> mediumAtomNums) {
		foreach (int atomNum in mediumAtomNums)
			AddAtomToLayerByNumber(1, atomNum);
	}

	public void SetLowAtoms( List<int> lowAtomNums) {
		foreach (int atomNum in lowAtomNums)
			AddAtomToLayerByNumber(2, atomNum);
	}

	public void AddLayer(string method="", string basis="", string options="", char layerName='H', int charge=0, int multiplicity=0 ) {
		Layer layer = new Layer (method, basis, options, layerName, charge, multiplicity);
		if (layerDict.ContainsKey (layerName)) {
			Debug.LogErrorFormat ("Layer {0} can't be added: layer already exists", layerName);
			return;
		}

		layerDict.Add (layerName, layer);
		if (!layerNames.Contains(layerName))
			layerNames.Add (layerName);
	}

	public void AddLink(int linkAtomConnection, int linkAtomHost, string linkAtomAmber) {
		if (!linkAtomAmberType.ContainsKey (linkAtomConnection)) {
			linkAtomAmberType.Add (linkAtomConnection, new Dictionary<int, string> ());
		}

		if (!linkAtomAmberType [linkAtomConnection].ContainsKey (linkAtomHost)) {
			linkAtomAmberType [linkAtomConnection].Add (linkAtomHost, linkAtomAmber);
		} 

	}

	public void AddAtomToLayerByName(char layerName, int atomNum) {
		AddAtomToLayerByNumber( nameToNum [layerName], atomNum);
	}

	public void AddAtomToLayerByNumber(int layerNum, int atomNum) {
		layerList [atomNum] = layerNum;
		atoms[atomNum].showSphere = atoms.globalSettings.layerSphereRenderDict[layerNum];
	}

	public IEnumerator AddAtomsToLayerByNumberAsync(int layerNum, List<int> atomNums) {
		foreach (int atomNum in atomNums) {
			AddAtomToLayerByNumber(layerNum, atomNum);
		}
		yield return null;
	}

	public IEnumerator AddAtomsToLayerByNameAsync(char layerName, List<int> atomNums) {
		foreach (int atomNum in atomNums) {
			AddAtomToLayerByName(layerName, atomNum);
		}
		yield return null;
	}

	public int GetAtomLayerFromAtomNum(int atomNum) {
		return layerList [atomNum];
	}

	public List<int> GetLayerAtomsByInt(int layerNum) {
		List<int> layerAtoms = new List<int> ();
		for (int atomNum = 0; atomNum < layerList.Length; atomNum++) {
			if (layerList [atomNum] == layerNum)
				layerAtoms.Add (atomNum);
		}
		return layerAtoms;
	}

	public List<int> GetLayerAtomsByName(char layerNum) {
		return GetLayerAtomsByInt( nameToNum[layerNum]);
	}

	public int GetLayerCharge(char layerName, bool cascade=true) {
		return GetLayerByName (layerName, cascade).charge;
	}

	public int GetLayerMultiplicity(char layerName, bool cascade=true) {
		return GetLayerByName (layerName, cascade).multiplicity;
	}

	public string GetLayerMethod(char layerName, bool cascade=true) {
		return GetLayerByName (layerName, cascade).method;
	}

	public string GetLayerBasis(char layerName, bool cascade=true) {
		return GetLayerByName (layerName, cascade).basis;
	}

	public string GetLayerOptions(char layerName, bool cascade=true) {
		return GetLayerByName (layerName, cascade).options;
	}

	public void SetLayerCharge(char layerName, int charge) {
		GetLayerByName (layerName).charge = charge;
	}

	public void SetLayerMultiplicity(char layerName, int multiplicity) {
		GetLayerByName (layerName).multiplicity = multiplicity;
	}

	public void SetLayerMethod(char layerName, string method) {
		GetLayerByName (layerName).method = method;
	}

	public void SetLayerBasis(char layerName, string basis) {
		GetLayerByName (layerName).basis = basis;
	}

	public void SetLayerOptions(char layerName, string options) {
		GetLayerByName (layerName).options = options;
	}

	//Hide structures behind methods
	private Layer GetLayerByName(char layerName, bool cascade=true) {
		if (!cascade)
			return layerDict [layerName];
		return layerDict [Cascade(layerName)];

	}

	//Use this to get the next layer up if this layerName doesn't exist
	private char Cascade(char layerName) {
		int layerNum = nameToNum [layerName];

		while (layerNum > 0) {
			layerName = numToName [layerNum];
			if (layerDict.ContainsKey (layerName))
				break;
			
			layerNum--;
		}
		return layerName;
	}

	public override string ToString ()
	{
		StringBuilder sb = new StringBuilder ();

		int layerNum = 0;
		while (true) {
			if (layerNum == layerNames.Count)
				break;
			sb.Append (string.Format("{0}: {1}", layerNames [layerNum], GetLayerAtomsByInt(layerNum).Count));
			sb.Append (", ");
			layerNum++;
		}
		return string.Format ("Layers: {0}", sb.ToString());
	}
}

public class Layer {
	public string method;
	public string basis;
	public string options;
	public int charge;
	public int multiplicity;

	public int layerOrder;
	private char _layerName;
	public char layerName {
		set {
			if (value == 'H') {
				_layerName = value;
				layerOrder = 0;
			} else if (value == 'M') {
				_layerName = value;
				layerOrder = 1;
			} else if (value == 'L') {
				_layerName = value;
				layerOrder = 2;
			} else {
				Debug.LogError(string.Format("Layer char '{0}' not recognised", value));
			}
		}

		get { return _layerName; }
	}

	public Layer(string method="", string basis="", string options="", char layerName='H', int charge=0, int multiplicity=0 ) {
		this.method = method;
		this.basis = basis;
		this.options = options;
		this.layerName = layerName;
		this.charge = charge;
		this.multiplicity = multiplicity;
	}

	public Atoms GenerateLayerAtoms() {
		//Create a new Atoms object from this layer.
		return new Atoms();
	}

	public override string ToString ()
	{
		return string.Format ("Layer(layerName = {0}, method = {1}, basis = {2}, options = {3}, charge = {4}, multiplicity = {5})", layerName, method, basis, options, charge, multiplicity);
	}
}

public class Link {
	public Atom linkAtom;
	public int linkHostIndex;
	public int linkConnectionIndex;
	public float scaleFactor;

	public Link(Atom linkAtom, int linkHostIndex, int linkConnectionIndex, float scaleFactor=0.7f) {
		this.linkAtom = linkAtom;
		this.linkHostIndex = linkHostIndex;
		this.linkConnectionIndex = linkConnectionIndex;
		this.scaleFactor = scaleFactor;
	}

	public void RecomputeLinkAtomPosition() {
		Atoms atoms = linkAtom.transform.parent.GetComponentInParent<Atoms> ();
		Vector3 vector = atoms.GetVector (linkConnectionIndex, linkHostIndex);
		linkAtom.transform.position = atoms.atomList [linkConnectionIndex].transform.position + vector * scaleFactor;
	}

}