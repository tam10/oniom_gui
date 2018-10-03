using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GaussianCalculator : MonoBehaviour {

	public int nProcessors;
	public int jobMemoryMB;

	public string checkpointPath;
	public string oldCheckpointPath;

	public string killJobLink;
	public int killJobAfter;

	//LAYERS
	public List<Layer> layers;
	
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

		layers = new List<Layer> ();

		oniomOptions = new List<string> ();
		guessOptions = new List<string> ();
		geomOptions = new List<string> ();
		freqOptions = new List<string> ();
		additionalKeywords = new List<string> ();

		title = "Title";

	}

	public void SortLayers() {
	}
}

public class Layer {
	public string method;
	public string basis;
	public string options;
	public int charge;
	public int multiplicity;
	public List<int> atomNums;

	public int layerOrder;
	private char _layer;
	public char layer {
		set {
			if (value == 'H') {
				_layer = value;
				layerOrder = 0;
			} else if (value == 'M') {
				_layer = value;
				layerOrder = 1;
			} else if (value == 'L') {
				_layer = value;
				layerOrder = 2;
			} else {
				Debug.LogError(string.Format("Layer char '{0}' not recognised", value));
			}
		}

		get { return _layer; }
	}

	public Layer(string method="", string basis="", string options="", char layer='H', int charge=0, int multiplicity=0 ) {
		this.method = method;
		this.basis = basis;
		this.options = options;
		this.layer = layer;
		this.charge = charge;
		this.multiplicity = multiplicity;
		this.atomNums = new List<int> ();
	}

	public Atoms GenerateLayerAtoms() {
		//Create a new Atoms object from this layer.
		return new Atoms();
	}

	public override string ToString ()
	{
		return string.Format ("Layer(layer = {0}, method = {1}, basis = {2}, options = {3}, charge = {4}, multiplicity = {5})", layer, method, basis, options, charge, multiplicity);
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
		Vector3 vector = atoms.getVector (linkConnectionIndex, linkHostIndex);
		linkAtom.transform.position = atoms.atomList [linkConnectionIndex].transform.position + vector * scaleFactor;
	}

}