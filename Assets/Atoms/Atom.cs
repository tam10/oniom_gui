using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Atom : MonoBehaviour {

	public Mesh mesh;
	public MeshFilter meshFilter;

	private string _element;
	public string element {
		get { return _element; }
		set { _element = value; }
	}

	private string _pdbName;
	public string pdbName {
		get { return _pdbName; }
		set { _pdbName = value; }
	}

	private string _amberName;
	public string amberName {
		get { return _amberName; }
		set { _amberName = value; }
	}


	private int _formalCharge;
	public int formalCharge {
		get { return _formalCharge; }
		set { _formalCharge = value; }
	}

	private float _partialCharge;
	public float partialCharge {
		get { return _partialCharge; }
		set { _partialCharge = value; }
	}

	private int _valency;
	public int valency {
		get { return _valency; }
		set { _valency = value; }
	}

	private float _vdwRadius;
	public float vdwRadius {
		get { return _vdwRadius; }
		set { _vdwRadius = value; }
	}

	private float _vdwToSphereRatio;
	public float vdwToSphereRatio {
		get { return _vdwToSphereRatio; }
		set { _vdwToSphereRatio = value; }
	}

	private string _chainID;
	public string chainID {
		get { return _chainID; }
		set { _chainID = value; }
	}

	private string _residueName;
	public string residueName {
		get { return _residueName; }
		set { _residueName = value; }
	}

	private int _residueNumber;
	public int residueNumber {
		get { return _residueNumber; }
		set { _residueNumber = value; }
	}

	private int _index;
	public int index {
		get { return _index; }
		set { _index = value; }
	}

	private int _frozen;
	public int frozen {
		get { return _frozen; }
		set { _frozen = value; }
	}

	private int _resolution;
	public int resolution {
		get { return _resolution; }
		set { _resolution = value; }
	}

	private Color _color;
	public Color color {
		get { return _color; }
		set { _color = value; }
	}



	// Use this for initialization
	void Awake () {
		
		_element = "X";
		_pdbName = "";
		_amberName = "";

		_formalCharge = 0;
		_partialCharge = 0f;

		_valency = 1;

		_vdwRadius = 1f;
		_vdwToSphereRatio = 0.4f;

		_chainID = "";
		_residueName = "";
		_residueNumber = 0;

		_index = 0;

		_frozen = 0;

		_color = new Color (1f, 0f, 1f);
		_resolution = 1;

		meshFilter = GetComponent<MeshFilter> ();
		mesh = meshFilter.mesh;
	}
	
	// Update is called once per frame
	void Update () {
		
	}

	public void Render() {

		//Generate sphere at (0,0,0)
		IcosphereGenerator.GenerateSphere(resolution, mesh, color);
		float r = vdwToSphereRatio * vdwRadius;
		transform.localScale = new Vector3 (r, r, r);


	}
}
