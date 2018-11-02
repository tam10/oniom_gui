using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Atom : MonoBehaviour {

	public Mesh mesh;
	public MeshFilter meshFilter;
	public SphereCollider connectivityCollider;
	public SphereCollider mouseCollider;
	public Transform sphereTransform;
	public Atoms parent;

	[SerializeField]
	private string _element;
	public string element {
		get { return _element; }
		set { _element = value; }
	}

	[SerializeField]
	private string _pdbName;
	public string pdbName {
		get { return _pdbName; }
		set { _pdbName = value; }
	}

	[SerializeField]
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

	[SerializeField]
	private float _partialCharge;
	public float partialCharge {
		get { return _partialCharge; }
		set { _partialCharge = value; }
	}

	[SerializeField]
	private double _maxValency;
	public double maxValency {
		get { return _maxValency; }
		set { _maxValency = value; }
	}

	[SerializeField]
	private double _valency;
	public double valency {
		get { return _valency; }
		set { _valency = value; }
	}

	[SerializeField]
	private float sphereScale {
		set {
			sphereTransform.localScale = new Vector3 (value, value, value);}
	}

	[SerializeField]
	private float _vdwRadius;
	public float vdwRadius {
		get { return _vdwRadius; }
		set { 
			_vdwRadius = value;
			transform.localScale = new Vector3 (value, value, value);
			sphereScale = value * _vdwToSphereRatio;
		}
	}

	[SerializeField]
	private float _vdwToSphereRatio;
	public float vdwToSphereRatio {
		get { return _vdwToSphereRatio; }
		set { 
			_vdwToSphereRatio = value; 
			sphereScale = value * _vdwToSphereRatio;
		}
	}

	[SerializeField]
	private string _chainID;
	public string chainID {
		get { return _chainID; }
		set { _chainID = value; }
	}

	[SerializeField]
	private string _residueName;
	public string residueName {
		get { return _residueName; }
		set { _residueName = value; }
	}

	[SerializeField]
	private int _residueNumber;
	public int residueNumber {
		get { return _residueNumber; }
		set { _residueNumber = value; }
	}

	[SerializeField]
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

	[SerializeField]
	private Color _color;
	public Color color {
		get { return _color; }
		set { _color = value; }
	}

	private AtomState _atomState;
	public AtomState atomState {
		get { return _atomState; }
		set { _atomState = value; }
	}

	[SerializeField]
	private char _layer;
	public char layer {
		get {return _layer;}
		set {
			renderSphere = parent.globalSettings.renderSphereDict[value];
			_layer = value;
		}
	}

	private bool _renderSphere;
	public bool renderSphere {
		get {return _renderSphere;}
		set {
			if (value) {
				Render();
			}
			_renderSphere = value;
		}
	}

	public Vector2 screenPosition {
		get { return parent.activeCamera.WorldToScreenPoint(transform.position); }
	}

	public void GetDefaults() {
		if (!parent.graph.amberEnvironment.elementEnvironments.ContainsKey(element)) {
			Debug.LogErrorFormat ("Element {0} not present in Element Environments", element);
			return;
		}
		ElementEnvironment defaultEnvironment = parent.graph.amberEnvironment.elementEnvironments [element];
		valency = 0.0;
		maxValency = defaultEnvironment.maxValency;
		color = defaultEnvironment.color;
		vdwRadius = defaultEnvironment.radius;
	}

	// Use this for initialization
	void Awake () {
		
		_element = "X";
		_pdbName = "";
		_amberName = "";

		_formalCharge = 0;
		_partialCharge = 0f;

		_valency = 1.0;
		_maxValency = 1.0;

		_vdwRadius = 1f;
		_vdwToSphereRatio = 0.2f;

		_chainID = "";
		_residueName = "";
		_residueNumber = 0;

		_index = 0;

		_frozen = 0;

		_color = new Color (0f, 0f, 0f);
		_resolution = 1;

		meshFilter = sphereTransform.GetComponent<MeshFilter> ();
		mesh = meshFilter.mesh;

		connectivityCollider.enabled = false;

	}

	public void Render() {
		IcosphereGenerator.GenerateSphere(resolution, mesh, color);

	}
}
