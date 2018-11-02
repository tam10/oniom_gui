using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Xml;
using System.Xml.Serialization;
using System.Xml.Linq;

public class Settings : MonoBehaviour {

	public string projectSettingsFilename = "ProjectSettings.xml";
	public string projectSettingsPath = "";
	public string graphicsSettingsPath = "";

	public int ballResolution = 1;
	public int stickResolution = 1;

	public List<string> baseSettingsDirectories;
	public List<string> baseDataDirectories;

	//PDB2PQR
	public string pdb2pqrCommand = "pdb2pqr";
	public List<string> pdb2pqrOptions;
	public float pH = 7.0f;


	//NOT YET IN FILE
	public int maxLinesToYield = 250;
	public float fogStartDistance = 1f;
	public float fogEndDistance = 10f;
	public float fogRatio = 0.2f;
	public float mdTimeStep = 0.1f;
	public float maxNonBondingCutoff = 15f;
	public float angstromToBohr = 1.8897259886f;
	public float mdDampingFactor = 0.5f;

	//Protonation
	private string _tempFolder;
	public string tempFolder {
		get {
			if (!Directory.Exists(_tempFolder))
				Directory.CreateDirectory(_tempFolder);
			return _tempFolder;
		}
		set {_tempFolder = value;}
	}
	public string standardResidueProtonationFilename;

	//Layers
	public Dictionary<int, bool> layerSphereRenderDict = new Dictionary<int, bool>() {{0, true}, {1, true}, {2, false}};
	public Dictionary<int, bool> layerWireRenderDict = new Dictionary<int, bool>() {{0, false}, {1, false}, {2, true}};

	// Use this for initialization
	void Awake () {

		standardResidueProtonationFilename = "SR_protonation";
		tempFolder = Path.Combine( Application.persistentDataPath, "Temp");
		baseSettingsDirectories = new List<string> { Path.Combine(Application.persistentDataPath, "Settings"), "Settings", Path.Combine("Assets", "Settings") };
		baseDataDirectories = new List<string> { Path.Combine(Application.persistentDataPath, Path.Combine("Resources", "Data")), "Data", Path.Combine(Path.Combine("Assets", "Resources"), "Data") };

		GetSettings ();
	}

	void GetSettings() {
		GetProjectSettings ();
		GetGraphicsSettings ();
	}

	void GetProjectSettings() {
		projectSettingsPath = GetPath (baseSettingsDirectories, projectSettingsFilename);
		if (projectSettingsPath == "")
			Debug.LogErrorFormat ("Could not find Projects Settings File. Using {0}", projectSettingsFilename);
		
		XDocument sX = FileIO.ReadXML (projectSettingsPath);
		XElement projectSettingsX = sX.Element ("projectSettings");


		//Graphics
		string graphicsSettingsFilename = projectSettingsX.Element("graphics").Element("graphicsSettingsFilename").Value;
		graphicsSettingsPath = GetPath(baseSettingsDirectories, graphicsSettingsFilename);
		if (graphicsSettingsPath == "")
			Debug.LogErrorFormat ("Could not find Graphics Settings File. Using {0}", graphicsSettingsFilename);

		//Protonation
		pdb2pqrOptions = new List<string> ();
		XElement pdb2pqrX = projectSettingsX.Element ("protonation").Element ("pdb2pqr");
		pdb2pqrCommand = pdb2pqrX.Element ("command").Value;
		pH = float.Parse (pdb2pqrX.Element ("pH").Value);
		foreach (XElement pdb2pqroption in pdb2pqrX.Elements("option")) {
			pdb2pqrOptions.Add (pdb2pqroption.Value);
		}
		
	}

	void GetGraphicsSettings() {
		int qualitySetting = QualitySettings.GetQualityLevel ();

		XDocument sX = FileIO.ReadXML (graphicsSettingsPath);
		foreach (XElement el in sX.Element("graphics").Elements("globalResolution")) {
			if (int.Parse (el.Attribute ("value").Value) == qualitySetting) {
				ballResolution = int.Parse (el.Element ("ballResolution").Value);
				stickResolution = int.Parse (el.Element ("stickResolution").Value);
			}
		}
		return;

	}

	public string GetPath(List<string> directories, string filename) {

		string fullpath = "";

		foreach (string directory in directories) {
			fullpath = Path.Combine (directory, filename);

			if (File.Exists (fullpath)) {
				return fullpath;
			}
		}

		return "";
	}

}
