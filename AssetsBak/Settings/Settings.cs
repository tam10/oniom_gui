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

	public List<string> baseDirectories;

	//PDB2PQR
	public string pdb2pqrCommand = "pdb2pqr";
	public List<string> pdb2pqrOptions;
	public float pH = 7.0f;


	//NOT YET IN FILE
	public int maxLinesToYield = 250;

	public Dictionary<char, bool> renderSphereDict = new Dictionary<char, bool> {{'H', true}, {'M', false}, {'L', false}};
	public Dictionary<char, bool> renderWireframeDict = new Dictionary<char, bool> {{'H', false}, {'M', true}, {'L', true}};

	// Use this for initialization
	void Awake () {

		baseDirectories = new List<string> { Application.persistentDataPath, "", "Assets/Settings" };

		GetSettings ();
	}

	void GetSettings() {
		GetProjectSettings ();
		GetGraphicsSettings ();
	}

	void GetProjectSettings() {
		projectSettingsPath = GetPath (baseDirectories, projectSettingsFilename);
		if (projectSettingsPath == "")
			Debug.LogErrorFormat ("Could not find Projects Settings File. Using {0}", projectSettingsFilename);
		
		XDocument sX = FileIO.ReadXML (projectSettingsPath);
		XElement projectSettingsX = sX.Element ("projectSettings");


		//Graphics
		string graphicsSettingsFilename = projectSettingsX.Element("graphics").Element("graphicsSettingsFilename").Value;
		graphicsSettingsPath = GetPath(baseDirectories, graphicsSettingsFilename);
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

	string GetPath(List<string> directories, string filename) {

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
