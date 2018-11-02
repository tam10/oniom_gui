using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Xml;
using System.Xml.Serialization;
using System.Xml.Linq;

public class AmberEnvironment : MonoBehaviour {

	public Dictionary<string, ElementEnvironment> elementEnvironments;

	public void GetEnvironmentFromXML () {

		elementEnvironments = new Dictionary<string, ElementEnvironment> ();

		XDocument envX = FileIO.ReadXML ("Assets/Resources/Data/Elements.xml");

		XElement atomsX = envX.Element ("elements");

		foreach (XElement atomEl in atomsX.Elements("atom")) {
			ElementEnvironment elementEnvironment = new ElementEnvironment();

			elementEnvironment.element = atomEl.Attribute ("element").Value;
			elementEnvironment.radius = float.Parse (atomEl.Element ("radius").Value);
			elementEnvironment.maxValency = double.Parse (atomEl.Element ("maxValency").Value);
			elementEnvironment.mass = float.Parse (atomEl.Element ("mass").Value);

			float r = float.Parse(atomEl.Element ("red").Value);
			float g = float.Parse(atomEl.Element ("green").Value);
			float b = float.Parse(atomEl.Element ("blue").Value);
			elementEnvironment.color = new Color (r, g, b);


			foreach (XElement stateEl in atomEl.Elements("state")) {
				elementEnvironment.atomStateDict[stateEl.Attribute("name").Value] = AtomStateFromXmlElement(stateEl);
			}

			foreach (XElement environmentEl in atomEl.Elements("environment")) {
				AtomEnvironment atomEnvironment = AtomEnvironmentFromXmlElement (environmentEl, elementEnvironment.atomStateDict);
				elementEnvironment.environmentDict[atomEnvironment.name] = atomEnvironment;
			}

			elementEnvironments [elementEnvironment.element] = elementEnvironment;
		}
	}

	AtomState AtomStateFromXmlElement(XElement stateEl) {
		AtomState atomState = new AtomState();
		atomState.name = stateEl.Attribute ("name").Value;
		atomState.singles = int.Parse(stateEl.Element ("single").Value);
		atomState.doubles = int.Parse(stateEl.Element ("double").Value);
		atomState.triples = int.Parse(stateEl.Element ("triple").Value);
		atomState.aromatics = int.Parse(stateEl.Element ("aromatic").Value);
		atomState.charge = int.Parse(stateEl.Element ("charge").Value);
		atomState.radical = int.Parse(stateEl.Element ("radical").Value);
		return atomState;
	}

	AtomEnvironment AtomEnvironmentFromXmlElement(XElement environmentEl, Dictionary<string, AtomState> atomStateDict) {
		AtomEnvironment atomEnvironment = new AtomEnvironment();
		atomEnvironment.name = environmentEl.Attribute ("name").Value;
		atomEnvironment.amber = environmentEl.Element ("amber").Value;
		atomEnvironment.radius = float.Parse(environmentEl.Element ("radius").Value);

		float r = float.Parse(environmentEl.Element ("red").Value);
		float g = float.Parse(environmentEl.Element ("green").Value);
		float b = float.Parse(environmentEl.Element ("blue").Value);
		atomEnvironment.color = new Color (r, g, b);

		string stateName = environmentEl.Element ("state").Value;
		if (!atomStateDict.ContainsKey (stateName)) {
			string[] keyArr = new string[atomStateDict.Keys.Count];
			atomStateDict.Keys.CopyTo (keyArr, 0);
			Debug.LogErrorFormat (
				"State {0} not available for Atom Environment {1}. Available states for element {2}: {3}", 
				stateName, 
				atomEnvironment.name, 
				environmentEl.Parent.Attribute("element").Value, 
				string.Join(", ", keyArr)
			);
		} else {
			atomEnvironment.atomState = atomStateDict [stateName];
		}

		foreach (XElement conditionEl in environmentEl.Elements("condition")) {
			atomEnvironment.environmentConditions.Add(EnvironmentConditionFromXmlElement (conditionEl, atomStateDict));
		}

		return atomEnvironment;
	}

	EnvironmentCondition EnvironmentConditionFromXmlElement(XElement conditionEl, Dictionary<string, AtomState> atomStateDict) {
		EnvironmentCondition environmentCondition = new EnvironmentCondition ();

		foreach (XElement elementCondition in conditionEl.Elements("element"))
			foreach (string elementString in elementCondition.Value.Split(new char[] {','}, System.StringSplitOptions.RemoveEmptyEntries)) {
				environmentCondition.elementConditions.Add (elementString);
			}

		foreach (XElement amberCondition in conditionEl.Elements("amber"))
			foreach (string amberString in amberCondition.Value.Split(new char[] {','}, System.StringSplitOptions.RemoveEmptyEntries)) {
				environmentCondition.amberConditions.Add (amberString);
			}

		foreach (XElement pdbCondition in conditionEl.Elements("pdb"))
			foreach (string pdbString in pdbCondition.Value.Split(new char[] {','}, System.StringSplitOptions.RemoveEmptyEntries)) {
				environmentCondition.pdbConditions.Add (pdbString);
			}

		foreach (XElement residueCondition in conditionEl.Elements("residue"))
			foreach (string residueString in residueCondition.Value.Split(new char[] {','}, System.StringSplitOptions.RemoveEmptyEntries)) {
				environmentCondition.residueConditions.Add (residueString);
			}

		foreach (XElement ringCondition in conditionEl.Elements("inRing"))
			foreach (string ringString in ringCondition.Value.Split(new char[] {','}, System.StringSplitOptions.RemoveEmptyEntries)) {
				environmentCondition.ringConditions.Add (int.Parse (ringString));
			}

		foreach (XElement neighbourCountCondition in conditionEl.Elements("neighbourCount"))
			foreach (string neighbourCountString in neighbourCountCondition.Value.Split(new char[] {','}, System.StringSplitOptions.RemoveEmptyEntries)) {
				environmentCondition.neighbourCountConditions.Add (int.Parse (neighbourCountString));
			}

		foreach (XElement neighboursEnvironmentCondition in conditionEl.Elements("neighbour")) 
			environmentCondition.neighboursEnvironmentConditions.Add(EnvironmentConditionFromXmlElement(neighboursEnvironmentCondition, atomStateDict));

		return environmentCondition;
	}


}

public class ElementEnvironment {
	public Dictionary<string, AtomState> atomStateDict = new Dictionary<string, AtomState>();
	public Dictionary<string, AtomEnvironment> environmentDict = new Dictionary<string, AtomEnvironment> ();
	public string element;
	public float radius;
	public Color color;
	public double maxValency;
	public float mass;
}

public class AtomState {
	public string name;
	public int singles;
	public int doubles;
	public int triples;
	public int aromatics;
	public int charge;
	public int radical;
}

public class EnvironmentCondition {
	public string name;
	public List<string> elementConditions = new List<string>();
	public List<string> amberConditions = new List<string>();
	public List<string> pdbConditions = new List<string>();
	public List<string> residueConditions = new List<string>();
	public List<int> neighbourCountConditions = new List<int>();

	public List<EnvironmentCondition> neighboursEnvironmentConditions = new List<EnvironmentCondition>();
	public List<int> ringConditions =  new List<int>();


	public string IntListToString(List<int> intList) {
		string[] ss = new string[intList.Count];
		for (int i = 0; i < intList.Count; i++)
			ss [i] = intList [i].ToString ();
		return string.Join (" ", ss);
	}
	public string StringListToString(List<string> stringList) {
		string[] ss = new string[stringList.Count];
		for (int i = 0; i < stringList.Count; i++)
			ss [i] = stringList [i].ToString ();
		return string.Join (" ", ss);
	}
	//Conditions:
	//0: True
	//1: False
	//2: Indeterminable (depends on neighbours)
	public int EvaluateConditions(Atoms atoms, int index, bool debug=false) {
		int elementResult = EvaluateElementCondition (atoms, index);
		if (debug)
			Debug.LogFormat ("Element condition: {0}. Element: {1}", StringListToString(elementConditions), atoms [index].element);
		if (elementResult > 0)
			return elementResult;

		int amberResult = EvaluateAMBERCondition (atoms, index);
		if (debug)
			Debug.LogFormat ("Amber condition: {0}. Amber: {1}", StringListToString(amberConditions), atoms [index].amberName);
		if (amberResult > 0)
			return amberResult;

		int pdbResult = EvaluatePDBCondition (atoms, index);
		if (debug)
			Debug.LogFormat ("PDB condition: {0}. PDB: {1}", StringListToString(pdbConditions), atoms [index].pdbName);
		if (pdbResult > 0)
			return pdbResult;

		int residueResult = EvaluateResidueCondition (atoms, index);
		if (debug)
			Debug.LogFormat ("Residue condition: {0}. Residue: {1}", StringListToString(residueConditions), atoms [index].residueName);
		if (residueResult > 0)
			return residueResult;

		int neighbourCountResult = EvaluateNeighbourCountCondition (atoms, index);
		if (debug)
			Debug.LogFormat ("Neighbour Count condition: {0}. Neighbour Count: {1}", IntListToString(neighbourCountConditions), atoms.graph.neighbourCount [index]);
		if (neighbourCountResult > 0)
			return neighbourCountResult;

		int ringResult = EvaluateRingCondition (atoms, index);
		if (debug)
			Debug.LogFormat ("Ring condition: {0}. 5?: {1}, 6?: {2}", StringListToString(elementConditions), atoms.graph.fiveMemberedRings [index], atoms.graph.sixMemberedRings [index]);
		if (ringResult > 0)
			return ringResult;

		//DEBUG START
		if (debug) {
			Debug.LogFormat("Index: {0}, PDB: {1}", index, atoms[index].pdbName);
			int neighbourEnvironmentResult2;
			List<int> neighbourPool2 = atoms.graph.GetNeighbours (index);
			foreach (EnvironmentCondition neighbourEnvironmentCondition in neighboursEnvironmentConditions) {
				foreach (int neighbourIndex in neighbourPool2) {
					Debug.LogFormat ("Evaluate neighbour: Index: {0}", neighbourIndex);
					neighbourEnvironmentResult2 = neighbourEnvironmentCondition.EvaluateConditions (atoms, neighbourIndex, debug);
					Debug.LogFormat ("Evaluate neighbour: Index: {0}, Result {1}", neighbourIndex, neighbourEnvironmentResult2);

					if (neighbourEnvironmentResult2 == 0) {
						neighbourPool2.Remove (neighbourIndex);
						goto NEXT;
					}
				}
				Debug.Log ("Failed");
				return 1;
				NEXT:
				Debug.Log ("Passed");
				;
			}

			return 0;
		}
		//DEBUG END

		int neighbourEnvironmentResult;
		List<int> neighbourPool = atoms.graph.GetNeighbours(index);
		foreach (EnvironmentCondition neighbourEnvironmentCondition in neighboursEnvironmentConditions) {
			foreach (int neighbourIndex in neighbourPool) {
				neighbourEnvironmentResult = neighbourEnvironmentCondition.EvaluateConditions (atoms, neighbourIndex);
				if (neighbourEnvironmentResult == 0) {
					neighbourPool.Remove (neighbourIndex);
					goto NEXT;
				}
			}
			return 1;
			NEXT:;
		}

		return 0;
	}

	public int EvaluateElementCondition(Atoms atoms, int index) {
		//Debug.LogFormat ("elementCondition: {0}. Element: {1}", elementCondition, atoms [index].element);
		if (elementConditions.Count == 0)
			return 0;

		foreach (string elementCondition in elementConditions) {
			if (atoms [index].element == elementCondition)
				return 0;
		}
		return 1;
	}

	public int EvaluateAMBERCondition(Atoms atoms, int index) {
		//Debug.LogFormat ("amberCondition: {0}. Amber: {1}", amberCondition, atoms [index].amberName);
		if (amberConditions.Count == 0)
			return 0;

		foreach (string amberCondition in amberConditions) {
			if (atoms [index].amberName == amberCondition)
				return 0;
		}
		return 1;
	}

	public int EvaluatePDBCondition(Atoms atoms, int index) {
		//Debug.LogFormat ("pdbCondition: {0}. PDB: {1}", pdbCondition, atoms [index].pdbName);
		if (pdbConditions.Count == 0)
			return 0;

		foreach (string pdbCondition in pdbConditions) {
			if (atoms [index].pdbName == pdbCondition)
				return 0;
		}
		return 1;
	}
		
	public int EvaluateResidueCondition(Atoms atoms, int index) {
		//Debug.LogFormat ("residueCondition: {0}. Residue: {1}", residueCondition, atoms [index].residueName);
		if (residueConditions.Count == 0)
			return 0;

		foreach (string residueCondition in residueConditions) {
			if (atoms [index].residueName == residueCondition)
				return 0;
		}
		return 1;
	}

	public int EvaluateNeighbourCountCondition(Atoms atoms, int index) {
		//Debug.LogFormat ("neighbourCountCondition: {0}. Neighbours: {1}", neighbourCountCondition, atoms.graph.neighbourCount[index]);
		if (neighbourCountConditions.Count == 0)
			return 0;

		foreach (int neighbourCountCondition in neighbourCountConditions) {
			if (atoms.graph.neighbourCount[index] == (double)neighbourCountCondition)
				return 0;
		}
		return 1;
	}

	public int EvaluateRingCondition(Atoms atoms, int index) {
		//Debug.LogFormat ("ringCondition: {0}. 5?: {1}. 6?: {2}", ringCondition, atoms.graph.fiveMemberedRings[index], atoms.graph.sixMemberedRings[index]);

		if (ringConditions.Count == 0)
			return 0;

		foreach (int ringCondition in ringConditions) {
			if (ringCondition == 5) {
				if (atoms.graph.fiveMemberedRings.Contains(atoms[index]))
					return 0;
			} else if (ringCondition == 6) {
				if (atoms.graph.sixMemberedRings.Contains(atoms[index]))
					return 0;
			} else {
				Debug.LogError ("Only 5- and 6-membered rings can be tested");
			}
		}
		return 1;

	}
}

public class AtomEnvironment {
	public string name;
	public AtomState atomState;
	public string amber;
	public float radius;
	public Color color;

	public List<EnvironmentCondition> environmentConditions = new List<EnvironmentCondition>();


}