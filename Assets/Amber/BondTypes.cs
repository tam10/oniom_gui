using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Xml;
using System.Xml.Serialization;
using System.Xml.Linq;

public class BondTypes : MonoBehaviour {

	public Dictionary<string, BondType> bondTypesDict;
	public float leeway = 0f;

	public void GetBondTypesFromXML () {

		string path = "Assets/Resources/Data/Bonds.xml";

		bondTypesDict = new Dictionary<string, BondType> ();

		XDocument bondsX = FileIO.ReadXML (path);

		XElement atomsX = bondsX.Element ("elements");

		string atomElement;
		string otherElement;
		string a0El;
		string a1El;

		foreach ( XElement atomEl in atomsX.Elements("atom")) {
			
			atomElement = atomEl.Attribute ("element").Value;

			foreach (XElement otherEl in atomEl.Elements("other")) {

				otherElement = otherEl.Attribute ("element").Value;
				BondType bondType = new BondType ();

				XElement singleEl = otherEl.Element ("single");
				if (singleEl != null)
					bondType.singleThresh = float.Parse(singleEl.Value);

				XElement aromaticEl = otherEl.Element ("aromatic");
				if (aromaticEl != null)
					bondType.aromaticThresh = float.Parse(aromaticEl.Value);

				XElement doubleEl = otherEl.Element ("double");
				if (doubleEl != null)
					bondType.doubleThresh = float.Parse(doubleEl.Value);

				XElement tripleEl = otherEl.Element ("triple");
				if (tripleEl != null)
					bondType.tripleThresh = float.Parse(tripleEl.Value);

				a0El = string.Compare (atomElement, otherElement) > 0 ? atomElement : otherElement;
				a1El = string.Compare (atomElement, otherElement) > 0 ? otherElement : atomElement;
				string key = string.Format("{0}-{1}", a0El, a1El);

				if (bondTypesDict.ContainsKey (key))
					Debug.LogErrorFormat ("Duplicate exists in XML file: {0} ({1})", path, key);
				else 
					bondTypesDict.Add(key, bondType);
			}

		}
	}

	public double GetBondType(Atoms atoms, int a0, int a1) {
		
		string a0El = string.Compare (atoms [a0].element, atoms [a1].element) > 0 ? atoms [a0].element : atoms [a1].element;
		string a1El = string.Compare (atoms [a0].element, atoms [a1].element) > 0 ? atoms [a1].element : atoms [a0].element;
		string key = string.Format("{0}-{1}", a0El, a1El);
		
		BondType bondType;

		if (bondTypesDict.TryGetValue (key, out bondType)) {
		
			float distance = atoms.GetDistance (a0, a1);
		
			if (distance < bondType.tripleThresh + leeway)
				return 3.0;
			if (distance < bondType.doubleThresh + leeway)
				return 2.0;
			if (distance < bondType.aromaticThresh + leeway)
				return 1.5;
			if (distance < bondType.singleThresh + leeway)
				return 1.0;
		
		}

		return 0.0;

	}
}

public class BondType {
	public float singleThresh = 0f;
	public float aromaticThresh = 0f;
	public float doubleThresh = 0f;
	public float tripleThresh = 0f;
}