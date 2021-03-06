﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Xml;
using System.Xml.Serialization;
using System.Xml.Linq;

public class Data: MonoBehaviour {

	//Filenames
	public string dataPath = "Assets/Resources/Data/";
	public string amberToElementFilename = "amberToElement.csv";
	public string pdbToElementFilename = "pdbToElement.csv";
	public string gaussianMethodsFilename = "GaussianMethods.xml";

	//Data lists
	public List<AmberToElement> amberToElementData;
	public Dictionary<string, string> pdbToElementDict;

	public List<string> gaussianMethods;

	//Data structures
	public struct AmberToElement {

		public string amberName;
		public string element;
		public int formalCharge;

		public AmberToElement(string amberName, string element, int formalCharge) {
			this.amberName = amberName;
			this.element = element;
			this.formalCharge = formalCharge;
		}

		public override string ToString ()
		{
			return "Amber: " + amberName + ". Element: " + element + ". Formal Charge: " + formalCharge.ToString();
		}
	}

	void Awake() {
		populateAmberData ();
		populatePDBToElementDict ();
		populateGaussianMethods ();
	}

	private void populateGaussianMethods() {

		XDocument gaussianMethodsX = FileIO.ReadXML (dataPath + gaussianMethodsFilename);

		gaussianMethods = new List<string> ();
		XElement methodsList = gaussianMethodsX.Element ("methods");

		foreach (XElement methodName in methodsList.Elements("method"))
			if (!gaussianMethods.Contains(methodName.Value)) {
				gaussianMethods.Add(methodName.Value);
			}

	}

	private void populateAmberData() {
		amberToElementData = new List<AmberToElement> ();
		string[,] data = FileIO.ReadCSV (dataPath + amberToElementFilename, 3);

		for (int lineNum = 0; lineNum < data.GetLength (0); lineNum++) {
			
			string amber = data [lineNum, 0];
			string element = data [lineNum, 1];
			int formalCharge = int.Parse( data [lineNum, 2]);

			amberToElementData.Add( new AmberToElement (amber, element, formalCharge));

		}
	}
	
	private void populatePDBToElementDict() {
		pdbToElementDict = new Dictionary<string, string> ();
		string[,] data = FileIO.ReadCSV (dataPath + pdbToElementFilename, 2);

		for (int lineNum = 0; lineNum < data.GetLength (0); lineNum++) {

			string pdb = data [lineNum, 0];
			string element = data [lineNum, 1];

			pdbToElementDict.Add (pdb, element);

		}

	}


}
