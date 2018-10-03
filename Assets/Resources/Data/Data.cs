using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class Data: MonoBehaviour {

	//Filenames
	public string dataPath = "Assets/Resources/Data/";
	public string amberToElementFilename = "amberToElement.csv";

	//Data lists
	public List<AmberToElement> amberToElementData;

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
	



}
