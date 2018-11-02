using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System.Xml;
using System.Xml.Serialization;
using System.Xml.Linq;

public static class FileIO {

	// Return array of strings for each line in file
	public static string[] Readlines(string filename, string commentString="") {
		if (!File.Exists (filename))
			Debug.LogError ("File " + filename + " does not exist.");

		string[] allLines = File.ReadAllLines (filename);

		List<string> lines = new List<string> ();

		for (int lineNum = 0; lineNum < allLines.Length; lineNum++) {
			string line = allLines [lineNum];

			line.Trim ();

			if (commentString != "") {
				if (!line.StartsWith (commentString))
					lines.Add (line);
			} else 
				lines.Add (line);
		}

		return lines.ToArray ();

	}

	public static string[,] ReadCSV(string filename, int columns) {

		string[] lines = Readlines (filename, "#");
		string[,] data = new string[lines.Length, columns];

		for (int lineNum = 0; lineNum < lines.Length; lineNum++) {

			string[] splitLine = lines [lineNum].Split (new []{ "," }, System.StringSplitOptions.RemoveEmptyEntries);

			for (int columnNum = 0; columnNum < columns; columnNum++) {
				if (columnNum < splitLine.Length) {
					data [lineNum, columnNum] = splitLine [columnNum];
				} else {
					data [lineNum, columnNum] = string.Empty;
				}
			}

		}

		return data;
	}

	public static XDocument ReadXML(string filename) {

		XDocument xmlObj = XDocument.Load (filename);
		return xmlObj;

	}
}
