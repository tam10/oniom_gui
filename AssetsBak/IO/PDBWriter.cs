using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Text;
using System.IO;

public class PDBWriter : MonoBehaviour {

	public Data data;

	public void WritePDBFile(Atoms atoms, string path, bool writeConnectivity) {

		StringBuilder sb = new StringBuilder ();

		string newLine = System.Environment.NewLine;
		string format = "ATOM  {0,5} {1,-4} {2,3} {3,1}{4,4}    {5,8:.000}{6,8:.000}{7,8:.000}                      {8,2}" + newLine;
		string pdbName;

		for (int atomNum = 0; atomNum < atoms.size; atomNum++) {

			Atom atom = atoms [atomNum];

			if (atom.pdbName.Length == 4 || data.pdbToElementDict.ContainsKey (atom.pdbName)) {
				pdbName = atom.pdbName;
			} else {
				pdbName = " " + atom.pdbName;
			}

			sb.Append (
				string.Format (
					format,
					atomNum + 1,
					atom.pdbName,
					atom.residueName,
					atom.chainID,
					atom.residueNumber,
					atom.transform.position.x,
					atom.transform.position.y,
					atom.transform.position.z,
					atom.element
				)
			);
		}

		if (writeConnectivity) {

			string cStr = "CONECT";
			format = "{0,5}";
			List<int> neighbours = new List<int> ();

			for (int atomNum = 0; atomNum < atoms.size; atomNum++) {
				sb.Append (cStr);
				sb.Append (string.Format(format, atomNum + 1));

				neighbours = atoms.graph.GetNeighbours (atomNum);
				for (int neighbourNum = 0; neighbourNum < neighbours.Count; neighbourNum++) {
					sb.Append (string.Format(format, neighbours[neighbourNum] + 1));
				}

				sb.Append (newLine);
			}
		}

		File.WriteAllText (path, sb.ToString ());

	}

}
