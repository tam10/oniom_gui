using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FileWriter : MonoBehaviour {

	public PDBWriter pdbWriter;

	public void WritePDBFile(Atoms atoms, string path, bool writeConnectivity) {
		pdbWriter.WritePDBFile (atoms, path, writeConnectivity);
	}

}
