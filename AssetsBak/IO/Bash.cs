using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Diagnostics;

public class Bash : MonoBehaviour {

	public Process GetBashProcess (string command) {

		//Escape escape characters
		command = command.Replace ("\"", "\"\"");

		ProcessStartInfo procInfo = new ProcessStartInfo {
			UseShellExecute = false,
			RedirectStandardOutput = true,
			RedirectStandardError = true,
			CreateNoWindow = false,
			FileName = "/bin/bash",
			Arguments = string.Format ("-c \"{0}\"", command)
		};

		Process proc = new Process { StartInfo = procInfo };

		return proc;
	}
}
