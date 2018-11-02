using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Diagnostics;
using System.IO;

public class Bash : MonoBehaviour {

	public Settings globalSettings;

	public Process GetBashProcess (string command, string directory="") {

		command = string.Format ("-c \"{0}\"", command);

		ProcessStartInfo procInfo = new ProcessStartInfo {
			UseShellExecute = false,
			RedirectStandardOutput = true,
			RedirectStandardError = true,
			CreateNoWindow = false,
			FileName = "/bin/bash",
			WorkingDirectory = Path.Combine(globalSettings.tempFolder, directory),
			Arguments = command
		};

		Process proc = new Process { StartInfo = procInfo };

		return proc;
	}
}
