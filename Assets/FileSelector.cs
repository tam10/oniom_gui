using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using TMPro;
using System.IO;

public class FileSelector : MonoBehaviour {

	public Camera activeCamera;
	public Canvas canvas;
	public Button confirmButton;
	public Button cancelButton;
	public TextMeshProUGUI tmpText;

	public string confirmedText;

	private string[] filetypes;

	private List<string> autocompleteList;

	// Use this for initialization
	void Awake () {
		activeCamera = Camera.main;
		canvas.worldCamera = activeCamera;

		confirmButton.onClick.AddListener(Confirm);
		cancelButton.onClick.AddListener(Cancel);
	}

	void SetFileTypes(string[] filetypes) {
		this.filetypes = filetypes;
	}

	void Cancel() {
		confirmedText = "";
		GameObject.Destroy(this);
	}
	
	void Confirm() {
		if (CheckText(tmpText.ToString())) {
			GameObject.Destroy(this);
		}
	}

	string Autocomplete(string text) {
		FileInfo fi = new FileInfo(text);
		string autocompletedText = text;
		string directory = fi.DirectoryName;
		string partial = fi.Name;

		string[] allFiles = Directory.GetFiles(directory);
		string[] allDirectories = Directory.GetDirectories(directory);

		List<string> matchedFilesDirectories = new List<string>();

		foreach (string f in allFiles) {
			if (f.Contains(partial)) {
				matchedFilesDirectories.Add(f);
			}
		}

		foreach (string d in allDirectories) {
			if (d.Contains(partial)) {
				matchedFilesDirectories.Add(d);
			}
		}

		if (matchedFilesDirectories.Count == 1) {
			autocompletedText = Path.Combine(directory,  matchedFilesDirectories[0]);
		}

		return autocompletedText;

	}

	bool CheckText(string text) {

		confirmedText = "";

		if (! File.Exists(text)) {
			return false;
		}

		if (filetypes.Length > 0) {
			string extension = Path.GetExtension(text);
			bool hasExtension = false;
			foreach (string filetype in filetypes) {
				if (extension == filetype) {
					hasExtension = true;
				}
			}
			if (! hasExtension) {
				return false;
			}
		}

		confirmedText = text;
		return true;
	}

	void Update() {
		if (Input.GetKeyDown(KeyCode.Tab)) {
			tmpText.SetText (Autocomplete(tmpText.ToString()));
		}
	}
}
