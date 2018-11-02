using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Main : MonoBehaviour {

	public Camera activeCamera;

	public Settings globalSettings;
	public FileReader fileReader;
	public FileWriter fileWriter;
	public Data data;
	public PostRender postRenderer;
	public Protonation protonation;
	public Tooltip tooltip;
	public Status status;

	public string activeJob;

	public Vector3 centre;

	Dictionary<string, Atoms> atomsDictionary;

	Queue<Job> jobQueue;

	// Use this for initialization
	void Start () {

		jobQueue = new Queue<Job> ();

		fileReader.SetGlobalSettings(globalSettings);
		fileReader.SetData(data);

		globalSettings.pdb2pqrCommand = "~/Utilities/pdb2pqr-osx-bin64-2.1.0/pdb2pqr";

		atomsDictionary = new Dictionary<string, Atoms> ();

		atomsDictionary["initial"] = Instantiate<Atoms> (fileReader.atomsPrefab);

		AddToQueue (LoadAtoms ("/Users/tristanmackenzie/Unity/ONIOM/1qxt.pdb", "initial"), "Load PDB...");
		AddToQueue (DeactivateAtoms ("initial"), "Deactivate initial atoms...");
		AddToQueue (ProtonateStandardResidues ("initial", "protonated"), "Protonate Standard Residues...");
		AddToQueue (CalculateConnectivity ("protonated", 5), "Get connectivity...");
		AddToQueue (ActivateAtoms ("protonated"), "Activate protonated atoms...");
		AddToQueue (SolveGraph ("protonated"), "Solve Graph...");
		AddToQueue (CalculateEnvironment ("protonated"), "Get environment...");

		StartCoroutine (RunJobs ());
	}

	void AddToQueue (IEnumerator job, string jobName) {
		jobQueue.Enqueue (new Job(job, jobName));
	}

	IEnumerator LoadAtoms(string filename, string atomsName) {
		Atoms atoms = Instantiate<Atoms> (fileReader.atomsPrefab);
		yield return fileReader.LoadAtomsFileAsync (filename, atoms);
		atoms.transform.parent = transform;
		atomsDictionary [atomsName] = atoms;
		yield return null;

	}

	IEnumerator ActivateAtoms(string atomsName) {

		Atoms atoms;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Atoms ({0}) could not be activated - object is null", atomsName);
			} else {
				atoms.gameObject.SetActive (true);
				atoms.graph.showLines = true;
				tooltip.activeAtoms = atoms;
				atoms.activeCamera = activeCamera;
				atoms.postRenderer = postRenderer;
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return null;
	}

	IEnumerator DeactivateAtoms(string atomsName) {

		Atoms atoms;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Atoms ({0}) could not be deactivated - object is null", atomsName);
			} else {
			atoms.gameObject.SetActive (false);
			atoms.graph.showLines = false;
			atoms.activeCamera = null;
			atoms.postRenderer = null;
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return null;
	}

	IEnumerator ProtonateStandardResidues(string initialAtomsName, string protonatedAtomsName) {

		Atoms atoms;
		if (atomsDictionary.TryGetValue (initialAtomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Atoms ({0}) could not be protonated - object is null", initialAtomsName);
			} else {
				atomsDictionary [protonatedAtomsName] = protonation.GetProtonatedStandardResidues (atoms);
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", initialAtomsName);
		}
		yield return null;
	}

	IEnumerator CalculateConnectivity(string atomsName, int maxIterations=5) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Connectivity of atoms ({0}) could not calculated - object is null", atomsName);
			} else {
				coroutine = StartCoroutine (atoms.graph.CalculateConnectivity (maxIterations));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator CalculateEnvironment(string atomsName, bool overwrite=true, int maxIterations=5) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("AMBER Environment of atoms ({0}) could not calculated - object is null", atomsName);
			} else {
				coroutine = StartCoroutine (atoms.graph.CalculateEnvironment (overwrite, maxIterations));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator SolveGraph(string atomsName, int maxDistance=3) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Graph of atoms ({0}) could not calculated - object is null", atomsName);
			} else {
				coroutine = StartCoroutine (atoms.graph.SolveGraph (maxDistance));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator RunJobs() {
		while (true) {
			if (jobQueue.Count > 0) {
				Job job = jobQueue.Dequeue ();
				activeJob = job.name;
				status.statusText.text = activeJob;
				yield return StartCoroutine (job.iEnumerator);
			} else {
				activeJob = "";
				status.statusText.text = "Ready";
				yield return null;
			}
		}
	}

	void Update() {

	}

	void FixedUpdate() {
		
		//if (ready)
			//atoms.transform.RotateAround (centre, Vector3.up, 20f * Time.deltaTime);
	}

	struct Job {
		public string name;
		public IEnumerator iEnumerator;

			public Job(IEnumerator iEnumerator, string name) {
			this.name = name;
			this.iEnumerator = iEnumerator;
		}
	}

	//public void FixedUpdate() {
	//	testAtom.transform.position += v * Time.deltaTime;
	//}
}
