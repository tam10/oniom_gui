using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Main : MonoBehaviour {

	public Camera activeCamera;
	public InputHandler inputHandler;

	public Settings globalSettings;
	public FileReader fileReader;
	public FileWriter fileWriter;
	public Data data;
	public PostRender postRenderer;
	public Protonation protonation;
	
	public Tooltip tooltip;
	public StretchBox stretchBox;
	public BendBox bendBox;
	public Status status;

	public string activeJob;

	public List<int> nonStandardResidueList;

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

		/* 
		PDB Protonation and MD

		AddToQueue (LoadAtoms ("/Users/tam10/Notebooks/GFPMutation/2b3p.pdb", "initial"), "Load PDB...");
		AddToQueue (DeactivateAtoms ("initial"), "Deactivate initial atoms...");
		AddToQueue (ProtonateStandardResidues ("initial", "protonated"), "Protonate Standard Residues...");

		//TEMP
		AddToQueue (GetNonStandardResidues("initial", "nonStandard"), "Getting Non-Standard Residues");
		AddToQueue (MergeAtoms("nonStandard", "protonated"), "Merging Non-Standard Residues with Standard Residues");
		AddToQueue (DeactivateAtoms ("nonStandard"), "Deactivate Non-Standard Residue atoms...");

		AddToQueue (CalculateConnectivity ("protonated", 3), "Get connectivity...");
		AddToQueue (ActivateAtoms ("protonated"), "Activate protonated atoms...");

		//AddToQueue (MoveResnumsToLayer("protonated", new List<int>{1,2}, 'H'), "Adding atoms to layer 'H'");
		AddToQueue (SolveGraph ("protonated"), "Solving Graph...");
		AddToQueue (CalculateEnvironment ("protonated"), "Get environment...");

		AddToQueue(AddParametersFromName("protonated", "amber.prm"), "Adding parameters...");

		AddToQueue(GetAmberEnergy("protonated", 0, true), "Getting Amber Energy...");

		AddToQueue(MD("protonated", 200, true), "Running MD...");
		*/

		/* Com file MD */
		

		//AddToQueue (LoadAtoms ("/Users/tam10/Notebooks/GFPMutation/geo3.com", "com"), "Load Gaussian Input...");
		//AddToQueue (LoadAtoms ("/Users/tristanmackenzie/Unity/ONIOM/model_amber2.com", "com"), "Load Gaussian Input...");
		//AddToQueue (LoadAtoms ("/Users/tristanmackenzie/Unity/ONIOM/geo4.com", "com"), "Load Gaussian Input...");
		//AddToQueue (LoadAtoms ("/Users/tristanmackenzie/Unity/ONIOM/nonbon_test.com", "com"), "Load Gaussian Input...");
		AddToQueue (LoadAtoms ("/Users/tristanmackenzie/Unity/ONIOM/hooh_test.com", "com"), "Load Gaussian Input...");
		//AddToQueue (LoadAtoms ("/Users/tristanmackenzie/Unity/ONIOM/o2_test.com", "com"), "Load Gaussian Input...");

		AddToQueue (ActivateAtoms ("com"), "Activate atoms...");

		AddToQueue (SolveGraph ("com"), "Solving Graph...");

		AddToQueue(AddParametersFromName("com", "amber.prm"), "Adding parameters...");

		//AddToQueue(MDVerlet("com", 500, 1, true), "Running MD...");

		AddToQueue(AMBEROpt("com", 250, 1, true), "Optimising...");


		StartCoroutine (RunJobs ());
	}

	void AddToQueue (IEnumerator job, string jobName) {
		jobQueue.Enqueue (new Job(job, jobName));
	}

	public void ClearQueue() {
		jobQueue.Clear();
	}

	IEnumerator RunJobs() {
		while (true) {
			if (jobQueue.Count > 0) {
				Job job = jobQueue.Dequeue ();
				activeJob = job.name;
				status.statusText.text = activeJob;
				yield return null;
				yield return StartCoroutine (job.iEnumerator);
			} else {
				activeJob = "";
				status.statusText.text = "Ready";
				yield return null;
			}
		}
	}

	IEnumerator LoadAtoms(string filename, string atomsName) {
		Atoms atoms = Instantiate<Atoms> (GameObject.FindObjectOfType<PrefabManager>().atomsPrefab);
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
				atoms.active = true;
				tooltip.activeAtoms = atoms;
				stretchBox.activeAtoms = atoms;
				bendBox.activeAtoms = atoms;
				inputHandler.activeAtoms = atoms;
				atoms.postRenderer = postRenderer;

				IEnumerator iEnumerator = GetBestView(atomsName);
				while (iEnumerator.MoveNext()) {}

			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}

		yield return null;

	}

	IEnumerator GetBestView(string atomsName) {

		Atoms atoms;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Could not focus on atoms ({0}) - object is null", atomsName);
			} else {

				int c;
				int size = atoms.size;
				double fov = activeCamera.fieldOfView * Mathf.Deg2Rad;
				double fill_amount = 0.9f;
				double min_distance = 10f;
				double[] positions = new double[size *  3];
				double[] new_camera_position = new double[3];
				for (int atomNum = 0; atomNum < atoms.size; atomNum++) {
					for (c=0; c<3; c++){
						positions[atomNum * 3 + c] = atoms[atomNum].p[c];
					}
				}
				Fortran.getBestCameraView(positions, ref size, new_camera_position, ref fov, ref fill_amount, ref min_distance);
				activeCamera.transform.position = new Vector3((float)new_camera_position[0], (float)new_camera_position[1], (float)new_camera_position[2]);
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
			atoms.active = false;

			if (tooltip.activeAtoms == atoms)
				tooltip.activeAtoms = null;

			if (stretchBox.activeAtoms == atoms)
				stretchBox.activeAtoms = null;

			if (bendBox.activeAtoms == atoms)
				bendBox.activeAtoms = null;

			if (inputHandler.activeAtoms == atoms)
				inputHandler.activeAtoms = null;

			atoms.postRenderer = null;
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return null;
	}

	IEnumerator GetNonStandardResidues(string originalAtomsName, string nonStandardAtomsName) {
		
		Atoms originalAtoms;
		Atoms nonStandardAtoms;
		if (atomsDictionary.TryGetValue (originalAtomsName, out originalAtoms)) {

			if (originalAtoms == null) {
				Debug.LogErrorFormat ("Atoms ({0}) object is null", originalAtomsName);
			} else {
				List<int> nonStandardAtomNums = originalAtoms.graph.SelectIf(resnums:nonStandardResidueList);

				nonStandardAtoms = originalAtoms[nonStandardAtomNums];
				nonStandardAtoms.name = nonStandardAtomsName;
				atomsDictionary[nonStandardAtomsName] = nonStandardAtoms;
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", originalAtomsName);
		}
		yield return null;
	}

	IEnumerator MergeAtoms(string mergeFromAtomsName, string mergeToAtomsName) {
		
		Atoms mergeFromAtoms;
		Atoms mergeToAtoms;

		if (atomsDictionary.TryGetValue (mergeFromAtomsName, out mergeFromAtoms)) {
			if (mergeFromAtoms == null) 
				Debug.LogErrorFormat ("Atoms ({0}) object is null", mergeFromAtomsName);
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", mergeFromAtomsName);
		}


		if (atomsDictionary.TryGetValue (mergeToAtomsName, out mergeToAtoms)) {
			if (mergeToAtoms == null) 
				Debug.LogErrorFormat ("Atoms ({0}) object is null", mergeToAtomsName);
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", mergeToAtomsName);
		}

		foreach(Atom atom in mergeFromAtoms) {
			mergeToAtoms.AddAtom(atom);
		}

		mergeToAtoms.graph.SetAtoms(mergeToAtoms);
		mergeToAtoms.gaussianCalculator.layers.SetAtoms(mergeToAtoms);

		yield return null;
		
	}

	IEnumerator ProtonateStandardResidues(string initialAtomsName, string protonatedAtomsName) {

		Atoms initialAtoms;
		Coroutine coroutine;
		if (atomsDictionary.TryGetValue (initialAtomsName, out initialAtoms)) {

			if (initialAtoms == null) {
				Debug.LogErrorFormat ("Atoms ({0}) could not be protonated - object is null", initialAtomsName);
			} else {
				Atoms protonatedAtoms = Instantiate<Atoms>(GameObject.FindObjectOfType<PrefabManager>().atomsPrefab, transform);
				coroutine = StartCoroutine(protonation.GetProtonatedStandardResidues (initialAtoms, protonatedAtoms, nonStandardResidueList));
				yield return coroutine;
				atomsDictionary [protonatedAtomsName] = protonatedAtoms;
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
				coroutine = StartCoroutine (atoms.graph.SolveGraph ());
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator AddParametersFromName(string atomsName, string paramsName, bool replace=false) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Parameters could not be added to atoms ({0}) - atoms object is null", atomsName);
			} else {
				if (replace) {
					coroutine = StartCoroutine (fileReader.LoadParametersFromNameAsync(paramsName, atoms.graph.parameters));
				} else {
					Parameters newParams = Instantiate<Parameters>(GameObject.FindObjectOfType<PrefabManager>().parametersPrefab, transform);
					IEnumerator iEnumerator = fileReader.LoadParametersFromNameAsync(paramsName, newParams);
					while (iEnumerator.MoveNext()) {}
					atoms.graph.parameters.UpdateParameters(newParams);
				}
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}
	IEnumerator GetAmberEnergy(string atomsName, int order=0, bool suppress=false) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Energy of atoms ({0}) could not calculated - object is null", atomsName);
			} else {
				coroutine = StartCoroutine (atoms.graph.GetAmberEGH ());
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator MD(string atomsName, int steps, bool suppress=false) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Could not initialise MD for ({0})- object is null", atomsName);
			} else {
				coroutine = StartCoroutine (atoms.graph.MD2(steps));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator MDVerlet(string atomsName, int steps, int updateAfterSteps, bool suppress=false) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Could not initialise MD for ({0})- object is null", atomsName);
			} else {
				coroutine = StartCoroutine (atoms.graph.MDVerlet(steps, updateAfterSteps));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator AMBEROpt(string atomsName, int steps, int updateAfterSteps, bool suppress=false) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Could not initialise AMBER Optimisation for ({0})- object is null", atomsName);
			} else {
				coroutine = StartCoroutine (atoms.graph.AMBEROpt(steps, updateAfterSteps));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator MoveAtomsToLayer(string atomsName, List<int> atomNums, char layerName) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Could not move atoms ({0}) to layer ({1}) - object is null", atomsName, layerName);
			} else {
				coroutine = StartCoroutine (atoms.gaussianCalculator.layers.AddAtomsToLayerByNameAsync(layerName, atomNums));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
	}

	IEnumerator MoveResnumsToLayer(string atomsName, List<int> resnums, char layerName) {
		Atoms atoms;
		Coroutine coroutine = null;
		if (atomsDictionary.TryGetValue (atomsName, out atoms)) {

			if (atoms == null) {
				Debug.LogErrorFormat ("Could not move atoms ({0}) to layer ({1}) - object is null", atomsName, layerName);
			} else {
				coroutine = StartCoroutine (atoms.gaussianCalculator.layers.AddAtomsToLayerByNameAsync(layerName, atoms.graph.SelectIf(resnums:resnums)));
			}
		} else {
			Debug.LogErrorFormat ("Atoms ({0}) are not in Atoms Dictionary", atomsName);
		}
		yield return coroutine;
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
