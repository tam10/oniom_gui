using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SelectionHalo : MonoBehaviour {

	public ParticleSystem hoverHalo;
	public ParticleSystem.EmissionModule emissionModule;
	private ParticleSystem.MainModule main;
	private ParticleSystem.ShapeModule shape;

	public Camera activeCamera;

	public float sizeRatio;

	private Atom atom;

	// Use this for initialization
	void Awake () {
		emissionModule = hoverHalo.emission;
		main = hoverHalo.main;
		shape = hoverHalo.shape;
		sizeRatio = 2f;
	}

	public void SetColor(Color newColor) {
		main.startColor = newColor;
	}

	public void SetSpeed(float newSpeed) {
		shape.arcSpeed = newSpeed;
	}

	public void SetRadius(float newRadius) {
		shape.radius = newRadius;
		main.startSize = newRadius;
	}

	public void SetAtom(Atom atom) {
		transform.position = atom.transform.position;
		SetRadius(atom.vdwRadius * atom.vdwToSphereRatio * sizeRatio);
		hoverHalo.Clear ();
		emissionModule.enabled = true;
		hoverHalo.Play ();
	}

	public void ClearAtom() {
		emissionModule.enabled = false;
	}
	
	// Update is called once per frame
	void Update () {
		transform.rotation = activeCamera.transform.rotation;
	}
}
