using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SelectionHalo : MonoBehaviour {

	public SpriteRenderer spriteRenderer;

	public Camera activeCamera;
	public Atoms parent;

	public float sizeRatio;
	private float rotationSpeed;
	private Vector3 rotation;
	public float brightness = 1.5f;

	// Use this for initialization
	public void Awake () {
		sizeRatio = 0.8f;
		activeCamera = Camera.main;
		rotationSpeed = 90f;
		rotation = new Vector3();
	}

	public void SetColor(Color newColor) {
		spriteRenderer.color = newColor * Mathf.LinearToGammaSpace(brightness);
	}

	public void SetSpeed(float newSpeed) {
		rotationSpeed = newSpeed;
	}

	public void SetRadius(float newRadius) {
		transform.localScale = new Vector3(newRadius, newRadius, 1f);
	}

	public void SetAtom(Atom atom) {
		transform.position = atom.transform.position;
		SetRadius(atom.vdwRadius * atom.vdwToSphereRatio * sizeRatio);
		spriteRenderer.enabled = true;
	}

	public void ClearAtom() {
		spriteRenderer.enabled = false;
	}
	
	// Update is called once per frame
	void Update () {
		rotation.x = activeCamera.transform.eulerAngles.x;
		rotation.y = activeCamera.transform.eulerAngles.y;
		rotation.z = - parent.haloZSyncTime * rotationSpeed;
		transform.eulerAngles = rotation;
	}
}
