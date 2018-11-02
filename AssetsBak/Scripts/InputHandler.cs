using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class InputHandler : MonoBehaviour {

	public ParticleSystem hoverHalo;
	public ParticleSystem.EmissionModule emissionModule;
	public Camera activeCamera;
	public Vector3 worldMousePosition;
	public Vector3 screenMousePosition;

	private Ray ray;
	private RaycastHit hitInfo;
	public int numHits;

	private MouseCollider mouseCollider;
	private Atom hitAtom;
	private Atom closestAtom;
	private float angle;
	private float closestAngle;

	private int closestAtomIndex;

	private int frameNum;
	public int framesUntilRefresh = 10;

	void Awake() {
		frameNum = 0;
		emissionModule = hoverHalo.emission;
	}

	// Update is called once per frame
	void Update () {

		//Keep track of mouse position here
		screenMousePosition = Input.mousePosition;
		worldMousePosition = activeCamera.ScreenToWorldPoint(new Vector3(screenMousePosition.x, screenMousePosition.y, activeCamera.nearClipPlane));

		//Handle 6-axis motion
		activeCamera.transform.localPosition += Input.GetAxis ("Pan") * activeCamera.transform.right + Input.GetAxis ("Vertical") * activeCamera.transform.up+ Input.GetAxis ("Forward") * activeCamera.transform.forward;
		activeCamera.transform.Rotate (Input.GetAxis ("Mouse Y") * Input.GetAxis ("Left Click"), Input.GetAxis ("Mouse X") * Input.GetAxis ("Left Click"), Input.GetAxis ("Mouse X") * Input.GetAxis ("Right Click"));
	
		//Get hovered atom
		frameNum++;
		if (frameNum == framesUntilRefresh) {
			frameNum = 0;
		} else {

			ray = activeCamera.ScreenPointToRay (screenMousePosition);

			RaycastHit[] hitInfos;
			hitInfos = Physics.RaycastAll (ray);
			numHits = hitInfos.Length;

			if (numHits > 0) {

				closestAtom = null;
				hitAtom = null;

				closestAngle = 90f;

				for (int i = 0; i < numHits; i++) {

					hitInfo = hitInfos [i];
					mouseCollider = hitInfo.collider.gameObject.GetComponent<MouseCollider> ();

					if (mouseCollider == null)
						continue;

					hitAtom = mouseCollider.parent;

					//Make sure the atom that's closest to the cursor is selected
					angle = Vector3.Angle (ray.direction, hitAtom.transform.position - activeCamera.transform.position);
					if (angle > closestAngle)
						continue;

					closestAtom = hitAtom;
					closestAtomIndex = hitAtom.index;
					closestAngle = angle;
				}

				if (closestAtom != null) {
					if (closestAtomIndex != hitAtom.parent.hoveredAtom) {
						hitAtom.parent.hoveredAtom = closestAtomIndex;
						hoverHalo.transform.position = closestAtom.transform.position;
						hoverHalo.transform.rotation = activeCamera.transform.rotation;

						hoverHalo.Clear ();
						emissionModule.enabled = true;
						hoverHalo.Play ();
					}
				}

			} else {
				if (hitAtom != null)
					hitAtom.parent.hoveredAtom = -1;
				emissionModule.enabled = false;
			}
				

		}

		//if (Physics.Raycast (ray, out hitInfo)) {
		//	MouseCollider hoveredAtomCollider = hitInfo.collider.gameObject.GetComponent<MouseCollider> ();
		//
		//	if (hoveredAtomCollider != null) {
		//		hoveredAtomCollider.parent.parent.hoveredAtom = hoveredAtomCollider.parent.index;
		//	}
		//}
	}
}
