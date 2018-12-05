using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class InputHandler : MonoBehaviour {

	public Camera activeCamera;
	public Main main;
	public Vector3 worldMousePosition;
	public Vector3 screenMousePosition;

	private Ray ray;
	private RaycastHit hitInfo;
	public int numHits;

	private MouseCollider mouseCollider;
	public Atoms activeAtoms;
	private Atom hitAtom;
	private Atom closestAtom;
	private float angle;
	private float closestAngle;

	private int closestAtomIndex;

	private int frameNum;
	public int framesUntilRefresh = 10;

	public float leftClickDown;
	public float oldLeftClickDown;
	public float rightClickDown;
	public float mouseX;
	public float mouseY;
	public Vector3 rotation;

	private bool dragged;
	private bool clicked;

	public float dragThreshold = 0.5f;

	void Awake() {
		frameNum = 0;
		dragged = false;
		clicked = false;
	}

	// Update is called once per frame
	void Update () {

		if (Input.GetKeyDown(KeyCode.Q)) {
			
			Fortran.deallocate_all();
			main.ClearQueue();
			 #if UNITY_EDITOR
				UnityEditor.EditorApplication.isPlaying = false;
			#else
				Application.Quit();
			#endif
		}

		//Handle key presses
		if (Input.GetKeyDown(KeyCode.Escape)) {
			if (activeAtoms)
				activeAtoms.ClearSelection();
		}

		//Keep track of mouse position here
		screenMousePosition = Input.mousePosition;
		worldMousePosition = activeCamera.ScreenToWorldPoint(new Vector3(screenMousePosition.x, screenMousePosition.y, activeCamera.nearClipPlane));

		leftClickDown = Input.GetAxis("Left Click");
		rightClickDown = Input.GetAxis("Right Click");

		//Handle 6-axis motion
		mouseX = Input.GetAxis ("Mouse X");
		mouseY = Input.GetAxis ("Mouse Y");

		rotation.x = leftClickDown * mouseY;
		rotation.y = leftClickDown * mouseX;
		rotation.z = rightClickDown * mouseX;

		activeCamera.transform.localPosition += Input.GetAxis ("Pan") * activeCamera.transform.right + Input.GetAxis ("Vertical") * activeCamera.transform.up + Input.GetAxis ("Forward") * activeCamera.transform.forward;
		activeCamera.transform.Rotate (rotation);

		//Differentiate between click and drag
		clicked = false;
		if (leftClickDown == 0) {
			clicked = (oldLeftClickDown != 0 && !dragged);
			dragged = false;
		} else {
			if (!dragged)
				dragged = (mouseX * mouseX + mouseY * mouseY > dragThreshold);
		}	
		oldLeftClickDown = leftClickDown;
		
		//Get hovered atom
		if (activeAtoms != null) {
			frameNum++;
			if (frameNum == framesUntilRefresh) {
				frameNum = 0;
			} else {

				ray = activeCamera.ScreenPointToRay (screenMousePosition);

				//MOVE ALL HOVER HALOES TO ATOMS

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

						if (hitAtom.parent != activeAtoms)
							continue;

						//Make sure the atom that's closest to the cursor is selected
						angle = Vector3.Angle (ray.direction, hitAtom.transform.position - activeCamera.transform.position);
						if (angle > closestAngle)
							continue;

						closestAtom = hitAtom;
						closestAtomIndex = hitAtom.index;
						closestAngle = angle;
					}

					if (closestAtom != null) {
						if (closestAtomIndex != activeAtoms.hoveredAtom) {
							activeAtoms.hoveredAtom = closestAtomIndex;
						}
						
						if (clicked) {
							activeAtoms.ToggleSelect(closestAtomIndex);
						}
					}

				} else {
					if (hitAtom != null)
						activeAtoms.hoveredAtom = -1;
				}
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
