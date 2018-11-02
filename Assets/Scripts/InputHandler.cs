using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class InputHandler : MonoBehaviour {

	public GameObject haloHolder;
	public SelectionHalo hoverHalo;
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

	public float leftClickDown;
	public float oldLeftClickDown;
	public float rightClickDown;
	public float mouseX;
	public float mouseY;
	public Vector3 rotation;

	private bool dragged;
	private bool clicked;
	private Dictionary<int, SelectionHalo> selectionDict;
	public Color selectionColor = new Color(0.5f, 0f, 1f, 0.2f);

	public float dragThreshold = 0.5f;

	void Awake() {
		frameNum = 0;
		dragged = false;
		clicked = false;
		selectionDict = new Dictionary<int, SelectionHalo>();
	}

	// Update is called once per frame
	void Update () {

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
						hoverHalo.SetAtom(closestAtom);
						hoverHalo.sizeRatio = 2.4f;
						
					}
					
					if (clicked) {
							if (closestAtom.parent.selection.Contains(closestAtomIndex)) {
								closestAtom.parent.selection.Remove(closestAtomIndex);
								GameObject.Destroy(selectionDict[closestAtomIndex].gameObject);
								selectionDict.Remove(closestAtomIndex);
							} else {
								closestAtom.parent.selection.Add(closestAtomIndex);
								SelectionHalo newHalo = Instantiate<SelectionHalo>(hoverHalo, haloHolder.transform);
								newHalo.SetColor(selectionColor);
								newHalo.SetAtom(closestAtom);
								newHalo.sizeRatio = 2.2f;

								selectionDict.Add(closestAtomIndex, newHalo);
							}
						}
				}

			} else {
				if (hitAtom != null)
					hitAtom.parent.hoveredAtom = -1;
				hoverHalo.ClearAtom();
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
