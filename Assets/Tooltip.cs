using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Tooltip : MonoBehaviour {

	public Canvas canvas;
	public InputHandler inputHandler;
	public RectTransform rectTransform;
	public GameObject tooltipContainer;

	public Text element;
	public Text index;
	public Text amber;
	public Text pdb;
	public Text residue;
	public Text valency;

	public int hoveredIndex;
	public int oldHoveredIndex;

	public Atoms activeAtoms;
	private Atom hoverAtom;

	private float distanceFromMouse;
	private Vector3 vectorToCentre;
	public Vector3 offset;

	public bool track;

	void Awake() {
		rectTransform = GetComponent<RectTransform> ();
		distanceFromMouse = 1.1f * Mathf.Sqrt((rectTransform.rect.width * rectTransform.rect.width) + (rectTransform.rect.width * rectTransform.rect.width)); 
		offset = new Vector3 (0f, 0f, canvas.planeDistance);

		track = false;

		float pad = 20f;
		Vector2 initPosition = new Vector2(Screen.width - pad - rectTransform.rect.width * 0.5f, pad + rectTransform.rect.height * 0.5f);
		rectTransform.anchoredPosition = initPosition;

	}

	void Update () {

		if (track) {
			vectorToCentre.x = Screen.width * 0.5f - inputHandler.screenMousePosition.x;
			vectorToCentre.y = Screen.height * 0.5f - inputHandler.screenMousePosition.y;
			offset = vectorToCentre.normalized * distanceFromMouse;

			rectTransform.anchoredPosition = inputHandler.screenMousePosition + offset;
		}

		//Populate/Display tooltip
		if (activeAtoms == null)
			return;
		
		hoveredIndex = activeAtoms.hoveredAtom;
		if (hoveredIndex != oldHoveredIndex) {
			if (hoveredIndex == -1) {
				tooltipContainer.SetActive(false);
				goto CLEANUP;
			} else {
				tooltipContainer.SetActive (true);
			}
				
			hoverAtom = activeAtoms [hoveredIndex];

			element.text = hoverAtom.element;
			index.text = hoverAtom.index.ToString();
			amber.text = hoverAtom.amberName;
			pdb.text = hoverAtom.pdbName;
			residue.text = hoverAtom.residueName + hoverAtom.residueNumber.ToString();
			valency.text = hoverAtom.valency.ToString ();

		}

		CLEANUP:
		oldHoveredIndex = hoveredIndex;

	}
}
