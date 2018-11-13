using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class StretchBox : MonoBehaviour {

	public Canvas canvas;
	public RectTransform rectTransform;
	public GameObject stretchBoxContainer;

	public Text atom0Text;
	public Text atom1Text;
	public Text lengthText;
	public Text energyText;
	public Text gradientText;
	public Text forceConstantText;


	public Atoms activeAtoms;
	private int a0;
	private int a1;
	private Atom atom0;
	private Atom atom1;

	private Connection connection;
	private Stretch stretch;
	private float[] energies;

	public Vector3 offset;

	void Awake() {
		rectTransform = GetComponent<RectTransform> ();
		offset = new Vector3 (0f, 0f, canvas.planeDistance);

		energies = new float[3];

		float pad = 20f;
		Vector2 initPosition = new Vector2(pad + rectTransform.rect.width * 0.5f, Screen.height - pad - rectTransform.rect.height * 0.5f);
		rectTransform.anchoredPosition = initPosition;

	}

	void Update () {

		//Populate/Display tooltip
		if (activeAtoms == null)
			return;

		if (activeAtoms.selection.Count == 2) {
			
			if (stretchBoxContainer.activeSelf == false) {
				stretchBoxContainer.SetActive(true);
			}
				
			a0 = activeAtoms.selection[0];
			a1 = activeAtoms.selection[1];

			connection = activeAtoms.graph.GetConnection(a0, a1);

			if (connection == null) {
				return;
			}

			atom0 = connection.atom0;
			atom1 = connection.atom1;

			stretch = activeAtoms.graph.GetStretchParameter(connection);

			Mathematics.EStretch(connection.length, stretch.keq, stretch.req, energies);

			atom0Text.text = string.Format("{0}{1} ({2})", atom0.element, atom0.index, atom0.amberName);
			atom1Text.text = string.Format("{0}{1} ({2})", atom1.element, atom1.index, atom1.amberName);
			lengthText.text = string.Format("{0}", connection.length);
			energyText.text = string.Format("{0}", energies[0]);
			gradientText.text = string.Format("{0}", energies[1]);
			forceConstantText.text = string.Format("{0}", energies[2]);

		} else {
			stretchBoxContainer.SetActive(false);
		}
	}
}

