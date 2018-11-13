using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class BendBox : MonoBehaviour {

	public Canvas canvas;
	public RectTransform rectTransform;
	public GameObject bendBoxContainer;

	public Text atom0Text;
	public Text atom1Text;
	public Text atom2Text;
	public Text angleText;
	public Text energyText;
	public Text gradientText;
	public Text forceConstantText;


	public Atoms activeAtoms;
	private int a0;
	private int a1;
	private int a2;
	private Atom atom0;
	private Atom atom1;
	private Atom atom2;

	private float[] p0;
	private float[] p1;
	private float[] p2;

	private float[] p0Old;
	private float[] p1Old;
	private float[] p2Old;

	private float[] force0;
	private float[] force1;
	private float[] force2;

	private AngleConnection angleConnection;
	private Bend bend;
	private float[] energies;

	public Vector3 offset;

	void Awake() {
		rectTransform = GetComponent<RectTransform> ();
		offset = new Vector3 (0f, 0f, canvas.planeDistance);

		energies = new float[3];

		p0 = new float[3];
		p1 = new float[3];
		p2 = new float[3];

		p0Old = new float[3];
		p1Old = new float[3];
		p2Old = new float[3];

		force0 = new float[3];
		force1 = new float[3];
		force2 = new float[3];

		float pad = 20f;
		Vector2 initPosition = new Vector2(pad + rectTransform.rect.width * 0.5f, Screen.height - pad - rectTransform.rect.height * 0.5f);
		rectTransform.anchoredPosition = initPosition;

	}

	void Update () {

		//Populate/Display tooltip
		if (activeAtoms == null)
			return;
	

		if (activeAtoms.selection.Count == 3) {
			
			if (bendBoxContainer.activeSelf == false) {
				bendBoxContainer.SetActive(true);
			}
				
			a0 = activeAtoms.selection[0];
			a1 = activeAtoms.selection[1];
			a2 = activeAtoms.selection[2];

			angleConnection = activeAtoms.graph.GetAngleConnection(a0, a1, a2);

			if (angleConnection == null) {
				return;
			}

			atom0 = angleConnection.atom0;
			atom1 = angleConnection.atom1;
			atom2 = angleConnection.atom2;

			bend = activeAtoms.graph.GetBendParameter(angleConnection);

			Mathematics.VectorFromVector3(atom0.p, p0);
			Mathematics.VectorFromVector3(atom1.p, p1);
			Mathematics.VectorFromVector3(atom2.p, p2);

			if (Mathematics.IsEqual(p0, p0Old) && Mathematics.IsEqual(p1, p1Old) && Mathematics.IsEqual(p2, p2Old))
				return;

			Mathematics.GetBendForce(p0, p1, p2, energies, force0, force1, force2, Mathf.Deg2Rad * bend.req, bend.keq);

			atom0Text.text = string.Format("{0}{1} ({2})", atom0.element, atom0.index, atom0.amberName);
			atom1Text.text = string.Format("{0}{1} ({2})", atom1.element, atom1.index, atom1.amberName);
			atom2Text.text = string.Format("{0}{1} ({2})", atom2.element, atom2.index, atom2.amberName);
			angleText.text = string.Format("{0}", angleConnection.angle);
			energyText.text = string.Format("{0}", energies[0]);
			gradientText.text = string.Format("{0}", energies[1]);
			forceConstantText.text = string.Format("{0}", energies[2]);

			Mathematics.Copy(p0, p0Old);
			Mathematics.Copy(p1, p1Old);
			Mathematics.Copy(p2, p2Old);
		} else {
			bendBoxContainer.SetActive(false);
		}
	}
}


