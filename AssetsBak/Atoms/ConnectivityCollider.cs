using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ConnectivityCollider : MonoBehaviour {

	public Atom parent;
	ConnectivityCollider otherCollider;

	int otherIndex;

	double newBondOrder;
	double oldBondOrder;

	void OnTriggerEnter(Collider other) {

		//Wrong type
		otherCollider = other.transform.GetComponent<ConnectivityCollider> ();
		if (otherCollider == null)
			return;

		//Different atom group
		if (!parent.parent == otherCollider.parent.parent)
			return;

		//Other is fully occupied
		if (otherCollider.parent.valency >= otherCollider.parent.maxValency)
			return;

		otherIndex = otherCollider.parent.index;
		newBondOrder = parent.parent.graph.bondTypes.GetBondType (parent.parent, parent.index, otherIndex);
		//oldBondOrder = parent.parent.graph.GetBondOrder (parent.index, otherIndex);

		if (newBondOrder > 0.0) {

			parent.parent.graph.Connect (parent.index, otherIndex, newBondOrder);
		}
	}

}
