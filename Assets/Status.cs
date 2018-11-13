﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Status : MonoBehaviour {

	public Text statusText;
	public RectTransform rectTransform;
	public GameObject statusContainer;

	void Awake() {
		rectTransform = GetComponent<RectTransform> ();
		statusContainer.SetActive(true);

		float pad = 20f;
		Vector2 initPosition = new Vector2(pad + rectTransform.rect.width * 0.5f, pad + rectTransform.rect.height * 0.5f);
		rectTransform.anchoredPosition = initPosition;
	}

}
