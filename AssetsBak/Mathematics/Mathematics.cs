using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class Mathematics {

	public static int[,] IntDot(int[,] A, int[,] B) {

		int a0 = A.GetLength (0);
		int a1 = A.GetLength (1);
		int b0 = B.GetLength (0);
		int b1 = B.GetLength (1);

		if (a1 != b0)
			Debug.LogErrorFormat ("Incorrect dimensionality for dot product: A({0}, {1}), B({2}, {3})", a0, a1, b0, b1);


		int[,] C = new int[a0,b1];

		for (int i = 0; i < a0; i++) {
			for (int j = 0; j < b1; j++) {
				C [i, j] = 0;
				for (int k = 0; k < a1; k++) {
					C [i, j] = C [i, j] + A [i, k] * B [k, j];
				}
			}
		}
		return C;
	}
}
