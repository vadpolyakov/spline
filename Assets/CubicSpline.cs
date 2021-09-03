using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CubicSpline
{
	#region Fields

	// N-1 spline coefficients for N points
	private float[] a;
	private float[] b;

	// Save the original x and y for Eval
	//private float[] xOrig;
	//private float[] yOrig;
	private List<Vector3> Orig;

	#endregion

	#region Ctor

	/// <summary>
	/// Default ctor.
	/// </summary>
	public CubicSpline()
	{
	}

	/// <summary>
	/// Construct and call Fit.
	/// </summary>
	/// <param name="x">Input. X coordinates to fit.</param>
	/// <param name="y">Input. Y coordinates to fit.</param>
	/// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint.</param>
	/// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint.</param>
	/// <param name="debug">Turn on console output. Default is false.</param>
	public CubicSpline(List<Vector3> points, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
	{
		Fit(points, startSlope, endSlope, debug);
	}

	#endregion

	#region Private Methods

	/// <summary>
	/// Throws if Fit has not been called.
	/// </summary>
	private void CheckAlreadyFitted()
	{
		if (a == null) Debug.LogError("Fit must be called before you can evaluate.");
	}

	private int _lastIndex = 0;

	/// <summary>
	/// Find where in xOrig the specified x falls, by simultaneous traverse.
	/// This allows xs to be less than x[0] and/or greater than x[n-1]. So allows extrapolation.
	/// This keeps state, so requires that x be sorted and xs called in ascending order, and is not multi-thread safe.
	/// </summary>
	private int GetNextXIndex(float x)
	{
		if (x < Orig[_lastIndex].x)
		{
			Debug.LogError("The X values to evaluate must be sorted.");
		}

		while ((_lastIndex < Orig.Count - 2) && (x > Orig[_lastIndex + 1].x))
		{
			_lastIndex++;
		}

		return _lastIndex;
	}

	/// <summary>
	/// Evaluate the specified x value using the specified spline.
	/// </summary>
	/// <param name="x">The x value.</param>
	/// <param name="j">Which spline to use.</param>
	/// <param name="debug">Turn on console output. Default is false.</param>
	/// <returns>The y value.</returns>
	private float EvalSpline(float x, int j, bool debug = false)
	{
		float dx = Orig[j + 1].x - Orig[j].x;
		float t = (x - Orig[j].x) / dx;
		float y = (1 - t) * Orig[j].y + t * Orig[j + 1].y + t * (1 - t) * (a[j] * (1 - t) + b[j] * t); // equation 9
		//if (debug) Console.WriteLine("xs = {0}, j = {1}, t = {2}", x, j, t);
		return y;
	}

	#endregion

	#region Fit*

	/// <summary>
	/// Fit x,y and then eval at points xs and return the corresponding y's.
	/// This does the "natural spline" style for ends.
	/// This can extrapolate off the ends of the splines.
	/// You must provide points in X sort order.
	/// </summary>
	/// <param name="x">Input. X coordinates to fit.</param>
	/// <param name="y">Input. Y coordinates to fit.</param>
	/// <param name="xs">Input. X coordinates to evaluate the fitted curve at.</param>
	/// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint.</param>
	/// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint.</param>
	/// <param name="debug">Turn on console output. Default is false.</param>
	/// <returns>The computed y values for each xs.</returns>
	public float[] FitAndEval(List<Vector3> points, float[] xs, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
	{
		Fit(points, startSlope, endSlope, debug);
		return Eval(xs, debug);
	}

	/// <summary>
	/// Compute spline coefficients for the specified x,y points.
	/// This does the "natural spline" style for ends.
	/// This can extrapolate off the ends of the splines.
	/// You must provide points in X sort order.
	/// </summary>
	/// <param name="x">Input. X coordinates to fit.</param>
	/// <param name="y">Input. Y coordinates to fit.</param>
	/// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint.</param>
	/// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint.</param>
	/// <param name="debug">Turn on console output. Default is false.</param>
	public void Fit(List<Vector3> points, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
	{
		if (float.IsInfinity(startSlope) || float.IsInfinity(endSlope))
		{
			Debug.LogError("startSlope and endSlope cannot be infinity.");
		}

		// Save x and y for eval
		this.Orig = points;

		int n = points.Count;
		float[] r = new float[n]; // the right hand side numbers: wikipedia page overloads b

		TriDiagonalMatrixF m = new TriDiagonalMatrixF(n);
		float dx1, dx2, dy1, dy2;

		// First row is different (equation 16 from the article)
		if (float.IsNaN(startSlope))
		{
			dx1 = points[1].x - points[0].x;
			m.C[0] = 1.0f / dx1;
			m.B[0] = 2.0f * m.C[0];
			r[0] = 3 * (points[1].y - points[0].y) / (dx1 * dx1);
		}
		else
		{
			m.B[0] = 1;
			r[0] = startSlope;
		}

		// Body rows (equation 15 from the article)
		for (int i = 1; i < n - 1; i++)
		{
			dx1 = points[i].x - points[i - 1].x;
			dx2 = points[i + 1].x - points[i].x;

			m.A[i] = 1.0f / dx1;
			m.C[i] = 1.0f / dx2;
			m.B[i] = 2.0f * (m.A[i] + m.C[i]);

			dy1 = points[i].y - points[i - 1].y;
			dy2 = points[i + 1].y - points[i].y;
			r[i] = 3 * (dy1 / (dx1 * dx1) + dy2 / (dx2 * dx2));
		}

		// Last row also different (equation 17 from the article)
		if (float.IsNaN(endSlope))
		{
			dx1 = points[n - 1].x - points[n - 2].x;
			dy1 = points[n - 1].y - points[n - 2].y;
			m.A[n - 1] = 1.0f / dx1;
			m.B[n - 1] = 2.0f * m.A[n - 1];
			r[n - 1] = 3 * (dy1 / (dx1 * dx1));
		}
		else
		{
			m.B[n - 1] = 1;
			r[n - 1] = endSlope;
		}

		//if (debug) Console.WriteLine("Tri-diagonal matrix:\n{0}", m.ToDisplayString(":0.0000", "  "));
		//if (debug) Console.WriteLine("r: {0}", ArrayUtil.ToString<float>(r));

		// k is the solution to the matrix
		float[] k = m.Solve(r);
		//if (debug) Console.WriteLine("k = {0}", ArrayUtil.ToString<float>(k));

		// a and b are each spline's coefficients
		this.a = new float[n - 1];
		this.b = new float[n - 1];

		for (int i = 1; i < n; i++)
		{
			dx1 = points[i].x - points[i - 1].x;
			dy1 = points[i].y - points[i - 1].y;
			a[i - 1] = k[i - 1] * dx1 - dy1; // equation 10 from the article
			b[i - 1] = -k[i] * dx1 + dy1; // equation 11 from the article
		}

		//if (debug) Console.WriteLine("a: {0}", ArrayUtil.ToString<float>(a));
		//if (debug) Console.WriteLine("b: {0}", ArrayUtil.ToString<float>(b));
	}

	#endregion

	#region Eval*

	/// <summary>
	/// Evaluate the spline at the specified x coordinates.
	/// This can extrapolate off the ends of the splines.
	/// You must provide X's in ascending order.
	/// The spline must already be computed before calling this, meaning you must have already called Fit() or FitAndEval().
	/// </summary>
	/// <param name="x">Input. X coordinates to evaluate the fitted curve at.</param>
	/// <param name="debug">Turn on console output. Default is false.</param>
	/// <returns>The computed y values for each x.</returns>
	public float[] Eval(float[] x, bool debug = false)
	{
		CheckAlreadyFitted();

		int n = x.Length;
		float[] y = new float[n];
		_lastIndex = 0; // Reset simultaneous traversal in case there are multiple calls

		for (int i = 0; i < n; i++)
		{
			// Find which spline can be used to compute this x (by simultaneous traverse)
			int j = GetNextXIndex(x[i]);

			// Evaluate using j'th spline
			y[i] = EvalSpline(x[i], j, debug);
		}

		return y;
	}

	/// <summary>
	/// Evaluate (compute) the slope of the spline at the specified x coordinates.
	/// This can extrapolate off the ends of the splines.
	/// You must provide X's in ascending order.
	/// The spline must already be computed before calling this, meaning you must have already called Fit() or FitAndEval().
	/// </summary>
	/// <param name="x">Input. X coordinates to evaluate the fitted curve at.</param>
	/// <param name="debug">Turn on console output. Default is false.</param>
	/// <returns>The computed y values for each x.</returns>
	public float[] EvalSlope(float[] x, bool debug = false)
	{
		CheckAlreadyFitted();

		int n = x.Length;
		float[] qPrime = new float[n];
		_lastIndex = 0; // Reset simultaneous traversal in case there are multiple calls

		for (int i = 0; i < n; i++)
		{
			// Find which spline can be used to compute this x (by simultaneous traverse)
			int j = GetNextXIndex(x[i]);

			// Evaluate using j'th spline
			float dx = Orig[j + 1].x - Orig[j].x;
			float dy = Orig[j + 1].y - Orig[j].y;
			float t = (x[i] - Orig[j].x) / dx;

			// From equation 5 we could also compute q' (qp) which is the slope at this x
			qPrime[i] = dy / dx
				+ (1 - 2 * t) * (a[j] * (1 - t) + b[j] * t) / dx
				+ t * (1 - t) * (b[j] - a[j]) / dx;

			if (debug) Debug.Log(($"[{0}]: xs = {1}, j = {2}, t = {3}", i, x[i], j, t));
		}

		return qPrime;
	}

	#endregion

	#region Static Methods

	/// <summary>
	/// Static all-in-one method to fit the splines and evaluate at X coordinates.
	/// </summary>
	/// <param name="x">Input. X coordinates to fit.</param>
	/// <param name="y">Input. Y coordinates to fit.</param>
	/// <param name="xs">Input. X coordinates to evaluate the fitted curve at.</param>
	/// <param name="startSlope">Optional slope constraint for the first point. Single.NaN means no constraint.</param>
	/// <param name="endSlope">Optional slope constraint for the final point. Single.NaN means no constraint.</param>
	/// <param name="debug">Turn on console output. Default is false.</param>
	/// <returns>The computed y values for each xs.</returns>
	public static float[] Compute(List<Vector3> points, float[] xs, float startSlope = float.NaN, float endSlope = float.NaN, bool debug = false)
	{
		CubicSpline spline = new CubicSpline();
		return spline.FitAndEval(points, xs, startSlope, endSlope, debug);
	}

	/// <summary>
	/// Fit the input x,y points using the parametric approach, so that y does not have to be an explicit
	/// function of x, meaning there does not need to be a single value of y for each x.
	/// </summary>
	/// <param name="x">Input x coordinates.</param>
	/// <param name="y">Input y coordinates.</param>
	/// <param name="nOutputPoints">How many output points to create.</param>
	/// <param name="xs">Output (interpolated) x values.</param>
	/// <param name="ys">Output (interpolated) y values.</param>
	/// <param name="firstDx">Optionally specifies the first point's slope in combination with firstDy. Together they
	/// are a vector describing the direction of the parametric spline of the starting point. The vector does
	/// not need to be normalized. If either is NaN then neither is used.</param>
	/// <param name="firstDy">See description of dx0.</param>
	/// <param name="lastDx">Optionally specifies the last point's slope in combination with lastDy. Together they
	/// are a vector describing the direction of the parametric spline of the last point. The vector does
	/// not need to be normalized. If either is NaN then neither is used.</param>
	/// <param name="lastDy">See description of dxN.</param>
	public static void FitParametric(List<Vector3> points, int nOutputPoints, out float[] xs, out float[] ys,
		float firstDx = float.NaN, float firstDy = float.NaN, float lastDx = float.NaN, float lastDy = float.NaN)
	{
		// Compute distances
		int n = points.Count;
		float[] dists = new float[n]; // cumulative distance
		dists[0] = 0;
		float totalDist = 0;

		for (int i = 1; i < n; i++)
		{
			float dx = points[i].x - points[i - 1].x;
			float dy = points[i].y - points[i - 1].y;
			float dist = Mathf.Sqrt(dx * dx + dy * dy);
			totalDist += dist;
			dists[i] = totalDist;
		}

		// Create 'times' to interpolate to
		float dt = totalDist / (nOutputPoints - 1);
		float[] times = new float[nOutputPoints];
		times[0] = 0;

		for (int i = 1; i < nOutputPoints; i++)
		{
			times[i] = times[i - 1] + dt;
		}

		// Normalize the slopes, if specified
		NormalizeVector(ref firstDx, ref firstDy);
		NormalizeVector(ref lastDx, ref lastDy);

		// Spline fit both x and y to times
		CubicSpline xSpline = new CubicSpline();

		float[] _x = new float[points.Count];
		float[] _y = new float[points.Count];
		for(int i = 0; i < points.Count; i++)
		{
			_x[i] = points[i].x;
			_y[i] = points[i].y;
		}

		List<Vector3> dist_x = new List<Vector3>();
		List<Vector3> dist_y = new List<Vector3>();
		for(int i = 0; i < dists.Length; i++)
		{
			dist_x.Add(new Vector3(dists[i], points[i].x));
			dist_y.Add(new Vector3(dists[i], points[i].y));
		}
		xs = xSpline.FitAndEval(dist_x, times, firstDx / dt, lastDx / dt);

		CubicSpline ySpline = new CubicSpline();
		ys = ySpline.FitAndEval(dist_y, times, firstDy / dt, lastDy / dt);
	}

	private static void NormalizeVector(ref float dx, ref float dy)
	{
		if (!float.IsNaN(dx) && !float.IsNaN(dy))
		{
			float d = Mathf.Sqrt(dx * dx + dy * dy);

			if (d > float.Epsilon) // probably not conservative enough, but catches the (0,0) case at least
			{
				dx = dx / d;
				dy = dy / d;
			}
			else
			{
				Debug.LogError("The input vector is too small to be normalized.");
			}
		}
		else
		{
			// In case one is NaN and not the other
			dx = dy = float.NaN;
		}
	}

	#endregion
}
