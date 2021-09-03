using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SmoothFunctions : MonoBehaviour
{
    public static List<Vector3> GetSpline(List<Vector3> points, int segment_count, float max_segment_size)
    {
        List<Vector3> result = new List<Vector3>();
        for (int i = 1; i < points.Count - 1; i++)
        {
            List<Vector3> angle = new List<Vector3>();
            angle.AddRange(GetSegments(points[i], points[i - 1], segment_count, max_segment_size));
            //angle.Add(points[i - 1]);
            angle.Reverse();
            //angle.Add(points[i]);
            angle.AddRange(GetSegments(points[i], points[i + 1], segment_count, max_segment_size));
            //angle.Add(points[i + 1]);

            result.AddRange(MakeSmoothCurve(angle));
            
        }
        return result;
    }

    public static List<Vector3> GetCubicSpline(List<Vector3> points)
    {
        float[] x, y;
        CubicSpline.FitParametric(points, 100, out x, out y, 1, 1);
        List<Vector3> result = new List<Vector3>();
        for (int i = 0; i < x.Length; i++)
            result.Add(new Vector3(x[i], y[i]));
        return result;
    }

    public static List<Vector3> MakeSmoothCurve(List<Vector3> arrayToCurve)
    {
        if(arrayToCurve.Count == 0 || arrayToCurve.Count < 2)
            return new List<Vector3>();

        List<Vector3> points;
        List<Vector3> curvedPoints;
        int pointsLength = 0;
        int curvedLength = 0;

        pointsLength = arrayToCurve.Count;

        curvedLength = pointsLength * 3 - 1;
        curvedPoints = new List<Vector3>(curvedLength);

        float t = 0.0f;
        for (int pointInTimeOnCurve = 0; pointInTimeOnCurve < curvedLength + 1; pointInTimeOnCurve++)
        {
            t = Mathf.InverseLerp(0, curvedLength, pointInTimeOnCurve);

            points = new List<Vector3>(arrayToCurve);

            for (int j = pointsLength - 1; j > 0; j--)
            {
                for (int i = 0; i < j; i++)
                {
                    points[i] = (1 - t) * points[i] + t * points[i + 1];
                }
            }

            curvedPoints.Add(points[0]);
        }

        return curvedPoints;
    }

    public static List<Vector3> GetSegments(Vector3 from, Vector3 to, int segment_count, float max_segment_size)
    {
        List<Vector3> result = new List<Vector3>();
        for (int i = 0; i < segment_count; i++)
        {
            Vector3 new_segment;
            if (i == 0)
            {
                if (Vector3.Distance(from, to) < max_segment_size)
                    break;
                new_segment = Vector3.MoveTowards(from, to, max_segment_size);
            }
            else
            {
                if (Vector3.Distance(result[i - 1], to) < max_segment_size)
                    break;
                new_segment = Vector3.MoveTowards(result[i - 1], to, max_segment_size);
            }
            result.Add(new_segment);
        }
        return result;
    }
}
