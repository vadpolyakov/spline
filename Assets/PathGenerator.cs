using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PathGenerator : MonoBehaviour
{
    public float OBJSPEED = 3f;

    public LineRenderer Trajectory;
    public LineRenderer Snake;
    public Transform OBJ;

    public float max_segment_size = .1f;
    public int segment_count = 3;

    public float oblet_distance = 1f;
    public int oblet_segments_count = 5;

    // Start is called before the first frame update
    void Start()
    {
        DrawSmoothPath();
        
    }

    private void DrawSmoothPath()
    {
        List<Vector3> snake_points = new List<Vector3>();
        for (int i = 0; i < Snake.positionCount; i++)
            snake_points.Add(Snake.GetPosition(i));

        List<Vector3> trajectory_points = new List<Vector3>();
        for (int i = 0; i < Trajectory.positionCount; i++)
            trajectory_points.Add(Trajectory.GetPosition(i));

        LineRenderer new_original = Instantiate(Trajectory);
        new_original.startColor = Color.white;
        new_original.endColor = Color.white;

        List<Vector3> s_snake_points = GetSmoothSnake(snake_points);

        List<Vector3> points = GetSmoothPath(trajectory_points, OBJ.position, snake_points);
        new_original.positionCount = points.Count;
        new_original.SetPositions(points.ToArray());

        StartCoroutine(MoveObj(points));
    }

    IEnumerator MoveObj(List<Vector3> points)
    {
        foreach (Vector3 toPoint in points)
        {
            while (OBJ.position != toPoint)
            {
                OBJ.position = Vector3.MoveTowards(OBJ.position, toPoint, Time.deltaTime * OBJSPEED);
                yield return new WaitForEndOfFrame();
            }
        }
    }

    private List<Vector3> GetSmoothPath(List<Vector3> path_points, Vector3 OBJ_pos, List<Vector3> snake_points)
    {
        List<Vector3> result = new List<Vector3>();
        List<Vector3> buffer = new List<Vector3>() { OBJ_pos };
        buffer.AddRange(path_points);
        buffer.AddRange(GetSmoothPathToSnake(path_points, snake_points[0], snake_points[1]));
        //buffer.AddRange(snake_points);
        result.AddRange(SmoothFunctions.GetCubicSpline(buffer));
        //result.AddRange(GetSmoothPathToSnake(path_points, snake_points[0], snake_points[1]));
        result.AddRange(GetSmoothSnake(snake_points));
        return result;
    }

    private List<Vector3> GetSmoothPathToSnake(List<Vector3> path_points, Vector3 first_snake_pos, Vector3 second_snake_pos)
    {
        List<Vector3> t_to_s = new List<Vector3>();
        Vector3 path_dir = path_points[path_points.Count - 1] - path_points[path_points.Count - 2];
        Vector3 snake_dir = first_snake_pos - second_snake_pos;
        Vector3 snake_to_path_dir = path_points[path_points.Count - 1] - second_snake_pos;
        //t_to_s.Add(path_points[path_points.Count - 1]);
        t_to_s.Add(path_points[path_points.Count - 1] + path_dir * (segment_count * max_segment_size / 10));
        t_to_s.Add(path_points[path_points.Count - 1] + snake_to_path_dir * (segment_count * max_segment_size / 10));
        t_to_s.Add(first_snake_pos + snake_dir * (segment_count * max_segment_size / 10));
        t_to_s.Add(first_snake_pos);
        return t_to_s;
    }

    private List<Vector3> GetSmoothSnake(List<Vector3> snake_points)
    {
        List<Vector3> result = new List<Vector3>() { snake_points[0] };

        
        for(int i = 1; i < snake_points.Count - 1; i += 2)
        {
            float maxsegmentsize = Vector3.Distance(snake_points[i], snake_points[i + 1]) / (oblet_distance * 10);
            result.Add(snake_points[i]);
            Vector3 pd = snake_points[i] - snake_points[i - 1];
            result.AddRange(GetOblet(pd, snake_points[i], snake_points[i + 1], oblet_distance * .1f, oblet_segments_count, maxsegmentsize));
            result.Add(snake_points[i + 1]);
        }
        result.Add(snake_points[snake_points.Count - 1]);
        return result;
    }

    private List<Vector3> GetOblet(Vector3 pd, Vector3 p0, Vector3 p1, float distance, int points_count, float maxsegmentsize)
    {
        Vector3 newp = ((p0 + p1) / 2) + pd * distance;
        List<Vector3> true_points = new List<Vector3>() { p0, newp, p1 };
        return SmoothFunctions.GetSpline(true_points, points_count, maxsegmentsize);
    }
}
