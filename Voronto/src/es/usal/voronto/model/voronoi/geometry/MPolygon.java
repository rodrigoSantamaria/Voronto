package es.usal.voronto.model.voronoi.geometry;

import java.awt.Polygon;
import java.awt.geom.Point2D;
import java.util.ArrayList;

import processing.core.*;

public class MPolygon {

	float[][] coords;
	int count;
	
	public MPolygon(){
		this(0);
	}

	public MPolygon(int points){
		coords = new float[points][2];
		count = 0;
	}
	
	public MPolygon(ArrayList<Point2D.Double> points)
		{
		count=points.size();
		coords=new float[count][2];
		for(int i=0;i<count;i++)
			{
			Point2D.Double p=points.get(i);
			coords[i][0]=(float)p.x;
			coords[i][1]=(float)p.y;
			}
		}
	public ArrayList<Point2D.Double> getPoints()
		{
		ArrayList<Point2D.Double> points=new ArrayList<Point2D.Double>();
		for(int i=0;i<count;i++)
			{
			points.add(new Point2D.Double(coords[i][0], coords[i][1]));
			}
		return points;
		}
	
	public MPolygon(Polygon p)
		{
		count=p.npoints;
		coords=new float[count][2];
		for(int i=0;i<p.npoints;i++)	{coords[i][0]=p.xpoints[i];	coords[i][1]=p.ypoints[i];}
		}

	public void translate(float deltaX, float deltaY)
		{
		for(int i=0;i<count;i++)
			{
			coords[i][0]+=deltaX;
			coords[i][1]+=deltaY;
			}
		}
	public void add(float x, float y){
		coords[count][0] = x;
		coords[count++][1] = y;
	}

	public void draw(PApplet p){
		draw(p.g);
	}

	public void draw(PGraphics g){
		g.beginShape();
		for(int i=0; i<count; i++){
			g.vertex(coords[i][0], coords[i][1]);
		}
		g.endShape(PApplet.CLOSE);
	}

	public double area()
		{
		float area=0;
		for(int i = 0; i < count; i++)
			{
			int j = (i + 1) % count;
			area += coords[i][0] * coords[j][1];
			area -= coords[i][1] * coords[j][0];

			//area += (coords[i][0]*coords[i+1][1])-(coords[i+1][0]*coords[i][1]);
			}
		area /= 2;
 
		//if they enter points counterclockwise the area will be negative but correct.
		if(area < 0)
			area *= -1;
 
		return area;
		}
	
	public int count(){
		return count;
	}
	
	
	public float[][] getCoords(){
		return coords;
	}
	
	public Polygon getPolygon()
		{
		int [] x=new int[count];
		int [] y=new int[count];
		for(int j=0;j< count;j++)
			{
			x[j]=(int)coords[j][0];
			y[j]=(int)coords[j][1];
			}
		
		Polygon p=new Polygon(x,y, count);
		return p;
		}
	
	
}