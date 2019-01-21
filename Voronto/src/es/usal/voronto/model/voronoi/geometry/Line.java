package es.usal.voronto.model.voronoi.geometry;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;


public class Line extends Line2D.Double{
	
	static final long serialVersionUID=667488;
	double slope;
	
	public Line(double x1, double y1, double x2, double y2)
	{
	super(x1,y1,x2,y2);
	if(x1!=x2)	slope=(y1-y2)/(x1-x2);
	else 
		{
		if(y1<y2) slope=1000000000;
		else	  slope=-1000000000;
		}
	return;
	}
	
	public Line(int x1, int y1, int x2, int y2)
		{
		super(x1,y1,x2,y2);
		if(x1!=x2)	slope=(y1-y2)/(x1-x2);
		else 	
			{
			if(y1<y2) slope=1000000000;
			else	  slope=-1000000000;
			}
		return;
		}
	
	public Line(Point2D.Double p1, Point2D.Double p2)
		{
		super(p1.x,p1.y,p2.x,p2.y);
		if(x1!=x2)	slope=(y1-y2)/(x1-x2);
		else 	
			{
			if(y1<y2) slope=1000000000;
			else	  slope=-1000000000;
			}
		return;
		}
	/**
	 * Returns a line parallel to this one that passed by x,y
	 * @param x
	 * @param y
	 * @return
	 */
	public Line getParallel(double x0, double y0)
		{
		//using point-slope line formula: y-y0=m(x-x0)
		double xf=0;
		double yf=0;
		double m=slope;
		if(x0!=0)	
			{
			xf=0;
			yf=(-m*x0+y0);
			}
		else
			{
			xf=1;
			yf=(m*(1-x0)+y0);
			}
		return new Line(x0,y0,xf,yf);
		}
	
	/**
	 * Gets a line of the same length and position that this line (respect to its initial and ending points),
	 * at a distance d
	 * @return
	 */
	public Line getParallelSegment(double d)
		{
		double angle=Math.atan(slope)+Math.PI/2;
		double x3=x1+Math.cos(angle)*d;
		double y3=y1-Math.sin(angle)*d;//minus because of opposite java coordinate system
		double x4=x2+Math.cos(angle)*d;
		double y4=y2-Math.sin(angle)*d;
		return new Line(x3,y3,x4,y4);
		}
	
	/**
	 * Returns the point that is in the same line than this segment, but at a distance d of x2,y2
	 * @param d
	 * @return
	 */
	public Point2D.Double getProlongationPoint(double d)
		{
		Point2D.Double point=new Point2D.Double();
		point.x=x2;
		point.y=y2;
		
		if(x2>=x1)		
			{
			point.x+=Math.cos(getAngle())*Math.abs(d);
			point.y+=Math.sin(getAngle())*Math.abs(d);
			}
		else
			{
			point.x-=Math.cos(getAngle())*Math.abs(d);
			point.y-=Math.sin(getAngle())*Math.abs(d);
			}
		
		return point;
		}
	
	/**
	 * Same as getProlongationPoint, but from x1,y1
	 * @param d
	 * @return
	 */
	public Point2D.Double getProlongationPointFromP1(double d)
	{
	Point2D.Double point=new Point2D.Double();
	point.x=x1;
	point.y=y1;
	if(x2>=x1)		
		{
		point.x-=Math.cos(getAngle())*Math.abs(d);
		point.y-=Math.sin(getAngle())*Math.abs(d);
		}
	else if(x2<x1)
		{
		point.x+=Math.cos(getAngle())*Math.abs(d);
		point.y+=Math.sin(getAngle())*Math.abs(d);
		}
	/*else	//infinity slope
		{
		if(y1>y2)	point.y+=Math.abs(d);
		}*/

	
	return point;
	}
	
	public double getAngle()
		{
		return Math.atan(slope);
		}
	
	/**
	 * Returns a line perpendicular to this one that passes by x,y
	 * @param x
	 * @param y
	 * @return
	 */
	public Line getPerpendicular(double x0, double y0)
		{
		//using point-slope line formula: y-y0=m(x-x0)
		double xf=0;
		double yf=0;
		double m=-1.0/slope;
		if(slope==0)	m=1000000000;//To avoid infinity
		if(x0!=0)	
			{
			xf=0;
			//yf=(-m*x0+y0);
			yf=-(-m*x0-y0);	//we are in a neg y quadrant in java!
			}
		else
			{
			xf=1;
			//yf=m*(1-x0)+y0;
			yf=-(m*(1-x0)-y0); //we are in a neg y quadrant in java!
			}
		return new Line(x0,y0,xf,yf);
		}
	
	public Line getBisector(double x0, double y0)
		{
		Line l=getPerpendicular(x0,y0);
		Point2D.Double p=intersection(l);
		l.x2=p.x;
		l.y2=p.y;
		return l;
		}
	
	/**
	 * Returns a point that is in a perpendicular line crossing with this line at the start or the end,
	 * and at a determinate distance 
	 * @param start
	 * @param distance
	 * @return
	 */
	public Point2D.Double getPerpendicularPoint(double x0, double y0, double d)
		{
		double angle=(Math.atan(slope)+Math.PI);
		Point2D.Double point=new Point2D.Double();
			if(x2<x1)
			{
			point.x=x0+Math.sin(angle)*d;
			point.y=y0-Math.cos(angle)*d;
			}
		else
			{
			point.x=x0-Math.sin(angle)*d;
			point.y=y0+Math.cos(angle)*d;
			}
		return point;
		}
	/**
	 * Returns the intersection point with line l, or null if lines are parallel or coincident
	 * @param l
	 * @return
	 */
	public Point2D.Double intersection(Line l)
		{
		//line-line intersection point
		//http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
		double x3=l.x1;
		double y3=l.y1;
		double x4=l.x2;
		double y4=l.y2;
		
		double den=(y4-y3)*(x2-x1)-(x4-x3)*(y2-y1);
		double numa=(x4-x3)*(y1-y3)-(y4-y3)*(x1-x3);
		double numb=(x2-x1)*(y1-y3)-(y2-y1)*(x1-x3);
		double ua=numa/den;
		if(numa==0 && numb==0)	
			{
			System.err.println("Parallel/coincident lines");
			return null;
			}
		double xi=x1+ua*(x2-x1);
		double yi=y1+ua*(y2-y1);
		return new Point2D.Double(xi,yi);
		}
	
	// calculates intersection and checks for parallel lines.
		// also checks that the intersection point is actually on
		// the line segment p1-p2
	public static Point2D.Double findIntersection(Point2D.Double p1,Point2D.Double p2,
	  Point2D.Double p3,Point2D.Double p4) {
	  double xD1,yD1,xD2,yD2,xD3,yD3;
	  double dot,deg,len1,len2;
	  double segmentLen1,segmentLen2;
	  double ua,div;
	  //double ub; //unnecessary

	  // calculate differences
	  xD1=p2.x-p1.x;
	  xD2=p4.x-p3.x;
	  yD1=p2.y-p1.y;
	  yD2=p4.y-p3.y;
	  xD3=p1.x-p3.x;
	  yD3=p1.y-p3.y;  

	  // calculate the lengths of the two lines
	  len1=(float)Math.sqrt(xD1*xD1+yD1*yD1);
	  len2=(float)Math.sqrt(xD2*xD2+yD2*yD2);

	  // calculate angle between the two lines.
	  dot=(xD1*xD2+yD1*yD2); // dot product
	  deg=dot/(len1*len2);

	  // if abs(angle)==1 then the lines are parallell,
	  // so no intersection is possible
	  if(Math.abs(deg)==1) return null;

	  // find intersection Pt between two lines
	  Point2D.Double pt=new Point2D.Double(0,0);
	  div=yD2*xD1-xD2*yD1;
	  ua=(xD2*yD3-yD2*xD3)/div;
	 // ub=(xD1*yD3-yD1*xD3)/div;
	  pt.x=(int)(p1.x+ua*xD1);
	  pt.y=(int)(p1.y+ua*yD1);

	  // calculate the combined length of the two segments
	  // between Pt-p1 and Pt-p2
	  xD1=pt.x-p1.x;
	  xD2=pt.x-p2.x;
	  yD1=pt.y-p1.y;
	  yD2=pt.y-p2.y;
	  segmentLen1=(float)(Math.sqrt(xD1*xD1+yD1*yD1)+Math.sqrt(xD2*xD2+yD2*yD2));

	  // calculate the combined length of the two segments
	  // between Pt-p3 and Pt-p4
	  xD1=pt.x-p3.x;
	  xD2=pt.x-p4.x;
	  yD1=pt.y-p3.y;
	  yD2=pt.y-p4.y;
	  segmentLen2=(float)(Math.sqrt(xD1*xD1+yD1*yD1)+Math.sqrt(xD2*xD2+yD2*yD2));

	  // if the lengths of both sets of segments are the same as
	  // the lenghts of the two lines the point is actually
	  // on the line segment.

	  // if the point isn't on the line, return null
	  if(Math.abs(len1-segmentLen1)>0.01 || Math.abs(len2-segmentLen2)>0.01)
	    return null;

	  // return the valid intersection
	  return pt;
	}

	public double getSlope() {
		return slope;
	}
}
