package es.usal.voronto.model.voronoi;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.geom.Area;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Double;
import java.util.ArrayList;
import java.util.Vector;

import java.util.Collections;
import es.usal.voronto.model.voronoi.geometry.Line;
import es.usal.voronto.model.voronoi.geometry.MPolygon;






/**
 * Computes Additively Weighted Power (AWP) voronoi tessellations. 
 * It is adapted from Takashi Ohyama's algorithm (http://www.nirarebakun.com/voro/epwvoro.html)
 * 
 * @author rodri
 *
 */
public class AWPTesselation {
	  int numCells;
	  int height,width;
	  
	  int k,i,j;	//contadores
	  int cont;
	  int m,k2pw,xzpw,xz2pw,yzpw,yz2pw;
	  double eipw,ejpw,diffWeights,cpw,xdistance,ydistance,dipw,yspw;
	  double x0pw,y0pw,yypw,xa1pw,ya1pw;
	  
	  double kxpw[];
	  double kypw[];
	  double kzpw[];
	  
	  double ekpw,diffWeights2,c2pw,a2pw,b2pw,di3pw,ys3pw;
	  double y20pw,y21pw,sa0pw,sa1pw,xxpw,dist_ik,uspw;
	  int br2pw=0;
	private Polygon boundingPolygon;
	private Rectangle boundingBox;
	private int MAX_POINT=1000000000;
	//  private int MAX_POINT=100000;
	  
	  
	  public enum Boundary { TOP, BOTTOM, LEFT, RIGHT, CENTER };
		
	  
	  
	  public AWPTesselation(Cell[] cells, int w, int h)
	  	{
	    width=w;
	    height=h;
	    numCells=cells.length;
	    order=new int[cells.length];
	    for(k=1;k<=numCells;k++)
	    	{
	        xpoint[k-1]=cells[k-1].position[0];
	        ypoint[k-1]=cells[k-1].position[1];
	        weight[k-1]=cells[k-1].weight;
	        order[k-1]=k-1;
	        }
	    
	    //El algoritmo tiene problemas si dos posiciones tienen exactamente la misma y o x (divisiones por cero, etc.)
	    for(k=0;k<numCells;k++)
	    	{
	    	for(j=k+1;j<numCells;j++)
	    		{
	    		if(ypoint[k]==ypoint[j])
	    			ypoint[j]+=2;
	    		if(xpoint[k]==xpoint[j])
	    			xpoint[j]+=2;
	    		}
	    	}
	    
	    this.boundingPolygon=new Polygon();
	    boundingPolygon.addPoint(0,0);
	    boundingPolygon.addPoint(width,0);
	    boundingPolygon.addPoint(width,height);
	    boundingPolygon.addPoint(0,height);
	    }
	  
	  public AWPTesselation(Cell[] cells, Polygon boundingPolygon)
	  	{
		boundingBox=boundingPolygon.getBounds();
	    width=boundingBox.width;
	    height=boundingBox.height;
	   
	    numCells=cells.length;
	    kxpw=new double[numCells];
	    kypw=new double[numCells];
	    kzpw=new double[numCells];
	    xpoint=new double[numCells];
	    ypoint=new double[numCells];
	    weight=new double[numCells];
	    spw=new double[numCells];
	    
	    this.boundingPolygon=new Polygon();
	    order=new int[cells.length];
		for(int i=0;i<boundingPolygon.npoints;i++)
	    	{
	    	this.boundingPolygon.addPoint(boundingPolygon.xpoints[i], boundingPolygon.ypoints[i]);
	    	}
	    this.boundingPolygon.translate(-boundingBox.x, -boundingBox.y);
	    for(k=1;k<=numCells;k++)
	    	{
	        xpoint[k-1]=cells[k-1].position[0]-boundingBox.x;
	        ypoint[k-1]=cells[k-1].position[1]-boundingBox.y;
	   //     if(!boundingPolygon.contains(xpoint[k-1], ypoint[k-1]))
	    //    	System.err.println("Generator outside of bounding polygon");
	        weight[k-1]=cells[k-1].weight;
	        order[k-1]=k-1;
	     //   System.out.println((k-1)+"\t"+cells[k-1].weight);
	        }
	    //El algoritmo tiene problemas si dos posiciones tienen exactamente la misma y o x (divisiones por cero, etc.)
	    for(k=0;k<numCells;k++)
	    	{
	    	for(j=k+1;j<numCells;j++)
	    		{
	    		if(ypoint[k]==ypoint[j])
	    			ypoint[j]+=2;
	    		if(xpoint[k]==xpoint[j])
	    			xpoint[j]+=2;
	    		}
	    	}
	    return;
	    }
	  
	  double xpoint[];
	  double ypoint[];
	  double weight[];
	  double spw[];
	  private int[] order;	//necessary because the method may change the original order of the points
	
	
	  

	  
	  void compute(double weights[],double xpoints[],double ypoints[],int numCells){
	    int i,iipw,jjpw,mmpw;
	    double weight_i,x_i,y_i,c1pw,c2pw,c3pw;
	   // double ox[]=xpoint.clone();
	   // double oy[]=ypoint.clone();
	   // double ow[]=weight.clone();//original values to keep order 
	    
	    for(i=(int)(numCells/2);i>=1;i--)
	      {
	      iipw=i;
	      weight_i=weights[iipw-1];x_i=xpoints[iipw-1];y_i=ypoints[iipw-1];
	      while(2*iipw<=numCells)
	      	{
	        jjpw=2*iipw;
	        if(jjpw+1<=numCells)
	        	{
	        	if(weights[jjpw-1]<weights[jjpw])
	        		{
	        		jjpw++;
	        		}
	        	}
	        if(weights[jjpw-1]<=weight_i)
	        	{
	        	break;
	        	}
	        weights[iipw-1]=weights[jjpw-1];xpoints[iipw-1]=xpoints[jjpw-1];ypoints[iipw-1]=ypoints[jjpw-1];
	        iipw=jjpw;
	      	}//wend
	      weights[iipw-1]=weight_i;xpoints[iipw-1]=x_i;ypoints[iipw-1]=y_i;
	      }//next kk
	    
	    for(mmpw=numCells-1;mmpw>=1;mmpw--){
	      c1pw=weights[mmpw];c2pw=xpoints[mmpw];c3pw=ypoints[mmpw];
	      weights[mmpw]=weights[0];xpoints[mmpw]=xpoints[0];ypoints[mmpw]=ypoints[0];
	      iipw=1;
	      while(2*iipw<=mmpw){
	        i=2*iipw;
	        if(i+1<=mmpw){
	          if(weights[i-1]<=weights[i]){
	            i++;
	          }
	        }
	        if(weights[i-1]<=c1pw){
	          break;
	        }
	        weights[iipw-1]=weights[i-1];xpoints[iipw-1]=xpoints[i-1];ypoints[iipw-1]=ypoints[i-1];
	        iipw=i;
	        
	      }//wend
	      weights[iipw-1]=c1pw;xpoints[iipw-1]=c2pw;ypoints[iipw-1]=c3pw;
	      }//next mm
	    
	  return;
	  }
	
	  private void checkOrder(double ox[], double oy[], double ow[])
	  {	
		   //recompute order
		   for(int k=0;k<numCells;k++)
		   	{
			for(int m=0;m<numCells;m++)
				{
				if(ox[k]==xpoint[m] && oy[k]==ypoint[m] && ow[k]==weight[m])
					{
					order[k]=m;
					break;
					}
				}
		   	}	  
	  }
	  
	  public MPolygon[] computeEdges(){
		  
	//	  long t=System.currentTimeMillis();
		  //-------------------- 1) Do Japanese Magic ----------------------------
		  Vector<Vector<Line2D.Double>> vect=new Vector<Vector<Line2D.Double>>(numCells);
		  for(int i=0;i<numCells;i++)	vect.add(new Vector<Line2D.Double>());
		  Vector<Vector<Line2D.Double>> vect0=new Vector<Vector<Line2D.Double>>(numCells);
		  for(int i=0;i<numCells;i++)	vect0.add(new Vector<Line2D.Double>());
		  Vector<Vector<Line2D.Double>> vectNoSmallLines=new Vector<Vector<Line2D.Double>>(numCells);
		  for(int i=0;i<numCells;i++)	vectNoSmallLines.add(new Vector<Line2D.Double>());
		//  double avgDist=0;
		  
			 double []ox=xpoint.clone();
			 double []oy=ypoint.clone();//compute may resort points, so we keep the originals so we return areas in the same order they were introduced
			 double []ow=weight.clone();
	    compute(weight,xpoint,ypoint,numCells);
	    
	    for(i=1;i<=numCells-1;i++)
	    	{
	    	//System.out.println("ivoro  "+(i-1));
	        for(j=i+1;j<=numCells;j++)
		        {
	 //       	System.out.println(" jvoro "+(j-1));
		    	eipw=Math.pow(xpoint[i-1],2)+Math.pow(ypoint[i-1],2);
		        ejpw=Math.pow(xpoint[j-1],2)+Math.pow(ypoint[j-1],2);
		        diffWeights=weight[j-1]-weight[i-1];
		        cpw=0.5*(ejpw-eipw-diffWeights);		//a kind of weighted distance
		        xdistance=xpoint[j-1]-xpoint[i-1];
		        ydistance=ypoint[j-1]-ypoint[i-1];
		        
		       // if(ydistance!=0)//TODO: comprobaci�n m�a
		        	{
		        	dipw=-xdistance/ydistance;
		        	yspw=cpw/ydistance;
		        	}
		       // else//TODO: cambio m�o
		        //	{
		        	//dipw=1000000000;
		        	//yspw=1000000000;
		        	//dipw=height-1;
		        	//yspw=height-1;
		        	//}
		        if(yspw>0 && yspw<height)
		        	{
		        	x0pw=0;y0pw=yspw;
		        	}
		        else
		        	{
		        	if(dipw>0)
		        		{
		        		x0pw=-yspw/dipw;y0pw=0;
		        		}
		        	else
		        		{
			            x0pw=(height-yspw)/dipw;
			            y0pw=height;
		        		}
		        	}
		        yypw=dipw*width+yspw;
		        if(yypw>0 && yypw<height)
		        	{
		        	xa1pw=width;ya1pw=yypw;
		        	}
		        else
		        	{
			        if(dipw>0)
			          	{
			            xa1pw=(height-yspw)/dipw;ya1pw=height;
			          	}
			        else
			          	{
			            xa1pw=-yspw/dipw;
			            ya1pw=0;
			          	}
		        	}
		        /*if(new Float(xa1pw).isNaN() || new Float(ya1pw).isNaN())
		        	{
		        	System.err.println("NANs!");
		        	//return null;
		        	}*/
		    
		        cont=1;
		        kxpw[cont-1]=x0pw;kypw[cont-1]=y0pw;
		        for(k=1;k<=numCells;k++)	
		        	{
		        	if(k!=i && k!=j)
		        		{
			            ekpw=Math.pow(xpoint[k-1],2)+Math.pow(ypoint[k-1],2);
			            diffWeights2=weight[k-1]-weight[i-1];
			            c2pw=0.5*(ekpw-eipw-diffWeights2);
			            a2pw=xpoint[k-1]-xpoint[i-1];
			            b2pw=ypoint[k-1]-ypoint[i-1];
			            di3pw=-a2pw/b2pw;
			            ys3pw=c2pw/b2pw;
			            y20pw=di3pw*x0pw+ys3pw;
			            y21pw=di3pw*xa1pw+ys3pw;
			            sa0pw=y0pw-y20pw;
			            sa1pw=ya1pw-y21pw;
			            if(sa0pw*sa1pw<=0)
				            {
				            cont++;
				            kxpw[cont-1]=(ys3pw-yspw)/(dipw-di3pw);
				            kypw[cont-1]=dipw*kxpw[cont-1]+yspw;
			            	}
		        		}//if(k!=i && k!=j)
		        	}//next k
		        cont++;
		        kxpw[cont-1]=xa1pw;
		        kypw[cont-1]=ya1pw;
	//	        System.out.println("loop1");
		        for(m=1;m<=cont;m++)
		        	{
		        	kzpw[m-1]=0;
		        	}
	//	        System.out.println("compute");
		        
		        compute(kxpw,kypw,kzpw,cont);
		        
		        for(k=1;k<=cont-1;k++)
		        	  {
	///	        	  System.out.println("  vorok2: "+(k-1));
			          k2pw=k+1;
			          xxpw=(kxpw[k-1]+kxpw[k2pw-1])/2;
			          yypw=dipw*xxpw+yspw;
			          dist_ik=Math.pow(xxpw-xpoint[i-1],2)+Math.pow(yypw-ypoint[i-1],2)-weight[i-1];
			          for(m=1;m<=numCells;m++)
			          	{
	//		        	System.out.println("   loop2 "+upw);
			            if(m!=i && m!=j)
				            {
				            uspw=Math.pow(xxpw-xpoint[m-1],2)+Math.pow(yypw-ypoint[m-1],2)-weight[m-1];
				           // if(yypw==ypoint[m-1])
				           // 	System.out.println("Iguales!");
				            if(uspw<dist_ik)
				            	{
				                br2pw=1;
				                break;
				            	}
				            }
			          	}//next u
			          if(br2pw==0)
			          	{
			        	//this is the edge to add
				        xzpw=(int)(kxpw[k-1]+0.5);
			            xz2pw=(int)(kxpw[k2pw-1]+0.5);
			            yzpw=(int)(kypw[k-1]+0.5);
			            yz2pw=(int)(kypw[k2pw-1]+0.5);
			            
			      //      if( (xzpw==0 && yzpw==298) ||(xz2pw==0 && yz2pw==0))
			      //      	System.err.println("QUe nos vamooos");
			            Point2D.Double p1=new Point2D.Double((double)xzpw,(double)yzpw);
			            Point2D.Double p2=new Point2D.Double((double)xz2pw,(double)yz2pw);
			            //if(p1.distance(p2)>0)	//sometimes it gives puntual lines! and negative points!
					    if(p1.distance(p2)>0 && (p1.x>=0 || p1.y>=0) && (p2.x>=0 || p2.y>=0)	//sometimes it gives puntual lines! and negative points!
					    		)	
				            {
			            	//if(i-1==13 || j-1==13)
			            	//	System.out.println("13");
				            Line2D.Double l=new Line2D.Double(p1,p2);
				            vect.get(i-1).add(l);
				            vect.get(j-1).add(l);
				            vect0.get(i-1).add(l);
				            vect0.get(j-1).add(l);
				            }
			            }//if br2==0
			          else
			          	{//if br2==1
			            br2pw=0;
			          	}
			        }//next k
		        }//next jpw
		    }//next ipw
	    
	   // System.out.println("Time in AWP "+(System.currentTimeMillis()-t));
	   // t=System.currentTimeMillis();
	    
	    
	    
	    checkOrder(ox, oy,ow);
		
	    
	    
	    
	  //Construimos los pol�gonos, teniendo cuidado con el orden de los puntos y con el clipping
	  MPolygon p[]=new MPolygon[numCells];
	  //In the case of a polygon formed just by a line (that cuts the bounding polygon), we cannot decide if the corresponding polygon
	  //is at one side or the other of the line. We save both (clipped) and then compare them with the rest of voronoi regions, choosing
	  //the one without intersections.
	  Area ambiguousAreas[][]=new Area[numCells][2];
	  
	 // boolean noSmallLines=false;
	  
	   //2) ----------------------- Build continuous polylines (one or two)---------------------
	  for(int i=0;i<p.length;i++)
	  	{
		Vector<Line2D.Double> v=null;
		v=vect.get(i);
		
		ArrayList<Point2D.Double> points=new ArrayList<Point2D.Double>();
		ArrayList<Point2D.Double> points2=new ArrayList<Point2D.Double>();
		
		if(v.size()==0)
			{
			//System.err.println("No lines for polygon "+i+"!!");//this usually happens when weight constrains lead to squeezing a polygon from 4-6 points to 4,3 and finally 2, when it becomes a line and does not appear here at all 
			//We either consider this a failed tessellation and return null
			//return null;
			//or continue with a null area...
			}
		else
			{//Form the continuous line
			buildContinuousLine(points, v);
			//If the polygon is formed by two discontinuous polylines, such as "/ /" or "< /", or "< >", we must:
			//a) form the second polyline. 
			//b) Check which edges are closer to join into a single polyline
			while(v.size()>0)
				{
				//a)
				buildContinuousLine(points2,v);
				//b)
				this.joinPolylines(points, points2);
				}
			}//if 1+ lines
		
		//If no polygon lines, add an empty polygon (!) TODO: try to solve above so no empty polygons arrive here
		if(points.size()==0)
			{
			MPolygon p2=new MPolygon(points);

			//We put this polygon on its position on the original order.
			for(int m=0;m<numCells;m++)
				if(order[m]==i)	
					{
					ambiguousAreas[m][0]=new Area(); 
					p[m]=p2;
					}
			}
		else
		//3) -----------------------------Clip polygon if hitting the boundaries -------------------------
			{
			//buildPolygon(points);
			
			//3a) Grow fringe polygons so the clipping gets their whole area
			growAndIntersect(points);
				
			//3b) area clipping
			Path2D.Double poly=new Path2D.Double();
			poly.moveTo(points.get(0).x, points.get(0).y);
			for(int k=1;k<points.size();k++)
				poly.lineTo(points.get(k).x, points.get(k).y);	
			poly.lineTo(points.get(0).x, points.get(0).y);	
			
			
			if(points.size()==2)//In this case the area is only defined by a line, so we cannot get an intersection with the bounding polygon
				{//We need a third point, on the infinite, that will depend on the slope of the line
				Polygon polyA=new Polygon();
				Polygon polyB=new Polygon();
				for(int k=0;k<points.size();k++)
					{
					polyA.addPoint((int)points.get(k).x, (int)points.get(k).y);	
					polyB.addPoint((int)points.get(k).x, (int)points.get(k).y);
					}
				if((points.get(0).x<0 && points.get(0).y<0) || (points.get(0).x>0 && points.get(0).y>0))
					{
					//line with negative slope, we can either add this one...
					polyA.addPoint(MAX_POINT,-MAX_POINT);
					Area aA=new Area(polyA);
					aA.intersect(new Area(boundingPolygon));
					
					//...or this one
					polyB.addPoint(-MAX_POINT,MAX_POINT);
					Area aB=new Area(polyB);
					aB.intersect(new Area(boundingPolygon));
					
					//We will select the one with the smallest area, although that's not necessarily correct, for example if it's a region with very large weight...
					for(int m=0;m<numCells;m++)
						if(order[m]==i)	ambiguousAreas[m][0]=aA;
					for(int m=0;m<numCells;m++)
						if(order[m]==i)	ambiguousAreas[m][1]=aB;
					}
				else
					{
					//line with positive slope, we can either add this one...
					polyA.addPoint(MAX_POINT,MAX_POINT);
					Area aA=new Area(polyA);
					aA.intersect(new Area(boundingPolygon));
					
					//...or this one
					polyB.addPoint(-MAX_POINT,-MAX_POINT);
					Area aB=new Area(polyB);
					aB.intersect(new Area(boundingPolygon));
					
					//we add both and when every polygon is computed we choose the one that has no intersections
					for(int m=0;m<numCells;m++)
						if(order[m]==i)	ambiguousAreas[m][0]=aA;
					for(int m=0;m<numCells;m++)
						if(order[m]==i)	ambiguousAreas[m][1]=aB;
					}
				}//ambiguous polygon
			else
				{
				//Area a0=new Area(poly);
				Area a=new Area(poly);
				a.intersect(new Area(boundingPolygon));
				//if(!a.isSingular())
					//System.err.println(i+": Troceando �reas!!!");
		//		if(a.isEmpty())
		//			System.err.println("�rea vac�a!");
				//else
					{
	//				noSmallLines=false;
					MPolygon p2=buildPolygonFromArea(a);
					//MPolygon p2=new MPolygon(points);
					//MPolygon p2=buildPolygonFromArea(a);
					//We put this polygon on its position on the original order.
					for(int m=0;m<numCells;m++)
						if(order[m]==i)	
							{
							ambiguousAreas[m][0]=a;//for the decision on possible ambiguous polygons later on //TODO 
							p[m]=p2;
							}
					}
				}//unambiguous polygon
				
			}//if polygon has points
		}//for each polygon
	
	  //4) ------------------------ Check possible ambiguous polygons -----------------------
	 /* int kk=0;
	  int sumNull=0;
	  for(i=0;i<numCells;i++)
		 if(p[i]==null)	sumNull++;
	  if(sumNull==numCells)
	  	{
		//we have a  \  / ! / situation: all of the polygons are ambiguous!!
		//1) Search for the smallest area (of the two possibles on each polygon) and assign it
		MPolygon minPolygon=null;
		Area minAreaObject=null;
		double minArea=999999999;
		int minIndex=-1;
		for(i=0;i<numCells;i++)
			{
			for(int j=0;j<2;j++)
				{
				Area area=(Area)ambiguousAreas[i][j].clone();
				MPolygon pol=buildPolygonFromArea(area);
				if(pol.area()<minArea || minPolygon==null)
					{minArea=pol.area(); minPolygon=pol; minIndex=i; minAreaObject=area;}
				}
			}
		p[minIndex]=new MPolygon(minPolygon.getPoints());
		
		minArea=999999999;
		minPolygon=null;
		int minIndex2=-1;
		int numResolved=1;
		
		//2) Search for the ambiguous polygon with the smallest intersection to the minPolygon
		while(numResolved<numCells)
			{
			Area temp=null;
			
			for(i=0;i<numCells;i++)
				{
				if(p[i]==null)
					{
					Area aA=(Area)ambiguousAreas[i][0].clone();
	  				MPolygon pA=buildPolygonFromArea(aA);
	  				
	  				Area aAi=(Area)ambiguousAreas[i][0].clone();
	  				aAi.intersect(minAreaObject);
	  				MPolygon pAi=buildPolygonFromArea(aAi);
	  				if(pAi.area()>0.1 && pA.area()<minArea)		//if there's an intersection and it's the minimum polygon	
	  					{
	  					temp=(Area)aA.clone();
	  					aA.subtract(minAreaObject); 
	  					minPolygon=buildPolygonFromArea(aA);
	  					minIndex2=i;
	  					}
					
	  				Area aB=(Area)ambiguousAreas[i][1].clone();
					MPolygon pB=buildPolygonFromArea(aB);
					
					Area aBi=(Area)ambiguousAreas[i][1].clone();
					aBi.intersect(minAreaObject);
					MPolygon pBi=buildPolygonFromArea(aB);
					if(pBi.area()>0.1 && pB.area()<minArea)	
						{
						temp=(Area)aB.clone();
	  					aB.subtract(minAreaObject); 
						minPolygon=buildPolygonFromArea(aB); 
						minIndex2=i;
						}
					}
				}
			if(minIndex2!=-1)
				{
				p[minIndex2]=new MPolygon(minPolygon.getPoints());
				minAreaObject=(Area)temp.clone();
				numResolved++;
				}
			}
		}
	  else*/
		  {
		  for(i=0;i<numCells;i++)
		  	if(p[i]==null)	
		  		{
		  		boolean isA=true;
		  		for(int j=0;j<numCells;j++)
		  	  		{
		  			if(i!=j && p[j]!=null)	//if not itself or another ambiguous area.
		  				{
		  				Area aA=(Area)ambiguousAreas[i][0].clone();
		  				aA.intersect(ambiguousAreas[j][0]);
		  		//		System.out.println("Intersection of pA with "+j+" has area "+buildPolygonFromArea(aA).area());
		  		  		Area aB=(Area)ambiguousAreas[i][1].clone();
						aB.intersect(ambiguousAreas[j][0]);
				//		System.out.println("Intersection of pB with "+j+" has area "+buildPolygonFromArea(aB).area());
			  	
						if(buildPolygonFromArea(aA).area()>0.1)//We allow a small intersection of 0.1 due to decimal adjustments...
		  					{isA=false; break;}
		  				}
		  			//else sumNull++;
		  	  		}
		  	  //	if(sumNull<numCells-1)
			  	  	{
			  	  	if(!isA)
			  	  		{
			  	  		MPolygon pB=buildPolygonFromArea(ambiguousAreas[i][1]);
			  	  		
			    		
						ambiguousAreas[i][0]=ambiguousAreas[i][1];
			  	  		ambiguousAreas[i][1]=null;
			  	  		p[i]=pB;
			  	  		}
			  	  	else
			  	  		{
			  	  		MPolygon pA=buildPolygonFromArea(ambiguousAreas[i][0]);
			  	  		ambiguousAreas[i][1]=null;
			  	  		p[i]=pA;
			  	  		}
			  	  	}
		  	  	}
		  }
	//  System.out.println("Time in retrieving polygons and clipping "+(System.currentTimeMillis()-t));
	  return p;
	  }
	  
	  private MPolygon buildPolygonFromArea(Area a)
	  	{
		  PathIterator pi=a.getPathIterator(null);
			Vector<Point2D.Float> ants=new Vector<Point2D.Float>(); 
			while(!pi.isDone())
				{
				float[] array=new float[6];
				int type=pi.currentSegment(array);
				if(array[0]<0.00001)	array[0]=0;
				if(array[1]<0.00001)	array[1]=0;
				if((type==PathIterator.SEG_MOVETO || type==PathIterator.SEG_LINETO))
					{
					boolean add=true;
					for(Point2D pant:ants)
						{
						if(pant.equals(new Point2D.Float(array[0],array[1])))	add=false;
						}
					if(add)	ants.add(new Point2D.Float(array[0], array[1]));
					}
				pi.next();
				}
			MPolygon p2=new MPolygon(ants.size());
			for(Point2D pant:ants)
				p2.add((float)pant.getX(), (float)pant.getY());
			return p2;
	  	}

	  
	  private void growAndIntersect(ArrayList<Point2D.Double> points)
	  	{
		if(points.size()==0)	
			{
			//System.err.println("No points to grow"); 
			return;
			}
		Line lp=new Line(points.get(0), points.get(1));
		Line lf=new Line(points.get(points.size()-1), points.get(points.size()-2));
		
		Polygon boundingBoxPolygon=new Polygon();
		boundingBoxPolygon.addPoint(0, 0);
		boundingBoxPolygon.addPoint(0, boundingPolygon.getBounds().height);
		boundingBoxPolygon.addPoint(boundingPolygon.getBounds().width, boundingPolygon.getBounds().height);
		boundingBoxPolygon.addPoint(boundingPolygon.getBounds().width, 0);

		
		for(int m=0;m<boundingBoxPolygon.npoints;m++)
			{
			Line2D.Double lbb;
			if(m!=boundingBoxPolygon.npoints-1)	
				lbb=new Line2D.Double(boundingBoxPolygon.xpoints[m], boundingBoxPolygon.ypoints[m], boundingBoxPolygon.xpoints[m+1], boundingBoxPolygon.ypoints[m+1]);
			else
				lbb=new Line2D.Double(boundingBoxPolygon.xpoints[m], boundingBoxPolygon.ypoints[m], boundingBoxPolygon.xpoints[0], boundingBoxPolygon.ypoints[0]);
			if(Line.linesIntersect(lp.x1, lp.y1, lp.x2, lp.y2, lbb.x1, lbb.y1, lbb.x2, lbb.y2))
				{//intersects with the bounding box, we have to grow the polygon in order to get a proper clipping
				Point2D.Double pr=lp.getProlongationPointFromP1(MAX_POINT);
				points.set(0, new Point2D.Double(pr.x,pr.y));
				}
			if(Line.linesIntersect(lf.x1, lf.y1, lf.x2, lf.y2, lbb.x1, lbb.y1, lbb.x2, lbb.y2))
				{//intersects with the bounding box, we have to grow the polygon in order to get a proper clipping
				Point2D.Double pr=lf.getProlongationPointFromP1(MAX_POINT);
				points.set(points.size()-1, new Point2D.Double(pr.x,pr.y));
				}
			}
			
			//After growing, if it's an open polygon with 3+ lines, we should check if the prolongated lines at the extremes intersect, 
			// and if so, add the intersection point instead of the superextremes (the narrow up end of the U problem, or the / \ problem)
			if(points.size()>3)
				{
				Point2D.Double p0=points.get(0);
				Point2D.Double p1=points.get(1);
				Point2D.Double pf0=points.get(points.size()-2);
				Point2D.Double pf1=points.get(points.size()-1);
				
				if(Line.linesIntersect(p1.x, p1.y, p0.x, p0.y, pf0.x, pf0.y, pf1.x, pf1.y))
					{
					Line l1=new Line(p1, p0);//this one goes backwards
					Line l2=new Line(pf0, pf1);
					Point2D.Double pi=l2.intersection(l1);
					if(pi!=null)//if l1 and l2 are not parallel or the same
						{
						points.set(0, pi);
						points.remove(points.size()-1);
						}
					}
				}
			
		return;
	  	}
	  
	  /**
	   * Joins two lines by the closest end that has an intersection
	   * @param points
	   * @param points2
	   */
	  private void joinPolylines(ArrayList<Point2D.Double> points, ArrayList<Point2D.Double> points2)
	  	{
	    Point2D.Double p1_0=points.get(0);
		Point2D.Double p1_1=points.get(1);
		Point2D.Double p1_f0=points.get(points.size()-1);
		Point2D.Double p1_f1=points.get(points.size()-2);
		
		Point2D.Double p2_0=points2.get(0);
		Point2D.Double p2_1=points2.get(1);
		Point2D.Double p2_f0=points2.get(points2.size()-1);
		Point2D.Double p2_f1=points2.get(points2.size()-2);
		
		//initial vs initial
		double dist=p1_0.distance(p2_0);	
		Line l1=new Line(p1_1, p1_0);
		Line l2=new Line(p2_1, p2_0);
		int pos=0;//edge closer to the other line (l1)
		int pos2=0;//point closer to the other line (l2)
		
		if(p1_0.distance(p2_f0)<dist)
			{
			dist=p1_0.distance(p2_f0);
			l1=new Line(p1_1, p1_0);
			l2=new Line(p2_f1, p2_f0);
			pos=0;
			pos2=points2.size()-1;
			}
		if(p1_f0.distance(p2_0)<dist)
			{
			dist=p1_f0.distance(p2_0);
			l1=new Line(p1_f1, p1_f0);
			l2=new Line(p2_1, p2_0);
			pos=points.size()-1;
			pos2=0;
			}
		if(p1_f0.distance(p2_f0)<dist)
			{
			dist=p1_f0.distance(p2_f0);
			l1=new Line(p1_f1, p1_f0);
			l2=new Line(p2_f1, p2_f0);
			pos=points.size()-1;
			pos2=points2.size()-1;
			}
		
		//need to grow here so it will find an intersection
		Point2D.Double pg11=l1.getProlongationPoint(MAX_POINT);
		l1.x2=pg11.x;
		l1.y2=pg11.y;
		Point2D.Double pg21=l2.getProlongationPoint(MAX_POINT);
		l2.x2=pg21.x;
		l2.y2=pg21.y;
		
		Point2D.Double pi=null;
		if(Line.linesIntersect(l1.x1, l1.y1, l1.x2, l1.y2, l2.x1, l2.y1, l2.x2, l2.y2))
			{
			pi=l2.intersection(l1);
			if(pi!=null)//if l1 and l2 are not parallel or the same
				{
				points.remove(pos);
				points.add(pi);
				points2.remove(pos2);
				if(pos2!=0)
					Collections.reverse(points2);
				for(int i=0;i<points2.size();i++)
					points.add(points2.get(i));
				
				
				}
			else	// The |�| case (it happens!) --> must select as joining point the middle on the infinite
				System.err.println("Insterset but parallel?!");
			}
		else
			{
				pi=l1.getPerpendicularPoint(l1.x2, l1.y2, 100);
				points.set(pos, pi);
				points2.remove(pos2);
				if(pos2!=0)
					Collections.reverse(points2);
				for(int i=0;i<points2.size();i++)
					points.add(pos, points2.get(i));
				//System.err.println("Cannot join polylines");
			}
		return;
		}
	
	
	  /**
	   * Searches for elements to add to the line formed by points (a vector with 2 points), from the array of lines in v
	   */
	  private void buildContinuousLine(ArrayList<Point2D.Double> points, Vector<Line2D.Double> v)
	  	{
		Line2D.Double l=v.remove(0);
		points.add((Point2D.Double)l.getP1());
		points.add((Point2D.Double)l.getP2());
		Point2D.Double	lastPoint=new Point2D.Double(l.x2, l.y2);
		Point2D.Double firstPoint=new Point2D.Double(l.x1, l.y1);
			
		while(v.size()>0)	//keep adding to the line
			{
			int j=0;
			int sizeAnt=v.size();
			for(j=0;j<v.size();j++)
				{
				l=v.get(j);
				if(l.getP1().equals(lastPoint))		//by the end
					{
					v.remove(j);
					lastPoint=(Point2D.Double)l.getP2();
					points.add((Point2D.Double)l.getP2());
					break;
					}
				else if(l.getP2().equals(lastPoint))
					{
					v.remove(j);
					lastPoint=(Point2D.Double)l.getP1();
					points.add((Point2D.Double)l.getP1());
					break;
					}
				else if(l.getP2().equals(firstPoint))	//or by the beginning
					{
					v.remove(j);
					firstPoint=(Point2D.Double)l.getP1();
					points.add(0, (Point2D.Double)l.getP1());
					break;
					}
				else if(l.getP1().equals(firstPoint))	
					{
					v.remove(j);
					firstPoint=(Point2D.Double)l.getP2();
					points.add(0, (Point2D.Double)l.getP2());
					break;
					}
				}
			if(v.size()==sizeAnt)	//we added as much as possible, no more continous lines to add to the polyline
				{
				solveCrosses(points);
				return;
				}
			}
		solveCrosses(points);
		return;
	  	}
	  
	  /**
	   * Checks if the corresponding lines cross with other lines in the polygon, and solves it by removing the point that causes the cross.
	   * (It sometimes happens, the computeEdges is not perfect)
	   * @param points
	   * @return
	   */
	  public void solveCrosses(ArrayList<Point2D.Double> points)
	  	{
		//  ArrayList<Line2D.Double> lines=new ArrayList<Line2D.Double>();
		  if(points.size()<=3)	return;
		
		//Line2D.Double l1=new Line2D.Double(points.get(0), points.get(1));
		//Line2D.Double l2=new Line2D.Double(points.get(points.size()-1), points.get(points.size()-2));
		Line2D.Double lc=new Line2D.Double(points.get(0), points.get(points.size()-1));//The line that closes an open polygon, if existing
		for(int i=1;i<points.size()-2;i++)//the lines that are not last or first should not intersect with the closing line
		 	{
			Line2D.Double lnew=new Line2D.Double(points.get(i), points.get(i+1));
			if(lc.intersectsLine(lnew))
				{
				//System.err.println("Crossed polygon, removing latest point with distance "+points.get(points.size()-2).distance(points.get(points.size()-1)));
				points.remove(points.size()-1);
				return;	
				}
		 	}
		return; 
	  	}
}
