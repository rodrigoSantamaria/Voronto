package es.usal.voronto.model.voronoi;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Rectangle2D.Float;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;
import java.util.TreeMap;

import es.usal.voronto.control.Voronto;
import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.ontology.OntologyTerm;
import es.usal.voronto.model.voronoi.geometry.MPolygon;
import es.usal.voronto.view.VoronoiVisualization;

public class BernhardtTessellation {

	float[][] edges;//edges of the regions (hard to map to regions!)
	int width;
	int height;
	public double scale;
	private Cell[] currentCells;//current hierarchical cell structure (for computations and drawings)
	private Cell[] cells;		//total hierarchical cell structure
	int numLeaves=0;
	
	public float unit;
	private double maxDistortion=0.01;
	private boolean anyoneNeedsReweight=true;
	private Polygon boundingPolygon;
	private int contawp;
	private String currentName="root";
	
	public int maxLevel=0;//max depth ontology level found (or maxDepth if maxLevel>maxDepth)
	public int maxDepth=3;//max depth ontology level determined 
	
	private ExpressionData expData;
	
	private boolean tessellationFailed=false;
	private boolean zeroAreas=false;
	private boolean allZeroAreas=false;
	
	private int mode=0; //0-normal execution, 1-debug mode
	private int ontology;
	
	private boolean save=true;//save the mapping 
	private String mapFile="map.txt";
	private BufferedWriter out=null;
	private BufferedWriter bw;
	private int minAnnotations=0;
	
	/**
	 * 
	 * Generates a voronoi diagram for non hierarchical cells
	 * @param cells
	 * @param width
	 * @param height
	 */
	public BernhardtTessellation(Cell[] cells, int width, int height)
		{	
		this.cells=cells;
		this.width=width;
		this.height=height;
		boundingPolygon=getBoundingBox();
		computeVoronoi();
		}
	
	

	/**
	 * TODO: Generates a voronoi diagram from an ontology hierarchy
	 * @param map
	 * @param width
	 * @param height
	 */
	public BernhardtTessellation(TreeMap<OntologyTerm,TreeMap> m, ExpressionData md, int width, int height, int ontology, int maxDepth) throws Exception
		{
		this.expData=md;
		this.ontology=ontology;
		this.maxDepth=maxDepth;
		this.minAnnotations=5;
		if(ontology==VoronoiVisualization.KEGG || ontology==VoronoiVisualization.CUSTOM)
			minAnnotations=0;
		
	//	bw=new BufferedWriter(new FileWriter("taxons.txt"));
		
		//1) Recursively generate the cell hierarchy
		long time=System.currentTimeMillis();
		OntologyTerm[] terms=m.keySet().toArray(new OntologyTerm[0]);//TODO: this might lead to unnefficient computing.
		ArrayList<Cell> cellsTemp=new ArrayList<Cell>();
		numLeaves=0;
		if(save)	out=new BufferedWriter(new FileWriter(mapFile));
	    for(int i=0;i<terms.length;i++)	//Covering leaf children (only for major nodes)
	    	{
	    	Cell c=generateRecursiveCell(m.get(terms[i]), terms[i],1);
	    	if(c!=null)	cellsTemp.add(c);
	    	}
	    if(save)	out.close();;
	   // bw.close();
		
	    cells=cellsTemp.toArray(new Cell[0]);
	    if(cells.length==0)	throw new Exception("No mapped genes. Check that gene ids match with ontology's available ids");
	    System.out.println("Number of cells "+numLeaves);
	    System.out.println("Time in building cells "+(System.currentTimeMillis()-time)/1000.0);
		
	    //2) Sort top terms by weight (lower terms have been sorted during the recursive generation)
	    time=System.currentTimeMillis();

    	Comparator<Cell> byWeight=new Voronto.WeightComparator();
		Arrays.sort(cells, byWeight);
		
	    System.out.println("Time in sorting cells "+(System.currentTimeMillis()-time)/1000.0);
		
		//3) Set initial polygon width and height
		this.width=width;
		this.height=height;
		
		//4) Compute places recursively
		this.boundingPolygon=getBoundingBox();
		this.currentCells=cells;
		
		long t=System.currentTimeMillis();
		recursiveComputeVoronoi(cells, 0, getBoundingBox());
		System.out.println("Time in computing voronoi "+ (System.currentTimeMillis()-t)/1000.0+" seconds");
		
		return;
		}
	
	private Polygon getBoundingBox()
		{
		Polygon p=new Polygon();
		p.addPoint(0, 0);
		p.addPoint(width, 0);
		p.addPoint(width, height);
		p.addPoint(0, height);
		return p;
		}
	
	
	public Cell buildCell(OntologyTerm term, int level)
		{
		Cell c=new Cell(1, term, level);
		switch(ontology)
			{
			case VoronoiVisualization.KEGG:
				c.numLeaves=c.getNumMatchedKOTerms(expData);
				break;
			case VoronoiVisualization.GO:
			case VoronoiVisualization.BP:
			case VoronoiVisualization.CC:
			case VoronoiVisualization.MF:
			case VoronoiVisualization.SLIMBP:
			case VoronoiVisualization.SLIMCC:
				c.numLeaves=c.getNumMatchedGOTerms(expData);
				break;
			case VoronoiVisualization.REACTOME:
				c.numLeaves=c.getNumMatchedTerms(expData);
				break;
			default:
				if(expData!=null)	c.numLeaves=c.getNumMatchedTerms(expData);
				break;
			}
		c.weight=c.numLeaves;	//In general, the term has a size proportional to the number of genes annotated with it
		numLeaves++;
		return c;
		}
	
	public Cell generateRecursiveCell(Map<OntologyTerm, Map> m, OntologyTerm term, int level)
	{
	OntologyTerm[] names=null;
	if(m!=null)	names=m.keySet().toArray(new OntologyTerm[0]);
	else		return null;
	
	if(names==null || names.length==0)	
		{
		if(expData!=null)
			{
		//	try{bw.write("\t\t"+term.id+"\t\t"+term.name);bw.newLine();}catch(Exception e){e.printStackTrace();}
			
			Cell c=buildCell(term, level);
			if(c.weight<=minAnnotations)	return null; //no children and few annotations, so not shown
			else			
				{
				if(save)	try{for(String g:term.geneIds)	{out.write(term.id+"\t"+term.name+"\t"+g);out.newLine();}}catch(Exception e){e.printStackTrace();}
				return c;
				}
			}
		else 
			{
			numLeaves++;
			
			if(mode==1)	System.out.println("Creating cell "+term.name+" at level "+level);
			
			if(term.geneIds==null || term.geneIds.size()==0)	return new Cell(1, term, level);
			else											
				{
			//	if(save)	try{for(String g:term.geneIds)	{out.write(term.id+"\t"+term.name+"\t"+g);out.newLine();}}catch(Exception e){e.printStackTrace();}
				return new Cell(term.geneIds.size(), term, level);
				}
			}
		}
	else	//has children
		{
		if(level>maxLevel)
			{
			maxLevel=Math.min(maxDepth, level);
			}
		ArrayList<Cell> cs=new ArrayList<Cell>();
		for(OntologyTerm n:names)
			{
			Cell sc=generateRecursiveCell((Map<OntologyTerm,Map>)(m.get(n)), n, level+1);
			if(sc!=null)
				cs.add(sc);
			}
		
		Cell c=null;
		if(expData!=null)
			{
			c=buildCell(term, level);
			c.subcells=cs.toArray(new Cell[0]);
			}
		else
			c=new Cell(term, level, cs);
		//if(c!=null && c.weight>0)
		if(c!=null && c.weight>minAnnotations)
			{
		//	if(save)	try{for(String g:term.geneIds)	{out.write(term.id+"\t"+term.name+"\t"+g);out.newLine();}}catch(Exception e){e.printStackTrace();}
			return c;
			}
		else
			return null;
		}
	}

	public void recursiveComputeVoronoi(Cell[] c, int level, Polygon totalArea)
	{
//	System.out.println("Computing voronoi for level "+level);
	//if(level<1)	//we only compute voronoi tessellation till maxDepth (in the case of GO complete, with 14 levels, this is adviceable)
	if(level<maxDepth)	//we only compute voronoi tessellation till maxDepth (in the case of GO complete, with 14 levels, this is adviceable)
		{
		currentCells=c;
	    boundingPolygon=totalArea;
		contawp=0;
		
		if(totalArea.npoints>2)
			{
			if(c!=null && c.length>1)
			//if(c!=null && c.length>0)
				computeVoronoi();
			else if(c!=null && c.length==1 && level==0)	//In the rare but possible case of having a root term with only 1 child
				c[0].region=new MPolygon(boundingPolygon);
			for(Cell cc:c)
				{
				if(cc.subcells!=null && cc.subcells.length>1)	
					{
					if(mode==1)	
						System.out.println("Computing voronoi for subcells of "+cc.term.name+" (level "+cc.level+")");
					currentName=cc.term.name;
					recursiveComputeVoronoi(cc.subcells, level+1, cc.region.getPolygon());
					}
				else if(cc.subcells!=null && cc.subcells.length==1)//TODO: a term with only one children that has a term with only one children, etc. is not reflected (only 1st level)
					{
					cc.subcells[0].region=new MPolygon(cc.region.getPoints());
					cc.subcells[0].centroid=cc.centroid;
					cc.subcells[0].position=cc.position;
					recursiveComputeVoronoi(cc.subcells, level+1, cc.region.getPolygon());
					}
				}
			}
		else
			System.err.println(currentName+"at level "+level+" has no area");
		}
	}

	
	public void computeVoronoi()
		{
		if(currentCells==null || currentCells.length==0 ){System.err.println("No cells defined"); return;}
		
		Cell[] cellsAnt=new Cell[currentCells.length];
		for(int i=0;i<cellsAnt.length;i++)		cellsAnt[i]=new Cell(currentCells[i]);	//deep copy
			
		if(mode==1)		
			if(currentName.contains("Transport and Catabolism"))	
				System.out.println("tal");
				
		anyoneNeedsReweight=true;
		//0) Initial placement: trying to optimize space based on weights
		if(mode==1)	System.out.println("0) Initial placement for "+currentName);
		
		//gridPlacement();
		//weightedGridPlacement();
		weightedPolygonPlacement();
		
		if(mode==1)
			for(Cell c:currentCells)
				System.out.println(c.term.name+" at "+c.position[0]+", "+c.position[1]+" with weight "+c.weight+" and numLeaves "+c.numLeaves);
		
		for(int kk=0;kk<15000;kk++)
			{
			
			//1) Compute voronoi regions from initial positions: TODO: this is the point where it should be weighted voronoi regions!
			if(mode==1)	System.out.println("1) AWPVT");
			computeRegionsAWP();
		//TODO: might return with a failed tessellation, no region built!
		
	    
		  //3) Estimate error between desired area and obtained area
	      if(kk==0)//first iteration, we must also compute and save desired areas
			  {
		      float sum=0;
		      for(Cell c:currentCells)	      	sum+=c.numLeaves;
		      for(Cell c:currentCells)
		        {
		      	c.desiredArea=new MPolygon(boundingPolygon).area()*c.numLeaves/sum;
		 		}
			  }

		  	float asum=0;
		  	for(Cell c:currentCells)
				asum+=c.region.area();
		  	if(mode==1)
			  	System.out.println("Regions sum area: "+asum +"\t vs \t Bounding polygon area: "+new MPolygon(boundingPolygon).area());
		  	
			boolean unmatchingAreas=false;
		    if(new MPolygon(boundingPolygon).area()>asum+1000 || new MPolygon(boundingPolygon).area()<asum-1000)
		    	{
		    	if(mode==1)	System.err.println(currentName+": Areas do not match");
		    	unmatchingAreas=true;
		    	}
	      	
	      if(!anyoneNeedsReweight)
			  break;

	      //2) Reweight if the area is very different from the desired area
	      if(mode==1)	System.out.println("2) Reweight");
		  float ratioError=computeAreaRatios();
		  
		  if(mode==1)	 System.out.println(tessellationFailed+"\t"+ratioError+"\t"+computeRatioError(cellsAnt));
		  float errIncrease=ratioError/computeRatioError(cellsAnt);
		  
		  boolean tooMuchIncrease=errIncrease>2;
		  if(zeroAreas)	tooMuchIncrease=errIncrease>2.5;
		  
		  if(allZeroAreas)
			  System.err.println("All areas are zero in this term!!!");
		  if(kk>0 && (unmatchingAreas || tessellationFailed || tooMuchIncrease || allZeroAreas))
		  	  {
			  if(mode==1) 
				  {
				  System.out.println("unmatching areas: "+unmatchingAreas+"\ttessellationFailed"+tessellationFailed+"\ttooMuchIncrease"+tooMuchIncrease+"\tallZeroAreas: "+allZeroAreas);
				  System.out.println(currentName+": Stop & Rollback");
				  }
	    	  for(int i=0;i<cellsAnt.length;i++)		currentCells[i]=new Cell(cellsAnt[i]);	//deep copy
	    	  tessellationFailed=false;
	    	  break;
	    	  }
	      else
	      	{
	    	cellsAnt=new Cell[currentCells.length];
			for(int i=0;i<cellsAnt.length;i++)		cellsAnt[i]=new Cell(currentCells[i]);	//deep copy
	      	}
	      
	 		
	      reweighting();

	      //3) Centroidal voronoi tesselation (modifies initial positions in order to make region shapes more regular)
	      if(mode==1)	System.out.println("3) CVT");
	      lloydMethod(0.5,0);
	     // lloydMethod(0.1,0);
	     // lloydMethod(0.01,0);
			}//iterative thingy
			
		//F) Latest centroid computation  
		computeCentroids();
		
		if(mode==1)	System.out.println("Finishing");
		}
	

	
	public void computeRegionsAWP()
		{
		Rectangle bp=boundingPolygon.getBounds();
		if(mode==1)		System.out.println("AWP"+contawp);
		if(mode==1)		if(contawp==130)
			System.out.println("awp");
		contawp++;
		
		zeroAreas=false;
		allZeroAreas=false;
		int numZeros=0;
		//if(currentName.equals("Transport and Catabolism"))
		//	System.out.println("Here"+contawp);
		AWPTesselation v=new AWPTesselation(currentCells, boundingPolygon);
		
		MPolygon[] p=v.computeEdges();
		if(p==null)
			{
			//System.err.println(currentName+": Tessellation failed");
			tessellationFailed=true;
			return;
			}
		for(MPolygon pp:p)
			if(pp==null)
				{
				p=null;
				tessellationFailed=true;
				return;
				}
		for(MPolygon pp:p)
			{
			if(pp.getPoints()==null || pp.getPoints().size()==0)
				{
				zeroAreas=true; //return;
				numZeros++;
				}
			}
		if(numZeros==p.length)	allZeroAreas=true;
			
		for(int i=0;i<p.length;i++)
			{
			p[i].translate(bp.x, bp.y);
			currentCells[i].region=p[i];
			}
		return;
		}
	
	/**
	 * Method to compute a centroid voronoi tessellation from the given voronoi points
	 */

	public void lloydMethod(double d, int it)
		{
		//1) Compute their centroids, which are the new generators
		computeCentroids();
		
		//2) Finish if convergence acquired, compute recursively if not
		float error=avgDistCentroidToGenerator();
		if(error<d || it>20)		
		//if(error<d || it>200)		
			{
//			System.out.println("Terminado, error: "+error);	
			return;
			}
		else 					
			{
			for(Cell c:currentCells)
				{
				if(this.boundingPolygon.contains(c.centroid[0], c.centroid[1]))
					c.position=c.centroid;
				}
			computeRegionsAWP();

			lloydMethod(d, it+1);
			}
		}

/**
 * Computes area ratios as described by Bernhardt et al. (different to Balzer/Deussen)
 * It is basically a measure of how different are the desired and the obtained area
 * Return the mean error on the ratios
 */
public float computeAreaRatios()
	{
	float sum=0;
	for(Cell c:currentCells)
		{	
		if(c.region.area()>0)		
			c.areaRatio=c.desiredArea/c.region.area();
		else						//c.areaRatio=c.desiredArea;
			{
			//System.err.println("Regrowing "+c.term.name);
			c.areaRatio=5.0;//In the case some area is 0, we give it a chance to grow
			//c.areaRatio=c.desiredArea;
			}
		sum+=Math.abs(1-c.areaRatio);
			
		if(mode==1)	System.out.println(c.term.name+"\tdesired: "+c.desiredArea+", actual: "+c.region.area()+"\tratio:"+c.areaRatio+" weight:"+c.weight+" leaves:"+c.numLeaves);
		
		}
	return sum/currentCells.length;
	}

/**
 * Computes the ratio error average of cells
 * 
 * @param cells
 * @return
 */
public float computeRatioError(Cell[] cells)
	{
	float err=0;
	for(Cell c:cells)
		{
		err+=Math.abs(1-c.areaRatio);
		}
	return err/=cells.length;
	}

public void reweighting()
	{
	anyoneNeedsReweight=false;
	for(Cell c:currentCells)
		{
		if((c.areaRatio-1) > maxDistortion)
			{
			c.weight=(float)(c.weight*(1+c.convergenceFactor*(c.desiredArea-c.region.area())/c.desiredArea));
			anyoneNeedsReweight=true;
			//c.convergenceFactor-=0.1;	//NOTE: not used so far
			}
		}
	}

	public Polygon getPolygon(MPolygon mp)
		{
		Polygon p=new Polygon();
		for(int i=0;i<mp.getCoords().length;i++)
			p.addPoint((int)mp.getCoords()[i][0], (int)mp.getCoords()[i][1]);
		return p;
		}

	public Cell[] getCells()
		{
		return cells;
		}

	public void computeCentroids()
		{
		for(int i=0;i<currentCells.length;i++)
		 	{
			currentCells[i].centroid=new float[2];
			currentCells[i].centroid[0]=0;
			currentCells[i].centroid[1]=0;
			float[][] cc=currentCells[i].region.getCoords();
			for(int j=0;j<cc.length;j++)
				{
				currentCells[i].centroid[0]+=cc[j][0];
				currentCells[i].centroid[1]+=cc[j][1];
				};
			currentCells[i].centroid[0]/=cc.length;
			currentCells[i].centroid[1]/=cc.length;
			}
		}

	protected boolean isEdgeShared(int face1[], int face2[]){
		for(int i = 0; i < face1.length; i++){
			int cur = face1[i];
			int next = face1[(i + 1) % face1.length];
			for(int j = 0; j < face2.length; j++){
				int from = face2[j];
				int to = face2[(j + 1) % face2.length];
				if(cur == from && next == to || cur == to && next == from)
					return true;
			}
		}
		return false;
	}

	public float avgDistCentroidToGenerator()
		{
		float dist=0;
		if(currentCells==null || currentCells.length==0)
			{
			System.err.println(currentName+": avgDistCentroidToGenerator: no generators or length different to centroids");
			return -1;
			}
		for(int i=0;i<currentCells.length;i++)
			{
			dist+=Math.sqrt((currentCells[i].position[0]-currentCells[i].centroid[0])*(currentCells[i].position[0]-currentCells[i].centroid[0])+(currentCells[i].position[1]-currentCells[i].centroid[1])*(currentCells[i].position[1]-currentCells[i].centroid[1]));
			}
		return dist/currentCells.length;
		}
	
	/**
	 * We make a grid and place each cell on a point that is 
	 * 1) Inside the bounding polygon
	 * 2) 
	 */
	public void gridPlacement()
		{
		double unitWidth=boundingPolygon.getBounds().width/Math.sqrt(currentCells.length);
		double unitHeight=boundingPolygon.getBounds().height/Math.sqrt(currentCells.length);
		int x0=boundingPolygon.getBounds().x;
		int y0=boundingPolygon.getBounds().y;
		int x=x0;
		int y=y0;
		boolean fitting=false;
		while(!fitting)
			{
			fitting=true;
			x=x0;
			y=y0;
			for(int i=0;i<currentCells.length;i++)	//For each cell
				{
				currentCells[i].position[0]=-1000000;
				currentCells[i].position[1]=-1000000;
				
				while(!boundingPolygon.contains(currentCells[i].position[0], currentCells[i].position[1]))//If it's not inside the polygon, go to next cross in the grid
					{
					if(x+unitWidth*0.5>x0+boundingPolygon.getBounds().width)	//If we reach the end of a row, go to next
						{
						y+=unitHeight;
						if(y>y0+boundingPolygon.getBounds().height-unitHeight*0.5)	//If we reach the last cross in the grid, reduce scale and retry
							{
							unitWidth=0.9*unitWidth;
							unitHeight=0.9*unitHeight;
							fitting=false; 
							break;
							}
						x=x0;
						}
					currentCells[i].position[0]=(float)(x+unitWidth*0.5);
					currentCells[i].position[1]=(float)(y+unitHeight*0.5);
					x+=unitWidth;
					}
				}
			}
		return;
		}
	
	//As grid placement, but grid is divided in as many points as the sum of the weights (not the sum of the cells), and 
	//places are assigned accordingly
	public void weightedGridPlacement()
	{
	int sumWeight=0;
	for(Cell c:currentCells)	sumWeight+=c.weight;
	
	double unitWidth=boundingPolygon.getBounds().width/Math.sqrt(sumWeight);
	double unitHeight=boundingPolygon.getBounds().height/Math.sqrt(sumWeight);
	int x0=boundingPolygon.getBounds().x;
	int y0=boundingPolygon.getBounds().y;
	int x=x0;
	int y=y0;
	boolean fitting=false;
	while(!fitting)
		{
		fitting=true;
		x=x0;
		y=y0;
		for(int i=0;i<currentCells.length;i++)	//For each cell
			{
			if(mode==1)	if(i>0)	System.out.println("cell "+(i-1)+"\t"+currentCells[i-1].weight+"\tposition "+currentCells[i-1].position[0]+", "+currentCells[i-1].position[1]);
			currentCells[i].position[0]=-1000000;
			currentCells[i].position[1]=-1000000;
			float w=currentCells[i].weight;
				
			while(w>0 || !boundingPolygon.contains(currentCells[i].position[0], currentCells[i].position[1]))//If it's not inside the polygon, go to next cross in the grid 
				{																						
				if(x+unitWidth*0.5>x0+boundingPolygon.getBounds().width)	//If we reach the end of a row, go to next
					{
					y+=unitHeight;
					if(y+unitHeight*0.5>y0+boundingPolygon.getBounds().height)	//If we reach the last cross in the grid, reduce scale and retry
						{
						if(mode==1)	System.out.println("cell "+i+" need to refactor, unitW and H"+unitWidth+", "+unitHeight);
						unitWidth=0.99*unitWidth;
						unitHeight=0.99*unitHeight;
						fitting=false; 
						break;
						}
					x=x0;
					}
				currentCells[i].position[0]=(float)(x+unitWidth*0.5);
				currentCells[i].position[1]=(float)(y+unitHeight*0.5);
				x+=unitWidth;
				w--;
				}
			 if(!boundingPolygon.contains(currentCells[i].position[0], currentCells[i].position[1]))//if the position falls outside the polygon
			 	{
				unitWidth=0.99*unitWidth;
				unitHeight=0.99*unitHeight;
				fitting=false; 
				break; 
			 	}
			}//for each cell
		}
	return;
	}

	public void weightedPolygonPlacement()
		{
		int sumWeight=0;
		for(Cell c:currentCells)	sumWeight+=c.weight;
		double unitWidth=boundingPolygon.getBounds().width/Math.sqrt(sumWeight);
		double unitHeight=boundingPolygon.getBounds().height/Math.sqrt(sumWeight);
		ArrayList<Rectangle2D.Float> polygons=new ArrayList<Rectangle2D.Float>(); //polygons reserved for each point
				
		if(unitWidth==0 || unitHeight==0)
			{
			//System.err.println("Empty starting area");
			return;
			}
		int x0=boundingPolygon.getBounds().x;
		int y0=boundingPolygon.getBounds().y;
		boolean fitting=false;
		ArrayList<Placement> inits=new ArrayList<Placement>();
		
		if(mode==1)	
			System.out.println("Placement for subterms of "+this.currentName);
		
		float sumArea=0;
		
		while(!fitting)
			{
			inits.clear();
			polygons.clear();
			//inits.add(new Point2D.Float(x0,y0));
			for(int j=0;j<boundingPolygon.npoints;j++)
				inits.addAll(getAllPlacements(new Point2D.Float(boundingPolygon.xpoints[j],boundingPolygon.ypoints[j])));
			
			fitting=true;
			for(Cell c: currentCells)
				sumArea+=Math.sqrt(c.weight)*unitWidth*Math.sqrt(c.weight)*unitHeight;
			if(mode==1)
				{
				System.out.println("Unit width and height: "+unitWidth+", "+unitHeight);
				System.out.println("Bounds: "+boundingPolygon.getBounds().width+", "+boundingPolygon.getBounds().height);
				System.out.println("Total area is "+boundingPolygon.getBounds().width*boundingPolygon.getBounds().height);
				System.out.println("Total polygonal area is "+new MPolygon(boundingPolygon).area());
				System.out.println("Total subcells area is "+sumArea);
				if(Double.isNaN(sumArea))
					System.out.println("sumArea is Nan!!!");
				}
			sumArea=0;
		
			//for(int i=0;i<currentCells.length;i++)	//For each cell
			for(int i=currentCells.length-1; i>=0;i--)
				{
				currentCells[i].position[0]=-1000000;
				currentCells[i].position[1]=-1000000;
				float w=currentCells[i].weight;
				Point2D.Float center=new Point2D.Float(-1000000,-1000000);
				//if(mode==1)					System.out.println("Placing "+i+" with weight "+w);
				
				//int cont=0;
				while(inits.size()>0)
					{
					Placement pinit=inits.remove(0);
					Point2D.Float init=pinit.point;
					
					float incx=(float)(Math.sqrt(w)*unitWidth);
					float incy=(float)(Math.sqrt(w)*unitHeight);
					Rectangle2D.Float r=null;
					
					switch(pinit.direction)
						{
						case Placement.DOWN_RIGHT:
							r=new Rectangle2D.Float(init.x,init.y,incx,incy);
							break;
						case Placement.DOWN_LEFT:
							r=new Rectangle2D.Float(init.x-incx,init.y,incx,incy);
							break;
						case Placement.UP_RIGHT:
							r=new Rectangle2D.Float(init.x,init.y-incy,incx,incy);
							break;
						case Placement.UP_LEFT:
							r=new Rectangle2D.Float(init.x-incx,init.y-incy,incx,incy);
							break;
						}
					center=new Point2D.Float((float)r.getCenterX(), (float)r.getCenterY());
				
					
					boolean intersected=false;
					for(Rectangle2D.Float p:polygons)	if(p.intersects(r))	{intersected=true;break;}
					
					if(boundingPolygon.contains(center) && !intersected)	//this one is good, put position	
						{
						//if fitting, because it's occupied, if not, because it's out of boundingPolygon, we add the square and inits
						//inits.remove(cont);
						
						//add up to three new inits --> TODO: the new three inits should depend on direction!
						//---
						/*
						Point2D.Float np=new Point2D.Float(r.x+r.width, r.y);
						if(np.x<x0+boundingPolygon.getBounds().width && !inits.contains(np))
							inits.addAll(getAllPlacements(np));	//some will be impossible, but it will be determined later on
							
						np=new Point2D.Float(r.x, r.y+r.height);
						if(np.y<y0+boundingPolygon.getBounds().height)
							inits.addAll(getAllPlacements(np));
						
						np=new Point2D.Float(r.x+r.width, r.y+r.height);
						if(np.y<y0+boundingPolygon.getBounds().height && np.x<x0+boundingPolygon.getBounds().width)
							inits.addAll(getAllPlacements(np));
							*/
						int xsign=1, ysign=1;
						
						switch(pinit.direction)
							{
							case Placement.DOWN_RIGHT:
								break;
							case Placement.DOWN_LEFT:
								xsign=-1;
								break;
							case Placement.UP_RIGHT:
								ysign=-1;
								break;
							case Placement.UP_LEFT:
								xsign=-1; ysign=-1; 
								break;
							}
						Point2D.Float np=new Point2D.Float(r.x+xsign*r.width, r.y);
						if(np.x<x0+boundingPolygon.getBounds().width && !inits.contains(np))
							inits.addAll(getAllPlacements(np));	//some will be impossible, but it will be determined later on
							
						np=new Point2D.Float(r.x, r.y+ysign*r.height);
						if(np.y<y0+boundingPolygon.getBounds().height)
							inits.addAll(getAllPlacements(np));
						
						np=new Point2D.Float(r.x+xsign*r.width, r.y+ysign*r.height);
						if(np.y<y0+boundingPolygon.getBounds().height && np.x<x0+boundingPolygon.getBounds().width)
							inits.addAll(getAllPlacements(np));
						//---
						currentCells[i].position[0]=center.x;
						currentCells[i].position[1]=center.y;
						if(mode==1)					//System.out.println("Placing "+i+" with area "+r.width+", "+r.height+ " at "+center.x+", "+center.y);
							System.out.println("Placing "+i+" with area "+r.width+", "+r.height+" and weight "+w+" at "+init.x+", "+init.y);
						
						sumArea+=r.width*r.height;
						//Check that the added polygon did not swallow some inits
						for(int k=0;k<inits.size();k++)
							{
							if(r.contains(inits.get(k).point))	
								{
								if(mode==1) System.out.println("point "+k+" swallowed, removing");
								inits.remove(k);
								}
							}
						
						polygons.add(new Rectangle2D.Float(init.x, init.y, r.width, r.height));
						
						break;
						}
					//cont++;
					}//for each init
				
				if(currentCells[i].position[0]==-1000000)
					{
					if(mode==1) System.out.println("Rescaling, not fitted at "+i+" with "+unitWidth+", "+unitHeight);
					unitWidth*=0.9;
					unitHeight*=0.9;
					
					fitting=false;
					break;
					}
				}//for each cell
			}
		return;
		}
	
	private class Placement
		{
		public Point2D.Float point;
		public int direction;
		public static final int UP_LEFT=0;
		public static final int UP_RIGHT=1;
		public static final int DOWN_LEFT=2;
		public static final int DOWN_RIGHT=3;
		public Placement(Point2D.Float p, int d)
			{
			point=p;
			direction=d;
			}
		}
	/**
	 * All possible placements of a polygon with a point p as initial point
	 * @param p
	 * @return
	 */
	private ArrayList<Placement> getAllPlacements(Point2D.Float p)
		{
		ArrayList<Placement> pl=new ArrayList<Placement>();
		for(int i=0;i<4;i++)
			pl.add(new Placement(p, i));
		return pl;
		}
	}