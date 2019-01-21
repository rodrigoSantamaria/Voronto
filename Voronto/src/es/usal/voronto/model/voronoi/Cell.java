package es.usal.voronto.model.voronoi;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import java.util.Arrays;
import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.ontology.OntologyTerm;
import es.usal.voronto.model.voronoi.geometry.MPolygon;
import es.usal.voronto.view.VoronoiVisualization;


public class Cell implements Cloneable {
	public float numLeaves;	//number of leaf nodes contained by this cell, indirectly determines weight and cell area
	public float weight;//weight assigned to the cell
	public float[] position;//location for voronoi diagram construction
	public float[] centroid;//centroid of its region polygon
	public OntologyTerm term;//term.geneIds must assure that contains every directly annotated gene plus every annotated gene in its subterms (without duplicates)
	public MPolygon region;
	public double desiredArea; //desired area based on initial weights
	public double areaRatio; //ro_i on Bernhardt et al. paper it's a measure of the difference between the desired and the actual area.
	public double convergenceFactor=1;
	public Cell[] subcells;//subcells in the hierarchy
	public int level=0;//level of the cell in the hierarchy
	public int id=-1;//TODO: by now, identifying by name
	public Color searchedColor=null;
	
	public ArrayList<Float> expressionLevel=new ArrayList<Float>(); //based on the expression of every gene in the term. If no gene is present on our expression data, it would be NaN and be drawn in grey
	public ArrayList<Float> expressionDeviation=new ArrayList<Float>();
	
	//public Color color=new Color(240,240,240);
	public ArrayList<Color> color=new ArrayList<Color>();
	
	//For label visualization
	public int labelX=-1;
	public int labelY=-1;
	public int labelSize=0;
	public String label;
	
	public Cell(float w, float[] p, String n)
		{
		numLeaves=w;
		position=p;
		term=new OntologyTerm(n, "");
		weight=numLeaves;
		}
	
	public Cell(float w, String n)
		{
		numLeaves=w;
		position=new float[2];
		centroid=new float[2];
		term=new OntologyTerm(n, "");
		weight=numLeaves;
		}	
	
	public Cell(float w, String n, int l)
	{
	numLeaves=w;
	position=new float[2];
	centroid=new float[2];
	term=new OntologyTerm(n, "");
	level=l;
	weight=numLeaves;
	}	
	
	public Cell(float w, OntologyTerm t, int l)
	{
	numLeaves=w;
	position=new float[2];
	centroid=new float[2];
	term=t;
	level=l;
	label=String.copyValueOf(term.name.toCharArray());
	weight=numLeaves;
	}	
	
	public Cell(OntologyTerm t, int l, ArrayList<Cell> sc)
	{
	term=t;
	level=l;
	position=new float[2];
	centroid=new float[2];
	label=String.copyValueOf(term.name.toCharArray());
	subcells=sc.toArray(new Cell[0]);
	int w=0;
	for(Cell cs:subcells)	w+=cs.numLeaves;
	weight=numLeaves=w;//Todo: change by total number of leaves!
	}

	
	/**
	 * Copy constructor
	 * @param c
	 */
	public Cell(Cell c)
		{
		numLeaves=c.numLeaves;
		position=new float[2];
		for(int i=0;i<position.length;i++)		position[i]=position[i];
		centroid=new float[2];
		for(int i=0;i<centroid.length;i++)		centroid[i]=centroid[i];
		term=new OntologyTerm("", "");
		term.name=c.term.name;
		term.id=c.term.id;
		term.geneIds=c.term.geneIds;
		term.koIds=c.term.koIds;
		
		level=c.level;
		weight=numLeaves;
		label=String.copyValueOf(c.label.toCharArray());
		
		desiredArea=c.desiredArea;
		areaRatio=c.areaRatio;
		convergenceFactor=c.convergenceFactor;
		id=c.id;
		
		if(c.region!=null)
				{
				region=new MPolygon(c.region.count());
				for(int i=0;i<c.region.count();i++)
					region.add(c.region.getCoords()[i][0], c.region.getCoords()[i][1]);
				}
		
		if(c.subcells!=null)
			{
			subcells=new Cell[c.subcells.length];
			for(int i=0;i<c.subcells.length;i++)
				subcells[i]=new Cell(c.subcells[i]);
			}
		}	
	
	/**
	 * Returns the number of genes in a expression matrix (md) that are in the KEGG OntologyTerm related to this cell.
	 * @param md
	 * @return
	 */
	public int getNumMatchedKOTerms(ExpressionData md)	//too slow compared with GO now
		{
		int cont=0;
		ArrayList<String> spcids=new ArrayList<String>();
		for(String id:term.geneIds)		
			if(id.startsWith(md.organismKegg+":"))	spcids.add(id.toLowerCase());
		for(String id:spcids)
			{
			if((md.koidList.indexOf(id))>=0)
				cont++;
			}
		
		return cont;
		}
	
	/**
	 * Returns the number of genes on this.term.geneIds that are in md.geneNamesHash
	 * @param md
	 * @return
	 */
	public int getNumMatchedTerms(ExpressionData md)
		{
		int cont=0;
		for(String id:this.term.geneIds)
			{
			if(Arrays.binarySearch(md.sortedGeneNames, id.toLowerCase())>=0)
				cont++;
			}
		
		return cont;
		}
	
	/**
	 * Returns the number of genes in a expression matrix (md) that are in the GO OntologyTerm related to this cell.
	 * @param md
	 * @return
	 */
	public int getNumMatchedGOTerms(ExpressionData md)
		{
		int cont=0;
		for(String id:term.geneIds)
			{
			if(Arrays.binarySearch(md.sortedGeneNames, id.toLowerCase())>=0)
				cont++;
			}
		
		return cont;
		}
	
	/**
	 * Gets the expression for a given column of a expression data matrix (md)
	 * @param column
	 */
	public float getExpression(int column)
		{
		return expressionLevel.get(column);
		}
	
	/**
	 * Links the ontology term to the expression of its corresponding genes
	 * @param md - expression data
	 * @param ontology - type of ontology (KEGG -0- or GO -1-)
	 */
	public void computeExpression(ExpressionData md, int ontology)
		{
		color=new ArrayList<Color>();
		if(md==null)
			{
			color.add(new Color(240,240,240));
			return;
			}
		for(int i=0;i<md.getNumConditions();i++)	
			color.add(new Color(240,240,240));

		//long t=System.currentTimeMillis();
		if(ontology==VoronoiVisualization.KEGG)
			computeKOExpression(md);
		else
			computeGOExpression(md);
		//System.out.println("Compute expression for "+term.name+" takes "+(System.currentTimeMillis()-t)/1000.0+" s");
		}
	
	public void computeExpressionFromChildren(ExpressionData md)
		{
		expressionLevel=new ArrayList<Float>();
		color=new ArrayList<Color>();
		for(int i=0;i<md.getNumConditions();i++)	color.add(new Color(240,240,240));
		
		int cont=0;
		for(int i=0;i<md.getNumConditions();i++)	//this might be in Cell.java
			{
			float temp=0;
			//float tempSD=0;
			cont=0;
			for(Cell cc:subcells)
				{
				temp+=cc.expressionLevel.get(i)*cc.numLeaves;
				cont+=cc.numLeaves;
				//tempSD+=cc.expressionDeviation.get(i);
				}
			if(cont==0)		
				{
				expressionLevel.add(Float.NaN);
				expressionDeviation.add(Float.NaN);
				}
			else			
				{
				expressionLevel.add(temp/cont);
				//cell.expressionDeviation.add(tempSD/cell.subcells.length);
				}
			//if(cell.expressionDeviation.size()==0)
			//	System.out.println(cell.term.name+" has deviation 0");
			}
		}
	
	public void computeExpressionProfile(ExpressionData md)
		{
		//0) Prepare cell expression data structure for gene in the cell and the expression data
		if(term.geneExs==null || term.geneExs.size()==0)
			{
			for(String id:term.geneIds)
				{
				if(Arrays.binarySearch(md.sortedGeneNames, id.toLowerCase())>=0)
					term.geneExs.put(id, new ArrayList<Float>());
				}
			}
	
		if(term.geneExs.size()>0)
			{
			Iterator<String> it=term.geneExs.keySet().iterator();
			while(it.hasNext())			
				{
				String i=it.next();
				double[] temp=md.getExpressionProfile(md.getGeneId(i));
				
				ArrayList<Float> profile=term.geneExs.get(i);
				for(int k=0;k<md.getNumConditions();k++)	
					{
					float level=(float)temp[k];
					profile.add(level);//NOTE: memory tests
					}
				}
			}
		
		return;
		}
	
	public void computeGOExpression(ExpressionData md)
		{
		expressionLevel=new ArrayList<Float>();
		
		
		//0) Prepare cell expression data structure for gene in the cell and the expression data
		for(String id:term.geneIds)
			{
			if(Arrays.binarySearch(md.sortedGeneNames, id.toLowerCase())>=0)
				term.geneExs.put(id.toLowerCase(), new ArrayList<Float>());
			}
	
		float[] e=new float[md.getNumConditions()];
		for(int j=0;j<e.length;j++)	e[j]=0;
		if(term.geneExs.size()>0)
			{
			//Computing average
			/*
			Iterator<String> it=term.geneExs.keySet().iterator();
			while(it.hasNext())			//NOTE: This loop can be slow (>1s)
				{
				String i=it.next();
				int id=md.getGeneId(i);
				double[] temp=md.getExpressionProfile(id);
				
				for(int k=0;k<md.getNumConditions();k++)	
					{
					float level=(float)temp[k];
					e[k]+=level;
					}
				}
			for(int k=0;k<md.getNumConditions();k++)	e[k]/=term.geneExs.size();
			*/
			//Computing median
			ArrayList<ArrayList<Float>> vals=new ArrayList<ArrayList<Float>>();
			for(int i=0;i<md.getNumConditions();i++)	vals.add(new ArrayList<Float>());
			
			Iterator<String> it=term.geneExs.keySet().iterator();
			while(it.hasNext())			
				{
				String i=it.next();
				int id=md.getGeneId(i);
				double[] temp=md.getExpressionProfile(id);
				
				for(int k=0;k<md.getNumConditions();k++)	
					vals.get(k).add((float)temp[k]);
				}
			
			for(int k=0;k<md.getNumConditions();k++)	
				{
				Float[] temp=vals.get(k).toArray(new Float[0]);
				Arrays.sort(temp);
				int middle = temp.length/2;
				if (temp.length%2 == 1) 	e[k]=temp[middle];
				else					e[k]= (float) ((temp[middle-1] + temp[middle]) / 2.0);
				}
			}
		else
			for(int j=0;j<e.length;j++)	e[j]=Float.NaN;
		
		for(int j=0;j<e.length;j++)	expressionLevel.add(e[j]);
		return;
		}
	/**
	 * Computes the expression for a given column of a expression data matrix (md)
	 * @param md
	 * @param column
	 */
	public void computeKOExpression(ExpressionData md)
	{
		
	//-------------------------- 1 FILL IN GENEEXS ------------------
	expressionLevel=new ArrayList<Float>();
	expressionDeviation=new ArrayList<Float>();
	float e=0;
	
	//long time=System.currentTimeMillis();
	//0) Do search only in the species
	ArrayList<String> spcids=new ArrayList<String>();
	for(String id:term.geneIds)		
		{
		if(id.startsWith(md.organismKegg) && !spcids.contains(id))	
			spcids.add(id.replace(md.organismKegg+":", ""));
		}
	//System.out.println("0)\t"+(System.currentTimeMillis()-time)/1000.0);
	//time=System.currentTimeMillis();
	//0b) And only take the elements that are in any KO term
	ArrayList<Integer> ids=new ArrayList<Integer>();
	ArrayList<String> idNames=new ArrayList<String>();
	for(String id:spcids)
		{
		int i=0;
		id=id.toLowerCase();
		if((i=Arrays.binarySearch(md.sortedGeneNames, id))>=0)
			{
			i=md.order.get(i);
			term.geneExs.put(id, new ArrayList<Float>());
			idNames.add(id);
			ids.add(i);
			}
		}
	
	//System.out.println("0b)\t"+(System.currentTimeMillis()-time)/1000.0);
	//0c) Compute the expression
	try{
		//expressionLevel=mean(ids, idNames, md);
		expressionLevel=median(md);
		}catch(Exception ex){ex.printStackTrace();}
	return;
	}
	
	public ArrayList<Float> mean(ArrayList<Integer> ids, ArrayList<String> idNames, ExpressionData md)
		{
		ArrayList<Float> el=new ArrayList<Float>();
 		for(int column=0;column<md.getNumConditions();column++)
			{
			float e=0;
			
			if(ids.size()>0)
				{
				for(int i=0;i<ids.size();i++)
					{
					float temp=(float)md.matrix[ids.get(i)][column];
					term.geneExs.get(idNames.get(i)).add(temp);
					e+=temp;
					}
				e/=ids.size();
				}
			else
				e=Float.NaN;
			
			el.add(e);
			}
		return el;
		}
	
	//public ArrayList<Float> median(ArrayList<Integer> ids, ArrayList<String> idNames, ExpressionData md)
	public ArrayList<Float> median(ExpressionData md)
		{
		ArrayList<Float> el=new ArrayList<Float>();
		float[] e=new float[md.getNumConditions()];
		ArrayList<ArrayList<Float>> vals=new ArrayList<ArrayList<Float>>();
		for(int i=0;i<md.getNumConditions();i++)	vals.add(new ArrayList<Float>());
		if(term.geneExs.size()>0)
			{
			Iterator<String> it=term.geneExs.keySet().iterator();
			while(it.hasNext())			
				{
				String i=it.next();
				int id=md.getGeneId(i);
				double[] temp=md.getExpressionProfile(id);
				
				for(int k=0;k<md.getNumConditions();k++)	
					vals.get(k).add((float)temp[k]);
				}
			
			for(int k=0;k<md.getNumConditions();k++)	
				{
				Float[] temp=vals.get(k).toArray(new Float[0]);
				Arrays.sort(temp);
				int middle = temp.length/2;
				if (temp.length%2 == 1) 	e[k]=temp[middle];
				else					e[k]= (float) ((temp[middle-1] + temp[middle]) / 2.0);
				}
			}
		else
			for(int j=0;j<e.length;j++)	e[j]=Float.NaN;
		
		for(int j=0;j<e.length;j++)	el.add(e[j]);
		
		return el;
		}
	
	
	/**
	 * Similar to computeKOExpression, but for KO terms instead of genes.
	 * This is necessary in the case that the KO pathway exists, but the specific species pathway does not.
	 * In that cases, we need to fill KO expression, not gene expression
	 * @param md
	 */
	public void computeKOIDExpression(ExpressionData md)
		{
		System.out.println("computeKOIDExpression");
		//-------------------------- 2 FILL IN KOEXS ------------------
		//0) Do search only in the KO terms that contain genes in the species
		ArrayList<String> spcids=new ArrayList<String>();
		for(String id:term.koIds.keySet().toArray(new String[0]))	//for each KO term
			{
			ArrayList<String> gids=term.koIds.get(id);				
			ArrayList<String> kospids=new ArrayList<String>();
			
			for(String gid:gids)
				if(gid.startsWith(md.organismKegg))	
					{
					if(!spcids.contains(id))	spcids.add(id);
					kospids.add(gid);
					}
			term.koIds.get(id).clear();
			term.koIds.get(id).addAll(kospids);
			//System.out.println(".");
			}

		for(String ko:spcids)	term.koExs.put(ko, new ArrayList<Float>());	//for KO terms in the pathway that have genes related to the expression species, but the expression data does not cover them

		
		//In the case of KEGG orthology, we also map expression to KOIDs
		for(int column=0;column<md.getNumConditions();column++)
			{
			for(String koid:term.koExs.keySet())
				{
				ArrayList<String> gids=term.koIds.get(koid);
				float koexp=0;
				int cont=0;
				for(String gid:gids)
					{
					int i=0;
					if((i=md.koidList.indexOf(gid.toLowerCase()))>=0)
						{
						koexp+=md.matrix[i][column];
						cont++;
						}
					}
				if(cont>0)	term.koExs.get(koid).add(koexp/cont);
				if(gids.size()==0)	term.koExs.get(koid).add(Float.POSITIVE_INFINITY);	//KO term in the pathway has no genes related to this species
				}
			}
		}

	
	/**
	 * Returns NaN if there are genes in the KO term for the species, but there are none of them in the expression data
	 * If there are 1 or more in the expression data, returns a number that is the average of their expressions
	 * If the KO term is in the pathway but has no assigned genes in the species, it leaves it infinite.
	 * 
	 * @param md
	 * @param column
	 * @return
	 */
	public Map<String, Float> computeKOExpression(ExpressionData md, int column)
		{
		Map<String,Float> koExpression=new TreeMap<String, Float>();//only for kegg orthology, we keep also track of gene expression for each KO term (that might include several genes)
		for(String ko:term.koIds.keySet())	koExpression.put(ko, Float.NaN);	//for KO terms in the pathway that have genes related to the expression species, but the expression data does not cover them

		if(term.koExs!=null && term.koExs.size()==term.koIds.size())
			{
			for(String koid:koExpression.keySet())
				{
				koExpression.put(koid, term.koExs.get(koid).get(column));
				}
			}
		else
			{
			//0) Do search only in the species
			ArrayList<String> spcids=new ArrayList<String>();
			for(String id:term.geneIds)		if(id.startsWith(md.organismKegg))	spcids.add(id);
	
			//In the case of KEGG orthology, we also map expression to KOIDs
			for(String koid:koExpression.keySet())
				{
				ArrayList<String> gids=term.koIds.get(koid);
					{
					float koexp=0;
					int cont=0;
					for(String gid:gids)
						{
						int i=0;
						if((i=md.koidList.indexOf(gid))>=0)
							{
							koexp+=md.matrix[i][column];
							cont++;
							}
						}
					if(cont>0)	koExpression.put(koid,koexp/cont);	//some genes in the expression, return its average
					if(gids.size()==0)	koExpression.put(koid, Float.POSITIVE_INFINITY);	//KO term in the pathway has no genes related to this species
					}
				}
			}
		return koExpression;
		}
	
	/*
	public void computeKOExpression(MicroarrayData md)
		{
		for(String ko:term.koIds.keySet())	term.koExs.put(ko, new ArrayList<Float>());	//for KO terms in the pathway that have genes related to the expression species, but the expression data does not cover them

		//0) Do search only in the species
		ArrayList<String> spcids=new ArrayList<String>();
		for(String id:term.geneIds)		if(id.startsWith(md.organismKegg))	spcids.add(id);

		//In the case of KEGG orthology, we also map expression to KOIDs
		for(int column=0;column<md.getNumConditions();column++)
			{
			for(String koid:term.koExs.keySet())
				{
				ArrayList<String> gids=term.koIds.get(koid);
				float koexp=0;
				int cont=0;
				for(String gid:gids)
					{
					int i=0;
					if((i=md.koidList.indexOf(gid))>=0)
						{
						koexp+=md.matrix[i][column];
						cont++;
						}
					}
				if(cont>0)	term.koExs.get(koid).add(koexp/cont);
				if(gids.size()==0)	term.koExs.get(koid).add(Float.POSITIVE_INFINITY);	//KO term in the pathway has no genes related to this species
				}
			}
		return;
		}*/
	
	public void translate(int x, int y)
		{
		if(region!=null)	region.translate(x, y);
		centroid[0]+=x;
		centroid[1]+=y;
		position[0]+=x;
		position[1]+=y;
		}
}
