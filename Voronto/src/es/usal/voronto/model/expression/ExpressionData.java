package es.usal.voronto.model.expression;



import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.StringTokenizer;

import javax.swing.JOptionPane;

import es.usal.voronto.model.ontology.KOparser;
import es.usal.voronto.model.voronoi.Cell;


/**
 * Class with data of Microarray expression levels, using Prefuse Tables
 * <p>
 *
 * MicroarrayData will contain:
 *<p>
 * For each gene/condition:	
 *<p>
 * 	"name" - name of the gene/condition (String)-, 
 *<p>
 * 	"id"   - unique identifier for the gene/condition (int)- and 
 *<p>
 * 	"rowRank"/"colRank" - order in which it is drawn (int). 
 *<p>
 * For each expression level:
 *<p>
 * 	"gene" - name of the gene for this expression level (String)-, 
 *<p>
 * 	"condition" - name of the condition for this expression level (String)-, 
 *<p>
 * 	"rowId" - id of the gene for this expression level (int)-, 
 *<p>
 * 	"colId" - id of the condition for this expression level (int)-, 
 *<p>
 * 	"rowRank" - order in which the gene is drawn (int)-, 
 *<p>
 *  "colRank" - order in which the condition is drawn (int)- and 
 *<p>
 *  "level" - expression level (double)
 *<p>
 *  
 *  In addition, an ancillary sparse matrix is built to avoid performance downgrades with very
 *  large matrices. It contains only numSparseGenes, typically a number around 200. The genes are
 *  randomly selected from the original matrix, and can also be selected by rebuilding the matrix
 *  with buildSparse(). It uses Prefuse Tables as described above, except it add a actualRowId that
 *  refers to the id in the original matrix, while rowId is used for the id in the sparse matrix.
 * @author Rodrigo Santamaria
 *
 */
public class ExpressionData 
	{
	/**
	 * matrix with expression levels replicated, to quicken some arithmetic operations
	 * TODO: remove, may lead to memory problems in large matrices
	 */
	public double matrix[][];
	int maxGenes=199;//A partir de 200
	
	//Stores the different R names of matrices loaded. Typically raw, merged, normalized, discretized, etc.
	//Typically, analysis and visualization will be done on the last matrix computed, but this is not always the case
	//for example, limma and gsea analysis cannot be done on merged matrices because they need several replicates of each efv
	public String filePath;//Absolute path with the file that corresponds to this data structure
	public String path;//Path for the folder containing the file that corresponds to this data structure
	public String name;//Name of the file that corresponds to this data structure
	public String rMatrixName;//Name in R for this data structure
	

	/**
	 * Names for the conditions.
	 */
	public String[] conditionNames;
	/**
	 * Gene names.
	 */
	public String[] geneNames;
	
	//Classes for quick searching in gene lists
	public HashMap<String, ArrayList<Integer>> geneNamesHash;
	public String[] sortedGeneNames;
	public HashMap<Integer, Integer> order;//maps ordered position to real position

	//Sorted gene names synonyms loaded from es/usal/voronto/data/synonym for the corresponding organismKegg, if available
	public String[] entrezgene;
	public String[] ensembl_gene_id;
	public String[] external_gene_id;
	//geneName as keys, synonym as values
	public HashMap<String,String> entrezgeneHash=null;
	public HashMap<String,String> ensembl_gene_idHash=null;
	public HashMap<String,String> external_gene_idHash=null;
	
	public static final int ENTREZ=0;
	public static final int SYMBOL=1;
	public static final int ENSEMBL=2;
	

	
	public String[] geneKOIDs;
	public ArrayList<String> koidList;
	
	/**
	 * Labels for the conditions. They are built initially from conditionNames and geneNames, but the user can change them
	 * to any combination of values present for gene ids or for experimental factors. 
	 */
	public String[] columnLabels;
	/**
	 * Gene names.
	 */
	public String[] rowLabels;
	
	public ArrayList<String> experimentFactors;
	/**
	 * Every key is a experimental factor, every array contains all the values for that factor (e.g. key="Time", values="0 min, 0 min, 20 min, 40 min") 
	 */
	public HashMap<String,String[]> experimentFactorValues;
	
	//For any other annotation present in the data in the first place (by now, not used)
	public ArrayList<String> geneFactors;
	public HashMap<String,String[]> geneFactorValues;
	
	
	/**
	 * Session class into which this microarray data is loaded
	 */
	//Session session=null;
	
	/**
	 * Type of gene names used. If GENENAME, annotations are searched with Entrez Gene and QuickGO
	 */
	public String rname="GENENAME";
	public String rdescription="DESCRIPTION";
	public String rgo="GO";
	/**
	 * Usually the R resource to search for gene annotations is an environment like for example illuminaHumanv1.db,
	 * but in some cases it is a database that cannot be accessed with get, mget; but with queries, and there is no
	 * environment counterpart. On these cases, in example Illumina Human HT12 for which there is no environment db
	 * but there is lumiHumanIDMapping package, isRDatabase is set to true and the search is done by dbGetQuery in R.
	 */
	public boolean isRDatabase=false;
	/**
	 * As isRDatabase, but with bioMaRt package, where getGene will be used for annotations
	 */
	public boolean isBioMaRt=false;
	
	/**
	 * As isRDatabase, but with GO.db package, where annotations of GO terms are in GOTERM
	 */
	public boolean isGO=false;
	
	
	int numGenes;
	int numSparseGenes;//for sparse matrices
	int numConditions;
	private int[] decimals;
	
	/**
	 * Average expression value of the whole expression matrix
	 */
	public double average=0;
	
	/**
	 * Maximum expression value of the whole expression matrix
	 */
	public double max;
	/**
	 * Minimum expression value of the whole expression matrix
	 */
	public double min;
	
	public String chip;//kind of microarray chip (any official for Affymetrix permitted, by now), or kind of name taken by genes (geneID and ORF permitted)
	public String organism;//Name of organism as registered in NCBI
	public boolean searchByR=false; //if there is information in the file about an available R package, gene annotations are taken from there, otherwise they're searched in NCBI
	boolean annotationsRetrieved=false;
	
	public int[] columnOrder;
	public double[] averageCols;
	public double[] minCols;
	public double[] maxCols;
	public double[] sdCols;
	public double meanSd;
	public double[] q75;//quantile 75 for each column
	public double[] q25;//quantile 25 for each column
	public double[] iqr;//interquantile (iq) range, computed as quantile75-quantile25
	public double whiskerRange=1.5;//this value multiplies iqr to get boxplot whiskers
	public double[] medianCols;
	public double median;
	public double[] quantiles=new double[100];//min. expression value to be in quantile i%
	public ArrayList<Double[]> quantileCols;
	
	
	public HashMap<Integer, int[]> outliers;//for each column, the gene above/below iqr*whiskerRange
	
	public String sortingFactor="Column ID";
	public String organismKegg;
	
	
	/**
	 * Constructor from a file
	 * @param path Path to the file with microarray information
	 * @param invert true if genes are columns (genes as rows are considered as the usual option)
	 * @param rowHeader Number of initial rows with row information (usually one)
	 * @param colHeader Number of initial columns with column information (usually one)
	 * @param nd	Number of decimals to be shown if numerically showing expression levels
	 */
	public ExpressionData(String path, String name, boolean invert, int rowHeader, int colHeader, int nd) throws Exception
		{
		//this.path=path.replace("\\", "/");
		this.path=path;
		experimentFactors=new ArrayList<String>();
		experimentFactorValues=new HashMap<String,String[]>();
		this.filePath=this.path;
		//this.name=path.substring(path.lastIndexOf("/")+1, path.lastIndexOf("."));
		this.name=name;
		this.rMatrixName=name.replace("-", "").replace(" ", ".");
		loadFromFile(rowHeader, colHeader);
		}
	
	public void simpleStats()
		{
		max=-1000000000;
		average=0;
		min=1000000000;
		maxCols=new double[numConditions];
		minCols=new double[numConditions];
		averageCols=new double[numConditions];
		sdCols=new double[numConditions];
		
		for(int i=0;i<this.numConditions;i++)
			{
			maxCols[i]=-1000000000;
			minCols[i]=1000000000;
			averageCols[i]=0;
			}
		for(int j=0;j<this.numConditions;j++)
			{
			for(int i=0;i<this.numGenes;i++)
				{
				if(maxCols[j]<matrix[i][j])	maxCols[j]=matrix[i][j];
				if(minCols[j]>matrix[i][j])	minCols[j]=matrix[i][j];
				if(max<matrix[i][j])	max=matrix[i][j];
				if(min>matrix[i][j])	min=matrix[i][j];
				averageCols[j]+=matrix[i][j];
				average+=matrix[i][j];
				}
			averageCols[j]/=numGenes;
			}
		average/=numGenes*numConditions;
	
		}

	/**
	 * Given an experimental factor (e.g. "Time") it returns its values for the conditions (e.g. "0 min", "0 min", "40 min", "40 min")
	 * @param experimentFactor
	 * @return
	 */
	public String[] getExperimentFactorValues(String experimentFactor)
		{
		return this.experimentFactorValues.get(experimentFactor);
		}

	
	/**
	 * Computes the standard deviation of each column. Mean should be previously computed (right now it's done in convert())
	 */
	private void computeSd()
		{
		long t=System.currentTimeMillis();
		meanSd=0;
		for(int i = 0; i < numConditions; i++)
			{
			sdCols[i]=0;
			for(int j = 0; j < numGenes; j++)
				sdCols[i]+=(averageCols[i]-matrix[j][i])*(averageCols[i]-matrix[j][i]);
			sdCols[i]/=numGenes;
			sdCols[i]=Math.sqrt(sdCols[i]);
			meanSd+=sdCols[i];
			}
		meanSd/=numConditions;
		System.out.println("Mean standard deviation: "+meanSd);
		System.out.println("Time to compute sd "+(System.currentTimeMillis()-t)/1000.0);
		}


	/**
	 * Builds a small form for the sample with id i
	 * @param i
	 * @return
	 */
	public String getDetailedSampleForm(int i)
		{
		String form="";
		form=form.concat("SampleId: "+conditionNames[i]+"\n");
		for(String s:experimentFactors)
			{
			form=form.concat(s+": "+experimentFactorValues.get(s)[i]+"\n");
			}
		form=form.concat("\n");
		return form;
		//TODO: Properly finish the form for samples
		}

	
	
	/**
	 * Returns the expression level of gene i under condition j
	 * @param i	position in the matrix of the gene
	 * @param j	position in the matrix of the condition
	 * @return	expression level in row i, column j
	 */
	public double getExpressionAt(int i, int j)
		{
		if(i>=0 && i<this.numGenes && j>=0 && j<this.numConditions)
			return matrix[i][j];
		else	return this.average;//TODO: this should be NA
		}

	public double[] getExpressionProfile(int i)
		{
		double[] profile=new double[numConditions];
		for(int j=0;j<numConditions;j++)
			profile[j]=matrix[i][j];
		return profile;
		}
	
	public float[] getAverageExpressionProfile(ArrayList<Integer> l)
		{
		float[] profile=new float[numConditions];
		for(int j=0;j<numConditions;j++)
			{
			profile[j]=0;
			for(int i:l)
				profile[j]+=matrix[i][j];
			profile[j]/=l.size();
			}
		return profile;
		}
	/**
	 * Returns the number of genes in the Microarray
	 * @return the number of genes in the Microarray
	 */
	public int getNumGenes() {
		return numGenes;
	}

	/**
	 * Returns the number of conditions in the Microarray
	 * @return the number of conditions in the Microarray
	 */
	public int getNumConditions() {
		return numConditions;
	}


	/**
	 * Returns the condition names
	 * @return	an array of strings with condition names, ordered by id
	 */
	public String[] getConditionNames() {
		return conditionNames;
	}
	
	/**
	 *	Returns the name of the condition at the specified position 
	 * @param pos	position of the condition in the matrix
	 * @return	name of the condition
	 */
	public String getConditionName(int pos)
		{
		return conditionNames[pos];
		}
	

	/**
	 *	Returns the name of the condition at the specified position 
	 * @param pos	position of the condition in the matrix
	 * @return	name of the condition
	 */
	public String getColumnLabel(int pos)
		{
		return columnLabels[pos];
		}

	/**
	 *	Returns the name of the condition at the specified position 
	 * @param pos	position of the condition in the matrix
	 * @return	name of the condition
	 */
	public String getRowLabel(int pos)
		{
		if(pos>=0 && pos<rowLabels.length)	return rowLabels[pos];
		else								return "NA";
		}

	/**
	 * Returns the gene names
	 * @return	an array of strings with gene names, ordered by id
	 */
	public String[] getGeneNames() {
		return geneNames;
	}
	
	public LinkedList<Integer> getGeneIds() {
		LinkedList<Integer> ret=new LinkedList<Integer>();
		for(int i=0;i<geneNames.length;i++)
			{
			ret.add(i);
			}
		return ret;
	}
	

	/**
	 * Returns the condition id of a condition name
	 * @param conditionName	condition name to know its id
	 * @return	id for the condition name, or -1 if it is not found
	 */
	public int getConditionId(String conditionName)
		{
		int id=-1;
		for(int i=0;i<conditionNames.length;i++)
			if(conditionNames[i].equals(conditionName))	
				return i;
		return id;
		}
	
	/**
	 * Returns the gene id of a gene name
	 * @deprecated: a gene name must correspond to more than one id
	 * @param geneName	gene name to know its id
	 * @return	id for the gene name, or -1 if it is not found
	 */
	public int getGeneId(String geneName)
		{
		int ret=-1;
		ArrayList<Integer> l=geneNamesHash.get(geneName.toLowerCase());
		if(l!=null)	ret=l.get(0);
		return ret;
		}
	
	public ArrayList<Integer> getGeneIds(String geneName)
		{
		return geneNamesHash.get(geneName.toLowerCase());
		}
	
	public ArrayList<Integer> getGeneIds(String[] geneNames)
		{
		ArrayList<Integer> ret=new ArrayList<Integer>();
		for(String g:geneNames)
			{
			ArrayList<Integer> id=geneNamesHash.get(g.toLowerCase());
			if(id!=null)	ret.addAll(id);
			}
		return ret;
		}
	
	
	/**
	 * Returns the constance by rows, columns or both, for a subset and genes and/or conditions
	 * The subset of genes and the subset of conditions may be null, but not both of them. If one 
	 * of them is null, all gene/condition profile is considered.
	 * @param genes - subgroup of genes
	 * @param conditions - subgroup of conditions 
	 * @param type - type of constance, by rows (0), by columns (1) or by both rows and columns (2)
	 * @return the constance value for the subgroup
	 */
	public float getConstance(ArrayList<String> genes, ArrayList<String> conditions, int type)
		{
		float constance=0;
		double matrixBic[][]=new double[genes.size()][conditions.size()];
		//1) Recuperamos la matriz sobre la que calcular la constancia
		for(int i=0;i<genes.size();i++)
			{
				String gene=genes.get(i);
			int row=getGeneId(gene);
			for(int j=0;j<conditions.size();j++)
				{
				String cond=conditions.get(j);
				int col=getConditionId(cond);
				matrixBic[i][j]=matrix[row][col];
				}
			}
		//2) Calculamos la constancia
		float sum=0;
		float sd=0;
		int num=0;
		switch(type)
			{
			case 0:
				num=genes.size();
				double[] midpoint=new double[conditions.size()];
				for(int i=0;i<conditions.size();i++)
					for(int j=0;j<genes.size();j++)
						midpoint[i]+=matrixBic[j][i];
				for(int i=0;i<conditions.size();i++)	midpoint[i]/=genes.size();
				for(int i=0;i<genes.size();i++)
					{
					double[] p0=new double[conditions.size()];
					for(int j=0;j<conditions.size();j++)	p0[j]=matrixBic[i][j];
					sd+=euclideanDistance(p0, midpoint);
					}
				constance=(float)Math.sqrt(sd/num);
					
				break;
			case 1:
				num=conditions.size();
				break;
			case 2:
				num=genes.size()*conditions.size();
				for(int i=0;i<genes.size();i++)
					for(int j=0;j<conditions.size();j++)
						sum+=matrixBic[i][j];
				for(int i=0;i<genes.size();i++)
					for(int j=0;j<conditions.size();j++)
						sd+=Math.abs(matrixBic[i][j]-sum);
				constance=sd/num;
				break;
			}
		
		return constance;
		}
	
	/**
	 * Returns the euclidean distance between two points
	 * @param p1 first point
	 * @param p2 second point
	 * @return euclidean distance
	 * TODO: move to utils
	 */
	private double euclideanDistance(double[] p1, double[] p2)
	{
	double ret=0;
	for(int i=0;i<p1.length;i++)
		ret+=(p1[i]-p2[i])*(p1[i]-p2[i]);
	return Math.sqrt(ret);
	}
    
	/**
	 * Returns a list of gene and condition names from two lists with genes and conditions ids
	 * @param lg	list of gene ids 
	 * @param lc	list of condition ids 
	 * @return	An ArrayList with Strings with corresponding gene and condition names
	 */
	public ArrayList<String> getNames(LinkedList<Integer> lg, LinkedList<Integer> lc)
		{
		ArrayList<String> ret=new ArrayList<String>();
		for(int i=0;i<lg.size();i++)		ret.add(geneNames[lg.get(i)]);
		for(int i=0;i<lc.size();i++)		ret.add(conditionNames[lc.get(i)]);
		return ret;
		}
	/**
	 * Returns a list of gene names from a list with genes ids
	 * @param lg	list of gene ids
	 * @return	An ArrayList with Strings with corresponding gene and condition names
	 */
	public HashSet<String> getGeneNames(LinkedList<Integer> lg)
		{
		HashSet<String> ret=new HashSet<String>();
		for(int i=0;i<lg.size();i++)		ret.add(geneNames[lg.get(i)].trim().toLowerCase());
		return ret;
		}
	
	
	public String getGeneName(int lg)
	{
	return geneNames[lg];
	}
	
	/**
	 * Returns a list of gene names from a list with genes ids
	 * @param lg	list of gene ids
	 * @return	An ArrayList with Strings with corresponding gene and condition names
	 */
	public ArrayList<String> getRowLabels(LinkedList<Integer> lg)
		{
		ArrayList<String> ret=new ArrayList<String>();
		for(int i=0;i<lg.size();i++)		ret.add(rowLabels[lg.get(i)]);
		return ret;
		}
	
	
	/**
	 * Returns a list of condition names a lists and conditions ids
	 * @param lc	list of condition ids 
	 * @return	An ArrayList with Strings with corresponding gene and condition names
	 */
	public ArrayList<String> getConditionNames(LinkedList<Integer> lc)
		{
		ArrayList<String> ret=new ArrayList<String>();
		if(conditionNames!=null)	
			for(int i=0;i<lc.size();i++)		ret.add(conditionNames[lc.get(i)]);
		return ret;
		}
	
	/**
	 * Returns a list of condition names a lists and conditions ids
	 * @param lc	list of condition ids 
	 * @return	An ArrayList with Strings with corresponding gene and condition names
	 */
	public ArrayList<String> getColumnLabels(LinkedList<Integer> lc)
		{
		ArrayList<String> ret=new ArrayList<String>();
		if(columnLabels!=null)	
			for(int i=0;i<lc.size();i++)		ret.add(columnLabels[lc.get(i)]);
		return ret;
		}
	
	public HashMap<String, String> invertedHash(HashMap<String, String> hm)
		{
		HashMap<String,String> ih=new HashMap<String, String>();
		Iterator<String> it=hm.keySet().iterator();
		while(it.hasNext())
			{
			String k=it.next(); 
			ih.put(hm.get(k), k);
			}
		return ih;
		}


	/**
	 * Returns the number of genes in the sparse matrix
	 * @return
	 */
	public int getNumSparseGenes() {
		return numSparseGenes;
	}
	/**
	 * Sets the number of genes in the sparse matrix
	 * @param numSparseGenes
	 */
	void setNumSparseGenes(int numSparseGenes) {
		this.numSparseGenes = numSparseGenes;
	}
	
	
	
	public LinkedList<Integer> searchGenes(String what, boolean exact)
		{
		LinkedList<Integer> genes=new LinkedList<Integer>();
		for(String s:geneNames)
			{
			if(!exact)
				{
				if(s.contains(what))
					{
					genes.add(getGeneId(s));
					//System.out.println("Adding gene "+s+"\t"+getGeneId(s)+"\t"+this.getGeneNames()[getGeneId(s)]);
					}
				}
			else
				{
				if(s.trim().equals(what.trim()))
					{
					genes.add(getGeneId(s));
					//System.out.println("Adding gene "+s+"\t"+getGeneId(s)+"\t"+this.getGeneNames()[getGeneId(s)]);
					}
				}
					
			}
		return genes;
		}
	
	public LinkedList<Integer> differentialExpressionGenes(int[] g1, int[] g2, float threshold, int type)
		{
		LinkedList<Integer> retList=new LinkedList<Integer>();
		for(int g=0;g<numGenes;g++)
			{
			double e1=0;
			for(int i:g1)			e1+=matrix[g][i];
			double e2=0;
			for(int i:g2)			e2+=matrix[g][i];
			e1/=g1.length; e2/=g2.length;
			
			if(type==0 && e1>e2+threshold)	
				retList.add(g);
			else if(type==1 && e2>e1+threshold)
				retList.add(g);
			}
		
		return retList;
		}
	
	public LinkedList<Integer> searchConditions(String what)
		{
		LinkedList<Integer> conditions=new LinkedList<Integer>();
		for(String s:conditionNames)
			{
			if(s.contains(what))
				conditions.add(getGeneId(s));
			}
		
		for(int i=0;i<numConditions;i++)
			for(String cad:this.experimentFactors)
				{
				if(experimentFactorValues.get(cad)[i].contains(what))
					conditions.add(i);
				}
			
		return conditions;
		}
/**
	 * Returns a string with the requested value as a text string, for a determinate dimension
	 * @param value	value to format to text
	 * @param dim	dimension to which the value pertains (it determines the number of decimals for that dimension)
	 * @return	a text string with the value, formatted to have the number of decimals set for its dimension
	 */
	public String format(double value, int dim)
		{
		int nc=decimals[dim];
		String cad;
		if(nc>0)	cad=new Double(value).toString();
		else		cad=new Double(Math.abs(Math.rint(value))).toString();
		int pos=0;
		if((pos=cad.indexOf("."))>=0)
			{
			int pos2=0;
			if((pos2=cad.indexOf("E"))>=0)
				{
				if(nc>0)	cad=cad.substring(0,pos)+cad.substring(pos, pos+1+nc)+cad.substring(pos2);
				else		cad=cad.substring(0,pos)+cad.substring(pos2);
				}
			else
				{
				if(nc>0)
					{
					if(pos+1+nc<=cad.length())	
						{
						double number;
						//double number = (double)(int)((value+0.005)*100.0)/100.0;
						if(nc>0) number = (double)(int)((value+Math.signum(value)*((0.1*(nc))/2))*(10.0*(nc)))/(10.0*(nc));
						else	 number = Math.rint(value);
						cad=new Double(number).toString();
						}
					else						;//En este caso ya esta, no tenemos que quitar decimales
					}
				else		cad=cad.substring(0,pos);
				}
			}
		
		return cad;
		}


	
		/**
	 * Returns the corresponding annotation package for the organism of the loaded microarray
	 * @return
	 */
	public String getAnnotationPackage()
		{
		String cad=null;
		StringTokenizer st=new StringTokenizer(organism);
		int numTokens=st.countTokens();
		String abbv=""+st.nextToken().toUpperCase().charAt(0)+st.nextToken().toLowerCase().charAt(0);
		if(numTokens==2)
			{
			//it works for everything except E. coli and Malaria
			cad="org."+abbv+".eg.db";
			
			if(abbv.equals("Pf"))//this makes malaria
				cad="org."+abbv+".plasmo.db";
			if(cad.contains("Ec"))//this makes E. coli (default K12)
				cad="org."+abbv+"K12.eg.db";
			}
		else
			{
			//this makes E. coli distinction
			if(organism.contains("K12"))
				cad="org."+abbv+"K12.eg.db";
			else if(organism.contains("Sakai"))
				cad="org."+abbv+"Sakai.eg.db";
			}
		return cad; 
		}
	
	public int loadFromFile(int rowHeader, int colHeader) throws Exception
	{
	String name=path.replace("\\","/");
	filePath=name;
	String description="";
	System.out.println("Loading "+path);
	
	try{
		BufferedReader in =	new BufferedReader(new FileReader(path));
		StringTokenizer sto=new StringTokenizer(in.readLine(),"\t");
		int numLines=1;
		numConditions=sto.countTokens();
		System.out.println("numColumns is "+numConditions);
		while(in.readLine()!=null)	numLines++;	

		System.out.println("numLines is "+numLines);
		in =	new BufferedReader(new FileReader(path));
		rowHeader=0;
		boolean stop=false;
		for(int i=0;i<numLines;i++)	
			{
			String st=in.readLine();
			Scanner scanner = new Scanner(st);
			scanner.useDelimiter("\\t");
			String scan="";
			if(i>0)
				{
				while(scanner.hasNext())
					{
				    if (scanner.hasNextDouble()) 
			        	{stop=true; break;}
				    scan=scanner.next();
				    if(scan.trim().startsWith("EF."))	break;	//To allow numeric values
				    double scand=-333;
				    try{scand=Double.parseDouble(scan.trim());}catch(NumberFormatException nfe){}
				    if(scand!=-333){stop=true; break;}
				    }
		        if(stop)
		        	{
		        	//System.out.println("Out because of "+scanner.next());
		        	System.out.println("Out because of "+scan);
		        	break;
		        	}
		        }
			rowHeader++;
	        }

		System.out.println("Number of row headers: "+rowHeader);
		numConditions-=colHeader;
		numGenes=numLines-rowHeader;
		
		}catch(Exception e){e.printStackTrace();}
			
	conditionNames=new String[numConditions];
	
	geneNames=new String[numGenes];
	entrezgene=new String[numGenes];
	ensembl_gene_id=new String[numGenes];
	external_gene_id=new String[numGenes];
	
	geneNamesHash=new HashMap<String, ArrayList<Integer>>();
	geneKOIDs=new String[numGenes];
	koidList=new ArrayList<String>();
	columnOrder=new int[numConditions];
	
	for(int i=0;i<numConditions;i++)	columnOrder[i]=i;
	System.out.println("Matrix with "+numGenes+" rows and "+numConditions+" columns");

	BufferedReader in =	new BufferedReader(new FileReader(path));
	StringTokenizer st=null;
		
	//Read sample related info
	if(rowHeader>=0)
		{
		String cad=in.readLine();
		st=new StringTokenizer(cad, "\t"); //First row contains sample names and the different gene names maybe
		for(int i=0;i<colHeader;i++)		//Depending on the colHeader, it can also contain additional information. TODO right now colHeader is always 1, with organism/geneID info
			{
			description=st.nextToken();//Pasamos los que no tienen q ver, que nombran las columnas de explicaci�n de genes
			
			if(description.contains("/"))
				{
				chip=description.substring(description.indexOf("/")+1);
				organism=description.substring(0, description.indexOf("/"));
				organismKegg=ExpressionData.organismKEGG(organism);
				System.out.println("Organism KEGG: "+organismKegg);
			
				if(chip.length()<2)	
					{
					JOptionPane.showMessageDialog(null,
							"Chip name is wrong: "+chip, 
							"Wrong format", JOptionPane.ERROR_MESSAGE);
					return 1;
					}
				if(organism.length()<2)
					{
					JOptionPane.showMessageDialog(null,
							"Organism name is wrong: "+organism, 
							"Wrong format", JOptionPane.ERROR_MESSAGE);
					return 1;
					}
				}
			else
				{
				throw new Exception("Organism/IDtype expected in the first cell (found: '"+description+"')");
				}
			}
		
		for(int i=0;i<numConditions;i++)	
			conditionNames[i]=st.nextToken().trim();
		for(int i=0;i<rowHeader-1;i++)	//Read experiment factors
			{
			st=new StringTokenizer(in.readLine(), "\t"); 
			String ef=st.nextToken().trim();
			for(int j=0;j<colHeader-1;j++)	st.nextToken();//Avoid the blanks due to column Headers
			String[] efvs=new String[numConditions];
			int cont=0;
			while(st.hasMoreTokens())		efvs[cont++]=st.nextToken().trim();//.replace("EF.", "")
			experimentFactors.add(ef);
			experimentFactorValues.put(ef, efvs);
			}
		}
		else	//Use integers as identifiers
			for(int i=0;i<numConditions;i++)			conditionNames[i]=new Integer(i).toString().trim();
			
		//Read gene related info (and expression levels)
		matrix=new double[numGenes][numConditions];
		if(colHeader>0)
			{
			for(int i=0;i<numGenes;i++)//Read gene names
				{
				String cad=in.readLine();
				//System.out.println(cad);
				if(cad==null)
					{
					throw new Exception("Possibly empty line at the end of the matrix, please remove.");
					}
				if(cad.contains("\t\t"))
					{
					throw new Exception("Empty sample name(s) found, remove additional tabs between fields.");
					}
				st=new StringTokenizer(cad,"\t");//El delimitador en Syntren es un tab.
				geneNames[i]=st.nextToken().trim();
				if(geneNamesHash.get(geneNames[i].toLowerCase())==null)	
					{
					ArrayList<Integer> l=new ArrayList<Integer>();
					l.add(i);
					geneNamesHash.put(geneNames[i].toLowerCase(),l);
					}
				else	
					geneNamesHash.get(geneNames[i].toLowerCase()).add(i);
				geneKOIDs[i]=organismKegg+":"+geneNames[i];//TODO: this might not be true for all kinds of ids and species
				koidList.add(geneKOIDs[i].toLowerCase());//TODO: must be changed not to geneNames but to the corresponding id for KEGG/species
				
				try{
					for(int j=0;j<numConditions;j++)
						{
						//System.out.println(geneNames[i]);
						matrix[i][j]=new Double(st.nextToken().replace(",", ".")).doubleValue();
						}
				}catch(Exception e)
					{
					throw new Exception("Number format error on expression levels ("+e.getMessage()+")");
					}
				//System.out.println(geneNames[i]);
				}
			}
		else	//Set up numbers as geneNames
			for(int i=0;i<numGenes;i++)	geneNames[i]=new Integer(i).toString().trim();
	
	//average non-unique gene ids
	//averageNonUniqueGeneNames();	//TODO: support non unique gene row ids 
	rowLabels=geneNames.clone();
	columnLabels=conditionNames.clone();
	
	simpleStats();
	computeSd();
	computeQuantiles();
	sortGeneNames();
	
	generateSynonyms();
	
	return 0;
	}

	
	
private void generateSynonyms()
	{
	try{
		BufferedReader in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/synonym/"+organismKegg+".syn")));
		entrezgene=Arrays.copyOf(geneNames, geneNames.length);
		external_gene_id=Arrays.copyOf(geneNames, geneNames.length);
		ensembl_gene_id=Arrays.copyOf(geneNames, geneNames.length);
		entrezgeneHash=new HashMap<String,String>();
		ensembl_gene_idHash=new HashMap<String,String>();
		external_gene_idHash=new HashMap<String,String>();
		
		String cad=null;
		int cont=0;
		while((cad=in.readLine())!=null)
			{
			String[] sp=cad.split("\t");
			String ens=null;
			String ent=sp[0];
			String ext=sp[1];
			if(sp.length>2)	ens=sp[2];
			int pos=-1;
			
			if(this.chip.equals("entrezgene"))
				pos=Arrays.binarySearch(sortedGeneNames, ent.toLowerCase());
			else if(this.chip.equals("ensembl_gene_id"))
				pos=Arrays.binarySearch(sortedGeneNames, ens.toLowerCase());
			else if(this.chip.equals("external_gene_id"))
				pos=Arrays.binarySearch(sortedGeneNames, ext.toLowerCase());
			if(pos>0)
				{
				entrezgene[pos]=ent;
				ensembl_gene_id[pos]=ens;
				external_gene_id[pos]=ext;
				if(ent!=null && ent.length()>0)	entrezgeneHash.put(sortedGeneNames[pos],ent.toLowerCase());
				if(ens!=null && ens.length()>0)	ensembl_gene_idHash.put(sortedGeneNames[pos],ens.toLowerCase());
				if(ext!=null && ext.length()>0)	external_gene_idHash.put(sortedGeneNames[pos],ext.toLowerCase());
				cont++;
				}
			}
		}
	catch(Exception ioe)
		{
		System.out.println("No synonym list recorded for organism "+organismKegg);
		}
	}

public void sortGeneNames()
	{
	long time=System.currentTimeMillis();
	sortedGeneNames=Arrays.copyOf(geneNames, geneNames.length);
	for(int i=0;i<sortedGeneNames.length;i++)	
		sortedGeneNames[i]=sortedGeneNames[i].toLowerCase();
	
	order=new HashMap<Integer,Integer>();
	Arrays.sort(sortedGeneNames);
	for(int i=0;i<geneNames.length;i++)
		{
		String gene=geneNames[i].toLowerCase();
		int oi=Arrays.binarySearch(sortedGeneNames, gene);
		order.put(oi, i);
		}
	System.out.println("Time in sorting: "+(time-System.currentTimeMillis())/1000.0);
	return;
	}
public void computeQuantiles()
	{
	medianCols=new double[numConditions];
	quantileCols=new ArrayList<Double[]>();
	
	long t=System.currentTimeMillis();
	double [] numbers=new double[numConditions*numGenes];
	int cont=0;
	double step=numGenes*0.01;
	
	for(int j=0;j<numConditions;j++)
		{
		Double[] qq=new Double[100];
		double []numberCols=new double[numGenes];
		for(int i=0;i<numGenes;i++)
			{
			numbers[cont++]=matrix[i][j];
			numberCols[i]=matrix[i][j];
			}
		Arrays.sort(numberCols);
		medianCols[j]=numberCols[(int)(numGenes*0.5)];
		for(int i=0;i<100;i++)
			{
			qq[i]=numberCols[(int)(step*i)];
			//System.out.println("quantile "+(i+1)+" for column "+j+" up to "+qq[i]);
			}
		
		quantileCols.add(qq);
		}
	Arrays.sort(numbers);
	median=numbers[(int)(numGenes*numConditions*0.5)];
	
	step=numGenes*numConditions*0.01;
	for(int i=0;i<100;i++)
		{
		quantiles[i]=numbers[(int)(step*i)];
		//System.out.println("quantile "+(i+1)+" up to "+quantiles[i]);
		}
	System.out.println("Median computation takes "+(System.currentTimeMillis()-t)/1000.0+" median="+median);
	}

/*
 * Median of the expression of a number of genes under a given column
 */
public float median(HashSet<String> genes, int column) 
	{
	float[] m=new float[genes.size()];
	Iterator<String> it=genes.iterator();
	for(int i=0;i<genes.size();i++)	m[i]=(float)matrix[getGeneId(it.next().replaceAll("\\(.*\\)", ""))][column];
	Arrays.sort(m);
	int middle = m.length/2;
	if (m.length%2 == 1) 	return m[middle];
	else					return (float) ((m[middle-1] + m[middle]) / 2.0);
	}
//public float median(float[] values)
/*
 * Median of the expression of a number of genes
 */
public float[] median(HashSet<String> genes)
	{
	float [] med=new float[numConditions];
	for(int i=0;i<med.length;i++)
		med[i]=median(genes, i);
	return med;
	}

public void loadMicroarray(String path, boolean invert, int rowHeader, int colHeader, int nd) 
	{
	
	String description="";
	System.out.println(path);
	
	conditionNames=new String[numConditions];
	geneNames=new String[numGenes];
	decimals=new int[numConditions];
	for(int i=0;i<numConditions;i++)	decimals[i]=nd;
	System.out.println("Microarray matrix with "+numGenes+" genes and "+numConditions+" conditions");
	
	try{
	//Leemos la primera fila, que tiene los nombres de los genes
	if(invert)
		{
		if(rowHeader>=1)
			{
			//BufferedReader in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));
			BufferedReader in =	new BufferedReader(new FileReader(path));
			StringTokenizer st=new StringTokenizer(in.readLine(),"\t");//El delimitador en Syntren es un tab.
			for(int i=0;i<numGenes;i++)	
				geneNames[i]=st.nextToken();
			}
		else	for(int i=0;i<numGenes;i++)	geneNames[i]=new Integer(i).toString();
		
		
		if(colHeader>=1)
			{
			}
		else	for(int i=0;i<numConditions;i++)			conditionNames[i]=new Integer(i).toString();
		}
	else	//Si no est�n invertidas las primeras columnas tienen informaci�n sobre los genes
		{	//Y las primeras filas sobre las condiciones
		if(colHeader>=0)
			{
			BufferedReader in =	new BufferedReader(new FileReader(path));
			for(int i=0;i<rowHeader;i++)	in.readLine(); //Pasamos las filas con explicaciones sobre los experimentos
			for(int i=0;i<numGenes;i++)
				{
				String cad=in.readLine();
				if(cad.contains("\t\t"))
					{
					JOptionPane.showMessageDialog(null,
							"Empty sample name(s) found, remove additional tabs between fields.", 
							"Wrong format", JOptionPane.ERROR_MESSAGE);
					}
				StringTokenizer st=new StringTokenizer(cad,"\t");//El delimitador en Syntren es un tab.
				geneNames[i]=st.nextToken();
				}
			}
		else	for(int i=0;i<numGenes;i++)	geneNames[i]=new Integer(i).toString();
		
		
		if(rowHeader>=1)
			{
			BufferedReader in =	new BufferedReader(new FileReader(path));
			StringTokenizer st=new StringTokenizer(in.readLine(),"\t");//El delimitador en Syntren es un tab.
			for(int i=0;i<colHeader;i++)		
				{
				description=st.nextToken();//Pasamos los que no tienen q ver, que nombran las columnas de explicaci�n de genes
				}
			for(int i=0;i<numConditions;i++)	conditionNames[i]=st.nextToken();
			}
		else	for(int i=0;i<numConditions;i++)			conditionNames[i]=new Integer(i).toString();
		}
	}catch(Exception e){System.err.println("loadMicroarray: Error reading file "+path); e.printStackTrace(); System.exit(1);}
	
	if(description.contains("/"))
		{
		chip=description.substring(description.indexOf("/")+1);
		organism=description.substring(0, description.indexOf("/"));
		
		
		if(chip.length()<2)	
			JOptionPane.showMessageDialog(null,
					"Chip name is wrong: "+chip, 
					"Wrong format", JOptionPane.ERROR_MESSAGE);
		if(organism.length()<2)	
			JOptionPane.showMessageDialog(null,
					"Organism name is wrong: "+organism, 
					"Wrong format", JOptionPane.ERROR_MESSAGE);
		}
	else
		{
		JOptionPane.showMessageDialog(null,
				"Chip/Organism description is wrong, use \"/\" as separator: "+description, 
				"Wrong format", JOptionPane.ERROR_MESSAGE);
		}

	}
	
	
	
	/**
	 * Returns the current columnOrder c, where c[i]=j means that the current position for column i is now j.
	 * @return
	 */
	public int[] getColumnOrder()
		{
		int[] co=new int[columnOrder.length];
		for(int j=0;j<columnOrder.length;j++)		co[columnOrder[j]]=j;  
		return co;
		}
	
	/**
	 * TODO: looks like not every organism follows this rule!!
	 * @param organism
	 * @return
	 */
	public static String organismKEGG(String organism)
		{
		HashMap<String,String> orgmap=KOparser.getKeggOrganismMap();
		
		String cad=organism;
		if(organism.indexOf(" ")>0)
			{
			cad=orgmap.get(organism);
			//deprecated, it's now more complicated than just this
			/*cad=organism.charAt(0)+organism.substring(organism.indexOf(" ")+1,organism.indexOf(" ")+3);
			if(organism.startsWith("Macaca mu"))
				cad="mcc";	//special case in kegg, if not it will take mmu!
			if(organism.startsWith("Macaca fa"))
				cad="mcf";	//special case in kegg, if not it will take mmu!*/
			}
		return cad.toLowerCase();
		}
	
	/**
	 * Selecciona los genes con valor de expresi�n por encima de sdsAbove desviaciones est�ndar de la media para todas las condiciones
	 * etiquetadas con el valor highEFV para el factor experimental highEF y por debajo de sdsBelow desviaciones est�ndar de la media para 
	 * todas las condiciones etiquetadas con el valor lowEFV para el factor experimental lowEF
	 * highEF or lowEF (but not both) can be "none333" if we only want to filter by one condition.
	 * Example: we want to get every gene above 2 sd for "diseased" but below mean for "control" for experimental factor "cancer" 
	 * selectHiLo("diseased", "cancer", 2, "control", "cancer", 0);
	 * @param highEFV
	 * @param highEF
	 * @param sdsAbove
	 * @param lowEFV
	 * @param lowEF
	 * @param sdsBelow
	 * @return
	 */
	public LinkedList<Integer> selectHiLo(String highEFV, String highEF, int sdsAbove, String lowEFV, String lowEF, int sdsBelow)
		{
		if(!experimentFactors.contains(highEF) && !highEF.equals("none333") && !highEF.equals("rest") )
			{System.err.println("Experimental factor "+highEF+" does not exist"); return null;}
		if(!experimentFactors.contains(lowEF) && !lowEF.equals("none333") && !lowEF.equals("rest"))
			{System.err.println("Experimental factor "+highEF+" does not exist"); return null;}
		
		LinkedList<Integer> ret=new LinkedList<Integer>(); 
		String[] efvsH=experimentFactorValues.get(highEF);
		String[] efvsL=experimentFactorValues.get(lowEF);
		if(highEF.equals("rest"))	efvsH=efvsL;
		if(lowEF.equals("rest"))	efvsL=efvsH;
		
		for(int i=0;i<numGenes;i++)
			{
			boolean add=true;
			for(int j=0;j<numConditions;j++)
				{
				if(!highEFV.equals("none333") && ((highEFV.equals("rest") && !efvsH[j].equals(lowEFV)) || efvsH[j].equals(highEFV))) //TODO:the first part of the if could be done less times.
					if(matrix[i][j]<averageCols[j]+sdsAbove*sdCols[j])	{add=false; break;}
				if(!lowEFV.equals("none333") && ((lowEFV.equals("rest") && !efvsL[j].equals(highEFV)) || efvsL[j].equals(lowEFV)))
					if(matrix[i][j]>averageCols[j]+sdsBelow*sdCols[j])	{add=false; break;}
				}
			if(add)	ret.add(i);
			}
		return ret;
		}
	/**
	 * Returns the ids of the conditions that have the corresponding experimental factor value for a given experimental factor
	 * @param ef - experimental factor to be checked
	 * @param efv - experimental factor value that have the conditions to be returned
	 * @param notEqual - if true, returns all the conditions except the ones corresponding to the efv for ef 
	 * @return
	 */
	public Integer[] getConditions(String ef, String efv, boolean notEqual)
		{
		if(ef==null || efv==null || ef.length()==0 || efv.length()==0 )
			{System.err.println("No ef or efv specified"); return null;}
		LinkedList<Integer> ret=new LinkedList<Integer>();
		String[] efvs=experimentFactorValues.get(ef);
		for(int i=0;i<efvs.length;i++)
			{
			if(notEqual)	{if(!efvs[i].equals(efv))	ret.add(i);}
			else			{if(efvs[i].equals(efv))		ret.add(i);}
			}
		return ret.toArray(new Integer[ret.size()]);
		}
	
	/**
	 * Selecciona los genes  con valor 
	 * @param highEFV
	 * @param highEF
	 * @param lowEFV
	 * @param lowEF
	 * @return
	 */
	public LinkedList<Integer> selectHiLo(String highEFV, String highEF, String lowEFV, String lowEF)
		{
		return selectHiLo(highEFV, highEF, 0, lowEFV, lowEF, 0);
		}

	public int getQuantile(float exp) 
		{
		for(int i=0;i<100;i++)
			if(exp < quantiles[i])	return i;
		return 100;
		}
	
	public int getQuantile(float exp, int column) 
	{
	for(int i=0;i<100;i++)
		if(exp < this.quantileCols.get(column)[i])	return i;
	return 100;
	}

	public String getSynonym(String name, int nameType) {
		switch(nameType)
			{
			case ENTREZ:
				if(entrezgeneHash==null)	return null;
				return entrezgeneHash.get(name.replaceAll("\\(.*\\)", ""));
			case SYMBOL:
				if(external_gene_idHash==null)	return null;
				return external_gene_idHash.get(name.replaceAll("\\(.*\\)", ""));
			case ENSEMBL:
				if(ensembl_gene_idHash==null)	return null;
				//return ensembl_gene_idHash.get(name.replaceAll("\\(.*\\)", ""));
				return ensembl_gene_idHash.get(name);
			default:
				return name;
			}
		}
	
	public String getIdFromSynonym(String name, int nameType) {
		switch(nameType)
			{
			case ENTREZ:
				if(entrezgeneHash==null)	return null;
				return invertedHash(entrezgeneHash).get(name.replaceAll("\\(.*\\)", ""));
			case SYMBOL:
				if(external_gene_idHash==null)	return null;
				return invertedHash(external_gene_idHash).get(name.replaceAll("\\(.*\\)", ""));
			case ENSEMBL:
				if(ensembl_gene_idHash==null)	return null;
				return invertedHash(ensembl_gene_idHash).get(name.replaceAll("\\(.*\\)", ""));
			default:
				return name;
			}
		}
	
	
	public HashSet<String> getSynonyms(HashSet<String> names, int nameType)
		{
		HashSet<String> syns=new HashSet<String>();
		for(String n:names)
			{
			String s=getSynonym(n.toLowerCase(),nameType);
			if(s!=null)	syns.add(s);
			}
		return syns;
		}
	
	public HashSet<String> getIdsFromSynonyms(HashSet<String> names, int nameType)
		{
		HashSet<String> syns=new HashSet<String>();
		for(String n:names)
			{
			String s=getIdFromSynonym(n.toLowerCase(),nameType);
			if(s!=null)	syns.add(s);
			}
		return syns;
		}
	}