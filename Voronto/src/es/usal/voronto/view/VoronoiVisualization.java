package es.usal.voronto.view;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Rectangle;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.WindowListener;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.apache.commons.collections15.CollectionUtils;
import org.philhosoft.p8g.svg.P8gGraphicsSVG;

import java.util.Arrays;

import es.usal.voronto.control.Voronto;
import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.filters.TIFFileFilter;
import es.usal.voronto.model.filters.SVGFileFilter;
import es.usal.voronto.model.filters.TextFileFilter;
import es.usal.voronto.model.ontology.OntologyTerm;
import es.usal.voronto.model.services.david.ChartReportClient;
import es.usal.voronto.model.stats.Stats;
import es.usal.voronto.model.voronoi.Cell;
import es.usal.voronto.model.voronoi.BernhardtTessellation;
import processing.awt.PSurfaceAWT.SmoothCanvas;
import processing.core.*;
import sample.session.service.xsd.SimpleChartRecord;

public class VoronoiVisualization extends PApplet{

	/**
	 * 
	 */
	private static final long serialVersionUID = 2515531690539870969L;
	BernhardtTessellation v;
	int hoveredRegion=-1;
	int selected=-1;
	int levelThreshold=3;
	int selectedCol=0;
	private List<Cell> hoveredCells=Collections.synchronizedList(new ArrayList<Cell>());
	private List<Cell> searchedCells=Collections.synchronizedList(new ArrayList<Cell>());
	private int minHoveredLevel;
	public String message=null;
	
	public double minExp=0;
	public double maxExp=0;
	public double avgExp=0;
	public ExpressionData expData;
	public String imagePath;//path where image captures are saved by default
	public String textPath;//path where text results are saved by default
	
	
	//scale constants to get minExp, maxExp and avgExp
	public static final int SCALE_MATRIX=0;//from the whole expression matrix
	public static final int SCALE_CONDITION=1;//from the current column
	public static final int SCALE_ONTOLOGY=2;//from the elements of the matrix in the ontology 
	//public int SCALE_MODE=SCALE_CONDITION;
	public int SCALE_MODE=SCALE_MATRIX;
	
	//coloring constants for
	public static final int COLOR_EXPRESSION=0;//max or min
	public static final int COLOR_DEVIATION=1;//overall deviation respect to average
	public static final int COLOR_INTERNAL_DEVIATION=2;//overall deviation respect to average (of each term)
	public int COLOR_MODE=COLOR_EXPRESSION;
	private boolean SHOW_LABELS=true;
	
	private int contAvgExp=0;
	
	
	private int START_Y=20;//Beginning and end of Voronoi tessellation
	private int END_Y=500+START_Y;
	private Cell[] cells;
	//Boxes for each element
	Rectangle2D.Float conditionBox;
	Rectangle2D.Float organismBox;
	Rectangle2D.Float ontologyBox;
	Rectangle2D.Float scaleBox;
	Rectangle2D.Float cellBox;
	Rectangle2D.Float noteBox;
	
	private int hoveredBox=-1;
	private final int ORGANISM=0;
	private final int ONTOLOGY=1;
	private final int CONDITION=2;
	private final int CELL=3;
	private final int SCALE=4;
	public static final int KEGG=5;
	public static final int GO=6;
	public static final int SLIMBP=7;
	public static final int SLIMCC=8;
	public static final int BP=9;
	public static final int CC=10;
	public static final int MF=13;
	public static final int REACTOME=11;
	public static final int CUSTOM=12;
	
	public static final int MIN=0;
	public static final int MEAN=1;
	public static final int MEDIAN=2;
	//public int whiteValue=MEDIAN;
	public int whiteValue=MEAN;
	
	private int type;
	private boolean computedLabels=false;
	
	private TreeMap<OntologyTerm,BernhardtTessellation> tessellations;//relation of tessellations we do have, in the case of redefining the tessellation by getting into it
	private OntologyTerm root; //terms that makes the root of the current tessellation
	private TreeMap<OntologyTerm, TreeMap> map;
	private int maxDepth;
	private boolean drawSearch=false;
//	private String searchString="";
	public float[] maxSd=null;
	private String helpMessage="press 'h' for help";
	private String noteLabel=helpMessage;
	public Voronto voronto=null;
	private boolean saving;
	private double medianExp;

	private String ontologyName;
	public String customOntologyName;
	
	public PFont font;
	public int fontSize=14;
	private boolean ctrlDown=false;
	private boolean altDown=false;
	
	private SearchFrame searchFrame=null;
	private DifExpSearch deSearchFrame;
	private JFrame colorFrame=null;
	
	//private SearchFrame searchPApplet;
	private Cell minHoveredCell;
	private CellHeatmap gh;
	private int textSize=22;
	private Cell selectedCell;
	private int lineSpacing=0;
	private boolean skewed=true;
	private final int AVERAGE_HEATMAP=0;
	private final int PROFILE_LINE=1;
	private final int PROFILE_PLOT=2;
	private int profileType=AVERAGE_HEATMAP;
	public boolean drawOnlyScale=false;
	public float scaleExpression=-999999;
	public int scaleSample=-1;
	public int nameType=1;
	private boolean test=false;
	//private boolean test=true;	//for printing tests
	
	//To keep degs/enriched terms in case of pressing enter
	private HashSet<String> degsi;
	private SimpleChartRecord[] records;
	private int termThreshold;
	private boolean entrezConv;
	private String searchedText;
	private int profileMargin;
	private HashMap<String, OntologyTerm> enrichedTerms=new HashMap<String,OntologyTerm>();
	private double enrichmentThreshold;
	
	public HashMap<String,PFont> fonts=new HashMap<String, PFont>();
	private Color[] palette;
	//public PFont[] fonts=new PFont[28];
	//public ArrayList<PFont> fonts=new ArrayList<PFont>();
	
	public VoronoiVisualization()
	 	{
		 super();
	 	}
	 
	public void setup()
		{
		font=loadFont("es/usal/voronto/font/AppleSymbols-14.vlw");
		for(int i=8;i<28;i++)
    		fonts.put(""+i, loadFont("es/usal/voronto/font/AppleSymbols-"+i+".vlw"));
    	smooth();
		noLoop();
	 	}
	
	public VoronoiVisualization(Cell[] c, int width, int height)
		{
		super();
		v=new BernhardtTessellation(c, width, height);//, this);
		}
	
	public VoronoiVisualization(TreeMap<OntologyTerm,TreeMap> m) throws Exception
		{
		super();
	 	size(width, height);
		v=new BernhardtTessellation(m, null, 700, 500,/* this,*/ VoronoiVisualization.KEGG, 3);
		cells=v.getCells(); 
		for(Cell c:cells)	recursiveRegionTranslation(c, 0, START_Y);
		}
	
	public VoronoiVisualization(TreeMap<OntologyTerm,TreeMap> m, String customName, ExpressionData md, int width, int height, int type, int maxDepth, Voronto parent) throws Exception
	{
	super();
	this.type=type;
	this.width=width;
	this.height=height;
	this.setSize(width, height);
	this.END_Y=height-100+START_Y;
	this.maxDepth=maxDepth;
	this.voronto=parent;
	this.imagePath=null;
	this.textPath=null;
	
	Color[] basePalette={Color.BLUE, Color.WHITE, Color.RED};
	changeColors(basePalette);
	File file=new File(parent.imgPath);
	if (!file.exists()) 
	       new FileOutputStream(file).close();
	file=new File(parent.txtPath);
   	if (!file.exists()) 
   	       new FileOutputStream(file).close();
   	//try{
	BufferedReader br=new BufferedReader(new FileReader(parent.imgPath));
	this.imagePath=br.readLine();
	br=new BufferedReader(new FileReader(parent.txtPath));
	this.textPath=br.readLine();
	//}catch(Exception ex){System.out.println("pathFile "+parent.txtPath+" or "+parent.imgPath+" non existing yet");}
	
	if(this.imagePath==null || imagePath.length()==0)
		imagePath=md.filePath;
	if(this.textPath==null || textPath.length()==0)
		textPath=md.filePath;

	if(customName!=null)	customOntologyName=customName;
	
 	v=new BernhardtTessellation(m, md, width, height-100, type, maxDepth);
 	this.expData=md;
 	if(type==VoronoiVisualization.KEGG)	this.levelThreshold=v.maxLevel+1;
	cells=v.getCells();
	for(Cell c:cells)	recursiveRegionTranslation(c, 0, START_Y);
	
	tessellations=new TreeMap<OntologyTerm, BernhardtTessellation>();
	root=new OntologyTerm("root", "root");
	tessellations.put(root, v);
	map=m;
	
	setOntologyName();
	
	if(md!=null)	
		{
		mapExpression(md);
		expression2color(md);
		System.out.println("Finished mapping colors");
		}
	else
		expression2color(null);
	}
	
	public void setOntologyName()
		{
		if(root.name.equals("root") && root.id.equals("root"))
			{
			if(customOntologyName!=null && customOntologyName.length()>0)
				ontologyName=customOntologyName+" ("+this.levelThreshold+")";
			else if(type==KEGG)	ontologyName="KEGG ("+this.levelThreshold+")";
			else if(type==REACTOME)	ontologyName="REACTOME ("+this.levelThreshold+")";
			else if(type==SLIMBP)	ontologyName="GO Slim - BP ("+this.levelThreshold+")";
			else if(type==SLIMCC)	ontologyName="GO Slim - CC ("+this.levelThreshold+")";
			else if(type==BP)	ontologyName="GO - Biological Process ("+this.levelThreshold+")";
			else if(type==CC)	ontologyName="GO - Cellular Component ("+this.levelThreshold+")";
			else if(type==MF)	ontologyName="GO - Molecular Function ("+this.levelThreshold+")";
			else
				ontologyName="Ontology ("+this.levelThreshold+")";
			}
		else
			ontologyName=root.name+" ("+this.levelThreshold+")";
		}
	
	public VoronoiVisualization(int width, int height)
		{
		super();
		}
	
	public void setOntology(TreeMap<OntologyTerm, TreeMap> m) throws Exception
		{
		v=new BernhardtTessellation(m, null, width, height, /*this, */VoronoiVisualization.KEGG, 3);
		}
	
	
	/**
	 * Returns the width in pixels of a text string (cad), drawn with the given font (f)
	 * @param cad
	 * @param f
	 * @return
	 */
	public float getWidth(String cad, PFont f)
		{
		float width=0;
		for(int i=0;i<cad.length();i++)
			{
			width+=f.width(cad.charAt(i));
			}
		return width*f.getSize();
		}

	public synchronized void draw() 
		{
		noStroke();
		fill(255);
		rect(0,0,width,height);
		noFill();
		
		//drawing labels
		textFont(font, fontSize);//TODO: check about fonts
		
		fill(0,0,0);
		if(expData!=null)
			{
			textAlign(LEFT, CENTER);
			text(expData.organism, 10, 8);
			organismBox=new Rectangle2D.Float(10,8,getWidth(expData.organism, font), fontSize);
			
			textAlign(RIGHT, CENTER);
			text(expData.getColumnLabel(selectedCol), width-10, 8);
			conditionBox=new Rectangle2D.Float(width-10-getWidth(expData.getColumnLabel(selectedCol), font),8, getWidth(expData.getColumnLabel(selectedCol), font), fontSize);
			
			if(!skewed)	drawScale();
			else		drawSkewedScale();
			}
		
		//Ontology label
		fill(0,0,0);
		textAlign(CENTER, CENTER);
		text(ontologyName, (int)(width*0.5), 8);	
		ontologyBox=new Rectangle2D.Float((int)(width*0.5-getWidth(ontologyName, font)*0.5),8, getWidth(ontologyName, font), font.getSize());
		
		//Note label
		if(!saving)
			{
			textAlign(LEFT, CENTER);
			fill(200,200,200);
			text(noteLabel, 10, this.END_Y+45);
			noFill();
			}
		
		textFont(font,fontSize);
		textAlign(LEFT, CENTER);
	
		//drawing regions
		if(!computedLabels)
			{
			long t0=System.currentTimeMillis();
			for(Cell c:cells)	
				recursiveLabelComputing(c,0);	
			System.out.println("Time in computing labels:" +(System.currentTimeMillis()-t0));
			computedLabels=true;
			}
		
		for(Cell c:cells)	recursiveRegionDrawing(c, 0);
		
		//drawing searched regions
		for(Cell s:searchedCells)
			{
			//color
			if(COLOR_MODE==COLOR_EXPRESSION)
				{
				if(s.searchedColor!=null)	stroke(s.searchedColor.getRed(), s.searchedColor.getGreen(), s.searchedColor.getBlue());
				else						stroke(0,200,0);
				}
			else								stroke(200,0,0);
			
			//stroke
			if(!s.term.significant && (s.term.degs==null || s.term.degs.size()==0))	
					strokeWeight(getWidth(s)+2+Math.min(10, s.term.numRelevantSubterms));
			else			//there are degs or enrichment
				{
				if(!s.term.significant && s.term.degs.size()>0)	//degs
						strokeWeight(getWidth(s)+2+(int)(10.0*s.term.degs.size()/s.term.geneExs.keySet().size()));
				else											//enrichment
					strokeWeight(getWidth(s)+2+Math.min(6,-(int)(Math.log10(s.term.pvalueBH))));
				}
				
				
			noFill();
			if(s.region!=null)	
				{
				if(s.level<=this.levelThreshold)
					{
					boolean draw=true;
					if(s.subcells!=null)
						for(Cell s2:s.subcells)
							{
							if(searchedCells.contains(s2) && levelThreshold>=s2.level && s2.region!=null)
								{
								draw=false;
								break;
								}
							}
						
					if(draw)	s.region.draw(this);
					}
				}
			strokeWeight(0);
			}
		
		//drawing hovered region
		synchronized(hoveredCells)
			{
			for(Cell h:hoveredCells)
				{
				if(h.level<=levelThreshold)
					{
					if(h.level!=min(levelThreshold, minHoveredLevel))
						{
						stroke(0,0,0);
						strokeWeight(getWidth(h)+1);
						
						noFill();
						h.region.draw(this);
						strokeWeight(0);
						}
					else
						{
						noFill();
						strokeWeight(getWidth(h)+1);
						this.stroke(0,0,0);
						h.region.draw(this);
						fill(0,0,0);
						strokeWeight(0);
						
						textAlign(LEFT, TOP);
						textFont(font,fontSize);
						
						String cad=h.term.name+" ("+(int)(h.numLeaves)+")";
						text(cad, 10, (int)(END_Y+font.getSize()*0.5));
						if(h.term.degs!=null && h.term.degs.size()>0)	
							{
							if(searchedCells.contains(h))	
								{
								if(h.searchedColor!=null)
									fill(h.searchedColor.getRed(), h.searchedColor.getGreen(), h.searchedColor.getBlue());
								else
									fill(0,200,0);
								}
							else							fill(0);
								
							if(h.term.significant || (h.term.degs!=null && h.term.degs.size()>0) || h.term.numRelevantSubterms>0) //if significant, it will be drawn in green, if not but has something, we show also it in black
								{
								textFont(font,12);
								DecimalFormat df=new DecimalFormat("#.#");
								DecimalFormat df2=new DecimalFormat("0.##E0");
								cad="  ";
								if(h.term.pvalue<1)
									cad=cad+"pval="+df2.format(h.term.pvalueBH)+"  ";
								if(h.term.degs.size()>0)
									cad=cad+"  "+h.term.degs.size()+" relevant genes ("+df.format((h.term.degs.size()*100.0)/h.numLeaves)+"%)";
								if(h.term.numRelevantSubterms>0)
									cad=cad+"  "+h.term.numRelevantSubterms+" relevant subterm(s)";
								text(cad, 10, (int)(END_Y+font.getSize()*1.5));
								}
							fill(255);
							}
						fill(255);
						
						switch(profileType)
							{
							case AVERAGE_HEATMAP:
								drawProfile(h);
								break;
							case PROFILE_LINE:
							default:
								drawProfileLine(h);
								break;
							case PROFILE_PLOT:
								drawProfilePlot(h);
								break;
							}
						}
					}
				}
			if(hoveredCells.size()==0 && selectedCell==null)
				{
				textFont(font,fontSize);
				textAlign(LEFT, TOP);
				fill(0,0,0);
				text("No term selected", 10, (int)(END_Y+font.getSize()*0.5));
				noFill();
				cellBox=new Rectangle2D.Float(10,this.END_Y+8, getWidth("No term selected", font), font.getSize());
				}
			}
		
		if(selectedCell!=null)
			{
			noFill();
			strokeWeight(getWidth(selectedCell)+3);
			this.stroke(254,227,13);
			selectedCell.region.draw(this);
			fill(0,0,0);
			strokeWeight(0);
			
			if(minHoveredCell==null)
				{
				textAlign(LEFT, TOP);
				textFont(font,fontSize);
			
				String cad=selectedCell.term.name+" ("+(int)(selectedCell.numLeaves)+")";
				text(cad, 10, (int)(END_Y+font.getSize()*0.5));
				if(selectedCell.term.degs!=null && selectedCell.term.degs.size()>0  && searchedCells.contains(selectedCell))	
					{
					fill(0,150,0);
					
					//if(selectedCell.term.pvalueBH!=-1)
					if(selectedCell.term.significant)
							{
						textFont(font,12);
						DecimalFormat df=new DecimalFormat("#.#");
						DecimalFormat df2=new DecimalFormat("0.##E0");
						cad="  BH="+df2.format(selectedCell.term.pvalueBH)+"  "+selectedCell.term.degs.size()+" relevant genes ("+df.format((selectedCell.term.degs.size()*100.0)/selectedCell.numLeaves)+"%)";
						text(cad, 10, (int)(END_Y+font.getSize()*1.2));
						}
					else
						{
						textFont(font,12);
						DecimalFormat df=new DecimalFormat("#.#");
						cad="  "+selectedCell.term.degs.size()+" relevant genes ("+df.format((selectedCell.term.degs.size()*100.0)/selectedCell.numLeaves)+"%)";
						text(cad, 10, (int)(END_Y+font.getSize()*1.2));
						}
					fill(255);
					}
			
				fill(255);				
				switch(profileType)
					{
					case AVERAGE_HEATMAP:
						drawProfile(selectedCell);
						break;
					case PROFILE_LINE:
					default:
						drawProfileLine(selectedCell);
						break;
					case PROFILE_PLOT:
						drawProfilePlot(selectedCell);
						break;
					}
				}
			}
				
		
		//drawing help messages
		switch(hoveredBox)
			{
			case ORGANISM:
				drawNote("Organism of the expression data", organismBox, font);
				break;
			case ONTOLOGY:
				drawNote("Ontology in the tessellation.\nEach cell represents an ontology term.\nCell size approximates to its number of genes\nHigh level terms have wider edges.\nChange depth level with up/down arrow keys.\nThe (number) is the max depth visualized", ontologyBox, font);
				break;
			case CONDITION:
				drawNote("Expression for this condition is color-mapped.\nBrowse conditions with left/right arrow keys", conditionBox, font);
				break;
			case CELL:
				if(type==KEGG)	drawNote("Term names appear here when hovered.\nIn parenthesis, the #�of genes in the term.\nCells' color is the mean gene expression level.\nA cell profile bar is shown close to this name.\nPress 'l' to show/hide fitting labels\nAlt-click on a cell to open the browser\nwith its colored KEGG pathway.\nDouble-click on a cell to show its heatmap.", cellBox, font);
				else			drawNote("Term names appear here when hovered.\nIn parenthesis, the #�of genes in the term.\nCells' color is the mean gene expression level.\nA cell profile bar is shown close to this name.\nPress 'l' to show/hide fitting labels\nAlt-click on a cell to open the browser\nwith its ontology web entry.\nDouble-click on a cell to show its heatmap.", cellBox, font);
				break;
			case SCALE:
				if(COLOR_MODE==COLOR_EXPRESSION)	
					{
					if(whiteValue==MEDIAN)	drawNote("Quantile color scale\nPress 's' to change scale type\n", scaleBox, font);
					if(whiteValue==MEAN)	drawNote("Numerical two-color scale\nPress 's' to change scale type\n", scaleBox, font);
					if(whiteValue==MIN)		drawNote("Numerical one-color scale\nPress 's' to change scale type\n", scaleBox, font);
					
					
					}
				else						drawNote("Color scale for the mapping:\n  -Bright green/white for max to no deviation.\nUse 'c' key to change the scale\n", scaleBox, font);
				break;
			default:
				break;
			}
		
		//System.out.println("Time in post-region drawing: "+(System.currentTimeMillis()-t0)/1000.0);
		}

public void drawNote(String text, Rectangle2D.Float bounds, PFont font)
	{
	int margin=5;
	int noteWidth=250;
	int noteHeight=fontSize;
	String temp=text;
	while(temp.indexOf("\n")!=-1)	
		{
		noteHeight+=font.getSize()*1.2+margin; 
		temp=temp.substring(temp.indexOf("\n")+1);
		}
	int noteX=(int)bounds.getX();
	if(bounds.x+noteWidth+margin>this.width) noteX=this.width-noteWidth-margin;	
	int noteY=(int)bounds.getY();
	if(bounds.y+noteHeight+margin+font.getSize()>this.height) noteY=this.height-noteHeight-margin-fontSize;	
	
	fill(250,250,200);
	stroke(100,100,100);
	rect(noteX, noteY+font.getSize(), noteWidth, noteHeight);
	fill(100,100,100);
	textFont(font, fontSize);
	textAlign(LEFT);
	text(text, noteX+margin, (int)(noteY+fontSize*1.5+margin));
	noFill();
	noStroke();
	return;
	}

/**
 * Draws profile for the given cell
 */
public void drawDetailedProfile(Cell c)
	{
	if(expData==null)	return;
	
	int squareSize=12;
	
	strokeWeight(1);
	stroke(150);
	for(int i=0;i<expData.getNumConditions();i++)
		{
		//draw empty square
		strokeWeight(1);
		stroke(150);
		fill(200);
		rect((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+i*squareSize, this.END_Y+20,squareSize,squareSize); 
		
		//fill with genes colors
		noStroke();
		int miniSize=(int)Math.max(1,((squareSize-1)/Math.ceil(Math.sqrt(c.term.geneExs.size()))));
		int contX=0;
		int contY=0;
		Iterator<String> it=c.term.geneExs.keySet().iterator();
		for(int j=0;j<Math.min((squareSize-1)*(squareSize-1), c.term.geneExs.size());j++)
			{
			String gene=it.next();
			Color co=getColor(c.term.geneExs.get(gene), i);
			fill(co.getRed(), co.getGreen(), co.getBlue());
			
			rect((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+i*squareSize+1+contX*miniSize, this.END_Y+20+1+contY*miniSize,miniSize,miniSize); 
			
			contX++;
			if((contX+1)*miniSize>squareSize-1)	{contX=0;contY++;}
			}
		}
	
	strokeWeight(2);
	stroke(0);
	line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize, this.END_Y+20+squareSize,(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize+squareSize,this.END_Y+20+squareSize);
	}

/**
 * Returns the color for the given expression array, corresponding to the given column
 * @param exs
 * @param col
 * @return
 */
public Color getColor(float ex, int col)
	{
	setScale(expData, col);
	int h=-1;
	switch(whiteValue)
		{
		case MEAN:	//raw coloring (two color scale)
			switch(COLOR_MODE)
				{
				/*case VoronoiVisualization.COLOR_EXPRESSION:
					if(ex>avgExp)
						{
						h=(int)Math.round(255-((ex-avgExp)/(maxExp-avgExp)*255));
						return new Color(255,h,h);
						}
					else
						{
						h=(int)Math.round(255-(Math.abs(ex-avgExp)/Math.abs(minExp-avgExp))*255);
						return new Color(h,h,255);
						}*/
				case VoronoiVisualization.COLOR_EXPRESSION:
						return palette[(int)Math.round((ex-minExp)/(maxExp-minExp)*510)];
				case VoronoiVisualization.COLOR_DEVIATION:
					h=(int)Math.round(255-(Math.abs(ex-avgExp)/Math.max(Math.abs(avgExp-minExp), Math.abs(avgExp-maxExp)))*255);
					return new Color(h,255,h);
				case VoronoiVisualization.COLOR_INTERNAL_DEVIATION://TODO: testing
					return null;
				}
			break;
		case MEDIAN:	//quantile coloring (two color scale) 
			switch(COLOR_MODE)
				{
				case VoronoiVisualization.COLOR_EXPRESSION:
					int q=expData.getQuantile(ex);
					if(q>=50)
						{
						h=(int)Math.round(255-((q-50.0)/50)*255);
						//return new Color(255,h,h);
						return palette[Math.abs(h-255)+254];
						//return palette[h+254];
						}
					else
						{
						h=(int)Math.round(255-((50.0-q)/50)*255);
						//return new Color(h,h,255);
						return palette[h];
						}
				case VoronoiVisualization.COLOR_DEVIATION:
					System.err.println("Option not supported for quantiles");
					return null;
				case VoronoiVisualization.COLOR_INTERNAL_DEVIATION://TODO: testing
					System.err.println("Option not supported for quantiles");
					return null;
				}
			break;
		case MIN:
		default:
			switch(COLOR_MODE)
				{
				case COLOR_EXPRESSION:
					h=(int)Math.round(255-((ex-minExp)/(maxExp-minExp)*255));
					return new Color(255, h, h);					
				case COLOR_DEVIATION:
					h=(int)Math.round(255-(Math.abs(ex-minExp)/(maxExp-minExp)*255));
					return new Color(255, h, h);
				case COLOR_INTERNAL_DEVIATION://TODO: testing
					/*if(cell.expressionDeviation.size()>0)	
						h=(int)(255-(cell.expressionDeviation.get(column)/maxSd[column])*255);
					else	h=255;
					cell.color.set(column, new Color(h,255,h));*/
					break;
				}break;
		}
	return null;
	}

public Color getColor(ArrayList<Float> exs, int col)
	{
	if(exs==null)	return null;
	return getColor(exs.get(col), col);
	}

public void drawProfile(Cell c)
	{
	if(expData==null)	return;
	
	int squareSize=12;
	profileMargin=25;
	strokeWeight(1);
	stroke(150);
	for(int i=0;i<expData.getNumConditions();i++)
		{
		if(selectedCol!=i)
			{
			Color co=c.color.get(i);
			fill(co.getRed(), co.getGreen(), co.getBlue());
			rect((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+i*squareSize, this.END_Y+profileMargin,squareSize,squareSize);
			}
		}
	
	if(selectedCol>=0)
		{
		Color co=c.color.get(selectedCol);
		fill(co.getRed(), co.getGreen(), co.getBlue());
		rect((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize, this.END_Y+profileMargin,squareSize,squareSize);
		if(test==false)
			{
			strokeWeight(2);
			stroke(0);
			line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize, this.END_Y+profileMargin+squareSize+3,(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize+squareSize,this.END_Y+profileMargin+squareSize+3);
			}
		}
	}

public void drawProfileLine(Cell c)
	{
	if(expData==null)	return;
	
	int squareSize=12;
	
	strokeWeight(1);
	stroke(220);
	int em=squareSize*2;
	int eM=0;
	int eAvg=(int)(squareSize*2-squareSize*2*(expData.average-expData.min)/(expData.max-expData.min));
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)em, 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)em);
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)eM, 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)eM);
	if(whiteValue==MEAN)
		{
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)eAvg, 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)eAvg);
		}
	else
		{
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)(em*0.5), 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)(em*0.5));
		}
	
	strokeWeight(2);
	stroke(100);
	for(int i=0;i<expData.getNumConditions()-1;i++)
		{
		double e1;
		double e2;
		if(whiteValue==MEAN)
			{
			e1=squareSize*2-squareSize*2*(c.expressionLevel.get(i)-expData.min)/(expData.max-expData.min);
			e2=squareSize*2-squareSize*2*(c.expressionLevel.get(i+1)-expData.min)/(expData.max-expData.min);
			}
		else
			{
			e1=squareSize*2-squareSize*2*(0.01*expData.getQuantile(c.expressionLevel.get(i)));
			e2=squareSize*2-squareSize*2*(0.01*expData.getQuantile(c.expressionLevel.get(i+1)));
			}
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+i*squareSize, this.END_Y+profileMargin+(int)e1, 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(i+1)*squareSize, this.END_Y+profileMargin+(int)e2);
		}
	
	stroke(0);
	fill(0);
	if(selectedCol>=0 && test==false)
		{
		double e;
		if(whiteValue==MEAN)
			e=squareSize*2-squareSize*2*(c.expressionLevel.get(selectedCol)-expData.min)/(expData.max-expData.min);
		else
			e=squareSize*2-squareSize*2*(0.01*expData.getQuantile(c.expressionLevel.get(selectedCol)));
		ellipse((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+selectedCol*squareSize, this.END_Y+profileMargin+(int)e,2,2);
		}
	}

public void drawProfilePlot(Cell c)
	{
	if(expData==null || c==null)	return;
	
	int squareSize=12;
	strokeWeight(1);
	stroke(220);
	int em=squareSize*2;
	int eM=0;
	int eAvg=(int)(squareSize*2-squareSize*2*(expData.average-expData.min)/(expData.max-expData.min));
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)em, 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)em);
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)eM, 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)eM);
	if(whiteValue==MEAN)
		{
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)eAvg, 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)eAvg);
		}
	else
		{
		line((int)(width*0.5-expData.getNumConditions()*squareSize*0.5), this.END_Y+profileMargin+(int)(em*0.5), 
				(int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+(expData.conditionNames.length-1)*squareSize, this.END_Y+profileMargin+(int)(em*0.5));
		}
	
	//if(minHoveredCell.term.geneExs==null || minHoveredCell.term.geneExs.values().iterator().next().size()!=expData.getNumConditions())
	if(c.term.geneExs==null || c.term.geneExs.values().iterator().next().size()!=expData.getNumConditions())
			c.computeExpressionProfile(expData);

	noStroke();
	fill(100);
	for(int i=0;i<expData.getNumConditions();i++)
		{
		int rnd=-2;
		Iterator<String> it=c.term.geneExs.keySet().iterator();
		for(int j=0;j<c.term.geneExs.size();j++)
			{
			double e;
			ArrayList<Float> p=c.term.geneExs.get(it.next());
			if(p==null || p.size()==0)
				{
				return;
				}
			if(whiteValue==MEAN)
				{
				e=squareSize*2-squareSize*2*(p.get(i)-expData.min)/(expData.max-expData.min);
				}
			else
				{
				e=squareSize*2-squareSize*2*(0.01*expData.getQuantile(p.get(i)));
				}
			
			if(selectedCol==i && !test)
				{
				fill(0);
				ellipse((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+i*squareSize+rnd, this.END_Y+profileMargin+(int)e,1,1);
				}
			else
				{
				fill(150);
				ellipse((int)(width*0.5-expData.getNumConditions()*squareSize*0.5)+i*squareSize+rnd, this.END_Y+profileMargin+(int)e,1,1);
				}
				
			rnd++;
			if(rnd>4)	rnd=-2;
			}
		}
	}

public void drawSkewedScale()
	{
	Cell hc=null;
	Cell finalhc=null;
	if(drawOnlyScale)
		{
		drawSkewedScale(this.scaleExpression, this.scaleSample);
		drawOnlyScale=false;
		}
	else
		{
		for(int i=0;i<hoveredCells.size();i++)
			{
			hc=hoveredCells.get(i);
			if(hc.level==min(levelThreshold, minHoveredLevel))
				finalhc=hc;
			}
		hc=finalhc;
		if(hc==null && this.selectedCell!=null)
			hc=selectedCell;
		
		if(hc!=null)
			drawSkewedScale(hc.getExpression(selectedCol), selectedCol);
		else
			drawSkewedScale((float)minExp-100, selectedCol);
		}

	}

public void drawSkewedScale(float expressionLevel, int sample)
	{
	textFont(font, 12);
	
	DecimalFormat df=new DecimalFormat("#.##");
	textAlign(CENTER, TOP);
	int scaleWidth=150;
	int scaleHeight=10;
	int scaleX=width-scaleWidth-15;
	int scaleY=(int)(END_Y+scaleHeight*1.7);
	scaleBox=new Rectangle2D.Float(scaleX, scaleY, scaleWidth, scaleHeight);
	strokeWeight(1);
	
	Color selCol=null;
	if(expressionLevel>=minExp) selCol=this.getColor(expressionLevel, sample);

	
	boolean refDrawn=false;
	double jump0=(255-(0/(float)scaleWidth)*255);
	double jump1=(255-(0/(float)scaleWidth)*255);
	double jump=Math.ceil(jump1-jump0);
	
	
	double mid=scaleWidth*0.5;//for median
	if(whiteValue==MEAN)	
		mid=scaleWidth*(expData.average-expData.min)/(expData.max-expData.min);
	else if (whiteValue==MIN)
		mid=minExp;
	
	double jumpBelow=Math.ceil(255/(255*mid/scaleWidth));
	double jumpAbove=Math.ceil(255/Math.abs(255*mid/scaleWidth-255));
	
	//Draw legend
	stroke(200,200,200);
	fill(200);
	line(scaleX,scaleY,scaleX,scaleY-3);
	if(whiteValue==MEAN || whiteValue==MIN)	
		{
		text(df.format(expData.min),  scaleX, scaleY-12);
		}
	else if(whiteValue==MEDIAN)	text("0", scaleX, scaleY-12);
	
	line((int)(scaleX+mid),scaleY+scaleHeight,(int)(scaleX+mid),scaleY-3);
	if(whiteValue==MEAN || whiteValue==MIN)		
		{
		text(df.format(expData.average), (int)(scaleX+mid), scaleY-12);
		}
	else if(whiteValue==MEDIAN)	text("50", (int)(scaleX+mid), scaleY-12);
	
	line(scaleX+scaleWidth-1,scaleY+scaleHeight,scaleX+scaleWidth-1,scaleY-3);
	if(whiteValue==MEAN || whiteValue==MIN)		
		{
		text(df.format(expData.max),  scaleX+scaleWidth-1, scaleY-12);
		}
	else if(whiteValue==MEDIAN)	text("100", scaleX+scaleWidth-1, scaleY-12);
	
	line((int)(scaleX+mid*0.5),scaleY+scaleHeight,(int)(scaleX+mid*0.5),scaleY-2);
	line((int)(scaleX+mid*2),scaleY+scaleHeight,(int)(scaleX+mid*2),scaleY-2);
	for(int i=0;i<scaleWidth;i++)
		{
		int h=0;
		boolean drawReference=false;
		switch(COLOR_MODE)
			{
			case COLOR_EXPRESSION:
				if(i>mid)
					{
					h=(int)Math.round(255-((i-mid)/Math.abs(scaleWidth-mid)*255));
					stroke(255,h,h);
					if(selCol!=null && selCol.getRed()==255 && Math.abs(selCol.getBlue()-h)<=jumpAbove)
						drawReference=true;
					}
				else	
					{
					h=(int)Math.round(255-(Math.abs(i-mid)/Math.abs(mid))*255);
					stroke(h,h,255);
					if(selCol!=null && selCol.getBlue()==255 && Math.abs(selCol.getRed()-h)<=jumpBelow)
							drawReference=true;
					}
				break;
			case COLOR_DEVIATION:
				h=(int)Math.round(255-(i/(float)scaleWidth)*255);
				stroke(h,255,h);
				if(selCol!=null && Math.abs(selCol.getBlue()-h)<=jump)
						drawReference=true;
				break;
			}
		line(scaleX+i,scaleY,scaleX+i,scaleY+scaleHeight);
		if(drawReference && !refDrawn)
			{
			stroke(0);
			fill(0);
			line(scaleX+i, scaleY-2, scaleX+i, scaleY+scaleHeight+2);
			fill(255,200);
			noStroke();
			rect(scaleX+i-10, scaleY+scaleHeight, 20,15);
			fill(0);
			stroke(0);
			if(whiteValue==MEAN || whiteValue==MIN)		text(df.format(expressionLevel),  scaleX+i, scaleY+scaleHeight+3);
			else if(whiteValue==MEDIAN)
				text(df.format(expData.getQuantile(expressionLevel)),  scaleX+i, scaleY+scaleHeight+3);
			refDrawn=true;
			}
		}
	stroke(200,200,200);
	noFill();
	rect(scaleX-1,scaleY-1,scaleWidth+1,scaleHeight+1);
	
	textFont(font, fontSize);
	
	noStroke();
	}


/**
 * Draws a scale bar
 */
public void drawScale()
	{
	int sw=150;
	int sh=10;
	int x=width-sw-15;
	int y=(int)(END_Y+sh*0.5);
	scaleBox=new Rectangle2D.Float(x, y, sw, sh);
	strokeWeight(1);
	
	Cell hc=null;
	Cell finalhc=null;
	for(int i=0;i<hoveredCells.size();i++)
		{
		hc=hoveredCells.get(i);
		if(hc.level==min(levelThreshold, minHoveredLevel))
			finalhc=hc;
		}
	hc=finalhc;
	if(hc==null && this.selectedCell!=null)
		hc=selectedCell;
	
	boolean refDrawn=false;
	double jump0=(255-(0/(float)sw)*255);
	double jump1=(255-(0/(float)sw)*255);
	double jump=Math.ceil(jump1-jump0)+2;
	
	for(int i=0;i<sw;i++)
		{
		int h=0;
		boolean drawReference=false;
		switch(COLOR_MODE)
			{
			case COLOR_EXPRESSION:
				if(i>sw*0.5)
					{
					h=(int)Math.round(255-((i-sw*0.5)/(sw*0.5)*255));
					stroke(255,h,h);
					if(hc!=null && hc.color.get(selectedCol).getRed()==255 && Math.abs(hc.color.get(selectedCol).getBlue()-h)<=jump)
						drawReference=true;
					}
				else	
					{
					h=(int)Math.round(255-(Math.abs(i-sw*0.5)/Math.abs(sw*0.5))*255);
					stroke(h,h,255);
					if(hc!=null && hc.color.get(selectedCol).getBlue()==255 && Math.abs(hc.color.get(selectedCol).getRed()-h)<=jump)
						drawReference=true;
					}
				break;
			case COLOR_DEVIATION:
				h=(int)Math.round(255-(i/(float)sw)*255);
				stroke(h,255,h);
				if(hc!=null && Math.abs(hc.color.get(selectedCol).getBlue()-h)<=jump)
					drawReference=true;
				break;
			}
		line(x+i,y,x+i,y+sh);
		if(drawReference && !refDrawn)
			{
			stroke(0);
			line(x+i, y-2, x+i, y+sh+2);
			refDrawn=true;
			}
		}
	stroke(200,200,200);
	noFill();
	rect(x-1,y-1,sw+1,sh+1);
	
	//Draw legend
	textAlign(CENTER, TOP);
	fill(200);
	line(x,y+sh,x,y+sh+2);
	if(whiteValue==MEAN)	text("min", x, y+sh+3);
	else if(whiteValue==MEDIAN)	text("0", x, y+sh+3);
	line((int)(x+sw*0.5),y+sh,(int)(x+sw*0.5),y+sh+2);
	if(whiteValue==MEAN)		
		{
		text("avg", (int)(x+sw*0.5), y+sh+3);
		text("("+expData.average+")", (int)(x+sw*0.5), y+sh+6);
		}
	else if(whiteValue==MEDIAN)	text("50", (int)(x+sw*0.5), y+sh+3);
	
	line(x+sw-1,y+sh,x+sw-1,y+sh+2);
	if(whiteValue==MEAN)		text("max", x+sw-1, y+sh+3);
	else if(whiteValue==MEDIAN)	text("100", x+sw-1, y+sh+3);
	
	line((int)(x+sw*0.1),y+sh,(int)(x+sw*0.1),y+sh+1);
	line((int)(x+sw*0.2),y+sh,(int)(x+sw*0.2),y+sh+1);
	line((int)(x+sw*0.3),y+sh,(int)(x+sw*0.3),y+sh+1);
	line((int)(x+sw*0.4),y+sh,(int)(x+sw*0.4),y+sh+1);
	line((int)(x+sw*0.6),y+sh,(int)(x+sw*0.6),y+sh+1);
	line((int)(x+sw*0.7),y+sh,(int)(x+sw*0.7),y+sh+1);
	line((int)(x+sw*0.8),y+sh,(int)(x+sw*0.8),y+sh+1);
	line((int)(x+sw*0.9),y+sh,(int)(x+sw*0.9),y+sh+1);
	noStroke();
	}


/**
 * Recursive browsing of the ontology hierarchy, drawing each cell following the tessellation placement and expression level
 * @param cell
 * @param level
 */
public void recursiveRegionDrawing(Cell cell, int level)
	{
	if(cell.region!=null && level<levelThreshold)
		{
		//1) Draw region
		this.stroke(120);
		strokeWeight(getWidth(cell));
		
		if(cell.color!=null && cell.color.size()>0)	fill(cell.color.get(selectedCol).getRGB());
		else										fill(240,240,240);
		cell.region.draw(this.g); // draw this shape
		
		if(SHOW_LABELS && cell.labelSize>0)
			{
			if(cell.color.get(this.selectedCol).getRed()<50 || cell.color.get(this.selectedCol).getBlue()<50 )	fill(255);
			else	fill(0);
			
			textAlign(LEFT, TOP);
			textFont(fonts.get(""+cell.labelSize));
			textLeading(cell.labelSize+lineSpacing);
			text(cell.label, cell.labelX, cell.labelY);
			textAlign(LEFT, BASELINE);
			}
		
		//2) Continue with recursion
		if(cell.subcells!=null && cell.subcells.length>0)	
			for(Cell cc:cell.subcells)
				recursiveRegionDrawing(cc, level+1);
		}
	}

public int getWidth(Cell cell)
	{
	//return Math.min(0, this.maxDepth*2-cell.level*2);
	return Math.max(0, 6-cell.level*2);
	}

/**
* Recursive browsing of the ontology hierarchy, computing label positions for each cell
* @param cell
* @param level
*/
public void recursiveLabelComputing(Cell cell, int level)
	{
	if(cell.region!=null)
		{
		//computeLabelPosition(cell);
		computeAndSplitLabelPosition(cell);
			
		//2) Continue with recursion
		if(cell.subcells!=null && cell.subcells.length>0)	
			for(Cell cc:cell.subcells)
				recursiveLabelComputing(cc, level+1);
			}
	}

public int numLines(String cad)
	{
	int numLines=1;
	String temp=cad;
	while(temp.indexOf("\n")>0)	{temp=temp.substring(temp.indexOf("\n")+1); numLines++;}	
	return numLines;
	}



public void computeAndSplitLabelPosition(Cell cell)
	{
	Rectangle r=cell.region.getPolygon().getBounds();
	int size=textSize;
	textFont(font, size);
	
	boolean fit=false;
	double diff=10000000;
	int bestSize=size;
	int bestX=-1, bestY=-1;
	
	while(!fit)
		{
		textSize(size);
		float tw=textWidth(cell.label);
		//NOTE: seems that textWidth and textFont have small discrepancies related to space widths. I make a small adjustment here based on that:
		float twu=tw/cell.label.length();//width per character
		tw+=twu*(cell.label.length()-cell.label.replaceAll(" ", "").length());//add the widht of spaces twice for computations here.
		
		int y=r.y+15;
		
		int ss=(int)(textAscent()+textDescent()+lineSpacing)*numLines(cell.label);
		
		while(y+ss<r.y+r.height)
			{
			int x=r.x+5;
			while(x+tw<r.x+r.width)
				{
				Rectangle rt=new Rectangle(x-2,y-2,(int)Math.ceil(tw)+4,ss+4);
				double newDist=Point2D.distance((double)cell.centroid[0], (double)cell.centroid[1], (double)x+tw*0.5, (double)y+ss*0.5);
				if(cell.region.getPolygon().contains(rt) &&	//here may fit
					diff>newDist)
					{
					diff=newDist;
					bestX=x; bestY=y; bestSize=size;
					fit=true;
					}
				else
					{
					x+=1;
					}
				}
			y+=1;
			}
		size--;
		
		textFont(font, size);
				
		if(fit==true)//return fitted label
			{
			cell.labelX=bestX;
			cell.labelY=bestY;
			cell.labelSize=bestSize;
			return;
			}
		else if(size<=8 && cell.label.indexOf(" ")>0 )//split
			{
			//instead of returning, try splitting the term by the mid point (if spaces are present)
			String temp=cell.label;
			int begin=0;
			ArrayList<Integer> spaces=new ArrayList<Integer>();
			while(temp.indexOf(" ")!=-1)	{spaces.add(begin+temp.indexOf(" ")); begin+=temp.indexOf(" ")+1; temp=temp.substring(temp.indexOf(" ")+1); }
			int middleSpace=-1;
			double minDist=99999999;
			for(int i=0;i<spaces.size();i++)
				{
				int space=spaces.get(i);
				double dist=Math.abs(space-cell.label.length()*0.5);
				if(dist<minDist)
					{minDist=dist; middleSpace=space;}
				}
			cell.label=cell.label.substring(0, middleSpace)+"\n"+cell.label.substring(middleSpace+1);
			computeAndSplitLabelPosition(cell);
			return;
			}
		else if(size<=8 && cell.label.indexOf(" ")==-1)//no fit and we cannote make it smaller (and readable)
			return;
		}
	}

/**
 * Displaces every region in the tessellation by (x,y)
 * @param cell
 * @param x
 * @param y
 */
public void recursiveRegionTranslation(Cell cell, int x, int y)
	{
	cell.translate(x,y);
	
	//2) Continue with recursion
	if(cell.subcells!=null && cell.subcells.length>=1)//TODO: not sure why but if I do >0 the sections move...
		{
		for(Cell cc:cell.subcells)
			recursiveRegionTranslation(cc, x,y);
		}
	}

/**
 * When the mouse id double-clicked, the web browser is opened with the corresponding KEGG pathway, colored.
 */
//by gene ids
public void mouseReleased() {
    // do something based on mouse movement
	// update the screen (run draw once)
	if(mouseEvent.getClickCount()==1 && hoveredCells.size()>0)
		{
		//if(mouseEvent.isAltDown())
		if(mouseEvent.isControlDown())
		//if(altDown)	//show webpage
			{
			selectedCell=null;
			if(type==KEGG)
				{
				//Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
				//setCursor(hourglassCursor);
				this.getSurface().setCursor(Cursor.WAIT_CURSOR);
				
				for(Cell c:hoveredCells)
					{
					if(c.subcells==null || c.subcells.length==0)
						{
						try{
							if(expData!=null)
								{//
								String cad="";
								
								//NOTE: this by-gene pathway mapping substitutes temporarily the KO mapping which is the right way to go (KO mapping is not working right now on KEGG)
								//If it persists not working (i sent a ticket) we may consider put 'proxy' gene expressions averaging the KO expression expected (which is probably not a good solution because sometimes a gene is in several KOs).
								/*c.computeExpressionProfile(this.expData);
								for(String g:c.term.geneExs.keySet())									{
									ArrayList<Float> ex=c.term.geneExs.get(g);
									if(ex!=null && ex.size()==this.expData.getNumConditions())
										{
										Color color=getColor(ex, this.selectedCol);
										cad=cad+g+"%09%23"+Color2Hex(color.getRed(), color.getGreen(), color.getBlue())+",%23446644/";
										}
									}*/
								
								if(c.term.koExs==null || c.term.koExs.size()==0)	
									c.computeKOIDExpression(this.expData);
								for(String ko:c.term.koExs.keySet())//generic KO pathway
									{
									ArrayList<Float> ex=c.term.koExs.get(ko);
									if(ex!=null && ex.size()==this.expData.getNumConditions())
										{
										Color color=getColor(ex, this.selectedCol);
										cad=cad+ko+"%09%23"+Color2Hex(color.getRed(), color.getGreen(), color.getBlue())+",%23007700/";
										}
									}
								
								//NOTE: New KEGG API version has made the coloring privative (shamefully)
								// We can try by modifying the html received with HTML5 area tag and style attribute,
								//<area shape="rect" coords="656,276,702,293" href="/dbget-bin/www_bget?K16030" title="K16030 (asmB)">  <-- Here, include the style attribute for each term... not sure if it will be possible.
								String url="https://www.kegg.jp/kegg-bin/show_pathway?"+this.expData.organismKegg.replace("kke", "map")+c.term.id+"/"+cad.toUpperCase();
								System.out.println(url);
								java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
								}
							else	//just browse to the pathway, no map
								{
								String url="http://www.kegg.jp/kegg-bin/show_pathway?map"+c.term.id;
								System.out.println(url);
								java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
								}
						}catch(Exception e){e.printStackTrace();}
						}
					}//for each hovered cell (usually one)
				redraw();
				
				//Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
				//setCursor(normalCursor);
				this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
				
				}//if KEGG
			else if(type==GO || type==BP ||type==CC || type==MF)
				{
				for(Cell c:hoveredCells)
					{
					if(c.level==this.levelThreshold)
						{
						try{
							String url="http://amigo.geneontology.org/cgi-bin/amigo/term_details?term="+c.term.id;
							System.out.println("Browsing to "+url);
							java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
							}catch(Exception e){e.printStackTrace();}
						}
	
					}
				}
			else if(type==REACTOME)
				{
				
				for(Cell c:hoveredCells)
					{
					if(c.level==minHoveredLevel)
						{
						try{
							String url="http://www.reactome.org/cgi-bin/eventbrowser?DB=gk_current&ID="+c.term.id;
							System.out.println("Browsing to "+url);
							java.awt.Desktop.getDesktop().browse(java.net.URI.create(url));
							}catch(Exception e){e.printStackTrace();}
						}
					}
				}
			
			//altDown=false;
			return;
			}
		else	//set selection
			{
//			if(isFocusOwner())
			if(focused)
				{
				if(selectedCell!=minHoveredCell)
					selectedCell=minHoveredCell;
				else
					selectedCell=null;
				redraw();
				}
			}
		}
		if(mouseEvent.getClickCount()==2 && hoveredCells.size()>0)						//draw heatmap
			{
			if(expData==null)	return;
			if(minHoveredCell!=null)
				{
				int x=voronto.getLocation().x+this.getGraphics().width;
				int y=voronto.getLocation().y;
				if(minHoveredCell.term.geneExs==null || minHoveredCell.term.geneExs.values().iterator().next().size()!=expData.getNumConditions())
					minHoveredCell.computeExpressionProfile(expData);
				
				
				if(gh!=null)
					((JFrame)((SmoothCanvas)gh.getSurface().getNative()).getFrame()).dispose();;
				gh=new CellHeatmap(minHoveredCell, font, this);
				gh.settings();
				
				PApplet.runSketch(new String[]{"--display=1", "--location="+x+","+y, "--sketch-path=" + gh.sketchPath(), minHoveredCell.term.name}, gh);
			    
				gh.getSurface().setTitle(minHoveredCell.term.name);
				
				//((JFrame)((SmoothCanvas)gh.getSurface().getNative()).getFrame()).setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
				WindowListener wl=((JFrame)((SmoothCanvas)gh.getSurface().getNative()).getFrame()).getWindowListeners()[0];
				((JFrame)((SmoothCanvas)gh.getSurface().getNative()).getFrame()).removeWindowListener(wl);
				//((JFrame)((SmoothCanvas)gh.getSurface().getNative()).getFrame()).setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
				//((JFrame)((SmoothCanvas)gh.getSurface().getNative()).getFrame()).addWindowListener(new SurfaceWindowListener());
				
				}
			
		return;
		}//if double click
    }
    
/**
 * Converts an RGB color to its corresponding hexadecimal value
 * @param r
 * @param g
 * @param b
 * @return
 */
public String Color2Hex(int r, int g, int b)
	{
      Color c = new Color(r,g,b);
      String s=Integer.toHexString( c.getRGB() & 0x00ffffff ).toUpperCase() ;
      while(s.length()<6)
    	  s="0"+s;
      return s;
   }

/**
 * Returns the map with term as the root
 * @param term
 * @return
 */
private TreeMap<OntologyTerm, TreeMap> getMap(TreeMap<OntologyTerm, TreeMap> m, OntologyTerm term)
	{
	if(m.containsKey(term))
		return m.get(term);
	else
		{
		for(OntologyTerm ot:m.keySet())
			{
			TreeMap<OntologyTerm, TreeMap> mret=getMap(m.get(ot), term);
			if(mret!=null)	return mret;
			}
		return null;
		}
	}

/**
 * Returns the parent in the hierarchy for a given term
 * @param term
 * @param root root of the current term
 * @return
 */
private OntologyTerm searchParent(OntologyTerm term, TreeMap<OntologyTerm, TreeMap> m, OntologyTerm root)
	{
	if(m==null)	return null;
	if(m.containsKey(term))
		return root;
	else
		{
		for(OntologyTerm ot:m.keySet())
			{
			OntologyTerm parent=searchParent(term, m.get(ot), ot);
			if(parent!=null)	return parent;
			}
		return null;
		}
	}

public void keyPressed()
	{
	switch(keyCode)
		{
		case 17:	//ctrl
			ctrlDown=true;
			break;
		case 18:	//alt
			altDown=true;
			break;
			
		}
	}


/**
 * Export gene list as an expression submatrix, with synonym ids if available
 * @param genes
 * @param saveFile
 */
public void export(HashSet<String> genes, String saveFile)
	{
	try{
		BufferedWriter bw=new BufferedWriter(new FileWriter(saveFile));
		if(expData.ensembl_gene_idHash!=null)
			bw.write("entrezgene\tensembl_gene_id\texternal_gene_id");
		else
			bw.write(expData.chip);
		
		for(String c:expData.conditionNames)
			bw.write("\t"+c);
		
		bw.newLine();
		for(String g:genes)
			{
			if(expData.ensembl_gene_idHash!=null)
				{
				String entrezgene=expData.getSynonym(g.toLowerCase(), ExpressionData.ENTREZ);
				String ensembl_gene_id=expData.getSynonym(g.toLowerCase(), ExpressionData.ENSEMBL);
				String external_gene_id=expData.getSynonym(g.toLowerCase(), ExpressionData.SYMBOL);
				
				if(expData.chip.equals("entrezgene") && entrezgene==null)	entrezgene=g;
				if(expData.chip.equals("ensembl_gene_id") && ensembl_gene_id==null)	ensembl_gene_id=g;
				if(expData.chip.equals("external_gene_id") && external_gene_id==null)	external_gene_id=g;
				
				if(entrezgene==null)	entrezgene="";
				if(ensembl_gene_id==null)	ensembl_gene_id="";
				if(external_gene_id==null)	external_gene_id="";
				
				bw.write(entrezgene+"\t"+ensembl_gene_id+"\t"+external_gene_id);
				}
			else
				bw.write(g);
			for(Double ex:expData.getExpressionProfile(expData.getGeneId(g)))
				bw.write("\t"+ex);
			bw.newLine();
			}
		bw.close();
	}catch(Exception e){e.printStackTrace();}
	}

/**
 * Export term list as an expression matrix
 * @param terms
 * @param saveFile
 */
public void exportCells(List<Cell> cells, String saveFile, boolean recursive)
	{
	try{
		List<String> names=new ArrayList<String>();
		BufferedWriter bw=new BufferedWriter(new FileWriter(saveFile));
		bw.write("name");
		
		for(String c:expData.conditionNames)
			bw.write("\t"+c);
		
		bw.newLine();
		int cont=0;
		for(Cell c:cells)
			{
			//TODO: temporal!
			if(!names.contains(c.term.id) && c.term.degs.size()/c.term.geneIds.size()>0.1 && c.term.pvalueBH==-1)
			//if(!names.contains(c.term.id))
					{
				bw.write(c.term.name);
				for(float ex:c.expressionLevel)
					{
					bw.write("\t"+ex);
					cont++;
					}
				bw.newLine();
				names.add(c.term.id);
				}
			if(recursive)
				recursiveExportCells(bw, names, c);
			}
		bw.close();
		System.out.println("Of which unique are "+cont+" occurrences");
		
	}catch(Exception e){e.printStackTrace();}
	}

public void recursiveExportCells(BufferedWriter bw, List<String> added, Cell c)
	{
	try{
		if(c.term.degs.size()>0)
			System.out.println(c.term.name+"\t"+c.term.degs.size()+"\t"+(float)(c.term.degs.size())/c.term.geneExs.size()+"\t"+c.term.pvalueBH);
		if(!added.contains(c.term.id) && (float)(c.term.degs.size())/c.term.geneIds.size()>0.1 && (c.term.pvalueBH==-1 || c.term.pvalueBH>0.1))
		//if(!added.contains(c.term.id))
			{
			System.out.println("THIS ONE!");
			bw.write(c.term.name);
			for(float ex:c.expressionLevel)
				bw.write("\t"+ex);
			bw.newLine();
			added.add(c.term.id);
			}
		if(c.subcells!=null)
			for(Cell sc:c.subcells)
				{
				if(sc==null || bw==null || added==null)
						System.out.println("null!");
				recursiveExportCells(bw, added, sc);
				}
	}catch(Exception e){e.printStackTrace();}
	}
/**
 * Export the genes in the term as an expression submatrix, with synonym ids if available
 * @param term
 */
public void export(OntologyTerm term)
	{
	JFileChooser selecFile = new JFileChooser();
	selecFile.addChoosableFileFilter(new TextFileFilter());
	
	if(expData!=null)	
		{
		selecFile.setCurrentDirectory(new File(expData.filePath));
		selecFile.setSelectedFile(new File(term.name.trim()+".txt"));
		}
	else
		selecFile.setSelectedFile(new File(term.name.trim()+".txt"));
	
//	int returnval = selecFile.showSaveDialog(this);
	int returnval = selecFile.showSaveDialog(this.voronto);

	try{
		if(returnval==JFileChooser.APPROVE_OPTION)
			{
			BufferedWriter bw=new BufferedWriter(new FileWriter(selecFile.getSelectedFile()));
			if(expData.ensembl_gene_idHash==null)
				bw.write(expData.chip);
			else
				bw.write("entrezgene\tensembl_gene_id\texternal_gene_id");
			for(String c:expData.conditionNames)
				bw.write("\t"+c);
			
			bw.newLine();
			for(String id:term.geneExs.keySet())
				{
				if(expData.ensembl_gene_idHash!=null)
					{
					String entrezgene=expData.getSynonym(id.toLowerCase(), ExpressionData.ENTREZ);
					String ensembl_gene_id=expData.getSynonym(id.toLowerCase(), ExpressionData.ENSEMBL);
					String external_gene_id=expData.getSynonym(id.toLowerCase(), ExpressionData.SYMBOL);
					
					bw.write(entrezgene+"\t"+ensembl_gene_id+"\t"+external_gene_id);
					}
				else
					bw.write(id);
				for(Float ex:term.geneExs.get(id))
					bw.write("\t"+ex);
				bw.newLine();
				}
			bw.close();
			}
	}catch(Exception e){e.printStackTrace();}
	}

public void keyReleased()
	{
	Cursor hourglassCursor, normalCursor;
	
	//System.out.println("key: "+keyCode);
	//TODO: check if keycodes for arrows are the same on any platform (should be)
	switch(keyCode)
		{
		case 39://right
			if(expData!=null)
				{
				selectedCol=Math.min(expData.getNumConditions()-1, selectedCol+1);
				expression2color(expData, selectedCol);
				}
			break;
		case 37://left
			selectedCol=Math.max(0, selectedCol-1);
			if(expData!=null)	expression2color(expData, selectedCol);
			break;
		case 38://up
			levelThreshold=Math.max(1, levelThreshold-1);
			selectedCell=null;
			mouseMoved();
			setOntologyName();
			break;
		case 40://down
			levelThreshold=Math.min(v.maxLevel+1, levelThreshold+1);
			selectedCell=null;
			mouseMoved();
			setOntologyName();
			break;
		case 10://enter
			//System.gc();
			searchedCells.clear();
			
			//dig deep into hierarchy
			//0) get hovered cell
			Cell cell=null;
			if(minHoveredCell!=null)	cell=minHoveredCell;
			
			//1) perform Benrhardt tessellation under it (to max 2 depth level again?)
			if(cell==null)	{System.err.println("No hovered cell"); return;}
			if(getMap(map, cell.term).size()<=1)	
				{
				System.err.println("Leaf cell, no further voronoi map");
				JOptionPane.showMessageDialog(null, "This term has zero or one subterms, tessellation will not be computed", "Warning", JOptionPane.WARNING_MESSAGE);
				return;
				}
			
			//hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
			//setCursor(hourglassCursor);
			this.getSurface().setCursor(Cursor.WAIT_CURSOR);

			
			try{tessellate(cell.term);}catch(Exception e){e.printStackTrace();}
			
			//renew search if any
			if(degsi!=null && degsi.size()>0)
				{
				if(enrichedTerms.size()>0)
					{
					for(Cell c:this.cells)
						searchedCells.addAll(recursiveEnrichmentRecovery(c,enrichmentThreshold,0));
					}
				else{
					if(records!=null)
						this.deEnrichment(records, entrezConv);
					else
						this.giSearch(degsi, termThreshold, null);
					}
				}
			else if(searchedText!=null && searchedText.length()>0)
				this.search(searchedText);
		
			
			this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);

			break;
			
		case 8://supr
			//return back into hierarchy
			this.getSurface().setCursor(Cursor.WAIT_CURSOR);

			if(root.id.equals("root") && root.name.equals("root"))	return;
			searchedCells.clear();
			OntologyTerm parent=new OntologyTerm("root","root");
			if(!altDown)	parent=searchParent(root, map, new OntologyTerm("root", "root"));
			v=tessellations.get(parent);
			if(v==null)	try{tessellate(parent);}catch(Exception e){e.printStackTrace();}
			
			root=parent;
		 	cells=v.getCells();
			hoveredCells.clear();
			selectedCell=null;
			computedLabels=false;
			setOntologyName();
			
			//renew search if any
			if(degsi!=null && degsi.size()>0)
				{
				if(records!=null)
					this.deEnrichment(records, entrezConv);
				else
					this.giSearch(degsi, termThreshold, null);
				}
			else if(searchedText!=null && searchedText.length()>0)
				this.search(searchedText);		
			this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);

			break;
		case 17:	//ctrl
			ctrlDown=false;
			break;
		case 18:	//alt
			altDown=false;
			break;
			
		default:
			//System.out.println(keyCode);
			break;
		}
	switch(key)
		{
		case KeyEvent.VK_ESCAPE:
			break;
		case 'n'://change gene ids shown on cell heatmap
			nameType=(nameType+1)%3;
			if(gh!=null)	gh.redraw();
			break;
		
		case 'l'://cell labels
			SHOW_LABELS=!SHOW_LABELS;
			break;
		case 'r': //pRofile type
			profileType=(profileType+1)%3;
			break;
		case 'e':
			if(this.selectedCell!=null)	export(this.selectedCell.term);
			else if(this.minHoveredCell!=null)	export(this.minHoveredCell.term);
			break;
		case 'c':
			Color[] palette= {Color.BLUE, Color.WHITE, Color.RED};
			JTextField[] samples= {new JTextField(10),new JTextField(10),new JTextField(10)};
			String[] labels= {"lowest expression", "median expression", "highest expression"};
			//JPanel jp=ColorChooser.getPanelPaleta(palette, labels, samples);
			if(colorFrame==null)
				{
				colorFrame=new ColorFrame(voronto, palette, labels, samples);
				colorFrame.setUndecorated(true);
				colorFrame.requestFocusInWindow();
				}
			else if(colorFrame.isVisible())	{colorFrame.setVisible(false); return;}
			
			colorFrame.setLocation((int)(voronto.getLocation().x+this.width*0.5-colorFrame.getWidth()*0.5), (int)(voronto.getLocation().y+this.height*0.5-colorFrame.getHeight()*0.5));
			colorFrame.setVisible(true);
			break;
		case 'm':
			this.exportCells(Arrays.asList(cells), "/Users/rodri/Documents/investigacion/sybaris/Voronto/plos/test.txt", true);
			break;
		case 'p':
			printImage();
			
			break;
		case 'f':	//searcher
			if(searchFrame==null)
				searchFrame=new SearchFrame(voronto);
			else if(searchFrame.isVisible())	{searchFrame.setVisible(false); return;}
			searchFrame.textField.setText("");
			searchFrame.setLocation((int)(voronto.getLocation().x+this.width*0.5-searchFrame.getWidth()*0.5), (int)(voronto.getLocation().y+this.height*0.5-searchFrame.getHeight()*0.5));
			searchFrame.setVisible(true);
			searchFrame.setAlwaysOnTop(true);
			drawSearch=false;
			break;
		case 'd':	//difexp searcher
			if(deSearchFrame==null)
				{
				deSearchFrame=new DifExpSearch(this);
				deSearchFrame.setUndecorated(true);
				deSearchFrame.requestFocusInWindow();
				}
			else if(deSearchFrame.isVisible())	{deSearchFrame.setVisible(false); return;}
			
			deSearchFrame.setLocation((int)(voronto.getLocation().x+this.width*0.5-deSearchFrame.getWidth()*0.5), (int)(voronto.getLocation().y+this.height*0.5-deSearchFrame.getHeight()*0.5));
			deSearchFrame.setVisible(true);
			break;
		case 'h':	//helper
			try{
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://vis.usal.es/~visusal/voronto/voronto/Help.html"));
			}catch(Exception e){JOptionPane.showMessageDialog(null, "Error", "Cannot open browser. Please visit http://vis.usal.es/~visusal/voronto/voronto/Help.html for help", JOptionPane.ERROR_MESSAGE);}
			break;
		case 'b': //go back to input interface
			if(voronto!=null)
				{
				gh.getSurface().setVisible(false);
				voronto.goBack();				
				}
			else
				System.err.println("Voronto is null, no chance to go back");
			break;
		case 's'://change medium value to mean or median or min
			whiteValue=(whiteValue+1)%3;	
			//whiteValue=(whiteValue+1)%2;	//JUST FOR THE VIDEO
			
			expression2color(this.expData);
			redraw();
			if(gh!=null)	gh.redraw();
			break;
		}
	redraw();
	}

private ArrayList<Cell> recursiveEnrichmentRecovery(Cell cell, double threshold, int level) {
	ArrayList<Cell> retList=new ArrayList<Cell>();
	if(cell.subcells==null || cell.subcells.length==0) //leaf node
		{
		cell.searchedColor=null;
		if(enrichedTerms.get(cell.term.id)!=null && enrichedTerms.get(cell.term.id).significant)
			{
			cell.term=enrichedTerms.get(cell.term.id);
			System.out.println("Adding\t"+cell.term.name);
			cell.searchedColor=new Color(0,200,0);
			retList.add(cell);
			}
		}
	else		//non leaf-node
		{
		boolean contains=false;
		if(enrichedTerms.get(cell.term.id)!=null && enrichedTerms.get(cell.term.id).significant)
			{
			contains=true;
			cell.term=enrichedTerms.get(cell.term.id);
			}
	
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveEnrichmentRecovery(c,threshold, level+1);
			retList.addAll(sublist);
			if(sublist.size()>0)	add=true;
			}
		
		//count subterms
		HashSet<String> set=new HashSet<String>();
		for(Cell c:retList)	
			if(c.term.significant)	
				set.add(c.term.id);
		cell.term.numRelevantSubterms=set.size();
		if(add || contains)
			{
			if(contains)
				cell.searchedColor=new Color(0,200,0);
			else
				{
				cell.searchedColor=new Color(0,100,0);
				cell.term=enrichedTerms.get(cell.term.id);
				}
			retList.add(cell); 
			}
		}
	return retList;
}

public void printImage()
	{
	JFileChooser selecFile = new JFileChooser();
	TIFFileFilter png=new TIFFileFilter();
	selecFile.addChoosableFileFilter(png);
	SVGFileFilter svg=new SVGFileFilter();
	selecFile.addChoosableFileFilter(svg);
	//selecFile.addChoosableFileFilter(new PDFFileFilter());
	
	if(imagePath!=null)
		selecFile.setCurrentDirectory(new File(this.imagePath));
	else if(expData!=null)
		selecFile.setCurrentDirectory(new File(expData.filePath));
	
	if(expData!=null)	
		selecFile.setSelectedFile(new File(expData.organismKegg+"-"+expData.conditionNames[this.selectedCol].replace("/", "-")));
	else
		selecFile.setSelectedFile(new File("voronto"));
	
	selecFile.setAcceptAllFileFilterUsed(false);

	int returnval = selecFile.showSaveDialog(this.voronto);


	if(returnval==JFileChooser.APPROVE_OPTION)
		{
		saving=true;
		try{Thread.sleep(100);}catch(Exception e){}
		if(selecFile.getFileFilter()==png)	//for TIFF
			{
			int scaleFactor=3;
			PGraphics hires = createGraphics(width*scaleFactor, height*scaleFactor, JAVA2D);
			PGraphics ant=this.g;
			this.g=hires;
			beginRecord(hires);
			hires.scale(scaleFactor);
			draw();
			endRecord();
			hires.save(selecFile.getSelectedFile().getAbsolutePath());
			this.g=ant;
			}
		else if(selecFile.getFileFilter()==svg)	//for SVG
			{
			font=loadFont("es/usal/voronto/font/AppleSymbols-14.vlw");
			computedLabels=false;
			
			P8gGraphicsSVG hires = (P8gGraphicsSVG)createGraphics(width, height, P8gGraphicsSVG.SVG);
			PGraphics ant=this.g;
			this.g=hires;
			beginRecord(hires);
			draw();
			hires.endRecord(selecFile.getSelectedFile().getAbsolutePath());
			
			font=loadFont("es/usal/voronto/font/AppleSymbols-14.vlw");
			//font=loadFont("es/usal/voronto/font/Arial-14.vlw");
			computedLabels=false;
			this.g=ant;
			}
		
		try{
		imagePath=selecFile.getSelectedFile().getParent();
		BufferedWriter bw=new BufferedWriter(new FileWriter(voronto.imgPath));
		bw.write(imagePath);
		bw.close();
		}catch(IOException e){e.printStackTrace();}
		
		saving=false;
		redraw();
		}
	}

public void exit()
	{
	//System.out.println("Exiting...");
	return;//i don't want esc to close the program
	}

/**
 * Computes a new tesselleation with cell as root node
 * @param cell
 */
public void tessellate(OntologyTerm term) throws Exception
	{
	boolean translate=false;
	if(tessellations.get(term)==null)	//search for the tessellation for this term as root
		{
		TreeMap<OntologyTerm, TreeMap> newMap=getMap(map,term);
		int numTerms=this.voronto.getTermList(newMap).size();
		if(numTerms<10000)
			maxDepth=20;	//GO ontologies now only curated, no IAE (much smaller)
		else 
			maxDepth=3;
	
		try
			{v=new BernhardtTessellation(getMap(map, term), expData, width, height-100, type, this.maxDepth);}
		catch(Exception e)
			{
			if(e.getMessage().startsWith("No mapped genes"))
				JOptionPane.showMessageDialog(null, "This term has "+getMap(map, term).size()+" subterms, but with no annotated genes, tessellation will not be computed.", "Warning", JOptionPane.WARNING_MESSAGE);
			else
				JOptionPane.showMessageDialog(null, "An error occurred during tessellation.", "Warning", JOptionPane.WARNING_MESSAGE);
			return;
			}
		tessellations.put(term, v);
	 	translate=true;
		}
	else
		v=tessellations.get(term);
	
	//2) substitute the visualization level (the root node)
	root=term;
	cells=v.getCells();
	
	if(type==VoronoiVisualization.KEGG)	this.levelThreshold=v.maxLevel+1;
 	if(type==VoronoiVisualization.SLIMBP || type==VoronoiVisualization.SLIMCC)	this.levelThreshold=v.maxLevel;
	
	if(translate)	for(Cell c:cells)	recursiveRegionTranslation(c, 0, START_Y);
	computedLabels=false;
	if(expData!=null)	
		{
		mapExpression(expData);
		expression2color(expData);
		}
	
	if(expData!=null)	
		{
		expression2color(expData);
		for(int i= 0;i<expData.getNumConditions();i++)
			expression2color(expData,i);
		System.out.println("Finished mapping colors");
		}
	
	setOntologyName();
	hoveredCells.clear();
	selectedCell=null;
	return;
	}

public void mouseExited(MouseEvent e)
	{
	minHoveredCell=null;
	minHoveredLevel=-1;
	hoveredCells.clear();
	redraw();
	}

public void mouseMoved() {
    if(v==null)	return;
	
	Cell[] cells = v.getCells();
	
	//1) Check hovered cells
	//long t0=System.currentTimeMillis();
	synchronized(hoveredCells)
		{
		minHoveredCell=null;
		minHoveredLevel=-1;
		hoveredCells.clear();
		
		for(int i=0; i<cells.length; i++)
			{
			hoveredCells.addAll(recursiveSearchHovered(cells[i]));
			if(hoveredCells.size()>0)	break;
			}
		}
	
	//2) If no hovered cells are selected, check if other element is selected to show help messages
	hoveredBox=-1;
	if(hoveredCells.size()==0)
		{
		if(organismBox!=null && organismBox.contains(mouseX, mouseY))
			hoveredBox=ORGANISM;
		else if(conditionBox!=null && conditionBox.contains(mouseX, mouseY))
			hoveredBox=CONDITION;
		else if(ontologyBox!=null && ontologyBox.contains(mouseX, mouseY))
			hoveredBox=ONTOLOGY;
		else if(cellBox!=null && cellBox.contains(mouseX, mouseY))
			hoveredBox=CELL;
		else if(scaleBox!=null && scaleBox.contains(mouseX, mouseY))
			hoveredBox=SCALE;
		}
	//System.out.println("Time in search: "+(System.currentTimeMillis()-t0)/1000.0);
	
    // update the screen (run draw once)
	redraw();
	}


/**
 * Returns a list of cells hovered by the mouse, that is, the cell hovered an all its parents.
 * @param cell
 * @return
 */
public ArrayList<Cell> recursiveSearchHovered(Cell cell)
	{
	ArrayList<Cell> ret=new ArrayList<Cell>();

	if(cell.region!=null)
		{
		//1) Check if region is hovered
		if(cell.region.getPolygon().contains(mouseX, mouseY))
			{
			if(minHoveredCell==null)	minHoveredCell=cell;
			if(minHoveredCell.level<cell.level && levelThreshold>=cell.level)	minHoveredCell=cell;
			if(cell.level>minHoveredLevel)	minHoveredLevel=cell.level;
			
			ret.add(cell);
			
			//2) Continue with recursion (only if the parent is hovered)
			if(cell.subcells!=null && cell.subcells.length>0)	
				for(Cell cc:cell.subcells)
					{
					ArrayList<Cell> toadd=recursiveSearchHovered(cc);
					if(toadd!=null && toadd.size()>0)
						{
						ret.addAll(toadd);
						break;
						}
					}
			}
		}
	return ret;
	}

/**Maps expression in md (all columns) to voronoi regions
 * 
 * @param md
 */
public void expression2color(ExpressionData md)
	{
	if(md==null)
		{
		for(Cell c:cells)	
			recursiveExpressionNormalization(c, -1);
		}
	else
		{	
		for(int i= 0;i<expData.getNumConditions();i++)
			expression2color(md,i);
		}
	return;
	}

public void mapExpression(ExpressionData md)
	{
	long t=System.currentTimeMillis();
	maxSd=new float[md.getNumConditions()];
	Cell[] cells=v.getCells();
	for(Cell c:cells)
		{
		//long t1=System.currentTimeMillis();
		//System.out.println("recursive mapping for "+c.term.name+" with "+c.term.geneIds.size()+" annot genes");
		recursiveExpressionMapping(c, md);
		//System.out.println("it took "+(System.currentTimeMillis()-t1)/1000.0);
		}
	System.out.println("Mapping expression to ontology takes "+(System.currentTimeMillis()-t)/1000.0);
	return;
	}

public void setScale(ExpressionData md, int col)
	{
	//Different options: 
		//a) against min and max of overall md
		//b) against min and max of overall md[col] -> we will start by this one, but not sure what's best, possibly last one
		//c) against min and max of mapped cell expressions
		switch(SCALE_MODE)
			{
			case SCALE_MATRIX:
			//	System.out.println("Changing scale to mode MATRIX");
				minExp=md.min;
				maxExp=md.max;
				avgExp=md.average;
				medianExp=md.median;
				break;
			case SCALE_CONDITION:
			//	System.out.println("Changing scale to mode CONDITION");
				minExp=md.minCols[col];
				maxExp=md.maxCols[col];
				avgExp=md.averageCols[col];
				medianExp=md.medianCols[col];
				break;
			case SCALE_ONTOLOGY:
				//System.out.println("Changing scale to mode ONTOLOGY");
				minExp=1000000000;
				maxExp=-1000000000;
				avgExp=0;
				contAvgExp=0;
				break;
			}
		//System.out.println("Mapping expression for "+md.getColumnLabel(col)+" in scale mode "+SCALE_MODE);
		
		
		if(SCALE_MODE==SCALE_ONTOLOGY)
			avgExp/=contAvgExp;
		}

/**Maps expression in md (in experiment col) to voronoi regions
 * 
 * @param md
 */
public void expression2color(ExpressionData md, int col)
	{
	
	if(expData==null)	{System.err.println("No expression data provided for mapping"); return;}
	setScale(md, col);
	Cell[] cells=v.getCells();
	for(Cell c:cells)	
		recursiveExpressionNormalization(c, col);
	return;
	}

public void recursiveExpressionMapping(Cell cell, ExpressionData md)
	{
	if(cell.expressionLevel.size()==md.getNumConditions())	
		return; //already computed (in GO, we can find inside loops!)

	cell.computeExpression(md, type);
	if(cell.subcells!=null && cell.subcells.length>0)	//for internal nodes
		for(Cell cc:cell.subcells)
			recursiveExpressionMapping(cc, md);
	return;
	}

public void recursiveExpressionNormalization(Cell cell, int column)
	{
	//if(cell.term.name.startsWith("Aminoacyl"))
	//	System.out.println("Aminoacyl");
	//1) Normalize color
	if(column>=0 && !Float.isNaN(cell.expressionLevel.get(column)) )
		{
		int h=-1;
		switch(whiteValue)
			{
			case MEAN:	//raw, two color scale
				switch(COLOR_MODE)
					{
				
					case COLOR_EXPRESSION:
						if(cell.expressionLevel.get(column)>avgExp)
							{
							h=(int)Math.round(255-((cell.expressionLevel.get(column)-avgExp)/(maxExp-avgExp)*255));
							//cell.color.set(column,new Color(255, h, h));
							cell.color.set(column, palette[Math.abs(h-255)+255]);
							}
						else
							{
							h=(int)Math.round(255-(Math.abs(cell.expressionLevel.get(column)-avgExp)/Math.abs(minExp-avgExp))*255);
							//cell.color.set(column, new Color(h,h, 255));
							cell.color.set(column, palette[h]);
							}
						
						
						//cell.color.set(column, palette[(int)Math.round((cell.expressionLevel.get(column)-avgExp)/(maxExp-avgExp)*510)]);
						
						break;
						
				
					case COLOR_DEVIATION:
						h=(int)Math.round(255-(Math.abs(cell.expressionLevel.get(column)-avgExp)/Math.max(Math.abs(avgExp-minExp), Math.abs(avgExp-maxExp)))*255);
						cell.color.set(column, new Color(h,255,h));
						break;
					case COLOR_INTERNAL_DEVIATION://TODO: testing
						System.out.println("Getting it for "+cell.term.name+" with length "+cell.expressionDeviation.size());
						if(cell.expressionDeviation.size()>0)	
							h=(int)(255-(cell.expressionDeviation.get(column)/maxSd[column])*255);
						else	h=255;
						cell.color.set(column, new Color(h,255,h));
						break;
					}
				break;
			case MEDIAN:	//quantile coloring --> note: color scaling 
				switch(COLOR_MODE)
					{
					case COLOR_EXPRESSION:
						int q=-1;
						
						if(SCALE_MODE==SCALE_CONDITION)		q=expData.getQuantile(cell.expressionLevel.get(column), column);
						else if(SCALE_MODE==SCALE_MATRIX)	q=expData.getQuantile(cell.expressionLevel.get(column));
						if(q>=50)
							{
							h=(int)Math.round(255-((q-50.0)/50)*255);
							//cell.color.set(column, new Color(255, h, h));
							cell.color.set(column, palette[Math.abs(h-255)+254]);
							}
						else
							{
							h=(int)Math.round(255-((50.0-q)/50)*255);
							//cell.color.set(column, new Color(h,h, 255));
							cell.color.set(column, palette[h]);
							}
						
						break;
					case COLOR_DEVIATION:
						System.err.println("Option not supported for quantiles");
						break;
					case COLOR_INTERNAL_DEVIATION://TODO: testing
						System.err.println("Option not supported for quantiles");
					}
				break;
			case MIN: 	//raw, one color scale
				switch(COLOR_MODE)
				{
				case COLOR_EXPRESSION:
					h=(int)Math.round(255-((cell.expressionLevel.get(column)-minExp)/(maxExp-minExp)*255));
					cell.color.set(column,new Color(255, h, h));
					break;
				case COLOR_DEVIATION:
					h=(int)Math.round(255-(Math.abs(cell.expressionLevel.get(column)-minExp)/(maxExp-minExp)*255));
					cell.color.set(column, new Color(255,h,h));
					break;
				case COLOR_INTERNAL_DEVIATION://TODO: testing
					System.out.println("Getting it for "+cell.term.name+" with length "+cell.expressionDeviation.size());
					if(cell.expressionDeviation.size()>0)	
						h=(int)(255-(cell.expressionDeviation.get(column)/maxSd[column])*255);
					else	h=255;
					cell.color.set(column, new Color(h,255,h));
					break;
				}
				break;
			default:
				if(column<0)	
					{
					cell.color=new ArrayList<Color>();
					cell.color.add(0, new Color(255,255,255));
					}
				else			
					cell.color.set(column, new Color(255,255,255));
				break;
		}
	}
	//2) Continue with recursion
	//if(cell.subcells!=null && cell.subcells.length>1)	
	if(cell.subcells!=null && cell.subcells.length>0)	
		for(Cell cc:cell.subcells)
			recursiveExpressionNormalization(cc, column);
	}


public void changeColors(Color[] newColors)
	{
	Color colorL=newColors[0];
	Color colorM=newColors[1];
	Color colorH=newColors[2];
	int grain=255;
	//int grain=123;
	palette = new Color[grain*2];
	
	float stepR=Math.abs(colorL.getRed()-colorM.getRed())/(float)grain;
	float stepG=Math.abs(colorL.getGreen()-colorM.getGreen())/(float)grain;
	float stepB=Math.abs(colorL.getBlue()-colorM.getBlue())/(float)grain;
	int senseR=(colorL.getRed()>colorM.getRed())?-1:1;
	int senseG=(colorL.getGreen()>colorM.getGreen())?-1:1;
	int senseB=(colorL.getBlue()>colorM.getBlue())?-1:1;
	
	
	for(int i=0;i<grain;i++)
		{
		palette[i]=new Color((int)(colorL.getRed()+senseR*stepR*i), 
				(int)(colorL.getGreen()+senseG*stepG*i),
				(int)(colorL.getBlue()+senseB*stepB*i));
		}
	
	stepR=Math.abs(colorM.getRed()-colorH.getRed())/grain;
	stepG=Math.abs(colorM.getGreen()-colorH.getGreen())/grain;
	stepB=Math.abs(colorM.getBlue()-colorH.getBlue())/grain;
	senseR=(colorM.getRed()>colorH.getRed())?-1:1;
	senseG=(colorM.getGreen()>colorH.getGreen())?-1:1;
	senseB=(colorM.getBlue()>colorH.getBlue())?-1:1;
	
	for(int i=0;i<grain;i++)
		{
		palette[i+grain]=new Color((int)(colorM.getRed()+senseR*stepR*i), 
				(int)(colorM.getGreen()+senseG*stepG*i),
				(int)(colorM.getBlue()+senseB*stepB*i));
		}
	
	if(expData!=null)
		//mapExpression(expData);
		expression2color(expData);
	redraw();
	}


public void search(String searchText) {
	this.getSurface().setCursor(Cursor.WAIT_CURSOR);
	
	searchedText=searchText;
	searchFrame.setVisible(false);
	//searchPApplet.searchText="";
	System.out.println("Searching...");
	searchedCells.clear();
	
	cleanSearch();
	degsi=new HashSet<String>();
	
	if(searchText==null || searchText.length()==0)	
		{
		redraw(); 
		return;
		}
	

	for(Cell c:cells)
		searchedCells.addAll(recursiveSearch(c, searchText, 0));
	System.out.println("Found "+searchedCells.size()+" occurrences");
	redraw();
	if(gh!=null)	gh.redraw();

	//setCursor(normalCursor);
	this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
	return;
	}

private void cleanSearch()
	{
	degsi=null;
	records=null;
	}

/**
 * Searches terms that have more than termThreshold genes differentially expressed on conditions g1 respect to conditions g2 (exp on g1 is above/below threshold from expression on g2)
 * @param g1	First group of conditions
 * @param g2	Second group of conditions
 * @param type	tells if they must be upreg on g1 respect to g2 (0) or the other way round (1)
 * @param threshold	differential threshold between the genes' average expression levels on g1 respect to g2
 * @param termThreshold	minimum number of differential genes that a term must contain to be highlighted
 */
public void deTermSearch(int[] g1, int[] g2, int type, float threshold, int termThreshold, String saveFile) {
	searchedText="";
	//Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
	//setCursor(hourglassCursor);
	//Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
	this.getSurface().setCursor(Cursor.WAIT_CURSOR);
	
	if(deSearchFrame!=null)	deSearchFrame.setVisible(false);
	System.out.println("Expression Search...");
	searchedCells.clear();
	
	if((g1==null || g1.length==0) && (g2==null || g2.length==0))
		{
		JOptionPane.showMessageDialog(null, "Please select at least one condition in order to perform a search", "No selected conditions", JOptionPane.INFORMATION_MESSAGE);
		redraw(); 
		//setCursor(normalCursor); 
		this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
		return;
		}
	
	//Highlight terms with more than termThreshold degs
	try{
		{
		for(Cell c:cells)
			searchedCells.addAll(recursiveDESearch2(c, g1, g2, type, threshold, 0));
		System.out.println("Found "+searchedCells.size()+" occurrences");
		}
	if(saveFile!=null)
		{
		List<Cell> real=new ArrayList<Cell>();	//we take only the terms that fulfill the criteria, not the terms that contain terms that fulfill the de criteria 
		for(Cell c:searchedCells)
			if(c.searchedColor.getGreen()==200)
				real.add(c);
		System.out.println("Of which non-containers are "+real.size()+" occurrences");
		exportCells(real, saveFile, false);
		}
	}catch(Exception e){System.err.println("Error during differentially expressed terms search:" +e.getMessage()); e.printStackTrace();}
	
	redraw();
	if(gh!=null)	gh.redraw();
	
	this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
	//setCursor(normalCursor);
	
	return;
	}

/**
 * 
 * Searches for occurrences in the ontology of a given list of genes of interest (GI)
 * @param degsi
 * @param termThreshold
 * @param saveFile
 */
public void giSearch(HashSet<String> degsi, int termThreshold, String saveFile) {
	searchedText="";
	//Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
	//setCursor(hourglassCursor);
	//Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
	this.getSurface().setCursor(Cursor.WAIT_CURSOR);
	
	if(deSearchFrame!=null)	deSearchFrame.setVisible(false);
	System.out.println("Expression Search...");
	searchedCells.clear();
	
	
	HashSet<String> degsi2=new HashSet<String>();
	for(String g:degsi)		{degsi2.add(g.toLowerCase());}
	this.degsi=degsi2;
	
	//Highlight terms with more than termThreshold degs
	try{
	if(degsi!=null && degsi.size()>0)
		{
		for(Cell c:cells)
			searchedCells.addAll(recursiveDESearch(c, degsi, termThreshold, 0, false));
		System.out.println("Found "+searchedCells.size()+" occurrences");
		}
	
	}catch(Exception e){System.err.println("Error during DAVID query:" +e.getMessage()); e.printStackTrace();}
	if(saveFile!=null)
		this.export(degsi, saveFile);
	
	redraw();
	if(gh!=null)	gh.redraw();
	
	//setCursor(normalCursor);
	this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
	
	return;
	}

/**
 * Performs a Fisher's Exact Test for each ontology term with the genes of interest (degsi).
 * TODO: Then performs a multiple hypothesis contrast correction (BH)
 * @param degsi Genes of interest
 * @param id_type
 * @param database
 * @param saveFile
 * @param pval Threshold for (BH corrected) p-values. Only terms lower than pval will be selected as enriched
 * @param count
 */
public void fisherEnrichment(HashSet<String> degsi, int id_type, String database, String saveFile, float pval, int count)
	{
	//Stats.fisherTest(degsi.size(), termGenes, relevantGenes, relevantGenesInTerm)
	}

public ArrayList<Cell> recursiveFisher(Cell cell, HashSet<String> degs, int threshold, int level)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	cell.term.resetEnrichment();
	if(cell.subcells!=null && cell.subcells.length>0)
		{
//		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveFisher(c, degs, threshold, level+1);
			retList.addAll(sublist);
	//		if(sublist.size()>0)	
		//		add=true;
			}
		}
	
	HashSet<String> intersect=new HashSet<String>();
	intersect.addAll(CollectionUtils.intersection(degs, cell.term.geneExs.keySet()));
	if(intersect.size()>threshold)
		{
		cell.searchedColor=new Color(0,200,0);
		cell.term.degs=intersect;
		retList.add(cell);
		}
	
	
	return retList;
	}

public void deEnrichment(HashSet<String> degsi, int count, float pval, String saveFile) {
	//Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
	//setCursor(hourglassCursor);
	//Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
	this.getSurface().setCursor(Cursor.WAIT_CURSOR);
	
	if(deSearchFrame!=null)	deSearchFrame.setVisible(false);
	System.out.println("Expression Enrichment...");
	searchedCells.clear();
	enrichedTerms.clear();
	this.degsi=degsi;

	
	if(degsi!=null && degsi.size()>2000)
		{
		System.err.println(degsi.size()+" differential genes (maximum is 2000). Narrow down your differential expression criteria");
		JOptionPane.showMessageDialog(null, degsi.size()+" differential genes (maximum is 2000). Narrow down your differential expression criteria", "Too many genes", JOptionPane.INFORMATION_MESSAGE);
		this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
		//setCursor(normalCursor);	
		return;
		}
	
	long t0=System.currentTimeMillis();
	
	//ENRICHMENT BY FISHER'S TEST
	ArrayList<Cell> cc=new ArrayList<Cell>();
	for(Cell c:cells)
		cc.addAll(recursiveDESearch(c, degsi, count, 0, true));
	Set<OntologyTerm> eterms=new LinkedHashSet<OntologyTerm>();
	for(Cell c:cc)
		eterms.add(c.term);
	System.out.println("Enrichment takes "+(System.currentTimeMillis()-t0)+" ms");
	System.out.println(cc.size()+" enriched cells corresponding to "+eterms.size()+" enriched terms");
	//Multiple hypothesis correction:
	System.out.println("Bonferroni: "+Stats.bonferroni(pval, eterms.size()));
	System.out.println("FWER: "+Stats.fwer(pval, eterms));
	System.out.println("FDR: "+Stats.fdr(pval, eterms));
	
//	enrichmentThreshold=Stats.bonferroni(pval, cc.size());
	enrichmentThreshold=Stats.fdr(pval, eterms);
	if(enrichmentThreshold==-1)//No enrichment found
		{	
		System.out.println("No enriched terms found");
		JOptionPane.showMessageDialog(null, " No enriched terms found for the "+this.degsi.size()+" selected genes for a corrected threshold of="+pval+"." , "No enrichment found", JOptionPane.INFORMATION_MESSAGE);
		this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
		//setCursor(normalCursor);
		return;
		}
	
	for(Cell c:cells)
		searchedCells.addAll(recursivePvalSearch(c, enrichmentThreshold, 0));
	for(Cell c:cc)
		enrichedTerms.put(c.term.id, c.term);
	
	saveEnrichment(saveFile, enrichedTerms);
	
	this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
	//setCursor(normalCursor);
	redraw();
	if(gh!=null)	gh.redraw();
	}

public void deEnrichmentDAVID(HashSet<String> degsi, int count, float pval, String saveFile) {
	//Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
	//setCursor(hourglassCursor);
	//Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
	this.getSurface().setCursor(Cursor.WAIT_CURSOR);
	
	if(deSearchFrame!=null)	deSearchFrame.setVisible(false);
	System.out.println("Expression Enrichment...");
	searchedCells.clear();

	
	if(degsi!=null && degsi.size()>2000)
		{
		System.err.println(degsi.size()+" differential genes (maximum is 2000). Narrow down your differential expression criteria");
		JOptionPane.showMessageDialog(null, degsi.size()+" differential genes (maximum is 2000). Narrow down your differential expression criteria", "Too many genes", JOptionPane.INFORMATION_MESSAGE);
		//setCursor(normalCursor);	
		this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
		return;
		}
	
	long t0=System.currentTimeMillis();
	
	// ENRICHMENT BY DAVID
	boolean entrezConv=false;
	if(expData.entrezgeneHash!=null && expData.entrezgeneHash.size()>0)
		{
		entrezConv=true;
		HashSet<String> entrezs=new HashSet<String>();
		for(String s: degsi)	
			{
			String e=expData.entrezgeneHash.get(s.toLowerCase());
			if(e!=null)	entrezs.add(e);
			}
		degsi=entrezs;
		}
	
	System.out.println(degsi.size()+" differential expressed genes found, checking for enrichment...");
	String degs="";
	for(String s:degsi) degs+=s+",";
	degs=degs.substring(0, degs.length()-1);
	
	//2) Determine DAVID database
	String database=null;
	switch(this.type)
		{
		case BP:
			database="GOTERM_BP_FAT";
			break;
		case MF:
			database="GOTERM_MF_FAT";
			break;
		case CC:
			database="GOTERM_CC_FAT";
			break;
		case KEGG:
			database="KEGG_PATHWAY";
			if(expData.organismKegg.equals("kke"))
				{
				JOptionPane.showMessageDialog(null, "DAVID does not support functional enrichment for KEGG's KO ids", "KO ids not supported", JOptionPane.INFORMATION_MESSAGE);
				System.err.println("Cannot make gene set enrichment of KEGG KO ids");
				this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
				//setCursor(normalCursor);	
				return;
				}
			break;
		case REACTOME:
			JOptionPane.showMessageDialog(null, "DAVID does not support functional enrichment for REACTOME", "Ontology not supported", JOptionPane.INFORMATION_MESSAGE);
			System.err.println("Cannot make gene set enrichment of Reactome");
			this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
			//setCursor(normalCursor);	
			return;
		default:
			JOptionPane.showMessageDialog(null, "DAVID does not support functional enrichment of custom ontologies", "Ontology not supported", JOptionPane.INFORMATION_MESSAGE);
			System.err.println("Cannot make gene set enrichment of custom databases");
			this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
			//setCursor(normalCursor);	
			return;
		}
	//3) Determine gene type id
	String id_type=null;
	if(this.expData.chip.equals("entrezgene") || entrezConv)
		id_type="ENTREZ_GENE_ID";
	else if(this.expData.chip.equals("ensembl_gene_id"))
		id_type="ENSEMBL_GENE_ID";
	else if(this.expData.chip.equals("external_gene_id"))
		id_type="GENE_SYMBOL";
	else
		{
		System.err.println("Cannot make gene set enrichment with gene ids of type "+expData.chip);
		JOptionPane.showMessageDialog(null, "Cannot make gene set enrichment with gene ids of type "+expData.chip, "Gene ids not supported", JOptionPane.INFORMATION_MESSAGE);
		this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
		//setCursor(normalCursor);	
		return;
		}
	
	
	try{
	ChartReportClient crc=new ChartReportClient();
	records=crc.invokeService(degs, id_type, database, saveFile, pval, count);
	this.entrezConv=entrezConv;
	this.degsi=degsi;
	deEnrichment(records, entrezConv);
	}catch(Exception e)
		{
		System.err.println("Error during DAVID query:" +e.getMessage()); e.printStackTrace();
		JOptionPane.showMessageDialog(null, "Error during DAVID enrichment: "+e.getMessage(),"Error during enrichment", JOptionPane.INFORMATION_MESSAGE);
		this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
		//setCursor(normalCursor); 	
		return;
		}
	
	this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
	//setCursor(normalCursor);
	redraw();
	if(gh!=null)	gh.redraw();
	}

public void saveEnrichment(String saveFile, HashMap<String, OntologyTerm> terms) {
	if(saveFile==null || saveFile=="")
		return;
	StringBuffer sb;
	sb=new StringBuffer();
	    
	sb.append("Category\tID\tName\tCount\t%\tPValue\tGenes\n");
	DecimalFormat df=new DecimalFormat("#.##");
	
	for (OntologyTerm t : terms.values())
		{
		if(t.significant)
			sb.append(customOntologyName+"\t"
				  +t.id+"\t"
		          +t.name+"\t"
		          +t.degs.size()+"\t"
		          +df.format(t.degs.size()*100.0/t.geneExs.size())+"\t"
		          +t.pvalue+"\t"
		          +t.degs.toString()+"\t"
		          //+simpleChartRecords[j].getListTotals()+"\t"
		          //+simpleChartRecords[j].getPopHits()+"\t"
		          //+simpleChartRecords[j].getPopTotals()+"\t"
		          //+simpleChartRecords[j].getBonferroni()+"\t"
		          //+simpleChartRecords[j].getBenjamini()+"\t"
		          //+simpleChartRecords[j].getAfdr()+"\t"
		          +"\n");
		}
	sb.append("\n");
	try{
		PrintWriter outfile = new PrintWriter(new FileOutputStream(saveFile));
		outfile.println(sb);
		outfile.close();
		}catch(Exception e){e.printStackTrace();}
	}


public void deEnrichment(SimpleChartRecord[] records, boolean entrezConv)
	{
	if(records!=null)
		{
		System.out.println("found "+records.length+" enriched terms");
		for(Cell c:cells)
			if(entrezConv)	searchedCells.addAll(recursiveEnrichedAnnotation(c, records, expData.invertedHash(expData.entrezgeneHash), 0));
			else			searchedCells.addAll(recursiveEnrichedAnnotation(c, records, null, 0));
		System.out.println("Found "+searchedCells.size()+" occurrences");
		}
	
	}

/**
 * Searches the ontology for occurrences of the genes differentially expressed in g1 conditions
 * respect to g2. Differential expression is very raw, no statistics involved. For more elaborate
 * differential expression lists, must be user provided and searched via giSearch
 * @param g1 Conditions of the expression matrix taken as group 1
 * @param g2 Conditions of the expression matrix taken as group 1
 * @param type
 * @param diffT differential threshold
 * @param pval		enrichment threshold 
 * @param saveFile
 */
public void deEnrichment(int[] g1, int[] g2, int type, float diffT, int count, float pval, String saveFile) {
	//Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
	//setCursor(hourglassCursor);
	//Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
	this.getSurface().setCursor(Cursor.WAIT_CURSOR);

	if(deSearchFrame!=null)	deSearchFrame.setVisible(false);
	System.out.println("Expression Enrichment...");
	searchedCells.clear();
	
	if(g1==null || g1.length==0 || g2==null || g2.length==0)	
		{
		JOptionPane.showMessageDialog(null, "No conditions selected, please select some on the list", "No conditions selectd", JOptionPane.INFORMATION_MESSAGE);
		System.err.println("No conditions selected");
		this.getSurface().setCursor(Cursor.DEFAULT_CURSOR);
		//setCursor(normalCursor);
		redraw(); 
		return;
		}
	
	//1) Get genes that are high in g1 respect to g2 (more than threshold)
	degsi=expData.getGeneNames(expData.differentialExpressionGenes(g1, g2, diffT, type));
	deEnrichment(degsi, count, pval, saveFile);
	
	return;
	}

/**
 * Sets degs in all the cells up to a threshold level
 * @param cell
 * @param degs
 * @param threshold
 * @param level
 * @return
 */
public ArrayList<Cell> recursiveDESearch(Cell cell, HashSet<String> degs, int threshold, int level, boolean computeEnrichment)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	cell.term.resetEnrichment();
	if(cell.subcells!=null && cell.subcells.length>0)
		{
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveDESearch(c, degs, threshold, level+1, computeEnrichment);
			retList.addAll(sublist);
			if(sublist.size()>0)	
				add=true;
			}
		if((add && level<=this.levelThreshold))
			{
			cell.term.setDegs(degs, computeEnrichment, expData.getNumGenes());
			if(computeEnrichment || cell.term.degs.size()>=threshold)//NEW
				retList.add(cell);
			}
		}
	
	cell.term.setDegs(degs, computeEnrichment, expData.getNumGenes());
	if(computeEnrichment || cell.term.degs.size()>=threshold)//N EW
		retList.add(cell);
	return retList;
	}

public ArrayList<Cell> recursivePvalSearch(Cell cell, double threshold, int level)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	if(cell.subcells==null || cell.subcells.length==0) //leaf node
		{
		cell.searchedColor=null;
		if(cell.term.pvalue<=threshold)
			{
			//System.out.println(cell.term.name+"\t"+cell.term.pvalue+"\t"+cell.term.degs.size());

			cell.searchedColor=new Color(0,200,0);
			retList.add(cell);
			cell.term.significant=true;
			}
		}
	else		//non leaf-node
		{
		boolean contains=false;
		if(cell.term.pvalue<=threshold)
			{
			contains=true;
			cell.term.significant=true;
			}
	
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursivePvalSearch(c,threshold, level+1);
			retList.addAll(sublist);
			if(sublist.size()>0)	add=true;
			}
		
		//count subterms
		HashSet<String> set=new HashSet<String>();
		for(Cell c:retList)	if(c.term.significant)	set.add(c.term.id);
		cell.term.numRelevantSubterms=set.size();
		
		
		if(add || contains)
			{
			if(contains)
				cell.searchedColor=new Color(0,200,0);
			else
				{
				cell.searchedColor=new Color(0,100,0);
				}
			retList.add(cell);
			}
		}
	return retList;
	}

/**
 * Searches terms that are differentially expressed on conditions g1 respect to conditions g2
 * Type tells if they must be upreg on g1 respect to g2 (0) or the other way round (1)
 * @param g1
 * @param g2
 * @param type
 */
public void deGeneSearch(int[] g1, int[] g2, int type, float threshold, int termThreshold, String saveFile) {
	searchedText="";

	if(deSearchFrame!=null)	deSearchFrame.setVisible(false);
	System.out.println("Expression Search...");
	searchedCells.clear();
	
	if(g1==null || g1.length==0 || g2==null || g2.length==0)	{redraw(); return;}
	
	//1) Get genes that are high in g1 respect to g2 (more than threshold)
	degsi=expData.getGeneNames(expData.differentialExpressionGenes(g1, g2, threshold, type));
	this.termThreshold=termThreshold;

	System.out.println(degsi.size()+" differential expressed genes found");
		
	for(Cell c:cells)
		searchedCells.addAll(recursiveDESearch(c, degsi, termThreshold, 0, false));
	
	System.out.println("Found "+searchedCells.size()+" difexp terms");
	
	if(saveFile!=null)
		this.exportCells(searchedCells, saveFile, false);
	
	redraw();
	if(gh!=null)	gh.redraw();

	return;
	}




public ArrayList<Cell> recursiveDESearch2(Cell cell, int[] g1, int[] g2, int type, float threshold, int level)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	cell.term.resetEnrichment();
	
	if(cell.subcells!=null && cell.subcells.length>0)
		{
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveDESearch2(c,g1, g2, type, threshold, level+1);
			retList.addAll(sublist);
			if(sublist.size()>0)	add=true;
			}
		if((add && level<=this.levelThreshold))
			{
			cell.searchedColor=new Color(0,100,0);
			retList.add(cell);
			}
		}
	double e1=0;
	for(int i:g1)			e1+=cell.expressionLevel.get(i);
	double e2=0;
	for(int i:g2)			e2+=cell.expressionLevel.get(i);
	e1/=g1.length; e2/=g2.length;
	if(g2.length==0)	//single group (1)
		{
		if(type==0 && e1>threshold)	
			{
			cell.searchedColor=new Color(0,200,0);
			retList.add(cell);
			}
		else if(type==1 && e1<threshold)
			{
			cell.searchedColor=new Color(0,200,0);
			retList.add(cell);
			}
		}
	else{
		if(g1.length==0)	//single group (2)
			{
			
			}
		else		//normal case: 2 groups
			{
			if(type==0 && e1>e2+threshold)	
				{
				cell.searchedColor=new Color(0,200,0);
				retList.add(cell);
				}
			else if(type==1 && e2>e1+threshold)
				{
				cell.searchedColor=new Color(0,200,0);
				retList.add(cell);
				}
			}
		}
	
	return retList;
	}

public ArrayList<Cell> recursiveSearch(Cell cell, String searchText, int level)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	cell.term.resetEnrichment();
	
	//leaf node
	if(cell.subcells==null || cell.subcells.length==0)
		{
		cell.searchedColor=null;
		if(cell.term.name.toLowerCase().contains(searchText.toLowerCase()) || cell.term.id.toLowerCase().contains(searchText.toLowerCase()))
			{
			cell.searchedColor=new Color(0,200,0);
			retList.add(cell);
			}
		for(String id:cell.term.geneIds)
			{
			if(id.toLowerCase().contains(searchText.toLowerCase()))	
				{
				cell.searchedColor=new Color(0,200,0);
				if(!retList.contains(cell))	retList.add(cell);
				cell.term.degs.add(id.substring(id.indexOf(":")+1).toLowerCase());
				degsi.add(id.substring(id.indexOf(":")+1).toLowerCase());
				}
			else
				{
				String syn0, syn1, syn2, id2=null;
				if(type==VoronoiVisualization.KEGG)
					{
					id2=id.substring(id.indexOf(":")+1).toLowerCase();
					syn0=expData.getSynonym(id2, 0);
					syn1=expData.getSynonym(id2, 1);
					syn2=expData.getSynonym(id2, 2);
					}
				else
					{
					syn0=expData.getSynonym(id, 0);
					syn1=expData.getSynonym(id, 1);
					syn2=expData.getSynonym(id, 2);
					}
				if((syn0!=null && syn0.toLowerCase().contains(searchText.toLowerCase())) ||
						(syn1!=null && syn1.toLowerCase().contains(searchText.toLowerCase())) || 
						(syn2!=null && syn2.toLowerCase().contains(searchText.toLowerCase())))
					{
					cell.searchedColor=new Color(0,200,0);
					if(!retList.contains(cell))	retList.add(cell); 
					if(type==VoronoiVisualization.KEGG)	{cell.term.degs.add(id2);degsi.add(id2);}
					else								{cell.term.degs.add(id);degsi.add(id);}
					}
				}
			}
 		}
	else	//internal node
		{
		boolean contains=false;
		if(cell.term.name.toLowerCase().contains(searchText.toLowerCase()) || cell.term.id.toLowerCase().contains(searchText.toLowerCase()))
			contains=true;
		
		for(String id:cell.term.geneIds)
			{
			if(id.toLowerCase().contains(searchText.toLowerCase()))	
				{
				contains=true;
				cell.term.degs.add(id.substring(id.indexOf(":")+1).toLowerCase());
				}
			else
				{
				String syn0, syn1, syn2, id2=null;
				if(type==VoronoiVisualization.KEGG)
					{
					id2=id.substring(id.indexOf(":")+1).toLowerCase();
					syn0=expData.getSynonym(id2, 0);
					syn1=expData.getSynonym(id2, 1);
					syn2=expData.getSynonym(id2, 2);
					}
				else
					{
					syn0=expData.getSynonym(id, 0);
					syn1=expData.getSynonym(id, 1);
					syn2=expData.getSynonym(id, 2);
					}
				if((syn0!=null && syn0.toLowerCase().contains(searchText.toLowerCase())) ||
						(syn1!=null && syn1.toLowerCase().contains(searchText.toLowerCase())) || 
						(syn2!=null && syn2.toLowerCase().contains(searchText.toLowerCase())))
					{
					cell.searchedColor=new Color(0,200,0);
					if(!retList.contains(cell))	retList.add(cell); 
					if(type==VoronoiVisualization.KEGG)	{cell.term.degs.add(id2);degsi.add(id2);}
					else								{cell.term.degs.add(id);degsi.add(id);}
					}
				}
			}
		
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveSearch(c,searchText, level+1);
			retList.addAll(sublist);
			if(sublist.size()>0)	add=true;
			}
	
		if((add && level<=this.levelThreshold) || contains)
			{
			if(contains)
				cell.searchedColor=new Color(0,200,0);
			else
				cell.searchedColor=new Color(0,100,0);
			retList.add(cell);
			}
		}
	return retList;
	}

public ArrayList<Cell> recursiveEnrichedAnnotation(Cell cell, SimpleChartRecord[] records, HashMap<String, String> conversionTable, int level)
	{
	ArrayList<Cell> retList=new ArrayList<Cell>();
	cell.term.resetEnrichment();
	if(cell.subcells==null || cell.subcells.length==0)
		{
		cell.searchedColor=null;
		boolean found=false;
		for(SimpleChartRecord r:records)
			{
			String recordName=r.getTermName().substring(r.getTermName().indexOf(":")+1);
			if(recordName.indexOf("~")>-1)	recordName=recordName.substring(recordName.indexOf("~")+1);
			if(recordName.equalsIgnoreCase(cell.term.name.trim()))
				{
				found=true;
				cell.term.pvalueBH=r.getBenjamini();
				
				ArrayList<String> genes=new ArrayList<String>();
				genes.addAll(Arrays.asList(r.getGeneIds().toLowerCase().split(", ")));
				if(conversionTable!=null)
					{
					ArrayList<String> conv=new ArrayList<String>();
					for(String s:genes)	
						{
						String c=conversionTable.get(s); 
						if(c!=null)	conv.add(c);
						}
					genes=conv;
					}
				
				cell.term.degs=new HashSet<String>();
				cell.term.degs.addAll(CollectionUtils.intersection(genes, cell.term.geneExs.keySet()));
				break;
				}
			}
		if(found)
			{
			cell.searchedColor=new Color(0,200,0);
			retList.add(cell);
			}
		else
			cell.term.resetEnrichment();
		}
	else		//non leaf-node
		{
		boolean contains=false;
		
		for(SimpleChartRecord r:records)
			{
			String recordName=r.getTermName().substring(r.getTermName().indexOf(":")+1);
			if(recordName.indexOf("~")>-1)	recordName=recordName.substring(recordName.indexOf("~")+1);
			if(recordName.equalsIgnoreCase(cell.term.name.trim()))
				{
				contains=true;
				cell.term.pvalueBH=r.getBenjamini();
				
				ArrayList<String> genes=new ArrayList<String>();
				genes.addAll(Arrays.asList(r.getGeneIds().toLowerCase().split(", ")));
				if(conversionTable!=null)
					{
					ArrayList<String> conv=new ArrayList<String>();
					for(String s:genes)	
						{
						String c=conversionTable.get(s); 
						if(c!=null)	conv.add(c);
						}
					genes=conv;
					}
				
				cell.term.degs=new HashSet<String>();
				cell.term.degs.addAll(CollectionUtils.intersection(genes, cell.term.geneExs.keySet()));
				break;
				}
			else
				cell.term.resetEnrichment();
			}
		
		boolean add=false;
		for(Cell c:cell.subcells)
			{
			ArrayList<Cell> sublist=recursiveEnrichedAnnotation(c,records, conversionTable, level+1);
			retList.addAll(sublist);
			if(sublist.size()>0)	add=true;
			}
		
		//count subterms
		HashSet<String> set=new HashSet<String>();
		for(Cell c:retList)	//if(c.term.pvalueBH!=999)	set.add(c.term.id);
			if(c.term.significant)	set.add(c.term.id);
		cell.term.numRelevantSubterms=set.size();
		
		
		//if((add && level<=this.levelThreshold) || contains)
		if(add || contains)
			{
			if(contains)
				cell.searchedColor=new Color(0,200,0);
			else
				{
				cell.searchedColor=new Color(0,100,0);
				}
			retList.add(cell);
			}
		}
	return retList;
	}
}
