package es.usal.voronto.view;

import java.awt.Color;
import java.awt.Toolkit;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.philhosoft.p8g.svg.P8gGraphicsSVG;

import ch.usi.inf.sape.hac.HierarchicalAgglomerativeClusterer;
import ch.usi.inf.sape.hac.agglomeration.AgglomerationMethod;
import ch.usi.inf.sape.hac.agglomeration.CentroidLinkage;
import ch.usi.inf.sape.hac.dendrogram.Dendrogram;
import ch.usi.inf.sape.hac.dendrogram.DendrogramBuilder;
import ch.usi.inf.sape.hac.dendrogram.DendrogramNode;
import ch.usi.inf.sape.hac.dendrogram.ObservationNode;
import ch.usi.inf.sape.hac.experiment.DissimilarityMeasure;

import es.usal.voronto.model.clustering.ExpressionSubset;
import es.usal.voronto.model.filters.TIFFileFilter;
import es.usal.voronto.model.filters.SVGFileFilter;
import es.usal.voronto.model.voronoi.Cell;

import processing.core.PApplet;
import processing.core.PFont;
import processing.core.PGraphics;

/**
 * Class to display static, simple heatmaps for the genes in a given term.
 * It is invoked from VoronoiVisualization when a term is double clicked with Alt key down
 * TODO: make a previous hierarchical clustering on rows.
 * @author rodri
 *
 */
public class CellHeatmap extends PApplet
//public class CellHeatmap implements PSurface
	{
	/**
	 * 
	 */
	private static final long serialVersionUID = -4416622212260260249L;
	Cell cell=null;
	int margin=10;//margin at the sides of the Frame
	int marginRows=100;//space for gene names		
	int marginCols=100;//space for condition names	
	int numRows=-1;
	int numCols=-1;
	VoronoiVisualization vv=null;
	String hoveredGene=null;
	int hoveredCondition=-1;
	
	PFont font;
	int size=12;//size of each expression level TODO: check if it becomes superlarge...
	private int maxItemsHeight;
	private int basicWidth;
	private ArrayList<Integer> order;
	private static int titleHeight=22;
	
	private int scaleFactor=1;
	private ExpressionSubset subset;
	private int fontSize=14;
	
	
	public CellHeatmap(Cell c, PFont font, VoronoiVisualization v)
		{
		super();
		vv=v;
		cell=c;
		scaleFactor=1;
		if(cell.term.geneExs==null || cell.term.geneExs.size()==0)
			{
			System.err.println("No single gene expression for this term, please choose another one");
			return;
			}
		numRows=cell.term.geneExs.size();
		
		//cell.completeWithChildren();
		numCols=cell.term.geneExs.get(cell.term.geneExs.keySet().iterator().next()).size();
		
		this.font=font;
		this.fontSize=vv.fontSize;
		
		//init();
		computeSize();
		this.setSize(width, height);
		computeOrder();
		/*setup();*/
		//start();
		}
	
	
	public void computeSize()
		{
		marginRows=0;
		for(String g:cell.term.geneExs.keySet())
			{
		//	marginRows=Math.max((int)textWidth(g),marginRows);
			marginRows=Math.max(g.length()*8,marginRows);
			}
		marginCols=0;
		for(String cond:vv.expData.columnLabels)
			{
			//marginCols=Math.max((int)textWidth(cond),marginCols);
			marginCols=Math.max(cond.length()*8,marginCols);
			}
				

		basicWidth=width=marginRows+margin*2+size*numCols;
		this.height=marginCols+margin*2+size*numRows;
		int screenHeight=Toolkit.getDefaultToolkit().getScreenSize().height;
		maxItemsHeight=(int)Math.floor((screenHeight-marginCols-margin*2-(double)titleHeight)/size)-2;
		
		if(height>screenHeight-titleHeight)
			{
			height=marginCols+margin*2+size*maxItemsHeight;
					
			int numItems=numRows;
			while(numItems>maxItemsHeight)
				{
				width+=basicWidth;
				numItems-=maxItemsHeight;
				}
			
			if(width>Toolkit.getDefaultToolkit().getScreenSize().width*4)
				{
				System.err.println("Too many genes, choose a smaller term");
				//setVisible(false);
				JOptionPane.showMessageDialog(null, "This term has too many genes for the heatmap visualization. Please choose a smaller term.", "Too many genes", JOptionPane.INFORMATION_MESSAGE);
				stop();
				return;
				}
			if(width>Toolkit.getDefaultToolkit().getScreenSize().width*2)
				{
				/*this.setSize(vv.getWidth()-margin, vv.getHeight()-margin);
				this.setPreferredSize(new Dimension(vv.getWidth()-margin, vv.getHeight()-margin));
				width=vv.getWidth()-margin;
				height=vv.getHeight()-margin;*/
				this.getSurface().setSize(width, height);
				return;
				}
			}
		}
	
	public void setup()
		{
		size(width, height);
		setSize(width, height);
		/*setPreferredSize(new Dimension(width, height));
		setMinimumSize(new Dimension(width, height));
		setMaximumSize(new Dimension(width, height));*/
		this.getSurface().setSize(width, height);
		
		
		smooth();
		noLoop();
		}

	/* This method is not considering which side is better to add first in the 
	 * conext of other branches, it always goes first for left, then for right*/
	public void retrieveOrder(ArrayList<Integer> order, DendrogramNode dn)
		{
		if(dn.getLeft().getClass().toString().contains("ObservationNode"))
			order.add(((ObservationNode)dn.getLeft()).getObservation());
		else
			retrieveOrder(order, dn.getLeft());
		
		if(dn.getRight().getClass().toString().contains("ObservationNode"))
			order.add(((ObservationNode)dn.getRight()).getObservation());
		else
			retrieveOrder(order, dn.getRight());
		}
	
	public void sortDendrogram(DendrogramNode dn)
		{
		if(dn.getClass().toString().contains("ObservationNode"))	
			return;
		//If both are leaf nodes, nothing to do
		if(dn.getLeft().getClass().toString().contains("ObservationNode") && dn.getRight().getClass().toString().contains("ObservationNode"))
			return;
		else
			{
			
			//Check order of the two left branches depending on the right node
			if(!dn.getLeft().getClass().toString().contains("ObervationNode"))
				{
				if(dn.getLeft().getLeft()!=null && dn.getLeft().getRight()!=null)
					{
					ArrayList<Float> rp=getAverageProfile(dn.getRight());
					ArrayList<Float> llp=getAverageProfile(dn.getLeft().getLeft());
					ArrayList<Float> lrp=getAverageProfile(dn.getLeft().getRight());
					double dll=subset.computeDissimilarity(rp, llp);
					double dlr=subset.computeDissimilarity(rp, lrp);
					if(dll<dlr)	//switch the branches
						{
						DendrogramNode temp=dn.getLeft().getRight();
						dn.getLeft().setRight(dn.getLeft().getLeft());
						dn.getLeft().setLeft(temp);
						}
					sortDendrogram(dn.getLeft());
					if(!dn.getRight().getClass().toString().contains("ObervationNode"))
						sortDendrogram(dn.getRight());
					}
				}
			//Check order of the two right branches depending on left node
			else
				{
				ArrayList<Float> lp=getAverageProfile(dn.getLeft());
				ArrayList<Float> rlp=getAverageProfile(dn.getRight().getLeft());
				ArrayList<Float> rrp=getAverageProfile(dn.getRight().getRight());
				double drl=subset.computeDissimilarity(lp, rlp);
				double drr=subset.computeDissimilarity(lp, rrp);
				if(drr<drl)	//switch the branches
					{
					DendrogramNode temp=dn.getRight().getLeft();
					dn.getRight().setLeft(dn.getRight().getRight());
					dn.getRight().setRight(temp);
					}
				sortDendrogram(dn.getRight());
				}
			}
		}

	public ArrayList<Float> getAverageProfile(DendrogramNode dn)
		{
		String[] names=cell.term.geneExs.keySet().toArray(new String[0]);
		
		if(dn.getClass().toString().contains("ObservationNode"))
			return cell.term.geneExs.get(names[((ObservationNode)dn).getObservation()]);
		else
			{
			ArrayList<Float> lp=getAverageProfile(dn.getRight());
			ArrayList<Float> rp=getAverageProfile(dn.getLeft());
			ArrayList<Float> ap=new ArrayList<Float>();
			for(int i=0;i<lp.size();i++)
				ap.add((float)((lp.get(i)+rp.get(i))/2.0));
			return ap;
			}
		}
		
	
	public void computeOrder()
		{
		String[] names=cell.term.geneExs.keySet().toArray(new String[0]);
		order=new ArrayList<Integer>();
		
		if(names.length>1)
			{
			subset= new ExpressionSubset(cell.term.geneExs, names);
			DissimilarityMeasure dissimilarityMeasure = (DissimilarityMeasure)subset;
			AgglomerationMethod agglomerationMethod = new CentroidLinkage();
			//AgglomerationMethod agglomerationMethod = new WeightedAverageLinkage();
			//AgglomerationMethod agglomerationMethod = new CompleteLinkage();
			DendrogramBuilder dendrogramBuilder = new DendrogramBuilder(subset.getNumberOfObservations());
			HierarchicalAgglomerativeClusterer clusterer = new HierarchicalAgglomerativeClusterer(subset, dissimilarityMeasure, agglomerationMethod);
			clusterer.cluster(dendrogramBuilder);
			Dendrogram dendrogram = dendrogramBuilder.getDendrogram();
			
			//dendrogram.dump();
			
			DendrogramNode dn=dendrogram.getRoot();
			sortDendrogram(dn);
			retrieveOrder(order, dn);
			}
		else
			{
			order.add(0);
			}
		}
	
	
	public void mouseMoved() 
		{
		redraw();
		}
	
	public void mouseReleased()
		{
		if(mouseEvent.getClickCount()==2 && hoveredGene!=null)
			{
			try{
			String entrezLabel=vv.expData.getSynonym(hoveredGene, vv.expData.ENTREZ);
				
			if(vv.expData.chip.equals("entrezgene"))
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+hoveredGene));
			else if(entrezLabel!=null)
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+entrezLabel));
			else
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://www.ncbi.nlm.nih.gov/gene?term="+hoveredGene.toUpperCase()+"%20AND%20"+vv.expData.organism.replace(" ", "%20")+"%5BOrganism%5D"));
			}catch(IOException e){System.out.println("Error: cannot show webpage: "+e.getMessage()); e.printStackTrace();}
			}
			
		}
	
	public void keyReleased()
		{
		
		switch(key)
			{
			case 'e':	//export gene ids 
				vv.export(cell.term);
				break;
			case 'p':	//print hi-res image
				printImage();
				break;
			default:
				vv.key=this.key;
				vv.keyReleased();
				break;
			}
			
		}
	
	public void printImage()
		{
		JFileChooser selecFile = new JFileChooser();
		TIFFileFilter png=new TIFFileFilter();
		selecFile.addChoosableFileFilter(png);
		SVGFileFilter svg=new SVGFileFilter();
		selecFile.addChoosableFileFilter(svg);
		
		if(vv.imagePath!=null)
			selecFile.setCurrentDirectory(new File(vv.imagePath));
		else if(vv.expData!=null)
			selecFile.setCurrentDirectory(new File(vv.expData.filePath));

		if(vv.expData!=null)	
			selecFile.setSelectedFile(new File(cell.term.name.trim()));
		else
			selecFile.setSelectedFile(new File(cell.term.name.trim()));
		selecFile.setAcceptAllFileFilterUsed(false);
		
		int returnval = selecFile.showSaveDialog(this.vv.frame);

		if(returnval==JFileChooser.APPROVE_OPTION)
			{
			try{Thread.sleep(100);}catch(Exception e){}
			
			if(selecFile.getFileFilter()==png)	//for TIFF
				{
				scaleFactor=3;
				PGraphics hires = createGraphics(width*scaleFactor, height*scaleFactor, JAVA2D);
				beginRecord(hires);
				hires.scale(scaleFactor);
				draw();
				endRecord();
				hires.save(selecFile.getSelectedFile().getAbsolutePath());
				scaleFactor=1;
				}
			else if(selecFile.getFileFilter()==svg)	//for SVG (letter sizes should change to SansSerif to avoid losing them)
				{
				font=loadFont("es/usal/voronto/font/SansSerif-20.vlw");
				fontSize=11;
				P8gGraphicsSVG hires = (P8gGraphicsSVG)createGraphics(width, height, P8gGraphicsSVG.SVG);
				beginRecord(hires);
				draw();
				hires.endRecord(selecFile.getSelectedFile().getAbsolutePath());
				//font=loadFont("es/usal/voronto/font/AppleSymbols-14.vlw");
				font=loadFont("es/usal/voronto/font/Arial-14.vlw");
				fontSize=14;
				}
			
			try{
				vv.imagePath=selecFile.getSelectedFile().getParent();
				BufferedWriter bw=new BufferedWriter(new FileWriter(vv.voronto.imgPath));
				bw.write(vv.imagePath);
				bw.close();
				}catch(IOException e){e.printStackTrace();}
				
			redraw();
			}
		}
	
	public void draw()
		{
		synchronized(vv)
		{
		if(order==null)	return;
		fill(255,255,255);
		stroke(255);
		rect(0,0,width,height);
		noFill();
		
		hoveredCondition=-1;
		//hoveredGene=null;
		vv.drawOnlyScale=false;
		vv.scaleExpression=(float)vv.minExp-100;
		
		textFont(font, fontSize);
		
		//----------drawHeatmap
		stroke(154);
		
		String[] names=cell.term.geneExs.keySet().toArray(new String[0]);
		int cont=0;
		int xDisplacement=0;
		
		while(cont<numRows)
			{
			if(cont % maxItemsHeight==0)
				{
				if(cont>0)	xDisplacement+=basicWidth;
				
				//draw column labels
				fill(5);
				textAlign(LEFT, CENTER);
				for(int i=0;i<vv.expData.getNumConditions();i++)
					{
					if(mouseX>(margin+marginRows+xDisplacement+i*size) && mouseX<(margin+marginRows+xDisplacement+(i+1)*size))
						fill(0);
					else
						fill(154);
					
					pushMatrix();
					translate((float)(margin+marginRows+xDisplacement+(i+0.5)*size), marginCols);
					rotate((float)(1.5*PI));
					text(vv.expData.conditionNames[i],0, 0);
					
					popMatrix();
					}
				}
			
			textAlign(RIGHT, CENTER);
			
			//draw gene labels and expression levels
			int i=0;
			for(i=0;i<Math.min(numRows, maxItemsHeight) && i+cont<names.length;i++)//for each gene
				{
				String gene=null;
			
				gene=names[order.get(i+cont)];
				if(mouseY>(margin+marginCols+(i)*size) && mouseY<(margin+marginCols+(i+1)*size) && mouseX>xDisplacement && mouseX<xDisplacement+basicWidth)
					fill(0);
				else
					{
					if(cell.term.degs!=null && cell.term.degs.contains(gene))
						fill(0,150,0);
					else
						fill(154);
					}
				
				String geneLabel=vv.expData.getSynonym(gene.toLowerCase(), vv.nameType);
				if(geneLabel==null)	geneLabel=gene;
				//text(gene.substring(gene.indexOf(":")+1).toUpperCase(),marginRows+xDisplacement, (float)(margin+marginCols+(i+0.5)*size));
				text(geneLabel.substring(geneLabel.indexOf(":")+1).toUpperCase(),marginRows+xDisplacement, (float)(margin+marginCols+(i+0.5)*size));
				
				ArrayList<Float> exs=cell.term.geneExs.get(gene);
				for(int j=0;j<numCols;j++)
					{
					Color co=vv.getColor(exs, j);
					fill(co.getRed(), co.getGreen(), co.getBlue());
					
					float xcell=(float)(margin+marginRows+xDisplacement+j*size);
					float ycell=(float)(margin+marginCols+i*size);
					
					rect(xcell,ycell, size,size);
					}//for each col
				}//for each gene
			cont+=i;
			}
		
		//Draw hovered cell
		cont=0;
		xDisplacement=0;
		while(cont<numRows)
			{
			if(cont % maxItemsHeight==0)
				if(cont>0)	xDisplacement+=basicWidth;
		
			int i=0;
			for(i=0;i<Math.min(numRows, maxItemsHeight) && i+cont<names.length;i++)//for each gene
				{
				noFill();
				
				stroke(0);
				
				for(int j=0;j<numCols;j++)
					{
					float xcell=(float)(margin+marginRows+xDisplacement+j*size);
					float ycell=(float)(margin+marginCols+i*size);
					if(j==0 && (mouseY>ycell && mouseY<ycell+size) && mouseX<xcell)
						{
						strokeWeight(3);
						rect(xcell,ycell, size*numCols,size);
						strokeWeight(1);
						vv.redraw();
						hoveredGene=names[order.get(i+cont)];
						return;
						}
					if(i==0 && (mouseX>xcell && mouseX<xcell+size) && mouseY<ycell)
						{
						strokeWeight(3);
						rect(xcell,ycell, size,size*Math.min(maxItemsHeight, numRows-cont));
						strokeWeight(1);
						vv.redraw();
						return;
						}
					if((mouseX>xcell && mouseX<xcell+size) && (mouseY>ycell && mouseY<ycell+size))
						{	
						strokeWeight(3);
						rect(xcell,ycell, size,size);
						strokeWeight(1);
						
						hoveredCondition=j; hoveredGene=names[order.get(i+cont)];
						vv.scaleExpression=cell.term.geneExs.get(hoveredGene).get(hoveredCondition);
						vv.scaleSample=hoveredCondition;
						vv.drawOnlyScale=true;
						vv.redraw();
						return;
						}
					}
				}
			cont+=i;
			}
		
		return;
		}
	}
	}
