package es.usal.voronto.view;

import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;

import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;

import javax.swing.JLabel;

import java.util.Collections;
import es.usal.voronto.control.Voronto;
import es.usal.voronto.model.filters.ExpressionDataFilter;
import es.usal.voronto.model.filters.GAFDataFilter;
import es.usal.voronto.model.filters.OBODataFilter;
import es.usal.voronto.model.ontology.GOparser;
import es.usal.voronto.model.ontology.OntologyTerm;

import java.awt.Font;


/**
* This code was edited or generated using CloudGarden's Jigloo
* SWT/Swing GUI Builder, which is free for non-commercial
* use. If Jigloo is being used commercially (ie, by a corporation,
* company or business for any purpose whatever) then you
* should purchase a license for each developer using Jigloo.
* Please visit www.cloudgarden.com for details.
* Use of Jigloo implies acceptance of these licensing terms.
* A COMMERCIAL LICENSE HAS NOT BEEN PURCHASED FOR
* THIS MACHINE, SO JIGLOO OR THIS CODE CANNOT BE USED
* LEGALLY FOR ANY CORPORATE OR COMMERCIAL PURPOSE.
*/
public class Introduction extends javax.swing.JPanel implements MouseListener{
	/**
	 * 
	 */
	private static final long serialVersionUID = 2715114188900530489L;
	private JLabel jLabel1;
	private JLabel jLabel2;
	private JLabel jLabel3;
	private JLabel jLabel4;
	private JButton jButton1;
	private JLabel jLabel5;
	private JButton jButton2;
	private JLabel jLabel6;
	
	public JComboBox jComboBox1;
	public File expressionFile;
	public File annotationFile;
	private File ontologyFile;
	public String ontologyType;
	private Voronto parent;
	private JLabel jLabel8;
	public JLabel jLabel9;
	private JLabel jLabel10;

	private String onto;
	private TermList tl;

	/**
	* Auto-generated main method to display this 
	* JPanel inside a new JFrame.
	*/
		
	public Introduction(Voronto p) {
		super();
		parent=p;
		initGUI();
	}
	
	private void initGUI() {
		try {
			this.setPreferredSize(new Dimension(507, 334));
			this.setLayout(null);
			this.setBackground(new java.awt.Color(255,255,255));
			this.setSize(700, 400);
			{
				jLabel1 = new JLabel();
				this.add(jLabel1);
				jLabel1.setText("Mapping gene expression to ontologies with voronoi tessellations");
				jLabel1.setBounds(194, 72, 357, 17);
				jLabel1.setFont(new java.awt.Font("Dialog",2,9));
			}
			{
				jLabel2 = new JLabel();
				this.add(jLabel2);
				jLabel2.setIcon(new ImageIcon(getClass().getClassLoader().getResource("es/usal/voronto/img/voronto.png")));
				jLabel2.setBounds(283, 17, 143, 55);
			}
			{
				jLabel3 = new JLabel();
				this.add(jLabel3);
				jLabel3.setText("1) Select gene expression data");
				jLabel3.setBounds(124, 124, 241, 13);
			}
			{
				jLabel4 = new JLabel();
				this.add(jLabel4);
				jLabel4.setText("Click here for sample data and format specs");
				jLabel4.setBounds(136, 146, 267, 16);
				jLabel4.setFont(new java.awt.Font("Dialog",0,10));
				jLabel4.addMouseListener(this);
				
			}
			{
				jButton1 = new JButton();
				this.add(jButton1);
				jButton1.setText("Select...");
				jButton1.setBounds(390, 124, 100, 24);
				jButton1.setBackground(new java.awt.Color(255,255,255));
				jButton1.setSize(98, 22);
				jButton1.addActionListener(new java.awt.event.ActionListener()	{
					public void actionPerformed(java.awt.event.ActionEvent e) 
						{
						JFileChooser selecFile = new JFileChooser();
						ExpressionDataFilter edf=new ExpressionDataFilter();
						selecFile.addChoosableFileFilter(edf);
						selecFile.setFileFilter(edf);
						
						String path=null;
						try{
							BufferedReader br=new BufferedReader(new FileReader(parent.expressionPath));
							path=br.readLine();
							br.close();
							}catch(Exception ex){System.out.println("pathFile non existing yet");}
						if(path!=null)				
							selecFile.setCurrentDirectory(new File(path));
						else						
							selecFile.setCurrentDirectory(new File("."));
						
						int returnval = selecFile.showDialog((Component)e.getSource(), "Load expression data");
						
						if(returnval == JFileChooser.APPROVE_OPTION) 
							{
							expressionFile = selecFile.getSelectedFile();
							jLabel8.setText(expressionFile.getAbsolutePath());
							jButton2.setEnabled(true);
							}
					}
					});
			}
			
			{
				jLabel5 = new JLabel();
				this.add(jLabel5);
				jLabel5.setText("2) Select ontology");
				jLabel5.setBounds(124, 218, 162, 14);
			}
			{
				ComboBoxModel jComboBox1Model = 
						new DefaultComboBoxModel(
								//		new String[] { "KEGG", "GO slim BP","GO slim CC", "GO full BP", "REACTOME", "Custom..." });
										new String[] { "KEGG", "GO BP","GO CC", "GO MF", "REACTOME", "Custom..." });
								
				jComboBox1 = new JComboBox();
				this.add(jComboBox1);
				jComboBox1.setModel(jComboBox1Model);
				jComboBox1.setBounds(380, 214, 130, 22);
				jComboBox1.setBackground(new java.awt.Color(255,255,255));
				jComboBox1.setSelectedIndex(0);
				jComboBox1.addActionListener(new java.awt.event.ActionListener() {

				

					@Override
					public void actionPerformed(ActionEvent e) 
						{
						switch(jComboBox1.getSelectedIndex())
							{
							case 5:
								JFileChooser selecFile = new JFileChooser();
								
								OBODataFilter odf=new OBODataFilter();
								selecFile.addChoosableFileFilter(odf);
								
								selecFile.setFileFilter(odf);
								
								String path=null;
								try{
									BufferedReader br=new BufferedReader(new FileReader(parent.ontologyPath));
									path=br.readLine();
									br.close();
									}catch(Exception ex){System.out.println("pathFile non existing yet");}
								if(path!=null)				
									selecFile.setCurrentDirectory(new File(path));
								else						
									selecFile.setCurrentDirectory(new File("."));
								
								int returnval = selecFile.showDialog((Component)e.getSource(), "Load ontology");
								
								if(returnval == JFileChooser.APPROVE_OPTION) 
									{
									ontologyFile = selecFile.getSelectedFile();
									ontologyType = selecFile.getFileFilter().getDescription();
									jLabel9.setText(ontologyFile.getAbsolutePath());
									
									if(ontologyType.contains("obo"))
										{
										jLabel10.setText("Click here to load an annotation file for OBO ontology...");
										jLabel10.setVisible(true);
										}
									}
							break;
						case 1://select root
							parent.root="biological_process";
							onto=parent.root;
							tl=null;
							jLabel9.setText("Root (click to select): "+parent.root);
							jLabel10.setText("Click here to load a custom .GAF file from GOA...");
							jLabel10.setVisible(true);
							break;
						case 2://select root
							parent.root="cellular_component";
							onto=parent.root;
							tl=null;
							jLabel9.setText("Root (click to select): "+parent.root);
							jLabel10.setText("Click here to load a custom .GAF file from GOA...");
							jLabel10.setVisible(true);
							break;
						case 3://select root
							parent.root="molecular_function";
							onto=parent.root;
							tl=null;
							jLabel9.setText("Root (click to select): "+parent.root);
							jLabel10.setText("Click here to load a custom .GAF file from GOA...");
							jLabel10.setVisible(true);
							break;
						default:
							ontologyFile = null;
							jLabel9.setText("");
							jLabel10.setText("");
							break;
						}
					}
					});
			}
			{
				jLabel6 = new JLabel();
				this.add(jLabel6);
				jLabel6.setText("3) Launch application!");
				jLabel6.setBounds(124, 309, 162, 19);
				
			}
			{
				jButton2 = new JButton();
				this.add(jButton2);
				jButton2.setText("Launch");
				jButton2.setBackground(new java.awt.Color(255,255,255));
				jButton2.setBounds(390, 306, 100, 22);
				//jButton2.setEnabled(false);	//TODO: commented just for tests!
				jButton2.addActionListener(new java.awt.event.ActionListener() {
					public void actionPerformed(java.awt.event.ActionEvent e) 
						{
						if(expressionFile==null)	
							{
							if(ontologyFile==null)	parent.launch(null, null, jComboBox1.getSelectedIndex());
							else					parent.launch(null, ontologyFile.getPath(), jComboBox1.getSelectedIndex());
							}
						else
							{
							if(ontologyFile==null)	parent.launch(expressionFile, null, jComboBox1.getSelectedIndex());
							else					parent.launch(expressionFile, ontologyFile.getPath(), jComboBox1.getSelectedIndex());
							}
						}
					});
			}
			{
				jLabel8 = new JLabel();
				this.add(jLabel8);
				jLabel8.setText("No file selected");
				jLabel8.setBounds(136, 166, 431, 15);
				jLabel8.setFont(new java.awt.Font("Dialog",0,10));
				jLabel8.setForeground(new Color(150,150,255));
			}
			{
				jLabel9 = new JLabel();
				this.add(jLabel9);
				jLabel9.setText("");
				jLabel9.setBounds(136, 238, 431, 15);
				jLabel9.setFont(new java.awt.Font("Dialog",0,10));
				jLabel9.setForeground(new Color(150,150,255));
				jLabel9.addMouseListener(this);
			}
			
			jLabel10 = new JLabel("Click here to load an annotation file for OBO ontology...");
			jLabel10.setForeground(Color.BLUE);
			jLabel10.setFont(new Font("Lucida Grande", Font.PLAIN, 10));
			jLabel10.setBounds(136, 255, 425, 16);
			jLabel10.setVisible(false);
			jLabel10.addMouseListener(this);
			add(jLabel10);
			

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	
	
	
	
	
	


public void reset()
	{
	//jLabel8.setText("");
	//jLabel9.setText("");
	//if(this.expressionFile==null)	jButton2.setEnabled(false);
	}

@Override
public void mouseClicked(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mouseEntered(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mouseExited(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mousePressed(MouseEvent arg0) {
	// TODO Auto-generated method stub
	
}

@Override
public void mouseReleased(MouseEvent arg0) {
	if(arg0.getSource().getClass()==JLabel.class)
		{
		if((JLabel)(arg0.getSource())==jLabel4)
			{
			try{
				java.awt.Desktop.getDesktop().browse(java.net.URI.create("http://vis.usal.es/~visusal/voronto/voronto/Help.html"));
			}catch(IOException ioe){ioe.printStackTrace();}
			}
		if((JLabel)(arg0.getSource())==jLabel10)
		{
		try{
			JFileChooser selecFileGAF = new JFileChooser();
			GAFDataFilter gdf=new GAFDataFilter();
			selecFileGAF.addChoosableFileFilter(gdf);
			selecFileGAF.setFileFilter(gdf);
			
			String path=null;
			try{
				BufferedReader br=new BufferedReader(new FileReader(parent.ontologyPath));
				path=br.readLine();
				}catch(Exception ex){System.out.println("pathFile non existing yet");}
			if(path!=null)		selecFileGAF.setCurrentDirectory(new File(path));
			else				selecFileGAF.setCurrentDirectory(new File("."));
			
			int returnval = selecFileGAF.showDialog(this, "Load annotations");
			if(returnval == JFileChooser.APPROVE_OPTION) 
				{
				annotationFile = selecFileGAF.getSelectedFile();
				jLabel10.setText("Annotation file: "+annotationFile.getAbsolutePath());
				}
		
			}catch(Exception ioe){ioe.printStackTrace();}
		}
		if((JLabel)(arg0.getSource())==jLabel9 && jComboBox1.getSelectedIndex()>0 && jComboBox1.getSelectedIndex()<4)
		{
		try{
			if(tl==null)
    			{
				Cursor hourglassCursor = new Cursor(Cursor.WAIT_CURSOR);
				setCursor(hourglassCursor);
				TreeMap<OntologyTerm, TreeMap> m=GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", onto, false);
		    	ArrayList<String> lista=parent.getTermList(m);
		    	
		    	Set<String> set=new HashSet<String>(lista);	//remove duplicates
		    	lista=new ArrayList<String>(set);
		    	Collections.sort(lista);					//sort
		    	System.out.println("Number of unique terms: "+lista.size());
		    	
		    	tl=new TermList(lista, this, parent);
				tl.setLocation(parent.getX()+(int)((parent.getWidth()-tl.getWidth())*0.5),parent.getY()+(int)((parent.getHeight()-tl.getHeight())*0.5));
				
				Cursor normalCursor = new Cursor(Cursor.DEFAULT_CURSOR);
				setCursor(normalCursor);
				}
    		
	    	tl.setVisible(true);
	    	}catch(Exception ioe){ioe.printStackTrace();}
		}

		}
}
}