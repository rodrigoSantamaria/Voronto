package es.usal.voronto.view;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.ComboBoxModel;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;
import javax.swing.JLabel;
import javax.swing.JComboBox;
import javax.swing.JList;

import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.model.filters.TextFileFilter;
import es.usal.voronto.model.ontology.GeneParser;

import java.awt.Color;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.HashSet;

import javax.swing.JTextField;
import javax.swing.JRadioButton;
import javax.swing.JSeparator;
import javax.swing.JCheckBox;

import java.awt.Font;


public class DifExpSearch extends JDialog {

	/**
	 * 
	 */
	private static final long serialVersionUID = -896238342032232889L;
	private final JPanel contentPanel = new JPanel();
	private JScrollPane list2;
	private VoronoiVisualization v;
	private JComboBox comboBox;
	private JList group1;
	private JList group2;
	private JTextField textField;
	private JTextField textField_1;
	private JTextField textField_2;
	private JCheckBox enrichment;
	private JRadioButton searchGenes;
	private JButton btnSelect;
	private JCheckBox chckbxSaveToFile;
	private String saveFile;
	private JRadioButton searchTerms;
	private String geneFile;
	private JRadioButton loadGenes;
	private JButton btnLoad;

	
	/**
	 * Create the dialog.
	 */
	public DifExpSearch(VoronoiVisualization v) {
		setBounds(100, 100, 412, 600);
		this.v=v;
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBackground(Color.WHITE);
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		/*try {
		contentPanel.setFont(Font.createFont(Font.TRUETYPE_FONT, new File("es/usal/voronto/font/AppleSymbols-14.vlw")));
		}catch(Exception e) {e.printStackTrace();}
		*/
		Font font=new Font("Arial", Font.PLAIN,13);
		Font bfont=new Font("Arial", Font.BOLD,13);
		contentPanel.setFont(font);
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		contentPanel.setLayout(null);
		{
			JLabel lblSearchFor = new JLabel("Search for ");
			lblSearchFor.setBounds(68, 33, 89, 16);
			lblSearchFor.setFont(font);
			contentPanel.add(lblSearchFor);
		}
		
		ComboBoxModel regulationModel = 
				new DefaultComboBoxModel(
						new String[] {"up regulation", "down regulation" });
		comboBox = new JComboBox(regulationModel);
		comboBox.setBounds(144, 29, 165, 27);
		comboBox.setBackground(Color.WHITE);
			
		
		contentPanel.add(comboBox);
		
		JLabel lblOnConditions = new JLabel("on conditions");
		lblOnConditions.setFont(font);
		lblOnConditions.setBounds(38, 74, 165, 16);
		contentPanel.add(lblOnConditions);
		
		JLabel lblRespectToConditions = new JLabel("respect to conditions");
		lblRespectToConditions.setBounds(38, 206, 173, 16);
		lblRespectToConditions.setFont(font);
		contentPanel.add(lblRespectToConditions);
		
		list2 = new JScrollPane(); 
		list2.setBounds(188, 206, 196, 118);
		contentPanel.add(list2);
		
		DefaultComboBoxModel model2 = new DefaultComboBoxModel();
		for(String c:v.expData.conditionNames)		model2.addElement(c);
		
		group2 = new JList();
		group2.setModel(model2);
		
		list2.setViewportView(group2);

		
		JScrollPane list1 = new JScrollPane();
		list1.setBounds(188, 74, 198, 120);
		contentPanel.add(list1);
		
		DefaultComboBoxModel model1 = new DefaultComboBoxModel();
		for(String c:v.expData.conditionNames)		model1.addElement(c);
		
		group1 = new JList();
		group1.setModel(model1);
		
		list1.setViewportView(group1);
		
		JLabel lblThreshold = new JLabel("Differential");
		lblThreshold.setToolTipText("Difference between the average level along conditions in the first list respect to the second");
		lblThreshold.setBounds(38, 338, 89, 16);
		lblThreshold.setFont(font);
		contentPanel.add(lblThreshold);
		
		textField = new JTextField();
		textField.setBounds(198, 332, 67, 28);
		contentPanel.add(textField);
		textField.setColumns(10);
		textField.setText(1.0+"");
		if(v.expData!=null)
			textField.setText(new DecimalFormat("#.##").format(v.expData.meanSd).replace(",", "."));
		
		searchGenes = new JRadioButton("Highlight");
		searchGenes.setForeground(new Color(0, 128, 0));
		searchGenes.setBackground(Color.WHITE);
		searchGenes.setSelected(true);
		searchGenes.setFont(bfont);
		searchGenes.setBounds(6, 400, 100, 23);
		searchGenes.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent arg0) {
						if(searchGenes.isSelected())
							{
							enrichment.setEnabled(true);
							searchTerms.setSelected(false);
							loadGenes.setSelected(false);
							btnLoad.setEnabled(false);
							}
						}
					});
		contentPanel.add(searchGenes);
		
		textField_1 = new JTextField();
		textField_1.setText("2");
		textField_1.setColumns(10);
		textField_1.setBounds(194, 398, 32, 57);
		contentPanel.add(textField_1);
		
		JLabel lblDegs = new JLabel("genes above");
		lblDegs.setBounds(228, 404, 96, 16);
		lblDegs.setFont(font);
		contentPanel.add(lblDegs);
		
		enrichment = new JCheckBox("Highlight");
		enrichment.setForeground(new Color(0, 128, 0));
		enrichment.setBackground(Color.WHITE);
		enrichment.setFont(bfont);
		enrichment.setToolTipText("p-values are Benjamini-corrected, and computed via DAVID Funtional Chart Annotation API");
		enrichment.setBounds(30, 456, 100, 23);
		contentPanel.add(enrichment);
		
		textField_2 = new JTextField();
		textField_2.setText("0.01");
		textField_2.setColumns(10);
		textField_2.setBounds(306, 454, 50, 28);
		contentPanel.add(textField_2);
		
		JSeparator separator = new JSeparator();
		separator.setBounds(10, 486, 348, 12);
		contentPanel.add(separator);
		
		chckbxSaveToFile = new JCheckBox("Save results to file...");
		chckbxSaveToFile.setBounds(6, 499, 173, 23);
		chckbxSaveToFile.setBackground(Color.WHITE);
		
		chckbxSaveToFile.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0)
				{
				if(chckbxSaveToFile.isSelected())
					btnSelect.setEnabled(true);
				else
					btnSelect.setEnabled(false);
				}
		});
		contentPanel.add(chckbxSaveToFile);
		
		btnSelect = new JButton("Select");
		btnSelect.setEnabled(false);
		btnSelect.setBounds(186, 498, 89, 29);
		btnSelect.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0)
				{
				save();
				}
		});
		contentPanel.add(btnSelect);
		
		searchTerms = new JRadioButton("Highlight");
		searchTerms.setForeground(new Color(0, 128, 0));
		searchTerms.setBackground(Color.WHITE);
		searchTerms.setFont(bfont);
		searchTerms.setBounds(6, 362, 100, 23);
		searchTerms.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				if(searchTerms.isSelected())
					{
					searchGenes.setSelected(false);
					loadGenes.setSelected(false);
					enrichment.setEnabled(false);
					btnLoad.setEnabled(false);
					}
				}
			});

		contentPanel.add(searchTerms);
		
		JSeparator separator_1 = new JSeparator();
		separator_1.setBounds(26, 386, 348, 12);
		contentPanel.add(separator_1);
		
		loadGenes = new JRadioButton("Highlight");
		loadGenes.setForeground(new Color(0, 128, 0));
		loadGenes.setBackground(Color.WHITE);
		loadGenes.setFont(bfont);
		loadGenes.setBounds(6, 432, 100, 23);
		loadGenes.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				if(loadGenes.isSelected())
					{
					enrichment.setEnabled(true);
					searchTerms.setSelected(false);
					searchGenes.setSelected(false);
					btnLoad.setEnabled(true);
					}
				}
			});

		contentPanel.add(loadGenes);
		
		btnLoad = new JButton("Load");
		btnLoad.setEnabled(false);
		btnLoad.setBounds(323, 426, 80, 25);
		btnLoad.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0)
				{
				try{
					loadGenes();
				}catch(Exception e){e.printStackTrace();}
				}
		});
		contentPanel.add(btnLoad);
		
		JLabel lblOnList = new JLabel("genes on list...");
		lblOnList.setToolTipText("you can load a text file with gene names on different lines or on a single, tab separated line");
		lblOnList.setBounds(228, 432, 146, 16);
		lblOnList.setFont(font);
		contentPanel.add(lblOnList);
		
		JLabel lblTermsWith = new JLabel("terms with >");
		//lblTermsWith.setToolTipText("Difference between the average level along conditions in the first list respect to the second");
		lblTermsWith.setBounds(107, 419, 132, 16);
		lblTermsWith.setFont(font);
		contentPanel.add(lblTermsWith);
		
		JLabel lblTermsWithPvalues = new JLabel("terms with p-values <");
		//lblTermsWithPvalues.setToolTipText("Difference between the average level along conditions in the first list respect to the second");
		lblTermsWithPvalues.setBounds(129, 460, 165, 16);
		lblTermsWithPvalues.setFont(font);
		contentPanel.add(lblTermsWithPvalues);
		
		/*JLabel label_1 = new JLabel("terms with >");
		label_1.setToolTipText("Difference between the average level along conditions in the first list respect to the second");
		label_1.setBounds(107, 435, 132, 16);
		label_1.setFont(font);
		contentPanel.add(label_1);*/
		
		JLabel lblThreshold_1 = new JLabel("threshold");
		lblThreshold_1.setFont(new Font("Arial", Font.ITALIC, 13));
		lblThreshold_1.setToolTipText("Difference between the average level along conditions in the first list respect to the second");
		lblThreshold_1.setBounds(113, 338, 165, 16);
		contentPanel.add(lblThreshold_1);
		
		JLabel lblThreshold_2 = new JLabel("threshold");
		lblThreshold_2.setFont(new Font("Arial", Font.ITALIC, 13));
		lblThreshold_2.setBounds(320, 404, 91, 16);
		contentPanel.add(lblThreshold_2);
		
		JLabel label = new JLabel("threshold");
		label.setFont(new Font("Arial", Font.ITALIC, 13));
		label.setBounds(320, 366, 91, 16);
		contentPanel.add(label);
		
		JLabel lblTermsWithAvg = new JLabel("terms with avg expression above");
		lblTermsWithAvg.setBounds(107, 366, 227, 16);
		lblTermsWithAvg.setFont(font);
		contentPanel.add(lblTermsWithAvg);

		JPanel buttonPane = new JPanel();
		buttonPane.setBackground(Color.WHITE);
		
		getContentPane().add(buttonPane, BorderLayout.SOUTH);
		
		JButton okButton = new JButton("OK");
		okButton.setActionCommand("OK");
		okButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				deSearch();
				setVisible(false);
				}
			});
		
		JButton cancelButton = new JButton("Cancel");
		cancelButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				setVisible(false);
				}
			});
		buttonPane.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 5));
		cancelButton.setActionCommand("Cancel");
		buttonPane.add(cancelButton);
		
		buttonPane.add(okButton);
		getRootPane().setDefaultButton(okButton);
	}
	
	private void save()
		{
		JFileChooser selecFile = new JFileChooser();
		selecFile.addChoosableFileFilter(new TextFileFilter());
		
		if(enrichment.isSelected())//p-vals
			{	
			if(v.expData!=null)	
				{
				selecFile.setCurrentDirectory(new File(v.expData.filePath));
				selecFile.setSelectedFile(new File("enrichedTerms-"+v.expData.organismKegg+".txt"));
				}
			else
				selecFile.setSelectedFile(new File("enrichedTerms.txt"));
			}
		else
			{
			if(searchGenes.isSelected())	//degs
				{
				if(v.expData!=null)	
					{
					selecFile.setCurrentDirectory(new File(v.expData.filePath));
					selecFile.setSelectedFile(new File("degs-"+v.expData.organismKegg+".txt"));
					}
				else
					selecFile.setSelectedFile(new File("degs.txt"));
				}
			else			//avg expression
				{
				if(v.expData!=null)	
					{
					selecFile.setCurrentDirectory(new File(v.expData.filePath));
					selecFile.setSelectedFile(new File("terms-"+v.expData.organismKegg+".txt"));
					}
				else
					selecFile.setSelectedFile(new File("terms.txt"));
				}
			}
		
		int returnval = selecFile.showSaveDialog(this);

		if(returnval==JFileChooser.APPROVE_OPTION)
			saveFile=selecFile.getSelectedFile().getAbsolutePath();
		else
			saveFile=null;
		}
	
	private void loadGenes() throws Exception
		{
		JFileChooser selecFile = new JFileChooser();
		selecFile.addChoosableFileFilter(new TextFileFilter());
		
		
		if(v.voronto.degPath!=null)
			{
			BufferedReader br=new BufferedReader(new FileReader(v.voronto.degPath));
			String path=br.readLine();
			br.close();
			selecFile.setCurrentDirectory(new File(path));
			}
		else
			{
			if(v.expData!=null)
				{
				selecFile.setCurrentDirectory(new File(v.expData.filePath	));
				selecFile.setSelectedFile(new File("terms-"+v.expData.organismKegg+".txt"));
				}
			else
				selecFile.setSelectedFile(new File("terms.txt"));
			}
			
		int returnval = selecFile.showOpenDialog(this);
		//geneFile=null;
		if(returnval==JFileChooser.APPROVE_OPTION)
			{
			geneFile=selecFile.getSelectedFile().getAbsolutePath();
			BufferedWriter bw=new BufferedWriter(new FileWriter(v.voronto.degPath));
			bw.write(geneFile);
			bw.close();
			}
		}
	
	private void deSearch()
		{
		if(searchTerms.isSelected())			//1) search difexp terms
			v.deTermSearch(group1.getSelectedIndices(), group2.getSelectedIndices(), comboBox.getSelectedIndex(), new Float(textField.getText()).floatValue(), new Integer(textField_1.getText()).intValue(), saveFile);
		else if(this.searchGenes.isSelected())	//2) search difexp genes
			{
			if(enrichment.isSelected())	//+DAVID enrichment
				{
				//v.fisherEnrichment();
				v.deEnrichment(group1.getSelectedIndices(), group2.getSelectedIndices(), comboBox.getSelectedIndex(), new Float(textField.getText()).floatValue(), new Integer(textField_1.getText()).intValue(), new Float(textField_2.getText()).floatValue(), saveFile);
				}
			else
				v.deGeneSearch(group1.getSelectedIndices(), group2.getSelectedIndices(), comboBox.getSelectedIndex(), new Float(textField.getText()).floatValue(), new Integer(textField_1.getText()).intValue(), saveFile);
			}
		else									//3) load difexp genes
			{
			HashSet<String> genes=GeneParser.parse(geneFile);
			if(genes==null)
				JOptionPane.showMessageDialog(null, "No gene file selected", "No file selected", JOptionPane.ERROR_MESSAGE);
			
			int numMatches=v.expData.getGeneIds(genes.toArray(new String[0])).size();
			System.out.println(numMatches+" genes on the list ("+numMatches*1.0/genes.size()+"%) are present on the data matrix");
			HashSet<String> newGenes=null;
			if(numMatches<genes.size()*0.5)
				{
				//Should convert them to the right gene ids in case they come from other sources that the matrix's gene ids
				HashSet<String> entrezs=v.expData.getIdsFromSynonyms(genes,ExpressionData.ENTREZ);
				if(entrezs.size()>numMatches)
					{
					numMatches=entrezs.size();
					newGenes=entrezs;
					}
				HashSet<String> ensembls=v.expData.getIdsFromSynonyms(genes,ExpressionData.ENSEMBL);
				if(ensembls.size()>numMatches)
					{
					numMatches=ensembls.size();
					newGenes=ensembls;
					}
				HashSet<String> symbols=v.expData.getIdsFromSynonyms(genes,ExpressionData.SYMBOL);
				if(symbols.size()>numMatches)
					{
					numMatches=entrezs.size();
					newGenes=symbols;
					}
				}
			if(newGenes!=null)
				{
				genes=newGenes;
				System.out.println(numMatches+" found using synonyms ("+numMatches*1.0/genes.size()+"%)");
				}
			else
				{
				newGenes=new HashSet<String>();
				for(String g:genes)	newGenes.add(g.toLowerCase());
				genes=newGenes;
				}
			
			if(enrichment.isSelected())	//+DAVID enrichment
				{
				v.deEnrichment(genes, new Integer(textField_1.getText()).intValue(), new Float(textField_2.getText()).floatValue(), saveFile);
				}
			else
				v.giSearch(genes, new Integer(textField_1.getText()).intValue(), saveFile);
			}
		saveFile=null;
		chckbxSaveToFile.setSelected(false);
		btnSelect.setEnabled(false);
		}
}
