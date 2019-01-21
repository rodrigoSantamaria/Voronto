package es.usal.voronto.model.ontology;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import es.usal.voronto.model.expression.ExpressionData;
import es.usal.voronto.view.VoronoiVisualization;


public class GOparser {
	
	public static void main(String[] args)
		{
		/*File dir=new File("/Users/rodri/Desktop/vorontoParsing/goa/goa");
		for(File f:dir.listFiles())
			{
			String s=f.getAbsolutePath();
			if(s.contains("gene_association.goa"))
				GOparser.map(s, "/Users/rodri/Desktop/vorontoParsing/goa/new/go-"+convert2(s.substring(s.indexOf(".goa")+5))+"_gene_ensembl.map");
			}*/
		GOparser.map("/Users/rodri/Desktop/vorontoParsing/goa/gene_association.GeneDB_Spombe.genedb_spombe", "/Users/rodri/Documents/investigacion/distros/git/voronto/Voronto/src/es/usal/voronto/data/go/go-spombe_gene_ensembl.map");
		//GOparser.map("/Users/rodri/Desktop/goa/gene_association.cgd", "/Users/rodri/Documents/investigacion/distros/git/voronto/Voronto/src/es/usal/voronto/data/go/go-calbicans_gene_ensembl.map");
		//GOparser.map("/Users/rodri/Documents/investigacion/distros/git/voronto/Voronto/data/calbicans/annotations/gene_association.cgd", "/Users/rodri/Documents/investigacion/distros/git/voronto/Voronto/src/es/usal/voronto/data/go/gobiological_process-calbicans_gene_ensembl.map");
		
		//GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", "biological_process", true);
		//GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", "cellular_component", true);
		//GOparser.parse("es/usal/voronto/data/go/gene_ontology_ext.obo", "molecular_function", true);
    	}
	
	private static List<OntologyTerm> lot;

	/**
	 * Unlike KO parser, this method only builds the hierarchy, but does not specify the gene mapping
	 * Then, method annotate should be used for gene mapping, providing the species name
	 * @param path
	 * @param ontology - either biological_process, cellular_component or molecular_function
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> parse(String path, String ontology, boolean redo)
		{
		TreeMap<OntologyTerm, TreeMap> map=new TreeMap<OntologyTerm, TreeMap>();
		String serName=null;
	//	if(path.contains("slim"))	serName="goSLIM";
	//	else						serName="go";
		
	//	String goTermName="terms";
	//	if(path.contains("slim"))	goTermName="goSLIM"+ontology+".txt";
	//	else						goTermName="go"+ontology+".txt";
		
		if(!redo)
			try
			   	{
				System.out.println("Reading GO serialized map");
				long time=System.currentTimeMillis();
				InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/go/go-"+ontology+".ser");
				ObjectInputStream oin = new ObjectInputStream(is);
				map = (TreeMap<OntologyTerm, TreeMap>)oin.readObject();
			    oin.close();
			    System.out.println("GO hierarchy read in "+(System.currentTimeMillis()-time)/1000.0);
			    }catch(Exception e){e.printStackTrace(); redo=true;}
		
		if(map.size()==0)
			{
			//NOTE: IMPORTANT: remove gene_ontology_ext.obo from the project before export, this is only done if redo is set or there's no serialized map for this species
			map=OBOparser.parse(path, ontology, ontology);//root names equals namespace name
			}
		
		if(redo)
			{
			try{
				//B3) Serialize the mapping for optimization (still it takes time to read after, better divide by species
				//FileOutputStream fos = new FileOutputStream(serName+ontology+".ser");
				FileOutputStream fos = new FileOutputStream("go-"+ontology+".ser");
				ObjectOutputStream out = new ObjectOutputStream(fos);
				out.writeObject(map);
				out.close();
				
				//B4) Write Go terms to a file as a list:
				/*BufferedWriter fichero=new BufferedWriter(new FileWriter(goTermName));
				for(String go:goTerms)	fichero.write(go+"\n");
				fichero.close();*/
				}catch(Exception e){e.printStackTrace();}
			}
		
		TreeMap<OntologyTerm, TreeMap> onto=new TreeMap<OntologyTerm, TreeMap>();
		onto.put(new OntologyTerm("GO slim ontology", ""), map);
		
		return onto;
		}
	
	/**
	 * Takes a GAF 2.0 file and generates a map from it (basically the same file, but smaller: for Human from 57M to 5.8M)
	 * GAF interesting columns:
	 * 		http://www.geneontology.org/GO.format.gaf-2_0.shtml
	 * 	3 - DB Object symbol - a gene symbol (external_gene_id can be mapped to)
	 *  5 - Term ID (mandatory as the former one)
	 *  10 - Gene name
	 *  11 - Synonyms (non mandatory, as the former one, but strongly recommended
	 * @param gafFilePath
	 */
	public static void map(String gafFilePath, String mapFilePath)
		{
		BufferedReader in;
		BufferedWriter out;
		try {
			in = new BufferedReader(new FileReader(gafFilePath));
			out = new BufferedWriter(new FileWriter(mapFilePath));
		//	out.write("go_id\texternal_gene_id\tevidence_id");	out.newLine();
			out.write("go_id\texternal_gene_id");	out.newLine();
			String cad=null;
			HashMap<String, ArrayList<String>> set=new HashMap<String,ArrayList<String>>();
			while((cad=in.readLine())!=null)
				{	
				if(!cad.startsWith("!"))
					{
					String[] fields=cad.split("\t");
					String go=fields[4];
					ArrayList<String> orfs=new ArrayList<String>();
					orfs.add(fields[2]);
					
					//This could be good, but makes the files too large for distribution: better to use separate synonym lists
					/*if(fields[10]!=null && fields.length>0)
						{
						String[] syn=fields[10].split("\\|");//This can be non existing cause it is not mandatory, should include fields[2]
						if(syn!=null && syn.length>0)
							for(String orf:syn)
								{
								if(orf.indexOf("_")>0)	orf=orf.substring(0, orf.indexOf("_"));
								if(!orfs.contains(orf))	orfs.add(orf);
								}
						}*///too large, synonyms via the synonym data
					
					
					//String evidence=fields[6];	//Not using evidences by now
					if(set.get(go)==null)	set.put(go, new ArrayList<String>());
					for(String orf:orfs)
					//	if(!evidence.equals("IEA") && !evidence.equals("NR"))
							//{out.write(go+"\t"+orf+"\t"+evidence); out.newLine();}
							{
							if(!set.get(go).contains(orf))
								{
								out.write(go+"\t"+orf); out.newLine();
								set.get(go).add(orf);
								}
							}
					
						//String orf=fields[10];//This can be non existing cause it is not mandatory, should include fields[2]
						//if(orf.indexOf("|")>0)	orf=orf.substring(0, orf.indexOf("|"));	//these are several synonyms... might be interesting to get all of them
					
					//String orf=fields[2];
					//out.write(go+"\t"+orf); out.newLine();
					}
				}
			out.close();
		}catch(Exception e){e.printStackTrace();}
		}
	
	public static void map2(String gafFilePath, String mapFilePath)
		{
		BufferedReader in;
		BufferedWriter out;
		try {
			in = new BufferedReader(new FileReader(gafFilePath));
			out = new BufferedWriter(new FileWriter(mapFilePath));
			out.write("go_id\tgenes");	out.newLine();
			String cad=null;
			while((cad=in.readLine())!=null)
				{	
				if(!cad.startsWith("!"))
					{
					String[] fields=cad.split("\t");
					String go=fields[4];
					//String orf=fields[10];//This can be non existing cause it is not mandatory, should include fields[2]
					//if(orf.indexOf("|")>0)	orf=orf.substring(0, orf.indexOf("|"));	//these are several synonyms... might be interesting to get all of them
					String orf=fields[2];
					out.write(go+"\t"+orf); out.newLine();
					}
				}
			out.close();
		}catch(Exception e){e.printStackTrace();}
		}
	
	
	/**
	 * Returns the map m, but with its ontology terms mapped to gene ids on the corresponding format (entrezgene, external_gene_id or ensembl_gene_id)
	 * This method reads annotation files produced by biomaRt via R (see script es.usal.voronto.rcode.retrieveGOannotations.R)
	 * These files are pregenerated and updated periodically (set and automate period... six months?)
	 * @param m			GO hierarchy to annotate
	 * @param species	full standard species name (e.g. Homo sapiens, but not H sapiens, hsa or human)
	 * @param gene_id	either entrezgene, external_gene_id, ensembl_gene_id or (standard id names from biomaRt)
	 * TODO: either generate a file for each GO type (Slim or not; BP/MF/CC) or big annotation files with everything and select depending if they appear in m
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> annotate(TreeMap<OntologyTerm, TreeMap> m, String species, String gene_id, int type, ExpressionData ed) throws Exception 
		{
		String ontology="";
		switch(type)
			{
			case VoronoiVisualization.SLIMBP:
				ontology="goSLIMbiological_process";
				break;
			case VoronoiVisualization.SLIMCC:
				ontology="goSLIMcellular_component";
				break;
			case VoronoiVisualization.BP://TODO: change name to generic go_
			case VoronoiVisualization.CC:
			case VoronoiVisualization.MF:
					//ontology="gobiological_process";
					ontology="go";
				break;
			}
		
		BufferedReader in = null; 
		
		try{
		String path="es/usal/voronto/data/go/"+ontology+"-"+convert(species)+"_gene_ensembl.map";//by now testing with GO slim bp
		in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));//para applet/jws
		}catch(Exception e)
			{
			throw new Exception("No GO mapping for organism "+species);
			}
		
		String header=in.readLine();
		String []fields=header.split("\t");
		int id=-1;//column with the requested gene ids.
		for(int i=0;i<fields.length;i++)
			if(fields[i].equals(gene_id))	{id=i; break;}
		
		String[] compareTo=ed.sortedGeneNames;
		HashMap<String,String> compareHash=null;
		//in the case the id in the expression data does not match with any of the ids on the GO mapping, we check if we have synonyms recorded for the expr. data
		if(id==-1)
			{
			if(gene_id.equals("entrezgene") && ed.entrezgene!=null && ed.entrezgene.length>0)			
				{
				compareTo=ed.entrezgene;
				for(int i=0;i<fields.length;i++)
					if(fields[i].equals("ensembl_gene_id")  || fields[i].equals("external_gene_id"))	
						{
						id=i;
						if(fields[i].equals("ensembl_gene_id"))
							compareHash=ed.invertedHash(ed.ensembl_gene_idHash);
						if(fields[i].equals("external_gene_id"))
							compareHash=ed.invertedHash(ed.external_gene_idHash);
						break;
						}
				}
			if(gene_id.equals("ensembl_gene_id") && ed.ensembl_gene_id!=null && ed.ensembl_gene_id.length>0)			
				{
				compareTo=ed.ensembl_gene_id;
				for(int i=0;i<fields.length;i++)
					if(fields[i].equals("entrezgene") || fields[i].equals("external_gene_id"))	
						{
						id=i;
						if(fields[i].equals("entrezgene"))
							compareHash=ed.invertedHash(ed.entrezgeneHash);
						if(fields[i].equals("external_gene_id"))
							compareHash=ed.invertedHash(ed.external_gene_idHash);
						break;
						}
				}
			if(gene_id.equals("external_gene_id") && ed.external_gene_id!=null && ed.external_gene_id.length>0)			
				{
				compareTo=ed.external_gene_id;
				compareHash=ed.external_gene_idHash;
				for(int i=0;i<fields.length;i++)
					if(fields[i].equals("entrezgene") || fields[i].equals("ensembl_gene_id"))
						{
						id=i;
						if(fields[i].equals("entrezgene"))
							compareHash=ed.invertedHash(ed.entrezgeneHash);
						if(fields[i].equals("ensembl_gene_id"))
							compareHash=ed.invertedHash(ed.ensembl_gene_idHash);
						break;
						}
				}
			
			}
		if(id==-1)
			{
			String s="Gene id '"+gene_id+"' is not supported for this GO map, please use one of these:\n";
			for(int i=1;i<fields.length;i++)	s+="\t"+fields[i];
			
			throw new Exception(s);
			}
		
		String cad;
		TreeMap<String, ArrayList<String>> annotations=new TreeMap<String, ArrayList<String>>();
		TreeSet<String> mappedGenes=new TreeSet<String>();
		while((cad=in.readLine())!=null)
			{
			fields=cad.split("\t");
			if(fields!=null && fields.length>id && fields[id].length()>0 && !fields[id].equals("NA"))	//add the annotation to a data structure
				{
				if(annotations.get(fields[0])==null)
					annotations.put(fields[0], new ArrayList<String>());
				if(compareHash!=null)
					{
					String gene=compareHash.get(fields[id].toLowerCase());
					if(gene!=null)
						{
						annotations.get(fields[0]).add(gene);
						//NOTE: provisional only for NDR data:
						annotations.get(fields[0]).add(gene+"(1)");
						annotations.get(fields[0]).add(gene+"(2)");
						annotations.get(fields[0]).add(gene+"(3)");
						annotations.get(fields[0]).add(gene+"(4)");
						annotations.get(fields[0]).add(gene+"(5)");
						}
					}
				else	
					{
					annotations.get(fields[0]).add(fields[id].toLowerCase());
					//NOTE: provisional only for NDR data:
					annotations.get(fields[0]).add(fields[id].toLowerCase()+"(1)");
					annotations.get(fields[0]).add(fields[id].toLowerCase()+"(2)");
					annotations.get(fields[0]).add(fields[id].toLowerCase()+"(3)");
					annotations.get(fields[0]).add(fields[id].toLowerCase()+"(4)");
					annotations.get(fields[0]).add(fields[id].toLowerCase()+"(5)");
					}
				mappedGenes.add(fields[id].toLowerCase());
				}
			}
		
		
		
		//---
		if(ed!=null)
			{
			int map=0;
			
			HashMap<String,String> compareHashI=null;
			if(compareHash!=null)	compareHashI=ed.invertedHash(compareHash);
			for(String mg:compareTo)
				{
				if(compareHashI!=null)	
					{
					mg=compareHashI.get(mg.toLowerCase());
					//mg=compareHashI.get(mg.toLowerCase().replaceAll("\\(.*\\)", ""));
					}
				if(mg!=null)
					if(mappedGenes.contains(mg.toLowerCase()))
					//if(mappedGenes.contains(mg.toLowerCase().replaceAll("\\(.*\\)", "")))
						map++;
				}
				
			System.out.println(map+"/"+ed.geneNames.length+" genes annotated on one or more terms ("+map*1.0/ed.geneNames.length+"%)");
			}
		//---
		
		
		long time=System.currentTimeMillis();
		addRecursiveAnnotations(annotations,m,null,0);
		
		System.out.println("time in adding annotations "+(System.currentTimeMillis()-time)/1000.0);
		
		return m;
		}
	
	/**
	 * 
	 * @param m
	 * @param path	GAF path
	 * @param gene_id
	 * @param type
	 * @param ed
	 * @return
	 * @throws Exception
	 */
	public static TreeMap<OntologyTerm, TreeMap> annotateFromGAF(TreeMap<OntologyTerm, TreeMap> m, String path, ExpressionData ed) throws Exception 
	{
	BufferedReader in = null; 

	try{
		in = new BufferedReader(new FileReader(path));
		}catch(Exception e)		{	throw new Exception("File "+path+" not found");		}
	
	String header=in.readLine();
	String []fields=header.split("\t");
	
	int id=1;//column of the gene id
	int symbolColumn=2;//column of the gene symbol
	int termColumn=4;
	int synoColumn=10;
	
	int idColumn=id;
	if(ed.chip.equals("GAFsymbol"))
		{
		idColumn=symbolColumn;
		symbolColumn=id;
		}
	else if(ed.chip.equals("GAFsynonym"))
		{
		idColumn=synoColumn;
		}
		
	String cad;
	TreeMap<String, ArrayList<String>> annotations=new TreeMap<String, ArrayList<String>>();
	TreeSet<String> mappedGenes=new TreeSet<String>();
	while((cad=in.readLine())!=null)
		{
		if(!cad.startsWith("!"))
			{	
			fields=cad.split("\t");
			String fid=fields[idColumn].split("\\|")[0].trim();
			if(fields!=null && fid.length()>0 && !fid.equals("NA"))	//add the annotation to a data structure
				{
				if(annotations.get(fields[termColumn])==null)
					annotations.put(fields[termColumn], new ArrayList<String>());
				annotations.get(fields[termColumn]).add(fid.toLowerCase());
				System.out.println("Adding id: "+fid.toLowerCase()+" to "+fields[termColumn]);
				mappedGenes.add(fid.toLowerCase());
				
				String symbol=fields[symbolColumn].trim();	//TODO: these "annotation by symbol/synonym" will add new ids, they should just attach synonyms/symbols to previous ids
				if(symbol!=null && symbol.length()>0 && !symbol.equals("NA"))
					{
					symbol=symbol.toLowerCase();
					if(ed.external_gene_idHash==null)	ed.external_gene_idHash=new HashMap<String,String>();
					ed.external_gene_idHash.put(fields[id].toLowerCase(), symbol);
					if(!annotations.get(fields[termColumn]).contains(symbol))
						{
						annotations.get(fields[termColumn]).add(symbol);
						System.out.println("Adding symbol: "+symbol.toLowerCase()+" to "+fields[termColumn]);
						mappedGenes.add(symbol.toLowerCase());
						}
					}
				if(fields.length>=10)	
					{
					String[] synonyms=fields[synoColumn].split("\\|");
					for(String synonym:synonyms)	
						{
						synonym=synonym.toLowerCase().trim();
						if(!annotations.get(fields[termColumn]).contains(synonym))
							{
							annotations.get(fields[termColumn]).add(synonym.toLowerCase());
							System.out.println("Adding synonym: "+synonym.toLowerCase()+" to "+fields[termColumn]);
							mappedGenes.add(synonym.toLowerCase());
							}
						}
					}
				}
			}
		}
	
	System.out.println("Genes annotated vs total genes in expression data: "+mappedGenes.size()+"/"+ed.sortedGeneNames.length+" ("+mappedGenes.size()/ed.sortedGeneNames.length+"%)");
	//---
	if(ed!=null)
		{
		int nomap=0;
		
		for(String mg:ed.sortedGeneNames)
			if(!mappedGenes.contains(mg))
				{
				nomap++;
				}
		System.out.println("Genes in the expression data but not annotated: "+nomap);
		}
	//---
	
	long time=System.currentTimeMillis();
	addRecursiveAnnotations(annotations,m,null,0);
	in.close();
	System.out.println("time in adding annotations "+(System.currentTimeMillis()-time)/1000.0);
	
	return m;
	}

	/**
	 * Adds to each OntologyTerm in m (recursively) genes related to it, following annotations
	 * NOTE: if an OntologyTerm is annotated with some genes but its parent is not, the parent appears but without annotated genes
	 * @param annotations
	 * @param m
	 * @param level
	 */
	public static void addRecursiveAnnotations(TreeMap<String, ArrayList<String>> annotations, TreeMap<OntologyTerm, TreeMap> m, OntologyTerm parent, int level)
		{
		for(OntologyTerm ot:m.keySet())
			{
			ArrayList<String> annot=annotations.get(ot.id);
			if(annot!=null && ot.geneIds.size()==0)	
				ot.geneIds.addAll(annot);	//add annotations only once
			
			TreeMap<OntologyTerm,TreeMap> mc=m.get(ot);			//and keep recursing
			if(mc!=null)	
				addRecursiveAnnotations(annotations, mc, ot, level+1);
			
			if(parent!=null)			//and check if the parent has it (because some ontologies are inconsistent!)	NOTE: superSLOW in GO-full+mice (70s)
				parent.geneIds.addAll(ot.geneIds);
			}
		return;
		}
	
	/**
	 * Converts standard species names to OBO species names
	 * E.g. Homo sapiens --> hsapiens
	 * @param species
	 * @return
	 */
	public static String convert(String species)
	{
	return species.trim().toLowerCase().charAt(0)+species.substring(species.indexOf(" ")+1);
	}
	public static String convert2(String species)
		{
		if(species.equals("arabidopsis"))	return "athaliana";
		if(species.equals("chicken"))	return "ggallus";
		if(species.equals("cow"))	return "btaurus";
		if(species.equals("dicty"))	return "ddiscoideum";
		if(species.equals("dog"))	return "cfamiliaris";
		if(species.equals("fly"))	return "dmelanogaster";
		if(species.equals("human"))	return "hsapiens";
		if(species.equals("mouse"))	return "mmusculus";
		if(species.equals("pig"))	return "sscrofa";
		if(species.equals("rat"))	return "rnorvegicus";
		if(species.equals("worm"))	return "celegans";
		if(species.equals("yeast"))	return "scerevisiae";
		if(species.equals("zebrafish"))	return "drerio";
			
		return species;
		}

	
	/**
	 * RESTful access, get annotations for a given term
	 * http://amigo.geneontology.org/cgi-bin/amigo/term-assoc.cgi?term=GO:0007200&term_assocs=direct&format=go_assoc
	 * term=GO:007200
	 * term_assocs=direct
	 * format=go_assoc
	 * 
	 * Get annotations for a given gene (GB: for entrezgene, ENSEMBL for ensembl)
	 * http://amigo.geneontology.org/cgi-bin/amigo/gp-assoc.cgi?gp=FB:FBgn0015946&format=go_assoc
	 * 
	 * So, por GO BP (full) we do have 22000 terms.
	 * We can
	 * 1) REST Search for each term its annotated genes --> 22000 queries
	 * 2) REST Search for each gene its annotated terms --> as much queries as genes in our matrix (similar to above)
	 * 3) Offline BiomaRt search for each term --> one for each organism, store large files
	 * 
	 * Pros and cons
	 * 		Memory
	 * 1	No memory increase
	 * 2	No memory increase
	 * 3	Large files storage: increase the size of the program up to 100M, even if only some organisms (human, mouse, yeast) recorded
	 * 
	 * 		Time
	 * 1	22000xt, being t the time for a REST query
	 * 2	Nxt, being N the number of genes in the expression data. Expected to be similar that in 1
	 * 3	Offline time, no charge on client
	 * 
	 * 		Update
	 * 1	Always up to date
	 * 2	Always up to date
	 * 3	Need of periodical updates
	 * 
	 * 		Memory	Time	Update
	 * 1	Good	Bad		Good
	 * 2	Good	Bad		Good
	 * 3	Bad		Good	Bad
	 * 
	 * But time is much more important in a visualizaiton, and is already burdened by tessellation.
	 * A solution like 3 is the one that has been selected for KEGG and GO slim.
	 * 
	 * Let's quantize how much it takes each thing:
	 * Time (term annotation retrieval)
	 * 1
	 * 2
	 * 3	3h15' for all the BP terms for mus musculus (4.6M) --> which means it's unafordable for running time at all
	 * 
	 *  Time (term annotation file parsing and mapping)
	 *  1
	 *  2
	 *  3
	 */
}
