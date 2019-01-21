package es.usal.voronto.model.ontology;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.BiochemicalReaction;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.Entity;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Protein;
import org.biopax.paxtools.model.level3.Xref;
import org.biopax.paxtools.impl.level3.PathwayImpl;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.io.jena.JenaIOHandler;
import org.biopax.paxtools.model.level3.Process;


import java.util.Arrays;
import es.usal.voronto.model.expression.ExpressionData;

public class ReactomeParser 
	{
	//public SimpleIOHandler handler;	//Simple handle cannot work with large files, apparently
	public JenaIOHandler handler;

	public static void main(String args[])
		{
		long t=System.currentTimeMillis();
		//A) Parse a single file
		//ReactomeParser.parse("/Users/rodri/Desktop/biopax3/Homo sapiens.owl");
		//B) Parse every file in a folder
		ReactomeParser.buildReactomeHierarchy("/Users/rodri/Desktop/biopax3/reactomeSep2012/mapped");
		System.out.println("parsing takes "+(System.currentTimeMillis()-t)/1000.0);
		//ReactomeParser.parse("/Users/rodri/Desktop/biopax3/mapped/Caenorhabditis elegans.owl");
		}
	
	/**
	 * Generates serialized maps for every OWL file in the folder path.
	 * Each serialized map is a TreeMap<OntologyTerm, TreeMap>, where OntologyTerm in this case is a 
	 * ReactomePathway and its TreeeMap contains every ReactomePathway in it
	 * NOTE: it takes about 3/4 of an hour to parse and serialize every OWL file in Reactome
	 * NOTE: fails for Homo sapiens right now
	 * @param folderPath
	 */
	public static void buildReactomeHierarchy(String folderPath)
		{
		File dir=new File(folderPath);
		for(String file:dir.list())
			{
			if(file.endsWith("owl"))	
				{
				System.out.println(file);
				ReactomeParser.parse(folderPath+"/"+file);
				}
			}
		}
	
	/**
	 * Reads and return the TreeMap in a serialized file with BioPAXparser.parse
	 * @param file Serialized file with Reactome hierarchy and annotations for a given species
	 * @param ed Expression data. Required if expression data ids are not entrezgenes, to make
	 * 				a proper synonym mapping
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> readSer(String file, ExpressionData ed) throws Exception
		{
		TreeMap<OntologyTerm, TreeMap> map=null;
		try
		   {
			System.out.println("Reading "+file);
			InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream(file);
			ObjectInputStream oin = new ObjectInputStream(is);
	        map= (TreeMap<OntologyTerm, TreeMap>)oin.readObject();
		    oin.close();
		    System.out.println("reactome hierarchy read");
		   }catch(Exception e){throw new Exception("No map file for this organism in Reactome");}
		    
		    //We now remove root elements with no children (they arise from inconsistencies on the OWL format)
		    OntologyTerm[] set=map.keySet().toArray(new OntologyTerm[0]);
		    for(OntologyTerm ot:set)
		    	{
		    	if(map.get(ot)==null || map.get(ot).isEmpty())
		    		map.remove(ot);
		    	}
		   // if(ed!=null && !ed.chip.equals("entrezgene") && ed.entrezgeneHash!=null)
		    if(ed!=null)
				{
		    	if(!ed.chip.equals("ensembl_gene_id") && ed.ensembl_gene_idHash!=null)
		    		//proceed with the translation of annotations
		    		//recursiveTranslation(map, ed.invertedHash(ed.entrezgeneHash));
		    		recursiveTranslation(map, ed.invertedHash(ed.ensembl_gene_idHash));
		    	else
		    		recursiveTranslation(map, null);
		    	}
		
		return map;	
		}
	
	private static void recursiveTranslation(TreeMap<OntologyTerm, TreeMap> m, HashMap<String,String> geneMap)
		{
		if(m==null)	return;
		for(OntologyTerm ot:m.keySet())
			{
			HashSet<String> annot=ot.geneIds;
			HashSet<String> annot2=new HashSet();
			System.out.println(ot.id+" has "+annot.size());
			for(String s:annot)
				{
				String st=null;
				if(geneMap!=null)	st=geneMap.get(s.toLowerCase());
				else				st=s;
				if(st!=null)	
					annot2.add(st.toLowerCase());
				}
			ot.geneIds.clear();
			ot.geneIds.addAll(annot2);
			recursiveTranslation(m.get(ot), geneMap);
			}
		return;
		}
	/**
	 * Reads OWL file, saves the hierarchy on a TreeMap, then annotates the hierarchy with the information in a map file with the same name,
	 * finally returns and saves the TreeMap in a serialized file with the same name.
	 * @param file
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> parse(String file)
		{
		TreeMap<OntologyTerm, TreeMap> map=new TreeMap<OntologyTerm, TreeMap>();
		
		ReactomeParser parser=new ReactomeParser();
		parser.handler=new JenaIOHandler(BioPAXLevel.L3);
		try{
			//1) Get the model
			long t=System.currentTimeMillis();
			Model model=parser.handler.convertFromOWL(new FileInputStream(file)); 
			System.out.println("file parsed in "+(System.currentTimeMillis()-t)/1000.0);
			
			//2) Parse pathways
			ArrayList<OntologyTerm> roots=new ArrayList<OntologyTerm>();//root terms
		
			//children paths for each term
			HashMap<String, ArrayList<OntologyTerm>> paths=new HashMap<String, ArrayList<OntologyTerm>>();
			OntologyTerm otem;
			for(BioPAXElement element:model.getObjects().toArray(new BioPAXElement[0]))//for each element
				{
				if(element.getClass().toString().contains("PathwayImpl"))//if it is a pathway
					{
					PathwayImpl pathway=(PathwayImpl)element;
					
					System.out.println("Pathway "+pathway.getDisplayName());
					
					String id=null;
					Iterator<Xref> itr=pathway.getXref().iterator();
					while(itr.hasNext())
						{
						Xref xref=itr.next();
						if(xref!=null && xref.getDb()!=null && xref.getDb().equals("Reactome Database ID"))
							{
							id=xref.getId();
							break;
							}
						}
					
					if(id==null)
						System.out.println("Error: pathway "+pathway.getName()+" has no reactome id!!");
					
					OntologyTerm ot=new OntologyTerm(pathway.getDisplayName(), id);
					
					System.out.println("\tWith "+pathway.getPathwayComponentOf().size()+" parents");//parents
					if(pathway.getPathwayComponentOf().size()==0)
						roots.add(ot);
					
					System.out.println("\tWith "+pathway.getPathwayComponent().size()+" children");//children
					Iterator<Process> it=pathway.getPathwayComponent().iterator();
					ArrayList<OntologyTerm> children=new ArrayList<OntologyTerm>();
					while(it.hasNext())
						{
						Process p=(Process)it.next();
						if(p.getClass().toString().contains("Pathway"))
							{
							PathwayImpl path=(PathwayImpl)p;
							String idchild=null;
							itr=path.getXref().iterator();
							while(itr.hasNext())
								{
								Xref xref=itr.next();
								if(xref!=null && xref.getDb()!=null && xref.getDb().equals("Reactome Database ID"))
									{
									idchild=xref.getId();
									break;
									}
								}
							if(id==null)
								System.err.println("ID of children "+path.getDisplayName()+" not found");
							else	
								children.add(new OntologyTerm(path.getDisplayName(), idchild));
							}
						}
					paths.put(ot.id, children);
					}
				}
			
			//3) Build map recursively
			System.out.println("Building the map");
			for(OntologyTerm ot:roots)
				{
				TreeMap<OntologyTerm, TreeMap> temp=new TreeMap<OntologyTerm, TreeMap>();
				
				for(OntologyTerm ch:paths.get(ot.id))
					temp.put(ch, null);
				map.put(ot, temp);
				addChildrenRecursively(map.get(ot), paths);
				}
			printMap(map, 0);
			
			//4) Add annotations from additional file
			System.out.println("Annotation");
			BufferedReader br=new BufferedReader(new FileReader(file.substring(0, file.indexOf("."))+".map"));
			String id=null;
			while((id=br.readLine())!=null)
				{
				String[] genes=br.readLine().split(" ");
				ArrayList<String> list=new ArrayList<String>();
				list.addAll(Arrays.asList(genes));
				OntologyTerm ot=searchFor(map, id);
				if(ot!=null)	ot.geneIds.addAll(list);
				else			System.err.println("id "+id+" not found");
				}
			
			//5) Save serialized object
			FileOutputStream fos = new FileOutputStream(file.substring(0, file.indexOf("."))+".ser");
			ObjectOutputStream out = new ObjectOutputStream(fos);
			out.writeObject(map);
			out.close();	
			}catch(Exception e){e.printStackTrace();}
		return map;
		}

	/**
	 * Searches for an OntologyTerm with the given id within the hierarchy.
	 * @param map
	 * @param id
	 * @return
	 */
	private static OntologyTerm searchFor(TreeMap<OntologyTerm, TreeMap> map, String id)
		{
		if(map==null)	return null;
		for(OntologyTerm ot:map.keySet())
			{
			if(ot.id.equals(id))	return ot;
			else					
				{
				OntologyTerm ret=searchFor(map.get(ot), id);
				if(ret!=null)	return ret;
				}
			}
		return null;
		}
	
	private static void printMap(TreeMap<OntologyTerm, TreeMap> map, int level)
		{
		if(map==null)	return;
		for(OntologyTerm ot:map.keySet())
			{
			for(int i=0;i<level;i++)	System.out.print("\t");
			System.out.println(ot.name+" | "+ot.id);
			printMap(map.get(ot), level+1);
			}
		}
	
	private static void addChildrenRecursively(TreeMap<OntologyTerm, TreeMap> map, HashMap<String, ArrayList<OntologyTerm>> children)
		{
		for(OntologyTerm ot:map.keySet())
			{
			if(ot.name.equals("Carbohydrate metabolism"))
				System.out.println("meta");
			
			ArrayList<OntologyTerm> c=children.get(ot.id);
			if(c!=null && c.size()>0)
				{
				map.put(ot,new TreeMap<OntologyTerm, TreeMap>());
				for(OntologyTerm ot2:c)	map.get(ot).put(ot2,null);
				addChildrenRecursively(map.get(ot), children);
				}
			}
		}
	}
	