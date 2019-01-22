package es.usal.voronto.model.ontology;

import java.io.BufferedReader;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.HttpURLConnection;
import java.net.URI;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;


import javax.xml.rpc.ServiceException;
import es.usal.voronto.model.expression.ExpressionData;

/**
 * Parsers KEGG orthology file into a Map
 * @author rodri
 *
 */
public class KOparser {

	public static ArrayList<String> getPathways(String organism)
		{
		ArrayList<String> paths=new ArrayList<String>();
		
		try{
			URL url = new URL("http://rest.kegg.jp/list/pathway/"+organism); 
			HttpURLConnection conn = (HttpURLConnection) url.openConnection(); 
			conn.setRequestMethod("GET");
			BufferedReader br = new BufferedReader(new InputStreamReader((conn.getInputStream())));
			String output; 
			//System.out.println("Output from Server .... \n");
			while ((output = br.readLine()) != null) 
				{
				//System.out.println(output);
				String p=output.split("\t")[0].substring(8);
				while(!p.substring(0, 1).matches("[0-9]"))
						p=p.substring(1);
				paths.add(p);
				}
			 conn.disconnect();
		}
		catch(FileNotFoundException e){System.out.println("Kegg pathway list for "+organism+" not found"); return null;}
		catch(Exception e){e.printStackTrace(); return null;}
		
		return paths;
		}
	/**
	 * Parses KEGG orthology (format from dec 2018)
	 * 
	 * @param path KEGG hierarchy serialized file
	 * @param organism KEGG organism, so it selects only entries for that organism (e.g. if Saccharomyces cerevisiae --> sce)
	 * 			TODO: if null, try elements as KO ids directly
	 * @param pathIDs IDs of the paths that must be included in the ontology (if null, all of them will be included)
	 * @param ed  ExpressionData which will determine selected gene annotations and gene id (it can be 'entrezgene', 'ensembl_gene_id' or 'external_gene_id')
	 * 					By default is null, and gene annotations will be the ones coming directly from KEGG (usually entrezgene, but not always)
	 * @param update If true, it uses programmatic remaps elements and updates the serialized file (takes time)
	 * @param brite If true, only the brite hierarchy is parsed (A code 09180). If false, KEGG pathways are parsed (A codes A09100 to A09160)
	 * @return
	 */
	public static TreeMap<OntologyTerm, TreeMap> parse(String path, String organism, ArrayList<String> pathIDs, ExpressionData ed, boolean update, boolean brite)
	{
	TreeMap<OntologyTerm, TreeMap> map=new TreeMap<OntologyTerm, TreeMap>();
	TreeMap<String, ArrayList<String>> koterms=new TreeMap<String, ArrayList<String>>();//mapping of KO terms to gene ids
	
	
	//getGenesByKO("K01805");
	FileInputStream fis = null;
	ObjectInputStream oin = null;
	boolean addElements=true;
	HashMap<String, String> geneMap=null;
	
	try
	   {
		//if(!update)	//For an incremental update, we can comment this if and read the previous serialization
			{
			System.out.println("Reading komappingTOTAL");
			//InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/kegg/komappingTOTAL.ser");
			InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/kegg/komap-"+organism+".ser");
			oin = new ObjectInputStream(is);
	        koterms = (TreeMap<String, ArrayList<String>>)oin.readObject();
		    oin.close();
		    System.out.println("koterms read");
			}
	    if(ed!=null)
	    	{
	    	/*
	    	 * KEGG specials. Except otherwise, KEGG genes are ENTREZ for all the organisms
	    	 * SPECIES		KEGG	ID	Comments
	    	 * C albicans	cal	EXTERNAL
	    	 * S pombe		spo	ENSEMBL
	    	 * S cerevisiae	sce	ENSEMBL
	    	 * C elegans	cel	ENSEMBL	adds CELE_ before the id
	    	 * E coli		eco	EXTERNAL?
	    	 * 
	    	 * checked entrez organisms
	    	 * H sapiens	hsa	ENTREZ
	    	 * M musculus	mmu	ENTREZ
	    	 */
	    	
	    	if(ed.organismKegg.equals("sce") || ed.organismKegg.equals("spo") || ed.organismKegg.equals("cel"))	
	    		{if(!ed.chip.equals("ensembl_gene_id") && ed.ensembl_gene_idHash!=null)	
	    			geneMap=ed.invertedHash(ed.ensembl_gene_idHash);}
	    	else if(ed.organismKegg.equals("cal"))
    			{if(!ed.chip.equals("external_gene_id") && ed.external_gene_idHash!=null)	
    				geneMap=ed.invertedHash(ed.ensembl_gene_idHash);}
	    	else								
	    		{if(!ed.chip.equals("entrezgene") && ed.entrezgeneHash!=null)	
	    			geneMap=ed.invertedHash(ed.entrezgeneHash);}
	    	}
	    }catch(Exception e){e.printStackTrace();}
	
	BufferedReader in;
	int numPaths=0;
	int cont=0;
	
	try {
		in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));//para applet/jws
		String s=null;
		String name=null;
		String id=null;
		OntologyTerm currentA=null, currentB=null, currentC=null, currentD=null;
		boolean ignore=false;
		
		while ((s = in.readLine()) != null)
			{
			switch(s.charAt(0))
				{
				case 'A':
					id=s.substring(0, s.indexOf(" ")).trim();
					name=s.substring(s.indexOf(" ")+1).trim();
					if(!brite) //PATHWAYS
						{
						if(id.equals("A09180") || id.equals("A09190")) //ignore "brite hierarchies" and "not included in pathway or brite" (drop out 40%, which is key for jnlp production tool)
							{ignore=true; break;}
						else 
							ignore=false;
						}
					else //BRITE
						{
						if(!id.equals("A09180"))
							{ignore=true; break;}
						else 
							ignore=false;
						}		
					System.out.println(name+"\t"+id);
					currentA=new OntologyTerm(name, id);
					map.put(currentA, new TreeMap<OntologyTerm, TreeMap>());
					break;
				case 'B':
					if(s.length()==1 || ignore)	break;
					
					s=s.substring(s.indexOf("B ")+2).trim();
					id=s.substring(0, s.indexOf(" ")).trim();
					name=s.substring(s.indexOf(" ")+1).trim();
					
					System.out.println("\t"+name+"\t"+id);
					currentB=new OntologyTerm(name, id);
					map.get(currentA).put(currentB, new TreeMap<OntologyTerm, TreeMap>());//Cannot be done with ontologies!!!
					break;
				case 'C'://pathway
					if(ignore)	break;
					s=s.substring(2).trim();
					id=s.substring(0, s.indexOf(" ")).trim();
					name=s.substring(s.indexOf(" ")+1).trim();
					if(name.contains("["))	name=name.substring(0, name.indexOf("["));
					
					currentC=new OntologyTerm(name, id);
					//if(name.contains("Fructose"))
						System.out.println("\t\t"+name+"\t"+id);
					//TODO ignore if not the one and so on
					if(brite || pathIDs==null || pathIDs.contains(id))	
						{
						((Map<OntologyTerm, TreeMap>)((Map<OntologyTerm, TreeMap>)map.get(currentA)).get(currentB)).put(currentC, new TreeMap<OntologyTerm, Map>());
						System.out.println("path "+name+" available for this organism");
						numPaths++;
						addElements=true;
						}
					else
						{
						addElements=false;
						//System.err.println("path "+name+" not available for this organism");
						}
					//System.out.print("\n"+name+":"+id+":");
					break;
				case 'D'://element
					if(!addElements || ignore)	break;
					
					name=s.substring(s.indexOf("K")+6).trim();
					id=s.substring(s.indexOf("K"), s.indexOf("K")+6).trim();
					if(name.contains("["))	name=name.substring(0, name.indexOf("["));
					//System.out.println("\t\t\t"+name);
					
					ArrayList<String> gids=new ArrayList<String>();
					if(koterms.containsKey("ko:"+id))
						gids=koterms.get("ko:"+id);	
					//else
					//	System.err.println("ko:"+id+" not found, should be added");
					
					//For quicker update, we can put the if(update) block here, although it won't consider
					//already added KOs, which might have added/removed some genes...
					//System.out.println(id);
					//if(id.equals("K10361"))
					//	update=true;
					
					if(update)
						{
						ArrayList<String> results=getGenesByKO(id);
						if(gids.size()!=results.size() || results.size()==0)
							{	
							System.out.println("Updating "+id+" "+gids.size()+" -> "+results.size());
							gids.clear();
							gids.addAll(results);
							koterms.put("ko:"+id, gids);
							cont++;
							if(cont%500==0)
								{
								System.out.println("\nSaving..."+cont);
								if(!brite)
									save(koterms, "komappingTOTAL-"+cont+".ser");
								else
									save(koterms, "britemappingTOTAL-"+cont+".ser");
								}
							}
						//else
						//	System.out.println("No need to update "+id);
						}
					
					if(organism.equals("kke"))
						gids.add("kke:"+id);
					else if(organism.equals("kegg_pathway"))
						gids.add("kegg_pathway:"+id);
					
					
					//TRANSLATION to gene ids on the expression data
					if(geneMap!=null)
						{
						ArrayList<String> gids2=new ArrayList<String>();
						for(String gid:gids)
							gids2.add(gid.substring(0, gid.indexOf(":"))+":"+geneMap.get(gid.substring(gid.indexOf(":")+1).toLowerCase()));
						gids=gids2;
						}
					
					
					//This should be something recursive on deeper ontologies
					Iterator<OntologyTerm> it=map.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentA.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism+":"))
										{
										t.geneIds.add(gid);
										l.add(gid);
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					
					Map<OntologyTerm, TreeMap> ma=map.get(currentA);
					it=ma.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentB.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism+":"))
										{
										t.geneIds.add(gid);
										l.add(gid);
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					Map<OntologyTerm, TreeMap> mb=ma.get(currentB);
					it=mb.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentC.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism+":"))
										{
										t.geneIds.add(gid);
										l.add(gid);
										//System.out.print(gid+" ");
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					
					break;
				default:
					break;
				}
			}
		
		//update=true;
		if(update)
			{
		//	save(koterms, "komappingTOTAL-new.ser");
			if(!brite)
				saveByOrganism(koterms, "komap");
			else
				saveByOrganism(koterms, "brite");
			}
		
		
	//saveTabDelimited(koterms, "KEGGalbicans");
			
	} catch (FileNotFoundException e) {
		e.printStackTrace();
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	} catch (Exception e){e.printStackTrace();}
	
	TreeMap<OntologyTerm, TreeMap> onto=new TreeMap<OntologyTerm, TreeMap>();
	onto.put(new OntologyTerm("KEGG Orthology", ""), map);
	System.out.println("Number of pathways parsed: "+numPaths);
	return onto;
	}
	
	/**
	 * For KEGG ontolgy version from 2012
	 * @param path
	 * @param organism
	 * @param pathIDs
	 * @param ed
	 * @param update
	 * @return
	 *//*
	public static TreeMap<OntologyTerm, TreeMap> parseV0(String path, String organism, ArrayList<String> pathIDs, ExpressionData ed, boolean update)
	{
	TreeMap<OntologyTerm, TreeMap> map=new TreeMap<OntologyTerm, TreeMap>();
	TreeMap<String, ArrayList<String>> koterms=new TreeMap<String, ArrayList<String>>();//mapping of KO terms to gene ids
	
	FileInputStream fis = null;
	ObjectInputStream oin = null;
	boolean addElements=true;
	HashMap<String, String> geneMap=null;
	
	try
	   {
		//if(!update)	//For an incremental update, we can comment this if and read the previous serialization
			{
			System.out.println("Reading komappingTOTAL");
			InputStream is=Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/kegg/komappingTOTAL.ser");
			oin = new ObjectInputStream(is);
	        koterms = (TreeMap<String, ArrayList<String>>)oin.readObject();
		    oin.close();
		    System.out.println("koterms read");
			}
	    if(ed!=null)
	    	{
	    	/*
	    	 * KEGG specials. Except otherwise, KEGG genes are ENTREZ for all the organisms
	    	 * SPECIES		KEGG	ID	Comments
	    	 * C albicans	cal	EXTERNAL
	    	 * S pombe		spo	ENSEMBL
	    	 * S cerevisiae	sce	ENSEMBL
	    	 * C elegans	cel	ENSEMBL	adds CELE_ before the id
	    	 * E coli		eco	EXTERNAL?
	    	 * 
	    	 * checked entrez organisms
	    	 * H sapiens	hsa	ENTREZ
	    	 * M musculus	mmu	ENTREZ
	    	 *
	    	
	    	if(ed.organismKegg.equals("sce") || ed.organismKegg.equals("spo") || ed.organismKegg.equals("cel"))	
	    		{if(!ed.chip.equals("ensembl_gene_id") && ed.ensembl_gene_idHash!=null)	
	    			geneMap=ed.invertedHash(ed.ensembl_gene_idHash);}
	    	else if(ed.organismKegg.equals("cal"))
    			{if(!ed.chip.equals("external_gene_id") && ed.external_gene_idHash!=null)	
    				geneMap=ed.invertedHash(ed.ensembl_gene_idHash);}
	    	else								
	    		{if(!ed.chip.equals("entrezgene") && ed.entrezgeneHash!=null)	
	    			geneMap=ed.invertedHash(ed.entrezgeneHash);}
	    	}
	    }catch(Exception e){e.printStackTrace();}
	
	BufferedReader in;
	int numPaths=0;
	int cont=0;
	
	try {
		in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));//para applet/jws
		String s=null;
		String name=null;
		String id=null;
		OntologyTerm currentA=null, currentB=null, currentC=null, currentD=null;
		
		while ((s = in.readLine()) != null)
			{
			switch(s.charAt(0))
				{
				case 'A':
					name=s.substring(s.indexOf("<b>")+3, s.indexOf("</b>")).trim();
					//System.out.println(name);
					currentA=new OntologyTerm(name, "");
					map.put(currentA, new TreeMap<OntologyTerm, TreeMap>());
					break;
				case 'B':
					if(s.length()==1)	break;
					
					name=s.substring(s.indexOf("<b>")+3, s.indexOf("</b>")).trim();
		//			System.out.println("\t"+name);
					currentB=new OntologyTerm(name, name);
					map.get(currentA).put(currentB, new TreeMap<OntologyTerm, TreeMap>());//Cannot be done with ontologies!!!
					break;
				case 'C'://pathway
					name=s.substring(s.indexOf("0")+5).trim();
					id=s.substring(s.indexOf("0"),s.indexOf("0")+5);
					if(name.contains("["))	name=name.substring(0, name.indexOf("["));
					currentC=new OntologyTerm(name, id);
	//				System.out.println("\t\t"+name+"\t"+id);
					//TODO ignore if not the one and so on
					if(pathIDs==null || pathIDs.contains(id))	
						{
						((Map<OntologyTerm, TreeMap>)((Map<OntologyTerm, TreeMap>)map.get(currentA)).get(currentB)).put(currentC, new TreeMap<OntologyTerm, Map>());
						//System.out.println("path "+name+" available for this organism");
						numPaths++;
						addElements=true;
						}
					else
						{
						addElements=false;
						//System.err.println("path "+name+" not available for this organism");
						}
					//System.out.print("\n"+name+":"+id+":");
					break;
				case 'D'://element
					if(!addElements)	break;
					
					name=s.substring(s.indexOf("K")+6).trim();
					id=s.substring(s.indexOf("K"), s.indexOf("K")+6).trim();
					if(name.contains("["))	name=name.substring(0, name.indexOf("["));
					//System.out.println("\t\t\t"+name);
					
					ArrayList<String> gids=new ArrayList<String>();
					if(koterms.containsKey("ko:"+id))
						gids=koterms.get("ko:"+id);	
					else
						{
						System.err.println("ko:"+id+" not found, should be added");
						//For quicker update, we can put the if(update) block here, although it won't consider
						//already added KOs, which might have added/removed some genes...
						if(update)
							{
							ArrayList<String> results=getGenesByKO(id);
							if(gids.size()!=results.size() || results.size()==0)
								{	
							//	System.out.println("Updating "+id);
							//	System.out.print("�");
								gids.clear();
								gids.addAll(results);
								koterms.put("ko:"+id, gids);
								cont++;
								if(cont%1000==0)
									{
									System.out.println("\nSaving..."+cont);
									save(koterms, "komappingTOTAL-"+cont+".ser");
									}
								}
							//else
							//	System.out.println("No need to update "+id);
							}
						}
					if(organism.equals("kke"))
						gids.add("kke:"+id);
					else if(organism.equals("kegg_pathway"))
						gids.add("kegg_pathway:"+id);
					
					//TRANSLATION to gene ids on the expression data
					if(geneMap!=null)
						{
						ArrayList<String> gids2=new ArrayList<String>();
						for(String gid:gids)
							gids2.add(gid.substring(0, gid.indexOf(":"))+":"+geneMap.get(gid.substring(gid.indexOf(":")+1).toLowerCase()));
						gids=gids2;
						}
					
					
					//This should be something recursive on deeper ontologies
					Iterator<OntologyTerm> it=map.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentA.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism+":"))
										{
										t.geneIds.add(gid);
										l.add(gid);
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					
					Map<OntologyTerm, TreeMap> ma=map.get(currentA);
					it=ma.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentB.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism+":"))
										{
										t.geneIds.add(gid);
										l.add(gid);
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					Map<OntologyTerm, TreeMap> mb=ma.get(currentB);
					it=mb.keySet().iterator();
					while(it.hasNext())
						{
						OntologyTerm t=it.next();
						if(t.name.equals(currentC.name))
							{
							if(organism==null)
								{
								t.geneIds.addAll(gids);
								t.koIds.put(id, gids);
								}
							else
								{
								ArrayList<String> l=new ArrayList<String>();
								for(String gid:gids)
									if(gid.startsWith(organism+":"))
										{
										t.geneIds.add(gid);
										l.add(gid);
										//System.out.print(gid+" ");
										}
								t.koIds.put(id, l);
								}
							break;
							}
						}
					
					break;
				default:
					break;
				}
			}
		if(update)
			{
			save(koterms, "komappingTOTAL-new.ser");
			}
		
	//saveTabDelimited(koterms, "KEGGalbicans");
			
	} catch (FileNotFoundException e) {
		e.printStackTrace();
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	} catch (Exception e){e.printStackTrace();}
	
	TreeMap<OntologyTerm, TreeMap> onto=new TreeMap<OntologyTerm, TreeMap>();
	onto.put(new OntologyTerm("KEGG Orthology", ""), map);
	System.out.println("Number of pathways parsed: "+numPaths);
	return onto;
	}*/

	public static void save(TreeMap<String, ArrayList<String>> map, String fileName) throws FileNotFoundException, IOException
		{
		System.out.println("Saving the whole thing");
		//Serialize the mapping for optimization
		FileOutputStream fos = new FileOutputStream(fileName);
		ObjectOutputStream out = new ObjectOutputStream(fos);
		out.writeObject(map);
		out.close();
		}
	
	/** As save, but generates one file per organism (given the current size of KEGG -2018-, it is necessary)
	 * 
	 * @param map
	 * @param fileName
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static void saveByOrganism(TreeMap<String, ArrayList<String>> map, String fileName) throws FileNotFoundException, IOException
		{
		System.out.println("Saving KEGG annotations by organism");
		HashMap<String,String> orgs=KOparser.getKeggOrganismMap();
		//Serialize the mapping for optimization
		for(String org:orgs.values())
			{
			FileOutputStream fos = new FileOutputStream(fileName+"-"+org+".ser");
			ObjectOutputStream out = new ObjectOutputStream(fos);
			TreeMap<String,ArrayList<String>> omap=KOparser.getMapByOrganism(map, org);
			out.writeObject(omap);
			out.close();
			}
		}
	
	public static TreeMap<String,ArrayList<String>> getMapByOrganism(TreeMap<String,ArrayList<String>> map,String org)
		{
		TreeMap<String,ArrayList<String>> omap=new TreeMap<String,ArrayList<String>>();
		for(String ko:map.keySet())
			{
			ArrayList<String> gl=new ArrayList<String>();
			for(String gene:map.get(ko))
				{
				if(gene.startsWith(org))
					gl.add(gene);
				}
			if(gl.size()>0)
				omap.put(ko, gl);
			}
		return omap;
		}
	
	/**
	 * Returns a map where the keys are the standard organism names and the values are the corresponding KEGG ids
	 * Works with an internal file generated from the query http://rest.kegg.jp/list/organism (performed on dec2018)
	 * @return
	 */
	public static HashMap<String,String> getKeggOrganismMap()
		{
		String cad="";
		HashMap <String,String> map=new HashMap<String,String>();
		
		try {		
			
		BufferedReader file = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream("es/usal/voronto/data/kegg/organisms.tab")));//para applet/jws
			
		//BufferedReader file = new BufferedReader(new FileReader());//File retrieved from http://rest.kegg.jp/list/organism on dec2018
		while((cad=file.readLine())!=null)
			{
			String[] fields=cad.split("\t");
			//System.out.println(cad);
			String id=fields[0];
			String org=fields[1];
			org=org.replaceAll("\\(.*\\)", "").trim();
			map.put(org, id);
			}
		file.close();
		}catch(Exception e){e.printStackTrace();}
		return map;
		}

	/*public static void saveTabDelimited(TreeMap <OntologyTerm, TreeMap> map, String fileName) throws FileNotFoundException, IOException
		{
		BufferedWriter out= new BufferedWriter(new FileWriter(fileName));
		out.write("term\tgene");
		for(OntologyTerm term:map.keySet())
			{
			for(String gene:map.get(term))
				out.write(term+"\t"+gene);
			}
		out.close();
		}*/
	/*public static void saveTabDelimited(TreeMap<String, ArrayList<String>> map, String fileName) throws FileNotFoundException, IOException
		{
		BufferedWriter out= new BufferedWriter(new FileWriter(fileName));
		out.write("term\tgene");
		for(String term:map.keySet())
			{
			for(String gene:map.get(term))
				out.write(term+"\t"+gene);
			}
		out.close();
		}*/
		
	public static ArrayList<String> getGenesByKO(String ko)
		{
		ArrayList<String> genes=new ArrayList<String>();
		
		try{
			URL url = new URL("http://rest.kegg.jp/link/genes/ko:"+ko);
			HttpURLConnection conn = (HttpURLConnection) url.openConnection(); 
			conn.setRequestMethod("GET");
			BufferedReader br = new BufferedReader(new InputStreamReader((conn.getInputStream())));
			String output; 
			while ((output = br.readLine()) != null) 
				{
				if(output.indexOf("\t")>=0)
					{	
					String g=output.split("\t")[1].trim();
					//if(g.startsWith("sal"))
					//	System.out.println("here");
					genes.add(g);
					}
				}
			 conn.disconnect();
			}catch(Exception e){e.printStackTrace();}
		//TODO: add itself as genes for metagenome/non-organism specific analyses
		genes.add("KOID:"+ko);
		return genes;
		}
	
	public static ArrayList<String> getKOsByPath(String path)
	{
	ArrayList<String> kos=new ArrayList<String>();
	
	try{
		URL url = new URL("http://rest.kegg.jp/link/ko/map"+path);
		HttpURLConnection conn = (HttpURLConnection) url.openConnection(); 
		conn.setRequestMethod("GET");
		BufferedReader br = new BufferedReader(new InputStreamReader((conn.getInputStream())));
		String output; 
	
		while ((output = br.readLine()) != null) 
			{
			String ko=output.split("\t")[1].trim().split(":")[1];
			kos.add(ko);
			}
		 conn.disconnect();
		}catch(Exception e){e.printStackTrace();}
	
	return kos;
	}
}

