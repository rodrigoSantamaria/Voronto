package es.usal.voronto.model.ontology;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

//import edu.emory.mathcs.backport.java.util.Arrays;
import java.util.Arrays;

/**
 * This class provides parsing functionanlity for OBO. We prefer to use this simple parser that the complicated and sparsely documented
 * API ofr OBO-Edit.
 * @author rodri
 *
 */
public class OBOparser {

	/**
	 * 
	 * @param path to the OBO file (full path)
	 * @param root term
	 */
	public static TreeMap<OntologyTerm, TreeMap> parse(String path, String rootLabel, String namespaceLabel)
		{
		BufferedReader in;
		TreeMap<OntologyTerm, TreeMap> map=new TreeMap<OntologyTerm, TreeMap>();
		Map<OntologyTerm, List<String>> parents=Collections.synchronizedMap(new TreeMap<OntologyTerm, List<String>>());//NOTE: considering more than 1 parent per node
		Map<String, List<String>> children=Collections.synchronizedMap(new TreeMap<String, List<String>>());
		ArrayList<String> goTerms=new ArrayList<String>();
		
		try {
			if(rootLabel==null)	in = new BufferedReader(new FileReader(path));
			else				in = new BufferedReader(new InputStreamReader(Thread.currentThread().getContextClassLoader().getResourceAsStream(path)));			
			
			String s=null;
			String name=null;
			String id=null;
			String namespace=null;
			String def=null;
			
			//1) map ontology terms to their parents
			long time=System.currentTimeMillis();
			System.out.println("Building OBO custom ontology: "+path+" with root "+rootLabel+" and namespace "+namespaceLabel);
			while ((s = in.readLine()) != null)
				{
				if(s.equals("[Term]"))	
					{
					//B11) Take first three fixed fields (id, name and namespace)
					id=in.readLine();
					id=id.substring(id.indexOf(":")+1).trim();
					name=in.readLine();
					name=name.substring(name.indexOf(":")+1).trim();
					if(name.length()==0)	//the term is over, it only has id: field
						{
						name=id;
						goTerms.add(id);
						OntologyTerm ot=new OntologyTerm(name, id);
						parents.put(ot, null);
						}
					
					else
						{
						if(namespaceLabel!=null)
							{
							namespace=in.readLine();
							namespace=namespace.substring(namespace.indexOf(":")+1).trim();
							}
						
						OntologyTerm ot=new OntologyTerm(name, id);
						
						if(namespaceLabel==null || namespace.equals(namespaceLabel))	//only for the ontology we are populating
							{
							//B12) Then go to is_a relations
							String cad=null;
							boolean obsolete=false;
							
							//considering more than one parent
							List<String> p=new ArrayList<String>();
							while((cad=in.readLine())!=null && cad.length()>0)	
								{
								//System.out.println(cad);
								if(cad.startsWith("is_obsolete: true"))	//this is a non valid node
									{	obsolete=true; break;	}
								if(cad.startsWith("is_a:"))
									{
									if(cad.indexOf("!")>0)	p.add(cad.substring(cad.indexOf(":")+1, cad.indexOf("!")).trim());//this is the root, no is_a: relationship
									else					p.add(cad.substring(cad.indexOf(":")+1).trim());//this is the root, no is_a: relationship
									}
								else if(cad.startsWith("relationship: part_of"))
									{
									//System.out.println("adding "+cad);
									if(cad.indexOf("!")>0)	p.add(cad.substring(cad.indexOf("part_of")+1, cad.indexOf("!")).trim());//this is the root, no is_a: relationship
									else					p.add(cad.substring(cad.indexOf("part_of")+1).trim());//this is the root, no is_a: relationship
									}
								}
							
		
							if(!obsolete)
								{
								goTerms.add(id);
								parents.put(ot, p);
								
								//with children list here
								for(String pp:p)
									{
									List<String> c=children.get(pp);
									if(c!=null)	c.add(ot.id);
									else	{c=new ArrayList<String>(); c.add(ot.id); children.put(pp, c);}
									}
								}
							}
						}
					}//if [Term]
				}
			
			System.out.println("Ontology terms recovered in "+(System.currentTimeMillis()-time)/1000.0+"s, "+parents.size()+" elements");
			
			//2) get all the ones without a parent (root) and build from there on
			Iterator<OntologyTerm> it=parents.keySet().iterator();
			ArrayList<OntologyTerm> roots=new ArrayList<OntologyTerm>();
			while(it.hasNext())
				{   
				OntologyTerm ot=it.next();
				
				if(rootLabel!=null && ot.name.equals(rootLabel))
					{
					roots.add(ot);
					map.put(ot, new TreeMap<OntologyTerm, Map>());
					}
				else if(rootLabel==null)
					{
					if(parents.get(ot)==null || parents.get(ot).size()==0)
						{
						roots.add(ot);
						map.put(ot, new TreeMap<OntologyTerm, Map>());
						}
					}
				}
			System.out.println("There are "+map.size()+" roots");
			//only one should be in the map by now
			//System.out.println("Parent is "+root.id+", "+root.name);
			OntologyTerm[] ots=parents.keySet().toArray(new OntologyTerm[0]);
			ArrayList<OntologyTerm> lot=new ArrayList<OntologyTerm>(Arrays.asList(ots));
			HashMap<String, OntologyTerm> mlot=new HashMap<String, OntologyTerm>();
			for(OntologyTerm ot:ots)	mlot.put(ot.id, ot);
			
			for(OntologyTerm root:roots)	
				addChildren(root, children, map.get(root), mlot,0);
		} 
		catch (FileNotFoundException e) {	e.printStackTrace();	} 
		catch (IOException e) 			{	e.printStackTrace();	}
		catch (Exception e)				{	e.printStackTrace();	}
	
	return map;
	}
	
	/**
	 * Recursively builds the hierarchy, given the parents.
	 * It cannot be done from the beginning because of the GO file structure.
	 * @param root
	 * @param parents
	 * @param m
	 * @param lot
	 */
	private static synchronized void addChildren(OntologyTerm root, Map<String, List<String>> children, Map<OntologyTerm, Map> m, Map<String, OntologyTerm> lot, int level)
		{
		List<String> c=children.get(root.id);
		if(c!=null)
			{
			for(String cc:c)
				{
				OntologyTerm ochild=lot.get(cc);
				//for(int i=0;i<level;i++)	System.out.print("\t");
				//System.out.println(ochild.id+"\t"+ochild.name);
				if(ochild!=null)
					{
					m.put(ochild, new TreeMap<OntologyTerm, Map>());
					addChildren(ochild, children, m.get(ochild), lot, level+1);
					}
				}
			}
		return;
		}

}
