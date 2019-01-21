package es.usal.voronto.model.ontology;

import java.util.ArrayList;
import java.util.TreeMap;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


/**
 * Simple XML ontology parser.
 * This is a specific Voronto format, but quite easy to follow:
 * <Ontology>
 * 	<Term id="34" name="Metabolism">
 * 		<Term id="3432" name="Metabolism of carbohydrates" link="http://kegg/3432">
 * 			<Element id="YBL3453-C"/>
 * 			<Element id="YBL3433-C"/>
 * 			...
 * 		</Term>
 *   ...
 *  </Term>
 * </Ontology>
 * 
 * The recursion depth is arbitrary. The leaf Terms musth contain Elements (usually genes) with ids that correspond to the expression matrix that is to be mapped to.
 * link attribute on leaf Terms is optional.
 * @author rodri
 *
 */
public class XMLparser extends DefaultHandler{
	
	TreeMap<OntologyTerm, TreeMap> map;
	private String tempVal;
	private TreeMap<OntologyTerm, TreeMap> currentMap;
	private ArrayList previousMaps=new ArrayList();
	private OntologyTerm currentTerm;
	public OntologyTerm root;
	
	public TreeMap<OntologyTerm, TreeMap> parse(String path) throws Exception
		{
		SAXParserFactory spf = SAXParserFactory.newInstance();
		SAXParser sp = spf.newSAXParser();
		sp.parse(path, this);

		return map.get(root);
		}
	
	//Event Handlers
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException 
		{
		//reset
		if(qName.equalsIgnoreCase("Ontology")) 
			{
			root=new OntologyTerm(attributes.getValue("name"),"");
			map=new TreeMap<OntologyTerm, TreeMap>();
			previousMaps.add(map);
			currentMap= new TreeMap<OntologyTerm, TreeMap>();
			map.put(root,currentMap);
			}
		
		if(qName.equalsIgnoreCase("Term"))
			{
			OntologyTerm ot=new OntologyTerm(attributes.getValue("name"), attributes.getValue("id"));
			TreeMap<OntologyTerm, TreeMap> m=new TreeMap<OntologyTerm, TreeMap>();
			currentMap.put(ot, m);
			previousMaps.add(currentMap);
			currentMap=m;
			currentTerm=ot;
			}
		if(qName.equalsIgnoreCase("Element"))
			{
			System.out.println("Adding element "+attributes.getValue("id")+" to term "+currentTerm.name);
			currentTerm.geneIds.add(attributes.getValue("id"));
			}
	
		}
	public void characters(char[] ch, int start, int length) throws SAXException 
		{
		tempVal = new String(ch,start,length);
		}

	public void endElement(String uri, String localName, String qName) throws SAXException 
		{
		if(qName.equalsIgnoreCase("Term"))
			currentMap=(TreeMap<OntologyTerm,TreeMap>)previousMaps.remove(previousMaps.size()-1);
		//else if(qName.equalsIgnoreCase("Ontology"))
		//	map=(TreeMap<OntologyTerm,TreeMap>)previousMaps.remove(previousMaps.size()-1);

		/*
		else if (qName.equalsIgnoreCase("id")) 			c.id=Integer.parseInt(tempVal);
		else if (qName.equalsIgnoreCase("anno")) 		c.anno=Integer.parseInt(tempVal);
		else if (qName.equalsIgnoreCase("autor"))		c.autor=tempVal;
		else if (qName.equalsIgnoreCase("titulo"))		c.titulo=tempVal;
		else if (qName.equalsIgnoreCase("pais"))		c.pais=tempVal;
		else if (qName.equalsIgnoreCase("tema"))		c.tema=tempVal;
		else if (qName.equalsIgnoreCase("publicacion"))	c.publicacion=tempVal;
		*/
		}

}
