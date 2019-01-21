package es.usal.voronto.model.ontology;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.collections15.CollectionUtils;

import es.usal.voronto.model.stats.Stats;

public class OntologyTerm implements Serializable, Cloneable, Comparable<OntologyTerm> {
	public String name;
	public String id;
	
	public double pvalue=999;	//enrichment statistic if required
	public double pvalueBH=999;	//enrichment BH.corrected statistic if required
	public boolean significant=false; //true only if the p-value is below the given threshold, 
	public HashSet<String> degs=new HashSet<String>();		//differential expressed genes in the term, if required
	public int numRelevantSubterms=0;
	
	//public ArrayList<String> geneIds=new ArrayList<String>(); //gene ids annotated with this term;
	public HashSet<String> geneIds=new HashSet<String>(); //gene ids annotated with this term;
	public Map<String, ArrayList<Float>> geneExs=new TreeMap<String, ArrayList<Float>>();//expression for each gene term, on each condition
	public Map<String, ArrayList<String>> koIds=new TreeMap<String, ArrayList<String>>();//ko ids for kegg orthology, mapped to Kegg genes
	public Map<String, ArrayList<Float>> koExs=new TreeMap<String, ArrayList<Float>>();//expression for each ko term, on each condition
	
	public OntologyTerm(String name, String id)
		{
			this.name=name;
			this.id=id;
		}

	@Override
	public int compareTo(OntologyTerm arg0) {
		// TODO Auto-generated method stub
		return 10*name.compareTo(arg0.name)+id.compareTo(arg0.id);
		}
	public void resetEnrichment()
		{
		pvalue=999;
		pvalueBH=999;
		degs.clear();
		numRelevantSubterms=0;
		}
	
	public void setDegs(HashSet<String> degs, boolean computeEnrichment, int totalGenes)
		{
		HashSet<String> genesInTerm=new HashSet<String>();
		genesInTerm.addAll(CollectionUtils.intersection(degs, geneExs.keySet()));
		
		this.degs=genesInTerm;
		if(computeEnrichment)
			{
			if(name.contains("Meiosis"))
				System.out.println(geneExs.keySet().size()+", "+degs.size()+", "+genesInTerm.size()+", "+totalGenes);
			this.pvalue=Stats.fisherTest(geneExs.keySet().size(), degs.size(), genesInTerm.size(), totalGenes);
			if(this.pvalue==0)
				pvalue=Double.MIN_VALUE;
			this.pvalueBH=this.pvalue;
			}
		return;
		}
}
