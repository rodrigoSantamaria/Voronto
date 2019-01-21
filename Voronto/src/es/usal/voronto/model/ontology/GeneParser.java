package es.usal.voronto.model.ontology;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

/**
 * This class reads a txt file searching for a list of genes.
 * If more than one line, it takes each line as a gene name.
 * If just one line, it considers the file as tab-delimited and takes each column as a gene name.
 * @author rodri
 *
 */
public class GeneParser {
	public static HashSet<String> parse(String path)
		{
		HashSet<String> genes=new HashSet<String>();
		
		try{
			BufferedReader in = new BufferedReader(new FileReader(path));
			String cad=in.readLine();
			String cad2=in.readLine();
			if(cad2==null || cad2.length()==0)	//genes in one line
				{
				String g[]=cad.split("\t");
				genes.addAll(Arrays.asList(g));	
				}
			else								//one gene per line
				{
				in = new BufferedReader(new FileReader(path));
				while((cad=in.readLine())!=null)
					genes.add(cad.trim());
				}
			in.close();
			return genes;
			}
		catch (FileNotFoundException e) {	e.printStackTrace();	return null; } 
		catch (IOException e) 			{	e.printStackTrace();	return null; }
		catch (Exception e)				{	e.printStackTrace();	return null; }
		}
	
}
