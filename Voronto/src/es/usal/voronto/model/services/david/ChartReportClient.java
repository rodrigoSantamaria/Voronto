package es.usal.voronto.model.services.david;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.TreeMap;

import javax.swing.JOptionPane;
//import david.xsd.*;
//import org.apache.axis2.AxisFault;

import es.usal.voronto.model.ontology.OntologyTerm;
import sample.session.service.xsd.*;
import sample.session.client.stub.DAVIDWebServiceStub;

public class ChartReportClient {

    public SimpleChartRecord[] invokeService(String inputIds, String idType, String category_names, String saveFile, float thrd, int ct) throws Exception{
      System.err.println("This method is by now deprecated");
      
      /*try {
            DAVIDWebServiceStub stub =
                    new DAVIDWebServiceStub(
                    "http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/");
            stub._getServiceClient().getOptions().setManageSession(true);
        	//System.out.println();
         	//user authentication by email address
  			//For user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm
 	    	String userAuthentication= (String) stub.authenticate("rodri@usal.es");
 	    	//System.out.println("User authentication: " +  userAuthentication);
 	    	//System.out.println();
 	   	  if (userAuthentication.equals("true"))
 	   		 {
				//list conversion types
				String conversionTypes = stub.getConversionTypes();
	 			//System.out.println("Conversion Types:\n"+ conversionTypes + "\n");

				//addlist
	        	//String inputIds = new String("1112_g_at,1331_s_at,1355_g_at,1372_at,1391_s_at,1403_s_at,1419_g_at,1575_at,1645_at,1786_at,1855_at,1890_at,1901_s_at,1910_s_at,1937_at,1974_s_at,1983_at,2090_i_at,31506_s_at,31512_at,31525_s_at,31576_at,31621_s_at,31687_f_at,31715_at,31793_at,31987_at,32010_at,32073_at,32084_at,32148_at,32163_f_at,32250_at,32279_at,32407_f_at,32413_at,32418_at,32439_at,32469_at,32680_at,32717_at,33027_at,33077_at,33080_s_at,33246_at,33284_at,33293_at,33371_s_at,33516_at,33530_at,33684_at,33685_at,33922_at,33963_at,33979_at,34012_at,34233_i_at,34249_at,34436_at,34453_at,34467_g_at,34529_at,34539_at,34546_at,34577_at,34606_s_at,34618_at,34623_at,34629_at,34636_at,34702_f_at,34703_f_at,34720_at,34902_at,34972_s_at,35038_at,35069_at,35090_g_at,35091_at,35121_at,35169_at,35213_at,35367_at,35373_at,35439_at,35566_f_at,35595_at,35648_at,35896_at,35903_at,35915_at,35956_s_at,35996_at,36234_at,36317_at,36328_at,36378_at,36421_at,36436_at,36479_at,36696_at,36703_at,36713_at,36766_at,37061_at,37096_at,37097_at,37105_at,37166_at,37172_at,37408_at,37454_at,37711_at,37814_g_at,37898_r_at,37905_r_at,37953_s_at,37954_at,37968_at,37983_at,38103_at,38128_at,38201_at,38229_at,38236_at,38482_at,38508_s_at,38604_at,38646_s_at,38674_at,38691_s_at,38816_at,38926_at,38945_at,38948_at,39094_at,39187_at,39198_s_at,39469_s_at,39511_at,39698_at,39908_at,40058_s_at,40089_at,40186_at,40271_at,40294_at,40317_at,40350_at,40553_at,40735_at,40790_at,40959_at,41113_at,41280_r_at,41489_at,41703_r_at,606_at,679_at,822_s_at,919_at,936_s_at,966_at,967_g_at");
				//String idType = new String("AFFYMETRIX_3PRIME_IVT_ID");
				String listName = new String("make_up");
				int listType = 0;
	    		double addListOutput =0;

	    		//Set user defined categories
				//String category_names = new String("BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE");
				//to use DAVID dafault categories, set category_names to be empty
				//String category_names = new String("");
				String validCategoryString = stub.setCategories(category_names);


				try{
					addListOutput = (double) stub.addList(inputIds, idType, listName, listType);
					String allListName_str = (String)stub.getAllListNames();
					String[] allListNames = allListName_str.split(",");
					String allSpecies = stub.getSpecies();
					String[] allSpecies_name = allSpecies.split(",");
					//System.out.println("All " + allSpecies_name.length + " species:\n");
					String currentSpecies = stub.getCurrentSpecies();
					//System.out.println("Current species: "+ currentSpecies + "\n");
					String[] currentSpecies_name = currentSpecies.split(",");
					//System.out.println("Total " + currentSpecies_name.length+" Current species\n");
					//String newCurrentSpecies = stub.setCurrentSpecies("1");
					//System.out.println("Set new current species to: "+ newCurrentSpecies + "\n");
					 
					//double thrd = 0.1;
					//int ct = 2;
					SimpleChartRecord[] simpleChartRecords=(SimpleChartRecord[])stub.getChartReport(thrd,ct);
					if (simpleChartRecords!=null && simpleChartRecords.length >0)
						{
						if(saveFile!=null)
							{
							StringBuffer sb;
							sb=new StringBuffer();
					 	    
							sb.append("Category\tTerm\tCount\t%\tPValue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n");
							writeSimpleChartRecords(simpleChartRecords,sb);
							
							PrintWriter outfile = new PrintWriter(new FileOutputStream(saveFile));
							outfile.println(sb);
							outfile.close();
							System.out.println(saveFile+" generated.");
							}
						return simpleChartRecords;
						}
					else
						{ 
						JOptionPane.showMessageDialog(null, "No enriched terms found by DAVID for the selected parameters", "No enriched terms", JOptionPane.INFORMATION_MESSAGE);
						System.out.println("0 chart record output in the report."); 
						return null;
						}
					
				//}catch (java.io.IOException e) {e.printStackTrace(); return null;
		 		}catch (Exception e) {e.printStackTrace(); throw e; }
			}
		}catch (Exception e) {e.printStackTrace(); throw e; }
		// }catch (AxisFault axisFault) { axisFault.printStackTrace(); return null;
        //} catch (java.rmi.RemoteException e) {e.printStackTrace(); return null;}
         * */
         
     return null;
    }

  protected void writeSimpleChartRecords(SimpleChartRecord[] simpleChartRecords, StringBuffer sb) {

			for (int j = 0; j < simpleChartRecords.length; j++)
			{
				sb.append(simpleChartRecords[j].getCategoryName()+"\t"
				          +simpleChartRecords[j].getTermName()+"\t"
				          +simpleChartRecords[j].getListHits()+"\t"
				          +simpleChartRecords[j].getPercent()+"\t"
				          +simpleChartRecords[j].getEase()+"\t"
				          +simpleChartRecords[j].getGeneIds()+"\t"
				          +simpleChartRecords[j].getListTotals()+"\t"
				          +simpleChartRecords[j].getPopHits()+"\t"
				          +simpleChartRecords[j].getPopTotals()+"\t"
				          +simpleChartRecords[j].getFoldEnrichment()+"\t"
				          +simpleChartRecords[j].getBonferroni()+"\t"
				          +simpleChartRecords[j].getBenjamini()+"\t"
				          +simpleChartRecords[j].getAfdr()+"\t"
				          +"\n");
			}
			sb.append("\n");
	}

    /*public static void main(String[] args) throws java.io.IOException {
        new ChartReportClient().invokeService();
    }*/
}
