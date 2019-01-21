package es.usal.voronto.model.filters;

import java.io.File;

import javax.swing.filechooser.FileFilter;


/**
 * OBO file format filter for JFileChooser
 * @author Rodrigo	Santamaria
 */
public class OBODataFilter extends FileFilter{
    static final String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
    }
    
  
    /**
     * Decides if a file is acceptable as Microarray file
     * @param	f	File to check
     * @return	true if the file extension is txt, false otherwise
     */
    public boolean accept(File f) 
    	{
        if (f.isDirectory()) 							       return true;

        String extension = getExtension(f);
        if (extension != null) 
        	{
            if (extension.equals("obo"))                  return true;
            else								          return false;
            }

        return false;
    	}

    //The description of this filter
    /**
     * Returns the description of Microarray files
     * @return	A brief String description of expected files for Microarray.
     */
    public String getDescription() {
        return "OBO format (.obo)";
    }
}