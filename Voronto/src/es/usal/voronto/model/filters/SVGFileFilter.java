package es.usal.voronto.model.filters;

import java.io.File;

import javax.swing.filechooser.FileFilter;


public class SVGFileFilter extends FileFilter {

	 /**
     * Returns the extension of a file (the three letters after the dot)
     * @param f	File to get the extension
     * @return	extension of f
     */
	public final String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 &&  i < s.length() - 1) {
            ext = s.substring(i+1).toLowerCase();
        }
        return ext;
    }
	

    
    /**
     * Decides if a file is acceptable as TRN file
     * @param	f	File to check
     * @return	true if the file extension is xml or gml, false otherwise
     */
    public boolean accept(File f) 
    	{
        if (f.isDirectory()) 							       return true;

        String extension = getExtension(f);
        if (extension != null) 
        	{
        	if (extension.equals("svg"))                  return true;
        	else							              return false;
            }

        return false;
    	}

    //The description of this filter
    /**
     * Returns the description of TRN files
     * @return	A brief description of expected files for TRN.
     */
    public String getDescription() {
        return "Scalable vector file (.svg)";
    }
}

