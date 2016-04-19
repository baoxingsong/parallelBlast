/**
 * 
 */
package me.songbx.parallelblast.util;

import java.io.File;

/**
 * @author Baoxing SONG
 * @mail songbaoxing168@163.com
 * @version 0.1
 * @since 04/10/2013
 */
public class DeleteFile {

	/**
	 * 
	 */
	public static synchronized boolean delete(String toBeDeleteFilePath){
		boolean flag = false;  
		File file = new File(toBeDeleteFilePath);  
	    if (file.isFile() && file.exists()) {  
	        file.delete();  
	        flag = true;  
	    }  
	    return flag; 
	}
}
