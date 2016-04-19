/**
 * 
 */
package me.songbx.parallelblast.util;

import java.nio.channels.FileChannel;
import java.io.FileInputStream;  
import java.io.FileOutputStream;  
import java.io.IOException;  
import java.nio.ByteBuffer;

/**
 * @author Baoxing SONG
 * @mail songbaoxing168@163.com
 * @version 0.1
 * @since 04/10/2013
 */
public class AddToFile {

	private static final int BUFSIZE = 1024 * 8;
	
	/**
	 * @param toBeCombinedfilePath
	 * @param poolFilePath
	 */
	
	/**
	 * @param 
	 */
	@SuppressWarnings("resource")
	public static synchronized void addTo(String toBeCombinedfilePath, String poolFilePath){

		FileChannel outChannel = null;
		try {
			
			FileOutputStream fo=new FileOutputStream(poolFilePath, true);
			outChannel = fo.getChannel();
            FileChannel fc = new FileInputStream(toBeCombinedfilePath).getChannel();   
            ByteBuffer bb = ByteBuffer.allocate(BUFSIZE);  
            while(fc.read(bb) != -1){  
                bb.flip();  
                outChannel.write(bb);  
                bb.clear();  
            }
            fc.close();  
        } catch (IOException ioe) {  
            ioe.printStackTrace();  
        } finally {  
            try {if (outChannel != null) {outChannel.close();}} catch (IOException ignore) {}  
        }
	}

	public static int getBufsize() {
		return BUFSIZE;
	}
	
	

}
