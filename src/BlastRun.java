import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;

import me.songbx.parallelblast.util.AddToFile;
import me.songbx.parallelblast.util.DeleteFile;
import Model.MyThreadCount;


public class BlastRun extends Thread {

	private MyThreadCount c;
	private String uuid;
	private String blastPath;
	private String tempPath;
	private String outPath;
	private ArrayList<String> cmdArrayList;
	
	public BlastRun(MyThreadCount c, String uuid, String blastPath, String tempPath, String outPath, ArrayList<String> cmdArrayList){
		this.c=c;
		this.uuid=uuid;
		this.blastPath=blastPath;
		this.tempPath=tempPath;
		this.outPath=outPath;
		this.cmdArrayList=cmdArrayList;
	}
	public void run(){
		
		String[] cmd=new String[cmdArrayList.size() + 5];
		cmd[0] = blastPath + File.separator + "blastall";
		for(int i=1; i<=cmdArrayList.size(); i++){
			cmd[i] = cmdArrayList.get(i-1);
		}
		cmd[cmdArrayList.size()+1]="-o";
		cmd[cmdArrayList.size()+2]=tempPath + uuid + ".blast";
		cmd[cmdArrayList.size()+3]="-i";
		cmd[cmdArrayList.size()+4]=tempPath + uuid;
		
		Runtime run = Runtime.getRuntime();
		try {
			//System.out.println(cmdArrayList.size());
			//System.out.println(cmd[0] + " " + cmd[1] + " " + cmd[2] + " " + cmd[3] + " " + cmd[4] + " " + cmd[5] + " " + cmd[6] + " " + cmd[7] + " " + cmd[8]);
			
            Process p = run.exec(cmd);
            BufferedInputStream in = new BufferedInputStream(p.getInputStream());  
            BufferedReader inBr = new BufferedReader(new InputStreamReader(in));  
            String lineStr;
            while ((lineStr = inBr.readLine()) != null){  
            	System.out.println(lineStr);
            }
            
            if (p.waitFor() != 0) {
                if (p.exitValue() == 1)  
                System.err.println("BLAST error!");  
            }
            inBr.close();  
            in.close();
            p.destroy();
            
        } catch (Exception e) {  
            e.printStackTrace();  
        }
		//System.out.println(tempPath + uuid);
		//System.out.println(DeleteFile.delete(tempPath + uuid));
		
		AddToFile.addTo(tempPath + uuid + ".blast", outPath);
		DeleteFile.delete(tempPath + uuid + ".blast");
		DeleteFile.delete(tempPath + uuid);
		c.countDown();
	}
}
