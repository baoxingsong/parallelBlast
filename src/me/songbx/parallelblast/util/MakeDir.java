package me.songbx.parallelblast.util;

import java.io.File;

public class MakeDir{
	private String path;
	public synchronized void makeDirs(String path){
		this.path=path;
		File dir = new File(path);
		if(dir!=null && !dir.exists()){
			dir.mkdirs();
		}
	}
	
	public synchronized void makeParentDirs(String path){
		this.path=path;
		File file = new File(path);
		File parent = file.getParentFile();
		if(parent!=null && !parent.exists()){
			parent.mkdirs();
		}
	}
	
	public static void main(String args[]){
		MakeDir mkd=new MakeDir();
		mkd.makeParentDirs("D:\\songs2\\BLASTTEST\\seperate");
	}
	
	public String getPath() {
		return path;
	}

	public void setPath(String path) {
		this.path = path;
	}
}
