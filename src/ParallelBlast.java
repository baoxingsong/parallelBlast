

import static java.lang.System.*; 

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Properties;
import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import Common.Constant;
import Model.MyThreadCount;
import me.songbx.parallelblast.error.NotFastaFormat;
import me.songbx.parallelblast.error.ParameterError;
import me.songbx.parallelblast.error.UnknownOs;
import me.songbx.parallelblast.util.MakeDir;

public class ParallelBlast {
	//parameters from BLAST
	private String program;
	private String database;
	private String inputFile;
	private Double expectationvalue;
	private int viewOptions;
	private String outPath;
	private String filterQuerySequence;
	private int costToOpenAGap;
	private int costToExtendAGap;
	private int xDropoffValueForGappedAlignment;
	private boolean showGIInDeflines;
	private int penaltyForANucleotideMismatch;
	private int rewardForANucleotideMatch;
	private int numberOfDatabaseSequencesToShowOneLineDescriptions;
	private int numberOfDatabaseSequenceToShowAlignmentsForB;
	private int thresholdForExtendingHits;
	private boolean performGappedAlignment;
	private int queryGeneticCodeToUse;
	private int dBGeneticCode;
	private int numberOfProcessorsToUse;
	private String seqAlignFile;
	private boolean believeTheQueryDefline;
	private String matrix;
	private int wordSize;
	private String effectiveLengthOfTheDatabase;
	private int numberOfBestHitsFromARegionToKeep;
	private int multipleHitOrsingleHit;
	private String  stringEffectiveLengthOfTheSearchSpace;
	private int queryStrandsToSearchAgainstDatabase;
	private boolean produceHTMLOutput;
	private String restrictSearchOfDatabaseToListOfGI;
	private boolean useLowerCaseFilteringOfFASTASequence;
	private String xDropoffValueForUngappedExtensionsInBits;
	private int xDropoffValueForFinalGappedAlignmentInBits;
	private String pSITBlastn;
	private boolean megaBlastSearch;
	private String locationOnQuerySequence;
	private int multipleHitsWindowSize;
	private int frameShiftPenalty;
	private int lengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments;
	private int numberOfConcatenatedQueries;
	private boolean forceUseOfTheLegacyBLASTEngine;
	private String useCompositionBasedScoreAdjustmentsForBlastpOrTblastn;
	private String computeLocallyOptimalSmithWatermanAlignments;
	
	private CommandLine cmd;
	
	private int threadNum=Runtime.getRuntime().availableProcessors();
	
	//other parameters
	
	private MyThreadCount c = new MyThreadCount(0);
	private String tempPath;
	private String blastPath;
	private MakeDir mkr=new MakeDir();
	private String osPlatform;
	private String os3264;
	
	private ArrayList<String> cmdArrayList = new ArrayList<String>();
	
	public ParallelBlast(CommandLine cmd) {
		
		this.cmd = cmd;
		program = cmd.getOptionValue("p"); 
		if(program != null) {
			cmdArrayList.add("-p");
			cmdArrayList.add(program);
		}
		
		database = cmd.getOptionValue("d");
		if(database != null) {
			cmdArrayList.add("-d");
			cmdArrayList.add(database);
		}

		inputFile = cmd.getOptionValue("i");
		if(inputFile != null) {
			//cmdArrayList.add("-i");
			//cmdArrayList.add(inputFile);
		}
		
		String expectationvalueS = cmd.getOptionValue("e");
		if(expectationvalueS != null) {
			cmdArrayList.add("-e");
			cmdArrayList.add(expectationvalueS);
		}
		
		
		if(cmd.getOptionValue("m")!= null) {
			viewOptions=Integer.parseInt(cmd.getOptionValue("m"));
			cmdArrayList.add("-m");
			cmdArrayList.add(String.valueOf(viewOptions));
		}
		
		outPath=cmd.getOptionValue("o");
		if(outPath!= null) {
			//cmdArrayList.add("-o");
			//cmdArrayList.add(outPath);
		}
		
		filterQuerySequence=cmd.getOptionValue("F");
		if(filterQuerySequence!= null) {
			cmdArrayList.add("-F");
			cmdArrayList.add(filterQuerySequence);
		}
		
		
		if(cmd.getOptionValue("G")!= null) {
			costToOpenAGap=Integer.parseInt(cmd.getOptionValue("G"));
			cmdArrayList.add("-G");
			cmdArrayList.add(String.valueOf(costToOpenAGap));
		}
		
		
		if(cmd.getOptionValue("E")!= null) {
			costToExtendAGap=Integer.parseInt(cmd.getOptionValue("E"));
			cmdArrayList.add("-E");
			cmdArrayList.add(String.valueOf(costToExtendAGap));
		}
		
		
		if(cmd.getOptionValue("X")!= null) {
			xDropoffValueForGappedAlignment=Integer.parseInt(cmd.getOptionValue("X"));
			cmdArrayList.add("-X");
			cmdArrayList.add(String.valueOf(xDropoffValueForGappedAlignment));
		}
		
		
		if(cmd.getOptionValue("I")!= null) {
			if(cmd.getOptionValue("I").equalsIgnoreCase("t")){
				showGIInDeflines=true;
			}else if(cmd.getOptionValue("I").equalsIgnoreCase("f")){
				showGIInDeflines=false;
			}else{
				try {
					throw new ParameterError("-I must be T or F");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			cmdArrayList.add("-I");
			if(showGIInDeflines){
				cmdArrayList.add("T");
			}else{
				cmdArrayList.add("F");
			}
		}
		
		if(cmd.getOptionValue("q")!= null) {
			penaltyForANucleotideMismatch=Integer.parseInt(cmd.getOptionValue("q"));
			cmdArrayList.add("-q");
			cmdArrayList.add(String.valueOf(penaltyForANucleotideMismatch));
		}
		
		if(cmd.getOptionValue("r")!= null) {
			rewardForANucleotideMatch=Integer.parseInt(cmd.getOptionValue("r"));
			cmdArrayList.add("-r");
			cmdArrayList.add(String.valueOf(rewardForANucleotideMatch));
		}
		
		if(cmd.getOptionValue("v")!= null) {
			numberOfDatabaseSequencesToShowOneLineDescriptions=Integer.parseInt(cmd.getOptionValue("v"));
			cmdArrayList.add("-v");
			cmdArrayList.add(String.valueOf(numberOfDatabaseSequencesToShowOneLineDescriptions));
		}
		
		if(cmd.getOptionValue("b")!= null) {
			numberOfDatabaseSequenceToShowAlignmentsForB=Integer.parseInt(cmd.getOptionValue("b"));
			cmdArrayList.add("-b");
			cmdArrayList.add(String.valueOf(numberOfDatabaseSequenceToShowAlignmentsForB));
		}
		
		if(cmd.getOptionValue("f")!= null) {
			thresholdForExtendingHits=Integer.parseInt(cmd.getOptionValue("f"));
			cmdArrayList.add("-f");
			cmdArrayList.add(String.valueOf(thresholdForExtendingHits));
		}
		
		if(cmd.getOptionValue("g")!= null) {
			if(cmd.getOptionValue("g").equalsIgnoreCase("t")){
				performGappedAlignment=true;
			}else if(cmd.getOptionValue("g").equalsIgnoreCase("f")){
				performGappedAlignment=false;
			}else{
				try {
					throw new ParameterError("-g must be T or F");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			cmdArrayList.add("-g");
			if(performGappedAlignment){
				cmdArrayList.add("T");
			}else{
				cmdArrayList.add("F");
			}
		}
		
		if(cmd.getOptionValue("Q")!= null) {
			queryGeneticCodeToUse=Integer.parseInt(cmd.getOptionValue("Q"));
			cmdArrayList.add("-Q");
			cmdArrayList.add(String.valueOf(queryGeneticCodeToUse));
		}
		
		if(cmd.getOptionValue("D")!= null) {
			dBGeneticCode=Integer.parseInt(cmd.getOptionValue("D"));
			cmdArrayList.add("-D");
			cmdArrayList.add(String.valueOf(dBGeneticCode));
		}
		
		if(cmd.getOptionValue("a")!= null) {
			numberOfProcessorsToUse=Integer.parseInt(cmd.getOptionValue("a"));
			cmdArrayList.add("-a");
			cmdArrayList.add(String.valueOf(numberOfProcessorsToUse));
		}
		
		seqAlignFile=cmd.getOptionValue("O");
		if(cmd.getOptionValue("O")!= null) {
			cmdArrayList.add("-O");
			cmdArrayList.add(seqAlignFile);
		}
		
		if(cmd.getOptionValue("J")!= null) {
			if(cmd.getOptionValue("J").equalsIgnoreCase("t")){
				believeTheQueryDefline=true;
			}else if(cmd.getOptionValue("J").equalsIgnoreCase("f")){
				believeTheQueryDefline=false;
			}else{
				try {
					throw new ParameterError("-J must be T or F");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			cmdArrayList.add("-J");
			if(believeTheQueryDefline){
				cmdArrayList.add("T");
			}else{
				cmdArrayList.add("F");
			}
		}
		
		matrix=cmd.getOptionValue("M");
		if(cmd.getOptionValue("M")!= null) {
			cmdArrayList.add("-M");
			cmdArrayList.add(matrix);
		}
		
		if(cmd.getOptionValue("W")!= null) {
			wordSize=Integer.parseInt(cmd.getOptionValue("W"));
			cmdArrayList.add("-W");
			cmdArrayList.add(String.valueOf(wordSize));
		}
		
		effectiveLengthOfTheDatabase=cmd.getOptionValue("z");
		if(cmd.getOptionValue("z")!= null) {
			cmdArrayList.add("-z");
			cmdArrayList.add(effectiveLengthOfTheDatabase);
		}
		
		if(cmd.getOptionValue("K")!= null) {
			numberOfBestHitsFromARegionToKeep=Integer.parseInt(cmd.getOptionValue("K"));
			cmdArrayList.add("-K");
			cmdArrayList.add(String.valueOf(numberOfBestHitsFromARegionToKeep));
		}
		
		if(cmd.getOptionValue("P")!= null) {
			multipleHitOrsingleHit=Integer.parseInt(cmd.getOptionValue("P"));
			cmdArrayList.add("-P");
			cmdArrayList.add(String.valueOf(multipleHitOrsingleHit));
		}
		
		stringEffectiveLengthOfTheSearchSpace=cmd.getOptionValue("Y");
		if(cmd.getOptionValue("Y")!= null) {
			cmdArrayList.add("-Y");
			cmdArrayList.add(stringEffectiveLengthOfTheSearchSpace);
		}
		
		if(cmd.getOptionValue("S")!= null) {
			queryStrandsToSearchAgainstDatabase=Integer.parseInt(cmd.getOptionValue("S"));
			cmdArrayList.add("-S");
			cmdArrayList.add(String.valueOf(queryStrandsToSearchAgainstDatabase));
		}
		
		if(cmd.getOptionValue("T")!= null) {
			if(cmd.getOptionValue("T").equalsIgnoreCase("t")){
				produceHTMLOutput=true;
			}else if(cmd.getOptionValue("T").equalsIgnoreCase("f")){
				produceHTMLOutput=false;
			}else{
				try {
					throw new ParameterError("-T must be T or F");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			cmdArrayList.add("-T");
			if(produceHTMLOutput){
				cmdArrayList.add("T");
			}else{
				cmdArrayList.add("F");
			}
		}
		
		restrictSearchOfDatabaseToListOfGI=cmd.getOptionValue("l");
		if(cmd.getOptionValue("l")!= null) {
			cmdArrayList.add("-l");
			cmdArrayList.add(restrictSearchOfDatabaseToListOfGI);
		}
		
		if(cmd.getOptionValue("U")!= null) {
			if(cmd.getOptionValue("U").equalsIgnoreCase("t")){
				useLowerCaseFilteringOfFASTASequence=true;
			}else if(cmd.getOptionValue("U").equalsIgnoreCase("f")){
				useLowerCaseFilteringOfFASTASequence=false;
			}else{
				try {
					throw new ParameterError("-U must be T or F");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			cmdArrayList.add("-U");
			if(useLowerCaseFilteringOfFASTASequence){
				cmdArrayList.add("T");
			}else{
				cmdArrayList.add("F");
			}
		}
		
		xDropoffValueForUngappedExtensionsInBits=cmd.getOptionValue("y");
		if(cmd.getOptionValue("y")!= null) {
			cmdArrayList.add("-y");
			cmdArrayList.add(xDropoffValueForUngappedExtensionsInBits);
		}
		
		if(cmd.getOptionValue("Z")!= null) {
			xDropoffValueForFinalGappedAlignmentInBits=Integer.parseInt(cmd.getOptionValue("Z"));
			cmdArrayList.add("-Z");
			cmdArrayList.add(String.valueOf(xDropoffValueForFinalGappedAlignmentInBits));
		}
		
		pSITBlastn=cmd.getOptionValue("R");
		if(cmd.getOptionValue("R")!= null) {
			cmdArrayList.add("-R");
			cmdArrayList.add(pSITBlastn);
		}
		
		if(cmd.getOptionValue("n")!= null) {
			if(cmd.getOptionValue("n").equalsIgnoreCase("t")){
				megaBlastSearch=true;
			}else if(cmd.getOptionValue("n").equalsIgnoreCase("f")){
				megaBlastSearch=false;
			}else{
				try {
					throw new ParameterError("-n must be T or F");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			cmdArrayList.add("-n");
			if(megaBlastSearch){
				cmdArrayList.add("T");
			}else{
				cmdArrayList.add("F");
			}
		}
		
		locationOnQuerySequence=cmd.getOptionValue("L");
		if(cmd.getOptionValue("L")!= null) {
			cmdArrayList.add("-L");
			cmdArrayList.add(locationOnQuerySequence);
		}
		
		if(cmd.getOptionValue("A")!= null) {
			multipleHitsWindowSize=Integer.parseInt(cmd.getOptionValue("A"));
			cmdArrayList.add("-A");
			cmdArrayList.add(String.valueOf(multipleHitsWindowSize));
		}
		
		if(cmd.getOptionValue("w")!= null) {
			frameShiftPenalty=Integer.parseInt(cmd.getOptionValue("w"));
			cmdArrayList.add("-w");
			cmdArrayList.add(String.valueOf(frameShiftPenalty));
		}
		
		if(cmd.getOptionValue("t")!= null) {
			lengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments=Integer.parseInt(cmd.getOptionValue("t"));
			cmdArrayList.add("-t");
			cmdArrayList.add(String.valueOf(lengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments));
		}
		
		if(cmd.getOptionValue("B")!= null) {
			numberOfConcatenatedQueries=Integer.parseInt(cmd.getOptionValue("B"));
			cmdArrayList.add("-B");
			cmdArrayList.add(String.valueOf(numberOfConcatenatedQueries));
		}
		
		if(cmd.getOptionValue("V")!= null) {
			if(cmd.getOptionValue("V").equalsIgnoreCase("t")){
				forceUseOfTheLegacyBLASTEngine=true;
			}else if(cmd.getOptionValue("V").equalsIgnoreCase("f")){
				forceUseOfTheLegacyBLASTEngine=false;
			}else{
				try {
					throw new ParameterError("-n must be T or F");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			cmdArrayList.add("-V");
			if(forceUseOfTheLegacyBLASTEngine){
				cmdArrayList.add("T");
			}else{
				cmdArrayList.add("F");
			}
		}
		
		useCompositionBasedScoreAdjustmentsForBlastpOrTblastn=cmd.getOptionValue("C");
		if(cmd.getOptionValue("C")!= null) {
			cmdArrayList.add("-C");
			cmdArrayList.add(useCompositionBasedScoreAdjustmentsForBlastpOrTblastn);
		}
		
		computeLocallyOptimalSmithWatermanAlignments=cmd.getOptionValue("s");
		if(cmd.getOptionValue("s")!= null) {
			cmdArrayList.add("-s");
			cmdArrayList.add(computeLocallyOptimalSmithWatermanAlignments);
		}
		
		this.init();
		this.doBlast();
	}
	
	public void init(){
		File outPutFile = new File(outPath);
		mkr.makeParentDirs(outPath);
		try {
			if (outPutFile.isFile() && outPutFile.exists()) {  
				outPutFile.delete();
			}
			outPutFile.createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
		this.setPathOfBlastAndTempPath();
		mkr.makeDirs(tempPath);
		
	}

	public synchronized void doBlast(){
		
		InputStreamReader read;
		try {
			File outFile=null;
			
			FileOutputStream fileOutputStream=null;
			OutputStreamWriter outputStreamWriter=null;
			BufferedWriter bufferedWriter=null;
			String uuid=null;
			
			read = new InputStreamReader(new FileInputStream(inputFile));
			BufferedReader bufferedReader = new BufferedReader(read);
			String line = null;
			int sequenceNumber=0;
			while((line=bufferedReader.readLine()) != null){
				Pattern fastaStartPattern = Pattern.compile("^>");
				Matcher fastaStartMatcher = fastaStartPattern.matcher(line);
				if(fastaStartMatcher.find()){
					sequenceNumber++;
					if(null != uuid){
						bufferedWriter.close();
						boolean isThisThreadUnrun=true;
						while(isThisThreadUnrun){
							if(c.getCount()<threadNum){
								BlastRun blastRun=new BlastRun(c, uuid, blastPath, tempPath, outPath, cmdArrayList);
								c.plusOne();
								blastRun.start();
			                    isThisThreadUnrun=false;
			                    break;
							}else{
								Thread.sleep(10);
							}
						}

					}
					uuid = UUID.randomUUID().toString();
					outFile=new File(tempPath + File.separator + uuid);
					outFile.createNewFile();
					fileOutputStream=new FileOutputStream(outFile);
					outputStreamWriter=new OutputStreamWriter(fileOutputStream);
					bufferedWriter=new BufferedWriter(outputStreamWriter);
				}
				bufferedWriter.write(line);
				bufferedWriter.newLine();
				bufferedWriter.flush();
			}
			bufferedWriter.close();
			read.close();
			boolean isThisThreadUnrun=true;
			while(isThisThreadUnrun){
				if(c.getCount()<threadNum){
					BlastRun blastRun=new BlastRun(c, uuid,blastPath, tempPath, outPath, cmdArrayList);
					c.plusOne();
					blastRun.start();
                    isThisThreadUnrun=false;
                    break;
				}else{
					Thread.sleep(10);
				}
			}
			while(true){// wait all thread over
	            if(!c.hasNext()) break; 
	        }

			
			if(sequenceNumber<1){
				throw new NotFastaFormat("The input file is not in FASTA format, please have a check!");
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NotFastaFormat e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
	
	private void setPathOfBlastAndTempPath(){
		File directory = new File("");
		Properties props=System.getProperties();
		String osName = props.getProperty("os.name");
		String osArch = props.getProperty("os.arch");
		if(osName.contains("Windows")){
			this.setOsPlatform(Constant.WINDOWSOS);
		}else{
			try {
				throw new UnknownOs("Sorry, your OS is unknown or unsupported, please contact us:" + Constant.MYEMAIL);
			} catch (UnknownOs e) {
				e.printStackTrace();
			}
			System.exit(1);
		}
		
		if(osArch.equalsIgnoreCase("x86")){
			this.setOs3264(Constant.X86OS);
		}else{
			try {
				throw new UnknownOs("Sorry, your OS Arch is unknown or unsupported, please contact us:" + Constant.MYEMAIL);
			} catch (UnknownOs e) {
				e.printStackTrace();
			}
			System.exit(1);
		}
		
		if(this.getOsPlatform().equals(Constant.WINDOWSOS) && this.getOs3264().equals(Constant.X86OS)){
			this.setBlastPath(directory.getAbsolutePath() + File.separator + Constant.WINDOWS32BLASTPATH);
		}
		tempPath = directory.getAbsolutePath() + File.separator + "temp" + File.separator ;
	}
	
	
	public static void main(String[] args) throws ParseException {
		Options options = new Options();
		options.addOption("h", false, "Lists short help");
		options.addOption("p", true, "Program Name [String]");
		options.addOption("d", true, "Database [String]\n"
				+ "\tdefault = nr");
		options.addOption("i", true, "Query File [File In]\n"
				+ "\tdefault = stdin");
		options.addOption("e", true, "Expectation value (E) [Real]\n"
				+ "\tdefault = 10.0");
		
		options.addOption("m", true, "alignment view options:\n"
				+ "0 = pairwise,\n"
				+ "1 = query-anchored showing identities,\n"
				+ "2 = query-anchored no identities,\n"
				+ "3 = flat query-anchored, show identities,\n"
				+ "4 = flat query-anchored, no identities,\n"
				+ "5 = query-anchored no identities and blunt ends,\n"
				+ "6 = flat query-anchored, no identities and blunt ends,\n"
				+ "7 = XML Blast output,\n"
				+ "8 = tabular,\n"
				+ "9 tabular with comment lines\n"
				+ "10 ASN, text\n"
				+ "11 ASN, binary [Integer]\n"
				+ "\tdefault = 0\n"
				+ "\trange from 0 to 11");
				
		options.addOption("o", true, "BLAST report Output File [File Out]  Optional\n"
				+ "\tdefault = stdout");
		options.addOption("F", true, "Filter query sequence (DUST with blastn, SEG with others) [String]\n"
				+ "\tdefault = T");
		options.addOption("G", true, "Cost to open a gap (-1 invokes default behavior) [Integer]\n"
				+ "\tdefault = -1");
		options.addOption("E", true, "Cost to extend a gap (-1 invokes default behavior) [Integer]\n"
				+ "\tdefault = -1");
		options.addOption("X", true, "X dropoff value for gapped alignment (in bits) (zero invokes default behavior)\n"
				+ "\tblastn 30, megablast 20, tblastx 0, all others 15 [Integer]\n"
				+ "\tdefault = 0");
		options.addOption("I", true, "Show GI's in deflines [T/F]\n"
				+ "\tdefault = F");
		options.addOption("q", true, "Penalty for a nucleotide mismatch (blastn only) [Integer]\n"
				+ "\tdefault = -3");
		options.addOption("r", true, "Reward for a nucleotide match (blastn only) [Integer]\n"
				+ "\tdefault = 1");
		options.addOption("v", true, "Number of database sequences to show one-line descriptions for (V) [Integer]\n"
				+ "\tdefault = 500");
		options.addOption("b", true, "Number of database sequence to show alignments for (B) [Integer]\n"
				+ "\tdefault = 250");
		
		options.addOption("f", true, "Threshold for extending hits, default if zero\n"
				+ "\tblastp 11, blastn 0, blastx 12, tblastn 13\n"
				+ "\ttblastx 13, megablast 0 [Real]\n"
				+ "\tdefault = 0");
		options.addOption("g", true, "Perform gapped alignment (not available with tblastx) [T/F]\n"
				+ "\tdefault = T");
		options.addOption("Q", true, "Query Genetic code to use [Integer]\n"
				+ "\tdefault = 1");
		options.addOption("D", true, "DB Genetic code (for tblast[nx] only) [Integer]\n"
				+ "\tdefault = 1");
		options.addOption("a", true, "Number of processors to use [Integer]\n"
				+ "\tdefault = 1");
		options.addOption("O", true, "SeqAlign file [File Out]  Optional\n");
		options.addOption("J", true, "Believe the query defline [T/F]\n"
				+ "\tdefault = F");
		options.addOption("M", true, "Matrix [String]\n"
				+ "\tdefault = BLOSUM62");
		options.addOption("W", true, "Word size, default if zero (blastn 11, megablast 28, all others 3) [Integer]\n"
				+ "\tdefault = 0");
		options.addOption("Z", true, "Effective length of the database (use zero for the real size) [Real]\n"
				+ "\tdefault = 0");
		options.addOption("K", true, "Number of best hits from a region to keep. Off by default.\n"
				+ "If used a value of 100 is recommended.  Very high values of -v or -b is also suggested [Integer]\n"
				+ "\tdefault = 0");
		options.addOption("P", true, "0 for multiple hit, 1 for single hit (does not apply to blastn) [Integer]\n"
				+ "\tdefault = 0");
		options.addOption("Y", true, "Effective length of the search space (use zero for the real size) [Real]\n"
				+ "\tdefault = 0");
		options.addOption("S", true, "Query strands to search against database (for blast[nx], and tblastx)\n"
				+ "\t3 is both, 1 is top, 2 is bottom [Integer]\n"
				+ "\tdefault = 3");
		options.addOption("T", true, "Produce HTML output [T/F]\n"
				+ "\tdefault = F");
		options.addOption("l", true, "Restrict search of database to list of GI's [String]  Optional\n");
		options.addOption("U", true, "\nUse lower case filtering of FASTA sequence [T/F]  Optional");
		options.addOption("y", true, "X dropoff value for ungapped extensions in bits (0.0 invokes default behavior)\n"
				+ "\tblastn 20, megablast 10, all others 7 [Real]\n"
				+ "\tdefault = 0.0");
		options.addOption("Z", true, "X dropoff value for final gapped alignment in bits (0.0 invokes default behavior)\n"
				+ "\tblastn/megablast 100, tblastx 0, all others 25 [Integer]\n"
				+ "\tdefault = 0");
		options.addOption("R", true, "PSI-TBLASTN checkpoint file [File In]  Optional\n");
		options.addOption("n", true, "MegaBlast search [T/F]\n"
				+ "\tdefault = F");
		options.addOption("L", true, "Location on query sequence [String]  Optional\n");
		options.addOption("A", true, "Multiple Hits window size, default if zero (blastn/megablast 0, all others 40 [Integer]\n"
				+ "\tdefault = 0");
		options.addOption("w", true, "Frame shift penalty (OOF algorithm for blastx) [Integer]\n"
				+ "\tdefault = 0");
		options.addOption("t", true, "Length of the largest intron allowed in a translated nucleotide sequence when linking multiple distinct alignments. (0 invokes default behavior; a negative value disables linking.) [Integer]\n"
				+ "default = 0\t");
		options.addOption("B", true, "Number of concatenated queries, for blastn and tblastn [Integer]  Optional\n"
				+ "\tdefault = 0");
		options.addOption("V", true, "Force use of the legacy BLAST engine [T/F]  Optional\n"
				+ "\tdefault = F");
		options.addOption("C", true, "Use composition-based score adjustments for blastp or tblastn:\n"
				+ "\tAs first character:\n"
				+ "\tD or d: default (equivalent to T)\n"
				+ "\t0 or F or f: no composition-based statistics\n"
				+ "\t2 or T or t: Composition-based score adjustments as in Bioinformatics 21:902-911,\n"
				+ "\t1: Composition-based statistics as in NAR 29:2994-3005, 2001\n"
				+ "\t\t2005, conditioned on sequence properties\n"
				+ "\t3: Composition-based score adjustment as in Bioinformatics 21:902-911,\n"
				+ "\t\t2005, unconditionally\n"
				+ "\tFor programs other than tblastn, must either be absent or be D, F or 0.\n"
				+ "\t\tAs second character, if first character is equivalent to 1, 2, or 3:\n"
				+ "\tU or u: unified p-value combining alignment p-value and compositional p-value in round 1 only\n"
				+ "\t[String]\n"
				+ "\tdefault = D");
		options.addOption("s", true, "Compute locally optimal Smith-Waterman alignments (This option is only\n"
				+ "\tavailable for gapped tblastn.) [T/F]\n"
				+ "\tdefault = F");
		
		CommandLineParser parser = new PosixParser(); 
		CommandLine cmd = parser.parse(options, args);
		HelpFormatter hf = new HelpFormatter(); 
		if (cmd.getOptions().length > 0) { 
			if(cmd.hasOption("h")) {
				hf.printHelp("Options", options);
				System.exit(0);
			}
			
			
			
			String programCli = cmd.getOptionValue("p");
			if(programCli == null) {
				try {
					throw new ParameterError("-p must be difinded");
				} catch (ParameterError e) {
					e.printStackTrace();
				}
				System.exit(0);
			}
			
			ParallelBlast pb=new ParallelBlast(cmd);
			//pb.doBlast();
		}else{
			hf.printHelp("Options", options);
		}
	}

	public String getInputFile() {
		return inputFile;
	}

	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}

	public String getTempPath() {
		return tempPath;
	}

	public void setTempPath(String tempPath) {
		this.tempPath = tempPath;
	}

	public String getOutPath() {
		return outPath;
	}

	public void setOutPath(String outPath) {
		this.outPath = outPath;
	}

	public String getDatabase() {
		return database;
	}

	public void setDatabase(String database) {
		this.database = database;
	}

	public MyThreadCount getC() {
		return c;
	}

	public void setC(MyThreadCount c) {
		this.c = c;
	}

	public int getThreadNum() {
		return threadNum;
	}

	public void setThreadNum(int threadNum) {
		this.threadNum = threadNum;
	}

	public String getBlastPath() {
		return blastPath;
	}

	public void setBlastPath(String blastPath) {
		this.blastPath = blastPath;
	}

	public MakeDir getMkr() {
		return mkr;
	}

	public void setMkr(MakeDir mkr) {
		this.mkr = mkr;
	}

	public String getOsPlatform() {
		return osPlatform;
	}

	public void setOsPlatform(String osPlatform) {
		this.osPlatform = osPlatform;
	}

	public String getOs3264() {
		return os3264;
	}

	public void setOs3264(String os3264) {
		this.os3264 = os3264;
	}

	public String getProgram() {
		return program;
	}

	public void setProgram(String program) {
		this.program = program;
	}

	public Double getExpectationvalue() {
		return expectationvalue;
	}

	public void setExpectationvalue(Double expectationvalue) {
		this.expectationvalue = expectationvalue;
	}

	public int getViewOptions() {
		return viewOptions;
	}

	public void setViewOptions(int viewOptions) {
		this.viewOptions = viewOptions;
	}

	public String getFilterQuerySequence() {
		return filterQuerySequence;
	}

	public void setFilterQuerySequence(String filterQuerySequence) {
		this.filterQuerySequence = filterQuerySequence;
	}

	public int getCostToOpenAGap() {
		return costToOpenAGap;
	}

	public void setCostToOpenAGap(int costToOpenAGap) {
		this.costToOpenAGap = costToOpenAGap;
	}

	public int getCostToExtendAGap() {
		return costToExtendAGap;
	}

	public void setCostToExtendAGap(int costToExtendAGap) {
		this.costToExtendAGap = costToExtendAGap;
	}

	public int getxDropoffValueForGappedAlignment() {
		return xDropoffValueForGappedAlignment;
	}

	public void setxDropoffValueForGappedAlignment(
			int xDropoffValueForGappedAlignment) {
		this.xDropoffValueForGappedAlignment = xDropoffValueForGappedAlignment;
	}

	

	public int getPenaltyForANucleotideMismatch() {
		return penaltyForANucleotideMismatch;
	}

	public void setPenaltyForANucleotideMismatch(int penaltyForANucleotideMismatch) {
		this.penaltyForANucleotideMismatch = penaltyForANucleotideMismatch;
	}

	public int getRewardForANucleotideMatch() {
		return rewardForANucleotideMatch;
	}

	public void setRewardForANucleotideMatch(int rewardForANucleotideMatch) {
		this.rewardForANucleotideMatch = rewardForANucleotideMatch;
	}

	public int getNumberOfDatabaseSequencesToShowOneLineDescriptions() {
		return numberOfDatabaseSequencesToShowOneLineDescriptions;
	}

	public void setNumberOfDatabaseSequencesToShowOneLineDescriptions(
			int numberOfDatabaseSequencesToShowOneLineDescriptions) {
		this.numberOfDatabaseSequencesToShowOneLineDescriptions = numberOfDatabaseSequencesToShowOneLineDescriptions;
	}

	public int getNumberOfDatabaseSequenceToShowAlignmentsForB() {
		return numberOfDatabaseSequenceToShowAlignmentsForB;
	}

	public void setNumberOfDatabaseSequenceToShowAlignmentsForB(
			int numberOfDatabaseSequenceToShowAlignmentsForB) {
		this.numberOfDatabaseSequenceToShowAlignmentsForB = numberOfDatabaseSequenceToShowAlignmentsForB;
	}

	public int getThresholdForExtendingHits() {
		return thresholdForExtendingHits;
	}

	public void setThresholdForExtendingHits(int thresholdForExtendingHits) {
		this.thresholdForExtendingHits = thresholdForExtendingHits;
	}

	public boolean isPerformGappedAlignment() {
		return performGappedAlignment;
	}

	public void setPerformGappedAlignment(boolean performGappedAlignment) {
		this.performGappedAlignment = performGappedAlignment;
	}

	public int getQueryGeneticCodeToUse() {
		return queryGeneticCodeToUse;
	}

	public void setQueryGeneticCodeToUse(int queryGeneticCodeToUse) {
		this.queryGeneticCodeToUse = queryGeneticCodeToUse;
	}

	public int getdBGeneticCode() {
		return dBGeneticCode;
	}

	public void setdBGeneticCode(int dBGeneticCode) {
		this.dBGeneticCode = dBGeneticCode;
	}

	public int getNumberOfProcessorsToUse() {
		return numberOfProcessorsToUse;
	}

	public void setNumberOfProcessorsToUse(int numberOfProcessorsToUse) {
		this.numberOfProcessorsToUse = numberOfProcessorsToUse;
	}

	public String getSeqAlignFile() {
		return seqAlignFile;
	}

	public void setSeqAlignFile(String seqAlignFile) {
		this.seqAlignFile = seqAlignFile;
	}

	public boolean isBelieveTheQueryDefline() {
		return believeTheQueryDefline;
	}

	public void setBelieveTheQueryDefline(boolean believeTheQueryDefline) {
		this.believeTheQueryDefline = believeTheQueryDefline;
	}

	public String getMatrix() {
		return matrix;
	}

	public void setMatrix(String matrix) {
		this.matrix = matrix;
	}

	public int getWordSize() {
		return wordSize;
	}

	public void setWordSize(int wordSize) {
		this.wordSize = wordSize;
	}

	public String getEffectiveLengthOfTheDatabase() {
		return effectiveLengthOfTheDatabase;
	}

	public void setEffectiveLengthOfTheDatabase(String effectiveLengthOfTheDatabase) {
		this.effectiveLengthOfTheDatabase = effectiveLengthOfTheDatabase;
	}

	public int getNumberOfBestHitsFromARegionToKeep() {
		return numberOfBestHitsFromARegionToKeep;
	}

	public void setNumberOfBestHitsFromARegionToKeep(
			int numberOfBestHitsFromARegionToKeep) {
		this.numberOfBestHitsFromARegionToKeep = numberOfBestHitsFromARegionToKeep;
	}

	public int getMultipleHitOrsingleHit() {
		return multipleHitOrsingleHit;
	}

	public void setMultipleHitOrsingleHit(int multipleHitOrsingleHit) {
		this.multipleHitOrsingleHit = multipleHitOrsingleHit;
	}

	public String getStringEffectiveLengthOfTheSearchSpace() {
		return stringEffectiveLengthOfTheSearchSpace;
	}

	public void setStringEffectiveLengthOfTheSearchSpace(
			String stringEffectiveLengthOfTheSearchSpace) {
		this.stringEffectiveLengthOfTheSearchSpace = stringEffectiveLengthOfTheSearchSpace;
	}

	public int getQueryStrandsToSearchAgainstDatabase() {
		return queryStrandsToSearchAgainstDatabase;
	}

	public void setQueryStrandsToSearchAgainstDatabase(
			int queryStrandsToSearchAgainstDatabase) {
		this.queryStrandsToSearchAgainstDatabase = queryStrandsToSearchAgainstDatabase;
	}

	public boolean isProduceHTMLOutput() {
		return produceHTMLOutput;
	}

	public void setProduceHTMLOutput(boolean produceHTMLOutput) {
		this.produceHTMLOutput = produceHTMLOutput;
	}

	public String getRestrictSearchOfDatabaseToListOfGI() {
		return restrictSearchOfDatabaseToListOfGI;
	}

	public void setRestrictSearchOfDatabaseToListOfGI(
			String restrictSearchOfDatabaseToListOfGI) {
		this.restrictSearchOfDatabaseToListOfGI = restrictSearchOfDatabaseToListOfGI;
	}

	public boolean isUseLowerCaseFilteringOfFASTASequence() {
		return useLowerCaseFilteringOfFASTASequence;
	}

	public void setUseLowerCaseFilteringOfFASTASequence(
			boolean useLowerCaseFilteringOfFASTASequence) {
		this.useLowerCaseFilteringOfFASTASequence = useLowerCaseFilteringOfFASTASequence;
	}

	public String getxDropoffValueForUngappedExtensionsInBits() {
		return xDropoffValueForUngappedExtensionsInBits;
	}

	public void setxDropoffValueForUngappedExtensionsInBits(
			String xDropoffValueForUngappedExtensionsInBits) {
		this.xDropoffValueForUngappedExtensionsInBits = xDropoffValueForUngappedExtensionsInBits;
	}

	public int getxDropoffValueForFinalGappedAlignmentInBits() {
		return xDropoffValueForFinalGappedAlignmentInBits;
	}

	public void setxDropoffValueForFinalGappedAlignmentInBits(
			int xDropoffValueForFinalGappedAlignmentInBits) {
		this.xDropoffValueForFinalGappedAlignmentInBits = xDropoffValueForFinalGappedAlignmentInBits;
	}

	public String getpSITBlastn() {
		return pSITBlastn;
	}

	public void setpSITBlastn(String pSITBlastn) {
		this.pSITBlastn = pSITBlastn;
	}

	public boolean isMegaBlastSearch() {
		return megaBlastSearch;
	}

	public void setMegaBlastSearch(boolean megaBlastSearch) {
		this.megaBlastSearch = megaBlastSearch;
	}

	public String getLocationOnQuerySequence() {
		return locationOnQuerySequence;
	}

	public void setLocationOnQuerySequence(String locationOnQuerySequence) {
		this.locationOnQuerySequence = locationOnQuerySequence;
	}


	public int getFrameShiftPenalty() {
		return frameShiftPenalty;
	}

	public void setFrameShiftPenalty(int frameShiftPenalty) {
		this.frameShiftPenalty = frameShiftPenalty;
	}

	public int getLengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments() {
		return lengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments;
	}

	public void setLengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments(
			int lengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments) {
		this.lengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments = lengthOfTheLargestIntronAllowedInATranslatedNucleotideSequenceWhenLinkingMultipleDistinctAlignments;
	}

	public int getNumberOfConcatenatedQueries() {
		return numberOfConcatenatedQueries;
	}

	public void setNumberOfConcatenatedQueries(int numberOfConcatenatedQueries) {
		this.numberOfConcatenatedQueries = numberOfConcatenatedQueries;
	}

	public boolean isForceUseOfTheLegacyBLASTEngine() {
		return forceUseOfTheLegacyBLASTEngine;
	}

	public void setForceUseOfTheLegacyBLASTEngine(
			boolean forceUseOfTheLegacyBLASTEngine) {
		this.forceUseOfTheLegacyBLASTEngine = forceUseOfTheLegacyBLASTEngine;
	}

	public String getUseCompositionBasedScoreAdjustmentsForBlastpOrTblastn() {
		return useCompositionBasedScoreAdjustmentsForBlastpOrTblastn;
	}

	public void setUseCompositionBasedScoreAdjustmentsForBlastpOrTblastn(
			String useCompositionBasedScoreAdjustmentsForBlastpOrTblastn) {
		this.useCompositionBasedScoreAdjustmentsForBlastpOrTblastn = useCompositionBasedScoreAdjustmentsForBlastpOrTblastn;
	}

	public String getComputeLocallyOptimalSmithWatermanAlignments() {
		return computeLocallyOptimalSmithWatermanAlignments;
	}

	public void setComputeLocallyOptimalSmithWatermanAlignments(
			String computeLocallyOptimalSmithWatermanAlignments) {
		this.computeLocallyOptimalSmithWatermanAlignments = computeLocallyOptimalSmithWatermanAlignments;
	}

	public ArrayList<String> getCmdArrayList() {
		return cmdArrayList;
	}

	public void setCmdArrayList(ArrayList<String> cmdArrayList) {
		this.cmdArrayList = cmdArrayList;
	}

	public boolean isShowGIInDeflines() {
		return showGIInDeflines;
	}

	public void setShowGIInDeflines(boolean showGIInDeflines) {
		this.showGIInDeflines = showGIInDeflines;
	}

	public int getMultipleHitsWindowSize() {
		return multipleHitsWindowSize;
	}

	public void setMultipleHitsWindowSize(int multipleHitsWindowSize) {
		this.multipleHitsWindowSize = multipleHitsWindowSize;
	}

	public CommandLine getCmd() {
		return cmd;
	}

	public void setCmd(CommandLine cmd) {
		this.cmd = cmd;
	}
}
