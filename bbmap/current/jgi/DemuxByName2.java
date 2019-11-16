package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map.Entry;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.BufferedMultiCros;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.MultiCros2;
import stream.MultiCros3;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;


/**
 * This class is designed to handle very large numbers of output files
 * with a small number of file handles.  An example use case is demultiplexing
 * Illumina Novaseq runs.
 * 
 * @author Brian Bushnell
 * @date May 1, 2019
 *
 */
public class DemuxByName2 {

	/** Code entrance from the command line */
	public static void main(String[] args){
		
		//Capture and subsequently restore static variables in case this is called by another class. 
		final int oldCap=Shared.numBuffers(), oldZipThreads=ReadWrite.MAX_ZIP_THREADS, oldZl=ReadWrite.ZIPLEVEL;
		final boolean oldPigz=ReadWrite.USE_PIGZ, oldUnpigz=ReadWrite.USE_UNPIGZ;
		
		Timer t=new Timer();
		DemuxByName2 x=new DemuxByName2(args);
		x.process(t);
		
		Shared.setBuffers(oldCap);
		ReadWrite.ZIPLEVEL=oldZl;
		ReadWrite.USE_PIGZ=oldPigz;
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		ReadWrite.MAX_ZIP_THREADS=oldZipThreads;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public DemuxByName2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		//Set some static variables
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=true;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		//Reduce default compression level to increase speed
		//This automatically gets bumped up to 4 if pigz is detected
		ReadWrite.ZIPLEVEL=2;
		
		Parser parser=new Parser();
		parser.overwrite=overwrite;
		
		for(int i=0; i<args.length; i++){//Parsing loop
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//Flag was captured by the parser; do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				BufferedMultiCros.verbose=verbose;
			}else if(a.equals("names") || a.equals("name") || a.equals("affixes")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						names.put(s, s);
					}
				}
			}else if(a.equals("threaded")){
				threaded=Tools.parseBoolean(b);
			}else if(a.equalsIgnoreCase("length") || a.equalsIgnoreCase("len") || a.equalsIgnoreCase("affixlength") || a.equalsIgnoreCase("affixlen")){
				fixedAffixLength=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("stats") || a.equalsIgnoreCase("report") || a.equalsIgnoreCase("results")){
				stats=b;
			}else if(a.equalsIgnoreCase("minreadstodump") || a.equalsIgnoreCase("minreads")){
				minReadsToDump=Tools.parseKMG(b);
			}else if(a.equalsIgnoreCase("mcrostype")){
				mcrosType=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("streams") || a.equalsIgnoreCase("maxstreams")){
				maxStreams=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("readsperbuffer") || a.equalsIgnoreCase("rpb")){
				readsPerBuffer=Tools.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("bytesPerBuffer") || a.equalsIgnoreCase("bpb")){
				bytesPerBuffer=Tools.parseIntKMG(b);
			}else if(a.equals("prefixmode") || a.equals("prefix") || a.equals("pm")){
				prefixMode=Tools.parseBoolean(b);
			}else if(a.equals("suffixmode") || a.equals("suffix") || a.equals("sm")){
				prefixMode=!Tools.parseBoolean(b);
			}else if(a.equals("column")){
				column=Integer.parseInt(b);
				assert(column>0 || column==-1) : "Column is 1-based; must be 1+ or else -1 to disable.";
				column--;
			}else if(a.equals("hdist") || a.equals("hamming") || a.equals("hammingdistance")){
				hdist=Integer.parseInt(b);
			}else if(a.equals("substringmode") || a.equals("substring")){
				substringMode=Tools.parseBoolean(b);
			}else if(a.equals("outu") || a.equals("outu1")){
				outu1=b;
			}else if(a.equals("outu2")){
				outu2=b;
			}else if(a.equals("pattern")){
				parser.out1=b;
			}else if(a.equals("perheader") || a.equals("persequence") || a.equals("pername")){
				perheader=Tools.parseBoolean(b);
			}else if(a.equals("delimiter")){
				if(b==null){delimiter=null;}
				
				//Convenience characters
				else if(b.equalsIgnoreCase("space")){
					delimiter=" ";
				}else if(b.equalsIgnoreCase("tab")){
					delimiter="\t";
				}else if(b.equalsIgnoreCase("whitespace")){
					delimiter="\\s+";
				}else if(b.equalsIgnoreCase("pound")){
					delimiter="#";
				}else if(b.equalsIgnoreCase("greaterthan")){
					delimiter=">";
				}else if(b.equalsIgnoreCase("lessthan")){
					delimiter="<";
				}else if(b.equalsIgnoreCase("equals")){
					delimiter="=";
				}else if(b.equalsIgnoreCase("colon")){
					delimiter=":";
				}else if(b.equalsIgnoreCase("semicolon")){
					delimiter=";";
				}else if(b.equalsIgnoreCase("bang")){
					delimiter="!";
				}else if(b.equalsIgnoreCase("and") || b.equalsIgnoreCase("ampersand")){
					delimiter="&";
				}else if(b.equalsIgnoreCase("quote") || b.equalsIgnoreCase("doublequote")){
					delimiter="\"";
				}else if(b.equalsIgnoreCase("singlequote") || b.equalsIgnoreCase("apostrophe")){
					delimiter="'";
				}
				
				//Java meta characters
				else if(b.equalsIgnoreCase("backslash")){
					delimiter="\\\\";
				}else if(b.equalsIgnoreCase("hat") || b.equalsIgnoreCase("caret")){
					delimiter="\\^";
				}else if(b.equalsIgnoreCase("dollar")){
					delimiter="\\$";
				}else if(b.equalsIgnoreCase("dot")){
					delimiter="\\.";
				}else if(b.equalsIgnoreCase("pipe") || b.equalsIgnoreCase("or")){
					delimiter="\\|";
				}else if(b.equalsIgnoreCase("questionmark")){
					delimiter="\\?";
				}else if(b.equalsIgnoreCase("star") || b.equalsIgnoreCase("asterisk")){
					delimiter="\\*";
				}else if(b.equalsIgnoreCase("plus")){
					delimiter="\\+";
				}else if(b.equalsIgnoreCase("openparen")){
					delimiter="\\(";
				}else if(b.equalsIgnoreCase("closeparen")){
					delimiter="\\)";
				}else if(b.equalsIgnoreCase("opensquare")){
					delimiter="\\[";
				}else if(b.equalsIgnoreCase("opencurly")){
					delimiter="\\{";
				}
				
				else{
					delimiter=b;
				}
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(perheader){
			//Do nothing
		}else{
			
			{//Process names
				String[] x=names.keySet().toArray(new String[names.size()]);
				names.clear();
				for(String s : x){//Process either literal names or files of names
					File f=new File(s);
					if(f.exists() && f.isFile()){
						TextFile tf=new TextFile(s);
						String[] lines=tf.toStringLines();
						for(String s2 : lines){
							names.put(s2, s2);
						}
					}else{
						names.put(s, s);
					}
				}
			}
			
			{//Find the lengths of all names to fill the affixLengths array
				BitSet bs=new BitSet();
				if(fixedAffixLength>0){
					bs.set(fixedAffixLength);
				}
				for(String s : names.keySet()){
					bs.set(s.length());
				}
				affixLengths=new int[bs.cardinality()];
				for(int i=0, bit=-1; i<affixLengths.length; i++){
					bit=bs.nextSetBit(bit+1);
					affixLengths[i]=bit;
				}
				Arrays.sort(affixLengths);
				Tools.reverseInPlace(affixLengths);
			}
			
			assert((affixLengths.length>0 && affixLengths[0]>0) || delimiter!=null || perheader) : 
				"Must include at least one name, an affix length, or a delimiter.";
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}

		//Validate patterns
		assert(out1==null || out1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout1==null || qfout1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout2==null || qfout2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		
		//Perform # replacement for twin files
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=outu1.replace("#", "2");
			outu1=outu1.replace("#", "1");
		}
		
		//Disable interleaving if in2 is specified
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		{//Perform various file validation
			assert(FastaReadInputStream.settingsOK());

			if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
			if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
				ByteFile.FORCE_MODE_BF2=true;
			}

			if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}

			if(!setInterleaved){
				assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
				if(in2!=null){ //If there are 2 input streams.
					FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}else{ //There is one input stream.
					if(out2!=null){
						FASTQ.FORCE_INTERLEAVED=true;
						FASTQ.TEST_INTERLEAVED=false;
						outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
					}
				}
			}

			if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
			if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}

			if(!Tools.testOutputFiles(overwrite, append, false, outu1, outu2, stats)){
				outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outu1+", "+outu2+", "+stats);
				throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outu1+", "+outu2+", "+stats+"\n");
			}

			ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
			ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);

			if(ffin1!=null && out1!=null && ffin1.samOrBam()){
				String ext=ReadWrite.rawExtension(out1);
				useSharedHeader=FileFormat.isSamOrBam(ext);
			}
			
			//This does not prevent expansion into a matching string
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
		}
		
		if(perheader){
			fixedAffixLength=-1;
			substringMode=false;
			affixLengths=null;
			delimiter=null;
		}
		
		assert(column<0 || delimiter!=null) : "Column may not be set if there is no delimiter.";
		if(delimiter!=null && delimiter.length()==1){
			delimiterChar=delimiter.charAt(0);
		}else if(delimiter!=null && delimiter.length()==2 && delimiter.charAt(0)=='\\'){
			delimiterChar=delimiter.charAt(1);
		}
		
		//Populate names table with mutants for hamming distance
		if(hdist>0 && names.size()>0){
			int mutants=mutate(names, hdist);
		}
		
		if(affixLengths.length>1){fixedAffixLength=-1;}
		
		//Populate name list
		for(Entry<String, String> e : names.entrySet()){
			nameList.add(e.getKey());
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary method.  Processes the data. */
	void process(Timer t){
		
		//Stream for input reads
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
		if(verbose){outstream.println("Started cris");}
		cris.start();
		
		final boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		//Streams for output reads other than unmatched
		final BufferedMultiCros mcros;
		if(out1!=null){
			
			if(mcrosType==2){
				mcros=new MultiCros2(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded);
			}else if(mcrosType==3){//Faster type
				mcros=new MultiCros3(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded, maxStreams);
			}else{
				throw new RuntimeException("Bad mcrosType: "+mcrosType);
			}
			mcros.readsPerBuffer=readsPerBuffer;
			mcros.bytesPerBuffer=bytesPerBuffer;
			mcros.minReadsToDump=minReadsToDump;
			
			if(paired && out2==null && (in1==null || !ffin1.samOrBam())){
				outstream.println("Writing interleaved.");
			}
			
			//Start the mcros thread if necessary
			if(threaded){mcros.start();}
		}else{
			mcros=null;
		}
		
		//Stream for unmatched output reads
		final ConcurrentReadOutputStream rosu;
		if(outu1!=null){
			//Number of output buffers; does not need to be high since access is singlethreaded
			final int buff=4;
			
			FileFormat ffout1=FileFormat.testOutput(outu1, FileFormat.FASTQ, extout, true, overwrite, append, false);
			FileFormat ffout2=(outu2==null ? null : FileFormat.testOutput(outu2, FileFormat.FASTQ, extout, true, overwrite, append, false));
			rosu=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, true);
			rosu.start();
		}else{
			rosu=null;
		}
		
		//Streams are set up, so process the reads
		processInner(cris, mcros, rosu);
		//At this point processing has finished.
		
		//Shut down mcros
		if(mcros!=null){
			mcros.close();
			errorState|=mcros.errorState();
		}
		
		if(minReadsToDump>0 && mcros!=null){
			//Dump the residual reads into the unmatched file
			mcros.dumpResidual(rosu);
			
			//Adjust statistics to reflect that residuals did not get demultiplexed
			readsOut-=mcros.residualReads;
			basesOut-=mcros.residualBases;
		}
		
//		errorState|=ReadStats.writeAll(); //Currently unused
		errorState|=ReadWrite.closeStreams(cris, rosu);
		
		//Report statistics to file
		if(stats!=null){printReport(mcros);}
		
		t.stop();
		
		{//Report statistics to screen
			double rpnano=readsProcessed/(double)(t.elapsed);
			double bpnano=basesProcessed/(double)(t.elapsed);

			outstream.println("Time:               "+t);
			outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", rpnano*1000000));
			outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", bpnano*1000));
			outstream.println("Reads Out:          "+readsOut);
			outstream.println("Bases Out:          "+basesOut);
		}
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	
	/** 
	 * Process all reads.
	 */
	void processInner(final ConcurrentReadInputStream cris, final BufferedMultiCros mcros, final ConcurrentReadOutputStream rosu) {
		
		//Fetch the first list
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		//Ensure pairing seems correct
		if(reads!=null && !reads.isEmpty()){
			Read r=reads.get(0);
			assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
		}
		
		//While there are more reads to process...
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			//List to accumulate unmatched reads
			ArrayList<Read> unmatched=new ArrayList<Read>();
			
			//For each read...
			for(Read r : reads){
				
				//Get the target file identifier
				String name=getValue(r);
				if(name!=null){
					//Set the name so that the mcros will send it to the right place
					r.obj=name;
					readsOut+=r.pairCount();
					basesOut+=r.pairLength();
				}else{
					//Send it to unmatched
					unmatched.add(r);
					readsUnmatched+=r.pairCount();
					basesUnmatched+=r.pairLength();
				}

				readsProcessed+=r.pairCount();
				basesProcessed+=r.pairLength();
			}
			if(rosu!=null){rosu.add(unmatched, ln.id);}
			if(mcros!=null){mcros.add(reads);}

			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		
		if(ln!=null){
			cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
		}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Generate mutants for all names according to the hamming distance */
	private static int mutate(HashMap<String, String> names, int hdist){
		ArrayList<String> list=new ArrayList<String>(names.size());
		list.addAll(names.values());
		int mutants=0;
		final HashSet<String> collisions=new HashSet<String>();
		for(String name : list){
			mutants+=addMutants(name.getBytes(), name, hdist, names, collisions);
		}
		System.err.println("Added "+mutants+" mutants to the initial "+list.size()+" names.");
		for(String key : list) {collisions.remove(key);}
		int removed=0;
		if(!collisions.isEmpty()){
			for(String key : collisions) {
				names.remove(key);
				removed++;
			}
			System.err.println("Removed "+removed+" collisions due to ambiguity.");
		}
		return mutants;
	}
	
	/** Generate mutants for one name, recursively */
	private static int addMutants(final byte[] keyArray, final String value, final int hdist, final HashMap<String, String> names, final HashSet<String> collisions){
		assert(hdist>0);
		int added=0;
		for(int i=0; i<keyArray.length; i++){
			final byte old=keyArray[i];
			if(AminoAcid.isACGTN(old)){//Only mutate bases
				for(byte b : symbols){
					if(b!=old){
						keyArray[i]=b;
						{
							String key=new String(keyArray);
							String oldValue=names.get(key);
							if(oldValue==null){
								names.put(key, value);
								added++;
//								System.err.println("Added "+key+"->"+value+"; oldValue="+oldValue);
							}else if(!oldValue.equals(value)){
//								assert(value.equals(names.get(s))) : "Collision between "+value+" and "+names.get(s)+" for mutant "+s;
								collisions.add(key);
//								System.err.println("Collision for "+key+"->"+value+"; oldValue="+oldValue);
							}
						}
						if(hdist>1){
							added+=addMutants(keyArray, value, hdist-1, names, collisions);
						}
					}
				}
				keyArray[i]=old;
			}
		}
		return added;
	}
	
	/** Returns the value this read maps to; meaning, the variable part of the output filename */
	private String getValue(Read r){
		String key=getKey(r);
		if(key==null){return null;}
		return nameList.isEmpty() ? key : names.get(key);
	}
	
	/** Generates a key from the read header */
	private String getKey(Read r){
		final String id=r.id;
		final int idlen=id.length();
		
		if(nameList.size()>0){
			if(substringMode){
				for(String s : nameList){
					if(id.contains(s)){return s;}
				}
				return null;
			}else if(affixLengths.length>0){
				for(int affixLen : affixLengths){
					final String sub=idlen>=affixLen ? prefixMode ? id.substring(0, affixLen) : id.substring(idlen-affixLen) : id;
					if(names.containsKey(sub)) {return sub;}
				}
				return null;
			}
		}
		
		final String name;
		if(fixedAffixLength>0){
			name=(id.length()<=fixedAffixLength ? id : prefixMode ? id.substring(0, fixedAffixLength) : id.substring(idlen-fixedAffixLength));
		}else if(delimiter!=null){
			if(column>-1){

				String[] split=id.split(delimiter);
				assert(split.length>1) : "Delimiter '"+delimiter+"' was not found in name '"+id+"'";
				
				int col=Tools.min(column, split.length-1);
				name=split[col];
				if(col!=column && !warned){
					System.err.println("*** WARNING! ***\n"
							+ "Only "+(col+1)+" columns for record "+id+"\n"
							+ "Further warnings will be suppressed.\n");
					warned=true;
					assert(errorState=true); //Change error state to true if assertions are enabled.
				}
			}else if(prefixMode){
				int idx=(delimiterChar>0 ? id.indexOf(delimiterChar) : id.indexOf(delimiter));
				assert(idx>=0) : "Delimiter '"+delimiter+"' was not found in name '"+id+"'";
				name=id.substring(0, idx);
			}else{
				int idx=(delimiterChar>0 ? id.lastIndexOf(delimiterChar) : id.lastIndexOf(delimiter));
				assert(idx>=0) : "Delimiter '"+delimiter+"' was not found in name '"+id+"'";
				name=id.substring(idx+delimiter.length());
			}
		}else{
			assert(perheader);
			name=id;
		}
		
		return name;
	}
	
	/** Statistics about demultiplexing */
	void printReport(BufferedMultiCros mcros){
		if(stats==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(stats, overwrite, append, true);
		bsw.start();

		ByteBuilder bb=new ByteBuilder();
		bb.append("#ReadsIn\t").append(readsProcessed).nl();
		bb.append("#BasesIn\t").append(basesProcessed).nl();
		bb.append("#ReadsOut\t").append(readsOut).nl();
		bb.append("#BasesOut\t").append(basesOut).nl();
		bb.append("#Name\tReads\tBases\n");
		if(names!=null && !names.isEmpty()){
			bb.append("Unmatched\t").append(readsUnmatched).tab().append(basesUnmatched).nl();
		}
		bsw.print(bb);
		
		if(mcros!=null){
			bb=mcros.report();
			bsw.print(bb);
		}
		
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Primary input file */
	private String in1=null;
	/** Optional input file for read 2 of paired reads in twin files */
	private String in2=null;
	
	/** Input qual file 1 */
	private String qfin1=null;
	/** Input qual file 2 */
	private String qfin2=null;

	/** Output file pattern 1 */
	private String out1=null;
	/** Output file pattern 2 */
	private String out2=null;

	/** Output qual file 1 */
	private String qfout1=null;
	/** Output qual file 2 */
	private String qfout2=null;

	/** Unmatched read output file 1 */
	private String outu1=null;
	/** Unmatched read output file 2 */
	private String outu2=null;
	
	/** File extension override for input files */
	private String extin=null;
	/** File extension override for output files */
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	/** Primary input file */
	private final FileFormat ffin1;
	
	/** Read 2 input file */
	private final FileFormat ffin2;
	
	/*--------------------------------------------------------------*/
	
	/** Input reads */
	long readsProcessed=0;
	long basesProcessed=0;
	
	/** Demultiplexed output reads */
	long readsOut=0;
	long basesOut=0;
	
	/** Output reads that did not get demultiplexed */
	long readsUnmatched=0;
	long basesUnmatched=0;

	/** Stop after this many input reads */
	private long maxReads=-1;

	/** File to print number of reads sent to each output file */
	private String stats=null;

	/** For splitting headers on a symbol */
	private String delimiter=null;
	
	/** For splitting headers faster if the delimiter is just one character */
	private char delimiterChar=0;
	
	/** If there is a delimiter, use this column after splitting.
	 * Column is 1-based. */
	private int column=-1;
	
	/** Use the prefix of a header.  If false, use the suffix. 
	 * This is ignored if column is set. */
	private boolean prefixMode=true;
	
	/** Allow names to be any substring of a header */
	private boolean substringMode=false;
	
	/** Demultiplex every sequence into its own file */
	private boolean perheader=false;
	
	/** Prevents issuing warnings multiple times */
	private boolean warned=false;

	/** Mcros buffers will be dumped after accumulating this many reads */
	private int readsPerBuffer=1600;

	/** Mcros buffers will be dumped after accumulating this many reads */
	private int bytesPerBuffer=3000000;
	
	/** Hamming distance for read indexes */
	private int hdist=0;
	
	/** 
	 * Do not create files with under this many reads.
	 * Increases memory usage since reads must be retained until processing is finished.
	 * Upon running out of memory the program may go very slowly.
	 */
	private long minReadsToDump=0;
	
	/** Affix length if fixed */
	private int fixedAffixLength=-1;
	
	/** All possible affix lengths if there are multiple.
	 * Basically, the lengths of everything in names. */
	private int[] affixLengths;

	/** 
	 * A set of recognized names and the output names they map to.
	 * The keys and values will be identical unless a Hamming distance is used.
	 * This is only filled if explicit names are provided.
	 */
	private HashMap<String, String> names=new HashMap<String, String>();
	
	/** Names in list form, for substring matching */
	private ArrayList<String> nameList=new ArrayList<String>();
	
	/** 
	 * Select MultiCros2 or MultiCros3.
	 * MultiCros2 only opens 1 stream at a time and does not do concurrent writes.
	 * MultiCros3 maintains multiple streams and does concurrent writes even with only 1 stream.
	 */
	private int mcrosType=3;
	
	/** Max concurrent streams for MultiCros3.  Does not seem to speed things up. */
	private int maxStreams=BufferedMultiCros.DEFAULT_MAX_STREAMS;
	
	/** 
	 * This puts the mcros handling in a second thread.
	 * Does not seem to speed things up.
	 */
	private boolean threaded=false;
	
	/*--------------------------------------------------------------*/
	
	/** Print messages here */
	private PrintStream outstream=System.err;
	
	/** True if errors were encountered */
	public boolean errorState=false;
	
	/** Permission to overwrite existing files. */
	private boolean overwrite=true;
	
	/** 
	 * Append to existing files rather than overwriting.
	 * It iss not advisable to set this flag for this class. 
	 * */
	private boolean append=false;
	
	/** Retain header of sam/bam files. */
	private boolean useSharedHeader=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Symbols allowed to substitute for Hamming distance */
	private static final byte[] symbols={'A', 'C', 'G', 'T', 'N'};
	
	/** Verbose messages for debugging */
	public static boolean verbose=false;
	
}
