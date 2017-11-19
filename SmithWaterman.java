import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * public domain snippet for SmithWaterman alignment
 * Smith, Temple F., and Michael S. Waterman. 
 * "Identification of common molecular subsequences."
 * Journal of molecular biology 147.1 (1981): 195-197.
 * 
 * author: yamule (https://github.com/yamule/smithwaterman)
 * usage ===
 *	SmithWaterman sw = new SmithWaterman();
 *	SWResult res = sw.align("EEEEMDQNNSLPPYAGGTWRYII","IIIIMDQNNSPPYAQGGTWRYEE");
 *	System.out.println("score: "+res.score);
 *	System.out.println(SmithWaterman.listToString(res.qseq));
 *	System.out.println(SmithWaterman.listToString(res.sseq));
 *		
 * output ===
 * score: 73
 * ----EEEEMDQNNSLPPYA-GGTWRYII--
 * IIII----MDQNNS-PPYAQGGTWRY--EE
 * 
 */
public class SmithWaterman {
	public static int TYPE_MATCH = 0;
	public static int TYPE_GAPINCOL = 1;
	public static int TYPE_GAPINROW = 2;
	SWCell[][] dpMat;
	ScoringMatrix smat = new BLOSUM62();
	double penalO = 10;
	double penalE = 0.5;
	
	
	public SWResult align(String qq,String ss){
		String qq2 = qq.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","");
		String ss2 = ss.replaceAll("^[\\s]*>[^\\r\\n]*[\\r\\n]+","");
		
		String name1 = "seq1";
		String name2= "seq2";
		Pattern npat = Pattern.compile(">([^\\s]+)");
		Matcher mat = npat.matcher(qq);
		if(mat.find()){
			name1 = mat.group(1);
		}
		 mat = npat.matcher(ss);
		if(mat.find()){
			name2 = mat.group(1);
		}
		
		Sequence q = new Sequence();
		Sequence s = new Sequence();
		q.name = name1;
		s.name = name2;
		q.add(qq2);
		s.add(ss2);
		return align(q,s);
		
		
	}
	
	public SWResult align(Sequence qq,Sequence ss){
		//Prepare letters for alignment.
		//Empty objects are removed and unknown letters are changed to 'X'.
		//Please see the filter method for exact process.
		ArrayList<Character> q = smat.filter(qq.seq);
		ArrayList<Character> s = smat.filter(ss.seq);
		int w = q.size()+1;
		int h = s.size()+1;
		
		dpMat = new SWCell[w][h];
		
		for(int xx = 0;xx < w;xx++){
			for(int yy = 0;yy < h;yy++){
				dpMat[xx][yy] = new SWCell();
			}
		}
		for(int xx = 0;xx < w;xx++){
			
			for(int ii = 0;ii < 3;ii++){
				dpMat[xx][0].setScoreAt(ii,0);
			}
		}
		for(int yy = 0;yy < h;yy++){
			for(int ii = 0;ii < 3;ii++){
				dpMat[0][yy].setScoreAt(ii,0);
			}
		}
		for(int xx = 1;xx <  w;xx++){
			for(int yy = 1;yy <  h;yy++){
				double fm = dpMat[xx-1][yy-1].getScoreAt(TYPE_MATCH) + smat.getScore(q.get(xx-1), s.get(yy-1));
				double fc = dpMat[xx-1][yy-1].getScoreAt(TYPE_GAPINCOL) + smat.getScore(q.get(xx-1), s.get(yy-1));
				double fr = dpMat[xx-1][yy-1].getScoreAt(TYPE_GAPINROW) + smat.getScore(q.get(xx-1), s.get(yy-1));
				SWCell currentcell = dpMat[xx][yy];
				fm = Math.max(0,fm);
				fc = Math.max(0,fc);
				fr = Math.max(0,fr);
				
				if(fm >= fr){
					if(fm >= fc){
						currentcell.setScoreAt(TYPE_MATCH,fm);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_MATCH);
					}else{
						currentcell.setScoreAt(TYPE_MATCH,fc);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINCOL);
					}
				}else{
					if(fr > fc){
						currentcell.setScoreAt(TYPE_MATCH,fr);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINROW);
			
					}else{
						currentcell.setScoreAt(TYPE_MATCH,fc);
						currentcell.setPrevTypeAt(TYPE_MATCH, TYPE_GAPINCOL);
					}
				}
				
				
				//Gap in col
				fm = dpMat[xx][yy-1].getScoreAt(TYPE_MATCH)-this.penalO;
				fc = dpMat[xx][yy-1].getScoreAt(TYPE_GAPINCOL)-this.penalE;
				fr = dpMat[xx][yy-1].getScoreAt(TYPE_GAPINROW)-this.penalO;// not used
				
				fm = Math.max(0,fm);
				fc = Math.max(0,fc);
				fr = Math.max(0,fr);
				
				
				if(fm >= fc){
					currentcell.setScoreAt(TYPE_GAPINCOL, fm);
					currentcell.setPrevTypeAt(TYPE_GAPINCOL, TYPE_MATCH);
				}else{
					currentcell.setScoreAt(TYPE_GAPINCOL, fc);
					currentcell.setPrevTypeAt(TYPE_GAPINCOL, TYPE_GAPINCOL);
				}
				
				//Gap in row
				fm = dpMat[xx-1][yy].getScoreAt(TYPE_MATCH)-this.penalO;
				fc = dpMat[xx-1][yy].getScoreAt(TYPE_GAPINCOL)-this.penalO;// not used
				fr = dpMat[xx-1][yy].getScoreAt(TYPE_GAPINROW)-this.penalE;
				
				fm = Math.max(0,fm);
				fc = Math.max(0,fc);
				fr = Math.max(0,fr);
				
				if(fm > fr){
					currentcell.setScoreAt(TYPE_GAPINROW, fm);
					currentcell.setPrevTypeAt(TYPE_GAPINROW, TYPE_MATCH);
				}else{
					currentcell.setScoreAt(TYPE_GAPINROW, fr);
					currentcell.setPrevTypeAt(TYPE_GAPINROW, TYPE_GAPINROW);
				}
			}
		}
		
		int maxpos[] = getMaxScorePos();
		int cx = maxpos[0];
		int cy = maxpos[1];
		int sx = cx;
		int sy = cy;
		int lx = cx;
		int ly = cy;
		if(cx == -1){//very short sequences which do not have positive value match between two.
			return new SWResult();
		}
		
		
		
		//backtrackpart
		double cs = dpMat[cx][cy].getScoreAt(TYPE_MATCH);
		double maxscore = cs;
		int ct = TYPE_MATCH;
		int cpt = dpMat[cx][cy].getPrevTypeAt(TYPE_MATCH);
		ArrayList<Character> qchar = new ArrayList<>();
		ArrayList<Character> schar = new ArrayList<>();
		qchar.add(q.get(cx-1));
		schar.add(s.get(cy-1));
		while(cs > 0){
			if(ct == TYPE_MATCH){
				cx--;
				cy--;
			}else if(ct == TYPE_GAPINCOL){
				cy--;
			}else if(ct == TYPE_GAPINROW){
				cx--;
			}
			ct = cpt;
			cs = dpMat[cx][cy].getScoreAt(ct);
			cpt = dpMat[cx][cy].getPrevTypeAt(ct);
			if(cs <= 0){
				break;
			}
			if(ct == TYPE_MATCH){
				qchar.add(q.get(cx-1));
				schar.add(s.get(cy-1));
				lx = cx;
				ly = cy;
			}else if(ct == TYPE_GAPINCOL){
				qchar.add('-');
				schar.add(s.get(cy-1));
				ly = cy;
			}else if(ct == TYPE_GAPINROW){
				qchar.add(q.get(cx-1));
				schar.add('-');
				lx = cx;
			}
		}
		
		
		//add unaligned nterminal part
		for(int xx = lx-1;xx > 0;xx--){
			qchar.add(q.get(xx-1));
			schar.add('-');
		}
		for(int yy = ly-1;yy > 0;yy--){
			schar.add(s.get(yy-1));
			qchar.add('-');
		}
		
		Collections.reverse(qchar);
		Collections.reverse(schar);
		
		
		//add unaligned cterminal part
		for(int xx = sx+1;xx < dpMat.length;xx++){
			qchar.add(q.get(xx-1));
			schar.add('-');
		}
		for(int yy = sy+1;yy < dpMat[0].length;yy++){
			schar.add(s.get(yy-1));
			qchar.add('-');
		}
		
		
		return new SWResult(qchar,schar,maxscore);
		
	}
	
	
	
	
	/**
	 * Returns x y position of the cell which has maximum scores.
	 * If there are multiple cells with the same score, the cell which has lowest x and lowest y is selected.
	 * @return 
	 */
	public int[] getMaxScorePos(){
		int w = dpMat.length;
		int h = dpMat[0].length;
		int ret[] = new int[2];
		ret[0] = -1;
		ret[1] = -1;
		double maxscore = 0;
		for(int xx = 0;xx < w;xx++){
			for(int yy = 0;yy < h;yy++){
				for(int i = 0;i < 3;i++){
					double sc = dpMat[xx][yy].getScoreAt(i);
					if(maxscore < sc){
						maxscore = sc;
						ret[0] = xx;
						ret[1] = yy;
					}
				}
			}	
		}
		return ret;
	}
	
	
	
	/**
	 * Changes ArrayList<Character> to String.
	 * @param al
	 * @return 
	 */
	public static String listToString(ArrayList<Character> al){
		StringBuffer ret = new StringBuffer();
		for(Character c:al){
			ret.append(c);
		}
		return ret.toString();
	}
	public static void main(String[] args){
		
		for(int ii = 0;ii < args.length;ii++){
			System.out.println(ii+" "+args[ii]);
		}
		ArrayList<Sequence> seq1 = Sequence.loadFasta(args[0]);
		ArrayList<Sequence> seq2 = Sequence.loadFasta(args[1]);
		
		SmithWaterman sw = new SmithWaterman();
		SWResult res = sw.align(seq1.get(0),seq2.get(0));
		System.out.println("# score: "+res.score);
		System.out.println(">"+seq1.get(0).name);
		System.out.println(listToString(res.qseq));
		System.out.println(">"+seq2.get(0).name);
		System.out.println(listToString(res.sseq));
		
		/*
		// for debug
		SmithWaterman sw = new SmithWaterman();
		SWResult res = sw.align("MTPPPPGRAAPSAPRARVPGPPARLGLPLRLRLLLLLWAAAASAQGHLRSGPRIFAVWKG","PGPGNPSPMSLSPAWPGHPDQPLPREQMTSPAPPRIITSATADPEGTETALAGDTSDGLA");
		System.out.println("score: "+res.score);
		System.out.println(SmithWaterman.listToString(res.qseq));
		System.out.println(SmithWaterman.listToString(res.sseq));
		*/
		/*
		BLOSUM62 bl = new BLOSUM62();
		SmithWaterman sw = new SmithWaterman();
		ArrayList<Sequence> seqs = Sequence.loadFasta("C:\\Users\\kimidori\\Desktop\\TBM\\sw\\testfas.txt");
		for(Sequence s:seqs){
			System.out.println(">"+s.name+" "+s.desc);
			for(Character c:s.seq){
				System.out.print(String.valueOf(c));
			}
			System.out.println();
		}
		SWResult res = sw.align(seqs.get(0),seqs.get(1));
		System.out.println(res.score);
		System.out.println(seqs.get(0).name);
		System.out.println(listToString(res.qseq));
		System.out.println(seqs.get(1).name);
		System.out.println(listToString(res.sseq));
		*/
	}
}



/**
 * Stores result of SmithWaterman Alignment.
 * @author kimidori
 */
class SWResult{
	ArrayList<Character> qseq = new ArrayList<>();
	ArrayList<Character> sseq = new ArrayList<>();
	double score = -1;
	SWResult(){
	}
	SWResult(ArrayList<Character> qq,ArrayList<Character> ss,double s){
		qseq.addAll(qq);
		sseq.addAll(ss);
		score = s;
	}
}



class SWCell{
	double score[] = new double[3];
	int prevType[] = new int[3];
	public void setScoreAt(int type,double sc){
		score[type] = sc;
	}
	public double getScoreAt(int type){
		return score[type];
	}
	public int getPrevTypeAt(int type){
		return prevType[type];
	}
	public void setPrevTypeAt(int ctype,int ptype){
		prevType[ctype] = ptype;
	}
	
	
}

class BLOSUM62 implements ScoringMatrix{
	
	
	String lines[] = {
		//https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *",
"A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 ",
"R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 ",
"N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 ",
"D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 ",
"C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 ",
"Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 ",
"E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
"G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 ",
"H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 ",
"I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 ",
"L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 ",
"K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 ",
"M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 ",
"F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 ",
"P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 ",
"S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 ",
"T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 ",
"W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 ",
"Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 ",
"V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 ",
"B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 ",
"Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ",
"X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 ",
"* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 "
	};
	
	
	HashSet<Character> acceptable = new HashSet<>();
	int[][] scoreMat = new int[128][128];
	
	BLOSUM62(){
		for(int ii = 0;ii < 128;ii++){
			for(int jj = 0;jj < 128;jj++){
				scoreMat[ii][jj] = 10000;
			}	
		}
		ArrayList<String> head = splitWithSpace(lines[0]);
		for(int ii = 1;ii < lines.length;ii++){
			ArrayList<String> pt = splitWithSpace(lines[ii]);
			String s = pt.get(0);
			char sc = s.charAt(0);
			acceptable.add(sc);
			for(int jj = 1;jj < pt.size();jj++){
				String q = head.get(jj-1);
				char qc = q.charAt(0);
				acceptable.add(qc);
				if(scoreMat[sc][qc] < 100){
					System.out.println(sc+"-"+qc+" duplicate pair?");
				}
				scoreMat[sc][qc] = Integer.parseInt(pt.get(jj));
			}
		}
		/*
		ArrayList<Character> aa = new ArrayList<>(acceptable);
		for(Character k:aa){
			for(Character q:aa){
				int sc3 = scoreMat[q][k];
				int sc4 = scoreMat[k][q];
				if(sc3 == sc4){
					System.out.println(k+"-"+q+" OK "+sc3);
				}else{
					System.out.println(k+q+" is different, "+sc3+","+sc4+",");
				}
			}
		}
		*/
		
	}
	
	
	/**
	 * Change letters which are not compatible with this matrix.
	 * @param al 
	 */
	public ArrayList<Character> filter(ArrayList<Character> al){
		ArrayList<Character> ret = new ArrayList<>();
		for(int ii = 0;ii < al.size();ii++){
			Character a = al.get(ii);
			if(a == null || a == Character.MIN_VALUE){
			}else{
				if(a.equals('-') || a.equals('.')){
				}else{
					if(acceptable.contains(a)){
						ret.add(a);
					}else{
						ret.add('X');
					}
				}
			}
		}
		return ret;
	}
	
	
	
	//For discrepancy between open java & oracle java
	public ArrayList<String> splitWithSpace(String str){
		ArrayList<String> ret = new ArrayList<>();
		String head[] = str.replaceAll("^[\\s]+","").replaceAll("[\\s]+$","").split("[\\s]+");
		for(String s:head){
			if(s.length() > 0){
				ret.add(s);
			}
		}
		return ret;
	}
	
	public int getScore(char q,char s){
		return scoreMat[q][s];
	}
}



interface ScoringMatrix{
	public int getScore(char q,char s);
	public ArrayList<Character> filter(ArrayList<Character> al);
}

class Sequence{
	String name;
	String desc;
	ArrayList<Character> seq = new ArrayList<Character>();
	public static ArrayList<Sequence> loadFasta(String filename){
		ArrayList<Sequence> ret = new ArrayList<>();
		try{
			BufferedReader br = new  BufferedReader(new FileReader(new File(filename)));
			String ln = null;
			Sequence currentseq = new Sequence();
			ret.add(currentseq);
			Pattern pat1 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]*)");
			Pattern pat2 = Pattern.compile("^[\\s]*>[\\s]*([^\\s]+)[\\s]+([^\\r\\n]+)");
			while((ln = br.readLine()) != null){
				Matcher mat = pat1.matcher(ln);
				if(mat.find()){
					String n = mat.group(1);
					String d = "";
					Matcher mat2 = pat2.matcher(ln);
					if(mat2.find()){
						d = mat2.group(2);
					}
					if(ret.size() == 1 && currentseq.seq.size() == 0){
					}else{
						currentseq = new Sequence();
						ret.add(currentseq);
					}
					currentseq.name = n;
					currentseq.desc = d;
				}else{
					currentseq.add(ln);
				}
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
		return ret;
	}
	public void add(String s){
		String[] pt = s.replaceAll("[\\r\\n]","").split("");
		for(String pp:pt){
			if(pp.length() == 1){
				seq.add(pp.charAt(0));
			}else{
				if(pp.length() == 0){
					
				}else{
					throw new RuntimeException("java implementation error.");
				}
			}
		}
	}
}
