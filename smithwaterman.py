

import re

MATCH = 0;
GAPINROW = 1;
GAPINCOL = 2;

blosumlines= ["A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *",
"A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4",
"R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4",
"N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4",
"D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4",
"C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4",
"Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4",
"E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4",
"G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4",
"H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4",
"I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4",
"L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4",
"K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4",
"M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4",
"F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4",
"P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4",
"S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4",
"T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4",
"W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4",
"Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4",
"V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4",
"B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4",
"Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4",
"X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4",
"* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1"
];
class SmithWaterman:
	def __init__(self):
		
		self.smat = {};
		self.go = 10.0;#must be positive
		self.ge = 0.5;
		heads = re.split("[\s]+",blosumlines[0]);
		for ii in range(1,len(blosumlines)):
			pt = re.split("[\s]+",blosumlines[ii]);
			for jj in range(1,len(pt)):
				self.smat[pt[0]+"_"+heads[jj-1]] = int(pt[jj]);
				
				
	def makeDPMatrix(self,seqA,seqB):
		self.aaA = list(re.sub("[^A-Z]","",seqA.upper()));
		self.aaB = list(re.sub("[^A-Z]","",seqB.upper()));	
		self.dpmat = [[0 for i in range(len(self.aaA)+1)] for j in range(len(self.aaB)+1)]
		
		for rr in range(len(self.aaB)+1):
			for cc in range(len(self.aaA)+1):
				sc = SWCell(rr,cc);
				sc.setMaxScores(self);
				self.dpmat[rr][cc] = sc;
	def backTrack(self):
		
		maxpos_r = 0;
		maxpos_c = 0;
		maxpath = 0;
		maxscore = 0;
		for rr in range(len(self.aaB)+1):
			for cc in range(len(self.aaA)+1):
				if(self.dpmat[rr][cc].score[MATCH] > maxscore):
					maxscore = self.dpmat[rr][cc].score[MATCH];
					maxpos_r = rr;
					maxpos_c = cc;
					maxpath = self.dpmat[rr][cc].path[MATCH];
		
		kr = maxpos_r;
		kc = maxpos_c;
		ks = maxscore;
		kp = maxpath;
		cstate = MATCH;
		
		saa = [];
		sbb = [];
		pcell = self.dpmat[kr][kc];
		
		while ks > 0:
			ccell = 0;
			pcell = self.dpmat[kr][kc];
			nstate = MATCH;
			if(cstate == MATCH):
				saa.append(self.aaA[kc-1]);
				sbb.append(self.aaB[kr-1]);
				ccell = self.dpmat[kr-1][kc-1];
				nstate = pcell.path[MATCH];
				kc -= 1;
				kr -= 1;
				
				
			if(cstate == GAPINCOL):
				saa.append("-");
				sbb.append(self.aaB[kr-1]);
				ccell = self.dpmat[kr-1][kc];
				nstate = pcell.path[GAPINCOL];
				kr -= 1;
				
				
			if(cstate == GAPINROW):
				saa.append(self.aaA[kc-1]);
				sbb.append("-");
				ccell = self.dpmat[kr][kc-1];
				nstate = pcell.path[GAPINROW];
				kc -= 1;
			
			cstate = nstate;
			ks = ccell.score[cstate];
		
		lastc = kc;
		lastr = kr;
		if(False):
			if(cstate == MATCH):
				lastr = kr+1;
				lastc = kc+1;
			
			if(cstate == GAPINCOL):
				lastr = kr+1;
			
			if(cstate == GAPINROW):
				lastc = kc+1;
		
		
		if(False):
		# add nterminals
			for ii in range(1,lastc+1):
				saa.append(self.aaA[lastc-ii]);
				sbb.append("-");
				
			for ii in range(1,lastr+1):
				sbb.append(self.aaB[lastr-ii]);
				saa.append("-");
			saa.reverse();
			sbb.reverse();
			# add cterminals
			for ii in range(maxpos_c+1,len(self.aaA)):
				saa.append(self.aaA[ii]);
				sbb.append("-");
				
			for ii in range(maxpos_r+1,len(self.aaB)):
				sbb.append(self.aaB[ii]);
				saa.append("-");
		else:
			saa.reverse();
			sbb.reverse();
		
		
		print("".join(saa));
		print("".join(sbb));
		
	

class SWCell:
	def __init__(self,r,c):
		self.r = r;
		self.c = c;
		self.score = [0]*3;
		self.path = [0]*3;
		
	def setMaxScores(self,sw):
		if(self.r == 0 or self.c == 0):
			self.score[0] = 0;
			self.score[1] = 0;
			self.score[2] = 0;
			return;
		r = self.r;
		c = self.c;
		
		#MATCH
		matchcode = sw.aaA[self.c-1]+"_"+sw.aaB[self.r-1];
		matchscore = -4;
		if(matchcode in sw.smat):
			matchscore = sw.smat[matchcode];
		#result was different from that of EMBOSS water when I was using = 
		
		if(sw.dpmat[r-1][c-1].score[MATCH] + matchscore > sw.dpmat[r-1][c-1].score[GAPINROW] + matchscore):
			if(sw.dpmat[r-1][c-1].score[MATCH] + matchscore > sw.dpmat[r-1][c-1].score[GAPINCOL] + matchscore):
				self.path[MATCH] = MATCH;
				self.score[MATCH] = sw.dpmat[r-1][c-1].score[MATCH] + matchscore;
			else:
				self.path[MATCH] = GAPINCOL;
				self.score[MATCH] = sw.dpmat[r-1][c-1].score[GAPINCOL] + matchscore;
		else:
			if(sw.dpmat[r-1][c-1].score[GAPINROW] + matchscore > sw.dpmat[r-1][c-1].score[GAPINCOL] + matchscore):
				self.path[MATCH] = GAPINROW;
				self.score[MATCH] = sw.dpmat[r-1][c-1].score[GAPINROW] + matchscore;
			else:
				self.path[MATCH] = GAPINCOL;
				self.score[MATCH] = sw.dpmat[r-1][c-1].score[GAPINCOL] + matchscore;
			
		#GAPINCOL
		if(sw.dpmat[r-1][c].score[MATCH] - sw.go > sw.dpmat[r-1][c].score[GAPINCOL] - sw.ge):
			self.path[GAPINCOL] = MATCH;
			self.score[GAPINCOL] = sw.dpmat[r-1][c].score[MATCH] - sw.go;
		else:
			self.path[GAPINCOL] = GAPINCOL;
			self.score[GAPINCOL] = sw.dpmat[r-1][c].score[GAPINCOL] - sw.ge;
		
		#GAPINROW
		if(sw.dpmat[r][c-1].score[MATCH] - sw.go >= sw.dpmat[r][c-1].score[GAPINROW] - sw.ge):
			self.path[GAPINROW] = MATCH;
			self.score[GAPINROW] = sw.dpmat[r][c-1].score[MATCH] - sw.go;
		else:
			self.path[GAPINROW] = GAPINROW;
			self.score[GAPINROW] = sw.dpmat[r][c-1].score[GAPINROW] - sw.ge;
		
		for ii in range(len(self.score)):
			if(self.score[ii] < 0):
				self.score[ii] = 0;
		
		
		


#for ii in range(0,len(blosumlines)):
#	blosumlines[ii]  = re.sub("[\s]+$","",re.sub("^[\s]+","",blosumlines[ii]));



sw = SmithWaterman();

str2 = "HTFSYRSHGPVIDLAEKLVSMAPVPMSKAYFTNSGSEANDTVVKL";
str1 = "LEIAKVSQDDIVYALGCGDGRIIITAAKDFNVK";
sw.makeDPMatrix(str1,str2);

sw.backTrack();