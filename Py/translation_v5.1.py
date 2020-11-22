#c:\Python27\python Translation_v5.1.py
import os,sys
import Tkinter as tk
from time import gmtime, strftime
#GGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATT
seq='GGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATTGGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATTGGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATTGGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATTGGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATTGGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATTGGGGGAGGTGGCGAGGAAGATGCTATATAAAAAGCTGTGCTCGTGATCGNNATCGNATTGGGGGAGGTGGCGAGGAAGATGCTATATAAAAAG'
out=['','','','','','','','',]

class Align_gui(object):
	def __init__(self):
		
		self.root = tk.Tk()

		self.root['bg'] = 'yellow'
		self.root.geometry('550x550+0+0')

		#Frames
		self.upper_frame = tk.Frame(self.root,height=10,width=400,bg='yellow')

		self.name_frame = tk.Frame(self.root,bd=1,highlightthickness=0,relief="sunken",padx=2,pady=0,width=0)

		self.intermediate0_frame = tk.Frame(self.root,height=15,width=400,bg='yellow')
		
		self.ORF_frame = tk.Frame(self.root,bd=1,highlightthickness=2,relief="sunken",padx=2,pady=0,width=400)
		self.intermediate1_frame= tk.Frame(self.root,height=15,width=400,bg='yellow')
	
		self.text_frame= tk.Frame(self.root,bd=1,highlightthickness=1,relief="sunken",padx=1,pady=0,width=400)
		self.intermediate2_frame= tk.Frame(self.root,height=10,width=400,bg='yellow')
	
		self.save_frame= tk.Frame(self.root,bd=1,highlightthickness=0,relief="sunken",padx=2,pady=0,height=2,width=100,bg='yellow')
		self.reset_frame= tk.Frame(self.root,bd=1,highlightthickness=0,relief="sunken",padx=2,pady=0,height=2,width=100,bg='yellow')
		self.radiovar = tk.IntVar()

		#Elements (labels, spinboxes, buttons)
		self.name_label = tk.Label(self.name_frame,text='DNA translation tool',bg ='white',fg='blue')
		self.name_label.configure(font=('Courier 10 bold'))
		self.ORF1 = tk.Spinbox(self.ORF_frame,values=('ORF1','ORF2','ORF3','ORF-1','ORF-2','ORF-3',)) 
	
		self.out1 = tk.Spinbox(self.ORF_frame,values=('DNA and Protein','Protein')) 
	
		self.ORF2 = tk.Button(self.ORF_frame,text='Translate',command=self.Check_seq)
		
		self.save_button = tk.Button(self.save_frame,text = 'Save in FASTA format',bg ='white',command=self.Save_FASTA)

		self.reset_button = tk.Button(self.reset_frame,text = 'RESET',bg ='white',command=self.Refresh)

		#Textbox
		self.textbox = tk.Text(self.text_frame)
		self.textbox.config(font='Arial 2')

		#for l in seq:										
		#self.textbox.insert('end',l)
		self.textbox.insert('end',seq)
		#self.textbox.insert('end','\n')

		#Packing frames
		self.upper_frame.pack()
		self.name_frame.pack()
		self.intermediate0_frame.pack()
		self.ORF_frame.pack()
		self.intermediate1_frame.pack()
		self.text_frame.pack()
		self.intermediate2_frame.pack()
		self.save_frame.pack(side=tk.LEFT)
		self.reset_frame.pack(side=tk.RIGHT)

		#Packing elements
		self.name_label.pack()
		self.ORF1.pack(side=tk.LEFT)
		self.out1.pack(side=tk.LEFT)
		self.ORF2.pack(side=tk.LEFT)
		self.save_button.pack()
		self.reset_button.pack()

		#Packing textbox
		self.textbox.pack()
		self.textbox.config(relief=tk.RAISED,font='Courier 10',fg='blue',width=60,height=25)
		self.root.mainloop()

	def Check_seq(self):
		out=['','','','','','','','',]
		seqD1=self.textbox.get('1.0','end')
		if seqD1[:1]=='\n':
			self.Results()
		else:
			seqD1=seqD1[:-1].upper()
			seqD2=''
			if '-' in seqD1:
				for i in xrange(len(seqD1)):
					if seqD1[i]!='-':
						seqD2+=seqD1[i]
				seqD1=seqD2
			counter=0
			for i in xrange(len(seqD1)):
				if seqD1[i]!='A' and seqD1[i]!='C' and seqD1[i]!='G' and seqD1[i]!='T' and seqD1[i]!='N':
					counter+=1
			if counter>0:
				self.Refresh()
				self.textbox.insert('1.0', 'Please, paste a DNA sequence')
			else:
				self.Transl_org(seqD1)	
		return

	def Inv_compl(self,seqD1):
		comp={'A':'T','T':'A','G':'C','C':'G','N':'N'}

		rev_compl = ''

		for i in xrange(-1,-len(seqD1)-1,-1)	:

			rev_compl = rev_compl + comp[seqD1[i]]

		return rev_compl


	def Transl(self,seqD1):
		gencode={'GGG':'G','GGA':'G','GGT':'G','GGC':'G','GAG':'E','GAA':'E','GAT':'D','GAC':'D','GTG':'V','GTA':'V','GTT':'V','GTC':'V','GCG':'A','GCA':'A','GCT':'A','GCC':'A','AGG':'R','AGA':'R','CGG':'R','CGA':'R','CGT':'R','CGC':'R','AAG':'K','AAA':'K','AAT':'N','AAC':'N','ATG':'M','ATA':'I','ATT':'I','ATC':'I','ACG':'T','ACA':'T','ACT':'T','ACC':'T','TGG':'W','TGA':'*','TAG':'*','TAA':'*','TGT':'C','TGC':'C','TAT':'Y','TAC':'Y','TTG':'L','TTA':'L','CTG':'L','CTA':'L','CTT':'L','CTC':'L','TTT':'F','TTC':'F','AGT':'S','AGC':'S','TCG':'S','TCA':'S','TCT':'S','TCC':'S','CAG':'Q','CAA':'Q','CAT':'H','CAC':'H','CCG':'P','CCA':'P','CCT':'P','CCC':'P','GGN':'G','GAN':'X','GTN':'V','GCN':'A','CAN':'X','CCN':'P','CGN':'R','CTN':'L','AGN':'X','AAN':'X','ATN':'X','ACN':'T','TGN':'X','TAN':'X','TTN':'X','TCN':'S','GNG':'X','GNA':'X','GNT':'X','GNC':'X','ANG':'X','ANA':'X','ANT':'X','ANC':'X','CNG':'X','CNA':'R','CNT':'X','CNC':'X','TNG':'X','TNA':'X','TNG':'X','TNC':'X','NAA':'X','NAT':'X','NAG':'X','NAC':'X','NGG':'X','NGA':'X','NGT':'X','NGC':'X','NCC':'X','NCA':'X','NCT':'X','NCG':'X','NTT':'X','NTA':'X','NTC':'X','NTG':'X','NNA':'X','NNG':'X','NNC':'X','NNT':'X','NAN':'X','NGN':'X','NCN':'X','NTN':'X','ANN':'X','GNN':'X','CNN':'X','TNN':'X','NNN':'X','A':'','G':'','C':'','T':'','AA':'','AC':'T','AG':'','AT':'','CA':'','CC':'P','CG':'R','CT':'L','GA':'','GC':'A','GG':'G','GT':'','TA':'','TC':'S','TG':'','TT':'','N':'','AN':'','CN':'','GN':'','TN':'','NN':'','NA':'','NG':'','NC':'','NT':'','TNT':'X'}
		seqP1=''
		for i in xrange(0,len(seqD1),3):

			seqP1+=gencode[seqD1[i:i+3]]
		return seqP1

	def Transl_org(self,seqD1):
		out[0]=seqD1
		out[1]=self.Inv_compl(seqD1)
		out[2]=self.Transl(out[0])
		out[3]=self.Transl(out[0][1:])
		out[4]=self.Transl(out[0][2:])
		out[5]=self.Transl(self.Inv_compl(out[0]))
		out[6]=self.Transl(self.Inv_compl(out[0])[1:])
		out[7]=self.Transl(self.Inv_compl(out[0])[2:])
		self.Results()
		return 

	def Refresh(self):
		self.textbox.delete('1.0', 'end')
		seqD1=''
		out=['','','','','','','','',]
		return 

	def Results(self):
		outD=[]
		outP=[]
		ORF_dic={1:out[2],2:out[3],3:out[4],-1:out[5],-2:out[6],-3:out[7]}
		ORF=int(self.ORF1.get()[3:])
		dir_rev1={1:out[0][abs(ORF)-1:],2:out[1][abs(ORF)-1:]}
		if ORF>0:
			dir_rev=1
		else:
			dir_rev=2
		if self.out1.get()=='DNA and Protein':
			total=60
			quantD=total
			digits=len(str(len(dir_rev1[dir_rev])))
			length=quantD+digits*2+6
			while length > total:
				quantD-=1
				length=quantD+len(str(len(dir_rev1[dir_rev])))*2+6
			quantD=int(quantD/3*3)
			lines=int(len(dir_rev1[dir_rev])/float(quantD))
			if lines<len(dir_rev1[dir_rev])/float(quantD):
				lines+=1
			out1=' '
			for i in xrange(len(ORF_dic[ORF])):
				out1+=ORF_dic[ORF][i]+2*' '
			ORF_dic[ORF]=out1
			for i in xrange(0,len(dir_rev1[dir_rev]),quantD):
				outD.append(' '+str(i+1)+(digits-len(str(i+1)))*' '+2*' '+dir_rev1[dir_rev][i:i+quantD]+2*' '+((quantD-len(dir_rev1[dir_rev][i:i+quantD]))+digits-len(str(i+quantD+1)))*' '+str(i+len(dir_rev1[dir_rev][i:i+quantD])))
				outP.append(' '+str(i/3+1)+(digits-len(str(i/3+1)))*' '+2*' '+ORF_dic[ORF][i:i+quantD]+2*' '+((quantD-len(ORF_dic[ORF][i:i+quantD]))+digits-len(str((i+quantD+1)/3)))*' '+str(i/3+len(dir_rev1[dir_rev][i:i+quantD])/3))
			self.textbox.delete('1.0', 'end')
			self.textbox.insert('end','\n')
			for i in xrange(len(outD)):
				self.textbox.insert('end',outD[i])
				self.textbox.insert('end','\n')
				self.textbox.insert('end',outP[i])
				self.textbox.insert('end','\n')
		else:
			outP=[]
			Plength=len(ORF_dic[ORF])
			total=60
			quantP=total
			digits=len(str(Plength))
			length=quantP+digits*2+6
			while length>total:
				quantP-=1
				length=quantP+len(str(Plength))*2+6
			quantP=int(quantP/3*3)
			lines=int(len(dir_rev1[dir_rev])/float(quantP))
			if lines<len(dir_rev1[dir_rev])/float(quantP):
				lines+=1
			for i in xrange(0,Plength,quantP):
				outP.append(' '+str(i+1)+(digits-len(str(i+1)))*' '+2*' '+ORF_dic[ORF][i:i+quantP]+2*' '+((quantP-len(ORF_dic[ORF][i:i+quantP]))+digits-len(str(i+quantP+1)))*' '+str(i+len(dir_rev1[dir_rev][i:i+quantP])))
			self.textbox.delete('1.0', 'end')
			self.textbox.insert('end','\n')
			for i in xrange(len(outP)):
				self.textbox.insert('end',outP[i])
				self.textbox.insert('end','\n')
		return

	def Save_FASTA(self):
		content=[]
		ORF_dic={1:out[2],2:out[3],3:out[4],-1:out[5],-2:out[6],-3:out[7]}
		ORF=int(self.ORF1.get()[3:])
		content.append('>Translation ORF ('+str(ORF)+')')

		for i in xrange(0,len(ORF_dic[ORF]),60):
			content.append(ORF_dic[ORF][i:i+60])
		self.textbox.delete('1.0', 'end')
		file=open('D:/Py/Translation/Results/Translation_ORF_'+str(ORF)+'_'+strftime("%d.%m.%y_%H.%M", gmtime())+'.FASTA','a')
		self.textbox.insert('end','\n')
		for l in content:										
			file.write("%s\n" % l)
			self.textbox.insert('end',l)
			self.textbox.insert('end','\n')
		file.close()
		return

Align_gui()































