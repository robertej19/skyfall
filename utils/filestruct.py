import pandas as pd
from icecream import ic

class fs:
	def __init__(self):
		self.phibins_test = [0,180,360]  
		self.tbins_test =  [0,1,12]
		self.xBbins_test = [0,.35,1]
		self.q2bins_test = [1,4,14]
		self.phibins =  [0,18,36,54,72,90,108,126,144,162,180,198,216,234,252,270,288,306,324,342,360]
		self.tbins = [0,0.2,0.4,0.6,0.8,1,2,4,7,12]
		self.xBbins = [0,0.2,.4,.6,.8, 1]
		self.q2bins = [0,1,2,2.5,3,3.5,4,4.5,5,6,10,14]
		self.f18_inbending_total_lumi = 5.511802214933477e+40
		self.f18_inbending_total_lumi_inv_nb = 5.511802214933477e+7
		self.Ebeam = 10.6 #GeV, CLAS12 beam energy
		self.m_p = 0.938 #GeV, proton mass
		self.alpha = 1/137 #, fine structure constant
		



		