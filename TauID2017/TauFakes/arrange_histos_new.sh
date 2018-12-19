mkdir -p Arranged_DY
hadd -f Arranged_DY/DY_0.root Out_DY*amc*0/Con*DY*.root 
hadd -f Arranged_DY/DY_1.root Out_DY*amc*1/Con*DY*.root 
hadd -f Arranged_DY/DY_2.root Out_DY*amc*2/Con*DY*.root 
hadd -f Arranged_DY/DY_3.root Out_DY*amc*3/Con*DY*.root 
hadd -f Arranged_DY/DY_4.root Out_DY*amc*4/Con*DY*.root 
hadd -f Arranged_DY/DY_5.root Out_DY*amc*5/Con*DY*.root 
hadd -f Arranged_DY/DY_6.root Out_DY*amc*6/Con*DY*.root 
hadd -f Arranged_DY/DY_7.root Out_DY*amc*7/Con*DY*.root 
hadd -f Arranged_DY/DY_8.root Out_DY*amc*8/Con*DY*.root 
hadd -f Arranged_DY/DY_9.root Out_DY*amc*9/Con*DY*.root 
hadd -f Arranged_DY/DY.root Arranged_DY/DY_*.root


mkdir -p Arranged_TT
hadd -f Arranged_TT/TT_had.root Out_TT_had*/Con*TT*.root 
hadd -f Arranged_TT/TT_semilep.root Out_TT_semilep*/Con*TT*.root 
hadd -f Arranged_TT/TT_2l2nu.root Out_TT_2l*/Con*TT*.root 


mkdir -p Arranged_WJets
hadd -f Arranged_WJets/WJets.root Out_WJets*/Con*WJets*.root 


#mkdir -p Arranged_QCD
#hadd -f Arranged_QCD/QCD_20to30.root Out_QCD_*20to30*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_30to50.root Out_QCD_*30to50*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_50to80.root Out_QCD_*50to80*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_80to120.root Out_QCD_*80to120*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_120to170.root Out_QCD_*120to170*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_170to300.root Out_QCD_*170to300*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_300to470.root Out_QCD_*300to470*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_470to600.root Out_QCD_*470to600*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_600to800.root Out_QCD_*600to800*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_800to1000.root Out_QCD_*800to1000*/Con*QCD*.root
#hadd -f Arranged_QCD/QCD_1000toInf.root Out_QCD_*1000toInf*/Con*QCD*.root

mkdir -p Arranged_WW
hadd -f Arranged_WW/WW.root Out_WW/Con*WW*.root 


mkdir -p Arranged_ZZ
hadd -f Arranged_ZZ/ZZ.root Out_ZZ/Con*ZZ*.root 


mkdir -p Arranged_WZ
hadd -f Arranged_WZ/WZ.root Out_WZ/Con*WZ*.root 


#mkdir -p Arranged_ST_top
#hadd -f Arranged_ST_top/ST_top.root Out_ST_top/Con*ST_top*.root 
#
#
#mkdir -p Arranged_ST_antitop
#hadd -f Arranged_ST_antitop/ST_antitop.root Out_ST_antitop/Con*ST_antitop*.root 


mkdir -p Arranged_data
hadd -f Arranged_data/dataB.root Out_SMu_*B*/Con*.root
hadd -f Arranged_data/dataC.root Out_SMu_*C*/Con*.root
hadd -f Arranged_data/dataD.root Out_SMu_*D*/Con*.root
hadd -f Arranged_data/dataE.root Out_SMu_*E*/Con*.root
hadd -f Arranged_data/dataF_0.root Out_SMu_*F_0*/Con*.root
hadd -f Arranged_data/dataF_1.root Out_SMu_*F_1*/Con*.root
hadd -f Arranged_data/dataF_2.root Out_SMu_*F_2*/Con*.root
hadd -f Arranged_data/dataF_3.root Out_SMu_*F_3*/Con*.root
hadd -f Arranged_data/dataF_4.root Out_SMu_*F_4*/Con*.root
hadd -f Arranged_data/dataF_5.root Out_SMu_*F_5*/Con*.root
hadd -f Arranged_data/dataF_6.root Out_SMu_*F_6*/Con*.root
hadd -f Arranged_data/dataF_7.root Out_SMu_*F_7*/Con*.root
hadd -f Arranged_data/dataF_8.root Out_SMu_*F_8*/Con*.root
hadd -f Arranged_data/data.root Arranged_data/dataB.root Arranged_data/dataC.root Arranged_data/dataD.root Arranged_data/dataE.root Arranged_data/dataF_*.root
