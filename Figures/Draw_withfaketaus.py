#!/usr/bin/env python
import ROOT
import re
from array import array

def add_lumi():
    lowX=0.58
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.06)
    lumi.SetTextFont (   42 )
    lumi.AddText("2017, 41.5 fb^{-1} (13 TeV)")
    return lumi

def add_CMS():
    lowX=0.21
    lowY=0.70
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.08)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi

def add_Preliminary():
    lowX=0.21
    lowY=0.63
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(52)
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("Preliminary")
    return lumi

def make_legend():
    output = ROOT.TLegend(0.65, 0.4, 0.92, 0.82, "", "brNDC")
    output.SetLineWidth(0)
    output.SetLineStyle(0)
    output.SetFillStyle(0)
    output.SetBorderSize(0)
    output.SetTextFont(62)
    return output

ROOT.gStyle.SetFrameLineWidth(3)
ROOT.gStyle.SetLineWidth(3)
ROOT.gStyle.SetOptStat(0)

c=ROOT.TCanvas("canvas","",0,0,600,600)
c.cd()

#file=ROOT.TFile("final.root","r")
file=ROOT.TFile("allhistos_control.root","r")

adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)

#categories=["passOS","failOS"] 
#ncat=2

var=[]
#var.append("tau_MVA")  
#var.append("ev_DRmutau")
#var.append("ev_Mt_raw")
tauID_cutoff=len(var)
var.append("ev_Mt")                
var.append("ev_Mvis")              
#var.append("ev_Mvis_TESup")        
#var.append("ev_Mvis_TESdown")      
#var.append("ev_Mvis_MESup")        
#var.append("ev_Mvis_MESdown")      
var.append("ev_Mvis_MinBiasdown") 
var.append("ev_Mvis_MinBiasup")   
#var.append("ev_Mvis_FRS_DM0_down")
#var.append("ev_Mvis_FRS_DM0_up")  
#var.append("ev_Mvis_FRS_DM1_down")
#var.append("ev_Mvis_FRS_DM1_up")  
#var.append("ev_Mvis_FRS_DM10_down")
#var.append("ev_Mvis_FRS_DM10_up") 
#var.append("ev_Mvis_antiisomu_down")
#var.append("ev_Mvis_antiisomu_up")
#var.append("ev_Mvis_antiisotau_down")
#var.append("ev_Mvis_antiisotau_up")
#var.append("tau_pt_TESup")         
#var.append("tau_pt_TESdown")       
#var.append("tau_pt_MESup")         
#var.append("tau_pt_MESdown")       
#var.append("tau_pt_FRS_DM0_down") 
#var.append("tau_pt_FRS_DM0_up")   
#var.append("tau_pt_FRS_DM1_down") 
#var.append("tau_pt_FRS_DM1_up")   
#var.append("tau_pt_FRS_DM10_down")
#var.append("tau_pt_FRS_DM10_up")  
#var.append("tau_pt_antiisomu_down")
#var.append("tau_pt_antiisomu_up") 
#var.append("tau_pt_antiisotau_down")
#var.append("tau_pt_antiisotau_up")
#var.append("tau_pt")               
#var.append("tau_eta")              
#var.append("tau_phi")              
#var.append("mu_pt")                
#var.append("mu_eta")               
#var.append("mu_phi")               
#var.append("ev_Nvertex")           
#var.append("ev_Nvertex_MinBiasdown")
#var.append("ev_Nvertex_Minbiasup")

nvar=len(var)
print tauID_cutoff
print "\n"
print nvar

photogenic_var={
"tau_MVA":  "#tau-ID MVA value",
"ev_DRmutau": "#DeltaR (#mu #tau)",
"ev_Mt_raw": "m_{T} (#mu E_{T}^{miss}) (GeV)",
"ev_Mt":                   "m_{T} (#mu E_{T}^{miss}) (GeV)",
"ev_Mvis":                 "m_{vis} (GeV)",
"ev_Mvis_TESup":           "m_{vis} (GeV)",
"ev_Mvis_TESdown":         "m_{vis} (GeV)",
"ev_Mvis_MESup":           "m_{vis} (GeV)",
"ev_Mvis_MESdown":         "m_{vis} (GeV)",
"ev_Mvis_MinBiasdown":     "m_{vis} (GeV)",
"ev_Mvis_MinBiasup":       "m_{vis} (GeV)",
"ev_Mvis_FRS_DM0_down":    "m_{vis} (GeV)",
"ev_Mvis_FRS_DM0_up":      "m_{vis} (GeV)",
"ev_Mvis_FRS_DM1_down":    "m_{vis} (GeV)",
"ev_Mvis_FRS_DM1_up":      "m_{vis} (GeV)",
"ev_Mvis_FRS_DM10_down":   "m_{vis} (GeV)",
"ev_Mvis_FRS_DM10_up":     "m_{vis} (GeV)",
"ev_Mvis_antiisomu_down":  "m_{vis} (GeV)",
"ev_Mvis_antiisomu_up":    "m_{vis} (GeV)",
"ev_Mvis_antiisotau_down": "m_{vis} (GeV)",
"ev_Mvis_antiisotau_up":   "m_{vis} (GeV)",
"tau_pt_TESup":            "#tau p_{T} (GeV)",
"tau_pt_TESdown":          "#tau p_{T} (GeV)",
"tau_pt_MESup":            "#tau p_{T} (GeV)",
"tau_pt_MESdown":          "#tau p_{T} (GeV)",
"tau_pt_FRS_DM0_down":     "#tau p_{T} (GeV)",
"tau_pt_FRS_DM0_up":       "#tau p_{T} (GeV)",
"tau_pt_FRS_DM1_down":     "#tau p_{T} (GeV)",
"tau_pt_FRS_DM1_up":       "#tau p_{T} (GeV)",
"tau_pt_FRS_DM10_down":    "#tau p_{T} (GeV)",
"tau_pt_FRS_DM10_up":      "#tau p_{T} (GeV)",
"tau_pt_antiisomu_down":   "#tau p_{T} (GeV)",
"tau_pt_antiisomu_up":     "#tau p_{T} (GeV)",
"tau_pt_antiisotau_down":  "#tau p_{T} (GeV)",
"tau_pt_antiisotau_up":    "#tau p_{T} (GeV)",
"tau_pt":                  "#tau p_{T} (GeV)",
"tau_eta":                 "#tau #eta",
"tau_phi":                 "#tau #phi",
"mu_pt":                   "#mu p_{T} (GeV)",
"mu_eta":                  "#mu #eta",
"mu_phi":                  "#mu #phi",
"ev_Nvertex":              "ev_Nvertex",
"ev_Nvertex_MinBiasdown":  "ev_Nvertex",
"ev_Nvertex_Minbiasup":    "ev_Nvertex",
}

tauIDs=[
#"cutbased_loose",   
#"cutbased_medium",  
#"cutbased_tight",   
#                    
#"MVA_2017v2vvloose",
#"MVA_2017v2vloose", 
#"MVA_2017v2loose",  
#"MVA_2017v2medium", 
"MVA_2017v2tight",  
#"MVA_2017v2vtight", 
#"MVA_2017v2vvtight",
]

photogenic_tauIDs={
"cutbased_loose":   "Cut-based loose",
"cutbased_medium":  "Cut-based medium",
"cutbased_tight":   "Cut-based tight",
                    
"MVA_2017v1vvloose":"MVA 2017v1 v. v. loose",
"MVA_2017v1vloose": "MVA 2017v1 v. loose",
"MVA_2017v1loose":  "MVA 2017v1 loose",
"MVA_2017v1medium": "MVA 2017v1 medium",
"MVA_2017v1tight":  "MVA 2017v1 tight",
"MVA_2017v1vtight": "MVA 2017v1 v. tight",
"MVA_2017v1vvtight":"MVA 2017v1 v. v. tight",
                    
"MVA_2017v2vvloose":"MVA 2017v2 v. v. loose",
"MVA_2017v2vloose": "MVA 2017v2 v. loose",
"MVA_2017v2loose":  "MVA 2017v2 loose",
"MVA_2017v2medium": "MVA 2017v2 medium",
"MVA_2017v2tight":  "MVA 2017v2 tight",
"MVA_2017v2vtight": "MVA 2017v2 v. tight",
"MVA_2017v2vvtight":"MVA 2017v2 v. v. tight",
}

dms=[
"allDMs",
"DM0",
"DM1",
"DM10",
]

eta=[
"alleta",
"barrel",
"endcap",
]

ptrange=[
"allpt",
"pt_20_40",
"pt_40_150",
]


ntauID=len(tauIDs)
categories=["pass",
#"fail"]#"pass_postfit","fail_postfit","zmm_prefit","zmm_postfit"
]

ncat=len(categories)                                                                                                                                                                                                                                                                                                                                                                         

for i in range (0,nvar):
    for j in range (0,len(dms)):
        for k in range (0,len(eta)):
            for l in range (0,ntauID):
                for m in range (0,ncat):
                    for p in range (0,len(ptrange)):
                        if i < tauID_cutoff:# or j == 0:
                            if j>0 or k>0 or l>0 or m>0:
                                continue
                            var_in = var[i]
                            ph_tauID_in = "No tau ID"
                            categ_in = ""
                        else:
                            var_in = var[i]+"_"+dms[j]+"_"+eta[k]+"_"+ptrange[p]+"_"+tauIDs[l]+"_"+categories[m]
                            ph_tauID_in = photogenic_tauIDs[tauIDs[l]]
                            categ_in = categories[m]
                        print var_in
                        ph_var_in = photogenic_var[var[i]]
                        
                        Data=file.Get(tauIDs[l]+"/data_"+var_in)
                        TT=file.Get(tauIDs[l]+"/TTB_"+var_in)
                        VV=file.Get(tauIDs[l]+"/VV_"+var_in)
                        DYS=file.Get(tauIDs[l]+"/DYS_"+var_in)
                        DYB=file.Get(tauIDs[l]+"/DYB_"+var_in)
                        Faketau=file.Get(tauIDs[l]+"/faketau_"+var_in)
                        
                        try:
                            Data.GetXaxis().SetTitle("")
                        except:
                            print "non existent file:  "+ tauIDs[l]+"/data_"+var_in
                            continue
                        Data.GetXaxis().SetTitleSize(0)
                        Data.GetXaxis().SetNdivisions(505)
                        Data.GetYaxis().SetLabelFont(42)
                        Data.GetYaxis().SetLabelOffset(0.01)
                        Data.GetYaxis().SetLabelSize(0.06)
                        Data.GetYaxis().SetTitleSize(0.075)
                        Data.GetYaxis().SetTitleOffset(1.04)
                        Data.SetTitle("")
                        Data.GetYaxis().SetTitle("Events/bin")
                        
                        
                        #Write number of events
                        if k == nvar-1:
                            nevents_file = open("nevents_"+var_in+".txt", 'w')
                            nevents_file.write("processus | nevents \n")
                            nevents_file.write("VV          " + str(VV.Integral()) + "\n")
                            nevents_file.write("TT          " + str(TT.Integral()) + "\n")
                            nevents_file.write("data        " + str(Data.Integral()) + "\n")
                            nevents_file.close()
                        
                        
                        
                        #QCD.SetFillColor(ROOT.TColor.GetColor("#ffccff"))
                        #W.SetFillColor(ROOT.TColor.GetColor("#de5a6a"))
                        VV.SetFillColor(ROOT.TColor.GetColor("#d89a6a"))
                        TT.SetFillColor(ROOT.TColor.GetColor("#9999cc"))
                        DYB.SetFillColor(ROOT.TColor.GetColor("#4496c8"))
                        DYS.SetFillColor(ROOT.TColor.GetColor("#ffcc66"))
                        Faketau.SetFillColor(ROOT.TColor.GetColor("#de5a6a"))
                        
                        Data.SetMarkerStyle(20)
                        Data.SetMarkerSize(1)
                        VV.SetLineColor(1)
                        TT.SetLineColor(1)
                        DYB.SetLineColor(1)
                        DYS.SetLineColor(1)
                        Faketau.SetLineColor(1)
                        Data.SetLineColor(1)
                        Data.SetLineWidth(2)
                    
                    
                        stack=ROOT.THStack("stack","stack")
                        #stack.Add(W)
                        stack.Add(VV)
                        stack.Add(TT)
                        stack.Add(Faketau)
                        stack.Add(DYB)
                        stack.Add(DYS)
                        
                        errorBand = TT.Clone()
                        #errorBand.Add(QCD)
                        errorBand.Add(DYB)
                        errorBand.Add(DYS)
                        errorBand.Add(VV)
                        errorBand.Add(Faketau)
                        errorBand.SetMarkerSize(0)
                        errorBand.SetFillColor(new_idx)
                        errorBand.SetFillStyle(3001)
                        errorBand.SetLineWidth(1)
                        
                        pad1 = ROOT.TPad("pad1","pad1",0,0.35,1,1)
                        pad1.Draw()
                        pad1.cd()
                        pad1.SetFillColor(0)
                        pad1.SetBorderMode(0)
                        pad1.SetBorderSize(10)
                        pad1.SetTickx(1)
                        pad1.SetTicky(1)
                        pad1.SetLeftMargin(0.18)
                        pad1.SetRightMargin(0.05)
                        pad1.SetTopMargin(0.122)
                        pad1.SetBottomMargin(0.026)
                        pad1.SetFrameFillStyle(0)
                        pad1.SetFrameLineStyle(0)
                        pad1.SetFrameLineWidth(3)
                        pad1.SetFrameBorderMode(0)
                        pad1.SetFrameBorderSize(10)
                        
                        Data.GetXaxis().SetLabelSize(0)
                        Data.SetMaximum(Data.GetMaximum()*2.5)#2.5)#FIXME
                        Data.SetMinimum(0.1)
                        Data.Draw("e")
                        stack.Draw("histsame")
                        errorBand.Draw("e2same")
                        Data.Draw("esame")
                    
                    
                        
                        legende=make_legend()
                        legende.AddEntry(Data,"Observed","elp")
                        legende.AddEntry(DYB,"Z#rightarrow#mu #mu","f")
                        legende.AddEntry(DYS,"Z#rightarrow#tau_#mu #tau_h","f")
                        legende.AddEntry(Faketau,"Fake #tau bg","f")
                        legende.AddEntry(TT,"t#bar{t}+jets","f")
                        #legende.AddEntry(W,"W+jets","f")
                        legende.AddEntry(VV,"Diboson","f")
                        #legende.AddEntry(QCD,"QCD multijet","f")
                        legende.AddEntry(errorBand,"Uncertainty","f")
                        legende.Draw()
                        
                        l1=add_lumi()
                        l1.Draw("same")
                        l2=add_CMS()
                        l2.Draw("same")
                        l3=add_Preliminary()
                        l3.Draw("same")
                        
                        pad1.RedrawAxis()
                        
                        finalstate  = ROOT.TLegend(0.21, 0.52+0.013, 0.43, 0.70+0.155)
                        finalstate.SetBorderSize(   0 )
                        finalstate.SetFillStyle(    0 )
                        finalstate.SetTextAlign(   12 )
                        finalstate.SetTextSize ( 0.06 )
                        finalstate.SetTextColor(    1 )
                        #finalstate.SetTextFont (   41 )
                        finalstate.SetHeader("#mu #tau")
                        finalstate.Draw("same")
                        
                        '''
                        categ  = ROOT.TPaveText(0.21, 0.45+0.013, 0.43, 0.65+0.155, "NDC")
                        categ.SetBorderSize(   0 )
                        categ.SetFillStyle(    0 )
                        categ.SetTextAlign(   12 )
                        categ.SetTextSize ( 0.06 )
                        categ.SetTextColor(    1 )
                        categ.SetTextFont (   41 )
                        categ.AddText("OS iso #mu anti-iso #tau")
                        #categ.AddText("Z#rightarrow#mu#mu CR")
                        categ.Draw("same")
                        '''
                        
                        tauID  = ROOT.TPaveText(0.21, 0.45+0.013, 0.43, 0.65+0.155, "NDC")
                        tauID.SetBorderSize(   0 )
                        tauID.SetFillStyle(    0 )
                        tauID.SetTextAlign(   12 )
                        tauID.SetTextSize ( 0.06 )
                        tauID.SetTextColor(    1 )
                        tauID.SetTextFont (   41 )
                        tauID.AddText(ph_tauID_in)
                        tauID.Draw("same")
                        
                        
                        c.cd()
                        pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.35)
                        pad2.SetTopMargin(0.05)
                        pad2.SetBottomMargin(0.35)
                        pad2.SetLeftMargin(0.18)
                        pad2.SetRightMargin(0.05)
                        pad2.SetTickx(1)
                        pad2.SetTicky(1)
                        pad2.SetFrameLineWidth(3)
                        pad2.SetGridx()
                        pad2.SetGridy()
                        pad2.Draw()
                        pad2.cd()
                        h1=Data.Clone()
                        h1.SetMaximum(1.3)#FIXME(1.5)
                        h1.SetMinimum(0.7)#FIXME(0.5)
                        h1.SetMarkerStyle(20)
                        h3=errorBand.Clone()
                        hwoE=errorBand.Clone()
                        for iii in range (1,hwoE.GetSize()-2):
                            hwoE.SetBinError(iii,0)
                        h3.Sumw2()
                        h1.Sumw2()
                        h1.SetStats(0)
                        h1.Divide(hwoE)
                        h3.Divide(hwoE)
                        h1.GetXaxis().SetTitle(ph_var_in)
                        h1.GetXaxis().SetLabelSize(0.08)
                        h1.GetYaxis().SetLabelSize(0.08)
                        h1.GetYaxis().SetTitle("Obs./Exp.")
                        h1.GetXaxis().SetNdivisions(505)
                        h1.GetYaxis().SetNdivisions(6)
                        
                        h1.GetXaxis().SetTitleSize(0.15)
                        h1.GetYaxis().SetTitleSize(0.15)
                        h1.GetYaxis().SetTitleOffset(0.56)
                        h1.GetXaxis().SetTitleOffset(1.04)
                        h1.GetXaxis().SetLabelSize(0.11)
                        h1.GetYaxis().SetLabelSize(0.11)
                        h1.GetXaxis().SetTitleFont(42)
                        h1.GetYaxis().SetTitleFont(42)
                        
                        h1.Draw("ep")
                        h3.Draw("e2same")
                        
                        c.cd()
                        pad1.Draw()
                        
                        ROOT.gPad.RedrawAxis()
                        
                        c.Modified()
                        c.SaveAs(var_in+".png")


