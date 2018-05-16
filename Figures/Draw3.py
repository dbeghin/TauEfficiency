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
    lumi.AddText("2017, 41.9 fb^{-1} (13 TeV)")
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
adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)

#categories=["passOS","failOS"] 
#ncat=2

var=[]
'''
var.append("tau_MVA")
var.append("mu_pt")
var.append("mu_eta")
var.append("mu_phi")
var.append("ev_DRmutau")
var.append("ev_Mt_raw")
'''
tauID_cutoff=len(var)

#var.append("ev_Mt")
var.append("ev_Mvis")
'''
var.append("ev_Mtot")
var.append("tau_pt")
var.append("tau_eta")
var.append("tau_phi")
var.append("tau_ptLead")
var.append("ev_METmumass")
var.append("ev_MET")
var.append("ev_METphi")
var.append("ev_Nvertex")
'''
nvar=len(var)
print tauID_cutoff 
print "\n"
print nvar

photogenic_var=[]
'''
photogenic_var.append("#tau-ID MVA value")
photogenic_var.append("#mu p_{T} (GeV)")
photogenic_var.append("#mu #eta")
photogenic_var.append("#mu #phi")
photogenic_var.append("#DeltaR (#mu #tau)")
photogenic_var.append("m_{T} (#mu E_{T}^{miss}) (GeV)")
photogenic_var.append("m_{T} (#mu E_{T}^{miss}) (GeV)")
'''
photogenic_var.append("m_{vis} (GeV)")
'''
photogenic_var.append("m_{tot} (GeV)")
photogenic_var.append("#tau p_{T} (GeV)")
photogenic_var.append("#tau #eta")
photogenic_var.append("#tau #phi")
photogenic_var.append("Leading Charged Cand. p_{T} (GeV)")
photogenic_var.append("m (#mu E_{T}^{miss}) (GeV)")
photogenic_var.append("E_{T}^{miss} (GeV)")
photogenic_var.append("#phi (E_{T}^{miss})")
photogenic_var.append("N_{vertex}")
'''

tauIDs=[
#"cutbased_loose",  
#"cutbased_medium",  
#"cutbased_tight",   
#"MVA_veryloose",   
#"MVA_loose",        
"MVA_medium",      
#"MVA_tight",       
#"MVA_verytight",   
#"MVA_veryverytight",
#"MVAnew_veryloose",   
#"MVAnew_loose",        
#"MVAnew_medium",      
#"MVAnew_tight",       
#"MVAnew_verytight",   
#"MVAnew_veryverytight"
]

photogenic_tauIDs={
"cutbased_loose" :   "Cut-based loose",
"cutbased_medium" :   "Cut-based medium",
"cutbased_tight" :    "Cut-based tight",
"MVA_veryloose" :    "MVA very loose",
"MVA_loose" :         "MVA loose",
"MVA_medium" :       "MVA medium",
"MVA_tight" :        "MVA tight",
"MVA_verytight" :    "MVA very tight",
"MVA_veryverytight" : "MVA very very tight",
"MVAnew_veryloose" :    "MVA very loose",
"MVAnew_loose" :         "MVA loose",
"MVAnew_medium" :       "MVA medium",
"MVAnew_tight" :        "MVA tight",
"MVAnew_verytight" :    "MVA very tight",
"MVAnew_veryverytight" : "MVA very very tight"
}

ntauID=len(tauIDs)
categories=["pass","fail","zmumu"]#"pass_postfit","fail_postfit","zmm_prefit","zmm_postfit"]
ncat=len(categories) #6

prepost=["prefit", "postfit"]



for j in range (0,ntauID):
    file=ROOT.TFile("ztt_mt_shapes_"+tauIDs[j]+".root","r")

    for k in range (0,len(prepost)):
        for i in range (0,ncat):
            folder_in = ""
            if categories[i] == "pass" :
                folder_in = "pass_"+prepost[k]
            elif categories[i] == "fail" :
                folder_in = "fail_"+prepost[k]
            elif categories[i] == "zmumu" :
                folder_in = "zmumu_"+prepost[k]
            else:
                folder_in = "error"
            #print var_in
            Data=file.Get(folder_in+"/data_obs")
            QCD=file.Get(folder_in+"/QCD")
            W=file.Get(folder_in+"/WJets")
            TTS=file.Get(folder_in+"/TTS")
            TTB=file.Get(folder_in+"/TTB")
            VV=file.Get(folder_in+"/VV")
            DYB=file.Get(folder_in+"/DYB")
            #if i<4:
            #DYJ=file.Get(categories[i]).Get("DYJ")
            #DYB.Add(DYJ)
            DYS=file.Get(folder_in+"/DYS")

            W.Add(VV)


            Data.GetXaxis().SetTitle("")
            Data.GetXaxis().SetTitleSize(0)
            Data.GetXaxis().SetNdivisions(505)
            Data.GetYaxis().SetLabelFont(42)
            Data.GetYaxis().SetLabelOffset(0.01)
            Data.GetYaxis().SetLabelSize(0.06)
            Data.GetYaxis().SetTitleSize(0.075)
            Data.GetYaxis().SetTitleOffset(1.04)
            Data.SetTitle("")
            Data.GetYaxis().SetTitle("Events/bin")
        
            if i<4:
                QCD.SetFillColor(ROOT.TColor.GetColor("#ffccff"))
            W.SetFillColor(ROOT.TColor.GetColor("#de5a6a"))
            if i<2: TTS.SetFillColor(ROOT.TColor.GetColor("#9999cc"))
            TTB.SetFillColor(ROOT.TColor.GetColor("#99992c"))
            DYB.SetFillColor(ROOT.TColor.GetColor("#4496c8"))
            if i<2:
                DYS.SetFillColor(ROOT.TColor.GetColor("#ffcc66"))

            Data.SetMarkerStyle(20)
            Data.SetMarkerSize(1)
            if i<4:
                QCD.SetLineColor(1)
            W.SetLineColor(1)
            if i<2: TTS.SetLineColor(1)
            TTB.SetLineColor(1)
            if i<2:
                DYS.SetLineColor(1)
            DYB.SetLineColor(1)
            Data.SetLineColor(1)
            Data.SetLineWidth(2)

        #DYS.Scale(0.84) #FIXME

            stack=ROOT.THStack("stack","stack")
            if i<4:
                stack.Add(QCD)
            stack.Add(W)
            if i<2: stack.Add(TTS)
            stack.Add(TTB)
            stack.Add(DYB)
            if i<2:
                stack.Add(DYS)
            
            errorBand = W.Clone()
            if i<4:
                errorBand.Add(QCD)
            if i<2: errorBand.Add(TTS)
            errorBand.Add(TTB)
            errorBand.Add(DYB)
            if i<2:
                errorBand.Add(DYS)
            errorBand.SetMarkerSize(0)
            print new_idx
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
            #pad1.SetLogy()
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
            Data.SetMaximum(Data.GetMaximum()*1.7)#1.5)#FIXME
            Data.SetMinimum(0.1)
            Data.Draw("e")
            stack.Draw("histsame")
            errorBand.Draw("e2same")
            Data.Draw("esame")

            legende=make_legend()
            legende.AddEntry(Data,"Observed","elp")
            if i<2: 
                legende.AddEntry(DYS,"Z#rightarrow#tau_{#mu}#tau_{h}","f")
            legende.AddEntry(DYB,"DY others","f")
            if i<2: legende.AddEntry(TTS,"t#bar{t}+jets (real #tau)","f")
            legende.AddEntry(TTB,"t#bar{t}+jets (fake #tau)","f")
            legende.AddEntry(W,"Electroweak","f")
            if i<4:
                legende.AddEntry(QCD,"QCD multijet","f")
            legende.AddEntry(errorBand,"Uncertainty","f")
            legende.Draw()

            l1=add_lumi()
            l1.Draw("same")
            l2=add_CMS()
            l2.Draw("same")
            l3=add_Preliminary()
            l3.Draw("same")
        
            pad1.RedrawAxis()

            tauID  = ROOT.TPaveText(0.21, 0.52+0.013, 0.43, 0.70+0.155, "NDC")
            tauID.SetBorderSize(   0 )
            tauID.SetFillStyle(    0 )
            tauID.SetTextAlign(   12 )
            tauID.SetTextSize ( 0.06 )
            tauID.SetTextColor(    1 )
            tauID.SetTextFont (   41 )
            tauID.AddText(photogenic_tauIDs[tauIDs[j]])
            tauID.Draw("same")

            categ  = ROOT.TPaveText(0.21, 0.45+0.013, 0.43, 0.65+0.155, "NDC")
            categ.SetBorderSize(   0 )
            categ.SetFillStyle(    0 )
            categ.SetTextAlign(   12 )
            categ.SetTextSize ( 0.06 )
            categ.SetTextColor(    1 )
            categ.SetTextFont (   41 )
            categ.AddText(categories[i])
            #categ.AddText("Z#rightarrow#mu#mu CR")
            categ.Draw("same")

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
            h1.SetMaximum(1.7)#FIXME(1.5)
            h1.SetMinimum(0.5)#FIXME(0.5)
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
            h1.GetXaxis().SetTitle(photogenic_var[0])
            h1.GetXaxis().SetLabelSize(0.08)
            h1.GetYaxis().SetLabelSize(0.08)
            h1.GetYaxis().SetTitle("Obs./Exp.")
            h1.GetXaxis().SetNdivisions(505)
            h1.GetYaxis().SetNdivisions(10)

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
            c.SaveAs("Mvis_"+tauIDs[j]+"_"+categories[i]+"_"+prepost[k]+".pdf")


