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
    #lumi.AddText("2017C, 9.8 fb^{-1} (13 TeV)")
    #lumi.AddText("2017E, 9.4 fb^{-1} (13 TeV)")
    #lumi.AddText("2017F, 13.6 fb^{-1} (13 TeV)")
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
file=ROOT.TFile("histos_mumu_80mb.root","r")
file_data=ROOT.TFile("histos_mumu_withOldData.root","r")


adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)

#categories=["passOS","failOS"] 
#ncat=2

var=[]
var.append("mu_pt")
var.append("mu_eta")
var.append("mu_phi")
#var.append("mu1_pt")
#var.append("mu1_eta")
#var.append("mu1_phi")
#var.append("mu2_pt")
#var.append("mu2_eta")
#var.append("mu2_phi")
#var.append("ev_DRmumu")
#var.append("ev_Mt_raw")
var.append("ev_Mt")
var.append("ev_Mvis")
#var.append("ev_Mvis_SS")
'''
var.append("ev_METmumass")
var.append("ev_MET")
var.append("ev_METphi")
'''
var.append("ev_Nvertex")
nvar=len(var)
print nvar

photogenic_var=[]
photogenic_var.append("#mu p_{T} (GeV)")
photogenic_var.append("#mu #eta")
photogenic_var.append("#mu #phi")
#photogenic_var.append("#mu_{1} p_{T} (GeV)")
#photogenic_var.append("#mu_{1} #eta")
#photogenic_var.append("#mu_{1} #phi")
#photogenic_var.append("#mu_{2} p_{T} (GeV)")
#photogenic_var.append("#mu_{2} #eta")
#photogenic_var.append("#mu_{2} #phi")
#photogenic_var.append("#DeltaR (#mu #mu)")
#photogenic_var.append("m_{T} (#mu E_{T}^{miss}) (GeV)")
photogenic_var.append("m_{T} (#mu E_{T}^{miss}) (GeV)")
photogenic_var.append("m_{vis} (GeV)")
#photogenic_var.append("m_{vis} (SS) (GeV)")
'''
photogenic_var.append("m (#mu E_{T}^{miss}) (GeV)")
photogenic_var.append("E_{T}^{miss} (GeV)")
photogenic_var.append("#phi (E_{T}^{miss})")
'''
photogenic_var.append("N_{vertex}")

for k in range (0,nvar):
    var_in = var[k]
    #print var_in
    Data=file_data.Get("data_"+var_in)
    W=file.Get("WJets_"+var_in)
    TT=file.Get("TTB_"+var_in)
    VV=file.Get("VV_"+var_in)#VV
    #W.Add(VV)
    DYB=file.Get("DYB_"+var_in)
    #if i<4:
    #DYJ=file.Get(categories[i]).Get("DYJ")
    #DYB.Add(DYJ)
    DYS=file.Get("DYS_"+var_in)

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
    
    
    #W.SetFillColor(ROOT.TColor.GetColor("#de5a6a"))
    TT.SetFillColor(ROOT.TColor.GetColor("#9999cc"))
    DYB.SetFillColor(ROOT.TColor.GetColor("#4496c8"))
    
    DYS.SetFillColor(ROOT.TColor.GetColor("#ffcc66"))

    Data.SetMarkerStyle(20)
    Data.SetMarkerSize(1)
    
    #W.SetLineColor(1)
    TT.SetLineColor(1)
    
    DYS.SetLineColor(1)
    DYB.SetLineColor(1)
    Data.SetLineColor(1)
    Data.SetLineWidth(2)

    #DYS.Scale(0.84) #FIXME

    stack=ROOT.THStack("stack","stack")
    
    #stack.Add(W)
    stack.Add(TT)
    stack.Add(DYS)
    stack.Add(DYB)
    
    
    errorBand = TT.Clone()
    
    #errorBand.Add(W)
    
    errorBand.Add(DYS)
    errorBand.Add(DYB)
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
    if ("Mvis" in var[k]): Data.GetXaxis().SetRangeUser(60,120)
    Data.SetMaximum(Data.GetMaximum()*2)
    Data.SetMinimum(0)
    Data.Draw("esame")
    stack.Draw("histsame")
    errorBand.Draw("e2same")
    Data.Draw("esame")

    legende=make_legend()
    legende.AddEntry(Data,"Observed","elp")
    
    legende.AddEntry(DYS,"Z#rightarrow#tau_{#mu}#tau_{h}","f")
    legende.AddEntry(DYB,"DY others","f")
    legende.AddEntry(TT,"t#bar{t}+jets","f")
    #legende.AddEntry(W,"Electroweak","f")
    
    legende.AddEntry(errorBand,"Uncertainty","f")
    legende.Draw()

    l1=add_lumi()
    l1.Draw("same")
    l2=add_CMS()
    l2.Draw("same")
    l3=add_Preliminary()
    l3.Draw("same")
    
    pad1.RedrawAxis()

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
    h1.GetXaxis().SetTitle(photogenic_var[k])
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
    c.SaveAs(var_in+".pdf")


