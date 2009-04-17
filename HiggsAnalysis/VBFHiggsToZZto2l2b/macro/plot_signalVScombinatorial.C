{

    // general root setting
  gROOT->Reset(); 
  //  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);

  TH1D * histo[33];
  histo[0]  = (TH1D*)Dphiqq->Clone();
  histo[1]  = (TH1D*)Detaqq->Clone();
  histo[2]  = (TH1D*)Ptminqq->Clone();
  histo[3]  = (TH1D*)Ptmaxqq->Clone();
  histo[4]  = (TH1D*)Etaminqq->Clone();
  histo[5]  = (TH1D*)Etamaxqq->Clone();
  histo[6]  = (TH1D*)Dphihh->Clone();
  histo[7]  = (TH1D*)Detahh->Clone();
  histo[8]  = (TH1D*)Ptminhh->Clone();
  histo[9]  = (TH1D*)Ptmaxhh->Clone();
  histo[10] = (TH1D*)Etaminhh->Clone();
  histo[11] = (TH1D*)Etamaxhh->Clone();
  histo[12] = (TH1D*)Ptqq->Clone();
  histo[13] = (TH1D*)Mqq->Clone();
  histo[14] = (TH1D*)Etaqq->Clone();
  histo[15] = (TH1D*)Pthh->Clone();
  histo[16] = (TH1D*)Mhh->Clone();
  histo[17] = (TH1D*)Etahh->Clone();
  histo[18] = (TH1D*)Ptll->Clone();
  histo[19] = (TH1D*)Mll->Clone();
  histo[20] = (TH1D*)Etall->Clone();
  histo[21] = (TH1D*)DphiTjetZjet->Clone();
  histo[22] = (TH1D*)DphiTjetZlep->Clone();
  histo[23] = (TH1D*)DphiminTZ->Clone();
  histo[24] = (TH1D*)DetaTjetZjet->Clone();
  histo[25] = (TH1D*)DetaTjetZlep->Clone();
  histo[26] = (TH1D*)DetaminTZ->Clone();
  histo[27] = (TH1D*)MassTjetZjet->Clone();
  histo[28] = (TH1D*)MassTjetZlep->Clone();
  histo[29] = (TH1D*)MassZjetZlep->Clone();
  histo[30] = (TH1D*)MassTZZ->Clone();
  histo[31] = (TH1D*)EtaTZZ->Clone();
  histo[32] = (TH1D*)PtTZZ->Clone();

  for (int index = 0; index < 33; index++ ) {
    histo[index]->SetLineColor(kRed);
    histo[index]->SetFillColor(kRed);
  }
  
  TH1D * f_histo[33];
  f_histo[0]  = (TH1D*)F_Dphiqq->Clone();
  f_histo[1]  = (TH1D*)F_Detaqq->Clone();
  f_histo[2]  = (TH1D*)F_Ptminqq->Clone();
  f_histo[3]  = (TH1D*)F_Ptmaxqq->Clone();
  f_histo[4]  = (TH1D*)F_Etaminqq->Clone();
  f_histo[5]  = (TH1D*)F_Etamaxqq->Clone();
  f_histo[6]  = (TH1D*)F_Dphihh->Clone();
  f_histo[7]  = (TH1D*)F_Detahh->Clone();
  f_histo[8]  = (TH1D*)F_Ptminhh->Clone();
  f_histo[9]  = (TH1D*)F_Ptmaxhh->Clone();
  f_histo[10] = (TH1D*)F_Etaminhh->Clone();
  f_histo[11] = (TH1D*)F_Etamaxhh->Clone();
  f_histo[12] = (TH1D*)F_Ptqq->Clone();
  f_histo[13] = (TH1D*)F_Mqq->Clone();
  f_histo[14] = (TH1D*)F_Etaqq->Clone();
  f_histo[15] = (TH1D*)F_Pthh->Clone();
  f_histo[16] = (TH1D*)F_Mhh->Clone();
  f_histo[17] = (TH1D*)F_Etahh->Clone();
  f_histo[18] = (TH1D*)F_Ptll->Clone();
  f_histo[19] = (TH1D*)F_Mll->Clone();
  f_histo[20] = (TH1D*)F_Etall->Clone();
  f_histo[21] = (TH1D*)F_DphiTjetZjet->Clone();
  f_histo[22] = (TH1D*)F_DphiTjetZlep->Clone();
  f_histo[23] = (TH1D*)F_DphiminTZ->Clone();
  f_histo[24] = (TH1D*)F_DetaTjetZjet->Clone();
  f_histo[25] = (TH1D*)F_DetaTjetZlep->Clone();
  f_histo[26] = (TH1D*)F_DetaminTZ->Clone();
  f_histo[27] = (TH1D*)F_MassTjetZjet->Clone();
  f_histo[28] = (TH1D*)F_MassTjetZlep->Clone();
  f_histo[29] = (TH1D*)F_MassZjetZlep->Clone();
  f_histo[30] = (TH1D*)F_MassTZZ->Clone();
  f_histo[31] = (TH1D*)F_EtaTZZ->Clone();
  f_histo[32] = (TH1D*)F_PtTZZ->Clone();


  TCanvas * canvas[33];
  for (int index =0; index < 33; index++ ) {
    //  for (int index =0; index < 2; index++ ) {
    char n[50];
    sprintf(n,"canvas_%d",index);
    canvas[index] = new TCanvas(n,"canvas",600,600);
    canvas[index]->Divide(2,2);
    canvas[index]->cd(1);
    histo[index]->Draw();
    canvas[index]->cd(2);
    f_histo[index]->Draw();
    canvas[index]->cd(3);
    histo[index]->SetNormFactor(1.);
    f_histo[index]->SetNormFactor(1.);
    f_histo[index]->Draw();
    histo[index]->Draw("same");
  }


}
