#ifndef __PDFTOOL_H__
#define __PDFTOOL_H__

//// use case
//
// 1. load in interactive
//   export LHAPATH=/afs/cern.ch/atlas/software/releases/16.0.2/sw/lcg/external/MCGenerators/lhapdf/5.8.2/share/PDFsets
//   root [0] gSystem->Load("libLHAPDF.so")
//   root [1] .L PDFTool.h+
//   // real name (SONAME) of shared library is libLHAPDF.so.0, so need the symlink and path(LD_LIBRARY_PATH) to it.
//
// 1'. or compile
//   g++ test.cc PDFTool.h libLHAPDF.a -lgfortran -m32 (or just make)
//   export LHAPATH=/afs/cern.ch/atlas/software/releases/16.0.2/sw/lcg/external/MCGenerators/lhapdf/5.8.2/share/PDFsets
//
// 2. scale beam energy from 3.5TeV to 4.0TeV
//   PDFTool pdfTool(3500000, 4/3.5);
//   in event loop:
//     use event_weight(pow(mcevt_pdf_scale,2), mcevt_pdf_x1, mcevt_pdf_x2, mcevt_pdf_id1, mcevt_pdf_id2, mcevt_pdf1);
//
// 2'. or change pdfset to 10042(cteq6ll)
//   PDFTool pdfTool(3500000, 1, 10042);
//   in event loop:
//     use event_weight(pow(mcevt_pdf_scale,2), mcevt_pdf_x1, mcevt_pdf_x2, mcevt_pdf_id1, mcevt_pdf_id2, mcevt_pdf1);
//
// 3. loop over many beam energy
//   PDFTool pdfTool(3500000,1,-1,mcevt_pdf1);// init_pdf can be set later by pdfTool.setPdfset(mcevt_pdf1)
//   in event loop:
//     pdfTool.setEventInfo(pow(mcevt_pdf_scale,2), mcevt_pdf_x1, mcevt_pdf_x2, mcevt_pdf_id1, mcevt_pdf_id2);
//     for(int i=0; i<8; i++) use pdfTool.scale((3.5+i*0.5)/3.5);//3.5TeV~7TeV
//
// 3'. or loop over many pdfsets
//   PDFTool pdfTool(3500000,1,-1,mcevt_pdf1);// init_pdf can be set later by pdfTool.setPdfset(mcevt_pdf1)
//   in event loop:
//     pdfTool.setEventInfo(pow(mcevt_pdf_scale,2), mcevt_pdf_x1, mcevt_pdf_x2, mcevt_pdf_id1, mcevt_pdf_id2);
//     for(int i=0; i<41; i++) use pdfTool.pdf(10100+i);//cteq61
//
// 0. simple test:
//   source setup.sh
//   root -x test.C
//   // get 1.20286
//
// 0' compile version:
//   source setup.sh
//   make;
//   ./a.out
//   // also get 1.20286
//
//// contact: Yousuke.Kataoka@cern.ch
//// wiki:  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/PDFReweight

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
// ~~~~~~~~~
#include "LHAPDF/LHAPDF.h"
// ~~~~~~~~~

class PDFTool {
 public:
  PDFTool(double beam_energy, double reweightScale=1, int reweightPdf=-1, int init_pdfset=-1, int nmxset=60)//i.e. beam_energy=3500000, init_pdfset=20650 (can be set later in event_weight)
   : m_newvalues(true), m_q2(0), m_x1(0), m_x2(0), m_f1(0), m_f2(0), m_pdf1org(0.), m_pdf2org(0.) {
    m_energy = beam_energy;
    setPdfset(init_pdfset);
    m_nmxset = nmxset;
    m_reweightPdf = reweightPdf;
    m_reweightScale = reweightScale;
  }
  virtual ~PDFTool(){};

  // change beam energy
  double scale(double scale){
    return weight(-1,scale);
  }

  // change pdf
  double pdf(int pdf){
    return weight(pdf,1.);
  }

  double event_weight(double pdf_scale2, double pdf_x1, double pdf_x2, int pdf_id1, int pdf_id2, int init_pdf=-1){//Q^2,x1,x2,f1,f2,init_pdf
    if(m_reweightPdf==-1 && m_reweightScale==1) return 1;

    if(init_pdf!=-1) setPdfset(init_pdf);

    setEventInfo(pdf_scale2,pdf_x1,pdf_x2,pdf_id1,pdf_id2);
    double reweight = weight(m_reweightPdf,m_reweightScale);
    return reweight;
  }

  // ratio of pdf, pdf==-1 to use original one
  double weight(int pdf=-1, double scale=1.){
    if(pdf==-1) {//-1 to use original one
      if(m_pdfset==-1){
        std::cout << "initial pdfset must be given" << std::endl;
        abort();
      }
      pdf=m_pdfset;
    }

    double q = sqrt(getQ2());

    //std::cout << "x1=" << getX1() << " x2=" << getX2() << " f1=" << getF1() << " f2=" << getF2() << " q=" << q << std::endl;
    if (m_f1==0||m_f2==0||m_q2==0) return 0;

    if (m_newvalues) {
      // pdf of particle1 by original pdf
      m_pdf1org = testPdf(m_pdfset,m_x1,q,m_f1)/m_x1;
      // pdf of particle2 by original pdf
      m_pdf2org = testPdf(m_pdfset,m_x2,q,m_f2)/m_x2;
      m_newvalues = false;
    }
    if (m_pdf1org == 0 || m_pdf2org == 0) return 0;

    // pdf of particle1 by given pdf
    double pdf1new = testPdf(pdf,m_x1/scale,q,m_f1)/(m_x1/scale);

    // pdf of particle2 by given pdf
    double pdf2new = testPdf(pdf,m_x2/scale,q,m_f2)/(m_x2/scale);


    double weight = (pdf1new*pdf2new)/(m_pdf1org*m_pdf2org)/pow(scale,2);
    //std::cout << "pdf(p1)=" << m_pdf1org << " pdf(p2)=" << m_pdf2org << " -> pdf(p1)=" << pdf1new << " pdf(p2)=" << pdf2new << " weight=" << weight << std::endl;

    return weight;
  }

  //////////////////////////////////////////// PDF

  double testPdf(int pdfset, double x, double q, int pdgid){
    std::string pdfname = pdfName(pdfset);
    int member = pdfMember(pdfset);
    return testPdf(pdfname, member, x, q, pdgid);
  }

  double testPdf(std::string pdfname, int member, double x, double q, int pdgid){
    int nset = setPdf(pdfname,member);

    double f[13];
    // evolvepdfm_(nset,x,q,f);
// ~~~~~~~~~
    LHAPDF::xfx(nset, x, q, f);
// ~~~~~~~~~

    if(pdgid==1) return f[7];
    else if(pdgid==-1) return f[5];
    else if(pdgid==2) return f[8];
    else if(pdgid==-2) return f[4];
    else if(pdgid==3) return f[9];
    else if(pdgid==-3) return f[3];
    else if(pdgid==4) return f[10];
    else if(pdgid==-4) return f[2];
    else if(pdgid==5) return f[11];
    else if(pdgid==-5) return f[1];
    else if(pdgid==6) return f[12];
    else if(pdgid==-6) return f[0];
    else if(pdgid==21) return f[6];
    else {
      std::cout << "pdgid=" << pdgid << " is not defined in testPdf()" << std::endl;
      return 0;
    }
  }

  std::string pdfName(int pdfset){
    if(pdfset>=10000&&pdfset<=10040) return "cteq6.LHpdf";
    else if(pdfset>=10041&&pdfset<=10041) return "cteq6l.LHpdf";
    else if(pdfset>=10042&&pdfset<=10042) return "cteq6ll.LHpdf";
    else if(pdfset>=10050&&pdfset<=10090) return "cteq6mE.LHgrid";
    else if(pdfset>=10100&&pdfset<=10140) return "cteq61.LHpdf";
    else if(pdfset>=10150&&pdfset<=10190) return "cteq6l.LHgrid";
    else if(pdfset>=10350&&pdfset<=10390) return "cteq65.LHgrid";
    else if(pdfset>=10550&&pdfset<=10594) return "cteq66.LHgrid";
    else if(pdfset>=10770&&pdfset<=10770) return "CT09MCS.LHgrid";
    else if(pdfset>=10771&&pdfset<=10771) return "CT09MC2.LHgrid";
    else if(pdfset>=10772&&pdfset<=10772) return "CT09MC2.LHgrid";
    else if(pdfset>=10800&&pdfset<=10852) return "CT10.LHgrid";
    else if(pdfset>=10860&&pdfset<=10870) return "CT10as.LHgrid";
    else if(pdfset>=11000&&pdfset<=11052) return "CT10nlo.LHgrid";
    else if(pdfset>=20000&&pdfset<=20004) return "MRST2001nlo.LHpdf";
    else if(pdfset>=20050&&pdfset<=20054) return "MRST2001nlo.LHgrid";
    else if(pdfset>=20060&&pdfset<=20061) return "MRST2001lo.LHgrid";
    else if(pdfset>=20070&&pdfset<=20074) return "MRST2001nnlo.LHgrid";
    else if(pdfset>=20100&&pdfset<=20130) return "MRST2001E.LHpdf";
    else if(pdfset>=20150&&pdfset<=20180) return "MRST2001E.LHgrid";
    else if(pdfset>=20550&&pdfset<=20580) return "MRST2006nnlo.LHgrid";
    else if(pdfset>=20650&&pdfset<=20650) return "MRST2007lomod.LHgrid";
    else if(pdfset>=20651&&pdfset<=20651) return "MRSTMCal.LHgrid";
    else if(pdfset>=21000&&pdfset<=21040) return "MSTW2008lo68cl.LHgrid";
    else if(pdfset>=21041&&pdfset<=21081) return "MSTW2008lo90cl.LHgrid";
    else if(pdfset>=21100&&pdfset<=21140) return "MSTW2008nlo68cl.LHgrid";
    else if(pdfset>=21141&&pdfset<=21181) return "MSTW2008nlo90cl.LHgrid";
    else if(pdfset>=21200&&pdfset<=21240) return "MSTW2008nnlo68cl.LHgrid";
    else if(pdfset>=21241&&pdfset<=21281) return "MSTW2008nnlo90cl.LHgrid";
    else if(pdfset>=22000&&pdfset<=22021) return "MSTW2008nlo_asmzrange.LHgrid";
    else if(pdfset>=60500&&pdfset<=60520) return "HERAPDF10_EIG.LHgrid";
    else if(pdfset>=60530&&pdfset<=60543) return "HERAPDF10_VAR.LHgrid";
    else if(pdfset>=60550&&pdfset<=60561) return "HERAPDF10_ALPHAS.LHgrid";
    else if(pdfset>=40850&&pdfset<=40850) return "abkm09_5_nlo.LHgrid";
    else if(pdfset>=80260&&pdfset<=80286) return "GJR08VFnloE.LHgrid";
    else if(pdfset>=192800&&pdfset<=192900) return "NNPDF21_100.LHgrid";
    else {
      std::cout << "pdfset=" << pdfset << " is not defined in pdfName()" << std::endl;
      abort();
    }
  }
  int pdfMember(int pdfset){
    if(pdfset>=10000&&pdfset<=10040) return pdfset-10000;
    else if(pdfset>=10041&&pdfset<=10041) return pdfset-10041;
    else if(pdfset>=10042&&pdfset<=10042) return pdfset-10042;
    else if(pdfset>=10050&&pdfset<=10090) return pdfset-10050;
    else if(pdfset>=10100&&pdfset<=10140) return pdfset-10100;
    else if(pdfset>=10150&&pdfset<=10190) return pdfset-10150;
    else if(pdfset>=10350&&pdfset<=10390) return pdfset-10350;
    else if(pdfset>=10550&&pdfset<=10594) return pdfset-10550;
    else if(pdfset>=10770&&pdfset<=10770) return pdfset-10770;
    else if(pdfset>=10771&&pdfset<=10771) return pdfset-10771;
    else if(pdfset>=10772&&pdfset<=10772) return pdfset-10772;
    else if(pdfset>=10800&&pdfset<=10852) return pdfset-10800;
    else if(pdfset>=10860&&pdfset<=10870) return pdfset-10860;
    else if(pdfset>=11000&&pdfset<=11052) return pdfset-11000;
    else if(pdfset>=20000&&pdfset<=20004) return pdfset-20000;
    else if(pdfset>=20050&&pdfset<=20054) return pdfset-20050;
    else if(pdfset>=20060&&pdfset<=20061) return pdfset-20060;
    else if(pdfset>=20070&&pdfset<=20074) return pdfset-20070;
    else if(pdfset>=20100&&pdfset<=20130) return pdfset-20100;
    else if(pdfset>=20150&&pdfset<=20180) return pdfset-20150;
    else if(pdfset>=20550&&pdfset<=20580) return pdfset-20550;
    else if(pdfset>=20650&&pdfset<=20650) return pdfset-20650;
    else if(pdfset>=20651&&pdfset<=20651) return pdfset-20651;
    else if(pdfset>=21000&&pdfset<=21040) return pdfset-21000;
    else if(pdfset>=21041&&pdfset<=21081) return pdfset-21041;
    else if(pdfset>=21100&&pdfset<=21140) return pdfset-21100;
    else if(pdfset>=21141&&pdfset<=21181) return pdfset-21141;
    else if(pdfset>=21200&&pdfset<=21240) return pdfset-21200;
    else if(pdfset>=21241&&pdfset<=21281) return pdfset-21241;
    else if(pdfset>=22000&&pdfset<=22021) return pdfset-22000;
    else if(pdfset>=60500&&pdfset<=60520) return pdfset-60500;
    else if(pdfset>=60530&&pdfset<=60543) return pdfset-60530;
    else if(pdfset>=60550&&pdfset<=60561) return pdfset-60550;
    else if(pdfset>=40850&&pdfset<=40850) return pdfset-40850;
    else if(pdfset>=80260&&pdfset<=80286) return pdfset-80260;
    else if(pdfset>=192800&&pdfset<=192900) return pdfset-192800;
    else {
      std::cout << "pdfset=" << pdfset << " is not defined in pdfMember()" << std::endl;
      abort();
    }
  }

  int setPdf(std::string pdfname, int member){
    for(unsigned int i=0; i<m_pdfsetNames.size(); i++){
      if(m_pdfsetNames[i]==pdfname) {
// ~~~~~~~~~
        LHAPDF::usePDFMember(i+1, member);
// ~~~~~~~~~
        return i+1;
      }
    }

    // max pdfset (defined in lhapdf library)
    if(m_pdfsetNames.size()==(unsigned int)m_nmxset) {
      m_pdfsetNames.pop_back();
      m_pdfsetMembers.pop_back();
    }

    // add pdfset
    m_pdfsetNames.push_back(pdfname);
    m_pdfsetMembers.push_back(member);
    int nset = m_pdfsetNames.size();

// ~~~~~~~~~
    std::cout << "PDFTool::initPDFSet " << nset << " "
              << pdfname.c_str() << " " << member << std::endl;
    LHAPDF::initPDFSet(nset, pdfname.c_str());
    LHAPDF::usePDFMember(nset, member);
// ~~~~~~~~~

    return nset;
  }

  /////////////////////////////////////////// q2, x, f

  double getQ2() { return m_q2; }
  double getX1() { return m_x1; }
  double getX2() { return m_x2; }
  int getF1() { return m_f1; }
  int getF2() { return m_f2; }
  int getPdfset() { return m_pdfset; }
  void setPdfset(int init_pdfset) { m_pdfset=init_pdfset; }

  // set q2,x,f manually
  void setEventInfo(double q2, double x1, double x2, int f1, int f2) {//i.e. pow(mcevt_pdf_scale,2), mcevt_pdf_x1, mcevt_pdf_x2, mcevt_pdf_id1, mcevt_pdf_id2
    m_q2 = q2;
    m_x1 = x1;
    m_x2 = x2;
    m_f1 = f1;
    m_f2 = f2;
    m_newvalues = true;
  }


 private:
  bool m_newvalues;
  double m_q2;
  double m_x1;
  double m_x2;
  int m_f1;
  int m_f2;
  double m_pdf1org;
  double m_pdf2org;
  std::vector<std::string> m_pdfsetNames;
  std::vector<int> m_pdfsetMembers;

  double m_energy;

  int m_pdfset;
  int m_nmxset;

  int m_reweightPdf;
  double m_reweightScale;
};

#endif

