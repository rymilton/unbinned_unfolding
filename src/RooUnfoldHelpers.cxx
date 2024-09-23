#include "RooUnfoldHelpers.h"
#include "TRandom.h"

#include <iostream>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

using std::cerr;

namespace RooUnfolding {
  void printVector(const char* name, const TVectorD& vec){
    std::cout << name << std::endl;
    for(int i=0; i<vec.GetNrows(); ++i){
      std::cout << vec[i] << std::endl;
    }
  }

  void add(TMatrixD& target, const TMatrixD& addition){
    for(int i=0; i<addition.GetNrows(); ++i){
      for(int j=0; j<addition.GetNcols(); ++j){
        target(i,j) = target(i,j) + addition(i,j);
      }
    }
  }

  TVectorD* resizeVector (const TVectorD& vec, Int_t n)
  {

    Int_t nbins;

    if (n < vec.GetNrows()){
      nbins = n;
    } else {
      cerr << "Warning: Requested vector size is smaller than the original vector. The initial size is maintained.";
      nbins = vec.GetNrows();
    }

    TVectorD* newvec = new TVectorD(nbins);

    // Add zeros to the new columns.
    for (Int_t i = 0; i < nbins; i++){
      if (i < vec.GetNrows()){
	(*newvec)[i] = vec[i];
      } else {
	(*newvec)[i] = 0.0;
      }
    }
  
    return newvec;
  }

  void resizeVector (TVectorD& vec, Int_t n)
  {

    if (n == vec.GetNrows()){
      return;
    }

    if (n < vec.GetNrows()){
      cerr << "Warning: Requested vector size is smaller than the original vector. The initial vector is maintained.";
      return;
    }

    // Add new columns with zero as value.
    vec.ResizeTo(n);
  }

  TMatrixD* squareMatrix (const TMatrixD& matrix, int nbins)
  {

    if(nbins < 0) nbins = std::max(matrix.GetNrows(),matrix.GetNcols());

    TMatrixD* newmatrix = new TMatrixD(matrix);

    newmatrix->ResizeTo(nbins,nbins);
  
    return newmatrix;
  }

  void squareMatrix (TMatrixD& matrix, int nbins)
  {

    if(nbins < 0) nbins = std::max(matrix.GetNrows(),matrix.GetNcols());    
    
    matrix.ResizeTo(nbins,nbins);
  }


  TH1D* getTH1(const TVectorD& vec, const TVectorD& errvec, const char* name, const char* title, bool overflow){
    
    Int_t i_start = 1;

    if(overflow){
      i_start = 0;
    }

    TH1D* hist = new TH1D(name, title, vec.GetNrows(), 0, 1);

    for (int i = 0; i < vec.GetNrows(); i++){
      hist->SetBinContent(i + i_start, vec[i]);
      hist->SetBinError(i + i_start, errvec[i]);
    }
    
    return hist;
  }

  // Add an empty bin on both sides of the vector.
  void addEmptyBins(TVectorD& v){
    
    Int_t bins = v.GetNrows();
    bins+=2;
    
    v.ResizeTo(bins);

    for (Int_t i = bins; i > 1; i--){
      v[i - 1] = v[i - 2];
    }
    v[0] = 0.0;
  }

  // Add one layer of empty bins on all sides of the matrix.
  void addEmptyBins(TMatrixD& m){
    
    Int_t bins_x = m.GetNrows();
    Int_t bins_y = m.GetNcols();
    bins_x+=2;
    bins_y+=2;

    m.ResizeTo(bins_x, bins_y);

    for (Int_t i = bins_x; i > 0; i--){
      for (Int_t j = bins_y; j > 0; j--){
	
	if (i == 1 || j == 1){
	  m[i - 1][j - 1] = 0.0;
	} else {
	  m[i - 1][j - 1] = m[i - 2][j - 2];
	}
      }
    }
  }


  TH2D* getTH2(const TMatrixD& matrix, const char* name, const char* title, bool overflow){
    Int_t i_start = 1;

    if(overflow){
      i_start = 0;
    }

    TH2D* hist = new TH2D(name, title, matrix.GetNrows(), 0, 1, matrix.GetNcols(), 0, 1);

    for (int i = 0; i < matrix.GetNrows(); i++){
      for (int j = 0; j < matrix.GetNcols(); j++){

	hist->SetBinContent(i + i_start, j + i_start, matrix[i][j]);
      }
    }
    
    return hist;
  }

  void printMatrix(const TMatrixD& m, const char* name, const char* format, Int_t cols_per_sheet){
  // Print the matrix as a table of elements.
  // Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
  // By default the format "%11.4g" is used to print one element.
  // One can specify an alternative format with eg
  //  format ="%6.2f  "

  if (!m.IsValid()) {
    m.Error("PrintMatrix","%s is invalid",name);
    return;
  }

  const Int_t ncols  = m.GetNcols();
  const Int_t nrows  = m.GetNrows();
  const Int_t collwb = m.GetColLwb();
  const Int_t rowlwb = m.GetRowLwb();

  if (!(format && format[0])) format= "%11.4g ";
  char topbar[1000];
  snprintf(topbar,1000,format,123.456789);
  Int_t nch = strlen(topbar)+1;
  if (nch > 18) nch = 18;
  char ftopbar[20];
  for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
  Int_t nk = 1 + Int_t(log10(ncols));
  snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
  Int_t nch2 = strlen(ftopbar);
  for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
  ftopbar[nch] = '|';
  ftopbar[nch+1] = 0;

  printf("\n%dx%d %s is as follows",nrows,ncols,name);

  if (cols_per_sheet <= 0) {
    cols_per_sheet = 5;
    if (nch <= 8) cols_per_sheet =10;
  }
  nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
  for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
  topbar[nk] = 0;
  for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
    printf("\n\n     |");
    for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
      printf(ftopbar,j+collwb-1);
    printf("\n%s\n",topbar);
    if (m.GetNoElements() <= 0) continue;
    for (Int_t i = 1; i <= nrows; i++) {
      printf("%4d |",i+rowlwb-1);
      for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
        printf(format,m(i+rowlwb-1,j+collwb-1));
      printf("\n");
    }
  }
  printf("\n");
  }


  TMatrixD& ABAT (const TMatrixD& a, const TMatrixD& b, TMatrixD& c){
    // Fills C such that C = A * B * A^T. Note that C cannot be the same object as A.
    TMatrixD d (b, TMatrixD::kMultTranspose, a);
    c.Mult (a, d);
    return c;
  }

  TMatrixD& ABAT (const TMatrixD& a, const TVectorD& b, TMatrixD& c){
    // Fills C such that C = A * B * A^T, where B is a diagonal matrix specified by the vector.
    // Note that C cannot be the same object as A.
    TMatrixD d (TMatrixD::kTransposed, a);
    d.NormByColumn (b, "M");
    c.Mult (a, d);
    return c;
  }

  void printTable (std::ostream& o, int dim, int ntxb, int ntyb,
                   const TVectorD& vTrainTrue, const TVectorD& vTrain, const TVectorD& vTrue,const TVectorD& vMeas, const TVectorD& vReco,
                   ErrorTreatment withError, const TVectorD& eTrue, const TVectorD& eReco, double chi_squ, bool overflow){
    // Prints entries from truth, measured, and reconstructed data for each bin.
    int _nm= vTrain    .GetNrows();
    int _nt= vTrainTrue.GetNrows();
    std::ostringstream fmt;
    fmt.copyfmt (o);
    Int_t iwid= (dim==3) ? 8 : (dim==2) ? 7 : 5;
    const char* xwid= (dim==3) ? "===" : (dim==2) ? "==" : "";
    o << "===============================================================================" << xwid << std::endl;
    o << std::setw(iwid) << ""      << std::setw(9) << "Train" << std::setw(9) << "Train"    << std::setw(9) << "Test"  << std::setw(9) << "Test"  << std::setw(9) << "Unfolded";
    if (withError)
      o << std::setw(10)<<"Error on"<<std::setw(9) << "Diff" << std::setw(9) << "Pull" << std::endl;
    else
      o << std::setw(9) << "Diff" << std::endl;
    o << std::setw(iwid) << "Bin"   << std::setw(9) << "Truth" << std::setw(9) << "Measured" << std::setw(9) << "Truth" << std::setw(9) << "Input" << std::setw(9) << "Output";
    if (withError)
      o << std::setw(10)<<"Unfolding";
    o << std::endl;
    o << "===============================================================================" << xwid << std::endl;
    Double_t true_train_tot=0;
    Double_t meas_train_tot=0;
    Double_t true_test_tot=0;
    Double_t meas_test_tot=0;
    Double_t unf_tot=0;
    Double_t chi2= 0.0;
    Int_t ndf= 0, first= (overflow ? 0 : 1);
    Int_t maxbin= _nt < _nm ? _nm : _nt;
    for (Int_t i = 0 ; i < maxbin; i++) {
      Int_t it= i%_nt;
      Int_t im= i%_nm;
      if (dim==2 || dim==3) {
        Int_t iw= (dim==2) ? 3 : 2;
        Int_t ix= it%ntxb;
        Int_t iy= ((it-ix)/ntxb)%ntyb;
        o << std::setw(iw) << ix+first << ',' << std::setw(iw) << iy + first;
        if (dim==3) o << ',' << std::setw(iw) << ((it-ix)/ntxb - iy)/ntyb + first;
      } else {
        o << std::setw(iwid) << i+first;
      }
      o << std::fixed << std::setprecision(0);
      true_train_tot+=vTrainTrue(it);
      meas_train_tot+=vTrain(im);
      if (vTrue.GetNrows()>0) true_test_tot+=vTrue(it);
      meas_test_tot+=vMeas(im);
      unf_tot+=vReco(it);
      if (i<_nt){
        o << ' ' << std::setw(8) << vTrainTrue(it);
      }
      else
        o << std::setw(9) << ' ';
      if (i<_nm)
        o << ' ' << std::setw(8) << vTrain(im);
      else
        o << std::setw(9) << ' ';
      if (vTrue.GetNrows()>0 && i<_nt)
        o << ' ' << std::setw(8) << vTrue(it);
      else
        o << std::setw(9) << ' ';
      if (i<_nm)
        o << ' ' << std::setw(8) << vMeas(im);
      else
        o << std::setw(9) << ' ';
      o << std::setprecision(1);
      if (i<_nt) {
        Double_t y= vReco(it), yerr = eReco(it);
        o << ' ' << std::setw(8) << y;
        if (withError)
          o << ' ' << std::setw(9) << yerr;
        if (vTrue.GetNrows()>0 &&
            (                       y!=0.0 || (withError &&                   yerr>0.0)) &&
            (vTrue(it)!=0.0 || (withError && eTrue(it)>0.0))) {
          Double_t ydiff= y - vTrue(it);
          o << ' ' << std::setw(8) << ydiff;
          if (withError && yerr>0.0) {
            ndf++;
            Double_t ypull = ydiff/yerr;
            chi2 += ypull*ypull;
            o << ' ' << std::setw(8) << ypull;
          }
        }
      }
      o << std::endl;
    }

    o << "===============================================================================" << xwid << std::endl
      << std::setw(iwid) << "" << std::fixed << std::setprecision(0)
      << ' ' << std::setw(8) << true_train_tot
      << ' ' << std::setw(8) << meas_train_tot;
    if (vTrue.GetNrows()>0)
      o << ' ' << std::setw(8) << true_test_tot;
    else
      o << std::setw(9) << ' ';
    Double_t toterr= 0.0;
    if (meas_test_tot>0.0 && meas_train_tot>0.0) toterr= sqrt(meas_test_tot)*true_train_tot/meas_train_tot;
    o << ' ' << std::setw(8) << meas_test_tot << std::setprecision(1)
      << ' ' << std::setw(8) << unf_tot;
    if (withError) 
      o << ' ' << std::setw(9) << toterr;
    o << ' ' << std::setw(8) << unf_tot-true_test_tot;
    if(withError && toterr>0.0)
      o << ' ' << std::setw(8) <<(unf_tot-true_test_tot)/toterr;
    o << std::endl
      << "===============================================================================" << xwid << std::endl;
    o.copyfmt (fmt);
    if (vTrue.GetNrows()>0) {
      if (chi_squ!=-999.0) {
        o << "Chi^2/NDF=" << chi_squ << "/" << ndf << " (bin-by-bin Chi^2=" << chi2 << ")" << std::endl;
      } else {
        o << "Bin-by-bin Chi^2/NDF=" << chi2 << "/" << ndf << std::endl;
        chi_squ= chi2;
      }
      if (chi_squ<=0.0) std::cerr << "Warning: Invalid Chi^2 Value" << std::endl;
    }
  }

  void printTable (std::ostream& o, const TVectorD& vTrainTrue, const TVectorD& vTrain,
                   const TVectorD& vMeas, const TVectorD& vReco){
    TVectorD vTrue(vTrainTrue.GetNrows());
    TVectorD eTrue(vTrainTrue.GetNrows());
    TVectorD eReco(vReco.GetNrows());        
    printTable (o, 1, vTrainTrue.GetNrows(), 0, vTrainTrue, vTrain, vTrue, vMeas, vReco, kNoError, eTrue, eReco, -999.0,false);
  }

  void randomize(TVectorD& values, const TVectorD& errors, TRandom& rnd){
    for(int i=0; i<values.GetNrows(); ++i){
      values[i] = rnd.Gaus(values[i],errors[i]);
    }
  }
  
  void randomize(TVectorD& values, TRandom& rnd){
    for(int i=0; i<values.GetNrows(); ++i){
      values[i] = rnd.Poisson(values[i]);
    }
  }
  
  void randomize(TMatrixD& values, const TMatrixD& errors, TRandom& rnd){
    for(int i=0; i<values.GetNrows(); ++i){
      for(int j=0; j<values.GetNcols(); ++j){
        values(i,j) = rnd.Gaus(values(i,j),errors(i,j));
      }
    }
  }

  bool sanitize(TMatrixD& mat,double t){
    bool retval = false;
    if(mat.GetNcols() == mat.GetNrows()){
      const size_t n = mat.GetNcols();
      for(size_t i=0; i<n; ++i){
        if(fabs(mat(i,i)) < t){
          bool ok = false;
          for(size_t j=0; j<n; ++j){
            if(i==j) continue;
            if(fabs(mat(i,j)) > t) ok = true;
            if(fabs(mat(j,i)) > t) ok = true;
          }
          if(!ok){
            retval = true;
            mat(i,i) = 1.;
          }
        }
      }
    } else {
      // no implementation for non-square matrices yet, TODO!
    }
    return retval;
  }
  
  void mNorm (TMatrixD& m, const TVectorD& norm){
    // normalize Matrix to values of vector
    for (Int_t j= 0; j < m.GetNcols(); ++j) {
      double fac = norm[j];
      if (fac != 0.0) fac= 1.0/fac;
      for (Int_t i= 0; i < m.GetNrows(); ++i) {
        m(i,j)= m(i,j) * fac;
      }
    }
  }
  
}
