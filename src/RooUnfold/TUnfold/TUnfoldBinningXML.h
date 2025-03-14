// Author: Stefan Schmitt
// DESY, 10/08/11

//  Version 17.9, parallel to changes in TUnfold
//
//  History:
//    Version 17.8, relaxed DTD definition
//    Version 17.7, in parallel to changes in TUnfold
//    Version 17.6, with updated doxygen comments
//    Version 17.5, in parallel to changes in TUnfold
//    Version 17.4, in parallel to changes in TUnfoldBinning
//    Version 17.3, support for repeated bins with the same width
//    Version 17.2, XML interface for class TUnfoldBinning

#ifndef ROOT_TUnfoldBinningXML
#define ROOT_TUnfoldBinningXML


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                                                                      //
//  TUnfoldBinningXML, an auxillary class to read and write             //
//  complex binning schemes in XML                                      //
//                                                                      //
//  Citation: S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

/*
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "RooUnfold/TUnfold/TUnfoldBinning.h"
#include <iostream>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TXMLNode.h>
#include <TXMLDocument.h>
#include <iostream>

class TXMLNode;
class TXMLDocument;

#define TUnfoldBinningXML TUnfoldBinningXMLV17

class TUnfoldBinningXMLV17 : public TUnfoldBinningV17 {
 public:
   /********************** XML interface to read binning schemes *************/
static TUnfoldBinningXMLV17 *ImportXML(const TXMLDocument *document,const char *name); // import binning scheme
   static Int_t ExportXML(const TUnfoldBinningV17 &binning,std::ostream &out,Bool_t writeHeader,Bool_t writeFooter,Int_t indent=0); // append binning scheme to file
   Int_t ExportXML(const char *fileName) const; // export this binning scheme
   static void WriteDTD(const char *fileName="tunfoldbinning.dtd"); // write DTD file
   static void WriteDTD(std::ostream &out); // write DTD to stream

   /// construct a new binning scheme, for use with the root streamer
   TUnfoldBinningXMLV17 (const char *name=0,Int_t nBins=0,const char *binNames=0) 
      : TUnfoldBinningV17 (name,nBins,binNames) { }
 protected:
   static TUnfoldBinningXMLV17 *ImportXMLNode(TXMLNode *node); // import the given node as binning scheme
   void AddAxisXML(TXMLNode *node); // import axis information
protected:

   ClassDef(TUnfoldBinningXMLV17, TUnfold_CLASS_VERSION) //Complex binning schemes for TUnfoldDensity
};

#endif
