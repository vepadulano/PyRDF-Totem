#ifndef __initialize_h_
#define __initialize_h_

#include "common_definitions.h"
#include "common_algorithms.h"
#include "parameters_global.h"
#include "parameters.h"
#include "common.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TH2D.h"

#include <unistd.h>
#include <cstdlib>

unsigned int timestamp_bins;
//#vector<Binning> binning_setup;
//vector<string> binnings;

void initialize() {
  printf("\n---------- RUNNING INITIALIZE ----------\n");
  Init("45b_56t");

  // default parameters
  auto input_n_si = 4.0;

  // select cuts
  anal.BuildCuts();
  anal.n_si = input_n_si;
  //
  // // alignment init
	// for (unsigned int i = 0; i < alignmentSources.size(); ++i)
	// {
	// 	printf("\n---------- alignment source %u ----------\n", i);
	// 	alignmentSources[i].Init();
  // }
  //
  timestamp_bins = timestamp_max - timestamp_min + 1;
  //
  // // binnings
	// binnings.push_back("ub");
	// //binnings.push_back("eb");
	// binnings.push_back("ob-1-10-0.2");
  // binnings.push_back("ob-1-30-0.2");
  //
  // for (unsigned int bi = 0; bi < binnings.size(); ++bi)
	// {
  //   unsigned int N_bins;
  //   double *bin_edges;
  //   BuildBinning(anal, binnings[bi], bin_edges, N_bins);
  //
  //   Binning b;
  //   b.N_bins = N_bins;
  //   b.bin_edges = bin_edges;
  //   binning_setup.push_back(b);
  // }
}

#endif
