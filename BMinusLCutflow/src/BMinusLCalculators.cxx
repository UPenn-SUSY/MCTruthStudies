#include "BMinusLCutflow/include/BMinusLCalculators.h"
#include "TruthNtupleLooper/include/Calculators.h"
#include <math.h>
#include <stdlib.h>

#include "TruthNtupleLooper/include/ObjectDefs.h"

// =============================================================================
static const double PI = 3.14159265359;

// -----------------------------------------------------------------------------
std::vector< std::pair< TruthNtuple::Particle*
                      , TruthNtuple::Particle*
                      >
           > BMinusL::doMblDiffPairing( const std::vector<TruthNtuple::Particle*>& b_list
                                      , const std::vector<TruthNtuple::Particle*>& l_list
                                      )
{
  std::vector<std::pair<TruthNtuple::Particle*, TruthNtuple::Particle*> > best_pairs;

  // For now, this will handle lists of exactly 2 b quarks and 2 leptons
  // TODO make this function more general
  size_t num_b = b_list.size();
  size_t num_l = l_list.size();
  if (num_b == 2 && num_l == 2) {
    double mbl_00 = TruthNtuple::invariantMass(b_list.at(0), l_list.at(0));
    double mbl_01 = TruthNtuple::invariantMass(b_list.at(0), l_list.at(1));
    double mbl_10 = TruthNtuple::invariantMass(b_list.at(1), l_list.at(0));
    double mbl_11 = TruthNtuple::invariantMass(b_list.at(1), l_list.at(1));

    double mbl_diff_00_11 = fabs(mbl_00 - mbl_11);
    double mbl_diff_01_10 = fabs(mbl_01 - mbl_10);

    if (mbl_diff_00_11 <= mbl_diff_01_10) {
      best_pairs.push_back(std::make_pair( b_list.at(0), l_list.at(0) ));
      best_pairs.push_back(std::make_pair( b_list.at(1), l_list.at(1) ));
    }
    else {
      best_pairs.push_back(std::make_pair( b_list.at(0), l_list.at(1) ));
      best_pairs.push_back(std::make_pair( b_list.at(1), l_list.at(0) ));
    }
  }

  return best_pairs;
}

// -----------------------------------------------------------------------------
std::vector< std::pair< TruthNtuple::Particle*
                      , TruthNtuple::Particle*
                      >
           > BMinusL::doMblRatioPairing( const std::vector<TruthNtuple::Particle*>& b_list
                                       , const std::vector<TruthNtuple::Particle*>& l_list
                                       )
{
  std::vector<std::pair<TruthNtuple::Particle*, TruthNtuple::Particle*> > best_pairs;

  // For now, this will handle lists of exactly 2 b quarks and 2 leptons
  // TODO make this function more general
  size_t num_b = b_list.size();
  size_t num_l = l_list.size();
  if (num_b == 2 && num_l == 2) {
    double mbl_00 = TruthNtuple::invariantMass(b_list.at(0), l_list.at(0));
    double mbl_01 = TruthNtuple::invariantMass(b_list.at(0), l_list.at(1));
    double mbl_10 = TruthNtuple::invariantMass(b_list.at(1), l_list.at(0));
    double mbl_11 = TruthNtuple::invariantMass(b_list.at(1), l_list.at(1));

    // compute the ratio (min mbl / max mbl)
    double mbl_ratio_00_11 = ( mbl_00 >= mbl_11
                             ? mbl_11/mbl_00
                             : mbl_00/mbl_11
                             );
    double mbl_ratio_01_10 = ( mbl_01 >= mbl_10
                             ? mbl_10/mbl_01
                             : mbl_01/mbl_10
                             );

    if (mbl_ratio_00_11 >= mbl_ratio_01_10) {
      // std::cout << "picking mbl_ratio_00_11 ("
      //           << mbl_ratio_00_11
      //           << ") instead of mbl_ratio_01_10("
      //           << mbl_ratio_01_10
      //           << ")\n";
      best_pairs.push_back(std::make_pair( b_list.at(0), l_list.at(0) ));
      best_pairs.push_back(std::make_pair( b_list.at(1), l_list.at(1) ));
    }
    else {
      // std::cout << "picking mbl_ratio_01_10 ("
      //           << mbl_ratio_01_10
      //           << ") instead of mbl_ratio_00_11("
      //           << mbl_ratio_00_11
      //           << ")\n";
      best_pairs.push_back(std::make_pair( b_list.at(0), l_list.at(1) ));
      best_pairs.push_back(std::make_pair( b_list.at(1), l_list.at(0) ));
    }
  }

  return best_pairs;
}

// -----------------------------------------------------------------------------
std::vector< std::pair< TruthNtuple::Particle*
                      , TruthNtuple::Particle*
                      >
           > BMinusL::doMblSqSumPairing( const std::vector<TruthNtuple::Particle*>& b_list
                                       , const std::vector<TruthNtuple::Particle*>& l_list
                                       )
{
  std::vector<std::pair<TruthNtuple::Particle*, TruthNtuple::Particle*> > best_pairs;

  // For now, this will handle lists of exactly 2 b quarks and 2 leptons
  // TODO make this function more general
  size_t num_b = b_list.size();
  size_t num_l = l_list.size();
  if (num_b == 2 && num_l == 2) {
    double mbl_00 = TruthNtuple::invariantMass(b_list.at(0), l_list.at(0));
    double mbl_01 = TruthNtuple::invariantMass(b_list.at(0), l_list.at(1));
    double mbl_10 = TruthNtuple::invariantMass(b_list.at(1), l_list.at(0));
    double mbl_11 = TruthNtuple::invariantMass(b_list.at(1), l_list.at(1));

    // compute square sum of the two mbl values
    double mbl_sq_sum_00_11 = (mbl_00*mbl_00 - mbl_11*mbl_11);
    double mbl_sq_sum_01_10 = (mbl_01*mbl_01 - mbl_10*mbl_10);

    if (mbl_sq_sum_00_11 >= mbl_sq_sum_01_10) {
      best_pairs.push_back(std::make_pair( b_list.at(0), l_list.at(0) ));
      best_pairs.push_back(std::make_pair( b_list.at(1), l_list.at(1) ));
    }
    else {
      best_pairs.push_back(std::make_pair( b_list.at(0), l_list.at(1) ));
      best_pairs.push_back(std::make_pair( b_list.at(1), l_list.at(0) ));
    }
  }

  return best_pairs;
}
