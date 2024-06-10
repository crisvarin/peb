#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_stats_rcpp(IntegerVector outcome,
                                 IntegerVector team1,
                                 IntegerVector team2,
                                 bool home_effect)
{
  int ngames = outcome.size();
  /* output */
  NumericMatrix stats(3, 3);
  
  for (int i = 0; i < (ngames - 1); i++) {
    for (int j = (i + 1); j < ngames; j++) {
      if ((team1[i] == team1[j]) || (team2[i] == team2[j]))
        stats(outcome[i] - 1, outcome[j] - 1) += 1;
      if (home_effect == false) {
        if ((team1[i] == team2[j]) || (team2[i] == team1[j]))
          stats(outcome[i] - 1, 3 - outcome[j]) += 1;
      }
    }
  }
  /* bye bye */
  return stats;
}




 
