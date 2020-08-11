#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
CharacterVector getLocalSubstr(CharacterVector x, NumericVector y, double cutoff) {
  // Find local substrings with total size less than cutoff.
  // The result contains same-length strings start from each position.
  // Substrings which have occurred in the previous substring have been omitted with "".
  // e.g. A B C D got substrings ABC BC C D, returns "ABC", "", "", "D"
  // getLocalSubstr(LETTERS[1:4], c(1, 3, 5, 6), 10)
  int size = x.size();
  int l = 0;
  CharacterVector out(size);
  for (int i=0; i<size; i++) {
    // Init
    double s = 0.0;
    std::string sC = "";
    int jj = i;
    // Get sequence with total size less than cutoff
    for (int j=i; j<size; j++) {
      if (s + y[j] >= cutoff) break;
      s += y[j];
      sC += x[j];
      jj++;
    }

    // std::cout << "i:" << i << std::endl;
    // std::cout << "l:" << l << std::endl;
    // std::cout << "jj:" << jj << std::endl;

    if (jj <= l) {
      sC = "";
    }
    l = jj;

    out[i] = sC;
  }

  return out;
}
