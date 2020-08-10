#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int LCS(std::string x, std::string y) {
  vector<vector<int>> f(x.size() + 1, vector<int>(y.size() + 1, 0));
  int max = -1;
  for (int i = 1; i <= x.size(); i++) {
    for (int j = 1; j <= y.size(); j++) {
      if (x[i - 1] != y[j - 1]) f[i][j] = 0;
      else if (x[i - 1] == y[j - 1]) f[i][j] = f[i - 1][j - 1] + 1;
      if (max < f[i][j]) {
        max = f[i][j];
      }
    }
  }
  return max;
}


// [[Rcpp::export]]
IntegerMatrix LCSMatrix(StringVector x, StringVector y, bool match) {
  int nrow = x.size(), ncol = y.size();
  IntegerMatrix out(nrow, ncol);
  std::string xv;
  std::string yv;

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      xv = std::string(x[i]);
      yv = std::string(y[j]);
      if (match) {
        out(i, j) = LCS(xv, yv);
      } else {
        out(i, j) = std::max(xv.size(), yv.size()) - LCS(xv, yv);
      }
    }
  }

  return out;
}
