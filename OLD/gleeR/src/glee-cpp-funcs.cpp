#include <Rcpp.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
NumericMatrix calcXbarStar(NumericMatrix A, double n_A, double n_B, double n_iter, double xbar0_replace) {
	int nA = (int) n_A;
	int nB = (int) n_B;
	int iter = (int) n_iter;
	
	boost::random::mt19937 gen;
	// boost::random::mt19937 gen(std::time(0));
	boost::random::uniform_int_distribution<> distA(0, nA-1);
	
	NumericMatrix xbar_star(iter*A.nrow(),2);
	for (int i=0; i<A.nrow(); i++) {
		for (int k=0;  k<iter; k++) {
			int idx = (i*iter) + k;
			xbar_star(idx,0) = 0;
			xbar_star(idx,1) = 0;
			// A
			for (int j=0; j<nA; ++j) { xbar_star(idx,0) += A(i, distA(gen)); }
			if (xbar_star(idx,0)==0) { xbar_star(idx,0) = xbar0_replace; } else { xbar_star(idx,0) /= nA; }
			// B
			for (int j=0; j<nB; ++j) { xbar_star(idx,1) += A(i, distA(gen)); }
			if (xbar_star(idx,1)==0) { xbar_star(idx,1) = xbar0_replace; } else { xbar_star(idx,1) /= nB; }
		}
	}
	return xbar_star;
}

/*
calculate p-values of numbers in x, using distribution in f
f need not be sorted
*/
// [[Rcpp::export]]
NumericVector calcPvals(NumericVector x, NumericVector f) {
	NumericVector pVal(x.size());
	for (int i = 0; i < x.size(); i++) {
		double p = (double) std::count_if(f.begin(), f.end(), [thresh = x[i]](double elem) { return (elem > thresh); }) / f.size();
		pVal[i] = 2*std::min(p,1-p);
	}
	return pVal;
}
