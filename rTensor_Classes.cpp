//
//  tensor.cpp
//  
//
//  Created by James Li on 8/16/13.
//
//

#include <Rcpp.h>
//#include <vector>
//#include <algorithm>


//[[Rcpp::export]]
Rcpp::CharacterVector ndTensor_validate(Rcpp::NumericVector x){
    //<need to implement checking>
    return Rcpp::wrap("Pass");
}












/*******************************************************************************
#include <marray.hxx>
namespace marray = andres;
 //[[Rcpp::export]]
 Rcpp::List marrayC (Rcpp::NumericVector x, Rcpp::IntegerVector modes_in){
 
 // Determine number of dimensions
 size_t num_modes = modes_in.size();
 
 //Copy IntegerVector modes_in into a vector
 //modes = Rcpp::as<std::vector<size_t> >(modes_in);
 
 //Initialize marray with dimensions from
 //std::vector<size_t>::iterator itr = modes.begin();
 Rcpp::IntegerVector::iterator itr = modes_in.begin();
 marray::Marray<double> a(marray::SkipInitialization, itr, itr+num_modes, marray::FirstMajorOrder);
 
 //Fill in marray with passed-in numeric vector
 std::copy(x.begin(), x.end(), a.begin());
 
 std::vector<size_t> modes (num_modes);
 for (size_t j=0; j != num_modes; ++j){
 modes[j] = a.shape(j);
 }
 
 return Rcpp::List::create(Rcpp::_["data"] = Rcpp::wrap(a), Rcpp::_["modes"] = Rcpp::wrap(modes), Rcpp::_["num_modes"] = a.dimension());
 }
*******************************************************************************/