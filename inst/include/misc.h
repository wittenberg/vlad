#ifndef MISC_H
#define MISC_H

inline arma::vec tapply3arma(arma::vec x, arma::ivec i) {
  std::map<int, std::vector<double> > groups;

  arma::vec::iterator x_it;
  arma::ivec::iterator i_it;

  for(x_it = x.begin(), i_it = i.begin(); x_it != x.end(); ++x_it, ++i_it) {
    groups[*i_it].push_back(*x_it);
  }
  arma::vec out(groups.size());

  std::map<int, std::vector<double> >::const_iterator g_it = groups.begin();
  arma::vec::iterator o_it = out.begin();
  for(; g_it != groups.end(); ++g_it, ++o_it) {
    *o_it = 0;
    for(std::vector<double>::const_iterator it = g_it->second.begin(); it != g_it->second.end(); ++it)
      *o_it += *it;
  }
  return out;
}

// armadillo vector has just one template type parameter
template< typename T, template <typename> class ARMA_VECTOR_TYPE >
ARMA_VECTOR_TYPE<T> vintersection( ARMA_VECTOR_TYPE<T> first, ARMA_VECTOR_TYPE<T> second )
{
  std::vector<T> output ;
  std::set_intersection( first.begin(), first.end(), second.begin(), second.end(),
                         std::back_inserter(output) ) ;
  std::reverse( output.begin(), output.end() ) ;

  ARMA_VECTOR_TYPE<T> result = arma::conv_to< ARMA_VECTOR_TYPE<T> >::from(output);
  return result ;
}
#endif
