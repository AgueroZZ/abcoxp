/// @file fixed.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj


template<class Type>
Type fixed(objective_function<Type>* obj)
{
  // Read in data
  // DATA_VECTOR(t); // ORDERED survival times (Don't need this; input the ranks directly)
  DATA_IVECTOR(cens); // censoring indicators, with 0 denotes right-censoring
  DATA_IVECTOR(ranks); // rank of each observation, correcting for ties by Breslow method
  int n = ranks.size(); // Sample size
  DATA_SPARSE_MATRIX(X); // Design matrix- fixed
  DATA_SPARSE_MATRIX(D); // Differencing matrix to compute delta;
  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)

  // Parameter
  PARAMETER_VECTOR(W); //
  
  int betadim = X.cols();
  vector<Type> beta(betadim);
  for (int i=0;i<betadim;i++) beta(i) = W(i);
  REPORT(beta);


  // Transformations
  vector<Type> eta = X * beta;
  REPORT(eta); // Check this works
  vector<Type> delta_red = D * eta;
  REPORT(delta_red);
  vector<Type> delta(n);
  delta(0) = 0;
  for (int i=1;i<n;i++) delta(i) = delta_red(i-1);
  REPORT(delta);
  
  // Log likelihood
  Type ll = 0;
  for (int i=0;i<n;i++){
    int nn = n-ranks(i)+1;
    vector<Type> delta_vec_i(nn); //rank starts at 1!!!
    for(int j=0;j<nn;j++) {
      delta_vec_i(j) = delta(n - nn + j);
      }
    vector<Type> diffvec(nn);
    for (int j=0;j<nn;j++) {
      diffvec(j) = exp(delta(i) - delta_vec_i(j));
      }
    ll += -cens(i) * log(diffvec.sum());
  }
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;
  // Cross product
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  
  // Final result!
  Type logpost = -1 * (ll + lpW);
  
  return logpost;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
