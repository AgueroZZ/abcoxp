/// @file smooth.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj


template<class Type>
Type smooth(objective_function<Type>* obj)
{
  // Read in data
  // DATA_VECTOR(t); // ORDERED survival times (Don't need this; input the ranks directly)
  DATA_IVECTOR(cens); // censoring indicators, with 0 denotes right-censoring
  DATA_IVECTOR(ranks); // rank of each observation, correcting for ties by Breslow method
  int n = ranks.size(); // Sample size
  DATA_SPARSE_MATRIX(B2); // Design matrix- B
  DATA_SPARSE_MATRIX(P2); // Penalty matrix
  DATA_SPARSE_MATRIX(D); // Differencing matrix to compute delta;
  
  int d = P2.cols(); // Number of B-Spline coefficients
  DATA_SCALAR(logP2det); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u); // pc prior, u param
  DATA_SCALAR(alpha); // pc prior, alpha param

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U), eta = B2 * U
  PARAMETER(theta); // theta = -2log(sigma)
  
  // Split the param into B-Spline coefficients and polynomial coefficients
  int Wdim = W.size();
  vector<Type> U(d);
  for (int i=0;i<d;i++) U(i) = W(i);
  REPORT(U);

  // Transformations
  vector<Type> eta = B2 * W;
  REPORT(eta); // Check this works
  vector<Type> delta_red = D * eta;
  REPORT(delta_red);
  vector<Type> delta(n);
  delta(0) = 0;
  for (int i=1;i<n;i++) delta(i) = delta_red(i-1);
  REPORT(delta);
  Type sigma = exp(-0.5*theta);
  REPORT(sigma);
  
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
  vector<Type> P2U = P2*U;
  Type UP2U = (U * P2U).sum();
  lpW += -0.5 * exp(theta) * UP2U; // U part

  
  // Log determinant
  Type logdet1 = d * theta + logP2det;
  lpW += 0.5 * logdet1; // P part
  REPORT(logdet1);
  REPORT(lpW);
  
  // Log prior for theta
  Type lpT = 0;
  Type phi = -log(alpha) / u;
  lpT += log(0.5 * phi) - phi*exp(-0.5*theta) - 0.5*theta;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}




#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this





