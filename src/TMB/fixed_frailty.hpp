// @file fixed_frailty.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj


template<class Type>
Type fixed_frailty(objective_function<Type>* obj)
{
  // Read in data
  // DATA_VECTOR(t); // ORDERED survival times (Don't need this; input the ranks directly)
  DATA_IVECTOR(cens); // censoring indicators, with 0 denotes right-censoring
  DATA_IVECTOR(ranks); // rank of each observation, correcting for ties by Breslow method
  int n = ranks.size(); // Sample size
  DATA_SPARSE_MATRIX(B1); // Design matrix- random intercepts
  DATA_SPARSE_MATRIX(X); // Design matrix- fixed
  DATA_SPARSE_MATRIX(D); // Difference matrix
  
  
  int K = B1.cols(); // Number of groups
  DATA_SCALAR(u); // pc prior, u param
  DATA_SCALAR(alpha); // pc prior, alpha param
  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)
  
  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta), eta = B1 * U + X * beta
  PARAMETER(theta); // theta = -2log(sigma)
  
  // Split the param into random effect coefficients and polynomial coefficients
  int Wdim = W.size();
  int betadim = Wdim - K;
  vector<Type> U(K);
  vector<Type> beta(betadim);
  for (int i=0;i<K;i++) U(i) = W(i);
  for (int i=0;i<betadim;i++) beta(i) = W(i+K);
  REPORT(U);
  REPORT(beta);
  
  // Transformations
  vector<Type> eta = B1 * U + X * beta;
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
  Type UU = (U * U).sum();
  lpW += -0.5 * exp(theta) * UU; // U part
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  
  // Log determinant
  Type logdet1 = K * theta;
  lpW += 0.5 * logdet1; // frailty part
  Type logdet2 = betadim * log(betaprec);
  lpW += 0.5 * logdet2; // beta part
  REPORT(logdet1);
  REPORT(logdet2);
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

