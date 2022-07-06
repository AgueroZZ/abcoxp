/// @file fixed_smooth_frailty.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj


template<class Type>
Type fixed_smooth_frailty(objective_function<Type>* obj)
{
  // Read in data
  // DATA_VECTOR(t); // ORDERED survival times (Don't need this; input the ranks directly)
  DATA_IVECTOR(cens); // censoring indicators, with 0 denotes right-censoring
  DATA_IVECTOR(ranks); // rank of each observation, correcting for ties by Breslow method
  int n = ranks.size(); // Sample size
  DATA_SPARSE_MATRIX(B1); // Design matrix- random
  DATA_SPARSE_MATRIX(B2); // Design matrix- random
  DATA_SPARSE_MATRIX(P2); // Precision matrix- random
  DATA_SPARSE_MATRIX(X); // Design matrix- fixed
  DATA_SPARSE_MATRIX(D); // Difference matrix
  
  
  int K1 = B1.cols(); // Number of groups
  int K2 = B2.cols(); // Number of knots

  DATA_SCALAR(u1); // pc prior, u param
  DATA_SCALAR(alpha1); // pc prior, alpha param
  DATA_SCALAR(u2); // pc prior, u param
  DATA_SCALAR(alpha2); // pc prior, alpha param

  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)
  
  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta), eta = B1 * U1 + B2 * U2 + X * beta
  PARAMETER(theta1); // theta = -2log(sigma1)
  PARAMETER(theta2); // theta = -2log(sigma2)

  // Split the param into random effect coefficients and polynomial coefficients
  int Wdim = W.size();
  int betadim = Wdim - K1 - K2;
  vector<Type> U1(K1);
  vector<Type> U2(K2);
  vector<Type> beta(betadim);
  for (int i=0;i<K1;i++) U1(i) = W(i);
  for (int i=0;i<K2;i++) U2(i) = W(i + K1);
  for (int i=0;i<betadim;i++) beta(i) = W(i + K1 + K2);
  REPORT(U1);
  REPORT(U2);
  REPORT(beta);
  
  // Transformations
  vector<Type> eta = B1 * U1 + B2 * U2 + X * beta;
  REPORT(eta); // Check this works
  vector<Type> delta_red = D * eta;
  REPORT(delta_red);
  vector<Type> delta(n);
  delta(0) = 0;
  for (int i=1;i<n;i++) delta(i) = delta_red(i-1);
  REPORT(delta);
  Type sigma1 = exp(-0.5*theta1);
  REPORT(sigma1);
  Type sigma2 = exp(-0.5*theta2);
  REPORT(sigma2);
  
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
  Type U1U1 = (U1 * U1).sum();
  lpW += -0.5 * exp(theta1) * U1U1; // U part
  Type U2U2 = (U2 * U2).sum();
  lpW += -0.5 * exp(theta2) * U2U2; // U part
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  
  // Log determinant
  Type logdet1 = K1 * theta1;
  lpW += 0.5 * logdet1; // frailty part
  Type logdet2 = K2 * theta2;
  lpW += 0.5 * logdet2; // IWP part
  Type logdet3 = betadim * log(betaprec);
  lpW += 0.5 * logdet3; // beta part
  REPORT(lpW);
  
  // Log prior for theta
  Type lpT = 0;
  Type phi1 = -log(alpha1) / u1;
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta1) - 0.5*theta1;
  Type phi2 = -log(alpha2) / u2;
  lpT += log(0.5 * phi2) - phi2*exp(-0.5*theta2) - 0.5*theta2;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);
  
  return logpost;
}



#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this





