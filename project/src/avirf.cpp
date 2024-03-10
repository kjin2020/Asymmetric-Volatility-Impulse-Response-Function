#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;  

int indicatorFunction(arma::mat r, arma::mat signs){
  r = r.t();
  
  int indicator = 1;
  int n = r.n_rows;
  for (int i = 0; i<n; i++){
    if(arma::as_scalar(signs.row(i)) * arma::as_scalar(r.row(i)) < 0){
      indicator = 0;
    }
  }
  return indicator;
}

// [[Rcpp::export]]
double expected_indicator_value(arma::mat r, arma::mat signs){
  int N =r.n_rows;
  double exp_indicator_value = 0;
  for (int i=0; i<N;i++){
    exp_indicator_value+=indicatorFunction(r.row(i),signs);
  }
  exp_indicator_value=exp_indicator_value/N;
  return exp_indicator_value;
}

// [[Rcpp::export]]
arma::mat elimination_mat(const int& n) {
  // Generates an elimination matrix for size 'n'
  int n1 = n * (n + 1) / 2;
  int n2 = pow(n, 2);
  
  arma::mat init = arma::eye(n1, n1);
  int oes = 1;
  
  arma::mat eli  = init.col(0);
  int block = n;
  
  while (eli.n_cols < n2) {
    if (eli.n_cols == 1) {
      eli = init.cols(0, block-1);
    } else {
      eli = arma::join_horiz(eli, init.cols(0, block-1));
    }
    
    if (init.n_cols > 1) {
      init = init.cols(block, init.n_cols-1);
    }
    
    eli = arma::join_horiz(eli, arma::zeros(eli.n_rows, oes));
    
    oes += 1;
    
    block -= 1;
  }
  
  return eli.cols(0, eli.n_cols-n-1);
}

// [[Rcpp::export]]
arma::mat commutation_mat(const int& n) {
  // generates a (square) commutation matrix for 'n'
  arma::mat K = arma::zeros(pow(n, 2), pow(n, 2));
  
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= n; j++) {
      K(i + n*(j - 1)-1, j + n*(i - 1)-1) = 1;
    }
  }
  return K;
}

// [[Rcpp::export]]
arma::mat duplication_mat(const int& n) {
  // Generates a duplication matrix for size 'n'
  int n2 = pow(n, 2);
  
  arma::mat el = elimination_mat(n);
  arma::mat co = commutation_mat(n);
  arma::mat m = arma::eye(n2, n2) + co;
  
  arma::mat dup = m*el.t()*arma::inv(el*m*el.t());
  
  return dup;
}

// [[Rcpp::export]]
arma::mat inv_gen(const arma::mat& m) {
  // Checks if a matrix is positive definit and calculates
  // the inverse or generalized inverse
  
  if (m.is_sympd() == TRUE) {
    return arma::inv_sympd(m);
  } else {
    return arma::pinv(m);
  }
}

// [[Rcpp::export]]
arma::mat eigen_value_decomposition(arma::mat& A){
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym( eigval, eigvec, A );
  
  
  arma::mat diag_mat_eigval = arma::diagmat(sqrt(eigval));
  return eigvec*diag_mat_eigval*eigvec.t();
  
}

// [[Rcpp::export]]

List updateH_A(const arma::mat &z, arma::mat &H, const arma::mat &A, 
               const arma::mat &tA,const arma::mat &B,const arma::mat &tB,
               const arma::mat &G,const arma::mat &tG,const arma::mat &CC, 
               int i,const arma::mat signs) {
  arma::mat e = z.row(i) * arma::chol(H) ;
  arma::mat Xi = indicatorFunction(e, signs) * e.t() * e;

  H = CC + tA * e.t() * e * A + indicatorFunction(e,signs)* tB * e.t() * e * B + tG * H * G;
  return Rcpp::List::create(Rcpp::Named("H") = Rcpp::wrap(H), 
                                   Rcpp::Named("Xi") = Rcpp::wrap(Xi),
                                   Rcpp::Named("e") = Rcpp::wrap(e));
}

// [[Rcpp::export]]
List monte_carlo_xi(const arma::mat H, const arma::mat z1, const arma::mat z2, const arma::mat A, const arma::mat tA, const arma::mat B, const arma::mat tB, const arma::mat G, const arma::mat tG,  const arma::mat CC, int v, int n, int N, int time, arma::mat& sign){
  arma::mat xi_1(N,N*v, arma::fill::zeros);
  arma::mat xi_2(N,N*v, arma::fill::zeros);
  for (int k = 0; k < n; ++k) {
    arma::mat ht1 = H;
    arma::mat ht2 = H;
    arma::uvec col_indices(N);
    for (size_t i = 0; i < col_indices.size(); ++i) {
      col_indices(i) = static_cast<unsigned int>(k + i * n);
    }
    arma::mat tmp1(N, N*v, arma::fill::zeros);
    arma::mat tmp2(N, N*v, arma::fill::zeros);
    for (int i = 0; i < v; ++i) {
      List lst_tmp1 = updateH_A(z1.cols(col_indices), ht1, A,tA,B,tB,G,tG,CC,i,sign);
      List lst_tmp2 = updateH_A(z2.cols(col_indices), ht2, A,tA,B,tB,G,tG,CC,i,sign);
      ht1  = as<arma::mat>(lst_tmp1["H"]);
      ht2  = as<arma::mat>(lst_tmp2["H"]);
      tmp1.submat(0, i*N, N-1, i*N+N-1) = as<arma::mat>(lst_tmp1["Xi"]);
      tmp2.submat(0, i*N, N-1, i*N+N-1) = as<arma::mat>(lst_tmp2["Xi"]);
    }
    xi_1 += tmp1;
    xi_2 += tmp2;
  }
  
  return Rcpp::List::create(Rcpp::Named("xi_1") = Rcpp::wrap(xi_1/n), 
                                     Rcpp::Named("xi_2") = Rcpp::wrap(xi_2/n));
}


// [[Rcpp::export]]
arma::mat virf_bekk_asymm(arma::mat H_t, arma::mat& A, arma::mat& B, arma::mat& G, arma::mat& CC, arma::mat& shocks, arma::mat& z1, arma::mat& z2, int& time, int& periods, int iterations, int& N, arma::mat& sign_1,arma::mat& sign_2,arma::mat& D_duplication, arma::mat& D_gen_inv){
  arma::mat H = H_t.submat(0,time*N-N,N-1,time*N-1);
  arma::mat virf = arma::zeros(periods, N*(N+1)/2);
  List xi_matrix = monte_carlo_xi(H,z1,z2,A,A.t(),B,B.t(),G,G.t(),CC,periods,iterations,N,time,sign_1);
  arma::mat Q_t = eigen_value_decomposition(H);
  arma::mat A_virf = D_gen_inv * kron(A,A).t() * D_duplication;
  arma::mat B_virf = D_gen_inv * kron(B,B).t() * D_duplication;
  arma::mat G_virf = D_gen_inv * kron(G,G).t() * D_duplication;
  arma::mat r = Q_t * shocks.row(0).t();
  arma::mat virf_temp = A_virf * D_gen_inv * arma::vectorise(r * r.t() - H)  + B_virf * D_gen_inv * arma::vectorise(indicatorFunction(r.t(),sign_2) * r * r.t() - as<arma::mat>(xi_matrix["xi_2"]).submat(0,0,N-1,N-1));
  virf.row(0) = virf_temp.t();
for (int i=1; i < periods; i++){
  virf_temp =  (A_virf+G_virf) * virf_temp + B_virf * D_gen_inv * arma::vectorise(as<arma::mat>(xi_matrix["xi_1"]).submat(0,i*N,N-1,(i+1)*N-1)-as<arma::mat>(xi_matrix["xi_2"]).submat(0,i*N,N-1,(i+1)*N-1));
  virf.row(i) = virf_temp.t();
}
  return virf;
}

// [[Rcpp::export]]

arma::mat virf_bekka(int start,int end, arma::mat& H_t, arma::mat& A, arma::mat& B, arma::mat& G, arma::mat& CC, arma::mat& shocks, arma::mat& z1, arma::mat& z2, int& periods, int iterations, int& N, arma::mat& sign_1,arma::mat& sign_2){
  arma::mat virf = arma::zeros(periods, N*(N+1)/2);
  arma::mat D_duplication = duplication_mat(N);
  arma::mat D_gen_inv = arma::inv(D_duplication.t() * D_duplication) * D_duplication.t();
  for (int time = start; time < end; ++time) {
    virf += virf_bekk_asymm(H_t, A, B, G, CC, shocks, z1, z2, time, periods, iterations, N, sign_1, sign_2,D_duplication,D_gen_inv);
    }
  return virf/(end-start);
}
