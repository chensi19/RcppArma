// [[Rcpp::depends(RcppArmadillo)]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat grad_desc(arma::vec Y, arma::mat X, double step, int dim, int itr, int n1, int n2, arma::vec& Qtrace){
  arma::mat beta_trace(dim,itr);
  arma::vec beta(dim);
  arma::mat a1(1,1);
  arma::vec y_x_beta(dim);
  for (int i=0; i<itr;i++){
    y_x_beta = Y-X*beta;
    a1 = ((y_x_beta.st())*y_x_beta);
    Qtrace(i)=a1(0,0);
    beta=beta-(step*(-2*X.st()*y_x_beta));
    beta_trace.col(i)=beta;
  }
  return beta_trace.cols(n1,n2);
}

/*** R
Y=c(1,1,1,1)
X=matrix(c(1,1,1,1,2,2,2,2),nrow=4)
step<-c(0.001,0.005,0.01)
dim<-2
itr<-100
n1=1
n2=10
split.screen(c(2, 2))
for (k in 1:3){
  Qtrace<-rep(0,itr)
  print(grad_desc(Y,X,step[k],dim,itr,n1,n2,Qtrace))
  screen(k)
  plot(c(1:itr),Qtrace,xlab='iterations',ylab='loss function',main=paste("step=",+step[k]),col=2,pch=5)
}
*/
