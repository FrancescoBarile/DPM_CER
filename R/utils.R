# Function to get the full adiacency matrix from its lower triangular representation
from_low_adj_to_full_adj<-function(low_adj){
  E<-length(low_adj)
  N<-(1+sqrt(1+8*E))/2
  m1 <- matrix(NA, N, N)
  m1<-`diag<-`(m1, 0)
  m1[lower.tri(m1, diag=FALSE)] <- low_adj
  m2 <- t(m1)
  m2[lower.tri(m2, diag=FALSE)] <- low_adj
  return(m2)
}
# Function to get the lower adiacency matrix from its full representation
from_full_adj_to_low_adj<-function(full_adj){
  low_adj<-t(full_adj)[lower.tri(t(full_adj))]
  return(low_adj)
}