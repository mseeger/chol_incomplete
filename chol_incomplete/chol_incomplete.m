%CHOL_INCOMPLETE Incomplete Cholesky decomposition for dense matrix
%  [LFACT,PIND] = CHOL_INCOMPLETE(N,ASEL,MAXD,PVTHRES,...)
%  Incomplete Cholesky decomposition for dense symmetric positive
%  definite matrix A (N-by-N).
%  Result is N-by-k lower triangular L, s.t. L*L' can be used as
%  low-rank approximation of P'*A*P, P a permutation matrix. The
%  algorithm is the usual Cholesky decomposition, but in each
%  iteration the next column is chosen to maximize the pivot (diagonal
%  element of L). L is returned in LFACT. L has at most MAXD
%  columns. The algorithm terminates with fewer columns if all
%  remaining pivots are below PVTHRES. P is returned in the index
%  vector PIND, meaning that P maps component i to PIND(i). PIND
%  converts pos. in the new ordering (of L*L') to pos. in the old
%  ordering (of A). The index I = [PIND(1),...,PIND(k)] contains
%  the selected active columns of A, so that L(1:k,:) is the
%  Cholesky factor of A(I,I). In fact:
%    P'*A(:,I) == L*L(1:k,:)'
%
%  The matrix A is defined via a function computing DIAG(A) and a
%  function computing a column of A. A is selected by the code ASEL
%  and an arbitrary number of further input arguments (the number
%  and types depend on the ASEL value). In the moment, the following
%  variants for A are implemented.
%
%  ASEL=0: Covariance matrix for RBF (Gaussian) covariance function
%    K(x,y) = C exp[ -(w/2) sum_{i=0}^{d-1} (x_i-y_i)^2 ] + v_b
%  Additional parameters to CHOL_INCOMPLETE:
%  - DATA:    Data matrix X (N-by-d)
%  - DVEC:    Vector diag(X*X') (N)
%  - HYPPARS: Hyperparameter vector [w; C; v_b], all positive
%
%  ASEL=1: Covariance matrix for squared-exponential covariance
%          function
%    K(x,y) = C \exp[ -(1/2) \sum_{i=0}^{d-1} w_i (x_i-y_i)^2 ] +
%    v_b
%  Additional parameters to CHOL_INCOMPLETE:
%  - DATA:    Data matrix X (N-by-d)
%  - DVEC:    Vector diag(X*W*X'), W=diag(w_i) (N)
%  - HYPPARS: Hyperparameter vector [w_1; ...; w_d; C; v_b], all
%             positive
