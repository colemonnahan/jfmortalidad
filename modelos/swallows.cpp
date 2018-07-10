/// See file swallows.stan for documentation 
template<class Type>
Type objective_function<Type>::operator() ()
{
  // str(data)  ## the real data
  // List of 11
  //  $ I     : int 322
  //  $ K     : int 18
  //  $ nfam  : num 72
  //  $ CH    : int [1:322, 1:18] 1 1 1 1 1 1 1 1 1 1 ...
  //   ..- attr(*, "dimnames")=List of 2
  //   .. ..$ : NULL
  //   .. ..$ : chr [1:18] "day0" "day1" "day2" "day3" ...
  //  $ carez : num [1:322] -0.584 -0.584 -0.584 -0.584 -0.584 ...
  //  $ year  : num [1:322] 1 1 1 1 1 1 1 1 1 1 ...
  //  $ agec  : num [1:18] -8.5 -7.5 -6.5 -5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5 ...
  //  $ family: num [1:322] 5 5 5 5 5 1 1 1 4 4 ...
  //  $ last  : int [1:322] 10 13 10 15 10 4 3 18 10 12 ...
  DATA_INTEGER(I);
  DATA_INTEGER(K);
  DATA_INTEGER(nfam);
  DATA_MATRIX(CH);
  DATA_VECTOR(carez);
  DATA_IVECTOR(year);
  DATA_VECTOR(agec);
  DATA_IVECTOR(family);
  DATA_IVECTOR(last);

  // efectos fijos
  PARAMETER(sigmayearphi);
  PARAMETER(sigmaphi);
  PARAMETER(sigmap);
  PARAMETER_VECTOR(a);
  PARAMETER(a1);
  PARAMETER_VECTOR(b0);
  PARAMETER_VECTOR(b1);
  // non-centered efectos aleatorios
  PARAMETER_VECTOR(fameffphi_raw);
  PARAMETER_VECTOR(fameffp_raw);
  PARAMETER_VECTOR(yeareffphi_raw);
  Type nll=0.0; // negativa log verosimiltud
  matrix<Type> p(I,K);
  matrix<Type> phi(I,K-1);
  matrix<Type> chi(I,K+1);

  // To bound below by zero
  Type sigmayearphi2=exp(sigmayearphi);
  Type sigmaphi2=exp(sigmaphi);
  Type sigmap2=exp(sigmap);

  p.setZero();
  phi.setZero();
  chi.setZero();

  int k;
  Type x;
  // TMB indexes from 0 not 1, so need to be careful to adjust that
  // below. I've added (-1) where needed.
  for(int i=0; i<I; i++){ // loop over each individual
    // calculate phi as a function of fixed and random effects
    for(int t=0; t<(K-1); t++) {
      x=a(t)+ a1*carez(i)+
      	sigmayearphi2*yeareffphi_raw(year(i)-1)+
      	sigmaphi2*fameffphi_raw(family(i)-1);
      phi(i,t) = inv_logit(x);
    }
    // calculate p as a function of fixed and random effects
    p(i,1-1) = 1;  // first occasion is marking occasion
    for(int t=1; t<K; t++){
      x=b0(year(i)-1)+ b1(year(i)-1)*agec(t)+
      	sigmap2*fameffp_raw(family(i)-1);
      p(i,t) = inv_logit(x);
    }
    // probabilitiy of never being seen after last observation. ind here is
    // a reverse index so this loop goes from K:2, recursively calculating
    // backward.
    chi(i,K+1-1) = 1.0;
    k = K;
    while (k > 1) {
      chi(i,k-1) = (1 - phi(i,k-1-1)) + phi(i,k-1-1) * (1 - p(i,k-1)) * chi(i,k+1-1);
      k = k - 1;
    }
    chi(i,1-1) = (1 - p(i,1-1)) * chi(i,2-1);
  }

  // random effects; non-centered
  nll-=dnorm(fameffphi_raw, Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(fameffp_raw,Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(yeareffphi_raw, Type(0.0), Type(1.0), true).sum();

  // // likelihood
  for(int i=0; i<I; i++){ // loop over each individual
    // probability of survival, known alive since k<last
    for (int t=1; t<last(i); t++) {
    	nll-= log(phi(i,t-1));
    }
    // // probability of observation given known alive
    for(int t=0; t< last(i); t++){
      // CH[i,t]~bernoulli(p[i,t]);
      if(CH(i,t)==1){
    	nll-= log(p(i,t));
      } else {
    	nll-= log(1-p(i,t));
      }
    }
    // probability of no observations after time period last
     nll-= log(chi(i,last(i)+1-1));
  }
  REPORT(fameffphi_raw);
  REPORT(fameffp_raw);
  REPORT(yeareffphi_raw);
  REPORT(p);
  REPORT(chi);
  REPORT(phi);
  return(nll);
}
