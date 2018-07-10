/// Un TMB modelo de tipo Cormack-Jolly-Seber de seccion 14.5 of
/// Korner-Nievergelt et al 2015. Actualizado para TMB y maximo
/// verosimiltud. Cole Monnahan, 7/2018.

#include <TMB.hpp>

// logit funcion
template<class Type>
Type inv_logit(Type x){
  Type y= 1/(1+exp(-x));
  return(y);
}

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
  DATA_INTEGER(I);              // numero de los individuos
  DATA_INTEGER(K);              // numero de los periodos
  DATA_INTEGER(nfam);		//  numero de las familias
  DATA_MATRIX(CH); 		// la historia de las capturas (0/1)
  DATA_VECTOR(carez);		// covariable
  DATA_IVECTOR(year);		// indice de los anos
  DATA_VECTOR(agec);		// covariable
  DATA_IVECTOR(family);		// indice de las familias
  DATA_IVECTOR(last);		// el ultimo periodo visto

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
  Type nll=0.0;			// negativa log verosimiltud
  matrix<Type> p(I,K); 		// probabilidad de las capturas
  matrix<Type> phi(I,K-1);	// prob. de la sobrevivencia
  matrix<Type> chi(I,K+1);	// prob. nunca de ver despues un periodo

  // para constrenir mas que 0
  Type sigmayearphi2=exp(sigmayearphi);
  Type sigmaphi2=exp(sigmaphi);
  Type sigmap2=exp(sigmap);

  p.setZero();
  phi.setZero();
  chi.setZero();

  int k;
  Type x;
  // TMB usa indices de 0, no de 1, entoces tenemos que estar cuidadoso, y
  // uso un "-1" para ser claro que pasa.
  for(int i=0; i<I; i++){ // loop over each individual
    // calcular prob. capturas como una funcion de efectos y covariables
    for(int t=0; t<(K-1); t++) {
      x=a(t)+ a1*carez(i)+
      	sigmayearphi2*yeareffphi_raw(year(i)-1)+
      	sigmaphi2*fameffphi_raw(family(i)-1);
      phi(i,t) = inv_logit(x);
    }
    // calcular prob. sobrevivencia como una funcion de efectos y
    // covariables
    p(i,1-1) = 1;  // periodo uno es la marca 
    for(int t=1; t<K; t++){
      x=b0(year(i)-1)+ b1(year(i)-1)*agec(t)+
      	sigmap2*fameffp_raw(family(i)-1);
      p(i,t) = inv_logit(x);
    }
    // la probabilidad de no ser visto nunca mas, usa un indice reverso
    // para calcular hacia atras usando recursion.
    chi(i,K+1-1) = 1.0;
    k = K;
    while (k > 1) {
      chi(i,k-1) = (1 - phi(i,k-1-1)) + phi(i,k-1-1) * (1 - p(i,k-1)) * chi(i,k+1-1);
      k = k - 1;
    }
    chi(i,1-1) = (1 - p(i,1-1)) * chi(i,2-1);
  }

  //// calcular la verosimiltud
  // efectos aleatorios
  nll-=dnorm(fameffphi_raw, Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(fameffp_raw,Type(0.0), Type(1.0), true).sum();
  nll-=dnorm(yeareffphi_raw, Type(0.0), Type(1.0), true).sum();
  // los datos
  for(int i=0; i<I; i++){ 
    // probabilidad de sobrevivencia, que es conocido porque k<last
    for (int t=1; t<last(i); t++) {
    	nll-= log(phi(i,t-1));
    }
    // probabilidad de captura, dado viva (como CH[i,t]~bernoulli(p[i,t]);)
    for(int t=0; t< last(i); t++){
      if(CH(i,t)==1){
    	nll-= log(p(i,t));
      } else {
    	nll-= log(1-p(i,t));
      }
    }
    // probabilidad de no ser capturado despues el proximo periodo fue visto 
    nll-= log(chi(i,last(i)+1-1));
  }

  // Informes sin incertidumbre
  REPORT(fameffphi_raw);
  REPORT(fameffp_raw);
  REPORT(yeareffphi_raw);
  REPORT(p);
  REPORT(chi);
  REPORT(phi);
  return(nll);
}
// final del archivo
