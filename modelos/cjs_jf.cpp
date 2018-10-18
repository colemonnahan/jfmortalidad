/// Un TMB modelo de tipo Cormack-Jolly-Seber de seccion 14.5 of
/// Korner-Nievergelt et al 2015. Actualizado para TMB y maximo
/// verosimiltud. Cole Monnahan, 7/2018.

#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// logit funcion
template<class Type>
Type inv_logit(Type x){
  Type y= 1/(1+exp(-x));
  return(y);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(I);              // numero de los individuos
  DATA_INTEGER(K);              // numero de los periodos
  DATA_MATRIX(CH);		// la historia de las capturas (0/1) [IxK]
  DATA_IVECTOR(last);		// el ultimo periodo visto
  DATA_MATRIX(counts);		// numeros de recapturas (individos x periodos)
  DATA_VECTOR(effort);		// la esfuerza por cada periodo
  DATA_VECTOR(lengths);		// longitudes de los individuos
  DATA_IVECTOR(first);		// el primero periodo
  DATA_FACTOR(sexo);		// sexo de los individuos
  DATA_FACTOR(evento)		// evento de las marcas
  DATA_VECTOR(lengths_pred);	// para calcular selectividad
  DATA_VECTOR(esfuerzo_pred); // para calcular catchability
  
  // efectos fijos
  PARAMETER(logM);		// mortalidad de naturleza
  PARAMETER(logr);		// effecto de los numeros de capturas
  //  PARAMETER(logk);		// effecto de esfuerza de probabilidad de captura
  PARAMETER(a); 		// efecto de la selectividad
  PARAMETER(b); 		// efecto de la selectividad
  PARAMETER_ARRAY(tau);		// efectos aleatorios por k por machos
  PARAMETER_MATRIX(mu_tau);	// los promedios de los efectos aleatorios
  PARAMETER(logsigma_tau);	// la varianza de todos los efectos aleatorios
  Type nll=0.0;			// negativa log verosimiltud
  matrix<Type> p(I,K); 		// probabilidad de las capturas
  matrix<Type> phi(I,K);	// prob. de la sobrevivencia
  matrix<Type> chi(I,K+1);	// prob. nunca de ver despues un periodo

  p.setZero();
  phi.setZero();
  chi.setZero();

  int kk;
  Type M=exp(logM);
  Type r=exp(logr);
  Type sigma_tau=exp(logsigma_tau);
  // la primero columna es por machos y la segunda por hembras
  array<Type> k(K,2,3); // periodos por generos por eventos
  k.setZero();
  for(int t=0;t<K;t++) { // periodos
    for(int ev=0; ev<3; ev++){ // eventos
    // non-centered efectos aleatorios
      k(t,0,ev)=exp(tau(t,0,ev)*sigma_tau+mu_tau(0,ev)); // machos
      k(t,1,ev)=exp(tau(t,1,ev)*sigma_tau+mu_tau(1,ev)); // hembras
    }
  }

  // TMB usa indices de 0, no de 1, entoces tenemos que estar cuidadoso, y
  // uso un "-1" para ser claro que pasa.
  for(int i=0; i<I; i++){ // iterando sobre cada individuos
    // inicializacion estan vivo en periodo uno
    phi(i,first(i)-1)=exp(-M-counts(i,first(i)-1)*r);
    p(i,first(i)-1)=1;
    for(int t=first(i); t<K; t++) {
      // calcular prob. capturas como una funcion de efectos y covariables
      p(i,t) = (1-exp(-k(t,sexo(i),evento(i))*effort(t)))/(1+exp(a+b*lengths(i)));
      // calcular prob. sobrevivencia como una funcion de efectos y
      // covariables
      phi(i,t) = exp(-M-counts(i,t-1)*r);
    }
    // la probabilidad de no ser visto nunca mas, usa un indice reverso
    // para calcular hacia atras usando recursion.
    chi(i,K+1-1) = 1.0;
    kk = K;
    while (kk > last(i)) {
      chi(i,kk-1) = (1 - phi(i,kk-1-1)) + phi(i,kk-1-1) * (1 - p(i,kk-1)) * chi(i,kk+1-1);
      kk = kk - 1;
    }
    chi(i,1-1) = (1 - p(i,1-1)) * chi(i,2-1);
  }

  //// calcular la verosimiltud
  // los datos
  for(int i=0; i<I; i++){ 
    // probabilidad de sobrevivencia, que es conocido porque first<k<last
    for (int t=first(i); t<last(i); t++) {
    	nll-= log(phi(i,t-1));
    }
    // probabilidad de captura, dado viva (como CH[i,t]~bernoulli(p[i,t]);)
    for(int t=0; t< last(i); t++){
      // NA significa que no hubo esfuerzo en este periodo
      if(!isNA(CH(i,t))){
	if(CH(i,t)>=1){
	  nll-= log(p(i,t));
	} else {
	  nll-= log(1-p(i,t));
	}
      }
    }
    // probabilidad de no ser capturado despues el proximo periodo fue visto 
    nll-= log(chi(i,last(i)+1-1));
  }

  // // Probabilidad de los efectos aleatorios
  for(int t=0;t<K;t++) { // periodos
    for(int ev=0; ev<3; ev++){ // eventos
      nll-=dnorm(tau(t,0,ev), Type(0), Type(1), true);
      nll-=dnorm(tau(t,1,ev), Type(0), Type(1), true);
    }
  }
  
  vector<Type> sel_pred(lengths_pred.size());
  for(int i=0; i<sel_pred.size(); i++){
    sel_pred(i)=1/(1+exp(a+b*lengths_pred(i)));
  }
  // vector<Type> catchabilityM_pred(esfuerzo_pred.size());
  // vector<Type> catchabilityH_pred(esfuerzo_pred.size());
  // for(int i=0; i<esfuerzo_pred.size(); i++){
  //   // dado el promedio de k
  //   catchabilityM_pred(i)=1-exp(-exp(mu_tauM)*esfuerzo_pred(i));
  //   catchabilityH_pred(i)=1-exp(-exp(mu_tauH)*esfuerzo_pred(i));
  // }
  
  // reportando
  vector<Type> kM1(K);
  vector<Type> kM2(K);
  vector<Type> kM3(K);  
  vector<Type> kH1(K);
  vector<Type> kH2(K);
  vector<Type> kH3(K);  
  kM1=k.col(0).col(0);
  kM2=k.col(1).col(0);
  kM3=k.col(2).col(0);
  kH1=k.col(0).col(1);
  kH2=k.col(1).col(1);
  kH3=k.col(2).col(1);
  ADREPORT(kM1);
  ADREPORT(kM2);
  ADREPORT(kM3);
  ADREPORT(kH1);
  ADREPORT(kH2);
  ADREPORT(kH3);
  ADREPORT(M);
  ADREPORT(r);
  // ADREPORT(a);
  // ADREPORT(b);
  //  ADREPORT(sel_pred);
  // ADREPORT(catchabilityM_pred);
  // ADREPORT(catchabilityH_pred);
  // REPORT(p);
  // REPORT(phi);
  // REPORT(CH);
  return(nll);
}
// final del archivo
