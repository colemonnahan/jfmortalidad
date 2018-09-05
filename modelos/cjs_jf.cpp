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
  DATA_INTEGER(I);              // numero de los individuos
  DATA_INTEGER(K);              // numero de los periodos
  DATA_MATRIX(CH);		// la historia de las capturas (0/1) [IxK]
  DATA_IVECTOR(last);		// el ultimo periodo visto
  DATA_MATRIX(counts);		// numeros de recapturas (individos x periodos)
  DATA_VECTOR(effort);		// la esfuerza por cada periodo
  
  // efectos fijos
  PARAMETER(logM);		// mortalidad de naturleza
  PARAMETER(logr);		// effecto de los numeros de capturas
  PARAMETER(logk);		// effecto de esfuerza de probabilidad de captura
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
  Type k=exp(logk);
  // TMB usa indices de 0, no de 1, entoces tenemos que estar cuidadoso, y
  // uso un "-1" para ser claro que pasa.
  for(int i=0; i<I; i++){ // iterando sobre cada individuos
    // inicializacion estan vivo en periodo uno
    phi(i,1-1)=exp(-M-counts(i,1-1)*r);
    p(i,1-1)=1;
    for(int t=1; t<K; t++) {
      // calcular prob. capturas como una funcion de efectos y covariables
      p(i,t) = 1*(1-exp(-k*effort(t)));
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
  return(nll);
}
// final del archivo
