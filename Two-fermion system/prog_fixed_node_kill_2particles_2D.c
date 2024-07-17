#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


//#define K_el 1
#define PI         3.14159265
#define printall 0
#define debug 1
#define maxNwf 100 //total allocated memory for the walkers' positions
#define stepsperblock 10000 //20
#define constraint 10  //it must be >=2 and <<Nw
//#define Tequil 5.0  // atomic units
//#define offset 1.0  // per terzo potenziale

void GeneraWalkerG(long int D, long int Nw, double **w);                    //inizializza le posizioni dei walker con distribuzione gaussiana
void GeneraWalkerU(long int D, long int Nw, double **w);                    //inizializza le posizioni dei walker con distribuzione uniforme
double V_pot(long int D, long int tipo, double K_el, double K_int, double *Vpar, double *x);                    //calcola il potenziale sentito dall'iesimo walker
void Propagazione(long int D,long int N, double **w, double np, double tau);            //ogni walker fa uno step di random walk
void Branching(long int D, long int *Nwt,long int Nw, double tau, double *A, double *B, double Et, double **w); //processo di nascita e di morte
double norm(double mean, double std_dev); //genera numeri gaussiani con la box-muller
double rand_val(int seed);         // funzione simile a rand, genera un double casuale compreso fra 0 e 1
void f_onda(long int Nbin, double Rmax,long int Nwt,long int D, double **w, double *); //serve a generare un istogramma normalizzato per rappresentare la f.onda
double binning(double *data , long int numMeas);                        //funzione utilizzata per il calcolo del reblocking
void f_pot(long int Nbin, double Rmax, long tipo, double K_el, double K_int, double* Vpar, double *pot); //serve a generare un istogramma normalizzato per rappresentare la f.onda
double psi_T(double *w, double np);

//per compilare: gcc -lm prog_fixed_node_kill_2particles_2D.c -o prog_fixed_node_kill_2particles_2D

	




int main(){


  double tau, Tequil, Et, Ebest, E;
  long int i, j, l, D, Nsteps, Nwt, Nw, tipo, n, distr, NstepsEq, nblock;
  double **w, **wnew, x_node[2]={0}, K_el=1.0, K_int=0.0, np=0.0;
  double *A, *B, alfa, beta;
  double *E_block, delta;
  double *phi, *phi_eq, Rmax=10.0;
  long int Nbin=1000;
  long int numMeas, binSize, flag;
  double *Vpar, *pot;

  //input
  do {
    printf("Potenziale: digita \n");  //dato il resto del programma, per D>1 funzionano solo potenziali radiali
    printf("1 per il potenziale armonico, \n");
    printf("2 per l'atomo di idrogeno, \n");  //circa 16 milioni di passi per l'atomo H, poi dà segmentation fault
    printf("3 per (x-a)^2*(x+a)^2, \n");
    printf("4 per (x-a)^n*(x+a)^n. \n");
    scanf("%ld",&tipo);
  }while((tipo<1)||(tipo>4));

  if (tipo==3){
    Vpar = (double*) malloc( 2*sizeof(double) );
    printf("Potenziale: b*(x-a)^2*(x+a)^2 \na=?");
    scanf("%lf",&Vpar[0]);
    printf("b=?");
    scanf("%lf",&Vpar[1]);
  }

  if (tipo==4){
    Vpar = (double*) malloc( 3*sizeof(double) );
    printf("Potenziale: b*(x-a)^n*(x+a)^n \nn=?");
    scanf("%lf",&Vpar[2]);
    printf("a=?");
    scanf("%lf",&Vpar[0]);
    printf("b=?");
    scanf("%lf",&Vpar[1]);
  }

  printf("Inserire dimensione\n");
  scanf("%ld",&D);
  

  if(tipo==1){
    if(D==1){
      printf("Inserisci la posizione del nodo \n");
      scanf("%lf", &x_node[0]);
	}
    else if(D==2){
      printf("Inserisci i parametri della superfice nodale x1=a*x2+b\n");
      printf("a=?");
      scanf("%lf", &x_node[0]);
      printf("b=?");
      scanf("%lf", &x_node[1]);
	}
    printf("Inserisci il valore della costante elastica K_el > 0 \n");
    scanf("%lf", &K_el);
    printf("Inserisci il valore della costante di interazione K_int >= 0 \n");
    scanf("%lf", &K_int);
    printf("Inserisci il parametro np della funzione d'onda di prova \n");
    scanf("%lf", &np);
  }
 
  do{
    printf("Scegli la distribuzione iniziale dei walker: digita 1 per uniforme, 2 per gaussiana \n");
    scanf("%ld", &distr);
  }while((distr!=1) && (distr!=2));

  printf("Inserire il time step tau in a.u.\n");
  scanf("%lf",&tau);

  printf("Inserire il numero di step\n");
  scanf("%ld",&Nsteps);
  
  printf("Inserire il tempo di equilibrazione in a.u.\n"); //dell'ordine di 10 a.u. o più
  scanf("%lf",&Tequil);

  printf("Inserire un numero di walker target\n");
  scanf("%ld",&Nw);
  Nwt = Nw;    // initial number of walkers = target

  //Et = 0.0;
  printf("Inserisci il valore di Et iniziale in a.u. \n");
  scanf("%lf",&Et);

  /*
  printf("Inserisci flag funzione d'onda (1: si (Rmax=10,Nbin=1000); 2: si input Rmax,Nbin;  else:altro) \n");
  scanf("%ld",&flag);
  if(flag==2){
    printf("inserisci distanza massima Rmax \n");
    scanf("%lf",&Rmax);
    printf("Inserisci numero di bin Nbin (pari) \n");
    scanf("%ld",&Nbin);
    flag = 1;
  }
  */


  // w: matrice con le posizioni dei walker
  w = (double**) malloc (maxNwf*Nw*sizeof(double*));
  for (i=0; i<Nw; i++) {
    w[i]=(double*) malloc (D*sizeof(double));
  }
  printf("memoria allocata \n");

  // inizializzazione funzioni random
  rand_val(time(NULL));
  srand(time(NULL));

  if(distr==1){
    GeneraWalkerU(D,Nw,w);
  }
  else if(distr==2){
    GeneraWalkerG(D,Nw,w);
  }
  else{
    printf("Distribuzione iniziale non definita correttamente \n");
    exit(EXIT_FAILURE);
  }

  printf("walker generati \n");


  Ebest = 0.0;
  E=0.0;

  //file di output
  FILE *fe;
  FILE *fb;
  FILE *fp;
  //FILE *ff;
  //FILE *fpot;

  fe=fopen("energie.txt","w");
  fprintf(fe,"#n \t  T  \t \t  Et \t \t  Ebest \n");

  fp=fopen("popolazione.txt","w");
  fprintf(fp, "#n \t T \t\t Nwt \n");

  fb=fopen("reblocking.txt","w");
  fprintf(fb, "binSize \t standard error \t numMeas \n");

  /*
  if (flag==1) {
    ff=fopen("f_onda_seq.txt","w");
    
    if (tipo==1 && D==2){
      fprintf(ff,"#x1 \t \t x2 \t \t Phi \n");
      phi = (double*) malloc( pow(Nbin, D) * sizeof(double) );
      phi_eq = (double*) malloc( pow(Nbin, D) * sizeof(double) );
	}
	else{
	  fprintf(ff,"#R \t \t Phi \n");
      phi = (double*) malloc( Nbin * sizeof(double) );
      phi_eq = (double*) malloc( Nbin * sizeof(double) );
      if (D==1) {
        fpot = fopen("pot.txt","w");
        fprintf(fpot, "# x \t \t pot(x)\n");
        pot = (double*) malloc( 2*Nbin * sizeof(double) );
        f_pot(Nbin,Rmax,tipo, K_el, K_int, Vpar,pot);  //fill "pot" with the Nbin values of the potential and then with Nbin associated positions
        for( i=0; i<Nbin; i++)
          fprintf(fpot,"%lf\t%lf\n",pot[i+Nbin],pot[i]);
        fprintf(fpot,"\n");
        fclose(fpot);
        free(pot); //release the allocated memory
      }
	}
  }
  */
  
  NstepsEq = (long)( Tequil / tau );
  numMeas = Nsteps - NstepsEq;
  if (numMeas>0){
    E_block = (double*) malloc( numMeas * sizeof(double) );  //per registrare il valore di Et ad ogni passo
  } else {
    E_block = (double*) malloc(0);
    numMeas = 0.0;
  }

  //inizio ciclo
  for( n=0; n<Nsteps; n++ ){

    // A: potenziali prima della diffusione
    // B: potenziali dopo la diffusione
    printf("%6ld  Nwt=%ld ", n, Nwt);
    A = (double*) malloc (Nwt*sizeof(double)); //l'allocazione va fatta ad ogni step perchè Nwt è variabile

    for(i=0;i<Nwt;i++){
      if(*(w+i)!=NULL){
        A[i]=V_pot(D, tipo, K_el, K_int, Vpar, w[i]);
      }
    }

    /*
    if(flag==1){
      delta=Rmax/(double)Nbin;
      if((n % stepsperblock)==0){  //costruisce la funzione d'onda durante l'intera simulazione
        fprintf(ff,"# funzione d'onda al passo numero %ld \n",n);
        fprintf(ff, "\n");
        f_onda(Nbin,Rmax,Nwt,D,w,phi); //costruisce l'istogramma nel vettore "phi"
        if(D>=3){
          for(i=0;i<Nbin;i++){
            fprintf(ff,"%f \t %f \n",delta*(i+0.5), phi[i]);
          }
          fprintf(ff, "\n");
        }
        else if(D==2){
          for(i=0;i<pow(Nbin, D);i++){
            fprintf(ff,"%f \t %f \t %f \n", delta*(i/Nbin+0.5-(Nbin/2)), delta*(i%Nbin+0.5-(Nbin/2)), phi[i]);
          }
          fprintf(ff, "\n");
		}
        else{
          for(i=0;i<Nbin;i++){
            fprintf(ff,"%f \t %f \n", delta*(i+0.5-(Nbin/2)), phi[i]);
          }
          fprintf(ff, "\n");
        }
      }
    }
    */

    fprintf(fp,"%ld\t%f\t%ld\n",n,n*tau, Nwt);		//population evolution

    Propagazione(D,Nwt,w,np,tau);          //Diffusion step

    B=(double*) malloc (Nwt*sizeof(double));
    for(i=0;i<Nwt;i++){
      if(*(w+i)!=NULL){
        B[i]=V_pot(D, tipo, K_el, K_int, Vpar, w[i]);
      }
    }

    //printf(">>> step propagazione n. %ld \n",n);

    Branching(D, &Nwt, Nw, tau, A,B, Et,w);       //Branching step

    //printf("step branching n. %ld eseguito \n",n);
    
    printf("newNwt=%ld  ", Nwt);

    free(A);
    free(B);


    // definire nuova energia Et
    if ( n<=NstepsEq ){
      nblock = n % stepsperblock;
      Ebest = ( Et + Ebest * nblock ) / (nblock+1);
    } else {
      nblock = n - NstepsEq;
      Ebest = ( Et + Ebest * nblock ) / (nblock+1);
    }

    // Aggiornare Et
    Et = Ebest-(1.0/tau)*log((double) Nwt / (double) Nw);       //ricalcolo di Et all'equilibrio

    // riempio array con valore Et dopo equilibrazione
    if ( n>=NstepsEq ) 
		E_block[n-NstepsEq]=Et;

    printf("Et= %f \t Ebest= %f \n", Et,Ebest);

    //fprintf(fe,"%f,%f,%f\n",n*tau, Et,Ebest);
    fprintf(fe,"%ld\t%f\t%f\t%f\n", n, n*tau, Et, Ebest);


    /*
    if(flag==1){
      if(n>=NstepsEq && D<=2){  //è corretto che ci sia >=, in modo che mediamo su Nsteps-NstepsEq campionamenti di phi (?)
        f_onda(Nbin,Rmax,Nwt,D,w,phi);
        for(i=0;i<pow(Nbin,D);i++){
          phi_eq[i]+=phi[i];
        }
      } else {
        f_onda(Nbin,Rmax,Nwt,D,w,phi);
        for(i=0;i<Nbin;i++){
          phi_eq[i]+=phi[i];
	    }
      }
    }
    */

  } // fine ciclo DMC
  
  free(w);
  if (tipo==3 || tipo==4){
    free(Vpar);
  }
  
  /*
  if(flag==1){
    fclose(ff);
    ff = fopen("f_onda.txt","w");
    printf("Scrivo funzione d'onda all'equilibrio \n" );
    // calcolo funzione d'onda all'equilibrio
    if(D<=2){
      for(i=0;i<pow(Nbin,D);i++){
        phi_eq[i]=phi_eq[i]/(Nsteps-NstepsEq);  //calcolata come la media delle funzioni d'onda ad ogni passo
      }
    } else {
      for(i=0;i<Nbin;i++){
        phi_eq[i]=phi_eq[i]/(Nsteps-NstepsEq);  //calcolata come la media delle funzioni d'onda ad ogni passo
      }
	}
    
    
    if(D>=3){
      for(i=0;i<Nbin;i++){
        fprintf(ff,"%f\t%f\n", delta*(i+0.5), phi_eq[i]);
      }
      fprintf(ff, "\n");
    }
    else if(D==2){
      for(i=0;i<pow(Nbin, D);i++){
        fprintf(ff,"%f \t %f \t %f \n", delta*(i/Nbin+0.5-(Nbin/2)), delta*(i%Nbin+0.5-(Nbin/2)), phi[i]);
      }
      fprintf(ff, "\n");
	}
    else{
      for(i=0;i<Nbin;i++){
        fprintf(ff,"%f\t%f\n", delta*(i+0.5-(Nbin/2)), phi_eq[i]);
      }
      fprintf(ff, "\n");
    }
    fclose(ff);
  }
  */

  fclose(fe);
  fclose(fp);
  
  //reblocking
  binSize = 1;
  while (numMeas >= 100) {

    fprintf(fb,"%ld \t\t %lf \t\t %ld \n", binSize, binning(E_block, numMeas), numMeas);
    binSize *= 2;
    numMeas /= 2;
  }
  fclose(fb);
  free(E_block);
  printf("Reblocking completato \n");
  //free(phi_eq);  //WARNING! Il programma non dealloca correttamente phi e phi_eq se i walker attraversano i bordi e si arresta in questo punto in maniera anomala, quindi lasciare queste due istruzioni sempre per ultime
  //free(phi);

  printf("Programma eseguito \n");

  return 0;
}






//functions

void GeneraWalkerG(long int D, long int Nw, double **w){
  long int i, j;
  for(i=0;i<Nw;i++){
    for(j=0;j<D;j++){
      w[i][j] = norm(0,1);		//norm(-2.0, 2.0/5.0);
	  /*if ((rand()/(double)(RAND_MAX)-0.5)>0)
        w[i][j] = norm(+4.0, 4.0/5.0);
      else
        w[i][j] = norm(-4.0, 4.0/5.0);		//doppia gaussiana */
	}
  }
}

void GeneraWalkerU(long int D, long int Nw, double **w){
  long int i, j;
  for(i=0;i<Nw;i++){
    for(j=0;j<D;j++){
      w[i][j]=5*(rand()/(double)(RAND_MAX)-0.5);		//(rand()/(double)(RAND_MAX)-0.5);
    }
  }
}

double V_pot(long int D, long int tipo, double K_el, double K_int, double *Vpar, double *x){
  double a=0;
  long int j;

  if(tipo == 1){
    a=0;
    for(j=0;j<D;j++){
      a+=(0.5)*K_el*(pow((x[j]),2));
    }
    a+=0.5*K_int*pow(x[0]-x[1],2);		//K_int>0
  }
  else if(tipo == 2){
    a=0;
    for(j=0;j<D;j++){
      a+=(pow((x[j]),2));
    }
    a=-1.0/sqrt(a);
  }
  else if(tipo == 3){
    a=0;
    for(j=0;j<D;j++){
      a+= Vpar[1] * (pow((x[j])-Vpar[0],2)) * (pow((x[j])+Vpar[0],2));
    }
  }
  else if(tipo == 4){
    a=0;
    for(j=0;j<D;j++){
      a+= Vpar[1] * (pow((x[j])-Vpar[0],Vpar[2])) * (pow((x[j])+Vpar[0],Vpar[2]));
    }
  }
  else{
    printf("potenziale non definito correttamente \n");
    exit(EXIT_FAILURE);
  }

  return a;
}

double psi_T(double *w, double np){		//trial wavefunction, remember to insert omega frequencies to get the correct dimensionality
	double p;
	//p=exp(-0.5*(pow(w[0],2)+pow(w[1],2)+pow(w[2],2)+pow(w[3],2)))*(w[0]-w[2]+(w[1]-w[3])*pow(sqrt(pow(w[1]-w[3],2)),np));
	p=exp(-0.5*(pow(w[0],2)+pow(w[1],2)+pow(w[2],2)+pow(w[3],2)))*(w[0]-w[2]+w[1]-w[3]+np*(pow(w[1],2)-pow(w[3],2)));
	return p;
}

void Propagazione(long int D,long int Nwt, double **w, double np, double tau){
  long int i,j;
  double w_temp[D];
  for(i=0;i<Nwt;i++){
    for(j=0;j<D;j++){
      w_temp[j]=w[i][j];		//fixed node with rejection only in 1D
      w[i][j]+= sqrt(tau)*norm(0,1);
    }
    if (psi_T(w[i],np)*psi_T(w_temp,np)<0)
      w[i][0]=1000;		//w[i]=NULL  forse dà problemi lasciare puntatori vuoti, si potrebbero accumulare (?)	//w[i][j]=10000; impractical way to almost surely kill a walker
  }
}

void Branching(long int D, long int *pNwt, long int Nw, double tau, double *A, double *B, double Et, double **w){
  long int i,j, k, Mi, Nnew=0.0, inew;
  double wi;
  double **wnew;
  long int *M;

  M = (long int *) malloc( *pNwt * sizeof(long int) ); //contiene il numero di copie generate da ogni walker

  for( i=0 ; i < *pNwt; i++ ){
    if(w[i]!=NULL && w[i][0]!=1000){
      wi = exp(-( tau * ( A[i] + B[i] - 2*Et )/2) );
    }
    else 
	  wi=0;  //non sempre cancella i walker correttamente
    if(wi>=constraint)
      wi=constraint;	//constraint for population "explosions", useful for potential with divergences
    Mi = (long int) floor(wi+(rand()/(double)RAND_MAX)); //il cast a long int tronca sempre per difetto il numero
    if (printall) printf("Mi= %ld \n",Mi);
    M[i] = Mi;
    Nnew += Mi;
  }
  //printf("Nwt=%ld\tNnew=%ld\n", (*pNwt), Nnew );

  if ( Nnew <= maxNwf * (*pNwt) ) {
    wnew = (double**) malloc( Nnew * sizeof(double*) );  //stesse dimensioni di w
  } else {
    printf("Errore, il programma ha utilizzato tutta la memoria allocata, aumentare maxNwf.\n");
    exit(EXIT_FAILURE);
  }

  inew = 0;
  for( i=0; i < *pNwt; i++ ){
    if ( M[i]==0 ) {
      // do not replicate walkers
      free( w[i] ); //serve perchè altrimenti si accumulano puntatori vuoti
    }
    else {
      wnew[inew] = w[i];  //è bene fare così perchè abbiamo già una porzione di memoria allocata per il primo walker
      inew++;
      for ( j=1; j<M[i]; j++ ) {  //registra la posizione delle copie dello stesso walker in wnew, in sequenza
        wnew[inew] = (double*) malloc( D*sizeof(double) );
        for (k=0;k<D;k++)
        	wnew[inew][k] = w[i][k];
        inew++;
      }
    }
  }

  if (inew != Nnew){
    printf("Errore: inew dopo ciclo diverso da Nnew");
    exit(EXIT_FAILURE);
  }

// libero memoria usata da wnew
  for( i=0; i<Nnew; i++ ){
    w[i] = wnew[i];
  }
  free(wnew);
  free(M);

  *pNwt = Nnew;
  if(*pNwt == 0){
    printf("Errore, tutti i walker sono morti \n");
    exit(EXIT_FAILURE);
  }
}



double norm(double mean, double std_dev)
{
  double   u, r, theta;           // Variables for Box-Muller method
  double   x;                     // Normal(0, 1) rv
  double   norm_rv;               // The adjusted normal rv

  // Generate u
  u = 0.0;
  while (u == 0.0)
  u = rand_val(0);

  // Compute r
  r = sqrt(-2.0 * log(u));

  // Generate theta
  theta = 0.0;
  while (theta == 0.0)
  theta = 2.0 * PI * rand_val(0);

  // Generate x value
  x = r * cos(theta);

  // Adjust x value for specified mean and variance
  norm_rv = (x * std_dev) + mean;

  // Return the normally distributed RV value
  return(norm_rv);
}


double rand_val(int seed) //Multiplicative LCG, Schrage’s algorithm
{
  const long  a =      16807;  // Multiplier
  const long  m = 2147483647;  // Modulus, 2^31-1, Mersenne prime
  const long  q =     127773;  // m div a
  const long  r =       2836;  // m mod a
  static long x;               // Random int value, "static" means its value is preserved even after the end of the function
  long        x_div_q;         // x divided by q
  long        x_mod_q;         // x modulo q
  long        x_new;           // New x value

  // Set the seed if argument is non-zero and then return zero
  if (seed > 0)
  {
    x = seed;
    return(0.0);
  }

  // RNG using integer arithmetic
  x_div_q = x / q;
  x_mod_q = x % q;
  x_new = (a * x_mod_q) - (r * x_div_q);
  if (x_new > 0)
  x = x_new;
  else
  x = x_new + m;

  // Return a random value between 0.0 and 1.0
  return((double) x / m);
}

void f_onda(long int Nbin, double Rmax,long int Nwt,long int D, double **w, double *histo){
  double deltaR,R_i=0.0,ausR=0.0, normalizzazione=0.0;
  long int i, j, bin;

  if(D==2){
    for(i=0;i<pow(Nbin, D);i++){
      histo[i]=0;
    }
  }
  else{
    for(i=0;i<Nbin;i++){
      histo[i]=0;
    }
  }

  if(D==1){
    for(i=0;i<Nwt;i++){
      bin=(long int)floor(w[i][0]*Nbin/Rmax+(Nbin/2.0)); //legge la posizione di un walker alla volta e incrementa di 1 la colonna corrispondente nell'istogramma. Probabilmente questo tipo di assegnazione (da posizione a bin) è asimmetrica, inducendo un bias nell'energia calcolata
      if (bin<0 || bin>Nbin-1){
        printf("bin=%ld, w=%lf, Aumenta Rmax \n",bin,w[i][0]);  //è possibile che bin=-1 o bin=Nbin perchè non stiamo rigettando gli step di Diffusione attraverso i bordi
        //exit(EXIT_FAILURE);
      }
      histo[bin]++;
    }
  }
  else if(D==2){
    for(i=0;i<Nwt;i++){
      bin=(long int)floor(w[i][0]*Nbin/Rmax+(Nbin/2.0))*Nbin+(long int)floor(w[i][1]*Nbin/Rmax+(Nbin/2.0)); //legge la posizione di un walker alla volta e incrementa di 1 la colonna corrispondente nell'istogramma. Probabilmente questo tipo di assegnazione (da posizione a bin) è asimmetrica, inducendo un bias nell'energia calcolata
      if (bin<0 || bin>Nbin*Nbin-1){
        printf("bin=%ld, w=(%lf,%lf), Aumenta Rmax \n",bin,w[i][0],w[i][1]);  //se un bordo viene attraversato, il walker finisce nella prossima riga, non è propriamente corretto
        //exit(EXIT_FAILURE);
      }
      histo[bin]++;
    }
  }
  else{
    for(i=0;i<Nwt;i++){
      for(j=0;j<D;j++){
        ausR+=pow(w[i][j],2);
      }
      R_i=sqrt(ausR);
      ausR=0.0;

      bin=(long int)floor(R_i*Nbin/Rmax);
      if(bin<Nbin){
        histo[bin]++;  //equivale ad integrare su tutte le coordinate (polari) tranne quella radiale
      } else {
        printf("Aumenta Rmax \n");
        //exit(EXIT_FAILURE);
      }
    }
  }
  
  if(D==1){
    for(i=0;i<Nbin;i++){
      normalizzazione+=pow(histo[i],2);
    }

    normalizzazione=sqrt(normalizzazione*Rmax/Nbin); //Rmax/Nbin=delta, appare quando si calcola l'integrale per la normalizzazione 

    for(i=0;i<Nbin;i++){
      histo[i]=histo[i]/normalizzazione;
    }
  }
  else if(D==2){
    for(i=0;i<pow(Nbin, D);i++){
      normalizzazione+=pow(histo[i],2);
    }

    normalizzazione=sqrt(normalizzazione*pow(Rmax/Nbin,D)); //Rmax/Nbin=delta, appare quando si calcola l'integrale per la normalizzazione 

    for(i=0;i<pow(Nbin, D);i++){
      histo[i]=histo[i]/normalizzazione;
    }
  }
  
  if(D==3){
    normalizzazione=Nwt*(Rmax/Nbin)/(8*sqrt(PI)); //Vale solo per l'atomo di idrogeno, in particolare normalizza 4pi*r^2*Phi(r)
    for(i=0;i<Nbin;i++){
      histo[i]=histo[i]/normalizzazione;
    }
  }
}

void f_pot(long int Nbin, double Rmax, long tipo, double K_el, double K_int, double* Vpar, double *pot){
  double delta, R;
  long int i, j, bin;

  delta = Rmax / Nbin;
  for(i=0;i<Nbin;i++){
    R = (i+0.5-Nbin/2.0) * delta;
    pot[i] = V_pot( 1, tipo, K_el, K_int, Vpar, &R );
    pot[i+Nbin] = R;
  }
}

double binning(double *data , long int numMeas){
  long int i, tmp = numMeas / 2;
  double mean = 0.0, variance = 0.0;

  for (i = 0; i < tmp; i++) {
    mean += data[2 * i] + data[2 * i + 1];
    variance += data[2 * i] * data[2 * i] +
    data[2 * i + 1] * data[2 * i + 1];
    data[i] = 0.5 * (data[2 * i] + data[2 * i + 1]);  //in tal modo preparo il vettore "data" per la prossima chiamata di binning, dimezzando di volta in volta la porzione da analizzare
  }
  if (2 * i < numMeas) {  //nel caso non venisse preso in considerazione l'ultimo elemento del vettore "data" nel ciclo for (?), es. numMeas dispari
    mean += data[2 * i];
    variance += data[2 * i] * data[2 * i];
  }
  mean /= numMeas; variance /= numMeas;
  //double error = sqrt((variance - mean * mean ) / (numMeas - 1));
  //printf("mean=%lf \t error= %lf \t numMeas=%ld \n", mean, error, );
  return sqrt((variance - mean * mean ) / (numMeas - 1));
}
