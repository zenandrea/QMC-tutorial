#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define K_el 1
#define PI         3.14159265
#define printall 0
#define debug 1
#define maxNwf 100 //total allocated memory for the walkers' positions
#define stepsperblock 5 // controls how often the walker distribution is saved in a text file
#define constraint 10  //it must be >=2 and <<Nw
//#define Tequil 5.0  // in atomic units

void GeneraWalkerG(long int D, long int Nw, double **w);                    //initialises walker positions with Gaussian distribution
void GeneraWalkerU(long int D, long int Nw, double **w);                    //initialises walker positions with uniform distribution
double V_pot(long int D, long int tipo, double *Vpar, double *x);                    //calculates the potential felt by the i-th walker
void Propagazione(long int D,long int N, double **w, double tau);            //each walker does a random walk step
void Branching(long int D, long int *Nwt,long int Nw, double tau, double *A, double *B, double Et, double **w); //birth-death process
double norm(double mean, double std_dev); //generates Gaussian numbers with the box-muller algorithm
double rand_val(int seed);         // function similar to rand, generates a random double between 0 and 1
void f_onda(long int Nbin, double Rmax,long int Nwt,long int D, double **w, double *); //it is used to generate a normalised histogram to represent the wavefunction
double binning(double *data , long int numMeas);                        //function used for reblocking calculation
void f_pot(long int Nbin, double Rmax, long tipo, double* Vpar, double *pot); //it is used to generate a normalised histogram to represent the wavefunction

//to compile: gcc -lm prog.c -o prog

//Alfonso Annarelli, alfonso.annarelli@gmail.com




int main(){


  double tau, Tequil, Et, Ebest, E;
  long int i, j, l, D, Nsteps, Nwt, Nw, tipo, n, distr, NstepsEq, nblock;
  double **w, **wnew;
  double *A, *B, alfa, beta;
  double *E_block, delta;
  double *phi, *phi_eq, Rmax=10.0;
  long int Nbin=1000;
  long int numMeas, binSize, flag;
  double *Vpar, *pot;

  //input
  do {
    printf("Potential: type \n");
    printf("1 for the harmonic potential, \n");
    printf("2 for the Hydrogen atom, \n");
    printf("3 for V(x)=(x-a)^2*(x+a)^2, \n");
    printf("4 for V(x)=(x-a)^n*(x+a)^n. \n");
    scanf("%ld",&tipo);
  }while((tipo<1)||(tipo>4));

  if (tipo==3){
    Vpar = (double*) malloc( 2*sizeof(double) );
    printf("Potential: b*(x-a)^2*(x+a)^2 \na=?");
    scanf("%lf",&Vpar[0]);
    printf("b=?");
    scanf("%lf",&Vpar[1]);
  }

  if (tipo==4){
    Vpar = (double*) malloc( 3*sizeof(double) );
    printf("Potential: b*(x-a)^n*(x+a)^n \nn=?");
    scanf("%lf",&Vpar[2]);
    printf("a=?");
    scanf("%lf",&Vpar[0]);
    printf("b=?");
    scanf("%lf",&Vpar[1]);
  }

  printf("Enter dimensionality\n");
  scanf("%ld",&D);

  do{
    printf("Choose initial walker distribution: enter 1 for uniform, 2 for Gaussian \n");
    scanf("%ld", &distr);
  }while((distr!=1) && (distr!=2));

  printf("Enter the time step dtau in a.u.\n");
  scanf("%lf",&tau);

  printf("Enter the number of simulation steps\n");
  scanf("%ld",&Nsteps);
  
  printf("Enter the equilibration time in a.u.\n"); //usually 10-20 a.u.
  scanf("%lf",&Tequil);

  printf("Enter a target number of walkers\n");
  scanf("%ld",&Nw);
  Nwt = Nw;    // initial number of walkers = target

  //Et = 0.0;
  printf("Enter the initial value of Et in a.u. \n");
  scanf("%lf",&Et);

  do {
    printf("Parameters of the wavefunction histogram (1: (Rmax=10,Nbin=1000); 2: user's input Rmax, Nbin) \n");
    scanf("%ld",&flag);
  } while((flag!=1)&&(flag!=2));
  
  if(flag==2){
    printf("Enter the maximum distance Rmax \n");  //choose a sufficiently large value for Rmax so that the walkers do not cross the system boundary
    scanf("%lf",&Rmax);
    printf("Enter the number of bins Nbin \n");
    scanf("%ld",&Nbin);
    flag = 1;
  }




  // w: matrix of walker positions
  w = (double**) malloc (maxNwf*Nw*sizeof(double*));
  for (i=0; i<Nw; i++) {
    w[i]=(double*) malloc (D*sizeof(double));
  }
  printf("memory has been allocated \n");

  // initialisation of random functions
  rand_val(time(NULL));
  srand(time(NULL));

  if(distr==1){
    GeneraWalkerU(D,Nw,w);
  }
  else if(distr==2){
    GeneraWalkerG(D,Nw,w);
  }
  else{
    printf("Initial distribution not correctly defined \n");
    exit(EXIT_FAILURE);
  }

  printf("Walkers has been generated\n");


  Ebest = 0.0;
  E=0.0;

  //output files
  FILE *fe;
  FILE *fb;
  FILE *fp;
  FILE *ff;
  FILE *fpot;

  fe=fopen("energie.txt","w"); // energy at each step
  fprintf(fe,"#n \t  T  \t \t  Et \t \t  Ebest \n");

  fp=fopen("popolazione.txt","w");
  fprintf(fp, "#n \t T \t\t Nwt \n");

  fb=fopen("reblocking.txt","w");
  fprintf(fb, "binSize \t standard error \t numMeas \n");

  if (flag==1) {
    ff=fopen("f_onda_equil.txt","w");  // equilibrium wavefunction
    fprintf(ff,"#R \t \t Phi \n");
    phi = (double*) malloc( Nbin * sizeof(double) );
    phi_eq = (double*) malloc( Nbin * sizeof(double) );
    if (D==1) {
      fpot = fopen("pot.txt","w");
      fprintf(fpot, "# x \t \t pot(x)\n");
      pot = (double*) malloc( 2*Nbin * sizeof(double) );
      f_pot(Nbin,Rmax,tipo,Vpar,pot);  //fill "pot" with the Nbin values of the potential and then with Nbin associated positions
      for( i=0; i<Nbin; i++)
        fprintf(fpot,"%lf\t%lf\n",pot[i+Nbin],pot[i]);
      fprintf(fpot,"\n");
      fclose(fpot);
      free(pot); //release the allocated memory
    }
  }

  NstepsEq = (long)( Tequil / tau );
  numMeas = Nsteps - NstepsEq;
  if (numMeas>0){
    E_block = (double*) malloc( numMeas * sizeof(double) );  //to record the Et value at each step
  } else {
    E_block = (double*) malloc(0);
    numMeas = 0.0;
  }

  //start of the main loop
  for( n=0; n<Nsteps; n++ ){

    // A: potentials before the diffusion
    // B: potentials after the diffusion
    printf("%6ld  Nwt=%ld ", n, Nwt);
    A = (double*) malloc (Nwt*sizeof(double)); //memory allocation must be done at each step because Nwt is variable

    for(i=0;i<Nwt;i++){
      if(*(w+i)!=NULL){
        A[i]=V_pot(D, tipo, Vpar, w[i]);
      }
    }

    if(flag==1){
      delta=Rmax/(double)Nbin;
      if(n<=Tequil/tau && ((n % stepsperblock)==0)){  //print the wavefunction during the equilibration phase
        fprintf(ff,"# wavefunction at step number %ld \n",n);
        fprintf(ff, "\n");
        f_onda(Nbin,Rmax,Nwt,D,w,phi); //creates the histogram in the ‘phi’ vector
        if(D>1){
          for(i=0;i<Nbin;i++){
            fprintf(ff,"%f \t %f \n",delta*(i+0.5), phi[i]);
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

    fprintf(fp,"%ld\t%f\t%ld\n",n,n*tau, Nwt);

    Propagazione(D,Nwt,w,tau);          //Diffusion step

    B=(double*) malloc (Nwt*sizeof(double));
    for(i=0;i<Nwt;i++){
      if(*(w+i)!=NULL){
        B[i]=V_pot(D, tipo, Vpar, w[i]);
      }
    }

    //printf(">>> diffusion step n. %ld \n",n);

    Branching(D, &Nwt, Nw, tau, A,B, Et,w);       //Branching step

    //printf("branching step n. %ld \n",n);
    
    printf("newNwt=%ld  ", Nwt);

    free(A);
    free(B);


    // update the value of the energy estimator <Et>
    if ( n<=NstepsEq ){
      nblock = n % stepsperblock;
      Ebest = ( Et + Ebest * nblock ) / (nblock+1);
    } else {
      nblock = n - NstepsEq;
      Ebest = ( Et + Ebest * nblock ) / (nblock+1);
    }

    // update Et
    Et = Ebest-(1.0/tau)*log((double) Nwt / (double) Nw);

    // fill the array with Et values after the equilibration phase
    if ( n>=NstepsEq ) 
		E_block[n-NstepsEq]=Et;

    printf("Et= %f \t Ebest= %f \n", Et,Ebest);

    fprintf(fe,"%ld\t%f\t%f\t%f\n", n, n*tau, Et, Ebest);



    if(flag==1){
      if(n>=NstepsEq){
        f_onda(Nbin,Rmax,Nwt,D,w,phi);
        for(i=0;i<Nbin;i++){
          phi_eq[i]+=phi[i];
        }
      }
    }

  } // end of the main loop
  free(w);
  if (tipo==3 || tipo==4){
    free(Vpar);
  }
  
  if(flag==1){
    fclose(ff);
    ff = fopen("f_onda.txt","w");
    printf("Saving the equilibrium wavefunction \n" );
    // wavefunction calculation at equilibrium
    for(i=0;i<Nbin;i++){
      phi_eq[i]=phi_eq[i]/(Nsteps-NstepsEq);  //calculated as the average of the wavefunctions at each step
    }
    
    if(D>1){
      for(i=0;i<Nbin;i++){
        fprintf(ff,"%f\t%f\n", delta*(i+0.5), phi_eq[i]);
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

  fclose(fe);
  fclose(fp);
  
  //reblocking, see J. Chem. Phys. 91, 461–466 (1989)
  binSize = 1;
  while (numMeas >= 100) {
    fprintf(fb,"%ld \t\t %lf \t\t %ld \n", binSize, binning(E_block, numMeas), numMeas);
    binSize *= 2;
    numMeas /= 2;
  }
  fclose(fb);
  free(E_block);
  printf("Reblocking completed \n");
  free(phi_eq);
  free(phi);

  printf("Programme successfully executed \n");

  return 0;
}






//functions

void GeneraWalkerG(long int D, long int Nw, double **w){
  long int i, j;
  for(i=0;i<Nw;i++){
    for(j=0;j<D;j++){
      w[i][j] = norm(0,1);		//the position of each walker is drawn from a gaussian distribution
	}
  }
}

void GeneraWalkerU(long int D, long int Nw, double **w){
  long int i, j;
  for(i=0;i<Nw;i++){
    for(j=0;j<D;j++){
      w[i][j]=(rand()/(double)(RAND_MAX)-0.5);		// drawn from an uniform distribution
    }
  }
}

double V_pot(long int D, long int tipo, double *Vpar, double *x){
  double a=0;
  long int j;

  if(tipo == 1){
    a=0;
    for(j=0;j<D;j++){
      a+=(0.5)*K_el*(pow((x[j]),2));
    }
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
    printf("potential not correctly defined \n");
    exit(EXIT_FAILURE);
  }

  return a;
}

void Propagazione(long int D,long int Nwt, double **w, double tau){  //diffusion function
  long int i,j;
  for(i=0;i<Nwt;i++){
    for(j=0;j<D;j++){
      w[i][j]+= sqrt(tau)*norm(0,1);
    }
  }
}

void Branching(long int D, long int *pNwt, long int Nw, double tau, double *A, double *B, double Et, double **w){
  long int i,j, k, Mi, Nnew=0.0, inew;
  double wi;
  double **wnew;
  long int *M;

  M = (long int *) malloc( *pNwt * sizeof(long int) ); //it contains the number of copies generated by each walker

  for( i=0 ; i < *pNwt; i++ ){
    wi = exp(-( tau * ( A[i] + B[i] - 2*Et )/2) );
    if(wi>=constraint)
    	wi=constraint;	//constraint for population "explosions", useful for potential with divergences
    Mi = (long int)(wi+(rand()/(double)RAND_MAX));
    if (printall) printf("Mi= %ld \n",Mi);
    M[i] = Mi;
    Nnew += Mi;
  }
  //printf("Nwt=%ld\tNnew=%ld\n", (*pNwt), Nnew );

  if ( Nnew <= maxNwf * (*pNwt) ) {
    wnew = (double**) malloc( Nnew * sizeof(double*) );  //same size of w
  } else {
    printf("Error: the programme has run out of allocated memory, please increase maxNwf.\n");
    exit(EXIT_FAILURE);
  }

  inew = 0;
  for( i=0; i < *pNwt; i++ ){
    if ( M[i]==0 ) {
      // do not replicate walkers
    }
    else {
      wnew[inew] = w[i];
      inew++;
      for ( j=1; j<M[i]; j++ ) {  //it records the position of copies generated from the same walker in wnew
        wnew[inew] = (double*) malloc( D*sizeof(double) );
        for (k=0;k<D;k++)
        	wnew[inew][k] = w[i][k];
        inew++;
      }
    }
  }

  if (inew != Nnew){
    printf("Error: inew after the cycle is different from Nnew");
    exit(EXIT_FAILURE);
  }

// free the memory used by wnew
  for( i=0; i<Nnew; i++ ){
    w[i] = wnew[i];
  }
  free(wnew);
  free(M);

  *pNwt = Nnew;
  if(*pNwt == 0){
    printf("Error: no walkers left \n");
    exit(EXIT_FAILURE);
  }
}



double norm(double mean, double std_dev)
{
  double   u, r, theta;           // Variables for Box-Muller method
  double   x;                     // Normal(0, 1) random variable (rv)
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

  for(i=0;i<Nbin;i++){		//initialize the histogram
    histo[i]=0;
  }

  if(D==1){
    for(i=0;i<Nwt;i++){
      bin=(long int)(w[i][0]*Nbin/Rmax+(Nbin/2.0)); 	//it reads the position of one walker at a time and increments the corresponding column in the histogram by 1.
      if (bin<0 || bin>Nbin-1){
        printf("bin=%ld, w=%lf, please set a larger Rmax. \n",bin,w[i][0]);
        exit(EXIT_FAILURE);
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

      bin=(long int)(R_i*Nbin/Rmax);
      if(bin<Nbin){
        histo[bin]++;  //it is equivalent to integrating over all (polar) coordinates except the radial one
      } else {
        printf("Please set a larger Rmax. \n");
        exit(EXIT_FAILURE);
      }
    }
  }
  
  if(D==1){
    for(i=0;i<Nbin;i++){
      normalizzazione+=pow(histo[i],2);
    }

    normalizzazione=sqrt(normalizzazione*Rmax/Nbin); //Rmax/Nbin=delta

    for(i=0;i<Nbin;i++){
      histo[i]=histo[i]/normalizzazione;
    }
  }
  if(D==3){
    normalizzazione=Nwt*(Rmax/Nbin)/(8*sqrt(PI)); //Only applies to the hydrogen atom, specifically normalising the distribution 4pi*r^2*Phi(r)
    for(i=0;i<Nbin;i++){
      histo[i]=histo[i]/normalizzazione;
    }
  }
}

void f_pot(long int Nbin, double Rmax, long tipo, double* Vpar, double *pot){
  double delta, R;
  long int i, j, bin;

  delta = Rmax / Nbin;
  for(i=0;i<Nbin;i++){
    R = (i+0.5-Nbin/2.0) * delta;
    pot[i] = V_pot( 1, tipo, Vpar, &R );
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
    data[i] = 0.5 * (data[2 * i] + data[2 * i + 1]);  //prepare the ‘data’ vector for the next binning iteration, each time halving the portion to be analysed
  }
  if (2 * i < numMeas) {
    mean += data[2 * i];
    variance += data[2 * i] * data[2 * i];
  }
  mean /= numMeas; variance /= numMeas;
  return sqrt((variance - mean * mean ) / (numMeas - 1));
}
