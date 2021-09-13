/*
  This program is a 3D FDTD simulation that will
  model the fields coming off of a microstrip, onto a patch antenna
  Absorbing Boundary Condition (ABC) is 1st order. 
  
  TO DO: Implement 2nd order accuracy for homogeneous regions.

*/


#include<stdio.h>
#include<math.h>

/* dimensions in the X, Y, and Z directions */
#define LIMX 60
#define LIMY 100
#define LIMZ 16

/* time at which source is switched off and ABC turned on */
#define SWITCH1 225 //405
#define DELAY 0

/* Total number of time steps */
#define totalT 3000

#define PI 3.14159265358979
#define MU0 1.25663706e-6
#define EPS0 8.854e-12
#define EPSR 2.2

/* globally declare fields */
double Ex[LIMX][LIMY][LIMZ], Ey[LIMX][LIMY][LIMZ], Ez[LIMX][LIMY][LIMZ];
double Hx[LIMX][LIMY][LIMZ], Hy[LIMX][LIMY][LIMZ], Hz[LIMX][LIMY][LIMZ];

/* globally declare stored field arrays for ABCs */
double HxABC1[LIMX][LIMZ], HzABC1[LIMX][LIMZ], HyABC2[LIMY][LIMZ], HzABC2[LIMY][LIMZ];
double HyABC3[LIMY][LIMZ], HzABC3[LIMY][LIMZ], HxABC4[LIMX][LIMZ], HzABC4[LIMX][LIMZ];
double HxABC5[LIMX][LIMY], HyABC5[LIMX][LIMY], ExABC6[LIMX][LIMZ], EzABC6[LIMX][LIMZ];
double ExABC5[LIMX][LIMZ], EzABC5[LIMX][LIMZ];

/* Storing the output to calculate S-parameters */
double EzOut[totalT];

/*  I want all variables declared globally */
int i, j, k, ntime, frame=0;


/*  Variables defining lattice and time steps, from Sheen, 1990 */
double delX, delY, delZ, delT;
double T, T0, temp;

/*  ABC Coefficients....and the FDTD coefficients */
double abcFSx, abcFSy, abcFSz, abcDIx, abcDIy, abcDIz, abcBx, abcBy, abcBz, cF, cB, cD;
double tMUX, tMUY, tMUZ, tEPX, tEPY, tEPZ, tERX, tERY, tERZ, tEBX, tEBY, tEBZ;

FILE *out;


/* declaration of functions */
void Initialize();
void UpdateEfields();
void Conductors();
void Source();
void FirstABC();
void UpdateHfields();
void SecondABC();

int main()
{

  FILE *in;
  char basename[80]="junk", filename[100];
  char outputF[20]="Incident.txt";

  out = fopen(outputF, "w");

  /* Define the Space */
  delX = 0.389e-3;
  delY = 0.400e-3;
  delZ = 0.265e-3;
  delT = 0.441e-12;


  /*  The source parameters */
  T = 15.e-12;
  T0 = 3.*T;

  /* Define Free Space ABC coefficients */
  cF = 1/sqrt(MU0*EPS0);
  abcFSx = (  delT*cF - delX  )/( delT*cF + delX );
  abcFSy = (  delT*cF - delY  )/( delT*cF + delY );
  abcFSz = (  delT*cF - delZ  )/( delT*cF + delZ );


  /* Define Dielectric ABC coefficients */
  cD = 1/sqrt(MU0*EPS0*EPSR);
  abcDIx = ( delT*cD - delX )/( delT*cD + delX );
  abcDIy = ( delT*cD - delY )/( delT*cD + delY );
  abcDIz = ( delT*cD - delZ )/( delT*cD + delZ );

  /* Define Boundary ABC coefficients */
  cB = 1/sqrt(MU0 * EPS0 * (EPSR+1.)/2. );
  abcBx = ( delT*cB - delX)/(delT*cB + delX);
  abcBy = ( delT*cB - delY)/(delT*cB + delY);
  abcBz = ( delT*cB - delZ)/(delT*cB + delZ);

  printf("abcBx = %lf, abcBy = %lf, abcBz = %lf\n", abcBx, abcBy, abcBz);

  /* Define H coefficients */
  tMUX = delT/MU0/delX;
  tMUY = delT/MU0/delY;
  tMUZ = delT/MU0/delZ;

  /* E coefficients (Free Space)*/
  tEPX = delT/EPS0/delX;
  tEPY = delT/EPS0/delY;
  tEPZ = delT/EPS0/delZ;

  /* E Coefficients (Dielectric) */
  tERX = delT/EPS0/EPSR/delX;
  tERY = delT/EPS0/EPSR/delY;
  tERZ = delT/EPS0/EPSR/delZ;

  /* E Coefficients (Boundary) */
  tEBX = delT/EPS0*2./(EPSR+1)/delX;
  tEBY = delT/EPS0*2./(EPSR+1)/delY;
  tEBZ = delT/EPS0*2./(EPSR+1)/delZ;



  /*  Zero Out the Fields */
  Initialize();


  printf("Pete Rules %lf %lf %lf\n", tEBX, tEBY, tEBZ);
  //printf("Enter the number of time steps\n");
  //scanf("%d", &totalT);



  /*Do time stepping */
  for(ntime=0; ntime<totalT; ntime++){

    printf("Doing time step %d\n", ntime);


	UpdateEfields();
	FirstABC();
	Conductors();
	Source();
	UpdateHfields();
	SecondABC();



    /* Write out E-field */
    k=2;
    if( ntime % 5 ==0){
      sprintf(filename, "%s.%d", basename, frame++);
      in=fopen(filename, "w");
      for(i=0; i<LIMX; i++)
		for(j=0; j<LIMY; j++)
		  fprintf(in, "%lf\n", Ez[i][j][k]);

	  fclose(in);
    }


  }/*End of time stepping*/

  fclose(out);
}


/*  Function:  Initialize Fields   */
/**********************************
 *  Zeros all fields and ABC storage arrays *
 *****************************************/
void Initialize(){

  /*Initializing fields to zero*/
  for(i=0; i<LIMX; i++)
    for(j=0; j<LIMY; j++)
      for(k=0; k<LIMZ; k++){
		Ex[i][j][k]=0.;
		Ey[i][j][k]=0.;
		Ez[i][j][k]=0.;
		Hx[i][j][k]=0.;
		Hy[i][j][k]=0.;
		Hz[i][j][k]=0.;
      }

    /* ABCs on wall Y = 0 and Y = LIMY */
	for(i=0; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
			HxABC1[i][k] = 0.;
			HzABC1[i][k] = 0.;
			HxABC4[i][k] = 0.;
			HzABC4[i][k] = 0.;
			ExABC6[i][k] = 0.;
			EzABC6[i][k] = 0.;
			ExABC5[i][k] = 0.;
			EzABC5[i][k] = 0.;
		}

	/* ABCs on wall X = 0 and X = LIMX */
	for(j=0; j<LIMY; j++)
		for(k=0; k<LIMZ; k++){
			HyABC2[j][k] = 0.;
			HzABC2[j][k] = 0.;
			HyABC3[j][k] = 0.;
			HzABC3[j][k] = 0.;
		}

	/* ABC on Z = LIMZ **/
	for(i=0; i<LIMX; i++)
		for(j=0; j<LIMY; j++){
			HxABC5[i][j] = 0.;
			HyABC5[i][j] = 0.;
		}
}
/*  End Initialize Function **********/




/*  Function:  UpdateEfields()    */
/**********************************
/*  Updates Ex, Ey, and Ez.
 *
 ***********************************/
void UpdateEfields(){

    /*Update Ex Field */
    for(i=0; i<LIMX; i++)
      for(j=1; j<LIMY-1; j++)
		for(k=1; k<LIMZ; k++){
		  if(k>3){
			Ex[i][j][k] += tEPY*(Hz[i][j][k] - Hz[i][j-1][k])
					  - tEPZ*(Hy[i][j][k] - Hy[i][j][k-1]);
			}else if(k==3){
				Ex[i][j][k] += tEBY*(Hz[i][j][k] - Hz[i][j-1][k])
				   - tEBZ*(Hy[i][j][k] - Hy[i][j][k-1]);
			}else{
				Ex[i][j][k] += tERY*(Hz[i][j][k] - Hz[i][j-1][k])
					   - tERZ*(Hy[i][j][k] - Hy[i][j][k-1]);
				  }
		}

    /*Special update for Ex because of PMC on wall y=0 */
	/*Simulated a PMC here.  See Sheen, 1990 for details.  */
	/*  Don't need to do this for Ez, because it is where the source is added */
    if(ntime < SWITCH1){
      j=0;
      for(i=0; i<LIMX; i++)
		for(k=1; k<LIMZ; k++){
			if(k>3){
			Ex[i][j][k] += tEPY* 2.*Hz[i][j][k]
				    - tEPZ*(Hy[i][j][k] - Hy[i][j][k-1]);
			}else if(k==3){
				Ex[i][j][k] += tEBY*2.*Hz[i][j][k]
				   - tEBZ*(Hy[i][j][k] - Hy[i][j][k-1]);
			}else{
				Ex[i][j][k] += tERY*2.*Hz[i][j][k]
					   - tERZ*(Hy[i][j][k] - Hy[i][j][k-1]);
			}
		}
    }

    /*Updating the Ey fields */
    for(i=1; i<LIMX; i++)
      for(j=0; j<LIMY-1; j++)
		for(k=1; k<LIMZ; k++){
		  if(k>3){
			Ey[i][j][k] += tEPZ * (Hx[i][j][k] - Hx[i][j][k-1])
						      - tEPX * (Hz[i][j][k] - Hz[i-1][j][k]);
			}else if(k==3){
				 Ey[i][j][k] += tEBZ*(Hx[i][j][k] - Hx[i][j][k-1])
						      - tEBX*(Hz[i][j][k] - Hz[i-1][j][k]);
			}else{
				Ey[i][j][k] += tERZ*(Hx[i][j][k] - Hx[i][j][k-1])
						      - tERX*(Hz[i][j][k] - Hz[i-1][j][k]);
			  }
		}

    /* Updating Ez fields */
    for(i=1; i<LIMX; i++)
      for(j=1; j<LIMY-1; j++)
		for(k=0; k<LIMZ; k++){
		if(k>2){
			Ez[i][j][k] += tEPX*(Hy[i][j][k] - Hy[i-1][j][k])
						      - tEPY*(Hx[i][j][k] - Hx[i][j-1][k]);
		}else if(k<3){
			Ez[i][j][k] += tERX * (Hy[i][j][k] - Hy[i-1][j][k])
						      - tERY*(Hx[i][j][k] - Hx[i][j-1][k]);
			}
		}

	/* Special Update for Ez fields on y=0 */
    if(ntime < SWITCH1){
      j=0;
      for(i=1; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
		  if(k>=3){
			Ez[i][j][k] += tEPX* (Hy[i][j][k] - Hy[i-1][j][k])
					    - tEPY*2.*Hx[i][j][k];
			}else if(k<3){
				Ez[i][j][k] += tERX * (Hy[i][j][k] - Hy[i-1][j][k])
					      - tERY*2.*Hx[i][j][k];
			  }
		}
    }
	/* Save Required E-field */
		/*now...22 is about the middle of the strip, 40 is arbitrary, 2 is just under the strip*/
	EzOut[ntime] = Ez[22][40][2] + Ez[22][40][1] + Ez[22][40][0];
	//fprintf(out, "%lf\n", EzOut[ntime]);
}
/* End UpdateEfields function ********************************/




/*  Function:  Conductors ()   */
/******************************
/*  Zeros the tangential (Ex, Ey) fields on the conductor
    surfaces (ground plane, microstrip, antenna)
/**************************************/
void Conductors(){

   /*Ground plane at z=0********************/
    k=0;
    for(i=0; i<LIMX; i++)
      for(j=0; j<LIMY; j++){
		Ex[i][j][k] = 0.;
		Ey[i][j][k] = 0.;
      }
	/*************************************/


    /* When you ONLY have the microstrip transmission line, you have this uncommented and
	   comment out the other uStrip and patch antenna sections *
    k=3;
    for(i=19; i<25; i++)
      for(j=0; j<LIMY; j++){
		Ex[i][j][k] = 0.;
		Ey[i][j][k] = 0.;
      }


    k=3;
    i=25;
    for(j=0; j<LIMY; j++){
      Ey[i][j][k] = 0.;
    }
    /* *******I add the above nodes to be zero to make the strip symmetric**********   */


    /* uStrip - zeroing tangential fields  **************************************/
    k=3;
    for(i=19; i<25; i++)
      for(j=0; j<50; j++){
		Ex[i][j][k] = 0.;
		Ey[i][j][k] = 0.;
      }


    k=3;
    i=25;
    for(j=0; j<50; j++){
      Ey[i][j][k] = 0.;
    }
	/**%%Above nodes zerod to make the strip symmetric, or fields look odd **/


	/*Patch antenna ********************************/
	k=3;
	for(i=14; i<46; i++)
		for(j=50; j<89; j++){
			Ex[i][j][k] = 0.;
			Ey[i][j][k] = 0.;
		}

	i=46;
	for(j=50; j<89; j++){
		Ey[i][j][k] = 0.;
	}

	//questionable here
	j=89;
	for(i=14; i<=46; i++){
		Ex[i][j][k] = 0.;
	}
	/*********************************************/
}

/* End function:  Conductors **************************/




/* Function:  Source ********************************/
/*
/*  Adds in the source *******************************/
void Source(){

    /* Source */
    //if(ntime < SWITCH1){
      j=0;
      for(i=19; i<=25; i++)
		  for(k=0; k<3; k++){
			  temp = (ntime*delT - T0)/T;
			  //Ez[i][j][k] = exp( -temp*temp );
			  //double coef = ntime/40 - 1;
			  //Ez[i][j][k] = 10*(1-2*PI*PI*coef*coef)*exp(-PI*PI*coef*coef);
			  Ez[i][j][k] = (1 - exp( - ntime/1000. ) )*cos( 2*PI* 6.32e9 * delT*ntime);
		  }

		  /*
	  for(i=19; i<25; i++)
		  for(k=0; k<3; k++){
			  Ex[i][j][k] = 0.;
		  }*/

    //}
    //fprintf(out, "%lf\n", exp( - (  (ntime*delT - T0 )/T ) * ( (ntime*delT-T0)/T ) ) );
    fprintf(out, "%lf\n", (1 - exp( - ntime/1000. ) )*cos( 2*PI* 6.32e9 * delT*ntime) );

}

/* End Function:   Source **********************************/



/* Function:  FirstABC() **********************************/
/* ************************************************       */
/* This first ABC is the only one applied to the E-fields.*/
/* Implementation details are in Scheen, 1990.  Performed */
/* after the source is turned off.  Also stores fields    */
/* needed for next round.                                 */
void FirstABC(){
	/* ABC on the wall y=0 */
	if(ntime >= SWITCH1 + DELAY){
		j=0;
		for(i=0; i<LIMX; i++)
			for(k=0; k<LIMZ; k++){
				if(k>3){
					Ex[i][j][k] = ExABC6[i][k] + abcFSy * (Ex[i][j+1][k] - Ex[i][j][k]);
					Ez[i][j][k] = EzABC6[i][k] + abcFSy * (Ez[i][j+1][k] - Ez[i][j][k]);
				}else if(k==3){
					Ex[i][j][k] = ExABC6[i][k] + abcBy * (Ex[i][j+1][k] - Ex[i][j][k]);
					Ez[i][j][k] = EzABC6[i][k] + abcFSy * (Ez[i][j+1][k] - Ez[i][j][k]);
				}else{
					Ex[i][j][k] = ExABC6[i][k] + abcDIy * (Ex[i][j+1][k] - Ex[i][j][k]);
					Ez[i][j][k] = EzABC6[i][k] + abcDIy * (Ez[i][j+1][k] - Ez[i][j][k]);
				}
			}
	}

	/*Store Fields for this ABC*/
	j=0;
	for(i=0; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
			ExABC6[i][k] = Ex[i][j+1][k];
			EzABC6[i][k] = Ez[i][j+1][k];
		}

	/* ABC on Wall Y = LIMY */
	j=LIMY-1;
	for(i=0; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
			if(k>3){
				Ex[i][j][k] = ExABC5[i][k] + abcFSy * (Ex[i][j-1][k] - Ex[i][j][k]);
				Ez[i][j][k] = EzABC5[i][k] + abcFSy * (Ez[i][j-1][k] - Ez[i][j][k]);
			}else if(k==3){
				Ex[i][j][k] = ExABC5[i][k] + abcBy * (Ex[i][j-1][k] - Ex[i][j][k]);
				Ez[i][j][k] = EzABC5[i][k] + abcFSy * (Ez[i][j-1][k] - Ez[i][j][k]);
			}else{
				Ex[i][j][k] = ExABC5[i][k] + abcDIy * (Ex[i][j-1][k] - Ex[i][j][k]);
				Ez[i][j][k] = EzABC5[i][k] + abcDIy * (Ez[i][j-1][k] - Ez[i][j][k]);
			}
		}

	/* Save Fields for this ABC */
	j=LIMY-1;
	for(i=0; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
			ExABC5[i][k] = Ex[i][j-1][k];
			EzABC5[i][k] = Ez[i][j-1][k];
		}

}
/* End Function:   FirstABC *******************************/




/* Function:  UpdateHfields() *****************************/
/* Updates H-fields.   Nothing special here.  *************/
void UpdateHfields(){

    /*Update Magnetic Fields */
    for(i=0; i<LIMX; i++)
      for(j=0; j<LIMY-1; j++)
		for(k=0; k<LIMZ-1; k++)
			Hx[i][j][k] += tMUZ*( Ey[i][j][k+1] - Ey[i][j][k] )
							- tMUY*(Ez[i][j+1][k] - Ez[i][j][k]);


    for(i=0; i<LIMX-1; i++)
      for(j=0; j<LIMY; j++)
		for(k=0; k<LIMZ-1; k++)
			Hy[i][j][k] += tMUX*(Ez[i+1][j][k] - Ez[i][j][k])
							- tMUZ*(Ex[i][j][k+1] - Ex[i][j][k]) ;

    for(i=0; i<LIMX-1; i++)
      for(j=0; j<LIMY-1; j++)
		for(k=0; k<LIMZ; k++)
		  Hz[i][j][k] += tMUY*(Ex[i][j+1][k] - Ex[i][j][k])
						      - tMUX*(Ey[i+1][j][k] - Ey[i][j][k]);

}
/* End Function:   UpdateHfields() ***********************/


/* Function:  SecondABC() *********************************/
/* Implements the remaining ABCs on the walls X = 0, LIMX */
/* and Y = LIMY, Z = LIMZ.   Also, the required fields are*/
/* then stored.                                           */
void SecondABC(){
	/*Implementation of ABCs */

	/* ABC 1 - on the wall y=LIMY *
	j=LIMY-1;
	for(i=0; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
			if(k>3){
				Hx[i][j][k] = HxABC1[i][k] + abcFSy * (Hx[i][j-1][k] - Hx[i][j][k]);
				Hz[i][j][k] = HzABC1[i][k] + abcFSy * (Hz[i][j-1][k] - Hz[i][j][k]);
			}else if(k==3){
				Hx[i][j][k] = HxABC1[i][k] + abcFSy * (Hx[i][j-1][k] - Hx[i][j][k]);
				Hz[i][j][k] = HzABC1[i][k] + abcBy  * (Hz[i][j-1][k] - Hz[i][j][k]);
			}else{
				Hx[i][j][k] = HxABC1[i][k] + abcDIy * (Hx[i][j-1][k] - Hx[i][j][k]);
				Hz[i][j][k] = HzABC1[i][k] + abcDIy * (Hz[i][j-1][k] - Hz[i][j][k]);
			}
		}

	/* ABC 2 - on the wall x=0 */
	i=0;
	for(j=0; j<LIMY; j++)
		for(k=0; k<LIMZ; k++){
			if(k>3){
				Hy[i][j][k] = HyABC2[j][k] + abcFSx * (Hy[i+1][j][k] - Hy[i][j][k]);
				Hz[i][j][k] = HzABC2[j][k] + abcFSx * (Hz[i+1][j][k] - Hz[i][j][k]);
			}else if(k==3){
				Hy[i][j][k] = HyABC2[j][k] + abcFSx * (Hy[i+1][j][k] - Hy[i][j][k]);
				Hz[i][j][k] = HzABC2[j][k] + abcBx  * (Hz[i+1][j][k] - Hz[i][j][k]);
			}else{
				Hy[i][j][k] = HyABC2[j][k] + abcDIx * (Hy[i+1][j][k] - Hy[i][j][k]);
				Hz[i][j][k] = HzABC2[j][k] + abcDIx * (Hz[i+1][j][k] - Hz[i][j][k]);
			}
		}

	/* ABC 3 - on the wall x=LIMX */
	i=LIMX-1;
	for(j=0; j<LIMY; j++)
		for(k=0; k<LIMZ; k++){
			if(k>3){
				Hy[i][j][k] = HyABC3[j][k] + abcFSx * (Hy[i-1][j][k] - Hy[i][j][k]);
				Hz[i][j][k] = HzABC3[j][k] + abcFSx * (Hz[i-1][j][k] - Hz[i][j][k]);
			}else if(k==3){
				Hy[i][j][k] = HyABC3[j][k] + abcFSx * (Hy[i-1][j][k] - Hy[i][j][k]);
				Hz[i][j][k] = HzABC3[j][k] + abcBx  * (Hz[i-1][j][k] - Hz[i][j][k]);
			}else{
				Hy[i][j][k] = HyABC3[j][k] + abcDIx * (Hy[i-1][j][k] - Hy[i][j][k]);
				Hz[i][j][k] = HzABC3[j][k] + abcDIx * (Hz[i-1][j][k] - Hz[i][j][k]);
		}
	}

	/* ABC 4 - the switched on one: y=0
	j=0;
	if(ntime >= SWITCH1){
		for(i=0; i<LIMX; i++)
			for(k=0; k<LIMZ; k++){
			if(k>3){
				Hx[i][j][k] = HxABC4[i][k] + ABCcoef * (Hx[i][j+1][k] - Hx[i][j][k]);
				Hz[i][j][k] = HzABC4[i][k] + ABCcoef * (Hz[i][j+1][k] - Hz[i][j][k]);
			}else if(k==3){
				Hx[i][j][k] = HxABC4[i][k] + ABCcoef  * (Hx[i][j+1][k] - Hx[i][j][k]);
				Hz[i][j][k] = HzABC4[i][k] + ABCcoefB * (Hz[i][j+1][k] - Hz[i][j][k]);
			}else{
				Hx[i][j][k] = HxABC4[i][k] + ABCcoefD * (Hx[i][j+1][k] - Hx[i][j][k]);
				Hz[i][j][k] = HzABC4[i][k] + ABCcoefD * (Hz[i][j+1][k] - Hz[i][j][k]);
			}
		}
	}

	/* ABC 5 - The wall z=LIMZ */
	k=LIMZ-1;
	for(i=0; i<LIMX; i++)
		for(j=0; j<LIMY; j++){
			if(k>3){
				Hx[i][j][k] = HxABC5[i][j] + abcFSz * (Hx[i][j][k-1] - Hx[i][j][k]);
				Hy[i][j][k] = HyABC5[i][j] + abcFSz * (Hy[i][j][k-1] - Hy[i][j][k]);
			}else if(k==3){
				Hx[i][j][k] = HxABC5[i][j] + abcFSz * (Hx[i][j][k-1] - Hx[i][j][k]);
				Hy[i][j][k] = HyABC5[i][j] + abcFSz * (Hy[i][j][k-1] - Hy[i][j][k]);
			}else{
				Hx[i][j][k] = HxABC5[i][j] + abcDIz * (Hx[i][j][k-1] - Hx[i][j][k]);
				Hy[i][j][k] = HyABC5[i][j] + abcDIz * (Hy[i][j][k-1] - Hy[i][j][k]);
			}
		}

	/* Saving */
	k=LIMZ-1;
	for(i=0; i<LIMX; i++)
		for(j=0; j<LIMY; j++){
			HxABC5[i][j] = Hx[i][j][k-1];
			HyABC5[i][j] = Hy[i][j][k-1];
		}

	/* Saving more fields */
	j=0;
	for(i=0; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
			HxABC4[i][k] = Hx[i][j+1][k];
			HzABC4[i][k] = Hz[i][j+1][k];
		}

	/* Save fields */
	i=0;
	for(j=0; j<LIMY; j++)
		for(k=0; k<LIMZ; k++){
			HyABC2[j][k] = Hy[i+1][j][k];
			HzABC2[j][k] = Hz[i+1][j][k];
		}

	/* Save fields */
	i=LIMX-1;
	for(j=0; j<LIMY; j++)
		for(k=0; k<LIMZ; k++){
			HyABC3[j][k] = Hy[i-1][j][k];
			HzABC3[j][k] = Hz[i-1][j][k];
		}


	/* Save required fields */
	j=LIMY-1;
	for(i=0; i<LIMX; i++)
		for(k=0; k<LIMZ; k++){
			HxABC1[i][k] = Hx[i][j-1][k];
			HzABC1[i][k] = Hz[i][j-1][k];
		}


}
/* End Function:   SecondABC() *****************************/
