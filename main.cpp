#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222
#define DATATYPE float
#define G 9.18
/*
vec3 *global_arr_MassSpeed; 

//note--> vec3 is X,Y and mass not X, Y and Z

void Initialize_Masses_and_Speeds(int light, int med, int heavy, vec3& global_arr_MassSpeed);
*/
int main(int argc, char* argv[]){
	
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
	}

	MPI_Init(&argc,&argv);

	int p, my_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = argv[1];
	int numParticleMedium = argv[2];
	int numParticleHeavy = argv[3];
	int numSteps = argv[4];
	int subSteps = argv[5];
	double timeSubStep = argv[6];
	int width = argv[7], height = argv[8];

	unsigned char* image;


	const int total_particles = numParticlesLight + numParticleMedium + numParticleHeavy;
	vec3 local_pos = new vec3[total_particles/p]; // assuming p evenly divides n.. 
	vec3 temp_pos = new vec3[total_particles/p];

	//root node stuff goes here
	if(my_rank == 0){
		//Initialize_Masses_and_Speeds(numParticlesLight, numParticleMedium, numParticleHeavy, global_arr_MassSpeed);


		//almost done, just save the image
		saveBMP(argv[9], image, width, height);
	}
	//all other nodes do this
	else{
		source = (my_rank+1)%p;
		dest = (my_rank-1+p)%p;
	}

	free(image);

	MPI_Finalize();
	return 0;
}
/*
void Initialize_Masses_and_Speeds(int light, int med, int heavy, vec3& global_arr_MassSpeed)
{
	global_arr_MassSpeed = new vec3[light+med+heavy];

	int iter = 0;
	for(int i = 0; i < light; i++){
		arr[iter++].x = velocityLightMin+drand48()*MAX_MIN_DIFF;
		arr[iter++].y = velocityLightMin+drand48()*MAX_MIN_DIFF;
		arr[iter++].z = massLightMin+drand48()*MAX_MIN_DIFF; //mass
	}
	for(int j = 0; j < med; j++){
		arr[iter++].x = velocityMediumMin+drand48()*MAX_MIN_DIFF;
		arr[iter++].y = velocityMediumMin+drand48()*MAX_MIN_DIFF;
		arr[iter++].z = massMediumMin+drand48()*MAX_MIN_DIFF; //mass
	}
	for(int k = 0; k < heavy; k++){
		arr[iter++].x = velocityHeavyMin+drand48()*MAX_MIN_DIFF;
		arr[iter++].y = velocityHeavyMin+drand48()*MAX_MIN_DIFF;
		arr[iter++].z = massHeavyMin+drand48()*MAX_MIN_DIFF; //mass
	}
}*/

double Calculate_Force(DATATYPE x1_pos, DATATYPE x2_pos, DATATYPE y1_pos, DATATYPE y2_pos,\
DATATYPE x1_vel, DATATYPE x2_vel, DATATYPE y1DATATYPE m1, DATATYPE m2, ){



}