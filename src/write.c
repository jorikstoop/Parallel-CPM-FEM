// file: write.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
void write_cells(VOX* pv, int increment)
{
	int v;
    int vx,vy;
    char filename[100];
    char astring[20];
	char bstring[20];
   	FILE *ofp;

   	myitostr(increment, astring);
   	myitostr(myRank, bstring);
	strcpy(filename, "./output/ctags");
   	strcat(filename, astring); //append timestep
	strcat(filename, "_r");
   	strcat(filename, bstring); //append rank id
   	strcat(filename, ".out");

	ofp = fopen(filename,"w");

	for(vx=0; vx<NVX; vx++) {
		for (vy=2; vy<my_nvy-2; vy++) { //only loop over owned
            v = vx + vy * NVX;
     		fprintf(ofp ,"%d ", pv[v].ctag);
        }
        fprintf(ofp, "\n");
    }
   	fflush(ofp);  fclose(ofp);
}

////////////////////////////////////////////////////////////////////////////////
void write_pstrain(VOX* pv, NOD* pn, int increment)
{
	int v;
    int vx,vy;
	double estrains[3],L1,L2,v1[2],v2[2];
    char filename[100];
    char astring[20];
	char bstring[20];
   	FILE *ofp;

   	myitostr(increment, astring);
   	myitostr(myRank, bstring);
	strcpy(filename, "./output/pstrain");
   	strcat(filename, astring); //append timestep
	strcat(filename, "_r");
   	strcat(filename, bstring); //append rank id
   	strcat(filename, ".out");

	ofp = fopen(filename,"w");
	for (vy=2; vy<my_nvy-2; vy++) { //only loop over owned voxels
		for(vx=0; vx<NVX; vx++) {

		get_estrains(pn,vx,vy,estrains);
		L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);
		if(L1>L2)
		{
			fprintf(ofp,"%d ", (int)(1000000*L1));
			fprintf(ofp,"%d ", (int)(1000*v1[0]));
			fprintf(ofp,"%d ", (int)(1000*v1[1]));
			fprintf(ofp,"%d\n",(int)(1000000*L2));
		}
		else
		{
			fprintf(ofp,"%d ", (int)(1000000*L2));
			fprintf(ofp,"%d ", (int)(1000*v2[0]));
			fprintf(ofp,"%d ", (int)(1000*v2[1]));
			fprintf(ofp,"%d\n",(int)(1000000*L1));
		}
		}
	}
	fflush(ofp); fclose(ofp);
}

void write_timing(
    double time_setup,
    double time_simloop,
    double time_output,
    double time_tracforce,
    double time_fem,
    double time_pcg,
	double time_fem_comm,
    double time_cpm,
	double time_cell_comm,
    int increment)
{
    char filename[100];
    char astring[20];
	char bstring[20];
    FILE *ofp;

    // create filename like ./output/timing_t3000_s300.csv (time and domain size)
	myitostr(increment, astring);
	myitostr(NVX,bstring);
	strcpy(filename, "./output/timing_t");
    strcat(filename, astring); //add time
	strcat(filename,"_s"); 
	strcat(filename,bstring); // add size
    strcat(filename, ".csv");

    ofp = fopen(filename, "w");
    if (ofp == NULL) {
        printf("Error opening timing file!\n");
        return;
    }

    // Write header
    fprintf(ofp, "Process,Time\n");

    // Write timing data
    fprintf(ofp, "setup,%f\n", time_setup);
    fprintf(ofp, "simloop,%f\n", time_simloop);
    fprintf(ofp, "out,%f\n", time_output);
    fprintf(ofp, "tf,%f\n", time_tracforce);
    fprintf(ofp, "fem,%f\n", time_fem);
    fprintf(ofp, "pcg,%f\n", time_pcg);
	fprintf(ofp, "node_comm,%f\n", time_fem_comm);
    fprintf(ofp, "cpm,%f\n", time_cpm);
	fprintf(ofp, "cell_comm,%f\n", time_cell_comm);

    fflush(ofp);
    fclose(ofp);
}