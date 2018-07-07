#include "fvCFD.H"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FILE *out, *out1, *rad, *U, *U_r ;

int ite,i=0, j;
float z = 0.0;
ite = 0;
int n_x_left_micro_1 = 60; // number of grid points in x direction
int n_y_left_micro_1 = 26; // number of grid points in y direction
int n_x_left_micro_2 = 60; // number of grid points in x direction
int n_y_left_micro_2 = 50; // number of grid points in y direction
int n_x_nano_1 = 100; // number of grid points in x direction
int n_y_nano_1 = 26; // number of grid points in y direction
int n_x_rht_micro_1 = 20; // number of grid points in x direction
int n_y_rht_micro_1 = 26; // number of grid points in y direction
int n_x_rht_micro_2 = 20; // number of grid points in x direction
int n_y_rht_micro_2 = 25; // number of grid points in y direction

int n_micro_end = (n_x_left_micro_1*n_y_left_micro_1)+(n_x_left_micro_2*n_y_left_micro_2);
int n_nano_end = (n_x_left_micro_1*n_y_left_micro_1)+(n_x_left_micro_2*n_y_left_micro_2)+(n_x_nano_1*n_y_nano_1);
int n_tot = (n_x_left_micro_1*n_y_left_micro_1)+(n_x_left_micro_2*n_y_left_micro_2)+(n_x_nano_1*n_y_nano_1)+(n_x_rht_micro_1*n_y_rht_micro_1)+(n_x_rht_micro_2*n_y_rht_micro_2); // total number of grid points

float x_new[n_tot],y_new[n_tot];
float U_x[n_tot], U_y[n_tot], U_z[n_tot] ;
char str[80];

out = fopen("mesh_points.txt", "w");
forAll(mesh.C(),counter)
{
ite ++;
fprintf(out,"%f  %f \n",mesh.C()[counter].x(),mesh.C()[counter].y());
 
}// end forAll
printf("number of mesh points %d \n" , ite);
fclose(out);

out = fopen("mesh_points.txt", "r");
U = fopen("U", "r");

// dumping the radial variation of individual current and concentration in these files
 
out1 = fopen("dump_xy", "w");
rad = fopen("r_variation", "w");
U_r = fopen("U_r", "w");

for(j=1;j<=n_tot;j++)
{
fscanf(out, "%f %f \n", &x_new[j],&y_new[j]);
fscanf(U, "%f %f %f \n",&U_x[j],&U_y[j],&U_z[j]);
}

// printing the radial variation of current and concentration at a particular axial distance

for(j=4561;j<=n_nano_end;j=j+n_x_nano_1) // j = 27301 is where x = 9005 nm (inside nanopore) starts; n_x is number of grid points
//for(j=27301;j<=n_micro_end;j=j+n_x_nano_1) // j = 27301 is where x = 9005 nm (inside nanopore) starts; n_x is number of grid points
{
fprintf(out1, "%f %f \n", x_new[j],y_new[j]);
fprintf(rad,"%f\n",y_new[j]);
fprintf(U_r, "%f \n", U_x[j]);
}
fclose(out);
fclose(out1);
fclose(U_r);






   Info<< "End\n" << endl;


    return(0);
}


// ************************************************************************* //
