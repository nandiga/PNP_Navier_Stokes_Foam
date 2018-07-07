#include "fvCFD.H"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FILE *points,*rad ;

int j=0;
int n_tot = 22957; // number of mesh points
double x_new[n_tot],y_new[n_tot],z_new[n_tot];

points = fopen("points", "r");
rad = fopen("r_variation", "a");
for(j=0;j<n_tot;j++)
{
fscanf(points, "%lf %lf %lf\n", &x_new[j],&y_new[j],&z_new[j]);
cout << x_new[j] << '\n' << endl;
if (x_new[j] == 500 && z_new[j] <=0)
fprintf(rad, "%c%lf %lf %lf%c \n",'(',x_new[j],y_new[j],z_new[j],')');

}

for(j=0;j<n_tot;j++)
{

if (x_new[j] == 10000.325 && z_new[j] <=0)
fprintf(rad, "%c%lf %lf %lf%c \n",'(',x_new[j],y_new[j],z_new[j],')');
}
for(j=0;j<n_tot;j++)
{
if (x_new[j] == 19500.65 && z_new[j] <=0)
fprintf(rad, "%c%lf %lf %lf%c \n",'(',x_new[j],y_new[j],z_new[j],')');

}

fclose(points);
fclose(rad);

    return(0);
}


// ************************************************************************* //
