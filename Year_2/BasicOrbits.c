#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    double k1x=0,k2x=0,k3x=0,k4x=0;
    double k1y=0,k2y=0,k3y=0,k4y=0;
    double k1vx=0,k2vx=0,k3vx=0,k4vx=0;
    double k1vy=0,k2vy=0,k3vy=0,k4vy=0;
    double x=0,y=0,vx=0,vy=0,t=0,r,h,M;
    int i,max;
    double denominator=0;
    FILE *myfile;

    double const G=6.67E-11;
    M=5.9736E+24; //earth mass
    y=0;

    x=35786000;
    vy=3074; //~35km out in the x axis, moving up in y at ~3km a second (GEOSTATIONARY ORBIT)

    h=100;

    printf("Enter the speed in the y direction and the number of data points.\n\n Recommended values are 3074 m/s and 20000 points\n\n Vy = ");
    scanf("%lf",&vy); //attribution to velocity with &
    printf(" Data points = ");
    scanf("%d",&max);

    myfile=fopen("BasicOrbitsData.txt","w"); //opens the file for output

    if (myfile!=NULL) //check to see if the file was actually opened
    {

        for(i=0; i<max; i++) // number of data points
        {
            t=i*h; //time steps with h

            k1x=vx;
            k1y=vy;
            k1vx= (-1*G*M*x)/pow(((x*x)+(y*y)),1.5);
            k1vy= (-1*G*M*y)/pow(((x*x)+(y*y)),1.5); //k's sequential. must be ordered as they are.

            r= pow((x+((h*k1x)/2)),2) + pow((y+((h*k1y)/2)),2);
            denominator = pow(r,1.5);
            k2x = vx + ((h*k1vx)/2);
            k2y = vy + ((h*k1vy)/2);
            k2vx = -1*((G*M*(x + ((h*k1x)/2))) / denominator); //tidied with denominator term above
            k2vy = -1*((G*M*(y + ((h*k1y)/2))) / denominator);

            r= pow((x+((h*k2x)/2)),2) + pow((y+((h*k2y)/2)),2);
            denominator = pow(r,1.5);
            k3x = vx + ((h*k2vx)/2);
            k3y = vy + ((h*k2vy)/2);
            k3vx = -1*((G*M*(x + ((h*k2x)/2))) / denominator);
            k3vy = -1*((G*M*(y + ((h*k2y)/2))) / denominator);

            r= pow((x+((h*k3x)/2)),2) + pow((y+((h*k3y)/2)),2);
            denominator = pow(r,1.5);
            k4x = vx + ((h*k3vx)/2);
            k4y = vy + ((h*k3vy)/2);
            k4vx = -1*((G*M*(x + ((h*k3x)/2))) / denominator);
            k4vy = -1*((G*M*(y + ((h*k3y)/2))) / denominator);


            y = y + h*((k1y/6) + (k2y/3) + (k3y/3) + (k4y/6)); // need to calculate k's above this clearly
            x = x + h*((k1x/6) + (k2x/3) + (k3x/3) + (k4x/6));

            vx = vx + h*((k1vx/6) + (k2vx/3) + (k3vx/3) + (k4vx/6));
            vy = vy + h*((k1vy/6) + (k2vy/3) + (k3vy/3) + (k4vy/6));

            fprintf(myfile,"%lf\t%lf\n",x,y);//print the position coords to a file

        }
    }

    fclose(myfile);//close the file
    printf("\nFile written\n");
    return 0;
}
