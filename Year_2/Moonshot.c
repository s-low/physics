/*
 * Moonshot.c
 *
 *  Created on: Dec 16, 2012
 *      Author: samlow
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    double k1x=0,k2x=0,k3x=0,k4x=0; //initialise all variables to zero. unnecessary in some cases but safe.
    double k1y=0,k2y=0,k3y=0,k4y=0;
    double k1vx=0,k2vx=0,k3vx=0,k4vx=0;
    double k1vy=0,k2vy=0,k3vy=0,k4vy=0;
    double x=0,y=0,vy=0,vx=0,t=0,h,d;
    int i=0;	//integer loopcount
    FILE *myfile; //pointer to file

    double const G=6.67E-11; //gravitational constant
    double const Me=5.9736E+24; //earth mass
    double const Mm=7.347E+24; //moon mass
    double const Xm=384400000; //position of moon in x (no y coordinate) 3.8E+8 m

    y=0;
    x=-7000000; //satellite starts in low orbit around earth, 7000km from CENTRE of earth. opposite to moon.
    vy=10381; // ~10380 m/s is good - figure of 8 'sweet spot'
    vx=0;

    h=1;

    myfile=fopen("MoonShot.txt","w"); //opens the file for output

    if (myfile!=NULL) //check to see if the file was actually opened
    {

        for(i=0; i<500000; i++) // number of data points
        {
            t=i*h; //time increments with step h in seconds

            k1x=vx;
            k1y=vy;
            k1vx= ((-Me*G*x)/pow((x*x)+(y*y),1.5))-((Mm*G*(x-Xm))/pow((pow((Xm-x),2)+(y*y)),1.5));
            k1vy= ((-Me*G*y)/pow((x*x)+(y*y),1.5))-((Mm*G*y)/pow(((pow((x-Xm),2))+(y*y)),1.5)); //k1's defined by eq. of motion

            k2x = vx + ((h*k1vx)/2); // sequence of k's important. interdependent.
            k2y = vy + ((h*k1vy)/2);
            k2vx = -1*((G*Me*(x + ((h*k1x)/2))) / pow(pow((x+((h*k1x)/2)),2) + pow((y+((h*k1y)/2)),2),1.5)) -1*((G*Mm*(x - Xm + ((h*k1x)/2))) / pow((pow(Xm-(x+((h*k1x)/2)),2)) + pow((y+((h*k1y)/2)),2),1.5));
            k2vy = -1*((G*Me*(y + ((h*k1y)/2))) / pow(pow((x+((h*k1x)/2)),2) + pow((y+((h*k1y)/2)),2),1.5)) -1*((G*Mm*(y + ((h*k1y)/2))) / pow((pow(Xm-(x+((h*k1x)/2)),2)) + pow((y+((h*k1y)/2)),2),1.5));

            k3x = vx + ((h*k2vx)/2);
            k3y = vy + ((h*k2vy)/2);
            k3vx = -1*((G*Me*(x + ((h*k2x)/2))) / pow(pow((x+((h*k2x)/2)),2) + pow((y+((h*k2y)/2)),2),1.5)) -1*((G*Mm*(x - Xm + ((h*k2x)/2))) / pow((pow(Xm-(x+((h*k2x)/2)),2)) + pow((y+((h*k2y)/2)),2),1.5));
            k3vy = -1*((G*Me*(y + ((h*k2y)/2))) / pow(pow((x+((h*k2x)/2)),2) + pow((y+((h*k2y)/2)),2),1.5)) -1*((G*Mm*(y + ((h*k2y)/2))) / pow((pow(Xm-(x+((h*k2x)/2)),2)) + pow((y+((h*k2y)/2)),2),1.5));

            k4x = vx + ((h*k3vx)/2);
            k4y = vy + ((h*k3vy)/2);
            k4vx = -1*((G*Me*(x + ((h*k3x)/2))) / pow(pow((x+((h*k3x)/2)),2) + pow((y+((h*k3y)/2)),2),1.5)) -1*((G*Mm*(x - Xm + ((h*k3x)/2))) / pow((pow(Xm-(x+((h*k3x)/2)),2)) + pow((y+((h*k3y)/2)),2),1.5));
            k4vy = -1*((G*Me*(y + ((h*k3y)/2))) / pow(pow((x+((h*k3x)/2)),2) + pow((y+((h*k3y)/2)),2),1.5)) -1*((G*Mm*(y + ((h*k3y)/2))) / pow((pow(Xm-(x+((h*k3x)/2)),2)) + pow((y+((h*k3y)/2)),2),1.5));

            y = y + h*((k1y/6) + (k2y/3) + (k3y/3) + (k4y/6)); //position values adjusted
            x = x + h*((k1x/6) + (k2x/3) + (k3x/3) + (k4x/6));

            vx = vx + h*((k1vx/6) + (k2vx/3) + (k3vx/3) + (k4vx/6));
            vy = vy + h*((k1vy/6) + (k2vy/3) + (k3vy/3) + (k4vy/6));

            fprintf(myfile,"%lf\t%lf\n",x,y); //print the distance of the satellite from lunar surface, when <500km

            if( pow(((pow((Xm-x),2))+((y)*(y))),.5) <(2.237E+6)) //only IF within ~500km + 1737km of moon's central coords.
            {
            	d=(pow(((pow((Xm-x),2))+pow(y,2)),0.5) - 1737000)/1000; //satellite to moon in km

            	printf("Success! When T = %lf s : Distance to lunar surface = %lf km\n",t,d);
            }

            if( pow(((pow((Xm-x),2))+((y)*(y))),0.5) <(2.237E+6)) //only IF within radius of moon
                        {
                        	d=(pow(((pow((Xm-x),2))+pow(y,2)),0.5) - 1737000)/1000; //satellite to moon in km

                        	printf("Crashed! When T = %lf s : Distance to lunar surface = %lf km\n",t,d);
                        }


        }
    }

    fclose(myfile);//close the file

    return 0;
}
