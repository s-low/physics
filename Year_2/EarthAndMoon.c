/*
 * EarthAndMoon.c
 *
 *  Created on: Dec 20, 2012
 *      Author: samlow
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double accx(double M,double x1,double x2,double y1,double y2); //acceleration in x and y dimensions written to subroutine
double accy(double M,double x1,double x2,double y1,double y2);

int main()
{
    double k1xe=0,k2xe=0,k3xe=0,k4xe=0; //all variables initialised to 0 as a precaution
    double k1xm=0,k2xm=0,k3xm=0,k4xm=0;
    double k1ye=0,k2ye=0,k3ye=0,k4ye=0;
    double k1ym=0,k2ym=0,k3ym=0,k4ym=0;
    double k1vxe=0,k2vxe=0,k3vxe=0,k4vxe=0;
    double k1vxm=0,k2vxm=0,k3vxm=0,k4vxm=0;
    double k1vye=0,k2vye=0,k3vye=0,k4vye=0;
    double k1vym=0,k2vym=0,k3vym=0,k4vym=0;
    double t=0,h;
    double xe,ye,vxe,vye;
    double xm,ym,vxm,vym;
    int i=0; 	//integer loop counter
    FILE *myfile;//pointer to file

    double const Me=5.97219E+24; //earth mass
    double const Mm=7.34767E+22; //moon mass

    xe=0;
    ye=0;
    vye=-12.58; //-12.58 from simple mv1 = mv2 momentum calculation
    vxe=0;

    xm=384400000; //384400000 earth-moon separation
    ym=0;
    vym=1023; // condition for stable lunar orbit
    vxm=0;

    h=10;
    printf("\n\nEARTH AND MOON\n");
    myfile=fopen("EarthAndMoon.txt","w"); //opens the file for output
    fprintf(myfile,"Xe\tYe\tXm\tYm\n");

    if (myfile!=NULL) //check to see if the file was actually opened
    {

        for(i=0; i<1000000; i++) // number of data points / time simulation runs for
        {
            t=i*h;// can be printed to the file - but not necessary for x vs y plots

            k1xe=vxe;//must calculate 4 k1's for earth and moon. 32 k's total.
            k1ye=vye;
            k1vxe= accx(Mm,xe,xm,ye,ym);//call acceleration subroutine
            k1vye= accy(Mm,xe,xm,ye,ym);

            k1xm=vxm;
            k1ym=vym;
            k1vxm= accx(Me,xm,xe,ym,ye);//swap m's and e's for moon expression
            k1vym= accy(Me,xm,xe,ym,ye);

            k2xe = vxe + ((h*k1vxe)/2);//all k's calculated sequentially earth,moon,earth,moon etc because of interdependence
            k2ye = vye + ((h*k1vye)/2);
            k2vxe = accx(Mm,xe+((h*k1xe)/2),xm+((h*k1xm)/2),ye+((h*k1ye)/2),ym+((h*k1ym)/2));
            k2vye = accy(Mm,xe+((h*k1xe)/2),xm+((h*k1xm)/2),ye+((h*k1ye)/2),ym+((h*k1ym)/2));//acceleration function called for adjusted values

            k2xm = vxm + ((h*k1vxm)/2);
            k2ym = vym + ((h*k1vym)/2);
            k2vxm = accx(Me,xm+((h*k1xm)/2),xe+((h*k1xe)/2),ym+((h*k1ym)/2),ye+((h*k1ye)/2));
            k2vym = accy(Me,xm+((h*k1xm)/2),xe+((h*k1xe)/2),ym+((h*k1ym)/2),ye+((h*k1ye)/2));

            k3xe = vxe + ((h*k2vxe)/2);
            k3ye = vye + ((h*k2vye)/2);
            k3vxe = accx(Mm,xe+((h*k2xe)/2),xm+((h*k2xm)/2),ye+((h*k2ye)/2),ym+((h*k2ym)/2));
            k3vye = accy(Mm,xe+((h*k2xe)/2),xm+((h*k2xm)/2),ye+((h*k2ye)/2),ym+((h*k2ym)/2));

            k3xm = vxm + ((h*k2vxm)/2);
            k3ym = vym + ((h*k2vym)/2);
            k3vxm = accx(Me,xm+((h*k2xm)/2),xe+((h*k2xe)/2),ym+((h*k2ym)/2),ye+((h*k2ye)/2));
            k3vym = accy(Me,xm+((h*k2xm)/2),xe+((h*k2xe)/2),ym+((h*k2ym)/2),ye+((h*k2ye)/2));

            k4xe = vxe + ((h*k3vxe)/2);
            k4ye = vye + ((h*k3vye)/2);
            k4vxe = accx(Mm,xe+(h*k3xe),xm+(h*k3xm),ye+(h*k3ye),ym+(h*k3ym));
            k4vye = accy(Mm,xe+(h*k3xe),xm+(h*k3xm),ye+(h*k3ye),ym+(h*k3ym));

            k4xm = vxm + ((h*k3vxm)/2);
            k4ym = vym + ((h*k3vym)/2);
            k4vxm = accx(Me,xm+(h*k3xm),xe+(h*k3xe),ym+(h*k3ym),ye+(h*k3ye));
            k4vym = accy(Me,xm+(h*k3xm),xe+(h*k3xe),ym+(h*k3ym),ye+(h*k3ye));//similar to k2 statement but no factor 1/2

            ye = ye + h*((k1ye/6) + (k2ye/3) + (k3ye/3) + (k4ye/6)); //final values below
            xe = xe + h*((k1xe/6) + (k2xe/3) + (k3xe/3) + (k4xe/6));

            vxe = vxe + h*((k1vxe/6) + (k2vxe/3) + (k3vxe/3) + (k4vxe/6));
            vye = vye + h*((k1vye/6) + (k2vye/3) + (k3vye/3) + (k4vye/6));

            ym = ym + h*((k1ym/6) + (k2ym/3) + (k3ym/3) + (k4ym/6));
            xm = xm + h*((k1xm/6) + (k2xm/3) + (k3xm/3) + (k4xm/6));

            vxm = vxm + h*((k1vxm/6) + (k2vxm/3) + (k3vxm/3) + (k4vxm/6));
            vym = vym + h*((k1vym/6) + (k2vym/3) + (k3vym/3) + (k4vym/6));

            fprintf(myfile,"%lf\t%lf\t%lf\t%lf\n",xe,ye,xm,ym);

        }
    }

    fclose(myfile);//close the file
    printf("\nFile written\n");

    return 0;
}




double accx(double M,double x1,double x2,double y1,double y2) //function to return acceleration in x direction for either body at any coord.
{
	double r,a;
	double const G=6.67E-11;//grav.constant

	r=pow((pow((x1-x2),2) + pow((y1-y2),2)),1.5);
	a=-(((G*M)*(x1-x2)))/r;
	return a;	//MUST RETURN SOMETHING as function is double type
}

double accy(double M,double x1,double x2,double y1,double y2) //function to return acceleration in y direction for either body at any coord.
{
	double r,a;
	double const G=6.67E-11;

	r=pow((pow((x1-x2),2) + pow((y1-y2),2)),1.5);
	a=-(((G*M)*(y1-y2)))/r;
	return a; //MUST RETURN as above
}






