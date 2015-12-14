#include <sstream>
#include <math.h>
#include <algorithm>
#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<iomanip>
#include <string>
#include "nr3.h"
//#include "liblbfgs/liblbfgs-1.10/lib/Makefile.in"


struct Base_interp


{
  int n,mm,jsav,cor,dj;
  const Doub *xx,*yy;
  Base_interp(VecDoub_I &x,const Doub *y, int m)
  : n(x.size()),mm(m),jsav(0),cor(0),xx(&x[0]),yy(y){
    dj = MIN(1,(int)pow((Doub)n,0.25));
  }
  
  Doub interp(Doub x){
    int j1o =cor ? hunt(x) : locate(x);
    return rawinterp(j1o,x);
  }
  
  int locate(const Doub x);
  int hunt(const Doub x);
  
  Doub virtual rawinterp (int j1o, Doub x) =0;
  
  
};


int Base_interp::locate(const Doub x)

{
  int ju,jm,j1;
  if (n<2 || mm < 2 || mm > n) throw ("locate size error");
  bool ascnd=(xx[n-1] >= xx[0]);
  j1=0;
  ju=n-1;
  while (ju-j1 > 1) {
    jm = (ju+j1) >> 1;
    if ((x >= xx[jm]) == ascnd)
      j1=jm;
    else
      ju=jm;
  }
  cor = abs(j1-jsav) > dj ? 0 :1;
  jsav = j1;
  return MAX(0,MIN(n-mm,j1-((mm-2)>>1)));
}

Int Base_interp::hunt(const Doub x)
{
  int j1=jsav,jm,ju,inc=1;
  if (n<2 || mm<2 || mm > n) throw ("hunt size error");
  bool ascnd=(xx[n-1] >= xx[0]);
  if (j1< 0|| j1>n-1) {
    j1=0;
    ju=n-1;
  }else{
    if((x >= xx[j1]) == ascnd) {
      for(;;) {
	ju = j1+inc;
	if (ju >= n-1) { ju = n-1; break;}
	else if ((x < xx[ju]) == ascnd) break;
		  else{
		    j1 =ju;
		    inc += inc;
		  }
      }
    } else {
      ju = j1;
      for (;;) {
	j1 =j1-inc;
	if (j1 <= 0) { j1 =0; break;}
	else if ((x >= xx[j1]) == ascnd) break;
		  else{
		    ju = j1;
		    inc += inc;
		  }
      }
    }
  }
  while (ju-j1 > 1 ) {
    jm = (ju+j1) >> 1;
    if ((x >= xx[jm]) == ascnd)
      j1=jm;
    else 
      ju=jm;
  }
  cor = abs(j1-jsav) < dj ? 0 :1;
  jsav =j1;
  return MAX(0,MIN(n-mm,j1-((mm-2)>>1)));
}










struct Poly_interp : Base_interp
{
  
Doub dy;
Poly_interp(VecDoub_I &xv,VecDoub_I &yv, int m)
      : Base_interp(xv,&yv[0],m),dy(0.){}
      Doub rawinterp (int j1, Doub x);
};




Doub Poly_interp::rawinterp(int j1, Doub x)


{






int i,m,ns=0;
Doub y,den,dif,dift,ho,hp,w;
const Doub *xa = &xx[j1], *ya = &yy[j1];
VecDoub c(mm),d(mm);
dif=abs(x-xa[0]);
for(i=0;i<mm;i++){
  if ((dift=abs(x-xa[i])) < dif){
    ns=i;
    dif=dift;
  }
  c[i]=ya[i];
  d[i]=ya[i];
}
y=ya[ns--];
for (m=1;m<mm;m++){
  for(i=0;i<mm-m;i++){
    ho=xa[i]-x;
    hp=xa[i+m]-x;
    w=c[i+1]-d[i];
    if((den=ho-hp) == 0.0) throw ("Poly_interp error");
      den=w/den;
    d[i]=hp*den;
    c[i]=ho*den;
  }
  y+= (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
}
 return y; 
}

struct Spline_interp : Base_interp
{
	VecDoub y2;
	
	Spline_interp(VecDoub_I &xv, VecDoub_I &yv, Doub yp1=1.e99, Doub ypn=1.e99)
	: Base_interp(xv,&yv[0],2), y2(xv.size())
	{sety2(&xv[0],&yv[0],yp1,ypn);}

	Spline_interp(VecDoub_I &xv, const Doub *yv, Doub yp1=1.e99, Doub ypn=1.e99)
	: Base_interp(xv,yv,2), y2(xv.size())
	{sety2(&xv[0],yv,yp1,ypn);}

	void sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn);
	Doub rawinterp(Int jl, Doub xv);
};



void Spline_interp::sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn)
{
	Int i,k;
	Doub p,qn,sig,un;
	Int n=y2.size();
	VecDoub u(n-1);
	if (yp1 > 0.99e99)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
		u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e99)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}
Doub Spline_interp::rawinterp(Int jl, Doub x)
{
	Int klo=jl,khi=jl+1;
	Doub y,h,b,a;
	h=xx[khi]-xx[klo];
	if (h == 0.0) throw("Bad input to routine splint");
	a=(xx[khi]-x)/h;
	b=(x-xx[klo])/h;
	y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo]
		+(b*b*b-b)*y2[khi])*(h*h)/6.0;
	return y;
}


struct Linear_interp : Base_interp
{
	Linear_interp(VecDoub_I &xv, VecDoub_I &yv)
		: Base_interp(xv,&yv[0],2)  {}
	Doub rawinterp(Int j, Doub x) {
		if (xx[j]==xx[j+1]) return yy[j];
		else return yy[j+1] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
	}
};

struct Rat_interp : Base_interp
{
	Doub dy;
	Rat_interp(VecDoub_I &xv, VecDoub_I &yv, Int m)
		: Base_interp(xv,&yv[0],m), dy(0.) {}
	Doub rawinterp(Int jl, Doub x);
};

Doub Rat_interp::rawinterp(Int jl, Doub x)
{
	const Doub TINY=1.0e-99;
	Int m,i,ns=0;
	Doub y,w,t,hh,h,dd;
	const Doub *xa = &xx[jl], *ya = &yy[jl];
	VecDoub c(mm),d(mm);
	hh=abs(x-xa[0]);
	for (i=0;i<mm;i++) {
		h=abs(x-xa[i]);
		if (h == 0.0) {
			dy=0.0;
			return ya[i];
		} else if (h < hh) {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;
	}
	y=ya[ns--];
	for (m=1;m<mm;m++) {
		for (i=0;i<mm-m;i++) {
			w=c[i+1]-d[i];
			h=xa[i+m]-x;
			t=(xa[i]-x)*d[i]/h;
			dd=t-c[i+1];
			if (dd == 0.0) throw("Error in routine ratint");
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
	}
	return y;
}












