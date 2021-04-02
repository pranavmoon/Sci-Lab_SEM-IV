// PDE: Heat equation: du/dt = D(d^2u/dx^2)=D u''
// Crank-Nicolson method: another version.
// Simplify from the average of the implicit and explicit formula.
// The central difference at time step j+1/2.
// Since u values are known only at integer value of j; so 
// u(n,j+1/2) replaced by 1/2[u(n,j+1) + u(n,j)];
// -const*u(n-1,j+1)+ (2+2const)*u(n,j+1)-const* u(n+1,j+1) =
// const*u(n-1,j)+ (2-2const)*u(n,j)+const* u(n+1,j); 
// const = D dt/dx^2;
// u(:,j+1) = Inv(M_j+1)[N_j *u(:,j)+ R_j -r_j];
//---------------------------------------------------------

clc;
clear;
stacksize("max");
gstacksize("max");
// Parameters defining range in space and time
L = 1;    //length of 1-d wire
T = 1.0;  // final time

// Parameters to solve the equation 
dx = 0.05    ;
dt = 4e-4;
x = [0:dx:L];
t = [0:dt:T];
D = 1/4;
const = D*dt/dx^2;
nn=length(x);
tt=length(t);
// Initial temperature of wire
for  n=1:nn
    u(n,1) =  sin(%pi*x(n));// 100*cos(%pi*x(n)/(2*L));// 
end
    
// Boundry condition
for j = 1:tt
    u(1,j) = 0; //100;   //0;
    u(nn,j) = 0;
end

// Solution through average of implicit and explicit method
// --column vector of b.c.(R_j and r_j) to be added in u(2:nn-1,j) --
for j = 1:tt 
    R = [const*u(1,j);zeros(nn-4,1);const*u(nn,j)];
    r = (-1)*R;
end
//---- generating tri-diagonal matrix and inverse---
M1(1:nn-3) = -const;
M0(1:nn-2) = 2+2*const;
M = diag(M0,0)+diag(M1,1)+diag(M1,-1);
M= inv(M);

N1(1:nn-3) = const;
N0(1:nn-2) = 2-2*const;
N = diag(N0,0)+diag(N1,1)+diag(N1,-1);

//--finding solution at various time j+1 -------
for j = 2:tt
    uu = N * u(2:nn-1,j-1)+R-r;
    u(2:nn-1,j) = M*uu;
end

// Plot
subplot(121)
plot(x,[u(:,1),u(:,200),u(:,600),u(:,1000)]);//['r-.','b','g','m']);
//plot(x',u(:,3),'r-.');
//xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = 100 cos(\frac{\pi x}{2L}), T(0,t)=100^oc, T(L,0) = 0^oc.$")// cos
xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = sin(\pi x), T(0,t)=0^oc, T(L,0) = 0^oc.$","$X$","$Temperature$") // sin

legend("t="+string(t(1))+"s","t="+string(t(200))+"s","t="+string(t(600))+"s","t="+string(t(1000))+"s")

subplot(122)
//[xx,yy,zz]=genfac3d(t,t,z);
//plot3d1(tt,xx,uu,theta=125,alpha=80,flag=[-2,2,4])
plot3d1(t,x,u',theta=135,alpha=60,flag=[-4,2,4])
//xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = 100 cos(\frac{\pi x}{2L}), T(0,t)=100^oc, T(L,0) = 0^oc.$") // cos
xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = sin(\pi x), T(0,t)=0^oc, T(L,0) = 0^oc.$","$Time$","$X $","$Temperature$") //sin

//
