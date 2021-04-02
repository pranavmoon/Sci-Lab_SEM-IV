// PDE:1-d Wave equation: d^2u/dt^2 = u_tt = C^2(d^2u/dx^2)=c^2 u_xx
// Explicit Method (FTCS: Forward time centered space)
// u_tt =(u(i,j-1)-2u(i,j)+u(i,j+1))/(dt^2) at all i i.e. position x
// u_xx = (u(i+1,j)-2u(i,j)+u(i-1,j))/(dx^2); at any time j
// 2*u(i,j+1) = const*u(i-1,j)+2(1-const)*u(i,j)+const*u(i+1,j)+2*dt*g(i,j)
// const = c^2*dt^2/dx^2;  // stability condition: const <1
// j=1 to n-1  ;  i = 2 to n-1 
//5-point formula; Forward time centered space
//-------------------------------------------------------
clc;
clear;
stacksize("max");
gstacksize("max");
// Parameters defining range in space and time
L = 1;    //length of 1-d wire
T = 1.0;  // final time

// Parameters to solve the equation 
dx = 0.02;
dt = 4e-4;
x = [0:dx:L];
t = [0:dt:T];
c2 =1; //c^2
xx = length(x);
tt = length(t);
const = c2*dt^2/dx^2;  

// Initial condition of string 
for  i=1:xx
    u(i,1) = 1/2*sin(%pi*x(i)); //100*cos(%pi*x(n)/(2*L));//
    g(i,1) = 0; // u'(x,0) = 0;
end
    
// Boundry condition
for j = 1:tt
    u(1,j) = 0; 
    u(xx,j) = 0;
end

// Solution............
    for j = 1:tt-1
     for i = 2:xx-1
     u(i,j+1) = const*u(i-1,j)+2*(1-const)*u(i,j)+const*u(i+1,j)+2*dt*g(i,1); 
     u(i,j+1) = (u(i,j+1))/2;
     end    
    end
// Plot
xset('window',0)
plot(x,[u(:,1),u(:,500),u(:,1000),u(:,2000)]);//['r-.','b','g','m']);
//plot(x',u(:,3),'r-.');
//xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = 100 cos(\frac{\pi x}{2L}), T(0,t)=100^oc, T(L,0) = 0^oc.$") //cos
xtitle("$\Large Variation\ of\ amplitude\ along\ the\ string\ at\ various\ time \\u(x,0) = 1/2 sin(\pi x), u(0,t)=0, u(L,t) = 0.$","$\Large X$","$\LARGE u(x,y)$") // sin
//xlabel("x");
//ylabel("u(x,y)")
legend("t="+string(t(1))+"s","t="+string(t(500))+"s","t="+string(t(1000))+"s","t="+string(t(2000))+"s")

xset('window',2)
//[xx,yy,zz]=genfac3d(t,t,z);
//plot3d1(tt,xx,uu,theta=125,alpha=80,flag=[-2,2,4])
plot3d1(t,x,u',theta=55,alpha=35,flag=[-4,2,4])
//xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = 100 cos(\frac{\pi x}{2L}), T(0,t)=100^oc, T(L,0) = 0^oc.$") // cos
xtitle("$\Large Variation\ of\ amplitude\ along\ the\ string\ at\ various\ time \\u(x,0) = 1/2 sin(\pi x), u(0,t)=0, u(L,t) = 0.$","$\Large X$","$\LARGE u(x,y)$") 
xlabel("Time");
ylabel("x")
zlabel("$\LARGE u(x,y)$")



