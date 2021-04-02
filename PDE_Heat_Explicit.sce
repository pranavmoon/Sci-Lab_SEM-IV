// PDE: Heat equation: du/dt = D(d^2u/dx^2)=D u''
// Explicit Method (FTCS: Forward time centered space)
// du/dt = (u(n,j+1)-u(n,j))/dt;  at all n i.e. position x
// u'' = (u(n+1,j)-2u(n,j)+u(n-1,j))/(dx^2); at any time j
// u(n,j+1) = u(n,j) + const*(u(n-1,j) -2*u(n,j)+u(n+1,j))
// const = D*dt/dx^2;  // stability condition: const <1
// 4-point formula; Forward time centered space
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
D = 1/4;
const = D*dt/dx^2;  // stability condition: const <1

// Initial temperature of wire
for  n=1:length(x)
    u(n,1) = sin(%pi*x(n)); //100*cos(%pi*x(n)/(2*L));//
end
    
// Boundry condition
for j = 1:length(t)-1
    u(1,j) = 0; //100;   //0;
    u(length(x),j) = 0;
end

// Solution through explicit method
if const<1 then
    for j = 1:length(t)-1
        for n = 2:length(x)-1
        u(n,j+1) = u(n,j) + const*(u(n-1,j) -2*u(n,j)+u(n+1,j)); 
        end    
    end
// Plot
xset('window',0)
plot(x,[u(:,1),u(:,200),u(:,600),u(:,1000)]);//['r-.','b','g','m']);
//plot(x',u(:,3),'r-.');
//xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = 100 cos(\frac{\pi x}{2L}), T(0,t)=100^oc, T(L,0) = 0^oc.$") //cos
xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = sin(\pi x), T(0,t)=0^oc, T(L,0) = 0^oc.$","$X$","$Temperature$") // sin
xlabel("x");
ylabel("Temperature")
legend("t="+string(t(1))+"s","t="+string(t(200))+"s","t="+string(t(600))+"s","t="+string(t(1000))+"s")

xset('window',2)
//[xx,yy,zz]=genfac3d(t,t,z);
//plot3d1(tt,xx,uu,theta=125,alpha=80,flag=[-2,2,4])
plot3d1(t,x,u',theta=55,alpha=35,flag=[-4,2,4])
//xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = 100 cos(\frac{\pi x}{2L}), T(0,t)=100^oc, T(L,0) = 0^oc.$") // cos
xtitle("$\LARGE Variation\ of\ Temperature\ along\ the\ wire\ at\ various\ time \\T(x,0) = sin(\pi x), T(0,t)=0^oc, T(L,0) = 0^oc.$","$\LARGE Time$","$\LARGE X$","$\LARGE Temperature$")  //sin
//xlabel("Time");
//ylabel("x")
//zlabel("Temperature")
else
    disp(const, "Stability condition is not valid as error =")
end



