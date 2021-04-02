// PDE: Helmholtz equation: u_xx+u_yy (= del^2u) +fu= g in 2-d
// 5-point formula; Finite difference method
// Gauss Seidel method
// u(i,j)[k+1] = 1/(4-h^2f)*(u(i+1,j)[k]+u(i,j+1)[k]
//              +u(i-1,j)[k+1]+u(i,j-1)[k+1]-h^2g)
// Lalace equation: f=0, g=0
// Poissons eq: f=0
//-------------------------------------------------------
clc;
clear;
stacksize("max");
gstacksize("max");
// Parameters defining range in space and time
L = 1;    //box of 1 sq m - 2-d (x,y)
// Parameters to solve the equation 
h = 0.01;
x = [0:h:L];
y = [0:h:L];

xn = length(x); 
yn = length(y);
// Boundry condition for x(i) and y(j)
// for x end points
    u(1,1:yn) = cosh(0.2*x(1))+cosh(0.2.*y);
    u(xn,1:yn) = cosh(0.2*x(xn))+cosh(0.2.*y);
// for y end point
    u(1:xn,1) = cosh(0.2.*x')+cosh(0.2*y(1));
    u(1:xn,yn) = cosh(0.2.*x')+cosh(0.2*y(yn));

// initial guess for gauss seidel

    u(2:xn-1,2:yn-1) = 1;
 
//--- Gauss Seidel method to solve Helmholtz equation
function [f,g] = f1(x,y)
        f = 0;   // for f(x,y)
        g = 0;  // for g(x,y)
endfunction

maxitr = input("  Enter maximum number of iteration : ");//40; // maximum iteration.
u1=u(2,2);
err = 1e-9;//input("  Enter desired accuracy  :")
disp("[  itr,    u1-u2]");
for k=1:maxitr
    u2 = u1;
    for j=2:yn-1
        for i=2:xn-1
        [f,g] = f1(x(i),y(j));
//        f = fg[1];
//        g = fg[2];
        c = 1/(4-h^2*f);
        u(i,j) =  c*(u(i+1,j)+u(i,j+1)+u(i-1,j)+u(i,j-1)-h^2*g);
        u1 = u(2,2);
        end
    end
    disp([k, abs(u1-u2)]);//,"[itr, u1-u2]");
    if abs(u1-u2)<err
        break;
    end
end

plot3d1(x,y,u',theta=55,alpha=35,flag=[-4,2,4])
xtitle("$\LARGE Laplace\ equation; \ Boundary\ condition: u(x,y) = cosh(x/5)+cosh(y/5)\\ Pranav Kumar$","$\LARGE X$","$\LARGE y$","$\LARGE u(x,y)$") 
//xlabel("X");
//ylabel("Y");
//zlabel("u(x,y");
