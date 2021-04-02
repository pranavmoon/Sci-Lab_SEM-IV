//Gauus -Laguerre Quadrature method of integration.
// int_ 0 to inf [a_0+a_1 x+a_2 x^2+..]dx = w_i * f(x_i)
// Pranav Kumar
clc;
clear
function d =f(x)
    d = x.^3 + sin(r);
    //d= 1 ./(1+x.^2);
endfunction

function LG()
//    ---------Finding roots/nodes and corresponding weights-----
    n = input("Enter order (n>0) of Laguerre polynomial : ");
s=poly(0,"x");

function Lar = Lr(n,s)
    if n==0 then
        Lar = 1;
    elseif n==1 then
        Lar = 1-s;
    else
        Lar = (1+2*(n-1)-s)*Lr(n-1,s)-(n-1)^2*Lr(n-2,s);
    end
endfunction
Lgr = Lr(n,s);
disp(Lgr, "Laguerre polynomial (recursive)")
L0=1;
L1=1-s;

function La = L(L0,L1,n)
    La = (1+2*n-s)*L1-n^2*L0;
endfunction

for i = 1:n-1
    Lg = L(L0,L1,i);
//    disp(Lg);
    L0=L1;
    L1=Lg;
end
disp(Lg," Laguerre polynomial ");
r = roots(Lg); // roots of nth order Laguerre polynomial
disp([r], "***-->Roots of the "+string(n)+" order Laguerre Polynomial are")
r=r';

A = [];
B=[];
//Generating matrice to solve for weight: A*W=B.....
// RHS column vector by integration of e^(-x^2) x^k for k=0,1,2...2n-1...
for k =0:2*n-1
    R = (r.^k);
    A = [A;R];
    Int = 'x^k*exp(-x)';
    II = integrate(Int,'x',0,100);
    B = [B;II];
end


disp([A,B], "**--> Augmanted Matrix to solve for weights ")
//
// A*W = B
w = (1/A)*B;
disp([r', w],["** Root/Node", "    Weight"])

//---------Integration--------------
I = f(r)*w;
disp(I, "-->Integration value from "+string(n)+"-point Gauss-Laguerre Quadrature:")
plot(r',w,"o")
xtitle("Plot of nodes and weights for Gauss-Laguerre Quad", "Nodes (x)", "Weights (w)")
endfunction

//------------
disp("  For n-point Gauss Laguerre Quadrature  ")

fl=1;
while fl>0
    disp([" 1. Enter 0 to exit program" ; " 2. Enter any none zero value to find Integration. ";])
    fl = input(" ");
    if fl ==0
        break;
    end    
    LG();           
end
disp(" Program stops ")
