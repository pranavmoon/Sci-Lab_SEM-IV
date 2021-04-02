//Gauus -Hermite Quadrature method of integration.
// int_-inf to inf [a_0+a_1 x+a_2 x^2+..]dx = w_i * f(x_i)
// Pranav Kumar
clc;
clear
function d =f(x)
    //d= cos(x)+x.^4
    d= 1 ./(1+x.^2);
endfunction

function GH()
//    ---------Finding roots/nodes and corresponding weights-----
n = input("Enter order (n>0) of Hermite polynomial : ");
//------ scilab does not derivat transcendental function 
//s=poly(0,"s"); // Define variable s. only for polynomials
//
//H = exp(-s^2); // poly does not define transcendental function.
//for i=1:n
//    H = derivat(H);
//end
//
//H = (-1)^n* exp(s^2)*H;
////----Use reccurrence relation to generate polynomial-----
s=poly(0,"x");
function hr = H(n,s)
    if n==0 then
     hr = 1;
    elseif n==1 then    
     hr = 2*s;   
    else
     hr = 2*s*H(n-1,s)-2*(n-1)*H(n-2,s);
    end
endfunction
hp = H(n,s);
disp(hp,"Hermite polynomial--")
//--another way to find H- polynomials...
//h0=1;
//h1=2*s;
//function hr = H(h0,h1,n)
//    hr=2*s*h1-2*n*h0;
//endfunction
////for i = 0:n-1
//    hp = H(h0,h1,i);
//    disp(hp);
//    h0=h1;
//    h1=hp;
//end
//-----------------------------
r = roots(hp); // roots of nth order Hermite polynomial
disp([r], "***-->Roots of the "+string(n)+" order Hermite Polynomial is")
r=r';

A = [];
B=[];
//Generating matrice to solve for weight: A*W=B.....
// RHS column vector by integration of e^(-x^2) x^k for k=0,1,2...2n-1...
for k =0:2*n-1
    R = (r.^k);
    A = [A;R];
    Int = 'x^k*exp(-x^2)';
    II = integrate(Int,'x',-100,100,1e-4);
    B = [B;II];
end


disp([A,B], "**--> Augmanted Matrix to solve for weights ")
//
// A*W = B
w = (1/A)*B; //
disp(linsolve(A,-B)'," Weights : "); 

disp([r', w],["** Root/Node", "    Weight"])

//---------Integration--------------
I = f(r)*w;
disp(I, "-->Integration value from "+string(n)+"-point Gauss-Hermite Quadrature:")
plot(r',w,"o")
xtitle("Plot of nodes and weights for Gauss-Hermite Quad", "Nodes (x)", "Weights (w)")
endfunction

//------------
disp("  For n-point Gauss Hermite Quadrature  ")

fl=1;
while fl>0
    disp([" 1. Enter 0 to exit program" ; " 2. Enter any none zero value to find Integration. ";])
    fl = input(" ");
    if fl ==0
        break;
    end    
    GH();           
end
disp(" Program stops ")
