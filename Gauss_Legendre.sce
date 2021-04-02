//Gauus -Legendre Quadrature method of integration.
// Pranav Kumar
clc;
clear
function d =f(t)
 //   g = input("enter function")
    x = (b-a)/2 .* t + (b+a)/2; // change of limit [a,b] to [-1,1]----
    d= 1 ./(1+x.^2);
    d = (b-a)/2 .*d;
endfunction


function GL()
//---------Finding roots/nodes and corresponding weights-----
s=poly(0,"s"); // Define variable s.
n = input("Enter order (n>0) of Legendre polynomial : ");
P = (s^2-1)^n; // the factor 1/(2^n n!) will go away when we find the root by equatin polynomial to zero.
for i=1:n
    P = derivat(P);
end
disp(P, " The "+string(n)+" order Legendre Polynomial:");
r = roots(P); // roots of nth order Legendre polynomial
disp([r], "***-->Roots of the "+string(n)+" order Legendre Polynomial are")
r=r';

A = [];
B=[];
//Generating matrice to solve for weight: A*W=B.....
//for k=0:2*n-1
//    R = (r.^k)
//    A = [A;R];
//    if modulo(k,2)==0
//        b = 2/(k+1);
//        B=[B;b];
//    else
//        b=0;
//        B=[B;b];
//        
//    end
//end

//or by integration on f(x) from -1 to 1
funcprot(0);
function y= Int(x,k), y = x^k, endfunction
for k =0:2*n-1
    R = (r.^k)
    A = [A;R];
    B = [B;intg(-1,1,Int,1e-4)]
end


disp([A,B], "**--> Augmanted Matrix to solve for weights ")
//
// A*W = B
w = (1/A)*B;
disp([r', w],["** Root/Node", "    Weight"])

//---------Integration--------------
a = input("** Enter lower limit: ");
b = input("** Enter upper limit: ");
I = f(r)*w;
disp(I, "-->Integration value from "+string(n)+"-point Gauss-Legendre Quadrature:")
plot(r',w,"o")
xtitle("Plot of nodes and weights for Gauss-Legendre Quad", "Nodes (x)", "Weights (w)")
endfunction

//------------
disp("  For n-point Gauss Legendre Quadrature  ")

fl=1;
while fl>0
    disp([" 1. Enter 0 to exit program" ; " 2. Enter any none zero value to find Integration. ";])
    fl = input(" ");
    if fl ==0
        break;
    end    
    GL();           
end
disp(" Program stops ")
