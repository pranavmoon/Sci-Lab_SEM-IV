//Least square fit for x^2+random error through lsq and Numerical solution.
//Pranav Kumar
clear;
//clc
clf
n=15;//input('Enter order of polynomial fit = ')

x=[1,2,3,4,5,6,7,8,9,10,11,12,13,14]';
//y=[1.5,4.6,8.7,15.9,25.1,35.8,49.6]';

rand("seed",10);
y = x.^2 + rand(x,"n"); //normal distribution random generator.
B=[y];

//------------with scilab function: lsq------------
yf= zeros(length(x),n);
er=zeros(n);
z=ones(x);

function T = coef(i)
    for i1=1:i
        z = [z, x.^i1];
    end
    T = lsq(z,B);
endfunction

for i =1:n
    c = coef(i)
    for j=1:length(c)
        yf(:,i) = yf(:,i)+c(j,:)*x.^(j-1);
        er(i)=sum((yf(:,i)-y).^2);
    end
xset("window",1)
plot2d(x,y,-i)
plot2d(x,yf(:,i),i)
end

[m,k]=min(er);

//-----Numerical---

function [X] = matsol(n)

//A=zeros(n+1,n+1);
for i=0:n
     for  j=i:n
        A(i+1,j+1)= sum(x.^(i+j));
        A(j+1,i+1)= A(i+1,j+1);
     end
     b(i+1) = sum(x.^i .*y); 
end

X=A\b;
endfunction

yfN= zeros(length(x),n);
erN=zeros(n);
for i =1:n
    c = matsol(i)
    for j=1:length(c)
        yfN(:,i) = yfN(:,i)+c(j,:)*x.^(j-1);
        erN(i)=sum((yfN(:,i)-y).^2);
    end
end
[q,p]=min(erN);

//+++++++++++++++plotting+++++++++++++++
xset("window",2)
//plot2d(x,[y,yf(:,(1:k)),yfN(:,(1:p))],[-k,(1:k),(1:p)+1])
plot2d(x,[y,yf(:,k),yfN(:,p)],[-2,k,-p+1])
legend('Data','Poly Order ' + string(k)+..
        ', error= '+ string(m),' Num_Poly order '+ string(p)+..
        ', error= '+ string(q),2)
xtitle('$Polynomial\ fit \ to\ data.\ y=x^2$','$\Large x \rightarrow$',..
        '$\Large y \rightarrow $');
        
xdel(0);  // close graphic window of specified number.
