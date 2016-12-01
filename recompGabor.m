%recompGabor
clc;
clear all;
disp('Início');
load('result_decomp.mat');
N = 32;
L = ceil(length(result_decomp)/N);
xcell = cell(L,1);
y = zeros(L*N,1);

for kBlock = 1:L
    
    disp(['Bloco: ' num2str(kBlock)]);
    structBook = result_decomp{1,kBlock};
    a = N*(kBlock-1)+1;
    b = a+N-1;
    
    %figure(1);
    %figure(2);
    x=zeros(1,N);
    for n = 1:length(structBook)
        alpha = structBook(n,1);
        s = structBook(n,2);
        u = structBook(n,3);
        xi = structBook(n,4);
        phi = structBook(n,5);
        
        %a = cos(xi*n + phi);
        %b = exp(-pi*((n-u)/s)^2);
        %c = 2^(1/4)*a*b;
        
        g = genGaborReal([s u xi],phi,N);
        
        %g = c;%/sqrt(c^2);
        %x1(n,:) = alpha*g;
        x = x + alpha*g;
        %figure(1);
        %hold on;
        %plot(g,'b');
        %plot(x,'r^:');
        %hold off;
        %pause
        
    end
    
    %x2 = sum(x1);
    %if kBlock == 1
    %    y(1,a:b) = xcell{1,kBlock};
    %else
    %a1 = a-32; b1 = b-32;
    y(a:b) = x;
    %figure(2);
    %plot(y);
    %pause
    %end
end
save 'result_recomp.mat' y;