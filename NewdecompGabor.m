% NewDecompGabor
tic;
clc;
clear all;
disp('Início');
N = 128;
NIter = N;
[z,Fs] = audioread('piano_midi_nota2.wav');
L = length(z)/N;
y = vertcat (z, zeros((ceil(L)-L)*N,1));
result_decomp = cell(1,length(y));
dicGamma = getDicParmDiscrete('std',N);
[row ,col] = size(dicGamma); 
numGamma = row;
dicGaborComplex = zeros(numGamma,N);

x = y';

innerProd = zeros(1,numGamma);
opt_phi = zeros(1,numGamma);

chosenParm = zeros(NIter,5);

for k = 1:numGamma                                                                                                                                                
        dicGaborComplex(k,:) = genGaborComplex( dicGamma(k,:),N);
end

for kBlock = 1:L
    %if ismember(kBlock,0:10:L)
    
    disp(['Bloco: ' num2str(kBlock)]);
    %end
    
    
    a = N*(kBlock-1)+1;
    b = a+N-1;
    
    residue=x(1,a:b);
    
    if norm(residue) == 0
        continue;
    end
    %figure(1);
    %figure(2);
    for kIter = 1: NIter
        %disp(['Iter: ' num2str(kIter)]);
        for k = 1:numGamma
            Pgamma = real(dicGaborComplex(k,:));
            Qgamma = imag(dicGaborComplex(k,:));
            
            innerProd_xP = residue*Pgamma';
            innerProd_xQ = residue*Qgamma';
            innerProd_PQ = Pgamma*Qgamma';
            
            normP2 = norm(Pgamma)^2;
            normQ2 = Qgamma*Qgamma';
            xi = dicGamma(k,3);
            
            [opt_phi(k),innerProd(k)] = computeOptPhaseInnerProd(innerProd_xP,innerProd_xQ,innerProd_PQ,normQ2,normP2,xi);
        end
        
        innerProdMax = max(abs(innerProd));
        indAtomMax = find(abs(innerProd)==innerProdMax);

        chosenAlpha = innerProd(indAtomMax);
        chosenOptPhi = opt_phi(indAtomMax);
        chosenGamma = dicGamma(indAtomMax,:);
    
        chosenAtom = genGaborReal(chosenGamma,chosenOptPhi,N);
        chosenParm(kIter,:) = [chosenAlpha chosenGamma chosenOptPhi];
        %figure(1);
        %plot(residue,'b');
        %hold on;
        %plot(chosenAlpha*chosenAtom,'r^:');
        %hold off;
        %pause
        %figure(2),plot(residue - chosenAlpha*chosenAtom,'m--');
        residue = residue - chosenAlpha*chosenAtom;
    end
    result_decomp{1,kBlock} = chosenParm;
 
end
save 'result_decomp.mat' result_decomp;

T = toc;
%result_decomp = cell(1,Nblock)
%result_decomp{kblock} = chosenParm;
%save 'result_decomp.mat' result_decomp