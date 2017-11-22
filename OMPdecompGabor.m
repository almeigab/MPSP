% OMPdecompGabor
tic;
clc;
clear all;
disp('Início');
N = 32;
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

for kBlock = 1:ceil(L)
    oldParm = zeros(N,5);
    
    disp(['Bloco: ' num2str(kBlock)]);
    
    ini = N*(kBlock-1)+1;
    fim = ini+N-1;
    
    residue=x(1,ini:fim);
    norm(residue)
    if norm(residue) == 0
        continue;
    end
    figure(1);
    v = zeros(0,1); b = zeros(0,1); Ainv = []; a = zeros(1,0); oldParm = []; aux = [];
    
    for kIter = 1: NIter
        disp(['Iter: ' num2str(kIter)]);
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

        chosenAlpha = innerProd(indAtomMax(1));
        chosenOptPhi = opt_phi(indAtomMax(1));
        chosenGamma = dicGamma(indAtomMax(1),:);
        chosenAtom = genGaborReal(chosenGamma,chosenOptPhi,N);
        chosenParm(kIter,:) = [chosenAlpha chosenGamma chosenOptPhi];
        
        oldParm(kIter,:)=[chosenGamma,chosenOptPhi,N];
        
        for i = 1:kIter-1
            v(i,1) = chosenAtom*genGaborReal(oldParm(i,1:3),oldParm(i,4),oldParm(i,5))';
        end
        
        for i=1:kIter
            aux(i,:) = genGaborReal(oldParm(i,1:3),oldParm(i,4),oldParm(i,5));
        end
        b = Ainv*v;
        beta = 1/(1-v'*b);
        Ainv = [Ainv + beta*b*b.', [-beta*b]; [-beta*b.', beta]];

        if kIter == 1
            u = chosenAtom;
        else
            u = chosenAtom - sum(b'*aux(1:kIter-1,:),1);
        end
        
        a(kIter) = residue*chosenAtom'/norm(u)^2;
        
        a(1:kIter-1) = a(1:kIter-1) - a(kIter)*b';
        
        plot(x(1,ini:fim),'b');
        
        residue = x(1,ini:fim) - sum(a*aux,1);
        hold on
        plot(sum(a*aux,1),'r^:');
        plot(residue,'k--');
        legend('signal', 'atoms', 'residue');
        hold off
        SNR(kIter) = 20*log10(norm(x(1,ini:fim))/norm(residue));
        disp(['SNR(db): ' num2str(SNR(kIter))]);
        pause
        if SNR(kIter) > 70, break, end
    end
    result_decomp{1,kBlock} = chosenParm;

 
end
save 'result_OMPdecomp.mat' result_decomp;

T = toc;
%result_decomp = cell(1,Nblock)
%result_decomp{kblock} = chosenParm;
%save 'result_decomp.mat' result_decomp
