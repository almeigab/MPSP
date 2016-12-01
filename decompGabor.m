% decompGabor
N = 32;
NIter = N;

dicGamma = getDicParmDiscrete('std',N);

[row ,col] = size(dicGamma); 
numGamma = row;

dicGaborComplex = zeros(numGamma,N);

for k = 1:numGamma
    dicGaborComplex(k,:) = genGaborComplex( dicGamma(k,:),N);
end

x1 = genGaborReal([32 0 0.589],pi/8,N); 
x2 = genGaborReal([16 16 2*0.589],pi/8,N); 
x3 = genGaborReal([8 2 3*0.589],pi/8,N); 

x = 3*x1-2*x2+x3;
%x = randn(1,N);
residue=x;

innerProd = zeros(1,numGamma);
opt_phi = zeros(1,numGamma);
chosenParm = zeros(NIter,5);
%figure(1);
%figure(2);
for kIter = 1: NIter
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


%result_decomp = cell(1,Nblock)
%result_decomp{kblock} = chosenParm;
%save 'result_decomp.mat' result_decomp




