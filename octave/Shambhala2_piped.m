% Shambhala2_piped.m
% Modified from Shambhala2.m with exactly 3 changes:
%   1. args.txt reading block removed (NH, NP, k injected via --eval preamble)
%   2. readExpressionData reads from /dev/stdin instead of P_prim.txt
%   3. Output written to stdout (fd 1) instead of Cu_bis.txt

inData=readExpressionData('/dev/stdin','log2');
Exp = inData.Samples;
SYMBOL = inData.GeneList;
SN = inData.SamplesName;
NS = length(SN);

for ( i = 1:NH )
    message = sprintf('Harmonizing sample %d out of %d',i,NH);
    disp(message);

    i0 = i;
    for ( jjj = 1:NP )
        i0 = [i0 (jjj+NH)];
    end

    EXP = Exp(:,i0);
    EXP = quantilenorm(EXP);

    dataN = CuBlock(real(EXP),[],k);

    log2e = log2(exp(1));
    DataN = dataN/log2e;

    EXPN = exp(DataN)-1;
    vecN = EXPN(:,1);

    if ( i == 1 )
        OUT = vecN;
    else
        OUT = [OUT vecN];
    end

end

nG=size(OUT,1);
nS=size(OUT,2);
for i=1:nG
    fprintf(1,'%s',SYMBOL{i,1});
    for j=1:nS
        fprintf(1,' %f',OUT(i,j));
    end
    fprintf(1, '\n');
end
