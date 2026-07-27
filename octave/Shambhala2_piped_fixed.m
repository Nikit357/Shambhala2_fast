% Shambhala2_piped_fixed.m
% Variant of Shambhala2_piped.m that uses CuBlock_fixed instead of CuBlock.
% Requires FIXED_CLUSTERS (G×1 integer vector, labels 1..k) to be set in the
% --eval preamble before sourcing this script.
%
% Additional --eval variable: FIXED_CLUSTERS (G×1 int32 vector of cluster labels)
%
% All other behaviour is identical to Shambhala2_piped.m.

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

    dataN = CuBlock_fixed(real(EXP),[],k,FIXED_CLUSTERS);

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
