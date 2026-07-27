% Shambhala2_piped_preqn.m
% Variant of Shambhala2_piped.m for use when Python has already applied quantile
% normalization (--precompute-qn-reference flag). Identical in every way except
% that the quantilenorm(EXP) call is removed, preventing double quantile normalization.
%
% Use this script instead of Shambhala2_piped.m when the input data has been
% pre-QN'd by the Python pipeline (conftest run_pipeline with precompute_qn_reference=True,
% or harmonize_parallel with skip_qn=True).

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
    % quantilenorm deliberately omitted: input was pre-QN'd by Python

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
