
function dataN = CuBlock_fixed(data, N, k, fixed_clusters)
%CuBlock_fixed: CuBlock normalization with precomputed gene cluster assignments.
%
%   Like CuBlock.m but skips k-means entirely. fixed_clusters (G×1 integer
%   vector, labels 1..k) are reused in every repetition.
%   Precompute once from P using precompute_clusters_from_p() in Python.
%
%   INPUTS:
%       - data           -> log2 transform of microarray. Rows=probes, cols=samples.
%       - N (OPTIONAL)   -> repetition count. Default 30.
%       - k (OPTIONAL)   -> number of clusters (must match fixed_clusters range). Default 5.
%       - fixed_clusters -> G×1 integer vector of precomputed cluster labels 1..k.
%
%   OUTPUTS:
%       - dataN -> Normalized G×1 vector for the first (input) sample.
%
if nargin<3 || isempty(k)
    k = 5;
end
if nargin<2 || isempty(N)
    N = 30;
end
data = double(data);
[nbProbes,~] = size(data);
dataN = zeros(nbProbes,1);
count = zeros(nbProbes,1);
for nRep = 1:N
    indProbes = fixed_clusters;                                                            %use precomputed clusters; skip k-means
    for j=1:1                                                                              % only the input sample; P columns calibrate k-means only
        for i=1:k
            if sum(indProbes==i)>100
                dataCurr = data(indProbes==i,j);

                % --- Primitive-based STD and MEAN ---
                valid_mask = ~isnan(dataCurr);
                valid_data = dataCurr(valid_mask);
                n_valid = sum(valid_mask);

                if n_valid > 1
                    curr_mean = sum(valid_data) / n_valid;
                    dataCurrStd = sqrt(sum((valid_data - curr_mean).^2) / (n_valid - 1));
                else
                    if n_valid == 1
                        curr_mean = valid_data(1);
                    else
                        curr_mean = 0;
                    end
                    dataCurrStd = 0;
                end
                % ------------------------------------

                if dataCurrStd>0
                    dataCurr = (dataCurr - curr_mean)/dataCurrStd;

                    [dataCurrS,indS] = sort(dataCurr);
                    p = 3:2:21;
                    tol = 1e-1;
                    X = (linspace(-1,1,numel(dataCurr))').^p;
                    [~,indStdUp] = min(abs(dataCurrS-1));
                    [~,indStdDown] = min(abs(dataCurrS+1));

                    X_subset = abs(X(indStdDown:indStdUp,:));
                    S = sum(X_subset, 1) / size(X_subset, 1);

                    indP = min([numel(p),find(S<tol,1)]);

                    V = [dataCurrS.^3, dataCurrS.^2, dataCurrS, ones(size(dataCurrS))];
                    pol = V \ X(:,indP);

                    currDataN = ModPol(dataCurr,indS,pol);
                    dataN(indProbes==i,1) = dataN(indProbes==i,1) + currDataN;
                    count(indProbes==i,1) = count(indProbes==i,1) + 1;
                end
            end
        end
    end
end
dataN = dataN./count;
end

function dataN = ModPol(data,indS,pol)

% --- Primitive-based cubic polyval ---
x_val = data(indS);
dataNS = pol(1)*(x_val.^3) + pol(2)*(x_val.^2) + pol(3)*x_val + pol(4);
% -------------------------------------

diff=dataNS(2:end)-dataNS(1:(end-1));
changeInDirectionDown = diff<0;
indDown1 = find(changeInDirectionDown,1);
n = numel(data);
if ~isempty(indDown1)
    indDownL = find(changeInDirectionDown,1,'last')+1;
    changeInDirectionUp = diff>0;
    indUp1 = find(changeInDirectionUp,1);
    indUpL = find(changeInDirectionUp,1,'last')+1;
    if dataNS(indUp1)==dataNS(1) && dataNS(indUpL)==dataNS(n)
        M = dataNS(indDown1);
        m = dataNS(indDownL);
        if (n-indDownL+1)<=indDown1
            dataNS(indDown1:indDownL) = M;
            dataNS((indDownL+1):end) = M + dataNS((indDownL+1):end) - m;
        else
            dataNS(indDown1:indDownL) = m;
            dataNS(1:(indDown1-1)) = m + dataNS(1:(indDown1-1)) - M;
        end
    else
        M = dataNS(indUpL);
        m = dataNS(indUp1);
        dataNS(indUpL:end) = M;
        dataNS(1:indUp1) = m;
    end
end
dataN = NaN(n,1);
dataN(indS) = dataNS;
end
