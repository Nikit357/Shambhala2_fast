
function dataN = CuBlock(data,N,k)
%CuBlock: Cross-platform normalization method based on dividing the gene
%expression microarray into blocks, approximating a cubical polynomial
%to each of them and transforming them accordingly.
%
%   DESCRIPTION: Normalization of a microarray using CuBlock, a
%   cross-platform normalization method for gene expression microarray.
%   The micro array cannot contain NaN values.
%   The micro array raw intensities must be log2 transformed and at the
%   probes level.
%   Reference: Valentin Junet, Judith Farrés, José M Mas, Xavier Daura,
%   CuBlock: a cross-platform normalization method for gene-expression microarrays,
%   Bioinformatics, 2021;, btab105, https://doi.org/10.1093/bioinformatics/btab105
%
%   INPUTS:
%
%       - data -> log2 transform of microarray raw intensities.
%       Rows are probes and columns are samples.
%       - N (OPTIONAL) -> the number of times the algorithm is
%       repeated. By default it is 30.
%       - k (OPTIONAL) -> the number of probe clusters for the
%       application of k-means to find probe-cluster partitions.
%       By default it is 5.
%
%   OUTPUTS:
%
%       - dataN -> Normalized microarray.
%
%   EXAMPLES:
%
%       -  dataN = CuBlock(data)
%       -  dataN = CuBlock(data,20,5)
%       -  dataN = CuBlock(data,[],5)
%
%default arguments
if nargin<3 || isempty(k)
    k = 5;
end
if nargin<2 || isempty(N)
    N = 30;
end
%initialize
data = double(data);
[nbProbes,nbSamples] = size(data);
dataN = zeros(nbProbes,1);
count = zeros(nbProbes,1);
%beginning of the algorithm
for nRep = 1:N
    indProbes = kmeans(data,k,'maxiter',1000);                                          %find the probes cluster partitions 
    for j=1:1                                                                           % only the input sample; P columns calibrate k-means only
        for i=1:k
            if sum(indProbes==i)>100                                                    %don't transform if the block is too small
                dataCurr = data(indProbes==i,j);
                
                % --- Primitive-based STD and MEAN ---
                valid_mask = ~isnan(dataCurr);
                valid_data = dataCurr(valid_mask);
                n_valid = sum(valid_mask);
                
                if n_valid > 1
                    curr_mean = sum(valid_data) / n_valid;
                    % Sample standard deviation using primitive sqrt and sum
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
                
                if dataCurrStd>0                                                        %don't transform if the block has a unique value
                    dataCurr = (dataCurr - curr_mean)/dataCurrStd;       %Z-transform                    
                    
                    [dataCurrS,indS] = sort(dataCurr);                                  %sort the Z-transformed block and keep the indices for sorting back
                    %The GetTargetValues algorithm (implemented in the main
                    %function); obtain the values that the polynomial in
                    %the next step should best fit
                    p = 3:2:21;
                    tol = 1e-1;
                    X = (linspace(-1,1,numel(dataCurr))').^p;
                    [~,indStdUp] = min(abs(dataCurrS-1));
                    [~,indStdDown] = min(abs(dataCurrS+1));

                    % Primitive-based column mean
                    X_subset = abs(X(indStdDown:indStdUp,:));
                    S = sum(X_subset, 1) / size(X_subset, 1);

                    
                    indP = min([numel(p),find(S<tol,1)]);
                    
                    % --- Primitive-based cubic polyfit ---
                    % Creates a Vandermonde matrix and uses core LAPACK left-division
                    V = [dataCurrS.^3, dataCurrS.^2, dataCurrS, ones(size(dataCurrS))];
                    pol = V \ X(:,indP);
                    % -------------------------------------
             
                    %The ModPol Algorithm; modify the decreasing part of
                    %the polynomial
                    currDataN = ModPol(dataCurr,indS,pol);
                    %Store the values
                    dataN(indProbes==i,j) = dataN(indProbes==i,j) + currDataN;
                    count(indProbes==i,j) = count(indProbes==i,j) + 1;
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
