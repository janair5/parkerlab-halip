% Quick and dirty function to replace stat toolbox randsample
function Out = randomSample(Samp,nValues,varargin)
disp('**ALERT** randomSample used instead of statistics toolbox randSample(). Change this back before production!');
lenSamp = length(Samp);
if lenSamp == 1 % randomsample(Samp,nValues) returns nValues values sampled from the integers 1 to Samp.
    lenSamp = Samp;
    Samp = 1:Samp;
end
if nargin == 2 
    Out = Samp(ceil(rand(1,nValues)*lenSamp));
elseif nargin == 4
    weightMatrix = varargin{2};
    if length(weightMatrix) ~= lenSamp
        error('Error: weight matrix must equal length of sample');
    end
    Out = zeros(1,nValues);
    for i = 1:nValues
        thisRand = rand;
        cumProb = 0;
        found = 0;
        for j = 1:lenSamp
            if found == 0
                cumProb = cumProb + weightMatrix(j);
                if thisRand <= cumProb
                    Out(i) = Samp(j);
                    found = 1;
                end
            end
        end
        if found == 0
            error('Error: Weight vector must sum to 1')
        end
    end
    
else
    error('Error: Incorrect usage of randomSample');
end