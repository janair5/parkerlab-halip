function fit_results = psychFit(deltaTowers, choices, deltaBins)

if nargin < 3; deltaBins = -12:4:12; end


%% Compute trials where the animal went right vs. evidence strength
numbins             = numel(deltaBins)-1;
numRight            = zeros(numbins,1);
numTrials           = zeros(numbins,1);
binVals             = zeros(numbins,1);

%% number right per bin
% actual x axis value is average of delta values within a bin
for iBin = 1:numbins
  if deltaBins(iBin) < 0
    idx      = deltaTowers >= deltaBins(iBin) & deltaTowers < deltaBins(iBin+1);
  else
    idx      = deltaTowers > deltaBins(iBin) & deltaTowers <= deltaBins(iBin+1);
  end
  numRight(iBin)   = sum(choices(idx) == 1);
  numTrials(iBin)  = sum(idx);
  binVals(iBin)    = mean(deltaTowers(idx));
end

%% Binomial CIs
[phat, pci]         = binointerval(numRight, numTrials, normcdf(-1));

%% Logistic function fit
sigmoid             = @(O,A,lambda,x0,x) O + A ./ (1 + exp(-(x-x0)/lambda));
sel                 = numTrials > 0;

if sum(sel) < 4
  psychometric      = [];
else
  psychometric      = fit ( binVals(sel), phat(sel), sigmoid  ...
    , 'StartPoint'      , [0 1 8 0]                           ...
    , 'Weights'         , ((pci(sel,2) - pci(sel,1))/2).^-2   ...
    , 'MaxIter'         , 100                                 ...
    );
end
delta               = linspace(deltaBins(1)-2, deltaBins(end)+2, 50);

%% Draw a line with error bars for data
fit_results.delta_data     = binVals(sel)';
fit_results.pright_data    = 100*phat(sel)';
fit_results.pright_error   = 100*pci(sel,:)';

if ~isempty(psychometric)
  fit_results.delta_fit  = delta';
  fit_results.pright_fit = psychometric(delta)*100;
  fit_results.slope = psychometric.A ./ (4*psychometric.lambda);
  fit_results.bias = psychometric.O;  
  
else
  fit_results.delta_fit  = [];
  fit_results.pright_fit = [];
  fit_results.slope = []; 
  fit_results.bias = []; 
end

