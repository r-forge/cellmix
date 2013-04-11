function visualizeMCData(MCData,varargin)
%visualizeMCData(MCData,param1,param1Obs,param2,param2Obs,...)
%
%Visualizes Markov chains and draws their histograms.
%
%INPUTS:
%MCData         An S elements long struct list, i.e., the Markov chain,
%               containing all parameter samples. Possibly obtained from
%               the function 'DSection' (see 'help DSection' for more 
%               details).
%param1         The name (string) of the first (multi-dimensional) parameter you 
%               want to visualize; e.g., 'x', and x is n-dimensional, i.e., 
%               MCData(i).x(d_1,d_2,...,d_n). 
%param1Obs      A logical array of size 'size(MCData(i).(param1))' that has
%               1's whenever that sub-dimension is wanted to be visualized.
%               SUGGESTED USE: if the chain for 'MCData(1:S).x(t,i,c)'
%               needs to be visualized, you input pairs 'x' and 'xObs',
%               where 'xObs(t,i,c) = 1' and 0 otherwise.
%
%PROJECT WEB PAGE:
%http://www.cs.tut.fi/~erkkila2/software/dsection/index.html
%
%WRITTEN BY:
%Timo Erkkilä, 18.12.2009

MCData = MCData(:);

n = numel(MCData);

Names = fieldnames(MCData(1));

nPairs = numel(varargin(:))/2;
varargin = reshape(varargin(:),[2,nPairs])';

ObsNames = varargin(:,1);

for i = 1:nPairs
    ObsName = ObsNames{i};
    if any(strcmpi(ObsName,Names))
        I = varargin{i,2};
        nI = sum(I(:));
        D = zeros(nI,n);
        for j = 1:n
            D(:,j) = MCData(j).(ObsName)(varargin{i,2});
        end
        figure;
        subplot(2,1,1); plot((1:n)',D');
        xlabel(['Samples for parameter ''',ObsName,''''],'FontSize',14);
        grid on;
        ylabel('Sampled value','FontSize',14);
        subplot(2,1,2); hist(D');
        xlabel(['Histogram for parameter ''',ObsName,''''],'FontSize',14);
        ylabel('Histogram counts','FontSize',14);
        grid on;
    else
        error('Wrong input argument: "%s" not found in the Data struct!',ObsName);
        break;
    end
end

