function [Summary,postDensity,ACF] = MCMCSummary(MCData,nThin)
%[Summary,Density,ACF] = MCMCSummary(MCData,nThin)
%
%This function is a generic implementation for summarizing Markov chains
%that follows the struct notation introduced in 'DSection', i.e., MCData is
%valid if an element MCData(s).x(a,b,c,...) is the s'th value for parameter
%x(a,b,c,...), etc. This function also computes a posterior density 
%estimate for each parameter, and implementation for computing 
%auto-covariance functions (OPTIONAL) for each parameter in MCData exists. 
%
%INPUTS:
%MCData         An S elements long struct list, i.e., the Markov chain,
%               containing all parameter samples. Possibly obtained from
%               the function 'DSection' (see 'help DSection' for more 
%               details).
%nThin          Amount of thinning applied when computing the MCMC
%               estimates. SUGGESTED USE: 'nThin = 1', i.e., every sample
%               in the chain is included.
%
%OUTPUTS:
%Summary        A struct containing all parameters estimates. e.g.
%               Summary.x(t,i,c) is the estimated expression of cell 
%               type t, probe/gene/etc. i, under experimental condition c.
%postDensity    A struct containing x- and y-coordinates for the parameter
%               posterior density. e.g. postDensity.x.xcoord is the
%               x-coordinate of posterior density for parameter x. NOTE:
%               THIS IS NOT IMPLEMENTED YET, the function will only ouput 0
%               at the moment.
%
%OPTIONAL OUTPUT:
%ACF            Auto-covariance functions for each parameter, computed
%               empirically from the Markoc chain. NOTE: can be very slow!
%               SUGGESTED USE: do not compute ACF, if possible.
%
%PROJECT WEB PAGE:
%http://www.cs.tut.fi/~erkkila2/software/dsection/index.html
%
%WRITTEN BY:
%Timo Erkkilä, 18.12.2009

postDensity = 0;

%Apply thinning to the Markov chain
MCData = MCData(:,1:nThin:end);
MCData = MCData(:);

%Number of samples in the chain after thinning
n = numel(MCData);

%Collect field names (yes, this is generic)
Names = fieldnames(MCData(1));
nNames = numel(Names);

%Go through each field name and calculate the MCMC estimate
for nameIdx = 1:nNames
    
    %Pick one name and print out a notification
    name = Names{nameIdx};
    fprintf('Now summarizing "%s"\n',name);
    Summary.(name) = 0;
    
    %If ACF must be computed, let's do it then
    if nargout == 3
        stackDimension = numel(size(MCData(1).(name)))+1;
        D = zeros(size(MCData(1).(name)));
        for i = 1:n
            D = cat(stackDimension,D,MCData(i).(name));
            
            %Calculate the MCMC estimate iteratively
            Summary.(name) = Summary.(name) + MCData(i).(name)/n;
        end
        
        [ACF.(name),Summary.(name)] = computeACF(D,stackDimension);
        
    else
        
        for i = 1:n
            
            %Calculate the MCMC estimate iteratively
            Summary.(name) = Summary.(name) + MCData(i).(name)/n;
        end
    end
end

end %[RG] FOR COMPATIBILITY WITH OCTAVE




%Function for calculating the ACF
function [ACF,mu] = computeACF(D,stackDimension)

mu = mean(D,stackDimension);

n = size(D,stackDimension);

hMax = min([100,n-1]);

ACF = [];

for h = 0:hMax
    
    ACFh = 0;
    
   for i = 1:(n-h) 
    
        ACFh = ACFh + 1/n*(getFromStack(D,stackDimension,i)-mu).*(getFromStack(D,stackDimension,i+h)-mu);

   end
    
   ACF = cat(stackDimension,ACF,ACFh);
    
end

end %[RG] FOR COMPATIBILITY WITH OCTAVE

%Function for getting elements from a stack
function Ds = getFromStack(D,stackDimension,idx)

str = repmat(':,',[1,stackDimension-1]);

eval(['Ds = D(',str,'idx);']);

end %[RG] FOR COMPATIBILITY WITH OCTAVE
