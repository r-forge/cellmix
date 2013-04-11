function [MCData,x_LS,lambda_LS] = DSection(Y,p0,W0,W_proposal,Treatment,nBurnIn,nSamples,samplep,summarize,verbose)
%[MCData,x_LS,lambda_LS] = DSection(Data,p0,W0,Wp,Treatment,B,S,samplep,summarize)
%
%This function estimates cell type (and condition) specific gene expression 
%profiles from heterogeneous measurements, given prior predictions about 
%cell type proportions for each tissue sample.
%
%INPUTS:
%Data           <I-by-J> matrix of measurements from heterogeneous tissues.
%               I is the number of probes/genes/etc., and J is the number
%               of tissues.
%p0             <T-by-J> matrix of prior predictions on cell type
%               proportions. T is the number of cell types, and columns in
%               p0 must be positive and add up to one.
%W0             Prior prediction weight, i.e., degree of confidence, on p0.
%               Defines the peakedness of Dirichlet density around p0.
%               NOTE: keep W0 >= T.
%W_proposal     Transition kernel weight, defines the peakedness of
%               Dirichlet density around p*, the old value. The higher
%               W_proposal is, the smaller the proposal steps around p* are. 
%Treatment      <1-by-J> vector of treatment indices, so that 
%               unique(Treatment) = [1,2,...,C], where C is the number of
%               treatments including control, i.e., "no treatment", if
%               available.
%B              Amount of burn-in. NOTE: keep B > 0.
%S              Amount of sampling. NOTE: keep S > 0.
%samplep        binary variable, indicating whether to sample from the
%               posterior for cell type proportions (1) or not (0).
%               SUGGESTED USE: sample from the posterior (samplep = 1).
%
%OUTPUTS:
%MCData         An S elements long struct list, i.e., the Markov chain,
%               containing all parameter samples. 
%               SUGGESTED USE: 'Summary = MCMCSummary(MCData,nThin)' 
%               summarizes the chain by outputting MCMC estimates with 
%               thinning nThin (see 'help MCMCSummary' for more details).
%
%OPTIONAL OUTPUTS:
%x_LS           Least-squares solution for cell type specific expression 
%               profiles,
%lambda_LS      Corresponding precision estimates (inverses of variances) 
%               for each probe/gene/etc.
%
%PROJECT WEB PAGE:
%http://www.cs.tut.fi/~erkkila2/software/dsection/index.html
%
%WRITTEN BY: 
%Timo Erkkilä, 18.12.2009

if( !verbose )
PAGER('/dev/null', 'local');
page_screen_output(1, 'local');
page_output_immediately(1, 'local');
end

%Find unique Treatment indices
TrIdc = unique(Treatment);

%Compute the number of unique treatments
E = TrIdc(end);

if sum(abs(TrIdc-(1:E))) ~= 0
    error('uncorrect formatting in Treatment indices!');
end

%Data dimensions
[I,J] = size(Y);

%Number of cell types
T = size(p0,1);
Tvec = (1:T)';

%Compute a LS solution to the linear regression problem
x_LS = zeros(T,I,E);
for i = 1:I
    for e = 1:E
    x_LS(:,i,e) = p0(:,Treatment == e)'\Y(i,Treatment == e)';
    end
end

%Compute variance of the regression residuals
V = zeros(I,1);
for e = 1:E
    V = V + var(Y(:,Treatment == e)-x_LS(:,:,e)'*p0(:,Treatment == e),[],2)/E;
end

%Compute LS estimates for precision
lambda_LS = 1./V;

%Deduce parameter values for the Bayes prior on precision
params = gamfit(lambda_LS);
alpha0 = params(1);
beta0 = 1/params(2);

%Initialize true expressions with the LS solution
MCData.x = x_LS;

%Initialize precisions with the LS solution
MCData.lambda = lambda_LS;

%Prior precision (weight) for x_LS
nu_LS = alpha0/beta0;

%Initialize mixture proportions with the prior guesses
MCData.p = p0;

%Make sure there are no 0's and 1's, because if there are, the sampler may
%get jammed.
MCData.p(MCData.p > 1-eps) = 1-eps;
MCData.p(MCData.p < eps) = eps;

%Prepare the parameter struct for the complete Markov chains
if summarize
MCData(2:3) = MCData;
else
MCData(2:nSamples+1) = MCData;
end

%Start timer
tic
if( verbose )
  fprintf('Iteration ');
end

%Start sampling
for idx = 2:(nBurnIn+nSamples+1)
    
    %Print our iteration index
    if( verbose )
    	fprintf('%i/%i',idx-1,nBurnIn+nSamples);
    end
    
    %If it's burn-in, stay where you are on the chain
    if idx <= nBurnIn + 1
        iter = 1;
    elseif summarize && iter >= 3 
	% samples 1 and 2 have already been computed and stored in MCData(2) and MCData(3)
	iter = 3;
    else
        iter = idx-nBurnIn;
        MCData(iter) = MCData(iter-1);
    end	
    
    %If sampling from the posterior
    if samplep
        for j = 1:J
            %Sample mixing vector for each condition using M-H
            MCData(iter).p(:,j) = samplepvec(Y(:,j),MCData(iter).lambda,MCData(iter).x(:,:,Treatment(j)),MCData(iter).p(:,j),p0(:,j),W0,W_proposal,T);
        end
    end

    
    for i = 1:I
        
        %Prepare a row vector for sampling precision
        DJ = zeros(1,J);
        for e = 1:E
            DJ(Treatment == e) = MCData(iter).x(:,i,e)'*MCData(iter).p(:,Treatment == e);
        end
            
        %Sample precision
        MCData(iter).lambda(i) = gamrnd(alpha0+J/2,1/(beta0+1/2*(Y(i,:)-DJ)*(Y(i,:)-DJ)'));
        
        for e = 1:E
            for t = 1:T
                %Helper data for sampling cell type expression
                Z = MCData(iter).p(Tvec ~= t,Treatment == e)'*MCData(iter).x(Tvec ~= t,i,e);
                P = MCData(iter).lambda(i)*(Y(i,Treatment == e)*MCData(iter).p(t,Treatment == e)'-MCData(iter).p(t,Treatment == e)*Z)+x_LS(t,i,e)*nu_LS;
                Q = MCData(iter).lambda(i)*sum(MCData(iter).p(t,Treatment == e).^2)+nu_LS;
                
                %Sample cell type expression
                MCData(iter).x(t,i,e) = P/Q+1/sqrt(Q)*randn(1);

            end
            
        end
        
    end

    % if one is interested in summarized data only: sum up
    if summarize && iter >= 3
        MCData(2).x += MCData(3).x;
	MCData(2).p += MCData(3).p;
	MCData(2).lambda += MCData(3).lambda;
    end
    
    %Erase the tail of the print-out
    if( verbose )
 	b = 1 + numel(num2str(idx-1)) + numel(num2str(nBurnIn+nSamples));
	str = repmat('\b',[1,b]);
	fprintf(str);
    end
end

%Remove the first sample from the struct, it's not part of the chain
MCData = MCData(2:end);

% If summarized data if requested then only MCData(1) is meaningfull
if summarize 
MCData = MCData(1);
MCData.x = MCData.x / nSamples;
MCData.p = MCData.p / nSamples;
MCData.lambda = MCData.lambda / nSamples;
end

%Stop timer
if( verbose )
    fprintf(' done!\n');
end
t = toc;
if( verbose )
    fprintf('Time elapsed: %fs\n',t);
end

end %[RG] FOR COMPATIBILITY WITH OCTAVE

%Function for sampling mixing proportions using M-H
function p_new = samplepvec(Y,lambda,x,p_old,p0,W1,W2,T)

%Propose a new p vector from the current Dirichlet distribution
p_prop = zeros(T,1);
for t = 1:T
    p_prop(t) = gamrnd(W2*p_old(t),1);
end
p_prop = p_prop/sum(p_prop);

%Push the proposal away from the boundaries, if it touches them
p_prop(p_prop > 1-eps) = 1-eps;
p_prop(p_prop < eps) = eps;

%A fraction product of Gamma functions
gamprod_frac = 1;
for t = 1:T
   gamprod_frac = gamprod_frac*gamma(W2*p_prop(t))/gamma(W2*p_old(t));
end

%Calculate the jump factor
alpha = exp(-1/2*(lambda'*((x'*p_prop).^2-2*Y.*(x'*p_prop)))+(W2*p_old-1)'*log(p_prop)+(W1*p0-1)'*log(p_prop)+1/2*(lambda'*((x'*p_old).^2-2*Y.*(x'*p_old)))-(W2*p_prop-1)'*log(p_old)-(W1*p0-1)'*log(p_old))*gamprod_frac;

%Make a decision whether to jump to a new state or to stay at the old one
if rand(1) < min(1,alpha)
    p_new = p_prop;
else
    p_new = p_old;
end

end %[RG] FOR COMPATIBILITY WITH OCTAVE

# overwrite old gamfit to stop optim::nmsmax it from outputing lots of verbose messages!
old_gamfit = @gamfit;

# overwrite old gamfit to stop it from outputing messages!
function res = gamfit(R)

  if (nargin != 1)
    print_usage;
  endif

  avg = mean(R);

  # This can be just about any search function. I choose this because it
  # seemed to be the only one that might work in this situaition...
  a=nmsmax( @gamfit_search, 1, [1e-3 inf inf 0 0], [], avg, R );

  b=a/avg;      # **

  res=[a 1/b];
endfunction

# Helper function so we only have to minimize for one variable. Also to
# inverting the output of gamlike, incase the optimisation function wants to
# maximize rather than minimize.
function res = gamfit_search( a, avg, R )
  b=a/avg;      # **
  if( b==0 )
    res = NaN;
  else
  res = -gamlike([a 1/b], R);
  end
endfunction

