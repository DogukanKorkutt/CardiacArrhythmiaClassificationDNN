function [denoisedSig,debug] = NLM_1dDarbon(signal,lambda,P,PatchHW)

if length(P)==1,  % scalar has been entered; expand into patch sample index vector
    Pvec = -P:P;
else
   Pvec = P;  % use the vector that has been input  
end
debug=[];
N = length(signal);

denoisedSig = NaN*ones(size(signal));

% to simpify, don't bother denoising edges
iStart=1+PatchHW+1;
iEnd = N-PatchHW;
denoisedSig(iStart:iEnd) = 0;

debug.iStart = iStart;
debug.iEnd = iEnd;

% initialize weight normalization
Z = zeros(size(signal));
cnt = zeros(size(signal));    

% convert lambda value to 'h', denominator, as in original Buades papers
Npatch = 2*PatchHW+1;
h = 2*Npatch*lambda^2;

for idx = Pvec  % loop over all possible differences: s-t
    % do summation over p  - Eq. 3 in Darbon
    k=1:N;
    kplus = k+idx;
    igood = find(kplus>0 & kplus<=N);  % ignore OOB data; we could also handle it
    SSD=zeros(size(k));
    SSD(igood) = (signal(k(igood))-signal(kplus(igood))).^2;
    Sdx = cumsum(SSD);
   
    for ii=iStart:iEnd  % loop over all points 's'
        distance = Sdx(ii+PatchHW) - Sdx(ii-PatchHW-1); % Eq 4; this is in place of point-by-point MSE
        % but note the -1; we want to icnlude the point ii-iPatchHW

        w = exp(-distance/h);  %Eq 2 in Darbon
        t = ii+idx;  % in the papers, this is not made explicit
        
        if t>1 && t<=N
            denoisedSig(ii) = denoisedSig(ii) + w*signal(t);
            Z(ii) = Z(ii) + w;
            cnt(ii) = cnt(ii)+1;
        end

    end
end % loop over shifts

% now apply normalization
denoisedSig = denoisedSig./(Z+eps);
denoisedSig(1:PatchHW+1) =signal(1:PatchHW+1);
denoisedSig(end-PatchHW+1:end) =signal(end-PatchHW+1:end);
debug.Z = Z;

return