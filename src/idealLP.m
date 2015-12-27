function H0 = idealLP(N,kcutoff)
% This function creates an ideal low-pass filter H0[k], where k ranges from
% 0 to N-1, given a cutoff frequency kcutoff.  The function sets H0[k]=1
% for k < kcutoff = 1 and H0[k]=0 for k>=kcutoff.  The function takes into
% account that indexing in matlab starts at 1.  So k=0, which corresponds
% to H0[0], is actually the vector element H0(1).  The functio has inputs
%       N: the length of vector representing the signal
%       kcutoff: the first index at which H0[k] = 0
% and outputs
   %    H0: a vector of length N that represents the filter H0[k]
%
% Written by Joel Voldman, 2013, 2014
        H0 = ones(N,1);     % initialize variable to all ones
        
        % Set the correct indexes to zero.  H[kcutoff] = 0.  Because matlab
        % starts indexing at 1, we set H0(kcutoff+1) = 0.  Near N, we stop
        % at N-kcutoff+1 because the last H0[k] that will equal zero is
        % kcutoff away from N, H0[N - kcutoff].  Again due to the
        % indexing in matlab, we need to add 1, ending up with
        % H0(N-kcutoff+1)
        H0(kcutoff+1:(N-kcutoff+1)) = 1.0e-6;  
end
