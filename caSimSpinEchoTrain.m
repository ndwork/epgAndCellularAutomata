
function out = caSimSpinEchoTrain( alpha0, alpha, t1, t2, T1, T2, ...
  varargin )
  % out = run_caSimulation( alpha0, alpha, t1, t2, T1, T2, [ M, N ] )
  % out is an array of size MxNx3
  %   out(., 1,2,3, .) is F+,F-,Z
  % alpha0 = the first flip angle
  % alpha = all following flip angles
  % M - (optional) the maximum time step of interest
  % N - (optional) the maximum coefficient of interest

  defaultM = 30;
  defaultN = 50;
  p = inputParser;
  p.addOptional('M', defaultM, @isnumeric );
  p.addOptional('N', defaultN, @isnumeric );
  p.parse(varargin{:});
  M = p.Results.M;
  N = p.Results.N;
  phi0 = 0;
  phi = 0;

  Ralpha0 = epgMakeR( alpha0, phi0 );
  Q0 = Ralpha0 * [0;0;1];
  CQ0 = caProcess( Q0, t1, t2, T1, T2, M, N, alpha );

  M0 = 1;
  Qtilde = zeros(3,N);
  Qtilde(1,2) = exp(-t2/T2)*(1-exp(-t1/T1))*M0*(-1i*exp(1i*phi))*sin(alpha);
  Qtilde(3,1) = exp(-t2/T1)*cos(alpha)*(1-exp(-t1/T1))*M0 + (1-exp(-t2/T1))*M0;
  CQtilde = caProcess( Qtilde, t1, t2, T1, T2, M-1, N, alpha );

  out = CQ0;
  for i=2:M
    out(i:M,:,:) = out(i:M,:,:) + CQtilde(1:M-i+1,:,:);
  end
  
end


function out = caProcess( Q0, t1, t2, T1, T2, M, N, alpha )

  out = zeros( M, 3, N );
  if numel(Q0) == 3
    out(1,:,1) = Q0;
  else
    out(1,:,:) = Q0;
  end

  Ralpha = epgMakeR( alpha );
  
  subQ = zeros( 3, 5 );
  for m=2:M
    lastQ = squeeze( out(m-1,:,:) );

    for n=1:N
      subQ(:) = 0;

      if n==1
        subQ(:,3:end) = lastQ(:,1:3);
        subQ(1,1) = conj( lastQ(2,3) );
        subQ(1,2) = conj( lastQ(2,2) );
        subQ(2,1) = conj( lastQ(1,3) );
        subQ(2,2) = conj( lastQ(1,2) );
        subQ(3,1) = conj( lastQ(3,3) );
        subQ(3,2) = conj( lastQ(3,2) );

      elseif n==2
        subQ(:,2:end) = lastQ(:,1:4);
        subQ(1,1) = conj( lastQ(2,2) );
        subQ(2,1) = conj( lastQ(1,2) );
        subQ(3,1) = conj( lastQ(3,2) );

      elseif n==N-1
        subQ(:,1:4) = lastQ(:,n-2:N);

      elseif n==N
        subQ(:,1:3) = lastQ(:,n-2:N);

      else
        subQ(:,:) = lastQ(:,n-2:n+2);

      end

      out(m,:,n) = caRule( subQ, Ralpha, t1, t2, T1, T2 );
    end

  end

end


% function out = caRule( subQ, Ralpha, t1, t2, T1, T2 )
% 
%   out = zeros(3,1);
%   tmp = zeros(3,1);
% 
%   E11 = exp(-t1/T1);
%   E12 = exp(-t1/T2);
%   E21 = exp(-t2/T1);
%   E22 = exp(-t2/T2);
% 
%   n = 3;
%   tmp(1) = subQ(1,n-1);
%   tmp(2) = subQ(2,n-1);
%   tmp(3) = subQ(3,n-1);
%   out(1) = Ralpha(1,:) * tmp;
%   
%   tmp(1) = subQ(1,n+1);
%   tmp(2) = subQ(2,n+1);
%   tmp(3) = subQ(3,n+1);
%   out(2) = Ralpha(2,:) * tmp;
% 
%   tmp(1) = subQ(1,n);
%   tmp(2) = subQ(2,n);
%   tmp(3) = subQ(3,n);
%   out(3) = Ralpha(3,:) * tmp;
%   
% %   n = 3;
% %   tmp(1) = E12 * subQ(1,n-2);
% %   tmp(2) = E12 * subQ(2,n);
% %   tmp(3) = E11 * subQ(3,n-1);
% %   out(1) = E22 * Ralpha(1,:) * tmp;
% % 
% %   tmp(1) = E12 * subQ(1,n);
% %   tmp(2) = E12 * subQ(2,n+2);
% %   tmp(3) = E11 * subQ(3,n+1);
% %   out(2) = E22 * Ralpha(2,:) * tmp;
% % 
% %   tmp(1) = E12 * subQ(1,n-1);
% %   tmp(2) = E12 * subQ(2,n+1);
% %   tmp(3) = E11 * subQ(3,n);
% %   out(3) = E21 * Ralpha(3,:) * tmp;
% 
% end

function out = caRule( subQ, Ralpha, t1, t2, T1, T2 )

  out = zeros(3,1);
  tmp = zeros(3,1);

  E11 = exp(-t1/T1);
  E12 = exp(-t1/T2);
  E21 = exp(-t2/T1);
  E22 = exp(-t2/T2);

  n = 3;
  tmp(1) = E12 * subQ(1,n-2);
  tmp(2) = E12 * subQ(2,n);
  tmp(3) = E11 * subQ(3,n-1);
  out(1) = E22 * Ralpha(1,:) * tmp;

  tmp(1) = E12 * subQ(1,n);
  tmp(2) = E12 * subQ(2,n+2);
  tmp(3) = E11 * subQ(3,n+1);
  out(2) = E22 * Ralpha(2,:) * tmp;

  tmp(1) = E12 * subQ(1,n-1);
  tmp(2) = E12 * subQ(2,n+1);
  tmp(3) = E11 * subQ(3,n);
  out(3) = E21 * Ralpha(3,:) * tmp;

end

