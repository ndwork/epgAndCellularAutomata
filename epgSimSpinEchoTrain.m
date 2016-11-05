function out = epgSimSpinEchoTrain( alpha0, alpha, t1, t2, T1, T2, ...
  varargin )
  % out = run_epgSimulation( alpha0, alpha, t1, t2, T1, T2, [ M, N ] )

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

  out = zeros(M,3,N);
  Q = zeros(3,N);
  Q(3,1) = 1;
  Ralpha0 = epgMakeR( alpha0, phi0 );
  Q = Ralpha0 * Q;

  out(1,:,:) = Q;
  Ralpha = epgMakeR( alpha, phi );
  for m=2:M
    Q = epgRelax( Q, t1, T1, T2 );
    Q = epgGrad( Q, 1 );
    Q = Ralpha * Q;
    Q = epgRelax( Q, t2, T1, T2 );
    Q = epgGrad( Q, 1 );
    out(m,:,:) = Q;
  end

end
