
function run_caSimSpinEchoTrain
  close all; clear;
  addpath( genpath('.') );

  alpha0 = 35 * pi/180;
  alpha = 30 * pi/180;

  t1 = 100;
  t2 = 115;

  T1 = 1000;
  T2 = 800;


  epgOut = epgSimSpinEchoTrain( alpha0, alpha, t1, t2, T1, T2 );
  caOut = caSimSpinEchoTrain( alpha0, alpha, t1, t2, T1, T2 );

  epgOut = squeeze( abs( epgOut(:,1,:) ) );
  caOut = squeeze( abs( caOut(:,1,:) ) );

  scale = 20;
  maxImgScale = 0.2;
  figure; imshow( rot90( imresize( epgOut, scale, 'nearest') ), [0 maxImgScale] );
  title('EPG Sim');
  figure; imshow( rot90( imresize( caOut, scale, 'nearest') ), [0 maxImgScale] );
  title('CA Sim');

  absDiff = abs( caOut - epgOut );
  figure; imshow( rot90( imresize( absDiff, scale, 'nearest') ), [] );
  title(['Diff - max: ', num2str(max(absDiff(:)))]);
end


