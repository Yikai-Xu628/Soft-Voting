%% genprojs and cl
%Orientation Determination of Cryo-EM Projection Images Using Reliable Common Lines and Spherical Embeddings
clear all;
clc;
%load('EMD33026SNR32.mat', 'corrstack','clstack');
        load('EMD33026SNR32.mat')
p=45;
method = 'sigmoid';
sigma = 0.1;
a =10;

        K = size(clstack, 2);
        n = size(projs,1);
        est_rots5=ComputeAnglesFromProjs_SE_wvote_nonlinear(clstack, corrstack, p, method, sigma, a);



        est_rots5 = align_rots(est_rots5, inv_rots_ref);
        est_ang5 = rot2euler(est_rots5, K);
        E_rots5 = est_rots5 - inv_rots_ref;
        RMSE5 = zeros(1,4);
        RMSE5(1) = sqrt(sum(E_rots5(:).^2)/numel(E_rots5));
        RMSE5(2) = sqrt(mean((est_ang5(:,1) - ang_ref(:,1)).^2));
        RMSE5(3) = sqrt(mean((est_ang5(:,2) - ang_ref(:,2)).^2));
        RMSE5(4) = sqrt(mean((est_ang5(:,3) - ang_ref(:,3)).^2));
        fprintf('RMSE_rot5: %f  RMSE_ang5:  %f %f %f\n', RMSE5);

        vol5 = recon3d_firm(projs,est_rots5,[],1e-6,100,zeros(n,n,n));
        vol5 = real(vol5);
        WriteMRC(vol5, 1.08, 'EMD33026SNR32sigmoidnew.mrc');

        
        
        
        





