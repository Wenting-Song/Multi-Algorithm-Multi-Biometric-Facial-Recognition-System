close all;clc;
%% Read Probe images
probeNames = dir('ProbeSet/*.pgm');

M=100; % no of images in gallery
Mp=200; % no of images in probe
N=50;
probe = uint8(zeros(N,N,length(probeNames)));
for i = 1:length(probeNames)
    probe(:,:,i) = imread(strcat('ProbeSet/',probeNames(i).name));
end



%% Read Gallery images

galleryNames = dir('GallerySet/*.pgm');

gallery = uint8(zeros(N,N,length(galleryNames)));

for i = 1:length(galleryNames)
    gallery(:,:,i) = imread(strcat('GallerySet/',galleryNames(i).name));
end



%% lbp

lbp_gal = zeros(size(gallery));
lbp_pro = zeros(size(probe));

for i=1:length(gallery)
   lbp_gal(:,:,i) = lbp(gallery(:,:,i));
end 
for i=1:length(probe)
   lbp_pro(:,:,i) = lbp(probe(:,:,i));
end 
lbp_pro = reshape(lbp_pro, N*N, Mp);
lbp_gal = reshape(lbp_gal, N*N, M);

%% PCA
probe = reshape(probe, N*N, Mp);
gallery = reshape(gallery, N*N, M);
[wg, wp] = PCA(probe, gallery, 100);

%% feature level fusion 

f_gal = [wg ; gallery ; lbp_gal];
f_pro = [wp ; probe ; lbp_pro];

%% creatSim
[ simMatrix3, imposter3, genuine3 ] = createSim(f_pro, f_gal); % similar Matrix for concatenated features

% normalize using min/max, min = -1, max =1
simMatrix3 = (simMatrix3+1)./2;
imposter3 = (imposter3+1)./2;
genuine3 = (genuine3+1)./2;
plotdist( imposter3, genuine3, -0.2, 1.4, 'imposter and genuine distribution for concatenated features');
%d = 2.7284
[ CMC3 ] = calCMC( simMatrix3, genuine3, 0, 'CMC for concatenated features');
[ FMR3, FNMR3 ] = calROC( imposter3, genuine3,1, 0, 'ROC for concatenated features' );
[ EER3, thr3 ] = findEER( FMR3, FNMR3 );




[ simMatrix, imposter, genuine ] = createSim(probe, gallery); % similar Matrix for gallery and probe
simMatrix = (simMatrix+1)./2;
imposter = (imposter+1)./2;
genuine = (genuine+1)./2;
plotdist( imposter, genuine, -0.2, 1.4, 'imposter and genuine distribution for correlation');
%d = 2.3339
[ CMC ] = calCMC( simMatrix, genuine, 0, 'CMC for correlation');
[ FMR, FNMR ] = calROC( imposter, genuine,1, 0, 'ROC for correlation' );
[ EER, thr ] = findEER( FMR, FNMR );





[ simMatrixPCA, imposterPCA, genuinePCA ] = createSim(wp, wg); % similar Matrix for pca
simMatrixPCA = (simMatrixPCA+1)./2;
imposterPCA = (imposterPCA+1)./2;
genuinePCA = (genuinePCA+1)./2;
plotdist( imposterPCA, genuinePCA, -1.5, 2, 'imposter and genuine distribution for PCA');
%d = 3.2313
[ CMCPCA ] = calCMC( simMatrixPCA, genuinePCA, 0, 'CMC for PCA');
[ FMRPCA, FNMRPCA ] = calROC( imposterPCA, genuinePCA,1, 0, 'ROC for PCA' );
[ EERPCA, thrPCA ] = findEER( FMRPCA, FNMRPCA );



[ simMatrixLBP, imposterLBP, genuineLBP ] = createSim(lbp_pro, lbp_gal); % similar Matrix for lbp
simMatrixLBP = (simMatrixLBP+1)./2;
imposterLBP = (imposterLBP+1)./2;
genuineLBP = (genuineLBP+1)./2;
plotdist( imposterLBP, genuineLBP, -0.5, 1.5, 'imposter and genuine distribution for LBP');
%d = 2.6863
[ CMCLBP ] = calCMC( simMatrixLBP, genuineLBP, 0, 'CMC for LBP');
[ FMRLBP, FNMRLBP ] = calROC( imposterLBP, genuineLBP,1, 0, 'ROC for LBP' );
[ EERLBP, thrLBP ] = findEER( FMRLBP, FNMRLBP );
%% Part 2 ---Z-score Normalization

% normalize Z-score

% similar Matrix for gallery and probe
mu = mean(mean(simMatrix));
sig = std2(simMatrix);
zsimMatrix = (simMatrix-mu)./sig;
zimposter = (imposter-mu)./sig;
zgenuine = (genuine-mu)./sig;

% similar Matrix for PCA
mu = mean(mean(simMatrixPCA));
sig = std2(simMatrixPCA);
zsimMatrixPCA = (simMatrixPCA-mu)./sig;
zimposterPCA = (imposterPCA-mu)./sig;
zgenuinePCA = (genuinePCA-mu)./sig;
% similar Matrix for LBP
mu = mean(mean(simMatrixLBP));
sig = std2(simMatrixLBP);
zsimMatrixLBP = (simMatrixLBP-mu)./sig;
zimposterLBP = (imposterLBP-mu)./sig;
zgenuineLBP= (genuineLBP-mu)./sig;

zsimMatrix_mul = zsimMatrixPCA .* zsimMatrixLBP .* zsimMatrix;
zimposter_mul = zimposterPCA .* zimposterLBP .* zimposter;
zgenuine_mul = zgenuinePCA .* zgenuineLBP .* zgenuine;

% normalize between zero and one for mul
minval = min(min(zsimMatrix_mul));
maxval = max(max(zsimMatrix_mul));
zsimMatrix_mul = (zsimMatrix_mul-minval)./(maxval-minval);
zimposter_mul = (zimposter_mul-minval)./(maxval-minval);
zgenuine_mul = (zgenuine_mul-minval)./(maxval-minval);

plotdist( zimposter_mul, zgenuine_mul, -0.2, 1.4, 'score distribution with Z-score normalization for product rule');
[ CMCzMUL ] = calCMC( zsimMatrix_mul, zgenuine_mul, 0, 'CMC with Z-score normalization for product rule');
[ FMRzMUL, FNMRzMUL ] = calROC( zimposter_mul, zgenuine_mul ,1, 0, 'ROC for with Z-score normalization for product rule' );
[ EERzMUL, thrzMUL ] = findEER( FMRzMUL,FNMRzMUL);


zsimMatrix_sum = zsimMatrixPCA + zsimMatrixLBP + zsimMatrix;
zimposter_sum = zimposterPCA + zimposterLBP + zimposter;
zgenuine_sum = zgenuinePCA + zgenuineLBP + zgenuine;
% normalize between zero and one for sum
minval = min(min(zsimMatrix_sum));
maxval = max(max(zsimMatrix_sum));
zsimMatrix_sum = (zsimMatrix_sum-minval)./(maxval-minval);
zimposter_sum = (zimposter_sum-minval)./(maxval-minval);
zgenuine_sum = (zgenuine_sum-minval)./(maxval-minval);


plotdist( zimposter_sum, zgenuine_sum, -0.2, 1.4, 'score distribution with Z-score normalization for sum rule');
[ CMCzSUM ] = calCMC( zsimMatrix_sum, zgenuine_sum, 0, 'CMC with Z-score normalization for sum rule');
[ FMRzSUM, FNMRzSUM ] = calROC( zimposter_sum, zgenuine_sum ,1, 0, 'ROC for with Z-score normalization for sum rule' );
[ EERzSUM, thrzSUM ] = findEER( FMRzSUM,FNMRzSUM);

zsimMatrix_max = zeros(size(zsimMatrixPCA));

for i = 1:size(zsimMatrixPCA,1)
    for j = 1:size(zsimMatrixPCA,2)
        zsimMatrix_max(i,j) = max([zsimMatrixPCA(i,j) zsimMatrixLBP(i,j) zsimMatrix(i,j)]);
    end
end

zgenuine_max = zeros(size(zgenuinePCA));
zimposter_max = zeros(size(zimposterPCA));

for i = 1:length(zgenuinePCA)
    zgenuine_max(i) = max([zgenuinePCA(i) zgenuineLBP(i) zgenuine(i)]);
end

for i = 1:length(zimposterPCA)
    zimposter_max(i) = max([zimposterPCA(i) zimposterPCA(i) zimposterPCA(i)]);
end

% normalize between zero and one for max
minval = min(min(zsimMatrix_max));
maxval = max(max(zsimMatrix_max));
zsimMatrix_max = (zsimMatrix_max-minval)./(maxval-minval);
zimposter_max = (zimposter_max-minval)./(maxval-minval);
zgenuine_max = (zgenuine_max-minval)./(maxval-minval);

plotdist( zimposter_max, zgenuine_max, -0.2, 1.4, 'score distribution with Z-score normalization for max rule');
[ CMCzMAX ] = calCMC( zsimMatrix_max, zgenuine_max, 0, 'CMC with Z-score normalization for max rule');
[ FMRzMAX, FNMRzMAX ] = calROC( zimposter_max, zgenuine_max ,1, 0, 'ROC for with Z-score normalization for max rule' );
[ EERzMAX, thrzMAX ] = findEER( FMRzMAX,FNMRzMAX);


%% min max

mm_simMatrix_mul = simMatrixPCA .* simMatrixLBP .* simMatrix;
mm_imposter_mul = imposterPCA .* imposterLBP .* imposter;
mm_genuine_mul = genuinePCA .* genuineLBP .* genuine;
% normalize between zero and one for mul
minval = min(min(mm_simMatrix_mul));
maxval = max(max(mm_simMatrix_mul));
mm_simMatrix_mul = (mm_simMatrix_mul-minval)./(maxval-minval);
mm_imposter_mul = (mm_imposter_mul-minval)./(maxval-minval);
mm_genuine_mul = (mm_genuine_mul-minval)./(maxval-minval);

plotdist( mm_imposter_mul, mm_genuine_mul, -0.2, 1.4, 'score distribution with min max normalization for product rule');
[ CMCMUL ] = calCMC( mm_simMatrix_mul, mm_genuine_mul, 0, 'CMC with min max normalization for product rule');
[ FMRMUL, FNMRMUL ] = calROC( mm_imposter_mul, mm_genuine_mul ,1, 0, 'ROC for with min max normalization for product rule' );
[ EERMUL, thrMUL ] = findEER( FMRMUL,FNMRMUL);

%%
mm_simMatrix_sum = simMatrixPCA + simMatrixLBP + simMatrix;
mm_genuine_sum = genuinePCA + genuineLBP + genuine;
mm_imposter_sum =  imposterPCA + imposterLBP + imposter;
% normalize between zero and one for sum

minval = min(min(mm_simMatrix_sum));
maxval = max(max(mm_simMatrix_sum));
mm_simMatrix_sum = (mm_simMatrix_sum-minval)./(maxval-minval);
mm_imposter_sum = (mm_imposter_sum-minval)./(maxval-minval);
mm_genuine_sum = (mm_genuine_sum-minval)./(maxval-minval);

plotdist( mm_imposter_sum, mm_genuine_sum, -0.2, 1.4, 'score distribution with min max normalization for sum rule');
[ CMCSUM ] = calCMC( mm_simMatrix_sum, mm_genuine_sum, 0, 'CMC with min max normalization for sum rule');
[ FMRSUM,FNMRSUM ] = calROC( mm_imposter_sum, mm_genuine_sum ,1, 0, 'ROC for with min max normalization for sum rule' );
[ EERSUM, thrSUM ] = findEER( FMRSUM,FNMRSUM);

%%
mm_simMatrix_max = zeros(size(simMatrixPCA));

for i = 1:size(simMatrixPCA,1)
    for j = 1:size(simMatrixPCA,2)
        mm_simMatrix_max(i,j) = max([simMatrixPCA(i,j) simMatrixLBP(i,j) simMatrix(i,j)]);
    end
end


mm_genuine_max = zeros(size(genuinePCA));
mm_imposter_max = zeros(size(imposterPCA));

for i = 1:length(genuinePCA)
    mm_genuine_max(i) = max([genuinePCA(i) genuineLBP(i) genuine(i)]);
end

for i = 1:length(imposterPCA)
    mm_imposter_max(i) = max([imposterPCA(i) imposterLBP(i) imposter(i)]);
end
% normalize between zero and one for max
minval = min(min(mm_simMatrix_max));
maxval = max(max(mm_simMatrix_max));
mm_simMatrix_max = (mm_simMatrix_max-minval)./(maxval-minval);
mm_imposter_max = (mm_imposter_max-minval)./(maxval-minval);
mm_genuine_max = (mm_genuine_max-minval)./(maxval-minval);

plotdist( mm_imposter_max, mm_genuine_max, -1, 2, 'score distribution with min max normalization for max rule');
[ CMCMAX ] = calCMC( mm_simMatrix_max, mm_genuine_max, 0, 'CMC with min max normalization for max rule');
[ FMRMAX, FNMRMAX ] = calROC( mm_imposter_max, mm_genuine_max ,1, 0, 'ROC for with min max normalization for max rule' );
[ EERMAX, thrMAX ] = findEER( FMRMAX,FNMRMAX);
%% Part 3

dec_simMatrixPCA = simMatrixPCA > thrPCA;
dec_imposterPCA = imposterPCA > thrPCA;
dec_genuinePCA = genuinePCA > thrPCA;

dec_simMatrixLBP = simMatrixLBP > thrLBP;
dec_imposterLBP = imposterLBP > thrLBP;
dec_genuineLBP = genuineLBP> thrLBP;


dec_simMatrix = simMatrix > thr;
dec_imposter = imposter > thr;
dec_genuine = genuine> thr;
%% AND
simAND = dec_simMatrixPCA & dec_simMatrixLBP & dec_simMatrix;
imposterAND = dec_imposterPCA & dec_imposterLBP & dec_imposter;
genuineAND = dec_genuinePCA & dec_genuineLBP & dec_genuine;

GARand = sum(genuineAND)./length(genuineAND);
FARand = sum(imposterAND)./length(imposterAND);
%% OR
simOR = dec_simMatrixPCA | dec_simMatrixLBP | dec_simMatrix;
imposterOR = dec_imposterPCA | dec_imposterLBP | dec_imposter;
genuineOR = dec_genuinePCA | dec_genuineLBP | dec_genuine;

GARor = sum(genuineOR)./length(genuineOR);
FARor = sum(imposterOR)./length(imposterOR);

%% majority
simm = (dec_simMatrixPCA+dec_simMatrixLBP+dec_simMatrix)>1;
imposterm = (dec_imposterPCA+dec_imposterLBP+dec_imposter)>1;
genuinem = (dec_genuinePCA+dec_genuineLBP+dec_genuine)>1;

GARm = sum(genuinem)./length(genuinem);
FARm = sum(imposterm)./length(imposterm);
