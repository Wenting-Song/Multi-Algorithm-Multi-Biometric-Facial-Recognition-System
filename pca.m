function [ wg, wp ] = PCA( probe, gallery, K )
%PCA Summary of this function goes here
%   Detailed explanation goes here
%   M: number of faces in gallery
%   Mp: no. of faces in probe
%   N: width/length of an image (50 for a 50*50 image)
%   K: size of output vector

%% Start PCA, find eigenfaces and eigenweights for gallery images
%unroll gallery to vectors

%gallery = reshape(gallery,N*N,M);

% find mean vector

M = size(gallery,2);
Mp = size(probe,2);
N = size(gallery,1);

galleryMean = mean(double(gallery),2);

% find the A matrix

A = double(gallery)-galleryMean*ones(1,M);

% covariance matrix

C = A'*A;

% find eigenvalues

[v,~] = eig(C);

% find u, M best eigenvectors of AA'

e = eig(C);
u = zeros(N, M);

for i = 1:M
    u(:,i) = A*v(:,i);
    u(:,i) = u(:,i) ./ norm(u(:,i),2);
end

[~,e_indices] = sort(e,'descend');


% take top K eigenfaces

uK = zeros(N,K);

for i=1:K
    uK(:,i) = u(:,e_indices(i));
end

% find weights of training data

wg = zeros (K,M);

for i = 1:M
    wg(:,i) = uK'*A(:,i);
end

%% find weights of probe

%probe = reshape(probe,N*N,Mp);

wp = zeros (K,Mp);

Ap = double(probe)- galleryMean*ones(1,Mp);

for i = 1:Mp
    wp(:,i) = uK'*Ap(:,i);
end


end
