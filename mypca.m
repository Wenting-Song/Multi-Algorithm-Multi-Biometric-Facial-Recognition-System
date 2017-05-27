function [ egvalue,egvector ] = mypca( X )
% By default, each line in X is a sample
% First column in returned egvalue and egvector represent the first PC

[line,column] = size(X);

%subtract mean
for i = 1:column
    X(:,i) = X(:,i) - mean(X(:,i)); 
end

covmat = X' * X;

[egvector,egvalue] = eig(covmat);

% Reverse the matrix so that PC arranges in a order of high to low
egvalue = fliplr(egvalue);
egvector = fliplr(egvector);


end

