function result = tau(t,n,W0,m,M)

p=1-(1-t/M)^(n-1);
X0 = (-1*M/2)*floor(W0/M)^2 + (W0-M/2)*floor(W0/M);
for k = 1:m
%     W(k) = W0*2^k;
    W(k) = (W0+1)*2^k-1;   % correct to this statement from above.('20.07.24)
    X(k) = (-1*M/2)*floor(W(k)/M)^2 + (W(k)-M/2)*floor(W(k)/M);
end
% W1 = W0*2^1;
% W2 = W0*2^2;
% W3 = W0*2^3;
% W4 = W0*2^4;
% W5 = W0*2^5;
% W6 = W0*2^6;
% 
% X0 = (-1*M/2)*floor(W0/M)^2 + (W0-M/2)*floor(W0/M);
% X1 = (-1*M/2)*floor(W1/M)^2 + (W1-M/2)*floor(W1/M);
% X2 = (-1*M/2)*floor(W2/M)^2 + (W2-M/2)*floor(W2/M); 
% X3 = (-1*M/2)*floor(W3/M)^2 + (W3-M/2)*floor(W3/M);
% X4 = (-1*M/2)*floor(W4/M)^2 + (W4-M/2)*floor(W4/M);
% X5 = (-1*M/2)*floor(W5/M)^2 + (W5-M/2)*floor(W5/M); 
% % X6 = (-1*M/2)*floor(W6/M)^2 + (W6-M/2)*floor(W6/M);

switch m
    
    case 0
    
    case 1 % 1 or 0
        SUM = 0;
    case 2 % 2 or 1
        SUM = X(1)*p/2;  
    case 3 % 3 or 2
        SUM = (X(2)*p^2)/4 + (X(1)*p)/2;
    case 4 % 4 or 3
        SUM = (X(3)*p^3)/8 + (X(2)*p^2)/4 + (X(1)*p)/2;
    case 5 % 5 or 4
        SUM = (X(4)*p^4)/16 + (X(3)*p^3)/8 + (X(2)*p^2)/4 + (X(1)*p)/2;
    case 6 % 6 or 5
        SUM = (X(5)*p^5)/32 + (X(4)*p^4)/16 + (X(3)*p^3)/8 + (X(2)*p^2)/4 + (X(1)*p)/2;
    
    otherwise
        error('this backoff stage value is not supported.');
end        

    switch m
        case 0
%             result = t - (W0 +1)/( W0 + 1 + X0*(1-p) );
%             result = t - (W0 +1)/( (W0 + 1 + X0)*(1-p) );
%             result = t - (W0 +1)/( X0*(1-p) + X0 + W0 + 1);
            result = t - (W0 +1)/(W0 + 1 + X0*(1-p) + X0);
            result = t - (W0 +1)/(W0 + 1 + X0);
        otherwise
            result = t - (W0 +1)/( W0 + 1 + (1-p)*X0 + (1-p) * SUM + X(m)*(p/2)^m);
    end

% result = tau - 2*(1-2*p) ./ ( (1-2*p)*(W/Nra+1)+(p*W/Nra).*(1-(2*p).^m));
