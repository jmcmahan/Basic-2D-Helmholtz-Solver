function R = dos_quad(n, tn)
% Compute the quadrature weights for the integral with the logarithmic
% term as a weight derived from integrating trigonmetric interpolation
% polynomials and taken from Colton & Kress. 
%
% n - the interval (0, 2*pi) is divided into 2*n-1 intervals. This is
%     the interval along which the integration variable is discretized.
% tn - a column vector of points to evaluate the integral at. This is
%      the variable which is independent of the integration variable.

    s = size(tn);
    if s(2) > s(1)
        tn = tn';
    end
    tj = (0:(2*n-1))*pi/n;
    m = (1:(n-1))';
    R = zeros(length(tn), length(tj));
    for i = 1:length(tn)
        ti = tn(i);
        % This is just a faster way of writing the
        % sum over m. If the index is (i,j,m), M looks like 
        % 
        % M = [(i,1,1)    (i,2,1)  ...  (i,2*n,1) ;
        %      (i,1,2)    (i,2,2)  ...  (i,2*n,2) ;
        %        ...        ...    ...     ...    ;
        %      (i,1,n-1) (i,2,n-1) ... (i,2*n,n-1);]


        M = cos(m*(ti - tj)) ./ repmat(m, size(tj))
        R(i,:) = -2*pi/n*sum(M, 1) - pi/n^2 * cos(n*(ti - tj));
    end
    
end
