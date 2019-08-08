% dct:discrete cosine transform, we use in instead of FFT with Neumann
% boundary condition, dct uses FFT so it is still o(NlogN`) operations
% abs error is misleading unless the exact solution is nealy zero. Rel
% error gives a more reliable indication of accuracy.
function a = myidct(b)  
n = size(b,1);
m = size(b,2);

% Pad or truncate b if necessary
bb = b(1:n,:);

if rem(n,2)==1 || ~isreal(b), % odd case
    % Form intermediate even-symmetric matrix.
    ww = sqrt(2*n) * exp(1i*(0:n-1)*pi/(2*n)).';
    ww(1) = ww(1) * sqrt(2);
    W = ww(:,ones(1,m));
    yy = zeros(2*n,m);
    yy(1:n,:) = W.*bb;
    yy(n+2:n+n,:) = -1i*W(2:n,:).*flipud(bb(2:n,:));
    
    y = ifft(yy);
    
    % Extract inverse DCT
    a = y(1:n,:);
    
else % even case
    % Compute precorrection factor
    ww = sqrt(2*n) * exp(1i*pi*(0:n-1)/(2*n)).';
    ww(1) = ww(1)/sqrt(2);
    W = ww(:,ones(1,m));
    
    % Compute x tilde using equation (5.93) in Jain
    y = ifft(W.*bb);
 
    
    % Re-order elements of each column according to equations (5.93) and
    % (5.94) in Jain
    a = zeros(n,m);
    a(1:2:n,:) = y(1:n/2,:);
    a(2:2:n,:) = y(n:-1:n/2+1,:);
end

if isreal(b), a = real(a); end
