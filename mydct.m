function b=mydct(a)
n = size(a,1);
m = size(a,2);

% Pad or truncate a if necessary
aa = a(1:n,:);

if ((rem(n,2)==1) || (~isreal(a))) % odd case
  % Form intermediate even-symmetric matrix.
  y = zeros(2*n,m);
  y(1:n,:) = aa;
  y(n+1:n+n,:) = flipud(aa);

  % Perform FFT
  yy = fft(y);


  % Compute DCT coefficients
  ww = (exp(-1i*(0:n-1)*pi/(2*n))/sqrt(2*n)).';
  ww(1) = ww(1) / sqrt(2);
  b = ww(:,ones(1,m)).*yy(1:n,:);

else % even case

  % Re-order the elements of the columns of x
  y = [ aa(1:2:n,:); aa(n:-2:2,:) ];

  % Compute weights to multiply DFT coefficients
  ww = 2*exp(-1i*(0:n-1)'*pi/(2*n))/sqrt(2*n);
  ww(1) = ww(1) / sqrt(2);
  W = ww(:,ones(1,m));

  % Compute DCT using equation (5.92) in Jain
  b = W .* fft(y);

end

if isreal(a), b = real(b); end

