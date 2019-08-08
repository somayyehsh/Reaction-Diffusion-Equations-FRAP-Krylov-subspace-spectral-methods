function b=mydct2(arg1)

% Basic algorithm.
b = mydct(mydct(arg1).').';
