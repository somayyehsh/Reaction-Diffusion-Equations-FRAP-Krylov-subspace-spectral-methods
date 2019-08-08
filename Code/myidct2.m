function a = myidct2(arg1)

% Basic algorithm.
a = myidct(myidct(arg1).').';
