% Find the maximum in the p.
% Cumsum: cummulative sum of elements
% Rand: returns a random integer between 0 and 1
% eps:  returns the distance from 1.0 to the next largest double-precision number

function a = choose(p)

a = max(find([-eps cumsum(p)] < rand));
