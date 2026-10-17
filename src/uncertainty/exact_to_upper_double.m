function value = exact_to_upper_double(exactValue)
%EXACT_TO_UPPER_DOUBLE Convert nonnegative exact sym to an enclosing double.
% The returned binary64 number is checked symbolically and stepped upward until
% it is no smaller than the exact rational certificate.

assert(isa(exactValue,'sym'),'FFSM:ExactRequired','Input must be symbolic exact.');
assert(isAlways(exactValue>=0),'FFSM:NonnegativeRequired','Certificate must be nonnegative.');
value=double(vpa(exactValue,40));
while ~isAlways(sym(value,'f')>=exactValue)
    value=value+eps(value);
end
end
