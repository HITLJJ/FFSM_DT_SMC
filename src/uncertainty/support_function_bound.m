function out = support_function_bound(J, set)
%SUPPORT_FUNCTION_BOUND Box support-function linearization radius.
% J is n-channel by 6. With numeric J (normally finite-difference derived),
% this is a CANDIDATE despite the box support formula itself being exact.

assert(size(J,2)==6,'FFSM:Dimension','J must have six parameter columns.');
out.linear_radius=abs(J)*set.halfwidth;
out.status='CANDIDATE';
out.reason='Parameter Jacobian is numerically differentiated; no global remainder certificate is implied.';
end
