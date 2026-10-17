function status=twin_validity_check(twin,currentTime)
%TWIN_VALIDITY_CHECK Basic runtime validity; not a containment certificate.
if any(~isfinite(twin.theta_hat)) || any(twin.theta_hat<=0)
    status='FAILED'; return
end
if any(twin.parameter_set.lower<=0) || any(twin.parameter_set.lower>twin.parameter_set.upper)
    status='FAILED'; return
end
if ~isempty(twin.previous_time) && currentTime-twin.previous_time>2*twin.options.update_period
    status='STALE'; return
end
if strcmpi(twin.mode,'online_update')
    status='EXPLORATORY';
else
    status='VALID_FIXED_MODEL';
end
end
