function validate_parameters(p)
%VALIDATE_PARAMETERS Validate physical and geometric parameter values.

positive = [p.m0 p.m1 p.m2 p.J0 p.J1 p.J2 p.l1 p.l2 p.c1 p.c2];
assert(all(isfinite(positive)) && all(positive > 0), ...
    'FFSM:InvalidParameters', 'Masses, inertias, lengths and COM offsets must be positive.');
assert(isfinite(p.a0) && p.a0 >= 0, 'FFSM:InvalidParameters', ...
    'The shoulder offset a0 must be finite and nonnegative.');
assert(p.c1 <= p.l1 && p.c2 <= p.l2, 'FFSM:InvalidParameters', ...
    'Each link COM must lie no farther than the corresponding link length.');

requiredPayload = {'enabled','model','mass','J'};
assert(all(isfield(p.payload, requiredPayload)), 'FFSM:InvalidPayload', ...
    'Payload configuration is incomplete.');
assert(any(strcmp(p.payload.model, {'point','rigid'})), 'FFSM:InvalidPayload', ...
    'Payload model must be ''point'' or ''rigid''.');
assert(p.payload.mass >= 0 && p.payload.J >= 0, 'FFSM:InvalidPayload', ...
    'Payload mass and inertia must be nonnegative.');
if p.payload.enabled
    assert(p.payload.mass > 0, 'FFSM:InvalidPayload', ...
        'An enabled payload must have positive mass.');
    if strcmp(p.payload.model, 'point')
        assert(p.payload.J == 0, 'FFSM:InvalidPayload', ...
            'A point payload must have zero rotational inertia.');
    end
end
end
