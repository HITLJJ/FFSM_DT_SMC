function measurement=plant_measurement_packet(time,x,sequence)
%PLANT_MEASUREMENT_PACKET One-way immutable-value Plant measurement packet.
% Values are copied into a new struct. Twin estimates must never overwrite the
% raw_Q/raw_Qdot fields or write into the Plant integrator state.
measurement=struct('time',time,'sequence',sequence, ...
    'raw_Q',x(1:5),'raw_Qdot',x(6:10),'source','SIMULATED_PLANT');
end
