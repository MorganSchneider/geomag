function alpha = find_alpha(Br_synth, r, theta, phi, g, N)

[Br, ~, ~] = find_B(r, theta, phi, g, N);
alpha = zeros(length(r), 1);
for i = 1:length(r)
    alpha(i) = Br_synth(i) - Br(i); % once I start using actual data, Br_synth will be Br_EEJ (= 0).
end

return