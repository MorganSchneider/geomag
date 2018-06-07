function B = findzero(r, theta_init, phi, m_MF, t)

Btot = synth_values([r 6000], [theta_init 70], [phi 70], m_MF, [t 5082]);
B = Btot(1,1);

end