function [holo, apoA, apoB] = add_fractions_3state(Kd, Kpre, L0, RNA0)
% from Reining et al. (2013) Nature
% Supplementary Method for P[holo] (i.e. SE) derivation
%
% takes Kd, Kpre, [L0], [RNA0]
% returns [holo], [apoA], [apoB], not SE

temp = Kd * (1 + Kpre) + Kpre * (L0 + RNA0);
holo = 1 / (2 * Kpre) * (temp - sqrt(-4 * Kpre^2 * L0 * RNA0 + temp^2));
L = L0 - holo;
apoA = holo * Kd / L;
apoB = holo * Kd / (Kpre * L);

