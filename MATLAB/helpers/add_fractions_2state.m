function [holo, apoB] = add_fractions_2state(Kd, L0, RNA0)
% from Reining et al. (2013) Nature
% Supplementary Method for P[holo] (i.e. SE) derivation
%
% takes Kd, [L0], [RNA0]
% returns [holo], [apoB], not SE

temp = RNA0 + L0 + Kd;
holo = (temp - sqrt(temp^2 - 4 * RNA0 * L0)) / 2;
apoB = RNA0 - holo;

