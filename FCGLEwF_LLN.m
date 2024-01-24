%% Lugiato-Lefever nonlinear RHS
function dy = LLwF_LLN(y,yT,et,ph,Delta)
    dy = -(1 - 1i*Delta).*abs(y).^2.*y + et.*exp(-1i.*ph).*yT;
end