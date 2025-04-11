function V=viscor1(d,hd,vplas,cpar,viscpar)
dOff = 2.4;
dCrit = 10.5;
d50 = 100;
eAmp = 1.1;
eWidth = 0.03;
ePeak = 0.6;
eHD = 1.18;
wMax = 2.6;

wAs = 0.;
if (dOff < d)
    wAs = wMax*(d-dOff)/(d+d50-2*dOff);
end
    
wPeak = 0.;
if (d > dOff && d <= dCrit) 
    wPeak = eAmp*(d-dOff)/(dCrit-dOff);

elseif (dCrit < d)
    wPeak = eAmp*exp(-eWidth*(d-dCrit));
end


wPH = wAs + wPeak*ePeak;
wEFF = wAs + wPeak*(1 + hd*eHD);
dPH = d - 2*wPH;
    
hdref = 0.45;
C = (cpar(1) + exp(cpar(2)*dPH)) * (-1. + 1./(1. + 10^cpar(3) * dPH^cpar(4))) + 1./(1. + 10^cpar(3) * dPH^cpar(4));
eta45 = viscpar(1) * exp(viscpar(2)*dPH) + viscpar(3) + viscpar(4) * exp(viscpar(5) * dPH^viscpar(6)); 
hdfac = (1-hd^C - 1.)/(1-hdref^C - 1.); 
etaVitro = 1. + (eta45 - 1.) * hdfac;
    
etaVivo = etaVitro*((d/(d-2*wEFF))^4)*vplas;
V=etaVivo;

end