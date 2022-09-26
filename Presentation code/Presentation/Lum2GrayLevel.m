function y = Lum2GrayLevel(x,Cg,gam,b0)

y = ((x-b0)./Cg).^(1/gam);