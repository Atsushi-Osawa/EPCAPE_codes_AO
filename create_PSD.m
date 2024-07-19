% TITLE: create_PSD
%
% PURPOSE: creates number and mass size distributions based on lognormal
% mode fitting parameters (number, mean size, and width).
%
% INPUT: D = diameter (microns), N = number concentration (cm^-3), mu =
% geometric mean diameter (microns), sig = geometric standard deviation
% (unitless). 
%
% OUTPUT: [PNSD = particle number size distribution, PMSD = particle mass
% size distribution].
%
% AUTHOR: Jeramy Dedrick, Scripps Institution of Oceanography, La Jolla,
% CA.
%
% CREATED: March 30, 2021. EDITED: April 5, 2021. 



function [PNSD, ...
         PMSD] = create_PSD(D, N, mu, sig)
     
      
     for i = 1:length(N)
     
     % lognormal mode number size distribution (S&P 2006)    
     PNSD(:,i) = N(i) ./ ...
           (sqrt(2 .* pi) .* log10(sig(i))) .* ...
           exp(-((log10(D)-log10(mu(i))).^2 ) ./ (2 .* log10(sig(i)).^2));  
       
     % lognormal mode mass size distribution (spherical particle
     % homogeneity with unit density)
     PMSD(:,i) = PNSD(:,i) .* ...
            (pi / 6) .* ...
            (1) .* ...
            (D .* 1e-6) .^3 .* ...
            (1e6 .* 1e6) .* 1e6;  
        
     end
        
        
end
%PMSD is the mass size distribution
     