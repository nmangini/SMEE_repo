%*  James Healy 6/18/2012
%*  Calculates the spin weighted (ss) spherical harmonic mode                  
%*  (ll, mm) at (th,ph).                                                       
%*  See:                                                                       
%*  http://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics#Calculating 
%*  for the equation and a Mathematica notebook containing this expression.
%*  Works with th,ph as meshes
%
function number=calcSWSH( ss, ll, mm, th, ph )

    all_coeff = (-1)^(mm) * sqrt( factorial(ll+mm)*factorial(ll-mm)*(2*ll+1) / ...
                ( 4*pi*factorial(ll+ss)*factorial(ll-ss) ));

    summ = 0;
    bot = max( 0, mm-ss );
    top = min( ll+mm, ll-ss );

    for i=bot:top
        sum_coeff = nchoosek( ll-ss, i ) * nchoosek( ll+ss, i+ss-mm );
        summ = summ + sum_coeff*((-1)^(ll-i-ss)) * cos(th/2).^(2*i+ss-mm) .* sin(th/2).^(2*(ll-i)+mm-ss);
    end

    realY = all_coeff*summ.*cos(mm*ph);
    imagY = all_coeff*summ.*sin(mm*ph);

    number = realY + sqrt(-1)*imagY;

end 
