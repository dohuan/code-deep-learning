function out = shp1d(index, dev , pt, sh, n_elm, type)

% SHP1D compute shape functions and their derivertives 
% for Lagrange linear and quadratic test function.
% index----test function index 
% pt----calculating point
% type----linear or quadratic
% dev----derivertives : 0 or 1

switch(type)
case 1
        max_n = n_elm+1;
case 2
       max_n = 2*n_elm +1;
otherwise
    error('Unknown type of shape function');
end

xi=(index-1)*sh;

if (dev ~= 0 & dev ~= 1)
   error('Wrong value for (dev)!!');
end

switch(type)
case  1
   if (pt>=xi-sh &  pt < xi)
          if dev == 0
             out = (pt-xi+sh)/sh;
         else 
             out = 1/sh;
          end
   elseif (pt>=xi & pt < xi+sh)
          if dev == 0
             out = (xi+sh-pt)/sh;
         else 
             out = -1/sh;
          end
   else 
	  out=0.0;
   end
case  2
   % Quadratic functions and their derivertives    
   if (mod(index,2)==1)
       if( pt >= xi-2*sh & pt < xi )
           if dev == 0
               out = (pt-xi+sh)*(pt-xi+2.0*sh)/(2*sh^2);
           else 
               out = (2*pt-2*xi+3*sh)/(2*sh^2);
           end
       elseif (pt >=xi & pt < xi+2*sh)
           if dev == 0
               out = (pt-xi-2.0*sh)*(pt-xi-sh)/(2*sh^2) ;
           else
               out = (2*pt-2*xi-3.0*sh)/(2*sh^2);
           end
       else
           out=0.0;
       end
   elseif (mod(index,2)==0)
       if(pt>=xi-sh & pt < xi+sh)
           if dev==0
               out = -(pt-xi+sh)*(pt-xi-sh)/sh^2;
           else
               out = -2.0*(pt-xi)/sh^2 ;
           end
       else
           out=0.0;
       end
   else
       out = 0.0
   end    
otherwise
   error('Error: Unknown shape function!!!!');
end
