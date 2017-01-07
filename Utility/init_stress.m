function init_stress(filename)

a=load(filename);
n=size(a);

str1=strcat(filename,'.dat');
out_f=fopen(str1,'w');

for i=1:n(1)
    fprintf(out_f,'%10.5f\t%10.5f\t%10.5f\n',a(i,1), a(i,5)/a(i,7), a(i,6)/a(i,7));
end

fclose(out_f);
clear *