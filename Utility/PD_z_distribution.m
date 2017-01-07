function  PD_z_distribution(fluid_data, flg)
 
format short e;

n_elm=150;             %number of element

Data_t0 = get_fluid_data(fluid_data, 'z', n_elm);
figure

if flg==1  % radius
    a=mean(Data_t0(:,2))
   plot(Data_t0(:,1), Data_t0(:,2));  %(mm, mm)
elseif flg==2 % pressure
    a=mean(Data_t0(:,3))
    plot(Data_t0(:,1)*1000, Data_t0(:,3)/1000); %(mm, kPa)
else  % magnitude of shear
    a=mean(Data_t0(:,4))
    plot(Data_t0(:,1)*1000, Data_t0(:,4));  %(mm, Pa)
end
 
clear *