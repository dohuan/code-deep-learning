function plot_q(t1,t2)

n_step=240;
max_t=1.2;
t=0:max_t/n_step:max_t;

q_t=zeros(n_step+1,1);

for i=1:n_step+1
    t_val=t(i);
    if  t_val >= 0 && t_val< t1
        q_t(i)=1.0;
    elseif t_val>= t1 && t_val<t2
        theta=pi*(t_val-t1)/(t2-t1);
        q_t(i)=0.5*(cos(theta)+1);
    else
        q_t(i)=0.0
    end
end

plot(t, q_t)