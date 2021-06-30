% spDist_postHoc_powerAnalysis


%part one: estimate our power using effect size form Li et al 2021
% using Figure 4b, v3ab, strongest result
li_n= 14; % reported subj
v3ab_mu = .1086; %estimated from observable means
v3ab_sem = .018; %estimated from observable SEM
n= 14; % reported subj


% estimate effect size
% for cohen's d one sample 

t_v3ab = v3ab_mu /v3ab_sem;
d_z_v3ab = t_v3ab/sqrt(li_n);

v3_mu = .055; %estimated from observable means
v3_sem = .01; %estimated from observable SEM



% estimate effect size
% for cohen's d one sample 

t_v3 = v3_mu /v3_sem;
d_z_v3 = t_v3/sqrt(li_n);

%computing our effect size
v123_mu = .3151;
v123_sem = .0890;
my_n = 7;
t_v123 = v123_mu /v123_sem;
d_z_v123 = t_v123/sqrt(my_n);

% now plot power as a function of sample size, using fixed effect sizes
% from 0.6-1.4 
d_point6 = [ 12    15    18    26    34]; % sample sizes collectd by inputting effect size and changing power in g*power
d_point8 = [8     9    12    15    19];
d_one = [6     7     8    11    13];
d_onepoint2 = [ 5     5     6     8    10];
d_onepoint4 = [4     5     5     7     8];
my_power = [.6 .7 .8 .9 .95];

figure; 
plot(d_point6, my_power,'linewidth',1.5)
hold on;
plot(d_point8, my_power,'linewidth',1.5)
plot(d_one, my_power,'linewidth',1.5)
plot(d_onepoint2, my_power,'linewidth',1.5)
plot(d_onepoint4, my_power,'linewidth',1.5)
legend('d=0.6','d=0.8','d=1','d=1.2','d=1.4','Fontsize',14)
xlabel('Sample Size')
ylabel('Power (1-beta)')
xlim([4 20])


%% finer discrim

my_power_ax = [.6 .65 .7 .75 .8 .85 .9 .95];

dpoint5 = [16 18 21 23 27 31 36 45];
dpoint6 = [12 13 15 17 19 22 26 32];
dpoint7 = [9 10 12 13 15 17 19 24];
dpoint8 = [8 8 9 10 12 13 15 19];
dpoint9 = [ 7 7 8 9 10 11 13 15];
done = [6 6 7 7 8 9 11 13];
donepointone = [ 5 6 6 7 7 8 9 11];
donepointtwo = [5 5 5 6 6 7 8 10];
donepointthree = [4 5 5 5 6 6 7 8 ];
donepointfour = [4 4 5 5 5 6 7 8];
donepointfive= [ 4 4 4 5 5 5 6 7];
figure; hold on;
plot(dpoint5,my_power_ax,'linewidth',1.5)
hold on;
plot(dpoint6,my_power_ax,'linewidth',1.5)
plot(dpoint7,my_power_ax,'linewidth',1.5)
plot(dpoint8,my_power_ax,'linewidth',1.5)
plot(dpoint9,my_power_ax,'linewidth',1.5)
plot(done,my_power_ax,'linewidth',1.5)
plot(donepointone,my_power_ax,'linewidth',1.5)
plot(donepointtwo,my_power_ax,'linewidth',1.5)
plot(donepointthree,my_power_ax,'linewidth',1.5)
plot(donepointfour,my_power_ax,'linewidth',1.5)
plot(donepointfive,my_power_ax,'linewidth',1.5)
plot(7,.9, 'k+')
legend('d=0.5','d=0.6','d=0.7','d=0.8','d=0.9','d=1','d=1.1','d=1.2','d=1.3','d=1.4','d=1.5','Fontsize',14)
xlabel('Sample Size')
ylabel('Power (1-beta)')
xlim([4 20])

% figure; hold on;
% plot(my_power_ax,dpoint5,'linewidth',1.5)
% hold on;
% plot(my_power_ax,dpoint6,'linewidth',1.5)
% plot(my_power_ax,dpoint7,'linewidth',1.5)
% plot(my_power_ax,dpoint8,'linewidth',1.5)
% plot(my_power_ax,dpoint9,'linewidth',1.5)
% plot(my_power_ax,done,'linewidth',1.5)
% plot(my_power_ax,donepointone,'linewidth',1.5)
% plot(my_power_ax,donepointtwo,'linewidth',1.5)
% plot(my_power_ax,donepointthree,'linewidth',1.5)
% plot(my_power_ax,donepointfour,'linewidth',1.5)
% plot(my_power_ax,donepointfive,'linewidth',1.5)
% legend('d=0.5','d=0.6','d=0.7','d=0.8','d=0.9','d=1','d=1.1','d=1.2','d=1.3','d=1.4','d=1.5','Fontsize',14)