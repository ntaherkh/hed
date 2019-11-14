clc
clear
close all

SNR= [0 10 15 20 25 30 35 ];
loop=100000
N=128;
N_t=130;
sishma_orginal=0.2752;
M_qam=16;
M_PPM=4;
N_sub=N/4;

I_H=1;
lambda_h=2;
lambda_l=2;
delta_org=30;
sishma_tuned_org=(1/delta_org)*I_H;
I_OFF_org=(lambda_l*sishma_tuned_org);
% I_ON_org=(I_H*0.5)+(lambda_l*sishma_tuned_org);
I_ON_org=(I_H)-(lambda_h*sishma_tuned_org);
delta_last=0.15*floor((min((1/((3*lambda_l/delta_org)+(lambda_h/delta_org))),(1/((3*lambda_h/delta_org)+(lambda_l/delta_org)))))/0.15)
tau=M_PPM/N_sub;
Count_location=nchoosek([1:N_sub],M_PPM); %  can we optimize N_sub and M_PPM ??
b_MPPM=floor(log2(length(Count_location))); %% bit in MPPM
M_location=Count_location(1:2^b_MPPM,:);
b_dco=(N/2)*log2(M_qam); %%  bit in dco
% I_ON=0.75*I_H; % 50% percent
%  % 50% percent
% I_m=(I_ON+I_OFF)/2;
% %% 
% lambda_l=(I_OFF-I_L)/sishma;
% lambda_h=(I_m-I_OFF)/sishma;
% lambda_l=(I_ON-I_m)/sishma;
% lambda_h=(I_H-I_ON)/sishma;
% % lambda_l=10;
% lambda_h=10;
gamma_index=0;
for gamma=1:.15:delta_last;
    gamma_index=gamma_index+1;
    ratio(gamma_index)=delta_org/gamma
    sishma_tuned=sishma_tuned_org*gamma;
    normalize_coff=sishma_tuned/sishma_orginal;
    sishma=sishma_tuned;
    sishma_tuned_vect(gamma_index)=sishma_tuned;
    I_OFF=(lambda_l*sishma_tuned);
    I_ON=(I_H)-(lambda_h*sishma_tuned);
    I_m=(I_ON+I_OFF)*0.5
%     I_L=(1-delta)*I_H;
    p_1=(sishma^2)*(((lambda_l^2)*qfunc(lambda_l))+((lambda_h^2)*qfunc(lambda_h)));
    p_2=(sishma^2)*(1/sqrt(2*pi))*(-((exp(-(lambda_l^2)/2))*lambda_l)-((exp(-(lambda_h^2)/2))*lambda_h));
    p_3=(sishma^2)*(1-qfunc(lambda_l)-qfunc(lambda_h));
    mean_clipped =(sishma)*(((1/sqrt(2*pi))*((exp(-(lambda_l^2)/2))-(exp(-(lambda_h^2)/2))))-lambda_l*(qfunc(lambda_l))+lambda_h*qfunc(lambda_h))
    sq_clipped=p_1+p_2+p_3;
    sishma_clipped=sqrt(sq_clipped);
    sishma_clipped_vect(gamma_index)=sishma_clipped;
    mean_clipped_vect(gamma_index)=mean_clipped ;
    PP_index=0;
    PP_vect=0;
    dimm_level(gamma_index)=(tau*(I_OFF+mean_clipped))+((1-tau)*(I_OFF+mean_clipped ));
end
gamma=1:.15:delta_last
plot(gamma,dimm_level)