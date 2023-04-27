%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NAME OF FILE: ppwd_live_link_COMSOL_single_freq_all_param_non_intr_POD_4_Github.m
%
% PURPOSE:  
% Reference: IntroductionToLiveLinkForMATLAB.pdf
%          LiveLinkForMATLABUsersGuide.pdf
%          COMSOL_ProgrammingReferenceManual.pdf
%
% Written by Gladwin Jos
%
% Date : 16/06/2022 
% Date: 22/02/2023 :  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) This MATLAB code works only with COMSOL livelink for MATLAB Software
% 2) The reference paper for this MATLAB code is
%    "A Non Intrusive  POD-PCE Scheme for Uncertainty Quantification
%     of Stochastic EM problems", IEEE Microwave and Wireless Technology
%     letters, 2023. Corresponding Author: Gladwin Jos

clc;
clear all;
close all;
tic;
format longe;
import com.comsol.model.*
import com.comsol.model.util.*
% COMSOL EM model for parallel palte waveguide
temp= mphload('parallel_plate_waveguide_25cm_3_pt_5cm_multilayer_single_freq_all_parameters.mph');
% Mean values of complex permittivity
epsilon_r_1_dash=2;
epsilon_r_2_dash=4;
epsilon_r_3_dash=3;
epsilon_r_4_dash=8;
epsilon_r_5_dash=6;
epsilon_r_1_dd=0.8;
epsilon_r_2_dd=0.12;
epsilon_r_3_dd=0.3;
epsilon_r_4_dd=0.22;
epsilon_r_5_dd=0.1;
freq_constant=3.9*10^9; 
No_of_dielectric_layer=5;  
No_of_rv=2*No_of_dielectric_layer;
monte_carlo_iteration=10000;
gaussian_RV=zeros(monte_carlo_iteration,No_of_dielectric_layer);
for ii=1:No_of_rv
    seed_value=ii;
    gaussian_RV(:,ii)=random_number_seed_gaussian(monte_carlo_iteration,seed_value);
end
xi_value = gaussian_RV(:,1:No_of_rv); 
poly_order=2;
% poly_order=3; % change accordingly
no_of_polynomial=factorial(poly_order+No_of_rv)/...
    (factorial(poly_order)*factorial(No_of_rv));
if poly_order==2
    No_of_LHS_samples=100; % ~2p (2*66)
end
if poly_order==3
    No_of_LHS_samples=500; % ~2p (2*286)
end
% Obtain Latin hypercube samples
[lhs_samples,min_length_lhs]= ...
        obtain_lhs_samples_final(No_of_LHS_samples ,No_of_rv ,xi_value);

temp.param.set('f1', strcat(num2str(freq_constant),'')); 
geom1=temp.geom.get('geom1'); 
geom1.run;
mesh1=temp.mesh.get('mesh1'); 
mesh1.run;
[meshstats,meshdata] = mphmeshstats(temp); %  LiveLinkForMATLABUsersGuide.pdf
tri_elm = meshdata.elem{2};
vtx = meshdata.vertex;
n=transpose(tri_elm)+1;
node_list=[transpose(vtx) zeros(length(vtx),1)];
figure(1) 
colorstring = 'kbgrmwy';
patch('Faces', n, 'Vertices', node_list, ...
    'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'w');
hold on;
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
% for dielectric dimension
Y2=2.5/100;
Dw=1/100
L_1=10/100; % in meters
L_2=15/100; 
L_3=0/100; 
L_4=Y2; 
% plotting dielectric
line([L_1 L_1],[L_3 L_4],'LineWidth',2.5,'Color','r');
hold on;
line([L_1 L_2],[L_4 L_4],'LineWidth',2.5,'Color','r');
hold on;
line([L_2 L_2],[L_3 L_4],'LineWidth',2.5,'Color','r');
hold on;
line([L_1+Dw L_1+Dw],[L_3 L_4],'LineWidth',2.5,'Color','r');
hold on;
line([L_1+2*Dw L_1+2*Dw],[L_3 L_4],'LineWidth',2.5,'Color','r');
hold on;
line([L_1+3*Dw L_1+3*Dw],[L_3 L_4],'LineWidth',2.5,'Color','r');
hold on;
line([L_1+4*Dw L_1+4*Dw],[L_3 L_4],'LineWidth',2.5,'Color','r');
hold on;
%%%%%%%%%%%%%%%%%%%% Analysis of mesh data %%%%%%%%%%%%%%%%%%%%%%
Ne=length(tri_elm)
No_of_nodes=length(vtx)
tri_edge=[1 2;2 3;3 1];
local_edge= edge_array_local_2D(n,Ne,tri_edge);      
[global_edge,local_edge_array_elm_no ] = ...
    edge_array_global_2D_unique(local_edge);
No_of_edges=length(global_edge) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Study the EM model
std1=temp.study.get('std1'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electric field, Ein for Pin=1W for Global matrix for Gamma
Ein= -1.481670943678778e+02 + 3.936113124114703e+00i
sigma_i_eps_r_dash=[0.1;0.1;0.1;0.1;0.1]; 
sigma_i_eps_r_dd=[0.1;0.1;0.1;0.1;0.1];  
% Generate the data matrix as in equation (1)
U_matrix=zeros(No_of_nodes,min_length_lhs);
for kk=1:min_length_lhs
    kk
    lhs_samples_itr=lhs_samples(kk,:);
    epsilon_r_1_dash_stoc=epsilon_r_1_dash ...
             +sigma_i_eps_r_dash(1,1)*lhs_samples_itr(1,1); 
    epsilon_r_1_dd_stoc=epsilon_r_1_dd+sigma_i_eps_r_dd(1,1)* ...
               lhs_samples_itr(1,1+No_of_dielectric_layer);
   
    [n,k]=refractive_index_COMSOL_function(epsilon_r_1_dash_stoc,epsilon_r_1_dd_stoc);
    temp.param.set('n1', strcat(num2str(n),'')); 
    temp.param.set('ki1', strcat(num2str(k),'')); 
    epsilon_r_2_dash_stoc=epsilon_r_2_dash ...
               +sigma_i_eps_r_dash(2,1)*lhs_samples_itr(1,2); 
    epsilon_r_2_dd_stoc=epsilon_r_2_dd+sigma_i_eps_r_dd(2,1)* ...
                lhs_samples_itr(1,2+No_of_dielectric_layer);
        
    [n,k]=refractive_index_COMSOL_function(epsilon_r_2_dash_stoc,epsilon_r_2_dd_stoc);
    temp.param.set('n2', strcat(num2str(n),'')); 
    temp.param.set('ki2', strcat(num2str(k),''));
    epsilon_r_3_dash_stoc=epsilon_r_3_dash ...
                 +sigma_i_eps_r_dash(3,1)*lhs_samples_itr(1,3); 
    epsilon_r_3_dd_stoc=epsilon_r_3_dd+sigma_i_eps_r_dd(3,1)* ...
                       lhs_samples_itr(1,3+No_of_dielectric_layer); 
        
    [n,k]=refractive_index_COMSOL_function(epsilon_r_3_dash_stoc,epsilon_r_3_dd_stoc);
    temp.param.set('n3', strcat(num2str(n),'')); 
    temp.param.set('ki3', strcat(num2str(k),'')); 
    epsilon_r_4_dash_stoc=epsilon_r_4_dash ...
                     +sigma_i_eps_r_dash(4,1)*lhs_samples_itr(1,4); 
    epsilon_r_4_dd_stoc=epsilon_r_4_dd+sigma_i_eps_r_dd(4,1)* ...
                    lhs_samples_itr(1,4+No_of_dielectric_layer);
       
    [n,k]=refractive_index_COMSOL_function(epsilon_r_4_dash_stoc,epsilon_r_4_dd_stoc);
    temp.param.set('n4', strcat(num2str(n),'')); 
    temp.param.set('ki4', strcat(num2str(k),''));
    epsilon_r_5_dash_stoc=epsilon_r_5_dash ...
                   +sigma_i_eps_r_dash(5,1)*lhs_samples_itr(1,5); 
    epsilon_r_5_dd_stoc=epsilon_r_5_dd+sigma_i_eps_r_dd(5,1)* ...
                   lhs_samples_itr(1,5+No_of_dielectric_layer);
      
    [n,k]=refractive_index_COMSOL_function(epsilon_r_5_dash_stoc,epsilon_r_5_dd_stoc);
    temp.param.set('n5', strcat(num2str(n),'')); 
    temp.param.set('ki5', strcat(num2str(k),'')); 
    temp.sol('sol1').run; 
%     M = mphevalglobalmatrix(temp,'n1')
%     mphgetcoords(temp, 'geom1', 'domain', 2)
%     mphglobal(temp)
    % Edim    Refer: LiveLinkMATLABUsersGuide.pdf
    Ey_field_domain=mpheval(temp,'emw.Ey' ,'Edim','domain','selection',[1 2 3 4 5 6 ]); 
    if kk==1
        Ey_field_port1=mpheval(temp, 'emw.Ey' ,'Edim','boundary','selection',1 );
        
        % Getting coordinates from solution
        node_pt_domain=Ey_field_domain.p;
        node_pt_domain_x=node_pt_domain(1,:);
        node_pt_domain_y=node_pt_domain(2,:);
        %%%%%%%% finding duplicates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [uniqe_cord, unique_index,new_index]=unique([node_pt_domain_x' node_pt_domain_y'],'rows');
        new_cord_x=uniqe_cord(:,1);
        new_cord_y=uniqe_cord(: ,2); 
        figure(1)
        plot(new_cord_x,new_cord_y,'.k','MarkerSize',12);
        hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        node_pt_port1=Ey_field_port1.p;
        node_pt_port1_x=node_pt_port1(1,:);
        node_pt_port1_y=node_pt_port1(2,:); 
        figure(1)
        plot(node_pt_port1_x,node_pt_port1_y,'sr','MarkerSize',14);
        %%%%%%% Find index from global Electric field %%%%%%%%%%%%%%%
        old_index=zeros(length(node_pt_port1_x),1);
        for ll=1:length(node_pt_port1_x)
                x_pt_port=node_pt_port1_x(ll);
                y_pt_port=node_pt_port1_y(ll);
                for mm=1: No_of_nodes
                    if x_pt_port==node_pt_domain_x(mm) && ...
                            y_pt_port==node_pt_domain_y(mm)
                            old_index(ll)=mm;
                        break;
                    end
                end
        end
        
    end
    Ey_field_domain_value=Ey_field_domain.d1;
    index_mapping_port1=new_index(old_index);
    Ey_field_unique =Ey_field_domain_value(unique_index);
    U_matrix(:,kk)=transpose(Ey_field_unique);     
end 

[SVD_U,SVD_Sigma,SVD_V_T]=svd(U_matrix);
sing_values=diag(SVD_Sigma);
figure(2) 
colorstring = 'krbgmwy'; 
x_temp=linspace(1,length(sing_values),length(sing_values));
bar(x_temp,sing_values,'b' ); 
hold on;
for ii=1:length(sing_values)
    t_1(1)=plot(x_temp(ii),sing_values(ii),'o' ,'LineWidth',1,'Color',colorstring(1));
    hold on;
end
legend([t_1(1) ],'Singular values' ,'Interpreter','latex');


count_sing_value=1;
N_value=length(sing_values);
RIC_value=zeros(N_value,1);
for ii =1:N_value
    index=linspace(1,ii,ii);
    RIC_value(ii)=sum(sing_values(index).*sing_values(index)) ...
        /sum(sing_values.*sing_values);
end

err_value_1=0.9999999;
err_value_2=0.9999990;
err_value_3=0.9999900;
err_value_4=0.9999000;
err_value_5=0.9990000;
err_value_6=0.9900000;
err_value_7=0.9000000;
err_value_8=0.8000000;
err_value=[err_value_1;err_value_2;err_value_3;err_value_4;...
    err_value_5;err_value_6;err_value_7;err_value_8]; 
m_count_list=zeros(length(err_value),1);
for jj=1:length(err_value)
    break_loop=0;
    for ii=1:N_value 
         if RIC_value(ii,1) >=err_value(jj,1) && break_loop==0
             m_count_list(jj,1)=ii;
             break_loop=1;
         end
    end
end
m_count=m_count_list(1)
% POD basis
U_POD=SVD_U(:,1:m_count); 
%=========== non-intrusive POD ============================================
% Actual POD coefficients from equation (9)
alpha_i=zeros(m_count,min_length_lhs);
for k=1:min_length_lhs
    alpha_i(:,k)=U_POD'* U_matrix(:,k);
end
%===============  orthogonality ================================= 
round((U_POD)'*U_POD)
%==================================================================
% Generate stochastic orthogonal polynomials
Psi_matrix_lhs=get_orthogonal_polynomial_for_Nr_10(no_of_polynomial,min_length_lhs,lhs_samples);
Psi_matrix=get_orthogonal_polynomial_for_Nr_10(no_of_polynomial,monte_carlo_iteration,xi_value);
% Solving for optimal POD coefficients as in equation (13)
alpha_coeff=zeros(no_of_polynomial,m_count);
for kk=1:m_count
    kk
    alpha_cap=transpose(alpha_i(kk,:));
    alpha_coeff(:,kk)=(transpose(Psi_matrix_lhs)*Psi_matrix_lhs)\(transpose(Psi_matrix_lhs)*alpha_cap);
end

gamma_POD=zeros(monte_carlo_iteration,1);
for jj=1:monte_carlo_iteration
    alpha_final=zeros(m_count,1);
    for kk=1:m_count 
        for ii=1:no_of_polynomial
            alpha_final(kk,1)=alpha_final(kk,1)+alpha_coeff(ii,kk)*Psi_matrix(jj,ii);
        end
    end
    % Generation of stochastic electric field as in equation (4)
    % from the estimated POD coefficients
    stochastic_E_field_non_intrusive=U_POD*alpha_final; 
    % Derive reflection ceofficients from the stochastic electric field.
    Ey_initial_POD_non_intrusive=stochastic_E_field_non_intrusive(index_mapping_port1);
    E_port_frm_global=sum(Ey_initial_POD_non_intrusive)/length(Ey_initial_POD_non_intrusive);
    gamma_POD(jj,1)=abs((E_port_frm_global*Ein-(Ein*Ein))/(Ein*Ein));
end
 
 



