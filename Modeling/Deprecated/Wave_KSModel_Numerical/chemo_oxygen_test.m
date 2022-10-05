Asp=50:25:200;
Time=20:20/(length(Asp)-1):40;
tic
% parfor i=1:1%length(Asp)
%     [wavespeed(i) final_profile(i) dyn_profile(i) dx]=chemo_oxygen_fun(Asp(i));
% 
% end
% %%

%%

 for i=1:length(Asp)
     [wavespeed(i) final_profile(i) dyn_profile(i) dx]=chemo_oxygen_fun_v2_no_oxg_grad(Asp(i),Time(i)*60);
      

 end
 for i=1:length(Asp)
    cell_density(i,:)=sum(final_profile(i).cell_density);
    [denmax(i),denmaxind(i)]=max(cell_density(i,60:end));
    wavedensity(i)=sum(cell_density(i,denmaxind(i)+40:end));
    cell_number(i)=wavedensity(i)*(dx*6*8.4/10*1.4);
 end
