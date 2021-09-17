close all
clear all
%% 00 - Pre-Definitions
% Total Number of DoFs
  NDoFs=6;
% The DoF to Draw the Plots for
  DoFtoDraw=4;
% External Load Increments Vector
  ExternalLoadIncrements=[0; 0; 0; 0.1; 0; 0];
%% 01 - Calculations
SCAN='%f'; for i=1:NDoFs-1; SCAN=[SCAN '%f']; end
system('OpenSees.exe Example1OpenSees_A_Corotational.tcl');
fid=fopen('OSNodalDisplacements.txt','r');
a=textscan(fid,SCAN,'CollectOutput',1); Displacements=a{1}; clear a;
fclose(fid);
delete('OSNodalDisplacements.txt');
OpenSeesCorotationalForces=zeros(NDoFs,1);
OpenSeesCorotationalDisplacements=zeros(NDoFs,1);
for i=1:size(Displacements,1)
OpenSeesCorotationalDisplacements(:,i+1)=transpose(Displacements(i,:));
OpenSeesCorotationalForces(:,i+1)=(i)*ExternalLoadIncrements;
end; clear Displacements;
OpenSeesCorotationalDisplacements=abs(OpenSeesCorotationalDisplacements);
OpenSeesCorotationalForces=abs(OpenSeesCorotationalForces);
system('OpenSees.exe Example1OpenSees_B_PDelta.tcl');
fid=fopen('OSNodalDisplacements.txt','r');
a=textscan(fid,SCAN,'CollectOutput',1); Displacements=a{1}; clear a;
fclose(fid);
delete('OSNodalDisplacements.txt');
OpenSeesPDeltaForces=zeros(NDoFs,1);
OpenSeesPDeltalDisplacements=zeros(NDoFs,1);
for i=1:size(Displacements,1)
OpenSeesPDeltalDisplacements(:,i+1)=transpose(Displacements(i,:));
OpenSeesPDeltaForces(:,i+1)=(i)*ExternalLoadIncrements;
end; clear Displacements;
OpenSeesPDeltalDisplacements=abs(OpenSeesPDeltalDisplacements);
OpenSeesPDeltaForces=abs(OpenSeesPDeltaForces);
clear b; clear fid; clear i;
figure('units','normalized','outerposition',[0 0 1 1]); 
ylim([0, 1.05*max(OpenSeesPDeltaForces(DoFtoDraw,size(OpenSeesPDeltaForces,2)),OpenSeesCorotationalForces(DoFtoDraw,size(OpenSeesCorotationalDisplacements,2)))]); 
xlim([0, 1.05*max(OpenSeesPDeltalDisplacements(DoFtoDraw,size(OpenSeesPDeltalDisplacements,2)),OpenSeesCorotationalDisplacements(DoFtoDraw,size(OpenSeesCorotationalDisplacements,2)))]);
grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
Title=['Load-Dispacement Plot for DoF ' num2str(DoFtoDraw) '']; title(Title,'fontsize',15,'fontweight','bold'); hold on;
xlabel('Displacement [in]','fontsize',13,'fontweight','bold'); ylabel('Load [lb]','fontsize',13,'fontweight','bold'); 
P1=plot(OpenSeesCorotationalDisplacements(DoFtoDraw,:),OpenSeesCorotationalForces(DoFtoDraw,:),'b','LineWidth',3);
P2=plot(OpenSeesPDeltalDisplacements(DoFtoDraw,:),OpenSeesPDeltaForces(DoFtoDraw,:),'c','LineWidth',3);
leg=legend([P1, P2],{'Corotational Formulation', 'P-Delta Formulation'},'Location','northwest','fontsize',15); 
Title=['DoF ' num2str(DoFtoDraw) ' Load-Dispacement Plot']; hold on;
Plot_Name=['OpenSees Results'];
h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;
clear ans; clear ax; clear DoFtoDraw; clear ExternalLoadIncrements; clear h1; clear leg; clear NDoFs; clear P1; clear P2; clear Plot_Name; clear SCAN; clear Title;
%% 02 - Save the Results
save('Example1Results_B_GeomNL_OpenSees')