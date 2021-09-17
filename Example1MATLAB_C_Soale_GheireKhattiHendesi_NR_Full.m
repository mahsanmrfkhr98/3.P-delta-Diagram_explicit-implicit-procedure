clear all
%% units: lb, psi, in
%% 00 - Pre-Definitions
% Pause Durations when Plotting
  Pause1=0.15; Pause2=1; 
% Plot Axis Scale
  XScale=1.05; YScale=1.05;
% Total Number of Steps
 % Steps=10;
% Steps=20;
% Steps=50;
 Steps=1000;
% Tolerance for the Error of Calculations
  Tolerance=0.00001;
% The DoF to Draw the Plots for
  DoFtoDraw=4;
% Only Plot the Fist MATLAB Step
  ShowOnlyFirstStep=0; ConnectStartToEnd=0; ZoomSteps=10; ZoomPause=2; ZoomXScale=1.15; ZoomYScale=1.15;
% Plot All MATLAB Steps
  ShowAllSteps=0; NoAnimation=1; ConnectSteps=0; ShowTheJump=0; DisableTrialandFinalPoints=0; 
% Plot the Result Curve without Iterations
  ShowOnlyTheResultCurve=0;         
% Draw Stiffness Variation or Load-Displacement
  DrawStiffnessVariations=1;
%% 01 - Input
% Nodes Data Matrix:                               Node#       x[in]      y[in]
                                     Nodes=[
                                                    1,           0,         0
                                                    2,         160,         0
                                                    3,         320,         0];
                           
% Elements Data Matrix:                          Element#   Start_Node   End_Node  A[in^2]   E[psi]      N0[lb]                                
                                  Elements=[
                                                    1,           1          2,       5,     36000000     1500;
                                                    2,           2          3,       5,     36000000     1500];

% Restrains Vector:                                DoFs
                                 Restrains=[
                                                    1;
                                                    2;
                                                    5;
                                                    6];
                           
% External Loads Vector:                      Load_Value[lb]   DoF#  
                                     Loads=[
                                                    0           %1
                                                    0           %2
                                                    0           %3
                                                  200           %4
                                                    0           %5
                                                    0];         %6
%% 02 - Complie Inputs
% Nodes Data                                                                                               
  % Available Data                                                                        
    x=Nodes(:,2);                                                                                        
    y=Nodes(:,3);                                                                                       
  % New Data                                                                               
    % Number of Nodes                                                                
      NNodes=size(Nodes,1);                                                                   
      NDoFs=2*NNodes;                                                                   
    % Free DoFs                                                                           
      FreeDoFs=(1:NDoFs)';                                                                      
      FreeDoFs(Restrains)=[];                                                                   
% Elements Data
  % Available Data                                                                     
    Start=Elements(:,2);                                                                        
    End=Elements(:,3);                                                                          
    A=Elements(:,4);                                                                            
    E=Elements(:,5);
    N=Elements(:,6);                                                                                                                               % Nirooye mehvarie avalieye elemana, az setoone 6 matrice "Elements" rikhte mishe too bordare "N". | N = Ye brodar ke derayeye i, nirooye mehvarie avalieye elemane i hast.                                                                                                         
  % New Data                                                                           
    % Number of Elements                                                                   
      NElements=size(Elements,1);     
    % Geometric and Mechanical Info Based on Initial Outputs
      DtoV=zeros(NElements,NDoFs);  
      QtoP=zeros(NDoFs,NElements);
      for i=1:NElements                                                                                                                                                              
          DoFs(i,:)=[2*Start(i)-1  2*Start(i)  2*End(i)-1  2*End(i)];                                              
          invL(i)=1/((x(End(i))-x(Start(i)))^2+(y(End(i))-y(Start(i)))^2)^0.5;                                                                     % Avale kar, Az mokhtasate avalieye nodhaye sar o tahe har eleman estefade mishe.                                                             
          C(i)=(x(End(i))-x(Start(i)))*invL(i);                                                                                                    % Avale kar, Az mokhtasate avalieye nodhaye sar o tahe har eleman estefade mishe.                                                                                               
          S(i)=(y(End(i))-y(Start(i)))*invL(i);                                                                                                    % Avale kar, Az mokhtasate avalieye nodhaye sar o tahe har eleman estefade mishe.                                                                                             
          a(i,:)=[-C(i) -S(i) C(i) S(i)];                                                                                                                                                                                      
          Kel(i,i)=A(i)*E(i)*invL(i);                                                                                                              % Baraye matrice sakhtie locale elemana az toole avalie estefade mishe. Formulation, safheye 8 az file "Session 03.pdf".
      end 
      x0=x;                                                                                                                                        % Mokhtasate x avalieye nod'ha too bordare "x0" save mishe.        | x0    = Ye brodar ke derayeye i, mokhtasate x avalieye node i hast.                             
      y0=y;                                                                                                                                        % Mokhtasate y avalieye nod'ha too bordare "y0" save mishe.        | y0    = Ye brodar ke derayeye i, mokhtasate y avalieye node i hast.                                                                      
      N0=N;                                                                                                                                        % Nirooye mehvarie ebteda'iie elemana too bordare "N0" save mishe. | N0    = Ye brodar ke derayeye i, Nirooye mehvarie avalieye elemane i hast.         
%% 03 - Calculations                                                                                                                                  
 % 00 - Pre-Definitions                                                                               
        Pint=zeros(NDoFs,1);                                                                                
        deltaD=zeros(NDoFs,1);                                                                          
        D=zeros(NDoFs,1); 
        V=zeros(NElements,1);                                                                                                                      % Agar inja ta'rif nashe, khatte 112 error migirim.  
        k=0;                                                                                                                                                         
        Reached(1)=1;  
 % 01 - Solve
        tic;
        for s=1:Steps
            % Prepare the External Loads Vectors
              Pext=(s/Steps)*Loads;
              Pr=Pext-Pint; 
               while norm(Pr(FreeDoFs))>Tolerance
                     % Assemble the Global Stiffness Matrix
                       KL=zeros(NDoFs,NDoFs);  
                       KNL=zeros(NDoFs,NDoFs);  
                       for i=1:NElements                                                                                                                                                                                                                                                                             
                           KL(DoFs(i,:),DoFs(i,:))=KL(DoFs(i,:),DoFs(i,:))+(transpose(a(i,:))*Kel(i,i)*a(i,:));                                    % Az zaviehaye jadid va toole avalie estefade mishe. Formulation, safheye 14 az file "Session 03.pdf".                       
                           KNL(DoFs(i,:),DoFs(i,:))=KNL(DoFs(i,:),DoFs(i,:))+(N(i)*invL(i)*[1, 0, -1, 0; 0, 1, 0, -1; -1, 0, 1, 0; 0, -1, 0, 1]);  % Az niroo mehvarie fe'li va toole fe'li estefade mishe. Formulation, safheye 14 az file "Session 03.pdf".            
                       end 
                       K=KL+KNL;                                                                                                                   % Az jam'e sakhtie khatti "KL" va sakhtie gheire khati "KNL", matrice sakhtie kol "K" be dast miad. Formulation, safheye 14 az file "Session 03.pdf".                  
                     % Calculate the Nodal Displacements Vectors
                       deltaD(FreeDoFs)=K(FreeDoFs,FreeDoFs)\Pr(FreeDoFs);
                       D=D+deltaD;
                     % Find the New Mechanical and Geometric  Information                                                                                  
                       % The New Position of the Nodes 
                         for j=1:NNodes                                                                                                            % Braye har node, mokhtasate jadide node, az jam'e mokhtsate avalie ba meghdare jabejaiie mohasebe shode baraye node too iteratioe jadid be dsat miad.
                             x(j)=x0(j)+D(2*j-1);                                                                                   %%1            % Derayeye 2j-1 az bordare "D" = jabejaiie node j dar jahate x, az avale kar ta iteratione fe'li.
                             y(j)=y0(j)+D(2*j);                                                                                     %%2            % Derayeye 2j az bordare "D"   = jabejaiie node j dar jahate y, az avale kar ta iteratione fe'li.
                         end
%                            x=x0+D(1:2:((2*NNodes)-1));                                                                            %%1
%                            y=y0+D(2:2:(2*NNodes));                                                                                %%2
                       % The New Stiffness Matrix, DtoV and QtoP Matrix
                         DtoV=zeros(NElements,NDoFs);  
                         QtoP=zeros(NDoFs,NElements);
                         for i=1:NElements                                                                                                                                                                                                           
                             invL(i)=1/((x(End(i))-x(Start(i)))^2+(y(End(i))-y(Start(i)))^2)^0.5;                                                  % Az mokhtasate fe'li estefade mishe.                                                             
                             C(i)=(x(End(i))-x(Start(i)))*invL(i);                                                                                 % Az mokhtasate fe'li estefade mishe.                                                                                             
                             S(i)=(y(End(i))-y(Start(i)))*invL(i);                                                                                 % Az mokhtasate fe'li estefade mishe.                                                                                             
                             a(i,:)=[-C(i) -S(i) C(i) S(i)];                                                                                       % Az mokhtasate fe'li estefade mishe. Formulation, safheye 5 az file "Session 03.pdf".                                                                                                                       
                             DtoV(i,DoFs(i,:))=a(i,:);                                                                                             % Az mokhtasate fe'li estefade mishe. Formulation, safheye 6 az file "Session 03.pdf".     
                             QtoP(DoFs(i,:),i)=transpose(a(i,:));                                                                                  % Az mokhtasate fe'li estefade mishe. Formulation, safheye 10 az file "Session 03.pdf".             
                         end                
                     % Calculate the Elements Excessive Deformations Vectors with the Old DtoV 
                       deltaV=DtoV*deltaD; 
                       V=V+deltaV;                                                                                                                 % Choon too har iteration, zavieye ozv taghiir karde va matrice "DtoV" avaz shode, dige nemishe bordare "V" ro be tore mostaghim ba [V=DtoV*D] be dast avord. Formulation, safheye 6 az file "Session 03.pdf".     
                     % Calculate the Elements Excessive Internal Forces Vector
                       Q=Kel*V;                                                                                                                    % Derayeye i az bordare "Q" = nirooye mehvarii ke az avale kar ta iteratione fe'li be nirooye mehvarie avalie [N0] ezafe shode. Formulation, safheye 8 az file "Session 03.pdf".     
                     % Update the Axial Force of the Elements 
                       N=N0+Q;                                                                                                                     % Ba tavajoh be taghiirate nirooye mehvarie elemana, meghdare kolle nirooye mehvarie elemana [N] update mishe. Formulation, safheye 8 az file "Session 03.pdf".     
                     % Calculate the Internal Nodal Loads Vector with the New QtoP
                       Pint=QtoP*N;                                                                                                 %%3            % Ba tavajoh be meghdare kolle nirooye dakhelie a'za, meghdare niroohaye be vojood oomade too har DoF hesap mishe. Formulation, safheye 10 az file "Session 03.pdf".                                                                                                                 
                       %Pint=Pint+K*deltaD;                                                                                         %%3    
                     % Update the Residual Loads Vector
                       Pr=Pext-Pint;
                     % Collect the Results
                       k=k+1;
                       Reached(s+1)=2*k+1;
                       ResultNodalDisplacements(:,2*k)=D;       ResultNodalDisplacements(:,2*k+1)=D;
                       ResultNodalForces(:,2*k)=Pext;           ResultNodalForces(:,2*k+1)=Pint;
                       ResultKlinear(:,2*k-1)=diag(KL);         ResultKlinear(:,2*k)=diag(KL); 
                       ResultKnonlinear(:,2*k-1)=diag(KNL);     ResultKnonlinear(:,2*k)=diag(KNL);
                       ResultKtotal(:,2*k-1)=diag(K);           ResultKtotal(:,2*k)=diag(K);
                       ResultKForce(:,2*k)=Pint;                ResultKForce(:,2*k+1)=Pint;
                       ResultMemberDeformations(:,k+1)=V; 
                       ResultMemberAxialForces(:,k+1)=Q;
              end
        end
        Duration=toc;
%% 04 - Visualization
close all;
if NoAnimation==1; Pause1=0; Pause2=0; end
ResultKForce(:,2*k+1)=[];%ResultKForce=ResultKForce-ResultKForce(:,1);
load('Example1Results_B_GeomNL_OpenSees')
figure('units','normalized','outerposition',[0 0 1 1]); 
   if DrawStiffnessVariations==0
       % Predefinition
         Y1=YScale*min(min(ResultNodalForces(DoFtoDraw,:)),min(min(OpenSeesCorotationalForces(DoFtoDraw,:)),min(OpenSeesPDeltaForces(DoFtoDraw,:))));                        Y2=YScale*max(max(ResultNodalForces(DoFtoDraw,:)),max(max(OpenSeesCorotationalForces(DoFtoDraw,:)),max(OpenSeesPDeltaForces(DoFtoDraw,:))));
         X1=XScale*min(min(ResultNodalDisplacements(DoFtoDraw,:)),min(min(OpenSeesCorotationalDisplacements(DoFtoDraw,:)),min(OpenSeesPDeltalDisplacements(DoFtoDraw,:))));  X2=XScale*max(max(ResultNodalDisplacements(DoFtoDraw,:)),max(max(OpenSeesCorotationalDisplacements(DoFtoDraw,:)),max(OpenSeesPDeltalDisplacements(DoFtoDraw,:))));
         Y3=ZoomYScale*min(ResultNodalForces(DoFtoDraw,Reached(1):Reached(2)));        Y4=ZoomYScale*max(ResultNodalForces(DoFtoDraw,Reached(1):Reached(2)));
         X3=ZoomXScale*min(ResultNodalDisplacements(DoFtoDraw,Reached(1):Reached(2))); X4=ZoomXScale*max(ResultNodalDisplacements(DoFtoDraw,Reached(1):Reached(2)));
         grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
         xlabel('Displacement [in]','fontsize',13,'fontweight','bold'); ylabel('Load [lb]','fontsize',13,'fontweight','bold');      
         Title=['DoF ' num2str(DoFtoDraw) ' Load-Dispacement Plot']; title(Title,'fontsize',15,'fontweight','bold'); hold on;
         if ShowOnlyFirstStep==1
            if ConnectStartToEnd==0; ylim([Y1, Y2]); xlim([X1, X2]); 
            else;                    ylim([Y3, Y4]); xlim([X3, X4]); end
            % Plotting - Real Behaviour
              P1=plot(OpenSeesCorotationalDisplacements(DoFtoDraw,:),OpenSeesCorotationalForces(DoFtoDraw,:),'b','LineWidth',2);
              P2=plot(OpenSeesPDeltalDisplacements(DoFtoDraw,:),OpenSeesPDeltaForces(DoFtoDraw,:),'c','LineWidth',2);
              if ConnectStartToEnd==0; pause(2*Pause2); end
            % Plotting - Target Load
              P3=plot([0,1*OpenSeesCorotationalForces(DoFtoDraw,size(OpenSeesCorotationalDisplacements,2))],[ResultNodalForces(DoFtoDraw,Reached(1+1)-1),ResultNodalForces(DoFtoDraw,Reached(1+1)-1)],'r','LineWidth',1.5);
              if ConnectStartToEnd==0; pause(ZoomPause); end
            % Zoom
            if ConnectStartToEnd==0
               plot([X3, X4, X4, X3, X3],[Y3, Y3, Y4, Y4, Y3],'color','k','LineWidth',2); hold on;
               pause(2*Pause1);
               for z=1:ZoomSteps
                   ylim([(Y1+((Y3-Y1)*(z/ZoomSteps)^2)), (Y2+((Y4-Y2)*(z/ZoomSteps)^2))]); 
                   xlim([(X1+((X3-X1)*(z/ZoomSteps)^2)), (X2+((X4-X2)*(z/ZoomSteps)^2))]);
                   pause(Pause1);
               end
               pause(Pause2);
            else
               plot([X3, X4, X4, X3, X3],[Y3, Y3, Y4, Y4, Y3],'color','k','LineWidth',2); hold on; 
            end
            if ConnectStartToEnd==0; pause(0.8*Pause2); end
            for i=1:(Reached(2)-1)/2
                % Plotting - Matlab Results
                  scatter(ResultNodalDisplacements(DoFtoDraw,Reached(1)),ResultNodalForces(DoFtoDraw,Reached(1)),'m','filled'); hold on;
                  if i==(Reached(2)-1)/2
                     P4=plot(ResultNodalDisplacements(DoFtoDraw,2*i-1:2*i+1),ResultNodalForces(DoFtoDraw,2*i-1:2*i+1),'g','LineWidth',1.5); hold on;
                     P5=scatter(ResultNodalDisplacements(DoFtoDraw,2*i+1),ResultNodalForces(DoFtoDraw,2*i+1),'m','filled'); hold on;
                  else
                     plot(ResultNodalDisplacements(DoFtoDraw,2*i-1:2*i+1),ResultNodalForces(DoFtoDraw,2*i-1:2*i+1),'g','LineWidth',1.5); hold on;
                     scatter(ResultNodalDisplacements(DoFtoDraw,2*i+1),ResultNodalForces(DoFtoDraw,2*i+1),'m','filled'); hold on;
                     if ConnectStartToEnd==0; pause(0.8*Pause2); end
                  end
            end
            if ConnectStartToEnd==1
               scatter(ResultNodalDisplacements(DoFtoDraw,1:2:Reached(2)),ResultNodalForces(DoFtoDraw,1:2:Reached(2)),'m','filled'); hold on;
               P7=scatter(ResultNodalDisplacements(DoFtoDraw,Reached(2)),ResultNodalForces(DoFtoDraw,Reached(2)),[],[0.4940 0.1840 0.5560],'filled'); hold on;
               pause(0.5*Pause2);
               P6=plot(ResultNodalDisplacements(DoFtoDraw,Reached(1:2)),ResultNodalForces(DoFtoDraw,Reached(1:2)),'color',[0.4660 0.6740 0.1880],'LineWidth',3); hold on;
               scatter(ResultNodalDisplacements(DoFtoDraw,Reached(1:2)),ResultNodalForces(DoFtoDraw,Reached(1:2)),[],[0.4940 0.1840 0.5560],'filled'); hold on;
               pause(0.5*Pause1);
               leg=legend([P1, P2, P3, P4, P5, P7, P6],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation', 'Target Loads', 'MATLAB Iterations in the First Step', 'MATLAB Trial Points in the First Step', 'MATLAB Final Points for the First Step', 'MATLAB Result Curve for the First Step'},'Location','northwest','fontsize',15); 
            else
               pause(0.8*Pause2);
               P7=scatter(ResultNodalDisplacements(DoFtoDraw,Reached(2)),ResultNodalForces(DoFtoDraw,Reached(2)),[],[0.4940 0.1840 0.5560],'filled'); hold on;
               pause(0.5*Pause1);
               leg=legend([P1, P2, P3, P4, P5, P7],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation', 'Target Loads', 'MATLAB Iterations in the First Step', 'MATLAB Trial Points in the First Step', 'MATLAB Final Point for the First Step'},'Location','northwest','fontsize',15); 
            end
         elseif ShowAllSteps==1
            ylim([Y1, Y2]); xlim([X1, X2]);
            % Plotting - Real Behaviour
              P1=plot(OpenSeesCorotationalDisplacements(DoFtoDraw,:),OpenSeesCorotationalForces(DoFtoDraw,:),'b','LineWidth',2);
              P2=plot(OpenSeesPDeltalDisplacements(DoFtoDraw,:),OpenSeesPDeltaForces(DoFtoDraw,:),'c','LineWidth',2);
              if ConnectSteps==0 && NoAnimation==0 && ShowTheJump==0 && DisableTrialandFinalPoints==0; pause(Pause2); end
            % Other
              if ShowTheJump==1; S=Steps:-1:1; else; S=1:Steps; end
              for s=S
                % Plotting - Target Load 
                  if ShowTheJump==0 && DisableTrialandFinalPoints==0
                     if s==1
                        P3=plot([0,1*OpenSeesCorotationalForces(DoFtoDraw,size(OpenSeesCorotationalDisplacements,2))],[ResultNodalForces(DoFtoDraw,Reached(s+1)-1),ResultNodalForces(DoFtoDraw,Reached(s+1)-1)],'r','LineWidth',1.5);
                     else
                        plot([0,1*OpenSeesCorotationalForces(DoFtoDraw,size(OpenSeesCorotationalDisplacements,2))],[ResultNodalForces(DoFtoDraw,Reached(s+1)-1),ResultNodalForces(DoFtoDraw,Reached(s+1)-1)],'r','LineWidth',1.5);
                     end
                     if ConnectSteps==0 && NoAnimation==0; pause(Pause1); end
                  end
                % Plotting - Matlab Results
                  if s==1 && ShowTheJump==0 && DisableTrialandFinalPoints==0
                     scatter(ResultNodalDisplacements(DoFtoDraw,Reached(s)),ResultNodalForces(DoFtoDraw,Reached(s)),'m','filled'); hold on;
                     if ConnectSteps==0 && NoAnimation==0; pause(Pause1); end
                  end
                  if s==Steps
                     P4=plot(ResultNodalDisplacements(DoFtoDraw,Reached(s):Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s):Reached(s+1)),'g','LineWidth',1.5); hold on;
                     if DisableTrialandFinalPoints==0 && ShowTheJump==0
                        P5=scatter(ResultNodalDisplacements(DoFtoDraw,Reached(s)+2:2:Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s)+2:2:Reached(s+1)),'m','filled'); hold on;
                     end
                     if DisableTrialandFinalPoints==0
                        P7=scatter(ResultNodalDisplacements(DoFtoDraw,Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s+1)),[],[0.4940 0.1840 0.5560],'filled'); hold on;
                     end
                  elseif s==1
                     if ShowTheJump==1
                        P8=plot(ResultNodalDisplacements(DoFtoDraw,Reached(s):Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s):Reached(s+1)),'color',[1 0.498 0.153],'LineWidth',2); hold on;  
                     else
                        plot(ResultNodalDisplacements(DoFtoDraw,Reached(s):Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s):Reached(s+1)),'g','LineWidth',1.5); hold on;
                     end
                     if DisableTrialandFinalPoints==0 && ShowTheJump==0 
                        scatter(ResultNodalDisplacements(DoFtoDraw,Reached(s)+2:2:Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s)+2:2:Reached(s+1)),'m','filled'); hold on;
                     end  
                     if DisableTrialandFinalPoints==0
                        scatter(ResultNodalDisplacements(DoFtoDraw,Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s+1)),[],[0.4940 0.1840 0.5560],'filled'); hold on;
                        if ConnectSteps==0 && NoAnimation==0 && ShowTheJump==0 ; pause(Pause1); end
                     end  
                  else
                     plot(ResultNodalDisplacements(DoFtoDraw,Reached(s):Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s):Reached(s+1)),'g','LineWidth',1.5); hold on;
                     if DisableTrialandFinalPoints==0 && ShowTheJump==0
                        scatter(ResultNodalDisplacements(DoFtoDraw,Reached(s)+2:2:Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s)+2:2:Reached(s+1)),'m','filled'); hold on;
                     end  
                     if DisableTrialandFinalPoints==0
                        scatter(ResultNodalDisplacements(DoFtoDraw,Reached(s+1)),ResultNodalForces(DoFtoDraw,Reached(s+1)),[],[0.4940 0.1840 0.5560],'filled'); hold on;
                        if ConnectSteps==0 && NoAnimation==0 && ShowTheJump==0; pause(Pause1); end
                     end
                     
                  end
              end
            if NoAnimation==1 && DisableTrialandFinalPoints==0
               scatter(ResultNodalDisplacements(DoFtoDraw,Reached(2:s+1)),ResultNodalForces(DoFtoDraw,Reached(2:s+1)),[],[0.4940 0.1840 0.5560],'filled'); hold on;
            end
            if ConnectSteps==1
               if DisableTrialandFinalPoints==0; scatter(ResultNodalDisplacements(DoFtoDraw,Reached),ResultNodalForces(DoFtoDraw,Reached),[],[0.4940 0.1840 0.5560],'filled'); hold on; end
               if NoAnimation==0; pause(Pause2); end
               P6=plot(ResultNodalDisplacements(DoFtoDraw,Reached),ResultNodalForces(DoFtoDraw,Reached),'color',[0.4660 0.6740 0.1880],'LineWidth',3); hold on;
               if DisableTrialandFinalPoints==0; scatter(ResultNodalDisplacements(DoFtoDraw,Reached),ResultNodalForces(DoFtoDraw,Reached),[],[0.4940 0.1840 0.5560],'filled'); hold on; end
            end
            pause(0.5*Pause1);
            if ConnectSteps==0 && DisableTrialandFinalPoints==0 && ShowTheJump==0 
               leg=legend([P1, P2, P3, P4, P5, P7],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation', 'Target Loads', 'MATLAB Iterations', 'MATLAB Trial Points', 'MATLAB Final Points'},'Location','northwest','fontsize',15); 
            elseif ConnectSteps==1 && DisableTrialandFinalPoints==0 && ShowTheJump==0 
               leg=legend([P1, P2, P3, P4, P5, P7, P6],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation', 'Target Loads', 'MATLAB Iterations', 'MATLAB Trial Points', 'MATLAB Final Points', 'MATLAB Result Curve'},'Location','northwest','fontsize',15); 
            elseif ConnectSteps==0 && DisableTrialandFinalPoints==1 && ShowTheJump==0 
               leg=legend([P1, P2, P4],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation', 'MATLAB Iterations'},'Location','northwest','fontsize',15); 
            elseif ConnectSteps==0 && DisableTrialandFinalPoints==0 && ShowTheJump==1             
               leg=legend([P1, P2, P4, P7, P8],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation','MATLAB Iterations', 'MATLAB Final Points', 'MATLAB Iterations in the First Step'},'Location','northwest','fontsize',15); 
            elseif ConnectSteps==1 && DisableTrialandFinalPoints==1 && ShowTheJump==0
               leg=legend([P1, P2, P4, P6],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation', 'MATLAB Iterations', 'MATLAB Result Curve'},'Location','northwest','fontsize',15); 
            elseif ConnectSteps==1 && DisableTrialandFinalPoints==0 && ShowTheJump==1                
               leg=legend([P1, P2, P4, P7, P8, P6],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation','MATLAB Iterations', 'MATLAB Final Points', 'MATLAB Iterations in the First Step', 'MATLAB Result Curve'},'Location','northwest','fontsize',15); 
            elseif ConnectSteps==0 && DisableTrialandFinalPoints==1 && ShowTheJump==1                
               leg=legend([P1, P2, P4, P8],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation','MATLAB Iterations', 'MATLAB Iterations in the First Step'},'Location','northwest','fontsize',15); 
            elseif ConnectSteps==1 && DisableTrialandFinalPoints==1 && ShowTheJump==1  
               leg=legend([P1, P2, P4, P8, P6],{'OpenSees Results: Corotational Formulation', 'OpenSees Results: P-Delta Formulation','MATLAB Iterations', 'MATLAB Iterations in the First Step', 'MATLAB Result Curve'},'Location','northwest','fontsize',15); 
            end
         elseif ShowOnlyTheResultCurve==1
            ylim([Y1, Y2]); xlim([X1, X2]);
              P1=plot(OpenSeesCorotationalDisplacements(DoFtoDraw,:),OpenSeesCorotationalForces(DoFtoDraw,:),'b','LineWidth',2);
              P2=plot(OpenSeesPDeltalDisplacements(DoFtoDraw,:),OpenSeesPDeltaForces(DoFtoDraw,:),'c','LineWidth',2);
              P4=plot(ResultNodalDisplacements(DoFtoDraw,Reached),ResultNodalForces(DoFtoDraw,Reached),'color',[0.4660 0.6740 0.1880],'LineWidth',3); hold on;
              leg=legend([P1, P2, P4],{'OpenSees Results: Corotationa Formulation', 'OpenSees Results: P-Delta Formulation', 'MATLAB Result Curve'},'Location','northwest','fontsize',15); 
         end
         Plot_Name=['Geomteric Non-Linear Structure - Full Newton-Raphson Method - ' num2str(Steps) ' Steps - Duration ' num2str(Duration) ' Seconds - ' num2str(k) ' Iterations']; hold on;
         h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;      
   else
       % Predefinition
         ylim([0, YScale*max(max(ResultKlinear(DoFtoDraw,:)),max(max(ResultKnonlinear(DoFtoDraw,:)),max(ResultKtotal(DoFtoDraw,:))))]); 
         xlim([min(ResultNodalForces(DoFtoDraw,:)), XScale*max(ResultNodalForces(DoFtoDraw,:))]);
         grid on; grid minor; ax=gca; ax.GridLineStyle='--'; ax.GridAlpha=0.6; ax.GridColor=['k']; ax.FontSize=12; ax.LineWidth=0.8; ax.TickLength=[0.01 0.01];    
         xlabel('Load [lb]','fontsize',13,'fontweight','bold'); ylabel('Stiffness [lb/in]','fontsize',13,'fontweight','bold');      
         Title=['DoF ' num2str(DoFtoDraw) ' Stiffness-Load Plot']; title(Title,'fontsize',15,'fontweight','bold'); hold on;
         P1=plot(ResultKForce(DoFtoDraw,:),ResultKlinear(DoFtoDraw,:),'r','LineWidth',1.5); hold on;
         P2=plot(ResultKForce(DoFtoDraw,:),ResultKnonlinear(DoFtoDraw,:),'b','LineWidth',1.5); hold on;
         P3=plot(ResultKForce(DoFtoDraw,:),ResultKtotal(DoFtoDraw,:),'g','LineWidth',1.5); hold on;
                     leg=legend([P1, P2, P3],{'Linear Component', 'Non-Linear Component', 'Total Stiffness'},'Location','northwest','fontsize',15); 
         Plot_Name=['Geomteric Non-Linear Structure - Full Newton-Raphson Method - ' num2str(Steps) ' Steps - Duration ' num2str(Duration) ' Seconds - ' num2str(k) ' Iterations']; hold on;
         h1=suptitle(Plot_Name); set(h1,'FontSize',17,'FontWeight','bold'); hold off;        
   end