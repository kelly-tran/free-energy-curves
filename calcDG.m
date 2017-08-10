% =========================================================================
% File: calcDG.m
% Author: Kelly Tran
% Purpose: This MATLAB script constructs a histogram of the potential
%     energy difference and calculates the free energy.
% =========================================================================

clear

% Specify the input file format and set variables
% =========================================================================
formatSpec = repmat('%f',[1 7]);
data_size = [7 Inf];

font_size = 12;
bin_size = 18;

% Read data for reference 2-/2- nuclear configuration
% =========================================================================
fileID = fopen('inp/estottime_22_het1.dat','r');
data22_1 = fscanf(fileID,formatSpec,data_size);
fclose(fileID);

fileID = fopen('inp/estottime_22_het2.dat','r');
data22_2 = fscanf(fileID,formatSpec,data_size);
fclose(fileID);

fileID = fopen('inp/estottime_23_het1.dat','r');
data23_1 = fscanf(fileID,formatSpec,data_size);
fclose(fileID);

fileID = fopen('inp/estottime_23_het2.dat','r');
data23_2 = fscanf(fileID,formatSpec,data_size);
fclose(fileID);

fileID = fopen('inp/estottime_32_het1.dat','r');
data32_1 = fscanf(fileID,formatSpec,data_size);
fclose(fileID);

fileID = fopen('inp/estottime_32_het2.dat','r');
data32_2 = fscanf(fileID,formatSpec,data_size);
fclose(fileID);

% Calculate potential energy difference (DV)
% =========================================================================
% Protein
V_0 = data22_1(3,:)+data22_2(3,:);
V_R = data32_1(3,:)+data32_2(3,:);
V_P = data23_1(3,:)+data23_2(3,:);
DV_0 = V_P - V_R;

% Construct histograms and calculate probability for DV
% =========================================================================
figure(1)
[n_0,x_0] = hist(DV_0,bin_size);
bar(x_0,(n_0/length(DV_0)),'k') 
pbaspect([1 1 1]); box on ; set(gca,'FontSize',font_size,'XDIR','Reverse')
xlabel('\DeltaV (kcal/mol)','FontSize',font_size)
ylabel('Probability','FontSize',font_size)
text(0.1,0.9,'Protein','FontSize',font_size,'Units','normalized')
ylim([0,0.2])

P_0 = n_0 / max(n_0);

% Calculate the free energy
% =========================================================================
G_0 = -(1/1.6774)*log(P_0);

figure(3)
plot(x_0,G_0,'k')
pbaspect([1 1 1]); box on ; set(gca,'FontSize',font_size,'XDIR','Reverse')
xlabel('Reaction Coordinate (kcal/mol)','FontSize',font_size)
ylabel('Free Energy (kcal/mol)','FontSize',font_size)
text(0.1,0.9,'Protein','FontSize',font_size,'Units','normalized')
set(gca,'XLim',[-100 100] , 'XTick',(-100 : 50 : 100))

% Free energy with best fit parabolas
% =========================================================================
[param_0,er_0] = polyfit(x_0,G_0,2)

xi = -120:120;
xmin = -120;
xmax =  120;
ymin = -1;
ymax =  10;

figure(4)
clf
hold on;
plot(xi,polyval(param_0,xi),'k');
plot(x_0,G_0,'k','Marker','o','LineStyle','none','MarkerSize',8)
pbaspect([1 1 1]); box on ; set(gca,'FontSize',font_size,'XDIR','Reverse')
hold off;
xlabel('Reaction Coordinate (kcal/mol)','FontSize',font_size)
ylabel('Free Energy (kcal/mol)','FontSize',font_size)
text(0.1,0.9,'Protein','FontSize',font_size,'Units','normalized')
set(gca,'XLim',[-100 100] , 'XTick',(-100 : 50 : 100))
set(gca,'YLim',[-1 10] , 'YTick',(0 : 2 : 10))

% Write to data files
% =========================================================================
fileID = fopen('inp/Hist_Protein.dat','w');
fprintf(fileID,'       n_0       x_0       P_0\n');
fprintf(fileID,'%10d %10.2f %10.5f \n', ... 
    [n_0(1,:) ; x_0 ; P_0 ]);
status=fclose(fileID);

