% Clean up

close all; clear all; clc;

%% Load in connectivity data
% Created Separate text files for the Target (inj_site) and Source columns
% of the Excel data sheet. Let's load them in. Copied SLNvec and FLNvec
% directly from Excel.

inj_site = dataread('file', '../raw_data/inj_site.txt', '%s', 'delimiter', '\n'); % injection sites (aka Targets)
cell_body_site = dataread('file', '../raw_data/source.txt', '%s', 'delimiter', '\n'); % The list of areas projecting to the Targets/Injection sites
injection_areas_list = unique(inj_site);
supragranular_vec =  dataread('file', '../raw_data/supragranular.txt', '%d', 'delimiter', '\n');
infragranular_vec =  dataread('file', '../raw_data/infragranular.txt', '%d', 'delimiter', '\n');

cell_body_areas_list = unique(cell_body_site);
num_inj_areas = length(injection_areas_list);
num_source_areas = length(cell_body_areas_list);
num_conns = length(inj_site);

%% Create connectivity matrices
%Initialise to make it fill quicker.
supragranular_mat = zeros(length(cell_body_areas_list),length(cell_body_areas_list));
infragranular_mat = zeros(length(cell_body_areas_list),length(cell_body_areas_list));


for i = 1:length(cell_body_areas_list);
   for j = 1:length(cell_body_areas_list);
       sprintf('Injection site %d , Source %d',i,j)
      if strcmp(cell_body_areas_list{i},cell_body_areas_list{j}); % note, as our strings (area names) are of different length we have to do strcmp instead of '=='
          %disp 'within area - filling with zeros'
          supragranular_mat(i,j) = 0;
          infragranular_mat(i,j) = 0;
      else
          %disp('Find the rows in Inj_Site associated with injection site i')
         [rnInj,~]=find(strcmp(cell_body_areas_list{i},inj_site)); % Find the rows in Inj_Site associated with injection site i
         %disp('Find the rows in Source associated with source region j')
         [rnSource,~]=find(strcmp(cell_body_areas_list{j},cell_body_site)); % Find the rows in Source associated with source region j
         %disp('Seeing where they match')
         [match idx] =ismember(rnSource,rnInj);
         if sum(idx) == 0;
             %disp('No match so filling in a zero')
          supragranular_mat(i,j) = 0;
          infragranular_mat(i,j) = 0;

         else
              %disp('Match! Looking up the connection strength in FLNvec')
              % If the same connections were tested across monkeys, take
              % the mean.
         supragranular_mat(i,j) = mean(supragranular_vec(rnInj(nonzeros(idx))));
         infragranular_mat(i,j) = mean(infragranular_vec(rnInj(nonzeros(idx))));
         end
      end
       
       
   end
end

% figure(1)
% imagesc(supragranular_mat)
% figure(2)
% imagesc(infragranular_mat)

%% Construct FLN
all_layer_conns_mat = supragranular_mat + infragranular_mat;
% Test number of non-zero connection
nonzeroconns = nnz(all_layer_conns_mat) 

summed_inputs = sum(all_layer_conns_mat,2); % sum the inputs (total in a row)
summed_inputs_repmat = repmat(summed_inputs,1,num_source_areas);

% Check the number of non-zero rows matches the number of injected areas
nnz(summed_inputs) == num_inj_areas;
injected_areas_index = find(summed_inputs);

FLNmat = zeros(num_source_areas,num_source_areas);
FLNmat(injected_areas_index,:) = all_layer_conns_mat(injected_areas_index,:)./summed_inputs_repmat(injected_areas_index,:);


sum(sum(FLNmat,2)) == num_inj_areas;

%% Construct SLN

SLNmat = supragranular_mat./(supragranular_mat + infragranular_mat);
% replace nans with zeros
SLNmat(isnan(SLNmat)) = 0;
% figure(3)
% imagesc(SLNmat)

%% Create design matrix X


X = zeros(num_conns,num_source_areas);
self_conns = zeros(num_conns,1);

for current_conn = 1:num_conns
    
    if strcmp(cell_body_site(current_conn),inj_site(current_conn))
    
        self_conns(current_conn) = 1;
    else
        
    % Fill in a -1 for the area containing the cell body
        X(current_conn,find(strcmp(cell_body_site(current_conn),cell_body_areas_list))) = -1;

        % Fill in a 1 for the area that was injected (containing the axon
        % terminals)
        X(current_conn,find(strcmp(inj_site(current_conn),cell_body_areas_list))) = 1;

    end
    
end
%%
% Delete rows corresponding to self-connections (e.g. V1 to V1)
X(logical(self_conns),:)=[];
supragranular_vec(logical(self_conns))=[];
infragranular_vec(logical(self_conns))=[];

% Delete zero connections
zero_conns = supragranular_vec + infragranular_vec ==0;

X(zero_conns,:)=[];
supragranular_vec(zero_conns)=[];
infragranular_vec(zero_conns)=[];


%%
% Performing a basic - non-exhaustive check of the legit-ness of the X matrix
if sum(sum(X'))==0
    disp('X matrix seems legit - good stuff')
else 
    disp('!!!!!! Ooops - not all rows in the X matrix add to zero - something is wrong!!!!!!! ')
end

% From Markov 2014b: "The matrix is singular (each row adds to zero); 
% in order to make the model identifiable, one column was therefore deleted, 
% fixing the hierarchical level for the corresponding area at 0. 
% The fitted hierarchy was arbitrarily normalized with respect to V1 in this fashion. 
% Delete V1 column from matrix to make matrix identifiable."
disp('deleting the column corresponding to V1 to make the X matrix non-singular')
X(:,find(strcmp('V1',cell_body_areas_list))) = [];

cell_body_areas_list_no_V1 = cell_body_areas_list;
cell_body_areas_list_no_V1(find(strcmp('V1',cell_body_areas_list))) = [];

% Delete areas with no connections
zero_columns = sum(X)==0;
X(:,zero_columns)=[];
cell_body_areas_list_no_V1_PIR_SUB = cell_body_areas_list_no_V1;
cell_body_areas_list_no_V1_PIR_SUB(zero_columns) = [];


save ../processed_data/beta_binomial_data.mat supragranular_vec infragranular_vec X

%% Calculate beta-binomial fit in R,(hierarchy_regressions.R)

%% annnnnd we're back

% loading in results of beta-binomial regression from R
load ../processed_data/betaBinHierVals.mat
hier_vals_unordered = hier_vals;
clear hier_vals;
% let's put V1 back in just below V2
V1_new_ind = find(strcmp('V2',cell_body_areas_list_no_V1_PIR_SUB));
% Fill V1 back in with value zero 
hier_vals_unordered = [ hier_vals_unordered(1:V1_new_ind-1); 0; hier_vals_unordered(V1_new_ind:end)];
cell_body_areas_list_no_PIR_SUB = [cell_body_areas_list_no_V1_PIR_SUB{1:V1_new_ind-1},{'V1'},cell_body_areas_list_no_V1_PIR_SUB{V1_new_ind:end}];
[hierarchy_vals_fullgraph hier_fullgraph_ind] = sort(hier_vals_unordered);
hierarchy_vals_fullgraph = hierarchy_vals_fullgraph/max(hierarchy_vals_fullgraph);
fullgraph_hierarchical_order = cell_body_areas_list_no_PIR_SUB(hier_fullgraph_ind)


save ../processed_data/hierarchy_89_areas.mat fullgraph_hierarchical_order hierarchy_vals_fullgraph
[~,subg_idx] = ismember(injection_areas_list,cell_body_areas_list_no_PIR_SUB);

hier_vals_subgraph_unordered = hier_vals_unordered(subg_idx);

hier_vals_subgraph_unordered = hier_vals_subgraph_unordered/max(hier_vals_subgraph_unordered);
[hierarchy_vals_subgraph hier_ind] = sort(hier_vals_subgraph_unordered);

subgraph_hierarchical_order = injection_areas_list(hier_ind);

subgraph_hierarchical_order(strmatch('7a',subgraph_hierarchical_order)) = {'7A'};
subgraph_hierarchical_order(strmatch('7b',subgraph_hierarchical_order)) = {'7B'};


%% plot hierarchy
% Beta-binomial model
disp('creating graph relating hierarchical rating to area')
myfig = figure('units','normalized','outerposition',[0 0 0.25 1])
set(gcf,'color','w');
scatter(hierarchy_vals_subgraph,1:num_inj_areas,1000,'.','k');
set(gca,'YTick', 1:num_inj_areas, 'YTickLabel', subgraph_hierarchical_order);
xlabel('Hierarchical Position', 'FontSize', 14)
% set(gca, 'FontSize', 20)
set(gca,'TickLabelInterpreter','none');
set(gca,'fontweight','bold','fontsize',20);

% myfig = figure('units','normalized','outerposition',[0 0 0.5 1])
% set(gcf,'color','w');
% scatter(hierarchy_vals_fullgraph,1:length(hierarchy_vals_fullgraph),200,'.','k');
% set(gca,'YTick', 1:length(hierarchy_vals_fullgraph), 'YTickLabel', fullgraph_hierarchical_order);
% xlabel('Hierarchical Position', 'FontSize', 14)


%% Create a version of the SLN and the FLN ordered by hierarchy

FLNmat_no_PIR_SUB = FLNmat;
FLNmat_no_PIR_SUB([strmatch('PIR',cell_body_areas_list), strmatch('SUB',cell_body_areas_list)],:) = [];
FLNmat_no_PIR_SUB(:,[strmatch('PIR',cell_body_areas_list), strmatch('SUB',cell_body_areas_list)]) = [];

SLNmat_no_PIR_SUB = SLNmat;
SLNmat_no_PIR_SUB([strmatch('PIR',cell_body_areas_list), strmatch('SUB',cell_body_areas_list)],:) = [];
SLNmat_no_PIR_SUB(:,[strmatch('PIR',cell_body_areas_list), strmatch('SUB',cell_body_areas_list)]) = [];

HierOrderedFLN89 = FLNmat_no_PIR_SUB(hier_fullgraph_ind,hier_fullgraph_ind)
HierOrderedSLN89 = SLNmat_no_PIR_SUB(hier_fullgraph_ind,hier_fullgraph_ind)

save ../processed_data/FLN_SLN_89_areas.mat fullgraph_hierarchical_order HierOrderedFLN89 HierOrderedSLN89



UnorderedFLNsubgraph      = FLNmat_no_PIR_SUB(subg_idx,subg_idx);
UnorderedSLNsubgraph      = SLNmat_no_PIR_SUB(subg_idx,subg_idx);

HierOrderedFLNsubgraph = UnorderedFLNsubgraph(hier_ind,hier_ind);
HierOrderedSLNsubgraph = UnorderedSLNsubgraph(hier_ind,hier_ind);


FigHandle = figure('units','normalized','outerposition',[0 0 0.5 0.8]);
set(gcf,'color','w');

imagesc(log10(HierOrderedFLNsubgraph))
colormap 'hot'
caxis([-7, 0])
c1 = colorbar('TickLabels',{'Absent','10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}',''});
xlabel('Projecting from');
ylabel('Projecting to');
set(gca, 'XTickLabel', subgraph_hierarchical_order);
set(gca,'XTick', 1:num_inj_areas, 'XTickLabel', subgraph_hierarchical_order);
set(gca,'YTick', 1:num_inj_areas, 'YTickLabel', subgraph_hierarchical_order);
xtickangle(90)
% c1.Label.String = 'Connection Strength';
% title('FLN')
% set(gca,'FontSize',20)
% c1.Label.FontSize = 30;
set(gca,'TickLabelInterpreter','none');
set(gca,'fontweight','bold','fontsize',20);

%%
FigHandle2 = figure('units','normalized','outerposition',[0.25 0 0.3 0.5]);
set(gcf,'color','w');
imagesc(HierOrderedSLNsubgraph)
colormap 'jet'
caxis([0, 1])
c2 = colorbar
xlabel('Projecting from');
ylabel('Projecting to');
set(gca,'XTick', 1:num_inj_areas, 'XTickLabel', subgraph_hierarchical_order);
set(gca,'YTick', 1:num_inj_areas, 'YTickLabel', subgraph_hierarchical_order);
c2.Label.String = 'Feedback                              Feedforward';
title('SLN')
%% save FLN, SLN, subgraph hierarchy
save ../processed_data/beta_bin_hierarchy_subgraph.mat subgraph_hierarchical_order hierarchy_vals_subgraph HierOrderedFLNsubgraph HierOrderedSLNsubgraph
