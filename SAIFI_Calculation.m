
function[SAIFI]=SAIFI_Calculation(gh)
clc
%to read the data file 
[~,~,CUST] = xlsread( gh, 'CONSUMER' );
[~,~,SW] = xlsread(gh, 'SW' );
[~,~,LINE] = xlsread( gh, 'LINE' );
[~,~,LINETYPE] = xlsread( gh, 'LINETYPE' );
[~,~,BACKUP] = xlsread( gh, 'SOURCENODES' );

% Removing header rows
CUST = CUST(2:end,:);
SW = SW(2:end,:);
LINE = LINE(2:end,:);
LINETYPE = LINETYPE(2:end,:);
BACKUP=BACKUP(2:end,:);
size(LINE);
% This is the list of all nodes in the network
node_names = unique( [LINE(:,[4 5]); SW(:,[2 3])] );
% This is the prototype G matrix, all zeros
G_matrics = zeros( length( node_names ) );
for line_ind = 1:size(LINE,1)
    % Reading begin and end node names from LINE row
    beg_node = LINE{line_ind,4};
    end_node = LINE{line_ind,5};
    
    % Determining the indices of said nodes in G
    beg_node_ind = ismember( node_names, beg_node );
    end_node_ind = ismember( node_names, end_node );
    
    % Modifying G so that line is represented in G (line == 1, in G)
    G_matrics(beg_node_ind,end_node_ind)=1;
    G_matrics(end_node_ind,beg_node_ind)=1;
end

for sw_ind = 1:size(SW,1)
    
    % Reading node names and determining indices for switch nodes
    % in a similar manner as was did for lines
    beg_node = SW{sw_ind,2};
    end_node = SW{sw_ind,3};
    beg_node_ind = ismember( node_names, beg_node );
    end_node_ind = ismember( node_names, end_node ); 
    % Modifying G so that the switch is represented in G 
    % (Open Sw == 3, Closed Sw == 2, Breaker == 4)
    if SW{sw_ind,5}==1 %Closed Switch
        if SW{sw_ind,4}==0 %Breaker
            G_matrics(beg_node_ind,end_node_ind)=4;
            G_matrics(end_node_ind,beg_node_ind)=4;
        else
            G_matrics(beg_node_ind,end_node_ind)=2;
            G_matrics(end_node_ind,beg_node_ind)=2;
        end
        
      elseif SW{sw_ind,5}==0 %Open Switch
        G_matrics(beg_node_ind,end_node_ind)=3;
        G_matrics(end_node_ind,beg_node_ind)=3;
   end
end
SAIFI=SAIFI_f(G_matrics,CUST,node_names,LINE,LINETYPE);

function [Calculate_SAIFI]=SAIFI_f(G_matrics,CUST,node_names,LINE,LINETYPE)
%% SAIFI Calculation
Calculate_SAIFI = 0;
% We modify G for the purposes for saifi calculation by removing
% open-switches and breakers from G
G_saifi = G_matrics;
G_saifi(G_saifi==4) = 0;
G_saifi(G_saifi==3) = 0;
G_saifi;
G_saifi = sparse( G_saifi );
TotalNCust = sum( cell2mat( CUST(:,2) ) );

for line = 1:size(LINE,1)
    % Reading line specific values for Line Type and Line length
    LineType = LINE{line,2};
    LineLength = LINE{line,3};  
    % Reading Line begin node name
    LineBeginNodeName = LINE{line,4}; 
    % graphtraverse requires the node index of the begin node in
    % 'node_names'. This can be calculated with:
    [ ~, LineBeginNodeIndex] = max( ismember( node_names, LineBeginNodeName ) );    
    % Determining indices of the nodes which experience the fault...
    FaultNodeIndices = graphtraverse( G_saifi, LineBeginNodeIndex );
    %h = view(biograph(G_saifi));
    % ...and extarcting the faulty node names from
    FaultNodeNames = node_names( FaultNodeIndices ); 
    % Find the total number of affected customers by referencing the
    % customers which are part of the previously defined faulty node set
    NFaultyCustomers = sum( cell2mat( CUST( ...
                     ismember(CUST(:,1),FaultNodeNames ),2)));             
   %Find out the effective BESS nodes,so that ineffective BESS nodes are not taken consideration.  
     Min_SW_time=min(cell2mat(SW(:,4)));             
     BESS_Nodes=cell(size(CUST(:,6)));
      for s=1:length(CUST(:,6))
        
          if cell2mat(CUST(s,6))>= Min_SW_time
          
           BESS_Nodes{s}=cell2mat(CUST(s,1));
          end
     
      end
   Efective_BESS_Nodes=BESS_Nodes(~cellfun('isempty',BESS_Nodes));
   BESS_supported_Cust=sum(cell2mat(CUST(ismember(Efective_BESS_Nodes,FaultNodeNames(:,1)),2)));
   %Remaining final intrrupted customer
   NFaultyCustomers =NFaultyCustomers- BESS_supported_Cust;
   % Line faut rate is read by reading the faultrate value from
   % LINETYPE from the row where Linetype mathces the type of this
   % specific line branch and multiplying it with line length.
   % NOTE: specific faultrates for linetypes are per km and line length is
   % in m, hence the division with 1000
   LineFaultRate = LineLength / 1000 * ...
        LINETYPE{ ismember(LINETYPE(:,1),LineType),7};  
   % Incrementing SAIFI:
   Calculate_SAIFI = Calculate_SAIFI + NFaultyCustomers/TotalNCust * LineFaultRate;
    
end

end

end
