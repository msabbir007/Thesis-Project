function[SAIDI_Source_Node]=SAIDI_Without_Bess(gh)
clc
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

SAIDI_Source_Node=0;

for line_S=  1:size(LINE,1)
   Fault_line=line_S  
   if line_S<size(LINE,1)
     Chk_target_nod= ismember(SW(:,[2 3]),LINE(line_S+1,4));
     target_nod=[SW( Chk_target_nod(:,1),3);SW( Chk_target_nod(:,2),3)];
     
  elseif line_S>=size(LINE,1)
    Chk_target_nod= ismember(SW(:,[2 3]),LINE(line_S-1,5));
    target_nod=[SW( Chk_target_nod(:,1),2);SW( Chk_target_nod(:,2),2)];
     
  end
  
G_saidi_islanding=G_matrics;
%to find out parmanent faulted nodes, all swithches, breaker should be kept open
G_saidi_islanding(G_saidi_islanding==4) = 0;
G_saidi_islanding(G_saidi_islanding==3) = 0;
G_saidi_islanding(G_saidi_islanding==2) = 0;
%for graph travers
G_saidi_islanding = sparse( G_saidi_islanding );

%h = view(biograph(G_saidi))
LineType = LINE{line_S,2};
%to find line lenth
LineLength = LINE{line_S,3};
LineBeginNodeName = LINE{line_S,4};
% graphtraverse requires the node index of the begin node in
    % 'node_names'. This can be calculated with:
[ ~, LineBeginNodeIndex] = max( ismember( node_names, LineBeginNodeName ) );
PrFaultNodeInd_islandng = graphtraverse( G_saidi_islanding, LineBeginNodeIndex );

%parmanent faulted nodes
PrFaultNodeName_islanding = node_names( PrFaultNodeInd_islandng );
% Find the total number of affected customers by referencing the
 % customers which are part of the previously defined faulty node set

PrFaultedCustomers = sum( cell2mat( CUST( ...
    ismember( CUST(:,1), PrFaultNodeName_islanding ), 2 ) ) );

Lamda = LineLength / 1000 * ...
    LINETYPE{ ismember(LINETYPE(:,1),LineType),7};
%Line fault repair time for specific faulted line can be calculated from
%linetype and if the time is given in secod unit hense division by 3600 for
%making in hour unit
FaultRepairTime= LINETYPE{ ismember(LINETYPE(:,1),LineType),8}/60;
SAIDI_PrFault=Lamda*FaultRepairTime*PrFaultedCustomers;


%Declearing empty cell array to store customer number for each soucenode
%connection.

  G_saidi_islanding=G_matrics;
  %Find switches that are needed to be open for faulty branch isolation
  Switch=ismember(SW(:,[2 3]),PrFaultNodeName_islanding);
  Swithch_involved=[SW(Switch(:,1),[2 3]);SW(Switch(:,2),[2 3])];
  Beg_required_sw=ismember(node_names, Swithch_involved(:,1));
  End_required_Sw = ismember(node_names, Swithch_involved(:,2));

  G_saidi_islanding(Beg_required_sw,End_required_Sw)=0;
  G_saidi_islanding(End_required_Sw,Beg_required_sw)=0;
  G_saidi_islanding(G_saidi_islanding==3) = 0;
  G_saidi_islanding = sparse( G_saidi_islanding );
  

  

  
    for Bacup_ind=1:size(BACKUP,1)
          Bacup_ind;
         %Finiding the connection path by sourcenodes while faulted branch
         %is isolated and open switches are open
         [ ~, BackupNodeIndex] = max( ismember( node_names,BACKUP(Bacup_ind,1) ) );
         [ ~, End_nodeIndex] = max( ismember( node_names, target_nod) );
         [dist, path, pred] = graphshortestpath(G_saidi_islanding,BackupNodeIndex,End_nodeIndex);
         short_path_Source=node_names( path );
    
        
        if isinf(dist)==1
         G_saidi_islanding=G_matrics;
         G_saidi_islanding(G_saidi_islanding==3) = 0;
         G_saidi_islanding = sparse( G_saidi_islanding );
         [disc, pred, closed] = graphtraverse( G_saidi_islanding, BackupNodeIndex );
         Total_NodeName = node_names(closed ) ; 
         Chk_Final_Path=ismember(Total_NodeName,LineBeginNodeName);
         Chk_Sw_Time=ismember(SW(:, [2 3]), Total_NodeName);
         Sw_Max_sec=max(cell2mat([SW(  Chk_Sw_Time(:,1),4);SW( Chk_Sw_Time(:,2),4)]))/3600;
            if max(Chk_Final_Path)==1      
                 G_saidi_islanding=G_matrics; 
                 G_saidi_islanding(Beg_required_sw,End_required_Sw)=0;
                 G_saidi_islanding(End_required_Sw,Beg_required_sw)=0;
                 G_saidi_islanding(G_saidi_islanding==3) = 0;
                 G_saidi_islanding = sparse( G_saidi_islanding );      
                 [disc, pred, closed] = graphtraverse( G_saidi_islanding, BackupNodeIndex );
                 Fast_Node_NodeName = node_names(closed );
                 Fast_SW_Time_Chk=ismember(SW(:,[2 3]),Fast_Node_NodeName);
                 Total_Path_Cust=sum(cell2mat(CUST(ismember(CUST(:,1),  Total_NodeName ),2)));
                 Fast_Node_Cust=sum(cell2mat(CUST(ismember(CUST(:,1),Fast_Node_NodeName ),2))) ;  
                 Fast_Time_to_reach=min(cell2mat([SW( Fast_SW_Time_Chk(:,1),4);SW(Fast_SW_Time_Chk(:,2),4)]))/3600;
                 Slow_SW_Cust=Total_Path_Cust- (Fast_Node_Cust+PrFaultedCustomers);
                 Slow_Time_to_reach=Sw_Max_sec;   

            end 
            
        Slow_Time_to_reach;
        Slow_SW_Cust;
        Fast_Time_to_reach;
        Fast_Node_Cust; 
        
        else
         G_saidi_islanding=G_matrics;
         G_saidi_islanding(G_saidi_islanding==3) = 0;
         G_saidi_islanding = sparse( G_saidi_islanding );
         [disc, pred, closed] = graphtraverse( G_saidi_islanding, BackupNodeIndex );
         Total_NodeName = node_names(closed ); 
         Chk_Final_Path=ismember(Total_NodeName,LineBeginNodeName);
         Chk_Sw_Time=ismember(SW(:, [2 3]), Total_NodeName);
         Sw_Max_sec=max(cell2mat([SW(  Chk_Sw_Time(:,1),4);SW( Chk_Sw_Time(:,2),4)]))/3600;
            if max(Chk_Final_Path)==1
                 G_saidi_islanding=G_matrics; 
                 G_saidi_islanding(Beg_required_sw,End_required_Sw)=0;
                 G_saidi_islanding(End_required_Sw,Beg_required_sw)=0;
                 G_saidi_islanding(G_saidi_islanding==3) = 0;
                 G_saidi_islanding = sparse( G_saidi_islanding );      
                 [disc, pred, closed] = graphtraverse( G_saidi_islanding, BackupNodeIndex );
                 Fast_Node_NodeName = node_names(closed );
                 Fast_SW_Time_Chk=ismember(SW(:,[2 3]),Fast_Node_NodeName);
                 Total_Path_Cust=sum(cell2mat(CUST(ismember(CUST(:,1),  Total_NodeName ),2)));
                 Fast_Node_Cust=sum(cell2mat(CUST(ismember(CUST(:,1),Fast_Node_NodeName ),2)));   
                 Fast_Time_to_reach=min(cell2mat([SW( Fast_SW_Time_Chk(:,1),4);SW(Fast_SW_Time_Chk(:,2),4)]))/3600;
                 Slow_SW_Cust=Total_Path_Cust- (Fast_Node_Cust+PrFaultedCustomers);
                 Slow_Time_to_reach=Sw_Max_sec;
           end     
            
        end  
        Slow_Time_to_reach;
        Slow_SW_Cust;
        Fast_Time_to_reach;
        Fast_Node_Cust;
        
      
    end
       
  
        Slow_Time_to_reach
        Slow_SW_Cust
        Fast_Time_to_reach
        Fast_Node_Cust 
TotalNCust = sum( cell2mat( CUST(:,2) ) );
SAIDI_BESS_SOURCE=(Slow_SW_Cust*Lamda*Slow_Time_to_reach)+(Fast_Node_Cust*Lamda*Fast_Time_to_reach )
SAIDI_Islan_Each_Line=(SAIDI_PrFault+SAIDI_BESS_SOURCE)/TotalNCust
SAIDI_Source_Node=SAIDI_Source_Node+SAIDI_Islan_Each_Line;
end
SAIDI_Source_Node

end