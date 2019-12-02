function[SAIDI_Final]=Nominal_MC_Support(gh)
clc
[~,~,CUST]=xlsread( gh,'CONSUMER');
[~,~,SW]=xlsread(gh,'SW');
[~,~,LINE]=xlsread(gh,'LINE');
[~,~,LINETYPE]=xlsread(gh,'LINETYPE');
[~,~,BACKUP]=xlsread(gh,'SOURCENODES');
% Removing header rows
%CUST, SW, LINE, LINETYPE, and BACKUP represent the customer data
%array,switching data array,Line data array, Line characteristics data
%array, and Sourcenode respectively.
CUST = CUST(2:end,:);
SW = SW(2:end,:);
LINE = LINE(2:end,:);
LINETYPE = LINETYPE(2:end,:);
BACKUP=BACKUP(2:end,:)
size(LINE);
% This is the list of all nodes in the network
node_names = unique([LINE(:,[4 5]); SW(:,[2 3])]);

% This is the prototype G matrix, all zeros
G_matrics = zeros(length(node_names));
for line_ind = 1:size(LINE,1)
    % Reading begin and end node names from LINE row
    beg_node = LINE{line_ind,4};
    end_node = LINE{line_ind,5};
    % Determining the indices of said nodes in G
    beg_node_ind = ismember(node_names, beg_node );
    end_node_ind = ismember(node_names, end_node );
    % Modifying G so that line is represented in G (line == 1, in G)
    G_matrics(beg_node_ind,end_node_ind)=1;
    G_matrics(end_node_ind,beg_node_ind)=1;
end

 for sw_ind = 1:size(SW,1)
    % Reading node names and determining indices for switch nodes
    % in a similar manner as was did for lines
    beg_node = SW{sw_ind,2};
    end_node = SW{sw_ind,3};
    beg_node_ind = ismember(node_names, beg_node);
    end_node_ind = ismember(node_names, end_node);
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
 
 G_SW = zeros( length(node_names));
 for line_ind = 1:size(LINE,1)
    % Reading begin and end node names from LINE row
    beg_node = LINE{line_ind,4};
    end_node = LINE{line_ind,5};
    % Determining the indices of said nodes in G
    beg_node_ind = ismember( node_names, beg_node );
    end_node_ind = ismember( node_names, end_node );
    % Modifying G so that line is represented in G (line == 1, in G)
    G_SW(beg_node_ind,end_node_ind)=eps;
    G_SW(end_node_ind,beg_node_ind)=eps;
end

 for sw_ind = 1:size(SW,1)
    % Reading node names and determining indices for switch nodes
    % in a similar manner as was did for lines
    beg_node = SW{sw_ind,2};
    end_node = SW{sw_ind,3};
    beg_node_ind = ismember( node_names, beg_node );
    end_node_ind = ismember( node_names, end_node );
    G_SW(beg_node_ind,end_node_ind)=(SW{sw_ind,4});
    G_SW(end_node_ind,beg_node_ind)=(SW{sw_ind,4});
    
 end
SAIDI=0;
 for line_S=  1:size(LINE,1)
%To show each index number during fault
Fault_line=line_S    
%G_Islanding is modified version of G matrics for using in different necessary condition 
G_Matrics_SW=G_SW;

G_saidi_islanding=G_matrics;
%to find out parmanent faulted nodes, all swithches, breaker should be kept open
G_saidi_islanding(G_saidi_islanding==4) = 0;
G_saidi_islanding(G_saidi_islanding==3) = 0;
G_saidi_islanding(G_saidi_islanding==2) = 0;
%Graph traveres required sparse matrics
G_saidi_islanding = sparse( G_saidi_islanding );
%h = view(biograph(G_saidi_islanding));
%Line_Type, Line_Length, and initial_fault_Node have been find out from
%data XL
Line_Type = LINE{line_S,2};
Line_Length = LINE{line_S,3};
Initial_Fault_Node = LINE{line_S,4};
% Graphtraverse requires starting index to start depth-first search algoritm 
[ ~, Fault_Node_Index] = max( ismember( node_names, Initial_Fault_Node ) );
%Pr_Fault_Node_indexs are parmanent faulted nodes search by Depth-first search algoritm.
Pr_Fault_Node_index = graphtraverse( G_saidi_islanding, Fault_Node_Index);

Pr_Fault_Node_Name = node_names( Pr_Fault_Node_index);
%Pr_Faulted_Cust_number denotes the total number of customar that sense the
%the parmanent fault
Pr_Faulted_Cust_initial=sum( cell2mat(CUST( ...
        ismember( CUST(:,1),Pr_Fault_Node_Name ),2)));
%Lamda and Fault_Repair_time is measured from different type of line 
%for SAIDI calculation issue    
Lamda = Line_Length / 1000 * ...
        LINETYPE{ ismember(LINETYPE(:,1),Line_Type),7};
Fault_Repair_Time= LINETYPE{ ismember(LINETYPE(:,1),Line_Type),8}/60; 
   
Min_SW_time=min(cell2mat(SW(:,4)));             
BESS_Nodes=cell(size(CUST(:,6)));
  for s=1:length(CUST(:,6))
        
    if cell2mat(CUST(s,6))>= Min_SW_time
          
       BESS_Nodes{s}=cell2mat(CUST(s,1));
    end
     
  end
Efective_BESS_Nodes=BESS_Nodes(~cellfun('isempty',BESS_Nodes));

for Bacup_ind=1:size(BACKUP,1)
        
       G_saidi_islanding=G_matrics;
       G_saidi_islanding(G_saidi_islanding==3)= 0;
       G_saidi_islanding = sparse( G_saidi_islanding );
       %h = view(biograph(G_saidi_islanding))
       [ ~, BackupNodeIndex] = max(ismember( node_names,BACKUP(Bacup_ind,1)));
       Node_In_Backup_path = graphtraverse(G_saidi_islanding,BackupNodeIndex);
       NodeName_In_Backup_path=node_names( Node_In_Backup_path);             
       check_fault_zone=ismember(NodeName_In_Backup_path(:,1),Initial_Fault_Node); 
       %Determine the fault zone hence the sourcenode which is responsible
       %for backup connection
       if max(check_fault_zone)==1
           Check_BESS_Availibility=ismember(Efective_BESS_Nodes,NodeName_In_Backup_path(:,1));
           BESS_Nodes_In_Backup_path=Efective_BESS_Nodes(Check_BESS_Availibility(:,1));
           Switch=ismember(SW(:,[2 3]),Pr_Fault_Node_Name);
           Swithch_involved=[SW(Switch(:,1),[2 3]);SW(Switch(:,2),[2 3])];
            for i=1:size(Swithch_involved)     
                  Beg_required_sw=ismember(node_names, Swithch_involved(i,1));
                  End_required_Sw = ismember(node_names, Swithch_involved(i,2));
                  G_saidi_islanding(Beg_required_sw,End_required_Sw)=0;
                  G_saidi_islanding(End_required_Sw,Beg_required_sw)=0;
            end
           G_saidi_islanding = sparse( G_saidi_islanding );
           %h = view(biograph(G_saidi_islanding))
           Short_time_faulted_node = graphtraverse(G_saidi_islanding,BackupNodeIndex);
           Short_time_faulted_node_name= node_names(Short_time_faulted_node);
           Reamaining_open_node=setxor(NodeName_In_Backup_path,Short_time_faulted_node_name);
           %Remaining_open_node represents the rest of the node which are not possible to rennect from sourcenode
           Open_node=setxor( Reamaining_open_node,Pr_Fault_Node_Name);
           %Open_node represents the nodes which are traeted as parmament fault
           BESS_In_OpenNode=ismember(Open_node,BESS_Nodes_In_Backup_path);
           Open_node_withoutBESS=Open_node(BESS_In_OpenNode);
           Check_short_Sw_Time=ismember(SW(:, [2 3]), Short_time_faulted_node_name);
           Total_Path_Cust=sum(cell2mat(CUST(ismember(CUST(:,1),NodeName_In_Backup_path),2)));
           Short_sw_time=min(cell2mat([SW(Check_short_Sw_Time(:,1),4);SW(Check_short_Sw_Time(:,2),4)]))/3600;
           Check_long_Sw_Time=ismember(SW(:, [2 3]), Short_time_faulted_node_name);
           %Finiding the longest switching perio
           Long_sw_time=max(cell2mat([SW(Check_long_Sw_Time(:,1),4);SW(Check_long_Sw_Time(:,2),4)]))/3600;
           %Finid out is there any BESS available in parmanent faulted nodes hence update the parmanent faulted customer
               
            if isempty(Open_node)==0
               Cust_treat_PRF=cell(length(BACKUP),1);
               Manual_SW_T=cell(length(BACKUP),1);
               %Try to reach parmanent fault nodes fron any possible sourcenodes.
               for B=1:size(BACKUP,1)  
                   Gmatrics_w=G_SW ;
                   Switch=ismember(SW(:,[2 3]),Pr_Fault_Node_Name);
                   Swithch_involved=[SW(Switch(:,1),[2 3]);SW(Switch(:,2),[2 3])];
                 for i=1:size(Swithch_involved)     
                     Beg_required_sw=ismember(node_names, Swithch_involved(i,1));
                     End_required_Sw = ismember(node_names, Swithch_involved(i,2));
                     Gmatrics_w(Beg_required_sw,End_required_Sw)=0;
                     Gmatrics_w(End_required_Sw,Beg_required_sw)=0;
                 end
                   Gmatrics_w= sparse(Gmatrics_w); 
                   %h = view(biograph(Gmatrics_w))
                   [ ~, BackupNodeIndex] =max(ismember(node_names,BACKUP(B,1)));
                   [ ~, TargetIndex]=max(ismember(node_names,Open_node(1,1)));
                   [dist, path, pred]= graphshortestpath(Gmatrics_w,BackupNodeIndex,TargetIndex);
                   distf=dist/3600;
                   %In case of possible connection through sourcenode
                      if isinf(distf)==0
                         Manual_SW_T{B}=fix(distf);
                      end 
                   if distf==inf
                     %In case of zero possibility connection through
                     %sourcenode,chk_bess represent the result of BESS availabilty in nodes 
                     chk_bess=ismember( BESS_Nodes_In_Backup_path,Open_node);
                     
                         if max(chk_bess)==1
                            Pr_Faulted_Custi =Pr_Faulted_Cust_initial;
                         else 
                            Pr_Faulted_Custi =(Pr_Faulted_Cust_initial+sum(cell2mat(CUST(ismember(CUST(:,1),Open_node),2))));   
                         end
                   else
                       Pr_Faulted_Custi =Pr_Faulted_Cust_initial ;
                   end 
              Cust_treat_PRF{B}=Pr_Faulted_Custi;
               end 
              %Manual_SW_Tf denotes the minimum switching duration among all N/O switches
              Manual_SW_Tf=min(cell2mat(Manual_SW_T(:,1)));
                 if isempty(Manual_SW_Tf)==1
                  Manual_SW_Tf=Fault_Repair_Time;   
                end
             Pr_Faulted_Cust= min(cell2mat(Cust_treat_PRF(:,1)));
            else
               Manual_SW_Tf=Long_sw_time;
               Pr_Faulted_Cust =Pr_Faulted_Cust_initial;
            end    
           %Update the short time faulted customer and longtime customer
           %based on updated parmanent faulted customer
           if max(Check_BESS_Availibility)==1
              BESS_Supported_customer=cell(size(BESS_Nodes_In_Backup_path(:,1)));
              for b=1:size(BESS_Nodes_In_Backup_path(:,1)) 
                [ ~, BESSIndex]=max( ismember(node_names,BESS_Nodes_In_Backup_path(b,1)));
                %BESSsuportedNode and BESSsuportedNode_name denotes BESS
                %connected node and node name respectively
                BESSsuportedNode=graphtraverse(G_saidi_islanding,BESSIndex);
                BESSsuportedNode_name=node_names(BESSsuportedNode); 
                %To check BESS option in parmanent foulted nodes and short
                %time faulted nodes
                Check_BESSIN_FaultedBrance=ismember(Pr_Fault_Node_Name,BESS_Nodes_In_Backup_path(b,1));
                Check_BESS_In_Shortpath=ismember( Short_time_faulted_node_name,BESS_Nodes_In_Backup_path(b,1));
                  if max(Check_BESSIN_FaultedBrance)==1    
                     BESS_Suported_Cust=0 ;   
                  elseif max(Check_BESS_In_Shortpath)==1   
                     BESS_Suported_Cust=0;                          
                  else  
                     BESS_Suported_Cust=sum(cell2mat(CUST(ismember(CUST(:,1),BESSsuportedNode_name),2)));    
                  end  
                %BESS_Supported_customer denotes the customer are backuped by BESS 
                BESS_Supported_customer{b} =BESS_Suported_Cust;          
             end
             Total_BESS_Supported_customer=max(cell2mat(BESS_Supported_customer));
             shortime_customer=sum(cell2mat(CUST(ismember(CUST(:,1),Short_time_faulted_node_name),2)))+Total_BESS_Supported_customer;
             longtime_customer= Total_Path_Cust-(shortime_customer+Pr_Faulted_Cust);             
           else                               
           shortime_customer=sum(cell2mat(CUST(ismember(CUST(:,1),Short_time_faulted_node_name),2)));
           longtime_customer= Total_Path_Cust-(shortime_customer+Pr_Faulted_Cust);         
           end

      end  

end
Total_Path_Cust;
Pr_Faulted_Cust;
shortime_customer;
longtime_customer;
Manual_SW_Tf;
Total_Cust_Number = sum( cell2mat( CUST(:,2))); 
%Total_Cust_Number represent the total customer connected the system
 if longtime_customer>=0
   if isempty(Manual_SW_Tf)==0
     SAIDI=SAIDI+((Lamda.*Pr_Faulted_Cust.*Fault_Repair_Time)+(Lamda.*longtime_customer.*Manual_SW_Tf)+(Lamda.* shortime_customer.*Short_sw_time))/Total_Cust_Number;
     SAIDI_Final=SAIDI ;
  end
 end
 end
SAIDI_Final;
end