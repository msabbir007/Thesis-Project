

function[G_matrics,SAIFI]=NetworkReliability(CUST,SW,LINE,LINETYPE)

[~,~,CUST] = xlsread( 'Lmaki4.xlsx', 'CONSUMER' );
[~,~,SW] = xlsread('Lmaki4.xlsx', 'SW' );
[~,~,LINE] = xlsread( 'Lmaki4.xlsx', 'LINE' );
[~,~,LINETYPE] = xlsread( 'Lmaki4.xlsx', 'LINETYPE' );
[~,~,BACKUP] = xlsread( 'Lmaki4.xlsx', 'SOURCENODES' );
% Removing header rows
CUST = CUST(2:end,:);
SW = SW(2:end,:);
LINE = LINE(2:end,:);
LINETYPE = LINETYPE(2:end,:);
BACKUP=BACKUP(2:end,:)

%% Network matrix, G

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

%Figure G
%figure;

spy( G_matrics );
set(gca,'YTick',1:length(G_matrics),'YTickLabel',node_names );
set(gca,'XTick',1:length(G_matrics),'XTickLabel',node_names );
end
