close all;
clear all;

tic;
% rng(1740);
N_SIM = 1;
List_Seed = 1731:(1731+N_SIM-1);


Sim_ID = '2100518_TEST_02'

for n_sim=1:N_SIM;
    n_sim
    rng(List_Seed(n_sim));
% Major parameters
    OCWmin = 7;  % (7 or 15)
    OCWmax = 31;  % (31 or 63 or 255)
    Mode_UORA = 4;
    
    
    % pattern for unassoc. STA
    N_unassoc_In = 5;  % number of unassociated STAs entering to BSS at each period.
    T_unassoc_In = 0;   % period (sec)
    N_unassoc_Out = 2;  % number of unassociated STAs entering to BSS at each period.
    T_unassoc_Out = 4;   % period (sec)
    % N_STAs_existing = 1;
    N_STAs_Assoc = 50; %[10 20 50 100];%[1 5:5:100];
%     N_Hetero_STAs = [10 20 30 40 50];
    
    Opt_Observation_TH = 1;     % set 1 to observe throughput per unit_time
    Opt_Observation_Delay = 1;  % set 1 to observe avg. delay per unit_time
    Opt_Save_Per_STA_TH = 1;
    
    if Opt_Observation_TH
       Observation_Period = 100000; %(slots), (note) 10,000 slot = 90 msec 
       log_slot_time_SUC = 0;
       log_SUC = 0;
       Next_Observation_Time_TH = Observation_Period;
    end

    if Opt_Observation_Delay
       Observation_Period_Delay = 1000; %(slots), (note) 10,000 slot = 90 msec 
       log_slot_time_SUC = 0;
       log_SUC = 0;
       Next_Observation_Time_Delay = Observation_Period_Delay;
    end

    if Opt_Save_Per_STA_TH
       Result_Per_TH = zeros(max(N_STAs_Assoc), length(N_STAs_Assoc))+NaN; 
    end

    % PHY parameters
    BW = 20; % (MHz)

    N_ss = 1; % number of spatial stream
    N_bp = [1 2 2 4 4 6 6 6 8 8 10 10]; % Number of Coded bit per subcarrier (depends on modulation order)
    R = [1/2 1/2 3/4 1/2 3/4 2/3 3/4 5/6 3/4 5/6 3/4 5/6]; % coding rate (depends on modulation order)
    T_dft = 12.8*10^-6;  % OFDM Symbol Duration (12.8 usec)
    T_gi = 1.6*10^-6;   % Guard Interval (0.8 usec)

    N_RU = 9;
    N_RA_RU_0 = 8;      % # of RA-RUs for associated STAs
    N_RA_RU_2045 = N_RU-N_RA_RU_0;   % # of RA-RUs for unassociated STAs
    Opt_Random_N_RA_RU_0 = 0;  % set to 1, to enalble random number of RA-RUs(AID=0) 
    if Opt_Random_N_RA_RU_0
        N_Random_RA_RU_0 = ceil(rand*N_RA_RU_0);
    end
    
    
    N_slot = 1111112*6;  % Simulation time in slot
    % N_Trigger = 10000; % number of the trigger frames

    for idx_Sim = 1:length(N_STAs_Assoc)

    idx_Sim
    if N_RA_RU_0 + N_RA_RU_2045 ~= N_RU
        error('Wrong RA-RU configuration!');
    else
            % set AID for each RU
        for i = 1:N_RA_RU_0
            RU.AID(i) = 0;
        end

        for i = N_RA_RU_0+1:N_RU
            RU.AID(i) = 2045;
        end

    end

        if N_RU == 9
            N_sd = [24 24 24 24 24 24 24 24 24]; % 9-RU: all RU has 26 tone
        elseif N_RU == 8
            N_sd = [24 24 24 24 24 24 24 48];% 8-RU: 7 RU has 26 tone, 1 RU has 56 tone
        elseif N_RU == 7
            N_sd = [24 24 24 24 24 48 48]; % 7-RU: 5 RU has 26 tone, 2 RU has 56 tone
        elseif N_RU == 6
            N_sd = [24 24 24 48 48 48]; % 6-RU: 3 RU has 26 tone, 3 RU has 56 tone
        elseif N_RU == 5
            N_sd = [24 48 48 48 48]; % 5-RU: 1 RU has 26 tone, 4 RU has 56 tone
        elseif N_RU == 4
            N_sd = [24 24 24 24]; %[48 48 48 48];
        elseif N_RU == 2
            N_sd = [102 102];
        elseif N_RU == 1
            N_sd = [234];
        else
            error('the value of N_RU is not supported yet!');
        end

    if N_RU ~= length(N_sd)
        error('configuration is ivalid!');
    end



    % MAC parameters
    % Set_N_STA = [1 5:5:100];
    % R_STA_Control = [0.2 0.8];
    % N_STA_Control1 = 0;

    N_STAs_existing = N_STAs_Assoc(idx_Sim);
    N_STA = N_STAs_existing + N_unassoc_In;  % initial number of STAs(unassoc.)

    %     N_STA_Control1 = ceil(R_STA_Control(1)*N_STA); 
    %     N_STA_Control2 = floor(R_STA_Control(2)*N_STA);
    % N_STA_Control1 = 0; 
    % N_STA_Control2 = N_STA-N_STA_Control1;
    %     N_STA_Control2 = 10; 
    %     N_STA_Control1 = N_STA-N_STA_Control2;   

    % if N_STA_Control1 + N_STA_Control2 ~= N_STA
    %     error('unexpected N_STA control!');
    % end
    % STA.OBO_Type(1:N_STA_Control1) = 1; % Non BEB STAs
    % STA.OBO_Type((N_STA_Control1+1):(N_STA_Control1+N_STA_Control2)) = 2;  % MIAD+OBO control STAs

    CWmin = OCWmin;
    CWmax = OCWmax;  %(31 or 63)
    Common_MCS = 5;  % common MCS order of STAs for data transmission.
    MCS = Common_MCS*ones(1,N_STA);% 0~11
    Basic_Data_Rate = (N_sd.*N_bp(0+1).*R(0+1)*N_ss)/(T_dft+T_gi);

    % Data_Rate = (N_sd.*N_bp(MCS+1).*R(MCS+1)*N_ss)/(T_dft+T_gi);% Data_Rate per RU
    U_slot = 9*10^-6; 


    SIFS = 16*10^-6;

    % Traffic
    L_MPDU = 2000*ones(1,N_STA); % (Byte)
    L_PHY = 40*10^-6;      % preamble length (40 usec)
    L_Trigger = 100*10^-6; % length of trigger frame
    L_BACK = 68*10^-6;     % length of block ack

    L_Assoc_Req = 38; % (Byte) Length of Association request frame 
    L_Assoc_Req_in_Slot = ceil( (L_PHY + 8*L_Assoc_Req/Basic_Data_Rate(1))/U_slot ); % in slot

    %L_PHY = ceil(L_PHY/U_slot);      % preamble length (40 usec)
    %L_Trigger = ceil((L_PHY+L_Trigger)/U_slot); % 
    %L_BACK = ceil((L_PHY+L_BACK)/U_slot);

    % Mode_UORA = 0;
    % 0: ax UORA, RU is randomly chosen.
    % 1: My UORA, RU is chosen based on the size of OBO.
    % 2: OBO decrease by alpha*RU, alpha is the value more than 0.
    % 3: OBO decrease by alpha*RU, alpha is adapted based on collision event (Centralized-Way).
    % 4: OBO decrease by alpha*RU, alpha is adapted based on collision event (Distributed-Way).
    % 5: [Prof.] MIMD
    % 6: go to init. alpha value when ?
    % 7: AIMD (combination 4 and 5)
    % 8: MIAD (another combination of 4 and 5)
    % 9: (test) use optimal OCW value according to nbr of contending STAs

    if Mode_UORA == 2
       alpha = 1.5;

    elseif Mode_UORA == 3
    %        Init_alpha = 1.0;
       alpha = 1.0; 
       Max_alpha = 2.0;
       Min_alpha = 0.1;
       Inc = 0.1;
       Dec = 0.1;
    elseif Mode_UORA == 4 || Mode_UORA == 5 || Mode_UORA == 6 || Mode_UORA == 7 || Mode_UORA == 8
       Init_alpha = 1.0;
       alpha = Init_alpha*ones(N_STA, 1); 
       Max_alpha = 2.0; %Inf
       Min_alpha = 0.1;
       Inc = 0.1%N_RA_RU_0*0.01;
       Dec = 0.1%N_RA_RU_0*0.01;
        
       On_Selfish_STA = 0;
       On_UORA_STD_STA = 0;
       
       if (On_Selfish_STA & On_UORA_STD_STA) == 1
           error('These options cannot coexist!');
       end
       
       if On_Selfish_STA   % the mode for the secenario where the selfish STAs exist.
          
           N_Selfish_STA = 40;
           Init_alpha_selfish = 2.0;
           alpha(1:N_Selfish_STA) = Init_alpha_selfish; % rule 1) the STAs from 1 to 'N_Selfish_STA' are defined as the selfish STA.
           
           if N_Selfish_STA > N_STA
              error('The number of selfish STAs cannot exceed the total number of STAs.');
           end
           
       elseif On_UORA_STD_STA
           N_UORA_STD_STA = 0;
           if N_UORA_STD_STA > N_STA
              error('The number of UORA standard STAs cannot exceed the total number of STAs.');
           end           
       end
           
       if Mode_UORA == 5 || Mode_UORA == 6 || Mode_UORA == 7 || Mode_UORA == 8
          Succ = zeros(N_STA, 2);  % 1st col: successive # of success
                                   % 2nd col: successive # of collisions 
    %        elseif Mode_UORA == 6
    %           
       end

       Opt_Observation_Alpha = 1;
       if Opt_Observation_Alpha 
          Target_STA = 1;
          Target_N_STAs = [10 20 100];
          log_alpha = [];
          log_slot_time = [];
       end

    elseif Mode_UORA == 9
        Tab_OCW = load('Tab_OCW.in');

    end

    Mode_OBO_BEB = 1; % 0: disable BEB, 1: enable BEB
    Mode_OBO_XIXD = 0; % 0: disable
    if Mode_OBO_XIXD
       Opt_Inc ='A';  % A: additive Increment, M: Multiple Increment
       Opt_Inc_val = CWmin*0.5; % N_RU or N_STA
       Opt_Dec ='M';  % A: additive Decrement, M: Multiple Decrement
       Opt_Dec_val = 2;
    end
    if Mode_OBO_BEB && Mode_OBO_XIXD > 0
        error('Mode is collided!');
    end


    % transformed to slot.

    T_unassoc_In = ceil(T_unassoc_In/U_slot);
    T_unassoc_Out = ceil(T_unassoc_Out/U_slot);
    T_Period_InOut = (T_unassoc_In + T_unassoc_Out);

    Next_T_unassoc_Out = T_unassoc_Out;  % next time when unassoc. STAs enter to BSS
    Next_T_unassoc_In = Next_T_unassoc_Out + T_unassoc_In;  % next time when unassoc. STAs enter to BSS



    % T_unassoc_Out = ceil(T_unassoc_Out/U_slot);
    % Next_T_unassoc_Out = T_unassoc_Out;  % next time when unassoc. STAs enter to BSS
    % 
    % T_unassoc_In = ceil(T_unassoc_In/U_slot);
    % Next_T_unassoc_In = Next_T_unassoc_Out+T_unassoc_In;  % next time when unassoc. STAs enter to BSS


    Sim_Rst = zeros(N_STA, 4); 
    % 1st col: # of Access
    % 2nd col: # of Success
    % 3rd col: # of Collisions
    % 4th col: (delay) cummulative number of trigger frames issued until transmission success

    Sim_Rst_RU = zeros(N_RU, 3); 
    % 1st col: # of IDLE
    % 2nd col: # of Success
    % 3rd col: # of collisions
    N_TF = 0; % number of trigger frames that the AP issued
    n_tf = zeros(N_STA, 1);

    Observation_N_assoc = zeros(1,3);
    % 1st col: time
    % 2nd col: # of N_assoc
    % 3rd col: # of N_data


    % begin simulation; Initialization
    STA.CWmin = CWmin*ones(1, N_STA);
    STA.CWmax = CWmax*ones(1, N_STA);
    STA.CW = STA.CWmin;
    if Mode_UORA == 9
        if Opt_Random_N_RA_RU_0
             STA.CW(1:N_STAs_existing) = Tab_OCW(knnsearch(Tab_OCW(:,1), N_STAs_existing), N_Random_RA_RU_0+1)*ones(1, N_STAs_existing);
        else
            STA.CW(1:N_STAs_existing) = Tab_OCW(knnsearch(Tab_OCW(:,1), N_STAs_existing), N_RA_RU_0+1)*ones(1, N_STAs_existing);
        end
    %     STA.CW = Tab_OCW(knnsearch(Tab_OCW(:,1), N_STA), 2)*ones(1, N_STA);
    end
    STA.Type = zeros(N_STA, 1); % 0: unassociated, 1: associated, -1:disabled
    STA.Type(1:N_STAs_existing) = ones(N_STAs_existing, 1); % 0: unassociated, 1: associated
    STA.Type(N_STAs_existing+1:N_STA) = zeros(N_unassoc_In, 1); % 0: unassociated, 1: associated

    STA.Assoc_Delay = zeros(N_STA, 3); % (in slot)
    STA.Assoc_Delay(1:N_STAs_existing,:) = zeros(N_STAs_existing, 3) + NaN; % (in slot)
    STA.Assoc_Delay(N_STAs_existing+1:N_STA,:) = zeros(N_unassoc_In, 3); % (in slot)
        % 1st col: entering time to BSS
        % 2nd col: Time when receiving Association Response(i.e., ACK)
        % 3rd col: association delay (2nd - 1st)

    STA.OBO = floor(STA.CW.*rand(1, N_STA)); % Initialize OBO,
    STA.Assigned_RU = zeros(N_STA, 1); % 1st col: RU #, 

    log_OCW_assoc = OCWmin*ones(1,2);
    log_OCW_unassoc = OCWmin*ones(1,2);

    i = 1; % ceil((L_PHY + L_Trigger + SIFS)/U_slot);

    while i <= N_slot

        if Next_T_unassoc_In <= i && N_unassoc_In > 0
            STA_ID_Disable = find(STA.Type == -1);
            if length(STA_ID_Disable) == N_unassoc_In
                STA.Type(STA_ID_Disable) = zeros(length(STA_ID_Disable),1);
                Next_T_unassoc_In = i + T_Period_InOut;
            else
                N_STA = N_STA + N_unassoc_In;
                Next_T_unassoc_In = i + T_Period_InOut;

                STA.Type = [STA.Type; zeros(N_unassoc_In, 1)];

                L_MPDU = [L_MPDU 2000*ones(1,N_unassoc_In)]; % (Byte)
                MCS = [MCS Common_MCS*ones(1,N_unassoc_In)];% 0~11
                STA.CWmin =[STA.CWmin CWmin*ones(1, N_unassoc_In)];
                STA.CWmax =[STA.CWmax CWmax*ones(1, N_unassoc_In)];
                STA.CW = [STA.CW CWmin*ones(1, N_unassoc_In)];
                STA.OBO =[STA.OBO floor(STA.CW((N_STA-N_unassoc_In+1):N_STA).*rand(1, N_unassoc_In))]; % Initialize OBO,
                STA.Assigned_RU = [STA.Assigned_RU; zeros(N_unassoc_In, 1)]; % 1st col: RU #, 

                STA.Assoc_Delay = [STA.Assoc_Delay; [ (Next_T_unassoc_In-T_unassoc_In)*ones(N_unassoc_In,1) zeros(N_unassoc_In, 2) ] ];
                Sim_Rst = [Sim_Rst; zeros(N_unassoc_In, 4)];
                n_tf = [n_tf; zeros(N_unassoc_In, 1)];

                if Mode_UORA > 1 && Mode_UORA ~= 9
                   alpha = [alpha; Init_alpha*ones(N_unassoc_In, 1)];
                   if Mode_UORA == 5 || Mode_UORA == 6 || Mode_UORA == 7 || Mode_UORA == 8
                      Succ = [ Succ; zeros(N_unassoc_In, 2)];  % 1st col: successive # of success  
                   end
                end
            end
        end

        if Next_T_unassoc_Out <= i && N_unassoc_Out > 0
           STA_ID_Assoc   = find(STA.Type == 1);
           Disable_STA_ID = [];
           while length(Disable_STA_ID) < N_unassoc_Out
                Temp_Disable_STA_ID = STA_ID_Assoc(ceil(rand*length(STA_ID_Assoc)));
                if ismember(Temp_Disable_STA_ID, Disable_STA_ID)
                else
                    Disable_STA_ID = [Disable_STA_ID Temp_Disable_STA_ID];
                    STA.Type(Temp_Disable_STA_ID) = -1;
                end
           end
           Next_T_unassoc_Out = i + T_Period_InOut;
        end

        STA_ID_Unassoc = find(STA.Type == 0);
        STA_ID_Assoc   = find(STA.Type == 1);

        if Mode_UORA == 0 || Mode_UORA == 9  % ax UORA: decrease OBO by the number of RU
            if Opt_Random_N_RA_RU_0
                STA.OBO(STA_ID_Unassoc) = STA.OBO(STA_ID_Unassoc) - N_RA_RU_2045; 
                STA.OBO(STA_ID_Assoc) = STA.OBO(STA_ID_Assoc) - N_Random_RA_RU_0;             
            else
                STA.OBO(STA_ID_Unassoc) = STA.OBO(STA_ID_Unassoc) - N_RA_RU_2045; 
                STA.OBO(STA_ID_Assoc) = STA.OBO(STA_ID_Assoc) - N_RA_RU_0; 
            end
        elseif Mode_UORA == 1 % My UORA
            STA_ID_TX = find(STA.OBO <= N_RU);
            STA.Assigned_RU(STA_ID_TX) = STA.OBO(STA_ID_TX);
            STA.OBO = STA.OBO - N_RU;
            error('not supported');
        elseif Mode_UORA == 2
            STA.OBO = STA.OBO - alpha*N_RU; 
            STA_ID_TX = find(STA.OBO <= 0);
            STA.Assigned_RU(STA_ID_TX) = ceil(N_RU*rand(length(STA_ID_TX),1));
            error('not supported');

        elseif Mode_UORA == 3
            STA.OBO = STA.OBO - alpha*N_RU; 
            STA_ID_TX = find(STA.OBO <= 0);
            STA.Assigned_RU(STA_ID_TX) = ceil(N_RU*rand(length(STA_ID_TX),1));
            error('not supported');

        elseif Mode_UORA == 4 || Mode_UORA == 5 || Mode_UORA == 6 || Mode_UORA == 7 || Mode_UORA == 8
            
            if On_Selfish_STA
                if ~isempty(STA_ID_Unassoc)
                    error('An unassociated STA is not considered in the selfish STA mode!');
                end
                
                if ~isempty(STA_ID_Assoc)
                    Idx_Selfish_STA = find(STA_ID_Assoc <= N_Selfish_STA);  % 
                    Idx_Non_Selfish_STA = STA_ID_Assoc; 
                    Idx_Non_Selfish_STA(Idx_Non_Selfish_STA <= N_Selfish_STA) = [];
                    
                    if ~isempty(Idx_Selfish_STA)
                        STA.OBO(Idx_Selfish_STA) = STA.OBO(Idx_Selfish_STA) - (alpha(Idx_Selfish_STA))'*N_RA_RU_0; 
                    end
                    
                    if ~isempty(Idx_Non_Selfish_STA)
                        STA.OBO(Idx_Non_Selfish_STA) = STA.OBO(Idx_Non_Selfish_STA) - (alpha(Idx_Non_Selfish_STA))'*N_RA_RU_0; 
                    end
                end
                
            elseif On_UORA_STD_STA
                
                if ~isempty(STA_ID_Unassoc)
                    error('An unassociated STA is not considered in the UORA standard STA mode!');
                end
                    Idx_UORA_STD_STA = find(STA_ID_Assoc <= N_UORA_STD_STA);  % 
                    Idx_Non_UORA_STD_STA = STA_ID_Assoc; 
                    Idx_Non_UORA_STD_STA(Idx_Non_UORA_STD_STA <= N_UORA_STD_STA) = [];
                    
                    if ~isempty(Idx_UORA_STD_STA)
                        STA.OBO(Idx_UORA_STD_STA) = STA.OBO(Idx_UORA_STD_STA) - N_RA_RU_0; 
                    end
                    
                    if ~isempty(Idx_Non_UORA_STD_STA)
                        STA.OBO(Idx_Non_UORA_STD_STA) = STA.OBO(Idx_Non_UORA_STD_STA) - (alpha(Idx_Non_UORA_STD_STA))'*N_RA_RU_0; 
                    end                
                
            else
                if ~isempty(STA_ID_Unassoc)
                    STA.OBO(STA_ID_Unassoc) = STA.OBO(STA_ID_Unassoc) - (alpha(STA_ID_Unassoc))'*N_RA_RU_2045; 
                end
                if ~isempty(STA_ID_Assoc)
                    if Opt_Random_N_RA_RU_0
                        STA.OBO(STA_ID_Assoc) = STA.OBO(STA_ID_Assoc) - (alpha(STA_ID_Assoc))'*N_Random_RA_RU_0; 
                    else
                        STA.OBO(STA_ID_Assoc) = STA.OBO(STA_ID_Assoc) - (alpha(STA_ID_Assoc))'*N_RA_RU_0; 
                    end
                end
                %         STA.OBO = STA.OBO - (alpha)'.*N_RU; 
        %         STA_ID_TX = find(STA.OBO <= 0);
        %         STA.Assigned_RU(STA_ID_TX) = ceil(N_RU*rand(length(STA_ID_TX),1));

                if Opt_Observation_Alpha && ismember(N_STA, Target_N_STAs);
                    log_slot_time = [log_slot_time i];
                    log_alpha = [log_alpha alpha(Target_STA)]; 
                end
            end

        end

        STA_ID_TX = find(STA.OBO <= 0);
        STA_ID_TX_Unassoc = setdiff(STA_ID_TX, STA_ID_Assoc);
        STA_ID_TX_Assoc   = setdiff(STA_ID_TX, STA_ID_Unassoc);

        if Opt_Random_N_RA_RU_0
            STA.Assigned_RU(STA_ID_TX_Assoc) = ceil(N_Random_RA_RU_0*rand(length(STA_ID_TX_Assoc),1));
        else
            STA.Assigned_RU(STA_ID_TX_Assoc) = ceil(N_RA_RU_0*rand(length(STA_ID_TX_Assoc),1));
        end
        STA.Assigned_RU(STA_ID_TX_Unassoc) = ceil(N_RA_RU_2045*rand(length(STA_ID_TX_Unassoc),1))+N_RA_RU_0;

        N_TF = N_TF + 1;
        n_tf = n_tf + 1;

    %     STA.OBO(STA_ID_TX) = floor(STA.CW*rand(length(STA_ID_TX), 1));

        STA_ID_SUC = [];
        STA_ID_COL = [];
        for j = 1:N_RU
            STA_ID_per_RU = find(STA.Assigned_RU == j);

            if isempty(STA_ID_per_RU) 
                Sim_Rst_RU(j, 1) = Sim_Rst_RU(j, 1) + 1; % count up # IDLE
            else
                Sim_Rst(STA_ID_per_RU, 1) = Sim_Rst(STA_ID_per_RU, 1) + 1; % count up # access

                % When collided Transmission 
                if length(STA_ID_per_RU) >1  % When collision occurs,
                    Sim_Rst(STA_ID_per_RU, 3) = Sim_Rst(STA_ID_per_RU, 3) + 1; % count up # collisions
                    Sim_Rst_RU(j, 3) = Sim_Rst_RU(j, 3) + 1;
                    STA_ID_COL = [STA_ID_COL; STA_ID_per_RU];

                    if Mode_OBO_BEB
                        if Mode_UORA == 9
                            if  RU.AID(j) == 0 && (N_unassoc_In > 0 || N_unassoc_Out > 0)
                                if Opt_Random_N_RA_RU_0
                                    STA.CW(STA_ID_per_RU) = Tab_OCW( knnsearch(Tab_OCW(:,1), length(STA_ID_Assoc)), N_Random_RA_RU_0+1)*ones(1, length(STA_ID_per_RU) );
                                else
                                    STA.CW(STA_ID_per_RU) = Tab_OCW( knnsearch(Tab_OCW(:,1), length(STA_ID_Assoc)), N_RA_RU_0+1)*ones(1, length(STA_ID_per_RU) );
                                end
                                if ismember(STA_ID_per_RU, 1)
                                    log_OCW_assoc = [log_OCW_assoc; [i STA.CW(STA_ID_per_RU(1))] ];
                                end
                                
                            elseif RU.AID(j) == 2045 && (N_unassoc_In > 0 || N_unassoc_Out > 0)
                                STA.CW(STA_ID_per_RU) = Tab_OCW( knnsearch(Tab_OCW(:,1), length(STA_ID_Unassoc)), N_RA_RU_2045+1)*ones(1, length(STA_ID_per_RU) );
                                log_OCW_unassoc = [log_OCW_unassoc; [i STA.CW(STA_ID_per_RU(1))] ];
                            end

                            % unassociated STA과 associated STA의 CW 갱신의 구분
                        else
                            STA.CW(STA_ID_per_RU) = (STA.CW(STA_ID_per_RU)+1)*2 - 1;
                            STA_ID_CWmax = find(STA.CW(STA_ID_per_RU) >= STA.CWmax(STA_ID_per_RU));
                            STA.CW(STA_ID_per_RU(STA_ID_CWmax)) = STA.CWmax(STA_ID_per_RU(STA_ID_CWmax));
                        end
                    elseif Mode_OBO_XIXD
                        if strcmp(Opt_Inc,'A')
                           STA.CW(STA_ID_per_RU) = round(STA.CW(STA_ID_per_RU) + Opt_Inc_val); 
                        elseif strcmp(Opt_Inc, 'M')
                           STA.CW(STA_ID_per_RU) = round(STA.CW(STA_ID_per_RU) * Opt_Inc_val);   
                        else
                            error('unexpected XIXD mode!');
                        end
                        if Mode_UORA == 9
                            % nothing to do(igonore OCWmin/max value)
                        else
                            STA_ID_CWmax = find(STA.CW(STA_ID_per_RU) >= STA.CWmax(STA_ID_per_RU));
                            STA.CW(STA_ID_per_RU(STA_ID_CWmax)) = STA.CWmax(STA_ID_per_RU(STA_ID_CWmax));
                        end
                    end

                    if Mode_UORA == 3
                        alpha = alpha - Dec;
                        alpha = max([alpha, Min_alpha]);
                    
                    elseif Mode_UORA == 4
                        if On_Selfish_STA
                            Idx_Selfish_STA = find(STA_ID_per_RU <= N_Selfish_STA);  % 
                            Idx_Non_Selfish_STA = STA_ID_per_RU; 
                            Idx_Non_Selfish_STA(Idx_Non_Selfish_STA <= N_Selfish_STA) = [];
                        
                            if ~isempty(Idx_Selfish_STA)
                                % alpha(Idx_Selfish_STA) = alpha(Idx_Selfish_STA)
                            end
                    
                            if ~isempty(Idx_Non_Selfish_STA)
                                alpha(Idx_Non_Selfish_STA) = alpha(Idx_Non_Selfish_STA) - Dec;
                                alpha(alpha < Min_alpha) = Min_alpha; % note) non-selfish STA만 alpha가 감소하므로.
                            end    
                            
                        elseif On_UORA_STD_STA
                            Idx_UORA_STD_STA = find(STA_ID_per_RU <= N_UORA_STD_STA);  % 
                            Idx_Non_UORA_STD_STA = STA_ID_per_RU; 
                            Idx_Non_UORA_STD_STA(Idx_Non_UORA_STD_STA <= N_UORA_STD_STA) = [];

                            if ~isempty(Idx_UORA_STD_STA)   % followw the Basic BEB operation
                               % do nothing because BEB already has done above. 
                            end

                            if ~isempty(Idx_Non_UORA_STD_STA)
                                alpha(Idx_Non_UORA_STD_STA) = alpha(Idx_Non_UORA_STD_STA) - Dec;
                                alpha(alpha < Min_alpha) = Min_alpha; % note) non-UORA STA만 alpha가 감소하므로, 구분하지 않음. 
                            end                              
                        else
                            alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Dec;
                            alpha(alpha < Min_alpha) = Min_alpha;
                        end
                        
                    elseif Mode_UORA == 5
                        Succ(STA_ID_per_RU, 2) = Succ(STA_ID_per_RU, 2) + 1;
                        Succ(STA_ID_per_RU, 1) = 0;
                        %alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Dec*log2(STA.CW(STA_ID_per_RU)'/CWmin);
                        alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Dec*2.^(Succ(STA_ID_per_RU, 2)-1);
                        alpha(alpha < Min_alpha) = Min_alpha;

                    elseif Mode_UORA == 6 
                        for idx = 1:length(STA_ID_per_RU)
                            Succ(STA_ID_per_RU(idx), 2) = Succ(STA_ID_per_RU(idx), 2) + 1;
                            if Succ(STA_ID_per_RU(idx), 1) > 0
                               alpha(STA_ID_per_RU(idx)) = Init_alpha; 
                            else
                               alpha(STA_ID_per_RU(idx)) = alpha(STA_ID_per_RU(idx)) - Dec;
                            end
                            Succ(STA_ID_per_RU(idx), 1) = 0;
                        end
                        alpha(alpha < Min_alpha) = Min_alpha;

                    elseif Mode_UORA == 7 || Mode_UORA == 8 % 7: AIMD, 8:MIAD
                        Succ(STA_ID_per_RU, 2) = Succ(STA_ID_per_RU, 2) + 1;
                        Succ(STA_ID_per_RU, 1) = 0;
                        %alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Dec*log2(STA.CW(STA_ID_per_RU)'/CWmin);
                        if Mode_UORA == 7
                            alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Dec;
                        elseif Mode_UORA == 8
                            alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Dec*2.^(Succ(STA_ID_per_RU, 2)-1);
                        else
                            error('unexpected Mode_UORA!');
                        end
                        alpha(alpha < Min_alpha) = Min_alpha;

                    end

                % When Successful Transmission 
                elseif length(STA_ID_per_RU) == 1
                    if STA.Type(STA_ID_per_RU) == 1 % TX is associated?
                    %if ismember(STA_ID_per_RU, STA_ID_Data_TX)  % is data transmission?
                        Sim_Rst(STA_ID_per_RU, 2) = Sim_Rst(STA_ID_per_RU, 2) + 1; % count up # success (only for data transmission)
                    else

                    end
                    Sim_Rst_RU(j, 2) = Sim_Rst_RU(j, 2) + 1;
                    Sim_Rst(STA_ID_per_RU, 4) = Sim_Rst(STA_ID_per_RU, 4) + n_tf(STA_ID_per_RU);
                    n_tf(STA_ID_per_RU) = 0;

                    STA_ID_SUC = [STA_ID_SUC; STA_ID_per_RU];

                    if Mode_OBO_BEB
                        if Mode_UORA == 9
                            if  RU.AID(j) == 0 && (N_unassoc_In > 0 || N_unassoc_Out > 0)
                                if Opt_Random_N_RA_RU_0
                                    STA.CW(STA_ID_per_RU) = Tab_OCW( knnsearch(Tab_OCW(:,1), length(STA_ID_Assoc)), N_Random_RA_RU_0+1)*ones(1, length(STA_ID_per_RU) );
                                else
                                    STA.CW(STA_ID_per_RU) = Tab_OCW( knnsearch(Tab_OCW(:,1), length(STA_ID_Assoc)), N_RA_RU_0+1)*ones(1, length(STA_ID_per_RU) );
                                end
                                if ismember(STA_ID_per_RU, 1)
                                    log_OCW_assoc = [log_OCW_assoc; [i STA.CW(STA_ID_per_RU(1))] ];
                                end
                            elseif RU.AID(j) == 2045 && (N_unassoc_In > 0 || N_unassoc_Out > 0)
                                STA.CW(STA_ID_per_RU) = Tab_OCW( knnsearch(Tab_OCW(:,1), length(STA_ID_Unassoc)), N_RA_RU_2045+1)*ones(1, length(STA_ID_per_RU) );
                                log_OCW_unassoc = [log_OCW_unassoc; [i STA.CW(STA_ID_per_RU(1))] ];


                            end

                            % unassociated STA과 associated STA의 CW 갱신의 구분
                        else
                            STA.CW(STA_ID_per_RU) = STA.CWmin(STA_ID_per_RU);
                        end
                    elseif Mode_OBO_XIXD
                        if strcmp(Opt_Dec,'A')
                           STA.CW(STA_ID_per_RU) = round(STA.CW(STA_ID_per_RU) - Opt_Dec_val); 
                        elseif strcmp(Opt_Dec, 'M')
                           STA.CW(STA_ID_per_RU) = round(STA.CW(STA_ID_per_RU) / Opt_Dec_val);  
                        else
                            error('unexpected XIXD mode!');
                        end
                        if Mode_UORA == 9
                            % nothing to do (ignore OCWmin/max values)
                        else
                            STA_ID_CWmin = find(STA.CW(STA_ID_per_RU) <= STA.CWmin(STA_ID_per_RU));
                            STA.CW(STA_ID_per_RU(STA_ID_CWmin)) = STA.CWmin(STA_ID_per_RU(STA_ID_CWmin));
                        end
                    end

                    if Mode_UORA == 3
                        alpha = alpha + Inc;
                        alpha = min([alpha, Max_alpha]);

                    elseif Mode_UORA == 4
                        if On_Selfish_STA
                            Idx_Selfish_STA = find(STA_ID_per_RU <= N_Selfish_STA);  % 
                            Idx_Non_Selfish_STA = STA_ID_per_RU; 
                            Idx_Non_Selfish_STA(Idx_Non_Selfish_STA <= N_Selfish_STA) = [];
                            
                            if ~isempty(Idx_Selfish_STA)
                                % alpha(Idx_Selfish_STA) = alpha(Idx_Selfish_STA)
                            end
                    
                            if ~isempty(Idx_Non_Selfish_STA)
                                alpha(Idx_Non_Selfish_STA) = alpha(Idx_Non_Selfish_STA) + Inc;
                                alpha(alpha > Max_alpha) = Max_alpha; % note) non-selfish STA만 alpha가 감소하므로.
                            end
                            
                        elseif On_UORA_STD_STA
                            Idx_UORA_STD_STA = find(STA_ID_per_RU <= N_UORA_STD_STA);  % 
                            Idx_Non_UORA_STD_STA = STA_ID_per_RU; 
                            Idx_Non_UORA_STD_STA(Idx_Non_UORA_STD_STA <= N_UORA_STD_STA) = [];

                            if ~isempty(Idx_UORA_STD_STA)   % followw the Basic BEB operation
                               % do nothing because BEB already has done above. 
                            end

                            if ~isempty(Idx_Non_UORA_STD_STA)
                                alpha(Idx_Non_UORA_STD_STA) = alpha(Idx_Non_UORA_STD_STA) + Inc;
                                alpha(alpha > Max_alpha) = Max_alpha; % note) non-UORA STA만 alpha가 감소하므로, 구분하지 않음. 
                            end                                           
                        else
                            
                            alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) + Inc;
                            alpha(alpha > Max_alpha) = Max_alpha;
                        end
                    elseif Mode_UORA == 5
                        Succ(STA_ID_per_RU, 1) = Succ(STA_ID_per_RU, 1) + 1;
                        Succ(STA_ID_per_RU, 2) = 0;
                        %alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Dec*log2(STA.CW(STA_ID_per_RU)'/CWmin);
                        alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) + Inc*2.^(Succ(STA_ID_per_RU, 1)-1);
                        %alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) - Inc*log2(STA.CW(STA_ID_per_RU)'/CWmin);
                        alpha(alpha > Max_alpha) = Max_alpha;
                    elseif Mode_UORA == 6
                        for idx = 1:length(STA_ID_per_RU)
                            Succ(STA_ID_per_RU(idx), 1) = Succ(STA_ID_per_RU(idx), 1) + 1;
                            if Succ(STA_ID_per_RU(idx), 2) > 0
                               alpha(STA_ID_per_RU(idx)) = Init_alpha+Inc;
                            else
                               alpha(STA_ID_per_RU(idx)) = alpha(STA_ID_per_RU(idx)) + Inc; 
                            end
                            Succ(STA_ID_per_RU(idx), 2) = 0;
                        end 
                        alpha(alpha > Max_alpha) = Max_alpha;
                    elseif Mode_UORA == 7 || Mode_UORA == 8 % 7: AIMD, 8:MIAD
                        Succ(STA_ID_per_RU, 1) = Succ(STA_ID_per_RU, 1) + 1;
                        Succ(STA_ID_per_RU, 2) = 0;
                        if Mode_UORA == 7
                            alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) + Inc*2.^(Succ(STA_ID_per_RU, 1)-1);
                        elseif Mode_UORA == 8
                            alpha(STA_ID_per_RU) = alpha(STA_ID_per_RU) + Inc;
                        else
                            error('unexpected UORA mode!');
                        end

                        alpha(alpha > Max_alpha) = Max_alpha;
                    end

                end
            end
        end

        if ~isempty(STA_ID_TX)
            STA.OBO(STA_ID_TX) = floor(STA.CW(STA_ID_TX).*rand(1,length(STA_ID_TX)));
            if N_RU == 1
                Data_Rate = (N_sd(STA.Assigned_RU(STA_ID_TX))'.*N_bp(MCS(STA_ID_TX)+1).*R(MCS(STA_ID_TX)+1)*N_ss)/(T_dft+T_gi);% Data_Rate per RU
            else
                Data_Rate = (N_sd(STA.Assigned_RU(STA_ID_TX)).*N_bp(MCS(STA_ID_TX)+1).*R(MCS(STA_ID_TX)+1)*N_ss)/(T_dft+T_gi);% Data_Rate per RU
            end
            STA_ID_assoc_req = STA_ID_TX(STA.Type(STA_ID_TX) == 0);
            STA_ID_Data_TX = STA_ID_TX(STA.Type(STA_ID_TX) == 1);

                % case 1) Both Data and Assoc. req. are exist
            if ~isempty(STA_ID_assoc_req) && ~isempty(STA_ID_Data_TX)
                i = i + max(ceil( (L_PHY + (8*L_MPDU(STA_ID_TX)./Data_Rate))/U_slot )) ...
                  + ceil((2*L_PHY + L_BACK + L_Trigger + 3*SIFS)/U_slot);
                % case 2) Only Assoc. req. is exist
            elseif ~isempty(STA_ID_assoc_req)
                i = i + L_Assoc_Req_in_Slot ...
                  + ceil((2*L_PHY + L_BACK + L_Trigger + 3*SIFS)/U_slot);
                % case 3) Only Data is exist
            elseif ~isempty(STA_ID_Data_TX)
                i = i + max(ceil( (L_PHY + (8*L_MPDU(STA_ID_TX)./Data_Rate))/U_slot )) ...
                  + ceil((2*L_PHY + L_BACK + L_Trigger + 3*SIFS)/U_slot);
            else
                error('unexpected situation!');
            end 

            for k = 1:length(STA_ID_assoc_req)
                if ismember(STA_ID_assoc_req(k), STA_ID_SUC)
                     STA.Type(STA_ID_assoc_req(k)) = 1;  % association success
                     STA.Assoc_Delay(STA_ID_assoc_req(k), 2) = i - ceil((2*L_PHY + L_BACK + L_Trigger + 3*SIFS)/U_slot);
                     STA.Assoc_Delay(STA_ID_assoc_req(k), 3) = STA.Assoc_Delay(STA_ID_assoc_req(k), 2) ...
                                                               - STA.Assoc_Delay(STA_ID_assoc_req(k), 1);
                end
            end

    %         i = i + max(ceil( (L_PHY + (8*L_MPDU(STA_ID_TX)./Data_Rate))/U_slot )) ...
    %               + ceil((2*L_PHY + L_BACK + L_Trigger + 3*SIFS)/U_slot);

            STA.Assigned_RU = zeros(N_STA, 1); % re-initialization for RU-assignment
            Observation_N_assoc(N_TF, 1) = i;
            Observation_N_assoc(N_TF, 2) = length(find(STA.Type == 0));
            Observation_N_assoc(N_TF, 3) = length(find(STA.Type == 1));
            Observation_N_assoc = [Observation_N_assoc; zeros(1,3)];
        else
            i = i + 1;
            Observation_N_assoc(N_TF, 1) = i;
            Observation_N_assoc(N_TF, 2) = length(find(STA.Type == 0));
            Observation_N_assoc(N_TF, 3) = length(find(STA.Type == 1));
            Observation_N_assoc = [Observation_N_assoc; zeros(1,3)];
        end

        if Opt_Observation_TH && Next_Observation_Time_TH <= i     
            log_slot_time_SUC = [log_slot_time_SUC (i-sum(log_slot_time_SUC))];
            log_SUC = [log_SUC  (sum(Sim_Rst(:,2)) - sum(log_SUC))];
            Next_Observation_Time_TH = Next_Observation_Time_TH + Observation_Period;
        end

        if Opt_Random_N_RA_RU_0
           N_Random_RA_RU_0 = ceil(rand*N_RA_RU_0); 
           
           if Mode_UORA == 9
               STA.CW(1:N_STAs_existing) = Tab_OCW(knnsearch(Tab_OCW(:,1), N_STAs_existing), N_Random_RA_RU_0+1)*ones(1, N_STAs_existing);
           end
           
        end
    end

        Idle_RU = mean(Sim_Rst_RU(:,1)'./sum(Sim_Rst_RU'))
        success_RU = mean(Sim_Rst_RU(:,2)'./sum(Sim_Rst_RU'))  
        collision_RU = mean(Sim_Rst_RU(:,3)'./sum(Sim_Rst_RU'))

        access_STA = mean(Sim_Rst(:,1)'./N_TF)
        success_STA = mean(Sim_Rst(:,2)./Sim_Rst(:,1))  
        collision_STA = mean(Sim_Rst(:,3)./Sim_Rst(:,1))  



    %     figure;
    %     hold on;
    %     bar(1:3, mean(Sim_Rst_RU));
    %     title('Average performance (RUs)');
    %     X_label = {'IDLE', 'Success', 'Collision'};
    %     grid on;
    %     set(gca, 'XTick',1:3,'XTickLabel',X_label);
    %     hold off;
    % 
    % 
    %     figure;
    %     hold on;
    %     bar(1:3, mean(Sim_Rst));
    %     title('Average performance (STAs)');
    %     X_label = {'Access', 'Success', 'Collision'};
    %     grid on;
    %     set(gca, 'XTick',1:3,'XTickLabel',X_label);
    %     hold off;

        Th_per_STA = Sim_Rst(:,2).*(8*L_MPDU')./(N_slot*U_slot) *10^-6; % (Mbps);
        
        Avg_Th = mean(Th_per_STA); % (Mbps);
        Total_Th = sum(Th_per_STA); % (Mbps);
        Fairness = (sum(Th_per_STA)^2)/(N_STA*sum(Th_per_STA.^2)); % Jain's Fairness index
    %     Avg_Th1 = mean(Th_per_STA(1:N_STA_Control1))
    %     Avg_Th2 = mean(Th_per_STA((N_STA_Control1+1):(N_STA_Control1+N_STA_Control2)))
        Avg_Delay1 = mean(Sim_Rst(:,4)./Sim_Rst(:,2));

        if N_RU == 1
            Data_Rate = (N_sd(1)'*N_bp(MCS(1)+1).*R(MCS(1)+1)*N_ss)/(T_dft+T_gi);% Data_Rate per RU
        else
            Data_Rate = (N_sd(1)*N_bp(MCS(1)+1).*R(MCS(1)+1)*N_ss)/(T_dft+T_gi);% Data_Rate per RU
        end
        Avg_Delay2 = mean(( Sim_Rst(:,4)*( ceil((L_PHY + L_Trigger + SIFS)/U_slot) ...
                                         + ceil((L_PHY + max(8*L_MPDU./Data_Rate))/U_slot ) ...
                                         + ceil((L_PHY + L_BACK + 2*SIFS)/U_slot)) ...
                                         - ceil((L_PHY + L_BACK + 2*SIFS)/U_slot))./Sim_Rst(:,2))*U_slot*10^3;   % (msec)
                                     % this delay calculation way does not consider different L_MPDU, Data_Rate among STAs.
    %     figure;
    %     hold on;
    %     bar(1:2, [Total_Th Avg_Th] );
    %     title('Throughput (Mbps)');
    %     X_label = {'Total', 'Average'};
    %     grid on;
    %     set(gca, 'XTick',1:2,'XTickLabel',X_label);
    %     hold off;
    if Mode_UORA == 4 
        if On_Selfish_STA
            
            Th_per_Selfish_STA = Th_per_STA(1:N_Selfish_STA);
            Th_per_Non_Selfish_STA = Th_per_STA(N_Selfish_STA+1:N_STA);
            Total_Th_Selfish = sum(Th_per_Selfish_STA);
            Total_Th_Non_Selfish = sum(Th_per_Non_Selfish_STA);
            AVG_Th_Selfish = mean(Th_per_Selfish_STA);
            AVG_Th_Non_Selfish = mean(Th_per_Non_Selfish_STA);
            Result_Set2(idx_Sim, :) = [Total_Th_Selfish Total_Th_Non_Selfish AVG_Th_Selfish AVG_Th_Non_Selfish];
            RESULTS_TH_Selfish(:,n_sim) = Result_Set2(:,1);
            RESULTS_TH_Non_Selfish(:,n_sim) = Result_Set2(:,2);
            RESULTS_AVG_TH_Selfish(:,n_sim) = Result_Set2(:,3);
            RESULTS_AVG_TH_Non_Selfish(:,n_sim) = Result_Set2(:,4);            
        elseif On_UORA_STD_STA
            Th_per_UORA_STD_STA = Th_per_STA(1:N_UORA_STD_STA);
            Th_per_Non_UORA_STD_STA = Th_per_STA(N_UORA_STD_STA+1:N_STA);
            Total_Th_UORA_STD = sum(Th_per_UORA_STD_STA);
            Total_Th_Non_UORA_STD = sum(Th_per_Non_UORA_STD_STA);
            AVG_Th_UORA_STD = mean(Th_per_UORA_STD_STA);
            AVG_Th_Non_UORA_STD = mean(Th_per_Non_UORA_STD_STA);
            Result_Set2(idx_Sim, :) = [Total_Th_UORA_STD Total_Th_Non_UORA_STD AVG_Th_UORA_STD AVG_Th_Non_UORA_STD];
            RESULTS_TH_UORA_STD(:,n_sim) = Result_Set2(:,1);
            RESULTS_TH_Non_UORA_STD(:,n_sim) = Result_Set2(:,2);
            RESULTS_AVG_TH_UORA_STD(:,n_sim) = Result_Set2(:,3);
            RESULTS_AVG_TH_Non_UORA_STD(:,n_sim) = Result_Set2(:,4);  
        end
    end
        
        Result_Set(idx_Sim, :) = [N_STA Idle_RU success_RU collision_RU access_STA success_STA collision_STA Avg_Th Total_Th Fairness Avg_Delay1 Avg_Delay2];
    %     Result_Set = [N_STA Idle_RU success_RU collision_RU access_STA success_STA collision_STA Avg_Th Total_Th Fairness Avg_Delay1 Avg_Delay2];
    if Opt_Save_Per_STA_TH
      Result_Per_TH(n_sim*N_STAs_Assoc(idx_Sim)-(N_STAs_Assoc(idx_Sim)-1): n_sim*N_STAs_Assoc(idx_Sim),idx_Sim) = Th_per_STA;
    end
    
    end
    RESULTS_TH(:,n_sim)   = Result_Set(:,9);
    RESULTS_ACC(:,n_sim)  = Result_Set(:,5);
    RESULTS_SUC(:,n_sim)  = Result_Set(:,6); 
    RESULTS_COL(:,n_sim)  = Result_Set(:,7);
    RESULTS_FAIR(:,n_sim) = Result_Set(:,10);
    
    
    if Opt_Save_Per_STA_TH
%        PATH = strcat(pwd,'\',
%        Result_Per_TH(1: N_STAs_Assoc(idx_Sim),idx_Sim) = Th_per_STA; 
    end
    
end

RESULTS_SUMMARY = [mean(RESULTS_ACC')' mean(RESULTS_SUC')' mean(RESULTS_COL')' mean(RESULTS_TH')' mean(RESULTS_FAIR')']; 
RESULTS_ALL = [RESULTS_ACC RESULTS_SUC RESULTS_COL RESULTS_TH RESULTS_FAIR]; 

if Mode_UORA == 4 
        if On_Selfish_STA
            
            RESULTS_SUMMARY = [ RESULTS_SUMMARY mean(RESULTS_TH_Selfish')' mean(RESULTS_TH_Non_Selfish')' mean(RESULTS_AVG_TH_Selfish')' mean(RESULTS_AVG_TH_Non_Selfish')'];

        elseif On_UORA_STD_STA
            
             RESULTS_SUMMARY = [ RESULTS_SUMMARY mean(RESULTS_TH_UORA_STD')' mean(RESULTS_TH_Non_UORA_STD')' mean(RESULTS_AVG_TH_UORA_STD')' mean(RESULTS_AVG_TH_Non_UORA_STD')'];

        end
end
    
%     if Mode_UORA == 4
%         if Opt_Observation_Alpha && ismember(N_STA, Target_N_STAs)
%             mkdir([pwd, '\',Sim_ID]);
%             save_alpha_log = [log_slot_time; log_alpha];
%             save([pwd, '\',Sim_ID,'\','Log_Alpha_', num2str(N_STA),'STAs','.mat'], 'save_alpha_log');
%         end
%     end
% %     close all;
% 
% 
if Opt_Observation_TH
    figure;
    plot(cumsum(log_slot_time_SUC)*(U_slot*10^3), (8*L_MPDU(1))*log_SUC./(log_slot_time_SUC*U_slot)*10^-6);
    grid on;
    xlabel('msec');ylabel('throughput(Mb/s)');
    hold on;
%     xlim([0 1000]);
    hold off;
    R1(:,1) = (cumsum(log_slot_time_SUC)*(U_slot*10^3))';
    R1(:,2) = ((8*L_MPDU(1))*log_SUC./(log_slot_time_SUC*U_slot) *10^-6)';
end
% 
% figure;
% hold on;
% plot(Observation_N_assoc(1:length(Observation_N_assoc)-1,1)*U_slot, Observation_N_assoc(1:length(Observation_N_assoc)-1,2), 'r-');       
% plot(Observation_N_assoc(1:length(Observation_N_assoc)-1,1)*U_slot, Observation_N_assoc(1:length(Observation_N_assoc)-1,3), 'b-');  
% grid on;
% xlabel('time(sec)');ylabel('number of STAs');
% legend('unassoc.', 'assoc.', 'Location' , 'Best');
% hold off;
% 
% Assoc_Delay = STA.Assoc_Delay;
% Assoc_Delay(Assoc_Delay(:,3)==0, :) = [];
% Assoc_Delay = sortrows(Assoc_Delay,2); % to sec
% Assoc_Delay = Assoc_Delay*U_slot*10^3; % to msec
% 
% figure;
% plot(Assoc_Delay(:, 2)*10^-3, Assoc_Delay(:, 3), 'kx-'); 
% xlabel('time(sec)');ylabel('Association delay(msec)');
% grid on;
% % ylim([0 1000]);
% 
% figure;
% hold on;
% plot(log_OCW_assoc(:,1), log_OCW_assoc(:,2));
% % plot(log_OCW_unassoc(:,1), log_OCW_unassoc(:,2));
% legend('assoc');
% grid on;
% hold off;
% 
% figure;
% hold on;
% % plot(log_OCW_assoc(:,1), log_OCW_assoc(:,2));
% plot(log_OCW_unassoc(:,1), log_OCW_unassoc(:,2));
% legend( 'unassoc');
% grid on;
% hold off;

toc;
beep;pause(1);beep;pause(1);beep;pause(1);beep;pause(1);beep;