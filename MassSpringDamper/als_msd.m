classdef als_msd < handle
%% Constrained ALS algorithm -- Mass Spring Damper
%
% Function Inputs: 
% data_yk:              measurements 
% data_uk:              inputs
% data_datapts:         number of data points considered 
% data_start:           data to be ignored in the beginning until initial condition is negligible - optional, default is 100
% N:                    autocovariance lags
% model.A:
% model.B:              
% model.C: 
% model.G:              
% estimator.Q:          Qw initial guess
% estimator.R:          Rv initial guess
% options_dt            Disctrete model time step
 
% Function outputs : 
% Qest_cell:            cell containing estimated Qw 
% Rest_cell:            cell containing estimated Rv
% Phi:                  vector containing least-squares portion of objective

    
    properties
        Ain                                         % Abar --> A-ALC
        data_yk                                     % measurement data
        data_uk                                     % input data
        data_xhat_                                  % calculated state estimate data
        data_datapts                                % number of points considered
        data_start                                  % data to be ignored in the beginning until initial condition is negligible
        datatrun                                    % data_datapts-data_start
        Eyy
        Eyy1
        Eyy2
        Eyy3
        Eyy4
        Eyy5
        Eyy6
        Eyy7
        estimator_Q                                 % Initial process noise matrix 
        estimator_R                                 % Initial measurement noise matrix
        estimator_L                                 % Initial Kalman gain
        estimator_errVar                            % variable noise structure
        ga                                          % number of columns of G
        gamma
        inntrun
        M1
        M2
        model_Aa                                    % A matrix of model
        model_Ba                                    % B matrix of model
        model_Ca                                    % C matrix of model
        model_Ga                                    % G matrix of model
        model_xhat0                                 % initial state
        model_c_scale
        N                                           % Number of lags
        na                                          % Number of rows of G
        OO                                          % Observability matrix where A is Ain
        options_dt                                  % Discrete time model time step
        pa                                          % Number of columns of B
        Phi                                         % Leas Square objective function
        Qest_cell                                   % Process noise estimate
        Rest_cell                                   % Measurement noise estimate
        A_scr
        Wm
       
    end
    
    methods
        function obj = als_msd(data,model,estimator,options)                                       % Constructor method
            obj.gamma = 0;
            obj.data_yk            = data.yk;
            obj.data_uk            = data.uk;
            obj.data_datapts       = data.datapts;
            obj.data_start         = data.start;
            obj.N                  = data.N;            
            obj.model_Aa           = model.A;
            obj.model_Ba           = model.B;
            obj.model_Ca           = model.C;
            obj.model_Ga           = model.G;
            obj.model_xhat0        = model.xhat0;
            obj.model_c_scale      = model.c_scale;
            obj.estimator_Q        = estimator.Q;
            obj.estimator_R        = estimator.R; 
            obj.estimator_errVar   = estimator.errVar;
            obj.options_dt         = options.dt;

            [~,~,~,obj.estimator_L{1}] = kalman(ss(obj.model_Aa{1},[obj.model_Ba{1} obj.model_Ga],obj.model_Ca{1}, 0, options.dt),obj.estimator_Q,obj.estimator_R,0);
            [~,~,~,obj.estimator_L{2}] = kalman(ss(obj.model_Aa{2},[obj.model_Ba{2} obj.model_Ga],obj.model_Ca{2}, 0, options.dt),obj.estimator_Q,obj.estimator_R,0);
            [~,~,~,obj.estimator_L{3}] = kalman(ss(obj.model_Aa{3},[obj.model_Ba{3} obj.model_Ga],obj.model_Ca{3}, 0, options.dt),obj.estimator_Q,obj.estimator_R,0);
            [~,~,~,obj.estimator_L{4}] = kalman(ss(obj.model_Aa{4},[obj.model_Ba{4} obj.model_Ga],obj.model_Ca{4}, 0, options.dt),obj.estimator_Q,obj.estimator_R,0);                 
            [~,~,~,obj.estimator_L{5}] = kalman(ss(obj.model_Aa{5},[obj.model_Ba{5} obj.model_Ga],obj.model_Ca{5}, 0, options.dt),obj.estimator_Q,obj.estimator_R,0);
            [~,~,~,obj.estimator_L{6}] = kalman(ss(obj.model_Aa{6},[obj.model_Ba{6} obj.model_Ga],obj.model_Ca{6}, 0, options.dt),obj.estimator_Q,obj.estimator_R,0);
            [~,~,~,obj.estimator_L{7}] = kalman(ss(obj.model_Aa{7},[obj.model_Ba{7} obj.model_Ga],obj.model_Ca{7}, 0, options.dt),obj.estimator_Q,obj.estimator_R,0);
            
            [obj.na,obj.ga]        = size(obj.model_Ga);
            obj.pa                 = size(obj.model_Ca{1},1);
                        
            if isfield(data,'xhatk')
                obj.data_xhat_     = data.xhatk;
            else
                xhat               = zeros(obj.na,obj.data_datapts);
                xhat_              = zeros(obj.na,obj.data_datapts);
                xhat_(1:n,1)       = model.xhat0;  
                for i = 1:obj.data_datapts
                    xhat(:,i) = xhat_(:,i) + obj.estimator_L*(obj.data_yk(:,i)-obj.model_Ca*xhat_(:,i));
                    xhat_(:,i+1) = obj.model_Aa*xhat(:,i) + obj.model_Ba*obj.data_uk(:,i);
                end
                obj.data_xhat_ = xhat_(:,1:end-1);
            end
            
        end
                 
        function obj = sdp_mrQ_diag_con(obj)
            obj.datatrun = obj.data_datapts - obj.data_start;
            obj.inntrun = [obj.data_yk(:,obj.data_start+1:obj.datatrun/7)-obj.model_Ca{1}*obj.data_xhat_(:,obj.data_start+1:obj.datatrun/7) ...                 % Initial innovations
                            obj.data_yk(:,1*obj.datatrun/7+1:2*obj.datatrun/7)-obj.model_Ca{2}*obj.data_xhat_(:,1*obj.datatrun/7+1:2*obj.datatrun/7)...
                            obj.data_yk(:,2*obj.datatrun/7+1:3*obj.datatrun/7)-obj.model_Ca{3}*obj.data_xhat_(:,2*obj.datatrun/7+1:3*obj.datatrun/7)...
                            obj.data_yk(:,3*obj.datatrun/7+1:4*obj.datatrun/7)-obj.model_Ca{4}*obj.data_xhat_(:,3*obj.datatrun/7+1:4*obj.datatrun/7)...      
                            obj.data_yk(:,4*obj.datatrun/7+1:5*obj.datatrun/7)-obj.model_Ca{5}*obj.data_xhat_(:,4*obj.datatrun/7+1:5*obj.datatrun/7)...      
                            obj.data_yk(:,5*obj.datatrun/7+1:6*obj.datatrun/7)-obj.model_Ca{6}*obj.data_xhat_(:,5*obj.datatrun/7+1:6*obj.datatrun/7)...
                            obj.data_yk(:,6*obj.datatrun/7+1:7*obj.datatrun/7)-obj.model_Ca{7}*obj.data_xhat_(:,6*obj.datatrun/7+1:7*obj.datatrun/7)];  
            
            obj.oneColumn_autocorrelation
            obj.Wm = eye(length(obj.Eyy));                                  % Use identity-based weighting
            obj.LSconstant_matrix_con
            obj.als_diag_con
         end
        
        function obj = oneColumn_autocorrelation(obj)                                           % Calculation of Autocorrelations for one column ALS
            
             obj.datatrun = obj.data_datapts - obj.data_start;
             obj.Eyy = [];
             obj.Eyy1 = [];
             obj.Eyy2 = [];
             obj.Eyy3 = [];
             obj.Eyy4 = [];
             obj.Eyy5 = [];
             obj.Eyy6 = [];
             obj.Eyy7 = [];

             for i = 0:obj.N-1
                 temp1 = obj.inntrun(:,i+1:obj.datatrun/7)*obj.inntrun(:,1:obj.datatrun/7-i)';                         % Rajamani thesis equation 3.14
                 temp1 = temp1./(obj.datatrun/7-i);
                 obj.Eyy1 = [obj.Eyy1; temp1];
                 
                 temp2 = obj.inntrun(:,i+obj.datatrun/7:2*obj.datatrun/7)*obj.inntrun(:,obj.datatrun/7:2*obj.datatrun/7-i)';                         
                 temp2 = temp2./(obj.datatrun/7-i);
                 obj.Eyy2 = [obj.Eyy2; temp2];
                 
                 temp3 = obj.inntrun(:,i+2*obj.datatrun/7:3*obj.datatrun/7)*obj.inntrun(:,2*obj.datatrun/7:3*obj.datatrun/7-i)';   
                 temp3 = temp3./(obj.datatrun/7-i);
                 obj.Eyy3 = [obj.Eyy3; temp3];
                 
                 temp4 = obj.inntrun(:,i+3*obj.datatrun/7:4*obj.datatrun/7)*obj.inntrun(:,3*obj.datatrun/7:4*obj.datatrun/7-i)';                         
                 temp4 = temp4./(obj.datatrun/7-i);
                 obj.Eyy4 = [obj.Eyy4; temp4];
                 
                 temp5 = obj.inntrun(:,i+4*obj.datatrun/7:5*obj.datatrun/7)*obj.inntrun(:,4*obj.datatrun/7:5*obj.datatrun/7-i)';                        
                 temp5 = temp5./(obj.datatrun/7-i);
                 obj.Eyy5 = [obj.Eyy5; temp5];
                 
                 temp6 = obj.inntrun(:,i+5*obj.datatrun/7:6*obj.datatrun/7)*obj.inntrun(:,5*obj.datatrun/7:6*obj.datatrun/7-i)';                         
                 temp6 = temp6./(obj.datatrun/7-i);
                 obj.Eyy6 = [obj.Eyy6; temp6];
                 
                 temp7 = obj.inntrun(:,i+6*obj.datatrun/7:end)*obj.inntrun(:,6*obj.datatrun/7:end-i)';                         
                 temp7 = temp7./(obj.datatrun/7-i);
                 obj.Eyy7 = [obj.Eyy7; temp7];
             end
             
             obj.Eyy1 = obj.Eyy1(:);
             obj.Eyy2 = obj.Eyy2(:);
             obj.Eyy3 = obj.Eyy3(:);
             obj.Eyy4 = obj.Eyy4(:);
             obj.Eyy5 = obj.Eyy5(:);
             obj.Eyy6 = obj.Eyy6(:);
             obj.Eyy7 = obj.Eyy7(:);
             obj.Eyy = obj.Eyy1;
        end
     
        function obj = LSconstant_matrix_con(obj)                                               % Building the constant matrix for the LS problem
           
        obj.Ain{1} = obj.model_Aa{1}-obj.model_Aa{1}*obj.estimator_L{1}*obj.model_Ca{1};                                    % Abar, Rajamani thesis equation 3.4
        obj.Ain{2} = obj.model_Aa{2}-obj.model_Aa{2}*obj.estimator_L{2}*obj.model_Ca{2};                                   
        obj.Ain{3} = obj.model_Aa{3}-obj.model_Aa{3}*obj.estimator_L{3}*obj.model_Ca{3};                                   
        obj.Ain{4} = obj.model_Aa{4}-obj.model_Aa{4}*obj.estimator_L{4}*obj.model_Ca{4};  
        obj.Ain{5} = obj.model_Aa{5}-obj.model_Aa{5}*obj.estimator_L{5}*obj.model_Ca{5};
        obj.Ain{6} = obj.model_Aa{6}-obj.model_Aa{6}*obj.estimator_L{6}*obj.model_Ca{6};
        obj.Ain{7} = obj.model_Aa{7}-obj.model_Aa{7}*obj.estimator_L{7}*obj.model_Ca{7};    

        obj.OO{1} = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO{1} = [obj.OO{1} ; obj.model_Ca{1}*temp];                                                                 % 1: Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain{1};
        end

        obj.M1{1} = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain{1},obj.model_Ga*II*obj.model_Ga');                                                       % 1: Qw part of lyapunov equation, Rajamani eq 3.8
                obj.M1{1}(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2{1} = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain{1},obj.model_Aa{1}*obj.estimator_L{1}*II*obj.estimator_L{1}'*obj.model_Aa{1}');          % 1: Rv part of lyapunov equation, Rajamani eq 3.8
                obj.M2{1}(:,i) = t2(:);
                i = i+1;
            end
        end
        
        obj.OO{2} = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO{2} = [obj.OO{2} ; obj.model_Ca{2}*temp];                                                                 % 2: Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain{2};
        end

        obj.M1{2} = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain{2},obj.model_Ga*II*obj.model_Ga');                                                       % 2: Qw part of lyapunov equation, Rajamani eq 3.8
                obj.M1{2}(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2{2} = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain{2},obj.model_Aa{2}*obj.estimator_L{2}*II*obj.estimator_L{2}'*obj.model_Aa{2}');          % 2: Rv part of lyapunov equation, Rajamani eq 3.8
                obj.M2{2}(:,i) = t2(:);
                i = i+1;
            end
        end
        
        obj.OO{3} = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO{3} = [obj.OO{3} ; obj.model_Ca{3}*temp];                                                                 % 3: Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain{3};
        end

        obj.M1{3} = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain{3},obj.model_Ga*II*obj.model_Ga');                                                       % 3: Qw part of lyapunov equation, Rrajamani eq 3.8
                obj.M1{3}(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2{3} = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain{3},obj.model_Aa{3}*obj.estimator_L{3}*II*obj.estimator_L{3}'*obj.model_Aa{3}');          % 3: Rv part of lyapunov equation, Rajamani eq 3.8
                obj.M2{3}(:,i) = t2(:);
                i = i+1;
            end
        end
        
        obj.OO{4} = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO{4} = [obj.OO{4} ; obj.model_Ca{4}*temp];                                                                 % 4: Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain{4};
        end

        obj.M1{4} = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain{4},obj.model_Ga*II*obj.model_Ga');                                                       % 4: Qw part of lyapunov equation, Rajamani eq 3.8
                obj.M1{4}(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2{4} = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain{4},obj.model_Aa{4}*obj.estimator_L{4}*II*obj.estimator_L{4}'*obj.model_Aa{4}');          % 4: Rv part of lyapunov equation, Rajamani eq 3.8
                obj.M2{4}(:,i) = t2(:);
                i = i+1;
            end
        end
        
        obj.OO{5} = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO{5} = [obj.OO{5} ; obj.model_Ca{5}*temp];                                                                 % 5: Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain{5};
        end

        obj.M1{5} = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain{5},obj.model_Ga*II*obj.model_Ga');                                                       % 5: Qw part of lyapunov equation, Rajamani eq 3.8
                obj.M1{5}(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2{5} = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain{5},obj.model_Aa{5}*obj.estimator_L{5}*II*obj.estimator_L{5}'*obj.model_Aa{5}');          % 5: Rv part of lyapunov equation, Rajamani eq 3.8
                obj.M2{5}(:,i) = t2(:);
                i = i+1;
            end
        end
        
        obj.OO{6} = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO{6} = [obj.OO{6} ; obj.model_Ca{6}*temp];                                                                 % 6: Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain{6};
        end

        obj.M1{6} = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain{6},obj.model_Ga*II*obj.model_Ga');                                                       % 6: Qw part of lyapunov equation, Rajamani eq 3.8
                obj.M1{6}(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2{6} = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain{6},obj.model_Aa{6}*obj.estimator_L{6}*II*obj.estimator_L{6}'*obj.model_Aa{6}');          % 6: Rv part of lyapunov equation, Rajamani eq 3.8
                obj.M2{6}(:,i) = t2(:);
                i = i+1;
            end
        end
        
        obj.OO{7} = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO{7} = [obj.OO{7} ; obj.model_Ca{7}*temp];                                                                 % 7: Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain{7};
        end

        obj.M1{7} = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain{7},obj.model_Ga*II*obj.model_Ga');                                                       % 7: Qw part of lyapunov equation, RRRrajamani eq 3.8
                obj.M1{7}(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2{7} = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain{7},obj.model_Aa{7}*obj.estimator_L{7}*II*obj.estimator_L{7}'*obj.model_Aa{7}');          % 7: Rv part of lyapunov equation, Rajamani eq 3.8
                obj.M2{7}(:,i) = t2(:);
                i = i+1;
            end
        end
        end
        
        function obj = als_diag_con(obj)                                                        % Diagonal ALS
 
            % Building block diagonal matrix A
            PSI = eye(obj.pa);                                                                                              % 1
            for  i= 1:obj.N-1
                PSI = [PSI; -obj.model_Ca{1}*obj.Ain{1}^(i-1)*obj.model_Aa{1}*obj.estimator_L{1}];
            end
            OOtemp = kron(obj.model_Ca{1},obj.OO{1});
            PSItemp = kron(eye(obj.pa),PSI);
            LH1 = OOtemp*obj.M1{1};   
            LH2 = OOtemp*obj.M2{1} + PSItemp;           
            As_diag1 = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];
             
            clear PSI OOtemp PSItemp LH1 LH2
             
            PSI = eye(obj.pa);                                                                                              % 2
            for  i= 1:obj.N-1
                PSI = [PSI; -obj.model_Ca{2}*obj.Ain{2}^(i-1)*obj.model_Aa{2}*obj.estimator_L{2}];
            end
            OOtemp = kron(obj.model_Ca{2},obj.OO{2});
            PSItemp = kron(eye(obj.pa),PSI);
            LH1 = OOtemp*obj.M1{2};   
            LH2 = OOtemp*obj.M2{2} + PSItemp;           
            As_diag2 = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];
             
            clear PSI OOtemp PSItemp LH1 LH2
             
            PSI = eye(obj.pa);                                                                                              % 3
            for  i= 1:obj.N-1
                PSI = [PSI; -obj.model_Ca{3}*obj.Ain{3}^(i-1)*obj.model_Aa{3}*obj.estimator_L{3}];
            end
            OOtemp = kron(obj.model_Ca{3},obj.OO{3});
            PSItemp = kron(eye(obj.pa),PSI);
            LH1 = OOtemp*obj.M1{3};   
            LH2 = OOtemp*obj.M2{3} + PSItemp;           
            As_diag3 = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];
             
            clear PSI OOtemp PSItemp LH1 LH2
             
            PSI = eye(obj.pa);                                                                                              % 4
            for  i= 1:obj.N-1
                PSI = [PSI; -obj.model_Ca{4}*obj.Ain{4}^(i-1)*obj.model_Aa{4}*obj.estimator_L{4}];
            end
            OOtemp = kron(obj.model_Ca{4},obj.OO{4});
            PSItemp = kron(eye(obj.pa),PSI);
            LH1 = OOtemp*obj.M1{4};   
            LH2 = OOtemp*obj.M2{4} + PSItemp;           
            As_diag4 = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];
             
            clear PSI OOtemp PSItemp LH1 LH2
             
            PSI = eye(obj.pa);                                                                                             % 5
            for  i= 1:obj.N-1
                PSI = [PSI; -obj.model_Ca{5}*obj.Ain{5}^(i-1)*obj.model_Aa{5}*obj.estimator_L{5}];
            end
            OOtemp = kron(obj.model_Ca{5},obj.OO{5});
            PSItemp = kron(eye(obj.pa),PSI);
            LH1 = OOtemp*obj.M1{5};   
            LH2 = OOtemp*obj.M2{5} + PSItemp;           
            As_diag5 = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];
             
            clear PSI OOtemp PSItemp LH1 LH2
             
            PSI = eye(obj.pa);                                                                                             % 6
            for  i= 1:obj.N-1
                PSI = [PSI; -obj.model_Ca{6}*obj.Ain{6}^(i-1)*obj.model_Aa{6}*obj.estimator_L{6}];
            end
            OOtemp = kron(obj.model_Ca{6},obj.OO{6});
            PSItemp = kron(eye(obj.pa),PSI);
            LH1 = OOtemp*obj.M1{6};   
            LH2 = OOtemp*obj.M2{6} + PSItemp;           
            As_diag6 = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];
             
            clear PSI OOtemp PSItemp LH1 LH2
             
            PSI = eye(obj.pa);                                                                                             % 7
            for  i= 1:obj.N-1
                PSI = [PSI; -obj.model_Ca{7}*obj.Ain{7}^(i-1)*obj.model_Aa{7}*obj.estimator_L{7}];
            end
            OOtemp = kron(obj.model_Ca{7},obj.OO{7});
            PSItemp = kron(eye(obj.pa),PSI);
            LH1 = OOtemp*obj.M1{7};   
            LH2 = OOtemp*obj.M2{7} + PSItemp;           
            As_diag7 = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];

            clear PSI OOtemp PSItemp LH1 LH2
             
            H = [blkdiag(As_diag1'*obj.Wm*As_diag1, As_diag2'*obj.Wm*As_diag2, As_diag3'*obj.Wm*As_diag3, As_diag4'*obj.Wm*As_diag4, As_diag5'*obj.Wm*As_diag5, As_diag6'*obj.Wm*As_diag6, As_diag7'*obj.Wm*As_diag7), zeros(21,5);zeros(5,26)];
            f = [-As_diag1'*obj.Wm*obj.Eyy1; -As_diag2'*obj.Wm*obj.Eyy2; -As_diag3'*obj.Wm*obj.Eyy3; -As_diag4'*obj.Wm*obj.Eyy4; -As_diag5'*obj.Wm*obj.Eyy5; -As_diag6'*obj.Wm*obj.Eyy6; -As_diag7'*obj.Wm*obj.Eyy7];

            % Noise Covariance Constraints
            Aeq = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1  0   -1     0  0; 
                   0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0 -1    0    -1  0;
                   0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0  0    0     0 -1;
                   0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1  0   -2     0  0;
                   0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0 -1    0    -2  0;
                   0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0  0    0     0 -1;
                   0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1  0   -3     0  0;
                   0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0  0 -1    0    -3  0;
                   0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0  0  0    0     0 -1;
                   0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 -1  0   -4     0  0;
                   0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0  0 -1    0    -4  0;
                   0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0  0  0    0     0 -1;
                   0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 -1  0   -5     0  0;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0  0 -1    0    -5  0;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0  0  0    0     0 -1;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -1  0   -6     0  0;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0  0 -1    0    -6  0;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0  0  0    0      0 -1;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 -1  0   -7     0  0;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0  0 -1    0    -7  0;
                   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1  0  0    0     0 -1];
            Beq = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

        	Xest_diag = quadprog(H, [f;zeros(5,1)], -eye(7*(obj.pa+obj.ga)+5), zeros(7*(obj.ga+obj.pa)+5,1), Aeq, Beq); 

            if prod(Xest_diag)==0 
                fprintf('Warning: Covariance estimate(s) is (are) at constraints! You may have bad data! \n')
            end
              
            obj.Qest_cell{1} = diag(Xest_diag([1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23 24 25])); 
            obj.Rest_cell{1} = diag(Xest_diag([3 6 9 12 15 18 21 26])); 

            obj.A_scr = blkdiag(As_diag1, As_diag2, As_diag3, As_diag4, As_diag5, As_diag6, As_diag7);
            phi = blkdiag(As_diag1,As_diag2,As_diag3,As_diag4,As_diag5,As_diag6,As_diag7)*Xest_diag(1:21)-[obj.Eyy1; obj.Eyy2; obj.Eyy3; obj.Eyy4; obj.Eyy5; obj.Eyy6; obj.Eyy7];
            obj.Phi = phi;
        end  
    end
end

