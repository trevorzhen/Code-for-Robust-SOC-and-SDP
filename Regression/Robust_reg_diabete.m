clc;clear;

dataDump = load('diabete.mat'); % http://www4.stat.ncsu.edu/~boos/var.select/diabetes.tab.txt
data = dataDump.data;

training=0.7;
data_model=data(end-floor(size(data,1)*training):end,:);
data_test=data(1:end-floor(size(data,1)*training)-1,:);

n=size(data,2)-1;
m=size(data_model,1);
rng(n);
A=data_model; 
b=A(:,end);
A(:,end)=[];

opt_gap =[];

for i =  1:9;
    Gam = i %G(i)
    
D=[eye(n) -eye(n); -eye(n) -eye(n); zeros(n,n) eye(n); zeros(1,n) ones(1,n)];
d=[zeros(2*n,1);ones(n,1);Gam];

Z=0.014*A;    Z(:,2)=0;

    % Upper Bound LDR approximation 
    model=xprog('Robust Regression');        % create a model, named 'Flow with uncertain demand'
    model.Param.Solver='MOSEK';
    x0=model.decision(1);                 % here-and-now decision
    x=model.decision(n);                   % here-and-now decision
    tau=model.decision(1);                 % here-and-now decision, epi
    y0=model.decision(size(d,1));          % LDR
    y1=model.decision(size(d,1),m);        % LDR
    z=model.random(m);                     % define random varaibles
    model.uncertain( norm(z) <= 1);        % uncertainty set                    
    model.min(  tau );                          % objective function with epi of second stage costs
    model.add(  sum(z)*x0 + z'*A*x + d'*(y0+y1*z)  - z'*b <= tau )  % constraints
    model.add( D(:,1:n)'*(y0+y1*z) == Z'*z.*x )         % constraints
    model.add( D(:,n+1:end)'*(y0+y1*z) == 0 ) 
    model.add( y0+y1*z >= 0 );                  % constraints
    model.solve;
    UB    =model.get% solve the problem
    UBtime=model.Solution.Time
    xldr  =x.get

    % Lower bound iterative procedure
    zeta_en = perms(1:n-1);
    zeta_en = zeta_en <= Gam;
    zeta_en = double(unique(zeta_en,'rows')) ;
    zeta_en = [zeta_en(:,1) zeros(size(zeta_en,1),1) zeta_en(:,2:end)] ;
    zeta_co = [];
    auximatrix = ones(2^Gam , Gam); count=1;
    
    for i = 1 : Gam ;
        minusone = nchoosek(1:Gam,i);
        for j = 1:size(minusone,1);
        auximatrix(count+j,minusone(j,:)) = -1; 
        end
        count= count+size(minusone,1);
    end

    for i = 1:size(zeta_en,1);
        ind1 = find(zeta_en(i,:)> 0);
        helpmatrix = zeros(2^Gam, n);
        helpmatrix(:,ind1) = auximatrix;
        zeta_co = [zeta_co;helpmatrix];
    end
    zeta_co = [zeta_co; zeros(1,10)];
    
    if Gam < 4
    
    model=xprog('Robust Regression');     % create a model, named 'Flow with uncertain demand'
    model.Param.Solver='MOSEK';
    x0=model.decision(1);                 % here-and-now decision
    x=model.decision(n);                   % here-and-now decision
    tau=model.decision(1);                 % here-and-now decision, epi                 
    model.min(  tau );                          % objective function with epi of second stage costs
    for i=1:size(zeta_co,1)
    model.add(  norm( ones(m,1)*x0 + A*x + Z*(zeta_co(i,:)'.*x)  - b ) <= tau )  % constraints
    end
    model.solve;
    Obj_opt= model.get
    x0_opt = x0.get
    x_opt  = x.get
    zeta_co_initial = zeta_co;
    
    else 
   
    violation_index = zeros(1, size(zeta_co,1)) ;
    for i = 1: size(zeta_co,1) ;
        if norm( ones(m,1)*x0_opt + A*x_opt + Z*(zeta_co(i,:)'.*x_opt)  - b ) > Obj_opt ;
           violation_index(i) = 1 ;
           violation_scenario(i,:) = zeta_co(i,:) ;
        end
    end
    sum(violation_index)
    zeta_vio=  zeta_co(find(violation_index==1),:);
    zeta_co = [zeta_co_initial; zeta_vio]; 
  
    model=xprog('Robust Regression');      
    model.Param.Solver='MOSEK';
    x0=model.decision(1);                  
    x=model.decision(n);                    
    tau=model.decision(1);                 
    model.min( tau );                      
    for i=1:size(zeta_co,1)
    model.add(  norm( ones(m,1)*x0 + A*x + Z*(zeta_co(i,:)'.*x)  - b ) <= tau )   
    end
    model.solve; 
    Obj_opt= model.get
    x0_opt = x0.get
    x_opt  = x.get;
    end
    
    opt_gap = [ opt_gap (UB - Obj_opt)/Obj_opt ]
end

