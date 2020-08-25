% Sensor data generator

clc;clear;

% Input Data

Results = [];
sizz = [10 20 30 40 50 60 70];

for siz =  1: size(sizz,2);
    
    for iter = 1 : 10;

      M = sizz(siz);                 % number of fixed points
      N = round(0.4*M);              % number of free points
      rng(iter)

      r = -1 + 2*rand(M,2);
      r = r./(max(abs(r),[],2)*[ 1 1]);
  
      R = 0.3;                        % a_head: abs dev. from fixed points ++

    % first N columns of A correspond to free points, last M columns correspond to fixed points
    C = zeros(N^2,M+N);

    for i = 1 : N 
        for j = 1 : N 
            C(j+(i-1)*N,i) = 1 ;
            C(j+(i-1)*N,randi([N+1,N+M])) = -1;
        end
    end

    Chelp = zeros(N, M+N);
    for i = 1: N-1;
         Chelp(i,i) = 1;
         Chelp(i,i+1) = -1;   
    end
         Chelp(N,1) = 1;
         Chelp(N,N) = -1;   
     C = [C ; Chelp];

    indexhelp = find(sum(abs(C))==0);
    if size(indexhelp,2) ~= 0 ;
        for i = 1 : size(indexhelp,2)
            Chelp = zeros(N,M+N);
            Chelp(:,1:N)= eye(N);
            Chelp(:,indexhelp(i)) = -1;
            C = [C;Chelp];
        end
    end

    C = unique(C,'rows');

    Cx = C(:,1:N);
    Cy = C(:,N+1:end);
    nolinks = size(C,1); 

    T = sqrt(2*M);
    D = [eye(2*M) -eye(2*M);-eye(2*M) -eye(2*M); zeros(2*M) eye(2*M); zeros(1,2*M) ones(1,2*M)];
    d = [zeros(4*M,1) ; ones(2*M,1);  T]; 

    % Convert to FME-type: Ax <= b

    A_fme = [d' zeros(1,M*2) ;  -D' [eye(M*2); zeros(size(D,2)-M*2,M*2)] ; D' -[eye(M*2);zeros(size(D,2)-M*2,M*2)] ;  -eye(size(D,1)) zeros(size(D,1), M*2)];
    b_fme = [1; zeros(size(D,1)+size(D,2)*2,1)];
    nadj = size(D,1);

    fprintf(1,'Computing the optimal locations of the free points...');

    model=xprog('Sensor LDR');        
    model.Param.Solver='MOSEK';
    x    =    model.decision(size(Cx,2),2);                       
    epi  =    model.decision(1,1);
    lambda0 = model.decision(nadj,1);
    lambda1 = model.decision(nadj,nolinks*2);

    V   = model.random(nolinks,2);

    for i = 1:nolinks 
    model.uncertain( norm(V(i,:)) <= 1);
    end

    model.min( epi  );   
    model.add( A_fme(:,1:nadj)*(lambda0+lambda1*[V(:,1);V(:,2)]) + A_fme(:,nadj+1:end)*reshape(R*(Cy'*V)', M*2 ,1) + b_fme*(sum(sum(V.*(Cx*x))) + sum(sum((Cy'*V).*r)) - epi) <= 0 )
    model.solve;
    epi = model.get
    UB = epi ;
    time = model.Solution.Time;
    x_help = x.get;
    xd = x;

    fprintf(1,'Done! \n');

    CBSV=[];nCBS=0;
    
for j= 1: size(A_fme,1);
       model=xprog('CBS');       
       model.Param.Solver='MOSEK';
       model.Param.IsPrint=0;
        V    =    model.decision(nolinks,2);       
       model.max( A_fme(j,1:nadj)*(lambda0.get+lambda1.get*[V(:,1);V(:,2)]) + A_fme(j,nadj+1:end)*reshape(R*(Cy'*V)', M*2 ,1) + b_fme(j)*(sum(sum((Cx*x_help).*V)) + sum(sum((Cy'*V).*r)) - epi)   );    % objective function          
       for i = 1:nolinks 
       model.add( norm(V(i,:)) <= 1);
       end
       model.solve;   
        if j/10 - round(j/10) ==0;
            j
        end
        if norm(model.get) >= -0.001 ;
            nCBS = nCBS +1;
        CBSV(:,:,nCBS) = V.get;
        end
end
    
        CBSz = []; 
        for j=1:nCBS;
            model=xprog('SensorPrimal LB');        
            model.Param.Solver='MOSEK';
            model.Param.IsPrint=0;
            z = model.decision(2*M,1);
            za = model.decision(2*M,1);
            model.max( sum(sum(R*(Cy'*CBSV(:,:,j)).*reshape(z,2,M)'))   );   
            model.add( D*[z;za] <= d )
            model.solve;                                    
            CBSz = [CBSz z.get];
        end
        nCBS = size(CBSz,2);

        model=xprog('Sensor Primal LB');       
        model.Param.Solver='MOSEK';
        x    =    model.decision(size(Cx,2),2);                       
        epi  =    model.decision(1,1);
        y    =    model.decision(size(Cy,2)*2,nCBS);
        z    =    model.decision(size(C,1),nCBS);
        model.min(  epi );  
        for j = 1:nCBS; 
            model.add( reshape(y(:,j),size(Cy,2),2) == r + R*reshape(CBSz(:,j),2,M)')
            for i=1:size(C,1);
                model.add( norm(C(i,:)*[x;reshape(y(:,j),size(Cy,2),2)])  <= z(i,j) )
            end
            model.add( sum(z(:,j)) <= epi ) 
        end
        model.solve;
        LBp = model.get
        Results(:,iter,siz) = [UB; LBp ; time]
    end
end

%% Plots

x = [xd;r]    
linewidth = 1;          % in points;  width of dotted lines
free_sum = x(1:N,:);
figure;
dots = plot(free_sum(:,1), free_sum(:,2), 'or', r(:,1), r(:,2), 'bs');
set(dots(1),'MarkerFaceColor','red');
hold on
legend('Free points','Fixed points','Location','Best');
for i=1:nolinks
  ind = find(C(i,:));
  line2 = plot(x(ind,1), x(ind,2), ':k');
  hold on
  set(line2,'LineWidth',linewidth);
end
axis([-1.2 1.2 -1.2 1.2]) ;
axis equal;

