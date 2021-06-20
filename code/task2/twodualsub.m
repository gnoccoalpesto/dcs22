%% dual subgradient algorithm for task assignation
% subcase of costraint coupled linear program

% this is a variation of the algorithm used on task1:
% uses inequality constraint (instead of equality ones)

function...
    [primCost,dualCost,primRA,dualRA,consErr,...
                lam,ZZ,ZRA]=...
                           twodualsub(MAXITERS,NN,AA,AANW,bb,cc,gg,GG,HH,LBB,UBB)

    % linprog options                  
    lpoptions = optimoptions('linprog','Display','none');
    
    msglnt2=fprintf('iterations: %d  \n',MAXITERS);
    msglnt2=fprintf('initializing...')+msglnt2;
    
    % cost init
    primCost=zeros(MAXITERS,1);
    dualCost=zeros(MAXITERS,1);
    primRA=zeros(MAXITERS,1);
    dualRA=zeros(MAXITERS,1);
    % consensus and agents init
    consErr=zeros(1,NN,MAXITERS);
    ZZ=zeros(NN,NN,MAXITERS);
    ZRA=zeros(NN,NN,MAXITERS);
    % multipliers init
    lam=zeros(NN,NN,MAXITERS);
    vv=zeros(NN,NN); 
    % b splitting
    bi=bb/NN;
    
    
    msglnt2=fprintf(" done!\n")+msglnt2;
    msglnt1=fprintf('progress: ');
    msglnt3=fprintf("0%%");   
    
    
    % algorithm loop
    for tt = 1:MAXITERS-1
        % algorithm progression
        if mod(tt,MAXITERS/100)==0
            fprintf((repmat('\b', 1,msglnt3)));
            msglnt3=fprintf('%d%%',tt*100/MAXITERS);
        end       
        % step size
        alpha_t = 1.6*(1/tt)^0.65;
        
        % update averaging
        for ii=1:NN
            N_ii = find(AANW(:,ii) == 1)';
            for kk=1:NN 
                
                vv(ii,kk) = AA(ii,ii)*lam(ii,kk,tt);    

                for jj = N_ii
                    vv(ii,kk) = vv(ii,kk) + AA(ii,jj)*lam(jj,kk,tt);
                end
            end
        end
        
        for ii=1:NN
            % data extraction
            ci=cc((ii-1)*NN+1 : ii*NN)';
            Hi=HH(:,NN*(ii-1)+1:NN*ii);% assignation submatrix
            UBi=UBB((ii-1)*NN+1:ii*NN);% bounds
            LBi=LBB((ii-1)*NN+1:ii*NN);
            Gi=GG(ii,NN*(ii-1)+1:NN*ii);% constraints
            gi=gg(ii);

            % local solution
            [ZZ(ii,:,tt),~,~]=linprog((ci+(vv(ii,:))*Hi)', [],[],...
                    Gi,gi,LBi,UBi,lpoptions);
            
            % running average
            if tt==1
                ZRA(ii,:,tt) = ZZ(ii,:,tt);
            else
                ZRA(ii,:,tt) = (1/tt)*((tt-1)*ZRA(ii,:,tt-1)...
                                                  +ZZ(ii,:,tt));
            end

        end
        % gradient and multiplier updating
        for ii=1:NN
            
            Hi=HH(:,NN*(ii-1)+1:NN*ii);
            
            grad_ii=Hi*ZZ(ii,:,tt)'-bi;
            lam(ii,:,tt+1)=vv(ii,:)'+alpha_t*grad_ii;
        end
    % multiplier averaging
	lamAvg = mean(lam(:,:,tt));
    
    % cost and consensus error update
        for ii=1:NN
            ci=cc((ii-1)*NN+1 : ii*NN)';
            
            ffii = ci * ZZ(ii,:,tt)';
            primCost(tt) = primCost(tt) + ffii;
            %
            ffii = ci * ZRA(ii,:,tt)';
            primRA(tt) = primRA(tt) + ffii;
            %
            qqi = ci * ZZ(ii,:,tt)'...
                +lam(ii,:,tt) * ( ZZ(ii,:,tt) - bi' )';
            dualCost(tt) = dualCost(tt) + qqi;
            %
            qqi = ci *ZRA(ii,:,tt)'...
                + lam(ii,:,tt) * ( ZRA(ii,:,tt) - bi' )';
            dualRA(tt) = dualRA(tt) + qqi;
                
            consErr(:,:,tt) = consErr(:,:,tt) + abs(lam(ii,:,tt) - lamAvg);

        end
        
    end
    
    
   % last iteration (as above)
    tt=MAXITERS;

    % algorithm progress printing
    fprintf((repmat('\b', 1,msglnt1+msglnt3)));%display cleaning
    msglnt1=fprintf('progress: 100%%\n');
%     fprintf('progress: %d%%\n',100);
    
    for ii=1:NN
        for kk=1:NN
            N_ii = find(AANW(:,ii) == 1)';

            vv(ii,kk) = AA(ii,ii)*lam(ii,kk,tt);

            for jj = N_ii
            	vv(ii,kk) = vv(ii,kk) + AA(ii,jj)*lam(jj,kk,tt);
            end
        end
    end
    
    lamAvg = mean(lam(:,:,tt)); 
    
    for ii=1:NN
        
        ci=cc((ii-1)*NN+1 : ii*NN)';
        Hi=HH(:,NN*(ii-1)+1:NN*ii);
        UBi=UBB((ii-1)*NN+1:ii*NN);
        LBi=LBB((ii-1)*NN+1:ii*NN);
        Gi=GG(ii,NN*(ii-1)+1:NN*ii);
        gi=gg(ii);

        [ZZ(ii,:,tt),~,~]=linprog((ci+vv(ii,:)*Hi)', [],[],...
                    Gi,gi,LBi,UBi,lpoptions);

        if tt==1
            ZRA(ii,:,tt) = ZZ(ii,:,tt);
        else
            ZRA(ii,:,tt) = (1/tt)*((tt-1)*ZRA(ii,:,tt-1)...
                                             +ZZ(ii,:,tt));
        end
        
        ffii = ci * ZZ(ii,:,tt)';
        primCost(tt) = primCost(tt) + ffii;
        %
        ffii = ci * ZRA(ii,:,tt)';
        primRA(tt) = primRA(tt) + ffii;
        %
        qqi = ci * ZZ(ii,:,tt)'...
            +lam(ii,:,tt) * ( ZZ(ii,:,tt) - bi' )';
        dualCost(tt) = dualCost(tt) + qqi;
        %
        qqi = ci *ZRA(ii,:,tt)'...
            + lam(ii,:,tt) * ( ZRA(ii,:,tt) - bi' )';
        dualRA(tt) = dualRA(tt) + qqi;
         
        consErr(:,:,tt) = consErr(:,:,tt) + abs(lam(ii,:,tt) - lamAvg);
    end
    fprintf((repmat('\b', 1,msglnt1+msglnt2)));
    fprintf("\n");
end
