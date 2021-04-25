function...
    [primCost,dualCost,primRA,dualRA,consErr,...
                lam,ZZ,ZRA,MAXITERS]=...
                           twodualsub(NN,AA,AANW,bb,cc,gg,GG,HH,LBB,UBB)
                       
    lpoptions = optimoptions('linprog','Display','none');

    msglnt=0;
      
    MAXITERS =1e3;
    fprintf('iterations: %d  \n',MAXITERS);

    fprintf('initializing...\n');
    primCost=zeros(MAXITERS,1);
    dualCost=zeros(MAXITERS,1);
    primRA=zeros(MAXITERS,1);
    dualRA=zeros(MAXITERS,1);
    
    consErr=zeros(1,NN,MAXITERS);
    ZZ=zeros(NN,NN,MAXITERS);
    ZRA=zeros(NN,NN,MAXITERS);

    %lam=0.2*ones(NN,NN,MAXITERS);
    %lam=rand(NN,NN,MAXITERS)
    lam=zeros(NN,NN,MAXITERS);
    vv=zeros(NN,NN); 

    bi=bb/NN;
    
    for tt = 1:MAXITERS-1
        if mod(tt,MAXITERS/100)==0
            msglnt=fprintf('progress: %d%%\n',tt*100/MAXITERS);
        end
        
        alpha_t = 1.6*(1/tt)^0.65;

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
            ci=cc((ii-1)*NN+1 : ii*NN)';
            Hi=HH(:,NN*(ii-1)+1:NN*ii);
            UBi=UBB((ii-1)*NN+1:ii*NN);
            LBi=LBB((ii-1)*NN+1:ii*NN);

%        Gi=GG(:,NN*(ii-1)+1:NN*ii);
%        gi=gg((ii-1)*NN+1 : ii*NN);
        Gi=GG(ii,NN*(ii-1)+1:NN*ii);
        gi=gg(ii);

         
            [ZZ(ii,:,tt),~,~]=linprog((ci+(vv(ii,:))*Hi)', [],[],...
                    Gi,gi,LBi,UBi,lpoptions);

            if tt==1
                ZRA(ii,:,tt) = ZZ(ii,:,tt);
            else
                ZRA(ii,:,tt) = (1/tt)*((tt-1)*ZRA(ii,:,tt-1)...
                                                  +ZZ(ii,:,tt));
            end

        end

        for ii=1:NN
            
            Hi=HH(:,NN*(ii-1)+1:NN*ii);
            
            grad_ii=Hi*ZZ(ii,:,tt)'-bi;
            lam(ii,:,tt+1)=vv(ii,:)'+alpha_t*grad_ii;
        end

	lamAvg = mean(lam(:,:,tt));
    
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
        
        if mod(tt-99,MAXITERS/100)==0
            fprintf((repmat('\b', 1,msglnt)));
        end
        
    end
    
    tt=MAXITERS;
    fprintf((repmat('\b', 1,msglnt)));
    fprintf('progress: %d%%\n',100);
    
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
%        Gi=GG(:,NN*(ii-1)+1:NN*ii);
%        gi=gg((ii-1)*NN+1 : ii*NN);
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
end
