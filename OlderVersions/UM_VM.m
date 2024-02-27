function [UM,VM] = UM_VM(AXH_J,AYH_J,q,BXH_J,BYH_J,BXH_W,BYH_W,TAU,TAUM1,Delta,L,V_INF,NM,CXH_M,CYH_M,CV,k)
    
    for h = 1:NM
        term1 = sum(AXH_J(h,:).*q');
        term2 = TAU*sum(BXH_J(h,:));
        term3 = dot(V_INF,[1 0]);
        term4 = L*(TAUM1-TAU)/Delta*BXH_W(h);
        term5 = 0;
        if NM >= 1
            for m = 1:NM
                if m ~= h
                    term5 = term5 + CXH_M(h,m)*CV(m);
                end
            end
        end
        UM(h) = term1 + term2 + term3 + term4 + term5;
    end

    for h = 1:NM
        term1 = sum(AYH_J(h,:).*q');
        term2 = TAU*sum(BYH_J(h,:));
        term3 = dot(V_INF,[0 1]);
        term4 = L*(TAUM1-TAU)/Delta*BYH_W(h);
        term5 = 0;
        if NM >= 1
            for m = 1:NM
                if m ~= h
                    term5 = term5 + CYH_M(h,m)*CV(m);
                end
            end
        end
        VM(h) = term1 + term2 + term3 + term4 + term5;
    end

end

