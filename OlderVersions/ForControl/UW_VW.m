function [UW,VW] = UW_VW(AXW_J,AYW_J,q,BXW_J,BYW_J,TAU,V_INF,NM,CXW_M,CYW_M,CV)

    term1 = sum(AXW_J.*q');
    term2 = TAU*sum(BXW_J);
    term3 = dot(V_INF,[1 0]);
    term4 = 0; % VERIFY IF THIS IS REALLY ZERO
    term5 = 0;
    if NM >= 1
        for m = 1:NM
            term5 = term5 + CXW_M(m)*CV(m);
        end
    end
    UW = term1 + term2 + term3 + term4 + term5;

    term1 = sum(AYW_J.*q');
    term2 = TAU*sum(BYW_J);
    term3 = dot(V_INF,[0 1]);
    term4 = 0; % VERIFY IF THIS IS REALLY ZERO
    term5 = 0;
    if NM >= 1
        for m = 1:NM
           term5 = term5 + CYW_M(m)*CV(m);
        end
    end
    VW = term1 + term2 + term3 + term4 + term5;
    
%     fprintf("%2.3f %2.3f %2.3f %2.3f %2.3f\n",term1,term2,term3,term4,term5);
%     pause();
    
end

