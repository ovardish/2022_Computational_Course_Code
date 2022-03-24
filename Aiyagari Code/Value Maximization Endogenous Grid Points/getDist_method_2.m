function [dist,trans] = getDist_continuous(param,grid,kpol,trans_only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to solve for the distribution of HHs with continuous policies
%% for reference see Young (2010, JEDC)
%% inputs:
%%       - grid: structure containing grids
%%       - kopt: indices of optimal capital policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding  a_i and a_j in between where my kpol(s, a) is lying!
% Function to find left and right points capital choice kpol(ial, iar)
  %[ial, iar, varphi]= locate(kpol,grid.k);
  tiny=0.1^(14);
% % %  Manual way to find upper and lower bounds of it.
% % % find index below capital policy
     ial = ones(param.nz,param.nkap);
     for i=2:param.nkap
         ial(kpol>=grid.k(i))=i;
     end
     ial = min(ial,param.nkap-1);
     iar = ial + 1;
% %     
% % % find weights attached to point below/above policy
     varphi = (kpol-grid.k(ial))./(grid.k(ial+1)-grid.k(ial));
     varphi = min(varphi,1); varphi = max(varphi,0); % just to be safe, should not be binding
% %     wbelow= 1-wabove;

   
    
    
    dist = NaN(param.nz,param.nkap);
   %if (trans_only~=1)

     lambda     = ones(param.nz,param.nkap)/(param.nz*param.nkap);
     lambda_new = zeros(param.nz,param.nkap);
     err=1;
 while err>tiny
         
     for j=1:param.nz                % Current productivity state (in t)
         for i=1:param.nkap       % Current asset state (in t)
          %  for j_p = 1:param.nz    % Future productivity state (in t+1)       
           %   Asset state (in t+1)           
           lambda_new(:,ial(j, i))=  lambda_new(:,ial(j, i))+grid.Pz(j,:)'.*        varphi(j, i)*lambda(j,i); 
           lambda_new(:,iar(j, i))=  lambda_new(:,iar(j, i))+grid.Pz(j,:)'.*(1.0d0-varphi(j, i))*lambda(j,i); 
           % end
         end
       end
   err=max(abs(lambda_new(:)-lambda(:)));
   lambda=lambda_new;
   lambda_new = lambda_new.*0;
   err
  end
   %end

    % reshape distribution to two dimensions
    dist=lambda;
    end  
    
    
 