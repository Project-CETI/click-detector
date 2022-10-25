function [mf,Sigl,C_inds,Reduced,Cut]=segment_IPIs(IPI)

      id_save={};
      Reduced=[0 0]; Sigl=[]; mf=[]; ms=[]; C_inds={}; Valid_flag=0;
      [input,I] = sort(1e3*IPI,'descend');
      
      q=0; DIF=1;
      while DIF>0 && length(input)>4
          L=length(input);
          [input, removed_ind] = rmoutliers(input);
          DIF=L-length(input);
          if DIF>0
              q=q+1;
              removed=find(removed_ind);
              id_save(q)={removed+ones(1,length(removed))*(q-1)};
              removed=[];
          end
      end
 
      if ~isempty(id_save)
          id_save_convert=cell2mat(id_save);
          I(id_save_convert)=[];
      end
      
      Original_input=input;     
      
       x=[1:length(input)];
       P = polyfit(x,input,1);
       yfit = polyval(P,x);
       if length(yfit)>1
           m=yfit(1)-yfit(2);
       else
            m=yfit(1);
       end
       
%       figure; plot(input,'o','Linewidth',2); hold on; plot(x,yfit,'r-.','Linewidth',2);
%       xlabel('Transient index','Fontsize',14); ylabel('IPI','Fontsize',14); ylim([0 18]); grid on;
% % % 
%       Rn=rand(1,56)*18e-3;
%       figure; histogram(input,10)
%       figure; histogram(Rn,10)
      
      
      
      
    figure; subplot(3,1,1);
    samples=3;
    input=[ input zeros(1,2*samples)];
    for q=1:length(input)-samples
      x=[q:q+samples];
      P = polyfit(x,input(q:q+samples),1);
      yfit = polyval(P,x);
      ms(q)=yfit(1)-yfit(2);
      plot(input,'x','Linewidth',2); grid on; xlabel('Click #'); ylabel('IPI [ms]');
%       xlim([0 length(IPI)]);
%       hold on; plot(x,yfit,'r-.','Linewidth',2);
    end


   subplot(3,1,2);
   plot(ms,'-x','Linewidth',2); grid on; xlabel('Click #'); ylabel('shape');
   
       Ds=diff(ms);
    subplot(3,1,3);
   plot(Ds,'-x','Linewidth',2); grid on; xlabel('Click #'); ylabel('shape');


      Count=0; Valid_flag=0; id_save=[];
      ms_orig=ms(2*samples+1:2*samples+length(Original_input));
      vali_inds=find(ms_orig<0.05);
      Follow=diff(vali_inds); 
      for i=1:length(Follow)
         if Follow(i)==1
             Count=Count+1;
             id_save=[id_save i];
         elseif Count>4
             Valid_flag=1;
             break
         else
             Count=0;
             id_save=[];
         end
      end
      
    if Valid_flag==1
      Ind_tmp=vali_inds(id_save);
      G=[Ind_tmp Ind_tmp(end)+1 Ind_tmp(end)+2 Ind_tmp(end)+3 Ind_tmp(end)+4 Ind_tmp(end)+5];
      C_inds={I(G)};
      Cut= C_inds;
       xval=[1:length(G)];
       P = polyfit(xval,input(G),1);
       yfit = polyval(P,xval);
       mval=yfit(1)-yfit(2);
       Sigval=std(Original_input(G));
       Reduced=[mval Sigval];
    else

    %%
    % subplot(3,1,3);
            T_inds={}; C_inds={}; count=0;
            [ps,los] =findpeaks(ms,'MinPeakHeight',0.15,'MinPeakDistance',5);  % Apply instantaneous energy detector (find peaks)          
            los=los(los>8);
            Padding=2*samples;
            start=Padding+1; 
            for i=1:length(los)
                neib=los(i)-4;
                if neib<samples
                    tmp=ms(1:neib+samples);
                elseif neib+samples>length(input)
                    tmp=ms(neib-samples:length(input));
                else
                    tmp=ms(neib-samples:neib+samples);
                end
                G(i)=sum(tmp<0.15);
                if G(i)>1
                    count=count+1;
    %                plot([start:los(i)+2],input([start:los(i)+2]),'x','Linewidth',2); grid on; xlabel('Click #'); ylabel('shape');
    %                hold on; xlim([0 length(IPI)]);
                    start=los(i)-samples-2-Padding;
                    End=los(i)-1-Padding;
                    indices=[los(i)-samples-G(i)-Padding:los(i)-1-Padding];                  
                    T_inds(count)={indices};
                    C_inds(count)={I(indices)};
                    Sigl(count)=std(Original_input(indices));
                end
            end

            input=Original_input;
            index=cellfun(@isempty, T_inds) == 0;
            T_inds=T_inds(index);
            for e=1:size(T_inds,2)
                Test=input(cell2mat(T_inds(e)));
                x=[1:length(Test)];
                P = polyfit(x,Test,1);
                yfit = polyval(P,x);

                mf(e)=yfit(1)-yfit(2);
                Test=[];
            end
    end   
         if ~isempty(mf)
             FI=find(mf==min(mf));
             Reduced=[mf(FI(1)) Sigl(FI(1))];
             Cut=T_inds(FI(1));
         else
             Reduced=[m std(input)];
             mf=m; Sigl=std(input);
             C_inds={I};
             Cut=C_inds;
         end
end


        
     