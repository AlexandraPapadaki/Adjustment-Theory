 clc
 clear all 
 close all
 
 %Load data
 e11=load('epoch2011.txt.');
 e12=load('epoch2012.txt.');
 
 %Epoch 2011
 x11=mean(e11);
 v11=x11-e11;
 s0_11=std(e11);
 r11=9;
 
 %Epoch 2012
 x12=mean(e12);
 v12=x12-e12;
 s0_12=std(e12);
 r12=24;
 
 %Plot residuals
 subplot(2,1,1);
 bar(v11); title('2011');
 subplot(2,1,2);
 bar(v12); title('2012');
 
 % Comments for the ploted residuals
 %No ordered residuals..although the look like so, because these are
 %measurements in time series and that happens...but not in that case
 %The problem would be if the graph would be according to time..then if the
 %show repeatness...that'sthe problem.
 %In 2012 there are 0...which is weard.
 %If we assume that std is around 0.0017 and we have residuals of 0.004, we
 %might have a blunder if the residuals were 3 to 4 times the std.
 %This is not the case but lets assume so and  check it with t-distribution
 %%
 
 %Plot F-distribution
 x=linspace(0,6,10^3);
 figure
 plot(x,fpdf(x,r11,r12)); title('F-distribution');
 
 %Confidence limits of F-distribution
 af_95=s0_11/s0_12*sqrt(1/finv(0.975,r11,r12))
 bf_95=s0_11/s0_12*sqrt(finv(0.975,r12,r11))
 
 %Chi^2 plot
 figure
 x=linspace(0,50,10^3);
 plot(x,chi2pdf(x,r11)),'b',x,chi2pdf(x,r12),'r'; 
 title('Chi^2 Dist');
 
 %Confidence limits of Chi^2-distribution 
 %!!!!!!!!!!!!!!!!!Now we will check if EPOCH belongs to the same theoretical reference std!!!!!!!!!!!!!
 achi2_11_95=s0_11*sqrt(r11/chi2inv(0.975,r11))
 bchi2_11_95=s0_11*sqrt(r11/chi2inv(0.025,r11))
 %Epoch 2011 is not so different. 
 %It belongs to the given theoritical reference std, it's inside the interval with a probability of 95%
 
 achi2_12_95=s0_12*sqrt(r12/chi2inv(0.975,r12))
 bchi2_12_95=s0_12*sqrt(r12/chi2inv(0.025,r12))
 %Epoch 2012 IS different. 
 %It DOESNOT belong to the given theoritical reference std, it's NOT inside the interval with a probability of 95%
 %So either we have a BLUNDER or we have done something wrong.
 %Soooo we check it by t-dist!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
%% chi2inv = critical value 
% Βρίσκει την critical value ώστε να ικανοποιείται η probability 95%, 
% δηλαδή να έχει αριστερά της 95% και δεξιά της 5%, στη γραφική παράσταση.
% Δες γραφική παράσταση στις σημειώσεις.
 
 

 
 
 