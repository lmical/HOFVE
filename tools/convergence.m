%test="../TestSteadyVortex";
%test="../STEADY_VORTEX";
%test="../SteadyVortex1";
%test="../SteadyVortex2";
%test="../SteadyVortex3";
%test="../TestSmoothSteadyVortex";%
%test="../TestSmoothSteadyVortex/nGPs4";%
%test="../TestSmoothSteadyVortex/WENO3";
%test="../TestSmoothVortex";
%test="../TestSmoothVortex/RK65";

test="../UnsteadyVortex/DeC5Test311";%
test="../UnsteadyVortex/stuMpDeCTest311";
%test = "../UnsteadyVortex/jacobiMPDeC5Test311";
test ="../TestLakeAtRest/mPDeC5NotWB";

Ns=2.^[2:10];
errors = zeros(length(Ns),3);
fid = fopen(sprintf("%s/convergence.tex",test),'w');
        
fprintf("  N    Error h    Order h  Error u    Order u  Error v    Order v\n")
fprintf(fid,"  N   & Error h  &  Order h & Error u  &  Order u & Error v  &  Order v \\ \n");
for i=1:length(Ns)
    N=Ns(i);
    try
        fileID = fopen(sprintf("%s/ErrorL1_%04d_%04d.dat",test,N,N),'r');
        mario = fscanf(fileID, ['%f %f %f']);
        errors(i,:)=mario;
        fclose(fileID);
    catch
        errors(i,:)=NaN;
    end
    if i>1
        order = -log(errors(i,:)./errors(i-1,:))/log(2);
    else
        order =zeros(3,1);
    end
    fprintf("%4d   %1.4e  %1.4f    %1.4e  %1.4f    %1.4e  %1.4f\n",...
        N,errors(i,1),order(1),errors(i,2),order(2),errors(i,3),order(3))
    fprintf(fid, " %4d  &   %1.3e  &  %1.3f  &  %1.3e & %1.3f  &  %1.3e & %1.3f \\ \n",...
        N,errors(i,1),order(1),errors(i,2),order(2),errors(i,3),order(3));
end
fclose(fid);

fig=figure();
alignedIdx=5;
alignedVar=1;
loglog(Ns,errors(:,1),"-*",'DisplayName',"Error h")
hold on
loglog(Ns,errors(:,2),"-+",'DisplayName',"Error u")
loglog(Ns,errors(:,3),"-x",'DisplayName',"Error v")
loglog(Ns,Ns.^-1*errors(alignedIdx,alignedVar)*Ns(alignedIdx),"--",'DisplayName',"1st order")
loglog(Ns,Ns.^-2*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^2,"--",'DisplayName',"2nd order")
loglog(Ns,Ns.^-3*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^3,"--",'DisplayName',"3rd order")
loglog(Ns,Ns.^-4*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^4,"--",'DisplayName',"4th order")
loglog(Ns,Ns.^-5*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^5,"--",'DisplayName',"5th order")
legend('location','best')
ylabel("Error")
xlabel("Grid points")
f = gcf;
exportgraphics(f,strcat(test,'/convergence.png'),'Resolution',600)
saveas(fig,strcat(test,"/convergence.pdf"))

