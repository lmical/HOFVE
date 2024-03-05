%test="../TestSteadyVortex";
%test="../STEADY_VORTEX";
%test="../TestSmoothSteadyVortex";
%test="../TestSmoothVortex";
%test="../TestSmoothWaves";
test="../TestOscillations/DeC5";
test="../TestOscillations/mpDeC5";
test="../UnsteadyVortex/DeC5Test311";
%test="../UnsteadyVortex/mPDeC5Test311";
%test = "../UnsteadyVortex/jacobiMPDeC5Test311";


delimiterIn   = ' ';
headerlinesIn = 1;
Ns=2.^[2:10];
%Ns=[25 50 100 200 300 400 500 600];
errors = zeros(length(Ns),3);

fprintf("  N    Error h    Order h  Error u    Order u\n")
for i=1:length(Ns)
    N=Ns(i);
    try
        filename = sprintf("%s/SOLUTION_%d.dat",test,N);
        mydata_solution = importdata(filename,delimiterIn,headerlinesIn);
        for nVar=1:3
            u = mydata_solution.data(:,nVar);
    
            U{i,nVar} = reshape(u,[N,N]);
            if i>1
                errors(i,nVar)=computeError(U{i,nVar},U{i-1,nVar});
            end
        end
    catch
        errors(i,:)=NaN;
    end
    if i>2
        order = -log(errors(i,:)./errors(i-1,:))/log(2);
    else
        order =zeros(3,1);
    end
    fprintf("%4d   %1.3e  %1.3f    %1.3e  %1.3f\n",N,errors(i,1),order(1),errors(i,2),order(2))
end

fig=figure();
alignedIdx=3;
alignedVar=2;
loglog(Ns,errors(:,1),"-*",'DisplayName',"Error h")
hold on
loglog(Ns,errors(:,2),"-+",'DisplayName',"Error u")
loglog(Ns,errors(:,3),"-x",'DisplayName',"Error v")
loglog(Ns,Ns.^-1*errors(alignedIdx,alignedVar)*Ns(alignedIdx),"--",'DisplayName',"1st order")
loglog(Ns,Ns.^-2*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^2,"--",'DisplayName',"2nd order")
loglog(Ns,Ns.^-3*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^3,"--",'DisplayName',"3rd order")
loglog(Ns,Ns.^-4*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^4,"--",'DisplayName',"4th order")
loglog(Ns,Ns.^-5*errors(alignedIdx,alignedVar)*Ns(alignedIdx)^4,"--",'DisplayName',"5th order")
legend('location','best')
ylabel("Error")
xlabel("Grid points")
f = gcf;
exportgraphics(f,strcat(test,'/convergence.png'),'Resolution',600)
saveas(fig,strcat(test,"/convergence.pdf"))

function err=computeError(U1,U2)
    if size(U1,1)~= 2*size(U2,1)
        error("dimension wrong")
    end
    UHelp = zeros(size(U2));
    for i=1:size(U2,1)
        for j=1:size(U2,2)
            UHelp(i,j)=mean(U1(2*i-1:2*i,2*j-1:2*j),'all');
        end
    end
%    figure()
%    surf(UHelp-U2)
    err=mean(abs(UHelp-U2),'all');
    %title(str(size(U2,1)))
end