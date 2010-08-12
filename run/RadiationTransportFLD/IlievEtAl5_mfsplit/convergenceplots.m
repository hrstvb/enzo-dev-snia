
% create plots showing convergence of errors

rdata_16 = load('nx16/radiusdata.txt');
rdata_32 = load('nx32/radiusdata.txt');
rdata_64 = load('nx64/radiusdata.txt');
rdata_128 = load('nx128/radiusdata.txt');

h = figure(1);
plot(rdata_16(:,1),rdata_16(:,2)./rdata_16(:,4),'k:',...
     rdata_32(:,1),rdata_32(:,2)./rdata_32(:,4),'b--',...
     rdata_64(:,1),rdata_64(:,2)./rdata_64(:,4),'r-',...
     rdata_128(:,1),rdata_128(:,2)./rdata_128(:,4),'g-.','LineWidth',2);
xlabel('t/t_{rec}','FontSize',14); ylabel('r_I/r_S','FontSize',14); grid on;
title('Convergence in I-front Position','FontSize',16);
hl = legend('16^3 mesh','32^3 mesh','64^3 mesh','128^3 mesh',2);
set(hl,'FontSize',14);  set(gca,'FontSize',14);  axis([0,4.5,0,2.5])
saveas(h,'rconv.fig');
saveas(h,'rconv.png');
saveas(h,'rconv.eps','epsc');

h = figure(4);
plot(rdata_16(:,1),abs(rdata_16(:,2)-rdata_128(:,2))./rdata_16(:,4),'k:',...
     rdata_32(:,1),abs(rdata_32(:,2)-rdata_128(:,2))./rdata_32(:,4),'b--',...
     rdata_64(:,1),abs(rdata_64(:,2)-rdata_128(:,2))./rdata_64(:,4),'r-',...
     'LineWidth',2);
xlabel('t/t_{rec}','FontSize',14); ylabel('|r_I-r_I^*|/r_S','FontSize',14); grid on;
title('Convergence in I-front Position','FontSize',16);
hl = legend('16^3 mesh','32^3 mesh','64^3 mesh');
set(hl,'FontSize',14);  set(gca,'FontSize',14);  %axis([0,4.5,0,3.0])
saveas(h,'rerr.fig');
saveas(h,'rerr.png');
saveas(h,'rerr.eps','epsc');


clear rdata_16 rdata_32 rdata_64 rdata_128;



vdata_16 = load('nx16/velocitydata.txt');
vdata_32 = load('nx32/velocitydata.txt');
vdata_64 = load('nx64/velocitydata.txt');
vdata_128 = load('nx128/velocitydata.txt');

h = figure(2);
plot(vdata_16(:,1),vdata_16(:,2),'k:',...
     vdata_32(:,1),vdata_32(:,2),'b--',...
     vdata_64(:,1),vdata_64(:,2),'r-',...
     vdata_128(:,1),vdata_128(:,2),'g-.','LineWidth',2);
xlabel('t/t_{rec}','FontSize',14); ylabel('v_I / (r_S/t_{rec})','FontSize',14); grid on;
title('Convergence in I-front Velocity','FontSize',16);
hl = legend('16^3 mesh','32^3 mesh','64^3 mesh','128^3 mesh');
set(hl,'FontSize',14);  set(gca,'FontSize',14);  axis([0,4.5,-0.2,1.2])
saveas(h,'vconv.fig');
saveas(h,'vconv.png');
saveas(h,'vconv.eps','epsc');

h = figure(5);
plot(vdata_16(:,1),abs(vdata_16(:,2)-vdata_128(:,2)),'k:',...
     vdata_32(:,1),abs(vdata_32(:,2)-vdata_128(:,2)),'b--',...
     vdata_64(:,1),abs(vdata_64(:,2)-vdata_128(:,2)),'r-','LineWidth',2);
xlabel('t/t_{rec}','FontSize',14); ylabel('|v_I-v_I^*| / (r_S/t_{rec})','FontSize',14); grid on;
title('Convergence in I-front Velocity','FontSize',16);
hl = legend('16^3 mesh','32^3 mesh','64^3 mesh');
set(hl,'FontSize',14);  set(gca,'FontSize',14);  %axis([0,4.5,-0.2,1.2])
saveas(h,'verr.fig');
saveas(h,'verr.png');
saveas(h,'verr.eps','epsc');

clear vdata_16 vdata_32 vdata_64 vdata_128;




ndata_16 = load('nx16/nfracdata.txt');
ndata_32 = load('nx32/nfracdata.txt');
ndata_64 = load('nx64/nfracdata.txt');
ndata_128 = load('nx128/nfracdata.txt');

h = figure(3);
plot(ndata_16(:,1),ndata_16(:,2),'k:',...
     ndata_32(:,1),ndata_32(:,2),'b--',...
     ndata_64(:,1),ndata_64(:,2),'r-',...
     ndata_128(:,1),ndata_128(:,2),'g-.','LineWidth',2);
xlabel('t/t_{rec}','FontSize',14); ylabel('1-x','FontSize',14); grid on;
title('Convergence in Total Neutral Fraction','FontSize',16);
hl = legend('16^3 mesh','32^3 mesh','64^3 mesh','128^3 mesh');
set(hl,'FontSize',14);  set(gca,'FontSize',14);  axis([0,4.5,0.5,1.1])
saveas(h,'nconv.fig');
saveas(h,'nconv.png');
saveas(h,'nconv.eps','epsc');

h = figure(6);
plot(ndata_16(:,1),abs(ndata_16(:,2)-ndata_128(:,2)),'k:',...
     ndata_32(:,1),abs(ndata_32(:,2)-ndata_128(:,2)),'b--',...
     ndata_64(:,1),abs(ndata_64(:,2)-ndata_128(:,2)),'r-','LineWidth',2);
xlabel('t/t_{rec}','FontSize',14); ylabel('|x^*-x|','FontSize',14); grid on;
title('Convergence in Total Neutral Fraction','FontSize',16);
hl = legend('16^3 mesh','32^3 mesh','64^3 mesh');
set(hl,'FontSize',14);  set(gca,'FontSize',14);  %axis([0,4.5,0.5,1.1])
saveas(h,'nerr.fig');
saveas(h,'nerr.png');
saveas(h,'nerr.eps','epsc');

clear ndata_16 ndata_32 ndata_64 ndata_128;





% end script