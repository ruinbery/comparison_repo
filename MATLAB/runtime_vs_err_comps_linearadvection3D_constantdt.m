[NDEerr, NDEdof, NDEdt, NDEsetup, NDErun, NDEM_vec, NDEp_vec] = parse_output_file("NodalDiagE_M12345_p12345678_CFL0p13Domega8dt0p001.txt");
[TPSSerr, TPSSdof, TPSSdt, TPSSsetup, TPSSrun, TPSSM_vec, TPSSp_vec] = parse_output_file("NodalTPSSLGL_M12345_p12345_CFL0p13Domega8dt0p001.txt");
[NOerr, NOdof, NOdt, NOsetup, NOrun, NOM_vec, NOp_vec] = parse_output_file("NodalOmega_M12345678_p12345_CFL0p13Domega8dt0p001.txt");
[MTerr, MTdof, MTdt, MTsetup, MTrun, MTM_vec, MTp_vec] = parse_output_file("ModalTensor_M12345678_p12345_CFL0p13Domega8dt0p001.txt");
%remove the M=2 data point, ie remove first row of the matrices
NDEerr(1,:) = [];
NDErun(1,:) = [];
TPSSerr(1,:) = [];
TPSSrun(1,:) = [];
NOerr(1,:) = [];
NOrun(1,:) = [];
MTerr(1,:) = [];
MTrun(1,:) = [];

figure(3)
subplot(2,3,1)
loglog(TPSSerr(:,1),TPSSrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,1),MTrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,1),NOrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,1),NDErun(:,1),'-o','LineWidth',2,'MarkerSize',8)
title("LinAdvec 3D, p=1 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,3,2)
loglog(TPSSerr(:,2),TPSSrun(:,2),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,2),MTrun(:,2),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,2),NOrun(:,2),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,2),NDErun(:,2),'-o','LineWidth',2,'MarkerSize',8)
title("LinAdvec 3D, p=2 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,3,3)
loglog(TPSSerr(:,3),TPSSrun(:,3),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,3),MTrun(:,3),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,3),NOrun(:,3),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,3),NDErun(:,3),'-o','LineWidth',2,'MarkerSize',8)
title("LinAdvec 3D, p=3 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,3,4)
loglog(TPSSerr(:,4),TPSSrun(:,4),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,4),MTrun(:,4),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,4),NOrun(:,4),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,4),NDErun(:,4),'-o','LineWidth',2,'MarkerSize',8)
title("LinAdvec 3D, p=4 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,3,5)
loglog(TPSSerr(:,5),TPSSrun(:,5),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,5),MTrun(:,5),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,5),NOrun(:,5),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,5),NDErun(:,5),'-o','LineWidth',2,'MarkerSize',8)
title("LinAdvec 3D, p=5 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

Lgnd=legend("TPSS-LGL", "SBP-CC", "SBP-Omega", "SBP-E",'Interpreter','latex','FontSize',16);
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.4;