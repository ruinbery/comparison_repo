[NDEerr, NDEdof, NDEdt, NDEsetup, NDErun, NDEM_vec, NDEp_vec] = parse_output_file("NodalDiagE_M2481216_p12345678_CFL0p1omega8dt0p005dp.txt");
[TPSSerr, TPSSdof, TPSSdt, TPSSsetup, TPSSrun, TPSSM_vec, TPSSp_vec] = parse_output_file("NodalTPSSLGL_M2481216_p12345678_CFL0p1omega8dt0p005dp.txt");
[NOerr, NOdof, NOdt, NOsetup, NOrun, NOM_vec, NOp_vec] = parse_output_file("NodalOmega_M2481216_p12345678_CFL0p1omega8dt0p005dp.txt");
[MTerr, MTdof, MTdt, MTsetup, MTrun, MTM_vec, MTp_vec] = parse_output_file("ModalTensor_M2481216_p12345678_CFL0p1omega8dt0p005dp.txt");
%remove the M=2 data point, ie remove first row of the matrices
NDEerr(1,:) = [];
NDErun(1,:) = [];
TPSSerr(1,:) = [];
TPSSrun(1,:) = [];
NOerr(1,:) = [];
NOrun(1,:) = [];
MTerr(1,:) = [];
MTrun(1,:) = [];

figure(1)
subplot(2,4,1)
loglog(TPSSerr(:,1),TPSSrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,1),MTrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,1),NOrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,1),NDErun(:,1),'-o','LineWidth',2,'MarkerSize',8)
title("p=1 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,4,2)
loglog(TPSSerr(:,2),TPSSrun(:,2),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,2),MTrun(:,2),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,2),NOrun(:,2),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,2),NDErun(:,2),'-o','LineWidth',2,'MarkerSize',8)
title("p=2 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)


subplot(2,4,3)
loglog(TPSSerr(:,3),TPSSrun(:,3),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,3),MTrun(:,3),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,3),NOrun(:,3),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,3),NDErun(:,3),'-o','LineWidth',2,'MarkerSize',8)
title("p=3 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,4,4)
loglog(TPSSerr(:,4),TPSSrun(:,4),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,4),MTrun(:,4),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,4),NOrun(:,4),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,4),NDErun(:,4),'-o','LineWidth',2,'MarkerSize',8)
title("p=4 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,4,5)
loglog(TPSSerr(:,5),TPSSrun(:,5),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,5),MTrun(:,5),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,5),NOrun(:,5),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,5),NDErun(:,5),'-o','LineWidth',2,'MarkerSize',8)
title("p=5 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,4,6)
loglog(TPSSerr(:,6),TPSSrun(:,6),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,6),MTrun(:,6),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,6),NOrun(:,6),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,6),NDErun(:,6),'-o','LineWidth',2,'MarkerSize',8)
title("p=6 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,4,7)
loglog(TPSSerr(:,7),TPSSrun(:,7),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,7),MTrun(:,7),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,7),NOrun(:,7),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,7),NDErun(:,7),'-o','LineWidth',2,'MarkerSize',8)
title("p=7 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

subplot(2,4,8)
loglog(TPSSerr(:,8),TPSSrun(:,8),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,8),MTrun(:,8),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,8),NOrun(:,8),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,8),NDErun(:,8),'-o','LineWidth',2,'MarkerSize',8)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

title("p=8 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)

Lgnd=legend("TPSS-LGL", "SBP-CC", "SBP-Omega", "SBP-E",'Interpreter','latex','FontSize',16);
Lgnd.Position(1) = 0.0;
Lgnd.Position(2) = 0.4;