[NDEerr, NDEdof, NDEdt, NDEsetup, NDErun, NDEM_vec, NDEp_vec] = parse_output_file("euler2d_dt0p0001NodalDiagE_M124812_p678_C_t0p1_rhomin0p9.txt");
[TPSSerr, TPSSdof, TPSSdt, TPSSsetup, TPSSrun, TPSSM_vec, TPSSp_vec] = parse_output_file("euler2d_dt0p0001NodalTPSSLGL_M1246_p678_C_t0p1_rhomin0p9.txt");
[NOerr, NOdof, NOdt, NOsetup, NOrun, NOM_vec, NOp_vec] = parse_output_file("euler2d_dt0p0001NodalOmega_M124812_p678_C_t0p1_rhomin0p9.txt");
[MTerr, MTdof, MTdt, MTsetup, MTrun, MTM_vec, MTp_vec] = parse_output_file("euler2d_dt0p0001ModalTensor_M124812_p678_C_t0p1_rhomin0p9.txt");
%remove the M=2 data point, ie remove first row of the matrices
NDEerr(1,:) = [];
NDErun(1,:) = [];
TPSSerr(1,:) = [];
TPSSrun(1,:) = [];
NOerr(1,:) = [];
NOrun(1,:) = [];
MTerr(1,:) = [];
MTrun(1,:) = [];

figure(4)
subplot(2,3,1)
loglog(TPSSerr(:,1),TPSSrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
hold on
loglog(MTerr(:,1),MTrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
loglog(NOerr(:,1),NOrun(:,1),'-o','LineWidth',2,'MarkerSize',8)
loglog(NDEerr(:,1),NDErun(:,1),'-o','LineWidth',2,'MarkerSize',8)
title("Isentropic Vortex, p=6 comparisons",'Interpreter','latex','FontSize',16)
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
title("Isentropic Vortex, p=7 comparisons",'Interpreter','latex','FontSize',16)
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
title("Isentropic Vortex, p=8 comparisons",'Interpreter','latex','FontSize',16)
xlabel("H-norm error",'Interpreter','latex','FontSize',16)
ylabel("run-time (s)",'Interpreter','latex','FontSize',16)
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',16)

Lgnd=legend("TPSS-LGL", "SBP-CC", "SBP-Omega", "SBP-E",'Interpreter','latex','FontSize',16);
Lgnd.Position(1) = -0.025;
Lgnd.Position(2) = 0.75;