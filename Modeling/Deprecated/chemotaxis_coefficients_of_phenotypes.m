dCWBias=0.001; CWBias=0.001:dCWBias:0.99;
MFTmodel = MFTmodelInitialize(CWBias);
figure
subplot(3,1,1)
plot(CWBias, MFTmodel.ki)
xlabel('TB')
ylabel('\chi (\mum^2/s)')
axis tight
subplot(3,1,2)
plot(CWBias, MFTmodel.Diff)
xlabel('TB')
ylabel('\mu (\mum^2/s)')
axis tight
subplot(3,1,3)
plot(CWBias, MFTmodel.ki ./ MFTmodel.Diff)
xlabel('TB')
ylabel('\chi/\mu')
axis tight
set(gcf, 'PaperPosition', [0 0 1 2] * 8);
print(gcf, 'chemotaxis_coefficients_of_phenotypes.jpg', '-djpeg', '-r300', '-opengl')
close(gcf)