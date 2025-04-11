function Param=input_Boas()
%% RheolParams reading file
ifp = fopen("RheolParams.dat", "r");
bifpar=fscanf(ifp,"%f");
fgetl(ifp);
cpar=fscanf(ifp,"%f");
fgetl(ifp);
viscpar=fscanf(ifp,"%f");
fgetl(ifp);
nitmax=fscanf(ifp,"%f",1);
tol=fscanf(ifp,"%f",1);
omega=fscanf(ifp,"%f",1);
fgetl(ifp);
nitmax1=fscanf(ifp,"%f",1);
qtol=fscanf(ifp,"%f",1);
hdtol=fscanf(ifp,"%f",1);
fgetl(ifp);
optw=fscanf(ifp,"%f",1);
optlam=fscanf(ifp,"%f",1);
fgetl(ifp);
constvisc=fscanf(ifp,"%f",1);
vplas=fscanf(ifp,"%f",1);
mcv=fscanf(ifp,"%f",1);
fgetl(ifp);
consthd=fscanf(ifp,"%f",1);
fgetl(ifp);
varyviscosity=fscanf(ifp,"%f",1);
fgetl(ifp);
phaseseparation=fscanf(ifp,"%f",1);
fclose(ifp);
%% conversiion
facfp = pi*1333./128./0.01*60./1.e6;
mcvcorr =(92./mcv)^0.33333;
%%
ParNames = who();
for i = 1:numel(ParNames)
    eval(['Param.' ParNames{i} '=' ParNames{i} ';' ]);
end
end