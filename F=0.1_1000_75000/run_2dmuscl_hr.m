% Time and Grid Parameters
tend=[5580];                                       % Total Simulation Time for two storms 5580+3600+6120?
dx=0.4;                                            % Grid Spacing
rand('seed',1)                                      % Define Seed for Random Number Generator

% HR Model Parameters
F=[0.1];                                         % Fraction of Excess Stream Power In Entrainment need to figure out
a=[1000];                                           % Detachability of Original Soil (Used to determine rainsplash detachment rate) need to figure out
ad=[75000];                                        % Detachability of Deposited Sediment (Used to determine rainsplash detachment rate for deposited sediment) need to figure out

h0=0.33*[0.002];
%h0=0.33*raindiameter;                                    % Critical flow depth in HR model original guess
mtstar0=[3];                                        % Critical mass of deposited sediment needed to shield original soil

% Parameters for Interception Model (Rutter, 1972)
vegcover=[0];                                       % Fraction vegetaion cover
pi=0.35;                                            % Free throughfall Coefficient 
Si=2.2/1000;                                        % Canopy capacity (m)
Ki=5e-8;                                            % Canopy Drainage Rate Coefficient (m/s)
gi=3.9*1000;                                        % Canopy Drainage Rate Exponent (m^{-1})

% Parameters Needed To Compute Critical Stream Power
d50=1/1000*[0.22];                                   % Median Grain Size [m]
g=9.81;                                             % Gravity
rhof=1000;                                          % Density of Fluid
rhos=2600;                                          % Density of Sediment

% Input Particle Size Distribution
DS=1/1000*[0.001,0.01,0.04,0.1,0.2,0.3,0.5,1.0,2.0,3.0];
PDS=[0.09,0.04,0.075,0.125,0.15,0.1,0.09,0.15,0.09,0.09];

% Input Rainfall need to change for fitting two storms
rainfallin=dlmread('LasLomas_1min_dec.txt');                % Rainfall data for 1 min interval
raintime=(0:60:60*(length(rainfallin)-1))/60;                % Time associated with each entry
rint=60;
rain=rainfallin;
for i=1:length(rain)
    if (rain(i)<0)
        rain(i)=0;
    end
end
rnum=length(rain);

% Rainfall Diameter need to change for fitting two storms
raindiameter=dlmread('LasLomas_1min_dec_diameter.txt');        % Rainfall  diameter data from disdrometer


% Rainfall Velocity need to change for fitting two storms
rainvelocity=dlmread('LasLomas_1min_dec_vel.txt');        % Rainfall fall velocity data from disdrometer

% Initial Topography
topo=dlmread('LasLomas04mDEMnov.txt'); %Changed to small area
itopo=topo;

% Masks Needed for Internal Boundaries and Spatially Variable Erodibiity; if don't need?
solid=dlmread('LasLomas04mMasknov.txt');
solidin=solid;
erodibilitymask=ones(size(topo));

[nx,ny]=size(topo);
xpts=(0:dx:dx*(nx-1));
ypts=(0:dx:dx*(ny-1));
[xgrid,ygrid]=meshgrid(ypts,xpts);
xgridfine=xgrid;
ygridfine=ygrid;

% Mask That Defines Extent of Channel Network
contributingarea=dlmread('LasLomas04mContributingArea.txt');           % Contributing Area (Could Be Used to Define a Channel Mask)
cmask=dlmread('LasLomas04mChannelMask.txt');                           % Channel Mask Based on Contributing Area Threshold of 1000 m^3

% Hydraulic Roughness
hillsloperoughness=[0.03];                         % Minimum Manning n on hillslopes need change for second storm
hcfriction=[0.003];                                 % Parameter in depth-dependent Manning roughness model
depthdependentexponent=[-0.33];                     % Parameter in depth-dependent Manning roughness model
roughness=hillsloperoughness.*ones(nx,ny);

% Infiltration
corrl=dx;                                           % Correlation Length [m] Associated with Ks Field
%hfin=0.001*ones(nx,ny);                              % Wetting Front Capillary Pressure Head [m]
vinfin=0*ones(nx,ny);                              % Volume Infiltrated [m]
theta0in=0.05*ones(nx,ny);                           % Initial Volumetric Water Content [-]
thetasin=0.4*ones(nx,ny);                           % Saturated Soil Volumetric Water Content [-]

% Input Critical Stream Power As Function of Bed Slope% how to calculate this one?
omegacstar=dlmread('omegastarcriticalcorrected.txt');
[sx,sy]=gradient(topo,dx,dx);
sgrad=(sx.^2+sy.^2).^(1/2);
sgrad(sgrad>0.6745)=0.6745;
for i=1:nx
    for j=1:ny
        % Find non dimensional critical stream power based on gradient
        nondimensionalcriticalshearstress=omegacstar(max(floor(sgrad(i,j)/0.001),1));
        
        % Convert to critical stream power
        omegac(i,j)=nondimensionalcriticalshearstress*rhof*(g*(rhos-rhof)/rhof*d50)^(3/2);
    end
end

saveinputvec=[];

% Main
idum=0;
for bb=1:length(F)
    for hr1=1:length(a)
        for hr2=1:length(ad)
            
            % Spatially Variable Saturated Hydraulic Conductivity with Given Correlation Length same for two storms
            ksmin=0;
            [kstemp,xx,yy]=rsgeng2DEM(max([nx,ny]),max([nx,ny])*dx,0.5*corrl,0.5*corrl);          % Create Random Field with normal distribution of Ks values
            kstemp=kstemp;
            ks=kstemp;                                                                                 % Exponentiate to end up with lognormal distribution of Ks values
            ks=ks(1:nx,1:ny);
            ks(ks<ksmin)=ksmin;
            ksin=ks/(3600*1000);                                                                            % Convert from [mm]/[h] to [m]/[s]
            %ksin=exprnd(27,nx,ny)/(3600*1000);                                                              % Exponential Distribution of Ks Values
            ksin(cmask==1)=0;

            % Spatially Variable wetting front capillary pressure head with Given Correlation Length same for two storms
            hfmin=0;
            [hftemp,xx,yy]=rsgeng2DEMhf(max([nx,ny]),max([nx,ny])*dx,0.5*corrl,0.5*corrl);          % Create Random Field with normal distribution of Ks values
            hftemp=hftemp;
            hf=hftemp;                                                                                 % Exponentiate to end up with lognormal distribution of Ks values
            hf=hf(1:nx,1:ny);
            hf(hf<hfmin)=hfmin;
            hfin=hf;                                                                            % Convert from [mm]/[h] to [m]/[s]
            hfin(cmask==1)=0;

            % Write Initial Conditions to Input Files
            topoin=reshape(itopo',1,nx*ny);
            ksin=reshape(ksin',1,nx*ny);
            theta0in=reshape(theta0in',1,nx*ny);
            thetasin=reshape(thetasin',1,nx*ny);
            hfin=reshape(hfin',1,nx*ny);
            vinfin=reshape(vinfin',1,nx*ny);
            solidin=reshape(solidin',1,nx*ny);
            erodibilitymaskin=reshape(erodibilitymask',1,nx*ny);
            cmaskin=reshape(cmask',1,nx*ny);
            roughnessin=reshape(roughness',1,nx*ny);
            omegacin=reshape(omegac',1,nx*ny);
            dsin=DS;
            pdsin=PDS;
            
            dlmwrite('topoin',topoin,'delimiter','\t');
            dlmwrite('solidin',solidin,'delimiter','\t');
            dlmwrite('erodibilitymaskin',erodibilitymaskin,'delimiter','\t');
            dlmwrite('cmaskin',cmaskin,'delimiter','\t');
            dlmwrite('roughnessin',roughnessin,'delimiter','\t');
            dlmwrite('rain',rain,'delimiter','\t');
            dlmwrite('raindiam',raindiameter,'delimiter','\t');
            dlmwrite('rainvel',rainvelocity,'delimiter','\t');
            dlmwrite('ksin',ksin,'delimiter','\t');
            dlmwrite('theta0in',theta0in,'delimiter','\t');
            dlmwrite('thetasin',thetasin,'delimiter','\t');
            dlmwrite('hfin',hfin,'delimiter','\t');
            dlmwrite('vinfin',vinfin,'delimiter','\t');
            dlmwrite('particlesizein',dsin,'delimiter','\t');
            dlmwrite('particlepercentin',pdsin,'delimiter','\t');
            dlmwrite('omegacin',omegacin,'delimiter','\t');
            
            Jentrain=[0.5*1000*(5.5)^2/a(hr1)];
            dlmwrite('input',[nx ny dx F(bb) tend a(hr1) ad(hr2) h0 mtstar0 rnum rint length(dsin) Jentrain hcfriction depthdependentexponent pi Si Ki gi vegcover],'delimiter','\t');
            
            disp(num2str(idum))
            
            % Compile and Run C Code
            tic
            %!gcc-4.9 -fopenmp -o 2dmuscl_hr_omp 2dmuscl_hr_omp.c -lm
            %!./2dmuscl_hr_omp
            toc
            
            % Read Output Data From C Code if No Simulation Finished Successfully
            cfl=dlmread('./cflstatus');
            if (cfl==0)
                
                idum=idum+1;
                
                topo=dlmread('./topo');
                topomovie=dlmread('./topomovie');
                initialtopo=dlmread('./initialtopo');
                depthmovie=dlmread('./depthmovie');
                velmovie=dlmread('./velocitymovie');
                cmovie=dlmread('./cmovie');
                %Emovie=dlmread('./Emovie');
                Dmovie=dlmread('./Dmovie');
                Mmovie=dlmread('./Mmovie');
                uh=dlmread('./uh');
                ch=dlmread('./ch');
                vel=dlmread('./vel');
                c=dlmread('./c');
                c_1=dlmread('./c1');
                c_2=dlmread('./c2');
                c_3=dlmread('./c3');
                m_1=dlmread('./m1');
                m_2=dlmread('./m2');
                m_3=dlmread('./m3');
                m_4=dlmread('./m4');
                m_5=dlmread('./m5');
                m_6=dlmread('./m6');
                m_7=dlmread('./m7');
                m_8=dlmread('./m8');
                m_9=dlmread('./m9');
                m_10=dlmread('./m10');
                m=dlmread('./m');
                vh=dlmread('./vh');
                depth=dlmread('depth');
                E1=dlmread('./E1');
                E2=dlmread('./E2');
                E3=dlmread('./E3');
                E4=dlmread('./E4');
                D=dlmread('./D');
                solid=dlmread('solid');
                stage=dlmread('stage');
                %stage_2=dlmread('stage2');
                infl=dlmread('infl');
                maxvel=dlmread('maxvel');
                saveflow=dlmread('saveflow');
                
                topo=topo';
                initialtopo=initialtopo';
                topomovie=topomovie';
                depthmovie=depthmovie';
                velmovie=velmovie';
                cmovie=cmovie';
                %Emovie=Emovie';
                Dmovie=Dmovie';
                Mmovie=Mmovie';
                uh=uh';
                ch=ch';
                vel=vel';
                c=c';
                c_1=c_1';
                c_2=c_2';
                c_3=c_3';
                m_1=m_1';
                m_2=m_2';
                m_3=m_3';
                m_4=m_4';
                m_5=m_5';
                m_6=m_6';
                m_7=m_7';
                m_8=m_8';
                m_9=m_9';
                m_10=m_10';
                vh=vh';
                m=m';
                E1=E1';
                E2=E2';
                E3=E3';
                E4=E4';
                D=D';
                depth=depth';
                solid=solid';
                infl=infl';
                maxvel=maxvel';
                
                eval(['topo',num2str(idum),'=topo;'])
                eval(['topomovie',num2str(idum),'=topomovie;'])
                eval(['depthmovie',num2str(idum),'=depthmovie;'])
                eval(['velmovie',num2str(idum),'=velmovie;'])
                eval(['cmovie',num2str(idum),'=cmovie;'])
                %eval(['Emovie',num2str(idum),'=Emovie;'])
                eval(['Dmovie',num2str(idum),'=Dmovie;'])
                eval(['Mmovie',num2str(idum),'=Mmovie;'])
                eval(['depth',num2str(idum),'=depth;'])
                eval(['uh',num2str(idum),'=uh;'])
                eval(['vel',num2str(idum),'=vel;'])
                eval(['ch',num2str(idum),'=ch;'])
                eval(['vh',num2str(idum),'=vh;'])
                eval(['c',num2str(idum),'=c;'])
                eval(['c_1',num2str(idum),'=c_1;'])
                eval(['c_2',num2str(idum),'=c_2;'])
                eval(['c_3',num2str(idum),'=c_3;'])
                eval(['m_1',num2str(idum),'=m_1;'])
                eval(['m_2',num2str(idum),'=m_2;'])
                eval(['m_3',num2str(idum),'=m_3;'])
                eval(['m',num2str(idum),'=m;'])
                eval(['E1',num2str(idum),'=E1;'])
                eval(['E2',num2str(idum),'=E2;'])
                eval(['E3',num2str(idum),'=E3;'])
                eval(['E4',num2str(idum),'=E4;'])
                eval(['D',num2str(idum),'=D;'])
                eval(['solid',num2str(idum),'=solid;'])
                eval(['infl',num2str(idum),'=infl;'])
                eval(['stage',num2str(idum),'=stage;'])
                %eval(['stage_2_',num2str(idum),'=stage_2;'])
                eval(['saveflow',num2str(idum),'=saveflow;'])
                eval(['maxvel',num2str(idum),'=maxvel;'])
                
                
                eval(['dlmwrite(''topo_',num2str(idum),''',topo',num2str(idum),')'])
                eval(['dlmwrite(''depth_',num2str(idum),''',depth',num2str(idum),')'])
                eval(['dlmwrite(''uh_',num2str(idum),''',uh',num2str(idum),')'])
                eval(['dlmwrite(''vh_',num2str(idum),''',vh',num2str(idum),')'])
                eval(['dlmwrite(''c_',num2str(idum),''',c',num2str(idum),')'])
                eval(['dlmwrite(''c_1_',num2str(idum),''',c_1',num2str(idum),')'])
                eval(['dlmwrite(''c_2_',num2str(idum),''',c_2',num2str(idum),')'])
                eval(['dlmwrite(''c_3_',num2str(idum),''',c_3',num2str(idum),')'])
                eval(['dlmwrite(''m_',num2str(idum),''',m',num2str(idum),')'])
                eval(['dlmwrite(''vel_',num2str(idum),''',vel',num2str(idum),')'])
                eval(['dlmwrite(''maxvel_',num2str(idum),''',maxvel',num2str(idum),')'])
                eval(['dlmwrite(''stage_',num2str(idum),''',stage',num2str(idum),')'])
                eval(['dlmwrite(''topomovie_',num2str(idum),''',topomovie',num2str(idum),')'])
                eval(['dlmwrite(''velmovie_',num2str(idum),''',velmovie',num2str(idum),')'])
                eval(['dlmwrite(''cmovie_',num2str(idum),''',cmovie',num2str(idum),')'])
                eval(['dlmwrite(''depthmovie_',num2str(idum),''',depthmovie',num2str(idum),')'])
                
                %saveinputvec(idum,:)=[nx ny dx F(bb) tend a(hr1) ad(hr2) h0 mtstar0 rnum rint length(dsin) normalmeanks normalvarks max(max(hf)) hillsloperoughness pi Si Ki gi vegcover];
                eval(['saveinputvec',num2str(idum),'=saveinputvec;'])
                
                savea(idum)=a(hr1);
                savead(idum)=ad(hr2);
                saveF(idum)=F(bb);
                
                % Plot A Few Results
                erosion=topo-initialtopo;
                erosion(erosion==0)=nan;
                figure(100+idum)
                set(gca,'FontName','Arial','FontSize',18)
                surf(xgrid,ygrid,erosion,erosion)
                shading interp
                axis tight
                xlabel('Distance (m)')
                ylabel('Distance (m)')
                title('Total Erosion (m)')
                view([90 90])
                caxis([-0.02 0.02])
                grid off
                box on
                colorbar
                set(gca,'DataAspectRatio',[1 1 1])
                set(gca,'FontName','Arial','FontSize',18)
                eval(['savefig(''TotalErosion',num2str(idum),''')'])
                
            end
            
        end
    end
end

erosion=topo-itopo;
erosion(erosion==0)=nan;
vel(maxvel==0)=nan;
depth(maxvel==0)=nan;
uh(maxvel==0)=nan;
vh(maxvel==0)=nan;
c(maxvel==0)=nan;
c_1(maxvel==0)=nan;
c_2(maxvel==0)=nan;
c_3(maxvel==0)=nan;
maxvel(maxvel==0)=nan;

% Grab A Slice of Data From the Movie File
j=49;
topotime1=topomovie(:,ny*j+1:ny*(j+1));

