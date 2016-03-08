% Overland Flow Model
% written by Emily Fairfax
% March 5th, 2016

clear all

%% Initialize
% Constants
    % Topography Constants
    topomax = 100; % max elevation of valley in meters
    m = .05; % slope of valley (for linear topography)
    
    % Flow Constants
    n = 0.05; % semi-weedy water
    
    % Rainfall Constants
    pwater = 100; %the period of the rainfall oscillations
    minutes = 15; %duration of the storm in minutes
    stormlength = 60*minutes; %length of storm in seconds
    stormrainheight = 0.02/stormlength; %meters of rain fallen over the length of the storm
    stormvelocity = .5; %m/s, velocity storm moves down slope
    
   
% Space Array
    dx = .5; % horizontal step in meters
    xmin = 0; % min horizontal space in meters
    xmax = 6*topomax; % max horizontal space in meters
    x = xmin:dx:xmax; % x array
    xedge = xmin+(dx/2):dx:xmax-(dx/2); % xedge array  
    
% Time Array
    dt = .1; % time step in seconds
    maxminutes = 90; % number of minutes to run the simulation
    tmax = 60*maxminutes; % tmax in seconds
    t = 0:dt:tmax; % time array
    imax = length(t); % loop length
    nplot = 2*maxminutes; % number of plots
    tplot = tmax/nplot; % times at which plots will be made
    
% Variable Arrays
    N = length(x); % nodes
    
    % Ground Arrays
    earthelevation = zeros(1,N); %create array for elevation profile of earth
    earthelevation(1:N) = topomax - m.*x; % linear valley profile
    slope = zeros(1,N); %preallocate slope array
    
    %Water Arrays
    Hinitial = 0; % initial water height
    H = Hinitial.*ones(1,N); %create array of water height in space based on initial height
    Hedge = zeros(1,N-1); %create array for height of water at edges of boxes
    jmax = length(Hedge); %loop for filling Hedge array
    for j=1:jmax
        Hedge(j) = (H(j)+H(j+1))/2; % calculate thickness of water at edges by averaging boxes 
    end 
    waterelevation = earthelevation(1:N) + H(1:N); %water elevation is elevation of earth plus water height
    
    %Flux Arrays
    Q = zeros(1,N); %preallocate Q array
    dQdx = zeros(1,N-2); %preallocate dQdx array
    BaseHydro = zeros(1,imax);
    MiddleHydro = zeros(1,imax);
    TopHydro = zeros(1,imax);
    
    %Storm Arrays
    stormfront = zeros(1,imax);
    stormtail = zeros(1,imax);
   
    %For Plotting
    xx = [x, fliplr(x)];
    bottomline=-200000000.*ones(1,N); %%Create a bottom line array for use in filling plots, needs to be same length as N
    
%% Run
for i = 2:imax
    % Slope Determination
    slope = abs(diff(waterelevation./dx));
    
    % Flux from Flow
    Q(2:N) = 1/n .* Hedge.^(5/3).*slope.^(1/2); % manning equation for mean water velocity depth averaged
    Q(1) = 0; % no water flux into the top, boundary condition
    
    %For tracking hydrographs, these are the locations of the monitoring
    %stations
    BaseHydro(i) = Q(N);
    MiddleHydro(i) = Q(ceil(N/2));
    TopHydro(i) = Q(ceil(N/10));
    
    %Net Fluxes
    dQdx(1:N-1) = diff(Q)./dx; % diff Q to get dQdx
    dQdx(N) = dQdx(N-1); % net flux out of last box is same as box before it, allows model to drain...boundary condition

    
    % Rainfall and Infiltration
    
    %Constant Rainfall Along Slope
    R = 0.000003.*ones(1,N); 

    %Infiltration
    I(1:N) = 0*ones(1,N); % constant infiltration along slope

    
    % Water Height Changes
    dHdt = -dQdx + R - I; % conservation statement: change in water height is a function of fluxes, rainfall, and infiltration
    H = H + dHdt.*dt; % update total height of water
    Hneg = find(H<0); % find negative water heights
    H(Hneg)=0; % set negative water heights to zero... no such thing as negative water.
    
    % Apply a Noise Filter to Smooth Instability at Wave Front
    windowsize = 11;
    b = (1/windowsize)*ones(1,windowsize);
    a = 1;
    smoothH = filter(b,a,H);
    
    % Recalculate Hedge
    jmax = length(Hedge);
    for j=1:jmax
        Hedge(j) = (smoothH(j)+smoothH(j+1))/2; % calculate thickness of water at edges by averaging boxes 
    end 
    
    %Update Elevation Profile
    waterelevation = earthelevation + H; %water = earth + waterheight
    
    %Plot the Results Each Time Step   
    if (rem(t(i),tplot)==0)
        figure(1)
        clf
        
        %Water Thickness Plot
            subplot('position',[.1 .55 .8 .4])
            plot(x,smoothH*100, 'k')
            yy = [bottomline, fliplr(smoothH*100)];
            fill(xx,yy,[.373 .557 .627]);

            %plot formatting
            title('Water Thickness Through Time');
            %xlabel('Distance Along Profile (m)');
            ylabel('Water Thickness (cm)');
            set(gca,'fontsize',16,'fontname','arial')
            ht=text(4/5*xmax,1.2,['  ',num2str(floor(round(t(i))/60)), ' minutes  '],'fontsize',18); %print time in animation
            axis([0 xmax 0 1.5])
        
        %Flux Plot
            subplot('position',[.1 .32 .8 .15])
            plot(x,Q*10000,'Color',[.627 .373 .494],'linewidth', 4)

            %plot formatting
            title('Discharge Through Time');
            %xlabel('Distance Along Profile (m)');
            str = {'Discharge per','Unit Width (cm^2/s)'};
            ylabel(str);
            set(gca,'fontsize',14,'fontname','arial')
            axis([0 xmax 0 25])
        
        %Rainfall Pattern Plot
            subplot('position',[.1 .1 .8 .15])
            zz = [bottomline, fliplr(R*3600*100)];
            fill(xx,zz,[.118 .18 .192]);
            hold all
            plot(x,R*3600*100,'Color', [.224 .38 .38],'linewidth',2)
            
            %plot formatting
            title('Rainfall Pattern Through Time');
            xlabel('Distance Along Profile (m)');
            ylabel('Rainfall (cm/hr)');
            set(gca,'fontsize',14,'fontname','arial')
            axis([0 xmax 0 2])
        
        
        pause(0.02)
        hold off
    end
end


%% Finalize
%plot hydrograph with time at different points along slope

figure(2)
clf

subplot(3,1,3)
plot(t/60,BaseHydro*10000,'Color',[.627 .373 .494],'linewidth', 2)
    title('Hydrograph Monitoring Station at Base of Slope')
    xlabel('Time (minutes)');
    ylabel(str);
    set(gca,'fontsize',14,'fontname','arial')
    axis([0 tmax/60 0 25])
subplot(3,1,2)
plot(t/60,MiddleHydro*10000,'Color',[.627 .373 .494],'linewidth', 2)
    title('Hydrograph Monitoring Station at Middle of Slope')
    %xlabel('Time (minutes)');
    ylabel(str);
    set(gca,'fontsize',14,'fontname','arial')
    axis([0 tmax/60 0 25])
subplot(3,1,1)
plot(t/60,TopHydro*10000,'Color',[.627 .373 .494],'linewidth', 2)
    title('Hydrograph Monitoring Station at Top of Slope')
    %xlabel('Time (minutes)');
    ylabel(str);
    set(gca,'fontsize',14,'fontname','arial')
    axis([0 tmax/60 0 25])


