%ALL CSV FILES PULLED FROM https://github.com/CSSEGISandData/COVID-19/
%specifically daily reports file. Collected as of 4/25/20
%attack rate source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7104686/
%efficacy equation: https://en.wikipedia.org/wiki/Vaccine_efficacy

location = mfilename('fullpath');
location = erase(location, 'project_ee_608');
location = append(location, '\*.csv');
directory = dir(location);
directory = {directory.name}; %import all csv file names to matlab for use
confirmed = []; % array of confirmed casses per date

%collect confirmed cases per day for all files
for i = 1:1:length(directory)
    x = readtable(char(directory(i))); %read a csv file by name
    data = table(x.Confirmed); 
    data{:,:}(isnan(data{:,:})) = 0;
    confirmed = [confirmed sum(data.Var1)]; %append confirmed cases to array
end

%variables
pop = 7594000000; %current population
enddate = 4000; %final day
vaccinated = 10000000; %amount vaccinated per day
startdate = 1000; %start date for vaccination introduction
ARU = .565; %infectivity of the disease, based of average of male and female AR
ARV = .3; %theroretical infectivity of vaccinated individuals

days = 1:1:enddate;

fit = polyfit(1:1:length(confirmed), confirmed, 2); %order is 2 because best fit
infecteq = polyval(fit,days); %view infection for the next 500 days

fprintf ('Infection Equation:\n%fx^2 + %fx + %f \n', fit(1), fit(2), fit(3))

%negative infection not possible
for i = 1:1:length(infecteq)
    if (infecteq(i) <= 0)
       infecteq(i) = 0;
    end
    if (infecteq >= pop)
       infecteq(i) = pop;
    end
end

%vaccine efficiancy
efficacy = (ARU-ARV)/ARU;

vacceq = [];
vaccdays = startdate:1:enddate;

for i = 1:startdate
    vacceq = [vacceq, pop];
end
vacceq = [vacceq, pop - (efficacy*vaccinated*(vaccdays-startdate))];

for i = 1:1:length(vacceq)
    if (vacceq(i) <= 0)
       vacceq(i) = 0;
    end
end

figure
plot (infecteq)
hold on
plot (vacceq)

stop = false;


for i = 1:length(infecteq)
    if (vacceq(i) - infecteq(i) < 0 && stop == false)
        dates = [i-1 i];
        stop = true;
    end
end

fprintf ('Intercept between days %d and %d \n', dates(1), dates(2))

%find exact number of people infected
syms x

%cap to determine what function we use
test = fit(1)*x.^2 + fit(2)*x + fit(3);
startmax = double(solve(test == pop));
for i = 1:length(startmax)
    if (startmax(i) > 0)
        startmax = startmax(i);
    end
end

if (startdate < startmax) %vaccine developed before hardcap
    funct = fit(1)*x.^2 + fit(2)*x + fit(3) - (pop - (efficacy*vaccinated*(x - startdate))) == 0;
else %vaccine developed post hardcap
   funct =  fit(1)*x.^2 + fit(2)*x + fit(3) - pop;
end

%calculat population infected at max
amount = double(solve(funct));
for i = 1:length(amount)
    if (amount(i) > 0)
        amount = amount(i);
    end
end
number_of_people = fit(1)*amount^2 + fit(2)*amount + fit(3);

fprintf('The total amount of people infected is %f people \n', number_of_people)