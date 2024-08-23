include("sonarForms.jl");

using DifferentialEquations;

kr = 1.0/CoppensFormulaBase(LayerData(2.0,34.7,0.0));

f(z,p,s) = kr / (sign(kr) * (sqrt(abs((1.0/CoppensFormula(LayerData(2.0,34.7,z)))^2 - kr^2)))) 

z0 = 500.0
sSpan = (0.0,6000.0)
prob = ODEProblem(f,z0,sSpan);

sol = solve(prob,Euler(),reltol = 0.001,abstol = 0.001,tstops = 1.0);

using GLMakie;

print(sol)
plot(sol)
