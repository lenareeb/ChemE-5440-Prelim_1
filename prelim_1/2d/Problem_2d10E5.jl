using DifferentialEquations
using Plots
gr(size=(500,500), show = true)  #use the gr backend for plotting

#Function for the lorenz equation
#u[1] = x; u[2] = y; u[3] = z
function lorenz!(du,u,p,t)
 ax = 0.15
 ay= 0.027
 Bx = 5.9
 By= 5.4
 dy= 1.07
 dz= 1.12
 zx = 0.000064
 yz= 0.011
 xz= 0.12
 xy= 0.00083
 nzx = 2.34
 nxz = 2.0
 nxy= 2.0
 nyz= 2.0
 S=100000
 du[1] = -u[1]+(ax+Bx*S)/(1+S+(u[3]/zx)^nzx)      #dx/dt
 du[2] = -dy*u[2]+(ay+By*S)/(1+S+(u[1]/xy)^nxy)
 du[3] = -dz*u[3]+1/(1+(u[2]/yz)^nyz+(u[1]/xz)^nxz)                   #dz/dt
end

u0 = [0.0;0.0;0.0]                      #intial conditions
tspan = (0.0,100.0)                     #start and end time
prob = ODEProblem(lorenz!,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                       #Solve the system

#Plot the results; the vars=(0,1) argument specifies to plot X (column 1 of sol)
#vs t (column 0 of sol)
plt1 = plot(sol,vars=(0,1), xaxis="t", yaxis = "X" )
display(plt1)