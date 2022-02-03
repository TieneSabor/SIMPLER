# CFD Final Project
## Program for SIMPLER Algorithm
### Revisit the Theory
#### Momentum Equation
- We take momentum equation for velocity in x direction(u) for example
$\frac{\partial(\rho u)}{\partial t}+\frac{\partial (\rho uu)}{\partial x}+\frac{\partial (\rho vu)}{\partial y}=-\frac{\partial P}{\partial x}+\mu(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 v}{\partial y^2})$
- Then define
$J_x=\rho uu-\mu\frac{\partial u}{\partial x}$ ; $J_y=\rho vu-\mu\frac{\partial v}{\partial y}$
$\implies\frac{\partial(\rho u)}{\partial t}+\frac{\partial (J_x)}{\partial x}+\frac{\partial (J_y)}{\partial y}=-\frac{\partial P}{\partial x}$
- We adopt following discretization
    - $\frac{\partial(\rho u)}{\partial t}\approx\frac{\rho u-\rho_0 u_0}{\Delta t}$
    - $\frac{\partial (J_x)}{\partial x}\approx\frac{J_{x,e}-J_{x,w}}{\Delta x}$
    - $\frac{\partial (J_y)}{\partial y}\approx\frac{J_{y,n}-J_{y,s}}{\Delta y}$
    - $\frac{\partial P}{\partial x}\approx\frac{P_e-P_w}{\Delta x}$
- From "Numerical Heat Transfer and Fluid Flow, Patankar", we can se from Equation (5.56)~(5.60), the discretize momentum equation is as follow:
    - $a_Pu_P=a_Eu_E+a_Wu_W+a_Nu_N+a_Su_S+b+Area(P_w-P_e)+WallFriction$
- We also adopted "Power Law" scheme for "Convection-Deffusion" equation and thus $A(|P|)=MAX(0,(1-0.1|P|)^5)$
#### Pressure Equation
- Consider the discretized continuity equation: $\frac{\rho_P-\rho_P^0\Delta x\Delta y}{\Delta t}+(\rho u_e-\rho u_w)\Delta y+(\rho u_n-\rho u_s)\Delta x=0$
- Combine the continuity equation and velocity correction equations, we get $a_PP_P=a_EP_E+a_WP_W+a_NP_N+a_SP_S+b$
### Programming
- Following is the procedure of SIMPLER
    - Calculate pseudo velocity: $\hat{u}=\frac{\sum a_{nb}u_{nb}+b}{a}$
    - Calculate the source terms of Pressure Equation with pseudo velocity 
    - Solve the Pressure Equation for the pressure field $P^*$
    - Calculate the source term of Momentum Equation with the pressure field $P^*$
    - Solve the Momentum Equation for the velocity field $u^*$
    - Calculate the source terms of Pressure Correction Equation with velocity $u^*$
    - Solve the Pressure Correction $P'$
    - Correct the velocity using $u^{n+1}=u^*+d(P_{w}'-P_{e}')$
    - Iterate above procedures.
- Here are the instructions to the project codes, the codes can be run in ubuntu 20.04. Also, to plot data correctly, we need python3, matplotlib and numpy packages of python.
```bash=
# get into the code folders
cd codes
# compile the program
make clean
make all
# To see Couette Flow, please follow the commands:
./bin/Final a > data1.txt
python3 script/show_flowfield.py -p data1.txt
# Then the flow field streamline and pressure contour will show
# Similarly, to see Cavity Flow, please type
./bin/Final b > data2.txt
python3 script/show_flowfield.py -p data2.txt
# Then the flow field streamline and pressure contour will show
```
## Solve Problems Using the Program
### Couette Flow
#### Analytical Solution
- From fundamental fluid mechanics, we knew that the solution for Couette Flow without external pressure field is simple:
$U(y)=U_b+(U_t-U_b)y/H$
    - $U(y)$: velocity in x variation along y axis
    - $U_b$: velocity at bottom plate
    - $U_t$: velocity at top plate
    - $H$: height between 2 plates
#### Boundary Conditions:
![](https://i.imgur.com/evNdqgR.png)
**(Fig 1.) Boundary Condition of the Couette Flow Case**
- At wall boundaries, we assume the flow is laminar and thus wall friction $\tau=Area*\mu \frac{\partial u}{\partial y}$
#### Flow field parameters:
- Computational domain: Height=0.025m; Width=1m
- $u_0=0.01\ (m/s)$; $Re=250$
- $\rho=998\ (kg/m^3)$
- $\mu=0.001\ (Pa\cdot s)$
#### Computational Parameters
- $\Delta t=0.1$
- Number of grids in y-direction: 10
- Number of grids in x-direction: 10
#### Results
- Stream line plot
![](https://i.imgur.com/CKlKuwx.png)
**(Fig 2.) Stream line of the Couette Flow Field**
- Pressure Contour
![](https://i.imgur.com/UXr6BMc.png)
**(Fig 3.) Pressure Contour of the Couette Flow Field**
- Variation in y axis of x-Velocity 
![](https://i.imgur.com/Ao08zH2.png)
**(Fig 4.) Comparison between analytical and numerical solution**

### Lid Driven Cavity FLow
#### Boundary Conditions:
![](https://i.imgur.com/clIzGMU.png)
**(Fig 5.) Boundary Condition of the Cavity Flow Case**
- At wall boundaries, we assume the flow is laminar and thus wall friction $\tau=Area*\mu \frac{\partial u}{\partial y}$
#### Flow field parameters:
- Computational domain: Height=0.01m; Width=0.01m
- $u_0=0.1\ (m/s)$; $Re=70$
- $\rho=1.23\ (kg/m^3)$
- $\mu=1.78*10^{-5}\ (Pa\cdot s)$
#### Computational Parameters
- $\Delta t=0.01$
- Number of grids in y-direction: 10
- Number of grids in x-direction: 10
#### Results
- Stream line plot
![](https://i.imgur.com/xHcgwjC.png)
**(Fig 6.) Stream line of the Cavity Flow Field**
- Pressure Contour
![](https://i.imgur.com/rpXnIMT.png)
**(Fig 7.) Pressure Contour of the Cavity Flow Field**

#### Validation
- I test 3 kind of grid size: $10*10,\ 15*15,\ 20*20$ to see if the value converges. Below are the variation of x-velocity on y axis at x=0.05m. 
![](https://i.imgur.com/Os2l6vU.png)
**(Fig 8.) Use the x-velocity along y axis to see if the result converges with finer grids**

### Conclusion
1. We can see there are 2 vortexes inside the cavity flow (Fig. 6) when Reynolds number is 70. One is larger and one is smaller.
2. The calculation result in Couette flow case align well with the analytical solution. (Fig. 4)
3. The validation test in Cavity flow case showed that the result converges to a certain curve as the grids become finer; that is, more grid points in same computational domain. (Fig. 8)
4. It shows that with given boundary condition and parameters, SIMPLER algorithm gives reasonable results in both cases.

### Reference
[1] "Numerical Heat Transfer and Fluid Flow", Suhas Patankar, 1980
[2] "Computational Fluid Dynamics", John Anderson, 1995
[3] "An Introduction to Computational Fluid Dynamics", Henk K. Versteeg, 2007
