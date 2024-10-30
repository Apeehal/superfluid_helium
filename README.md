# superfluid_helium
This is a 2D simulation of superfluid helium flowing down a pipe, which handles both the fluid flow and heat transfer via the FDM method.
The main source for the equations that are being solved can be found here:

https://cds.cern.ch/record/1479179/files/CERN-THESIS-2012-120.pdf

![image](https://github.com/user-attachments/assets/4fdba8fa-bdf1-47c9-b3db-5567d48c8b6b)

ρ_s = superfluid density
ρ_n = normal density
T = Temperature
ρ = superfluid density + normal density
P = pressure
s = entropy

![image](https://github.com/user-attachments/assets/13842639-e58e-4309-ba58-6ebc1eaede54)

q = heat input
d = diameter of pipe
β = 8 (constant from Hagen–Poiseuille equation)
η = viscosity

Additionally some data was taken to linearly interpolate the normal density and entropy values from the tables from this source:
https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1029.pdf


For now the simulation is only in the Landau regime, but in the future it would be interesting to add turbulence and Kapitza resistance.

Note that all units are in SI.
