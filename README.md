# NEOLAB

---

*Authors*: Mikaël Cugnet, Florian Gallois, Angel Kirchev, Denys Dutykh

*Last page update*: April 02, 2021

*Latest version*: 6.1.0

*Licence*: GNU Lesser General Public License (2.1)

---

Based on the existing knowledge in physics, chemistry and mathematics, 
NEOLAB is a lightweight library which provides the minimum working example to simulate the behaviour of the negative electrode of a lead-acid battery.
The model is etablished according to a minimal set of ordinary and partial differential equations describing the physics 
behind the lead electrode used as the negative electrode. A study to optimize the performance of the software has been carried out, in particular the impact 
of equations scaling. 
The toolbox consists of two Scilab routines, the main one in which are set up all the varibles, and a second one computating the residue of the equations.
The resolution of equations is done with a Scilab differential/algebraic system solver.

---

## Download and install
All the user needs to do is downloading the last version of Scilab (6.1.0) and open the main file. When compilling the program, 
the path to the residue function must be specified. The user has the choice to proceed the simulation with both dimensional and dimensionless methods. 

---

## Problems or questions
If you have any problems or questions, please contact the author: Mikaël Cugnet (mikael.cugnet@cea.fr)

---

## Files and folders
* README.md: Basic information (this file)
* CoA.txt: Authors' identity and contributions
* Licence.txt: Licence information
*data.sod: Experimental data
* main_dassl_dim.sce: Main routine using dassl solver and dimensional method 
* main_dassl_adim.sce: Main routine using dassl solver and dimensionless method
* main_dassl_dim_tol.sce: Main routine using dassl solver and dimensional method, with a tolerence study
* main_dassl_adim_tol.sce: Main routine using dassl solver and dimensionless method, with a tolerence study
* res_dassl_dim.sci: Residue function using dassl solver and dimensinal method
* res_dassl_adim.sci: Residue function using dassl solver and dimensinless method
* res_dassl_dim_tol.sci: Residue function using dassl solver and dimensinal method, with a tolerence study
* res_dassl_adim_tol.sci: Residue function using dassl solver and dimensinless method, with a tolerence study