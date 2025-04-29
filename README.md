# BasiliskRiverModel
To install, first download and install [basilisk C](http://basilisk.fr/src/INSTALL). 
Then clone this repository. You can compile by writing 
```
make 
```
in your terminal, then you can run the simulation by writing
```bash
./river n
```
with ```n``` being the resolution integer. Default is seven, and incrementing/decrementing the value with one will double/halve the amount of grid points. To try out different parameters, edit the parameter list in ```river.c```, then recompile by writing ```make``` over again in the terminal. 

The simulation will output all data to the ```plots```-directory, which is made during compilation. The outputed data comes in the form of ```.vts```-files which can be read with [paraview](https://www.paraview.org/). If you have any questions, open an issue here on github, or send me an email if you are taking the course MEK3800/4800 at UiO. 