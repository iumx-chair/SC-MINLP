# SC-MINLP

## Dependencies for Python

### CPLEX

SC-MINLP requires the full version of the CPLEX solver if you want to run models without restrictions.
If you already have a version installed on your system, you need to link it with Python with the command below:

```console
python3 <CPLEX-DIRECTORY>/python/setup.py install
```

###### Some default values for  <CPLEX-DIRECTORY>
| OS            | Path                                              |
| ------        | -----------                                       |
| Linux         | /opt/ibm/ILOG/CPLEX_Studio[edition]\<version>     |
| MacOS         | /Applications/CPLEX_Studio[edition]\<version>     |
  
  
### SymPy
 
SymPy is a library for symbolic mathematics. Use the command below to install this package for Python :
  
```console
python3 -m pip install sympy
```

## Other dependencies
  
To install and compile the ASL library, you can use the manual process or use the script below:
  
```console
cd aslpy
./GET_asl.sh
```
  
Then, you must compile the shared library responsible for communication between SC-MINLP and ASL:

```console
cd aslpy
./build.sh
```
  
### Test

To test the SC-MINLP, type the command:
```console
 ./SCMINLP
 ```

The solver will load the file "tmp/prob.nl".
  
## Contributors
  * [Renan Spencer Trindade](https://www.renan-st.com/) (École polytechnique, France)
  * [Claudia D'Ambrosio](https://www.lix.polytechnique.fr/~dambrosio/) (École polytechnique, France)
  * [Antonio Frangioni](http://www.di.unipi.it/~frangio/) (Università di Pisa, Italy)
  * [Claudio Gentile](http://www.iasi.cnr.it/~gentile/) (Istituto di Analisi dei Sistemi ed Informatica "Antonio Ruberti", Italy)
