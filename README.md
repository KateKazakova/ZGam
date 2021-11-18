# The Zgamma package
 Background estimation induced by jet->γ misidentification from electroweak Z(νν̄)γ process

 ## Example of checking out (CMake)

 ```
 mkdir MyAnalysis;
 cd MyAnalysis;
 mkdir source run build;
 cd source;
 git clone;
 ```

 For compiling the package:
 ```
 cd build;
 cmake ../source;
 cmake --build .;
```

### Functions included in the package
* RfactorMC() - allows to calculate the correlation factor on MC
* RfactorDataCounting() - allows to calculate the correlation factor with Data-Driven method
* Leakage() - allows to define the leakage parameters
* CentralValueCounter() - allows to find the number of background events in signal region induced by jet->γ misidentification
