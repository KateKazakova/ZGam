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

###Functions included in the package
1. RfactorMC() - allows to calculate the correlation factor on MC
2. RfactorDataCounting() - allows to calculate the correlation factor with Data-Driven method
3. Leakage() - allows to define the leakage parameters
4. CentralValueCounter() - allows to find the number of background events in signal region induced by jet->γ misidentification
