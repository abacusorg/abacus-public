## Tests

## Ewald
The Ewald test checks the accuracy of the forces against a pre-computed result. See [`Ewald/README.md`](Ewald/README.md).

<details>
<summary>Example output</summary>

```console
nid001044:~/abacus/Tests/Ewald$ ./run_ewald.py 
Beginning run at 2024-05-23 08:05:58, running for at most 230.0 minutes.

Using parameter file abacus.par and working directory /pscratch/sd/l/lgarriso/Ewald.
Copying derivatives from "/pscratch/sd/l/lgarriso/Derivatives" to "/dev/shm/lgarriso/Derivatives"
Ingesting IC from /pscratch/sd/l/lgarriso//Ewald/ic ... Skipping convolution
Beginning abacus steps:
Running singlestep for step 0
        Finished state 0. a = 1.0000, dlna = 0, rms velocity = 0
Performing convolution for step 1
         Finished convolution for state 1
Running singlestep for step 1
        Finished state 1. a = 1.0000, dlna = 0, rms velocity = 0
Final redshift of 0 reached; terminating normally after 0.01 hours.
Reached maxsteps limit at Thu May 23 08:06:19 2024; exiting job w/o requeue after 0.01 hours.

    (New) Ewald Test
    ================
    
    Fractional error
    ----------------
              Max: 2.2759e-03
              99%: 1.3711e-04
           Median: 1.0078e-05
              Min: 5.8454e-08
    Min (nonzero): 5.8454e-08
    

    Configuration
    -------------
      CPD: 11
    Order: 8
      ppc: 49.2
    dtype: float32
    
Plotting to ewald_storeforces.png ...
Plotting to ewald_storeforces.large.png ...
```
</details>


## CheckForces
CheckForces checks that the forces on a simple cubic lattice of particles are held to zero, or nearly zero.  The results are quoted as equivalent displacements, i.e. the displacements that would yield the measured force under the Zel'dovich approximation, in units of the lattice spacing.

This test is particularly useful for big MPI parallel runs because the IC lattice is generated on the fly. The force results are written to disk and can be analyzed with the provided helper script.

Powers of 2 are often useful lattice edge sizes, as they will be co-prime with CPD and thus maximally sample the unit cell.

### Usage
```console
$ cd CheckForces
$ ./check_forces.py [ppd] [cpd]
```

<details>
<summary>Example output</summary>

```console
nid001044:~/abacus/Tests/CheckForces$ ./check_forces.py 256 65
Running CheckForces with PPD 256, CPD 65
Reusing existing ICs
Beginning run at 2024-05-23 08:26:22, running for at most 230.0 minutes.

Using parameter file abacus.par and working directory /pscratch/sd/l/lgarriso/CheckForces.
Beginning abacus steps:
Running singlestep for step 0
        Finished state 0. a = 0.0100, dlna = 0, rms velocity = 0
Performing convolution for step 1
         Finished convolution for state 1
Running singlestep for step 1
        Finished state 1. a = 0.0100, dlna = 0, rms velocity = 0
Final redshift of 99 reached; terminating normally after 0.01 hours.
Reached maxsteps limit at Thu May 23 08:26:43 2024; exiting job w/o requeue after 0.01 hours.
Read acc time: 0.10 seconds (2.78e+03 MB/s)
Compute stats time: 0.91 seconds (18.4 Mpart/s)
Force statistics (equivalent ZA displacement in units of interparticle spacing):
Max: 2.6115e-05
Min: 9.90571e-10
Min (nonzero): 9.90571e-10
Mean: 1.42297e-06
RMS: 2.33235e-06
Plotting to checkforces_storeforces_absolute.png
```
</details>

## Spiral
Spiral compares the growth of a plane-wave perturbation against the analytic result. Particles actually move in this test, so it exercises the time integration as well as the forces.

The default Spiral is 256^3. Bigger versions can be run for exercising the MPI code, but it's usually not as convenient as CheckForces because it involves writing ICs to disk.

### Usage
```console
$ cd Spiral
$ ./SpiralTest
```

<details>
<summary>Example output</summary>

```console
./nid001044:~/abacus/Tests/Spiral$ ./SpiralTest.py 
   N1D: 256
     N: 16777216
  kvec:  1.000000e+00   0.000000e+00   0.000000e+00
 phase:  3.141593e+00   0.000000e+00   0.000000e+00
 Ainit:  9.000000e-02
Across:  1.000000e-01
fsmooth:  3.000000e-01
Beginning run at 2024-05-23 08:12:45, running for at most 230.0 minutes.

Using parameter file abacus.par and working directory /pscratch/sd/l/lgarriso/spiral.
Ingesting IC from /pscratch/sd/l/lgarriso//spiral/ic ... Skipping convolution
Beginning abacus steps:
Running singlestep for step 0
        Finished state 0. a = 0.0900, dlna = 0, rms velocity = 0
Performing convolution for step 1
         Finished convolution for state 1
Running singlestep for step 1
        Finished state 1. a = 0.0922, dlna = 0.025, rms velocity = 0.0245
Performing convolution for step 2
         Finished convolution for state 2
Running singlestep for step 2
        Finished state 2. a = 0.0932, dlna = 0.0108, rms velocity = 0.0253
Performing convolution for step 3
         Finished convolution for state 3
Running singlestep for step 3
        Finished state 3. a = 0.0943, dlna = 0.0108, rms velocity = 0.0256
Performing convolution for step 4
         Finished convolution for state 4
Running singlestep for step 4
        Finished state 4. a = 0.0953, dlna = 0.0108, rms velocity = 0.026
Performing convolution for step 5
         Finished convolution for state 5
Running singlestep for step 5
        Finished state 5. a = 0.0963, dlna = 0.0108, rms velocity = 0.0263
Performing convolution for step 6
         Finished convolution for state 6
Running singlestep for step 6
        Finished state 6. a = 0.0974, dlna = 0.0108, rms velocity = 0.0267
Performing convolution for step 7
         Finished convolution for state 7
Running singlestep for step 7
        Finished state 7. a = 0.0984, dlna = 0.0108, rms velocity = 0.0271
Performing convolution for step 8
         Finished convolution for state 8
Running singlestep for step 8
        Finished state 8. a = 0.0995, dlna = 0.0109, rms velocity = 0.0275
Performing convolution for step 9
         Finished convolution for state 9
Running singlestep for step 9
        Finished state 9. a = 0.1006, dlna = 0.0109, rms velocity = 0.0279
Performing convolution for step 10
         Finished convolution for state 10
Running singlestep for step 10
        Finished state 10. a = 0.1017, dlna = 0.0109, rms velocity = 0.0283
Performing convolution for step 11
         Finished convolution for state 11
Running singlestep for step 11
        Finished state 11. a = 0.1028, dlna = 0.0109, rms velocity = 0.0287
Performing convolution for step 12
         Finished convolution for state 12
Running singlestep for step 12
        Finished state 12. a = 0.1038, dlna = 0.0101, rms velocity = 0.0291
Performing convolution for step 13
         Finished convolution for state 13
Running singlestep for step 13
        Finished state 13. a = 0.1045, dlna = 0.00684, rms velocity = 0.0295
Performing convolution for step 14
         Finished convolution for state 14
Running singlestep for step 14
        Finished state 14. a = 0.1051, dlna = 0.006, rms velocity = 0.0298
Performing convolution for step 15
         Finished convolution for state 15
Running singlestep for step 15
        Finished state 15. a = 0.1058, dlna = 0.00627, rms velocity = 0.03
Performing convolution for step 16
         Finished convolution for state 16
Running singlestep for step 16
        Finished state 16. a = 0.1065, dlna = 0.00677, rms velocity = 0.0302
Performing convolution for step 17
         Finished convolution for state 17
Running singlestep for step 17
        Finished state 17. a = 0.1073, dlna = 0.00718, rms velocity = 0.0305
Performing convolution for step 18
         Finished convolution for state 18
Running singlestep for step 18
        Finished state 18. a = 0.1080, dlna = 0.00711, rms velocity = 0.0308
Performing convolution for step 19
         Finished convolution for state 19
Running singlestep for step 19
        Finished state 19. a = 0.1089, dlna = 0.00772, rms velocity = 0.031
Performing convolution for step 20
         Finished convolution for state 20
Running singlestep for step 20
        Finished state 20. a = 0.1097, dlna = 0.00766, rms velocity = 0.0313
Performing convolution for step 21
         Finished convolution for state 21
Running singlestep for step 21
        Finished state 21. a = 0.1100, dlna = 0.00277, rms velocity = 0.0316
Performing convolution for step 22
         Finished convolution for state 22
Running singlestep for step 22
        Finished state 22. a = 0.1100, dlna = 0, rms velocity = 0.0317
Final redshift of 8.09091 reached; terminating normally after 0.13 hours.
(16777216, 6)
(8192, 2)
Vx: rms 0.095616, max 0.135475
Vy: rms 3.862041e-09, max 1.438581e-07
Vz: rms 3.872732e-09, max 1.458242e-07
Ratio of max velocity (analytic/computed): 0.999320
All particles present with no duplicates.
```
</details>

## Galilean
Galilean advects a lattice with constant velocity, checking that we haven't lost any particles. It's particularly useful for exercising the particle exchange in the MPI code. The test itself is somewhat dependent on the internal state format.
