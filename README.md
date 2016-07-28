# Diffusion-analysis---Moving-boundary-condition
Following are the details of the machine I used to compile and run the code:
MacBook Air (13-inch, Early 2015)
Processor : 1.6 GHz Intel Core i5
Memory : 4 GB 1600 MHz DDR3



Example Input: 

In the file "ConservativeIFF.cpp"

In line 17/153

int lambda = 1;                         // Parameter describing geometry of system (=1 for p$

        double s_0 = 1.;                        // Initial interface position
        double R   = 5.;                        // Position of far boundary

        int nAlpha = 100;                       // Number of points in phase alpha
        double dAlpha = 1.22664E-11;            // Diffusion coefficient of phase alpha
        double initialAlpha = 0.63;             // Initial concentration in phase alpha
        double interAlpha   = 0.4545;           // Interfacial concentration in phase alpha


        int nBeta    = 100;
        double dBeta = 8.9348E-3;
        double initialBeta = 1.0;
        double interBeta   = 0.85;

        double time_step=0.1;
        int n_time_steps=3000;                  // As written, output written to maximum 100000000 s$

        double tol = 1.E-8;                     // Used to determine whether linearisation converges$
  


Command to compile:
g++ ConservativeIFF.cpp InOut.cpp subroutines.cpp trimatrix.cpp

Command to run:
./a.out

Command to view output:
nano results.txt
  

  
Example Output:

Time	Interface position
Time	Interface Position
0	1
10	0.841793
20	0.774097
30	0.722167
40	0.678392
...
