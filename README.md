<p align="center">
   <img src="https://github.com/amalbavera/pysco/blob/main/docs/pysco_logo.png" alt="pySCO logo" align="center" width="350px"/>
</p>

Please ensure that the following external libraries are installed beforehand:

- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)

We use the [International System of Units](https://www.nist.gov/pml/owm/metric-si/si-units) as follows:

| Property              | Symbol                    | Unit                          |
| :-------------------: | :-----------------------: | :---------------------------: |
| Temperature           | $T$                       | K                             |
| Spin-Crossover Energy | ${\Delta E}_\mathrm{sco}$ | kJ mol $^{-1}$                |
| Enthalpy Change       | $\Delta H$                | kJ mol $^{-1}$                |
| Entropy Change        | $\Delta S$                | J mol $^{-1}$ K $^{-1}$       |
| Gibbs Free Energy     | $G$                       | kJ mol $^{-1}$                |
| Interaction Parameter | $\Gamma$                  | kJ mol $^{-1}$                |
| High-Spin Population  | $n_\mathrm{HS}$           | $0 \leq n_\mathrm{HS} \leq 1$ |

---
---

# Table of Contents

1. [Reading Outputs from Electronic Structure Codes](#1-reading-outputs-from-electronic-structure-codes)
    1. [VASP](#11-vasp)
    2. [Quantum Espresso](#12-quantum-espresso)
    3. [Gaussian](#13-gaussian)
    4. [Orca](#14-orca)
    5. [NWChem](#15-nwchem)
    6. [Passing additional options](#16-passing-additional-options)
        1. [Magnetization per Atom](#161-magnetization-per-atom)
        2. [Jahn-Teller Distortion](#162-jahn-teller-distortion)
        3. [Spin-Orbit Coupling](#163-spin-orbit-coupling)
2. [Computing Properties](#2-computing-properties)
    1. [Spin-Crossover Energy](#21-spin-crossover-energy)
    2. [Transition Temperature](#22-transition-temperature)
    3. [High-Spin Population as function of Temperature](#23-high-spin-population-as-function-of-temperature)
    4. [Gibbs Free Energy as Function of Temperature](#24-gibbs-free-energy-as-function-of-temperature)
    5. [Fitting the Interaction Parameter](#25-fitting-the-interaction-parameter)

---
---

# 1. Reading Outputs from Electronic Structure Codes

Currently, the library is compatible with [VASP](#11-vasp), [Gaussian](#13-gaussian), [Orca](#14-orca), and [NWChem](#15-nwchem). We are working on adding compatiblity with [Quantum Espresso](#12-quantum-espresso).

---

## 1.1 VASP

Since VASP produces several output files, each collecting different data, we need to locate the folder with the outputs and pass it to pySCO. Before proceeding any further, please keep in mind that this is a three-step process, namely,

1. relax the geometry using a $k$-point mesh,
2. run a single-point energy calculation with a $k$-point mesh, and
3. compute the vibrational frequencies with a $q$-point mesh.

Here we are interested in the EIGENVAL, POSCAR, and OUTCAR files from step (2), and in the OUTCAR file from step (3). There are two additional files we need to prepare for every spin state, namely, MAGCAR and VIBCAR. The first one simply contains the analogous magnetization as reported in the OUTCAR from step (2), but with the correct rounded magnetic moments, e.g.

```text

 magnetization (x)

# of ion       s       p       d       tot
------------------------------------------
    1        0.000   0.000   4.000   4.000
    2        0.000   0.000   4.000   4.000
    3        0.000   0.000   0.000   0.000
    ...
  256        0.000   0.000   0.000   0.000
--------------------------------------------------
tot          0.000   0.000   8.000   8.000
```

where the first two atoms are assigned four unpair electrons each, whereas the rest show no magnetic moment. The second file, VIBCAR, is the OUTCAR that contains the harmonic vibrational frequencies computed with the chosen set of $q$-points in step (3). We can either rename the original OUTCAR file from step (3) to VIBCAR, or save time by running the following command from the console

```bash
grep 2PiTHz OUTCAR > VIBCAR
```

Now we can pass the location to the directory that contains the OUTCAR, EIGENVAL, POSCAR, MAGCAR, and VIBCAR files. Remember that we need to do this for each spin state.

```python
from pysco import read

low_spin  = read.vasp( "path/to/the/low/spin/state/directory" )
high_spin = read.vasp( "path/to/the/high/spin/state/directory" )
```

---

## 1.2 Quantum Espresso

We currently are working on making pySCO compatible with Quantum Espresso files.

---

## 1.3 Gaussian

There are several approaches to generating the output files from Gaussian. The most straightforward is to run a geometry relaxation followed by a frequency analysis in the same job. However, this usually is discuraged by most exprienced users and, instead, a two-step approach is preferred, where we first do the geometry relaxation and then, on a separate run, compute the frequency analysis. Here we will favor this second approach and hence the output of interest is that for the frequency analysis. We can pass its location as

```python
from pysco import read

low_spin  = read.gaussian( "path/to/the/low/spin/state/gaussian/output.log" )
high_spin = read.gaussian( "path/to/the/high/spin/state/gaussian/output.log" )
```

Remember that you need to do this for each spin state. The file extension is irrelevant for as long as you include the complete file name.

---

## 1.4 Orca

There are several approaches to generating the output files from Orca. The most straightforward is to run a geometry relaxation followed by a frequency analysis in the same job. However, this usually is discuraged by most exprienced users and, instead, a two-step approach is preferred, where we first do the geometry relaxation and then, on a separate run, compute the frequency analysis. Here we will favor this second approach and hence the output of interest is that for the frequency analysis. We can pass its location as

```python
from pysco import read

low_spin  = read.orca( "path/to/the/low/spin/state/orca/output.out" )
high_spin = read.orca( "path/to/the/high/spin/state/orca/output.out" )
```

Remember that you need to do this for each spin state. The file extension is irrelevant for as long as you include the complete file name.

**Advise to the user**: Orca sometimes truncates the number of molecular orbitals written in the output, which is not helpful in our case, therefore we recommend including `! PrintMOs` in your input to ensure that all orbitals are printed.

---

## 1.5 NWChem

There are several approaches to generating the output files from NWChem. The most straightforward is to run a geometry relaxation followed by a frequency analysis in the same job. However, this usually is discuraged by most exprienced users and, instead, a two-step approach is preferred, where we first do the geometry relaxation and then, on a separate run, compute the frequency analysis. Here we will favor this second approach and hence the output of interest is that for the frequency analysis. We can pass its location as

```python
from pysco import read

low_spin  = read.nwchem( "path/to/the/low/spin/state/nwchem/output.out" )
high_spin = read.nwchem( "path/to/the/high/spin/state/nwchem/output.out" )
```

Remember that you need to do this for each spin state. The file extension is irrelevant for as long as you include the complete file name.

---

## 1.6 Passing Additional Options
The [read](#1-reading-outputs-from-electronic-structure-codes) module also allows for the inclusion of additional options when reading output files. These are specified as keyword arguments, namely, [magnetization](#161-magnetization-per-atom), [jahnteller](#162-jahn-teller-distortion), and [orbit](#163-spin-orbit-coupling). These keywords are not mutually exclusive and hence may be included simultaneously.

### 1.6.1 Magnetization per Atom

`magnetization: list[int]`. We refer to magnetization, $\zeta$, as the difference between spin-up electrons, $`N_\uparrow`$, and spin-down electrons, $`N_\downarrow`$. Hence, $`\zeta = N_\uparrow - N_\downarrow`$. By default, pySCO will search for the multiplicity to determine the magnetization and assign it to any Cr, Mn, Fe, or Co in the structure. This choice certainly is convenient for mononuclear metal complexes, but needs to be handled properly for multicenter systems. The `magnetization` argument can therefore be used to override the default behavior and specify the desired magnetic moment for each metal center, or any other atom in particular.

- It is a list of integers with the same size as the total number of atoms in the metal-complex. 
- The order of the elements in this list must be consistent with the ordering of the atoms in the geometry.

The following is an example using the Orca reader,

```python
from pysco import read

low_spin_mag  = 2*[5.0] + 98*[0.0]
high_spin_mag = 2*[1.0] + 98*[0.0]

low_spin  = read.orca( "path/to/the/low/spin/state/orca/output.out",  magnetization=low_spin_mag )
high_spin = read.orca( "path/to/the/high/spin/state/orca/output.out", magnetization=high_spin_mag )
```

It is important to notice that the `magnetization` argument is not needed for [VASP](#11-vasp) because we already provide the MAGCAR. Keep in mind that pySCO will use `magnetization` instead of the MAGCAR file should it be specified.

### 1.6.2 Jahn-Teller Distortion

`jahnteller: int`. By default pySCO assumes no Jahn-Teller distortions. The argument `jahnteller` may be used to specify the number of distortions in the structure. The following is an example using the Gaussian reader,

```python
from pysco import read

low_spin  = read.gaussian( "path/to/the/low/spin/state/gaussian/output.log",  jahnteller=1 )
high_spin = read.gaussian( "path/to/the/high/spin/state/gaussian/output.log", jahnteller=3 )
```

### 1.6.3 Spin-Orbit Coupling

`orbit: int`. By default pySCO assumes no spin-orbit coupling. The argument `orbit` may be used to specify the orbital angular momentum, $L$, to compute the entropic contribution $`{\Delta S}_\mathrm{orb}`$ in the form

```math
{\Delta S}_\mathrm{orb} = k_B\,N_A\,\ln\left[ \frac{2\,L_\mathrm{HS} + 1}{2\,L_\mathrm{LS} + 1} \right]
```

where LS and HS label the low- and high-spin state, respectively. The following is an example using the NWChem reader,

```python
from pysco import read

low_spin  = read.nwchem( "path/to/the/low/spin/state/nwchem/output.out",  orbit=1 )
high_spin = read.nwchem( "path/to/the/high/spin/state/nwchem/output.out", orbit=3 )
```
---
---

# 2. Computing Properties

---

## 2.1 Spin-Crossover Energy

The energy difference between the high-spin and low-spin states, $`{E}_\mathrm{HS}`$ and $`{E}_\mathrm{LS}`$, respectively, plus the zero-point energy correction, $`{\Delta E}_\mathrm{zpe}`$, defines the spin-crossover energy as

```math
{\Delta E}_\mathrm{sco} = {E}_\mathrm{HS} - {E}_\mathrm{LS} + {\Delta E}_\mathrm{zpe}
```

There are two optional arguments that modify the behavior of the function:

- `metals: tuple(str)`. By default $`{\Delta E}_\mathrm{sco}`$ is normalized to the total number of spin conversion centers that is determined by comparing which transition metals undergo a spin transition. This list of metals defaults to `metals=('Cr', 'Mn', 'Fe', 'Co)` and can be modified by the user.
- `centers: float`. The automatic normalization can also be overruled by means of the optional keyword `centers` that fixes the total number of spin-crossover centers.

```python
from pysco import thermo

# Default, in kJ/mol
Esco = thermo.spin_crossover_energy( ls=low_spin, hs=high_spin )

# Modify the metals of interest for automatic normalization
Esco = thermo.spin_crossover_energy( ls=low_spin, hs=high_spin, metals=('Fe', 'Mn') )

# Specify the number of spin conversion centers
Esco = thermo.spin_crossover_energy( ls=low_spin, hs=high_spin, centers=2.0 )
```

---

## 2.2 Transition Temperature

The equilibrium temperature at which the populations of the high-spin and low-spin states are equal defines the transition temperature, $`T_{1/2}`$ as 

```math
T_{1/2} = \frac{\Delta H}{\Delta S}
```

where $`\Delta H`$ and $`\Delta S`$ are the changes in enthalpy and entropy between, respectively.

There are four optional arguments that modify the behavior of the function:

- `metals: tuple(str)`. By default the contribution $`{\Delta E}_\mathrm{sco}`$ is normalized to the total number of spin conversion centers that is determined by comparing which transition metals undergo a spin transition. This list of metals defaults to `metals=('Cr', 'Mn', 'Fe', 'Co)` and can be modified by the user.
- `centers: float`. The automatic normalization can also be overruled by means of the optional keyword `centers` that fixes the total number of spin-crossover centers.
- `dH_and_dS: bool`. If set to `True`, the function will return both $`\Delta H`$ and $`\Delta S`$ in addition to $`T_{1/2}`$.
- `guess: float`. A rough estimate of the transition temperature can be provided via the optional keyword `guess` to accelerate the numerical procedure. It defaults to `guess=298.15` Kelvin.

```python
from pysco import thermo

# Default, in Kelvin
Thalf = thermo.transition_temperature( ls=low_spin, hs=high_spin )

# Additionally return both ΔH and ΔS, in kJ/mol and J/(mol K), respectively
Thalf, dH, dS = thermo.transition_temperature( ls=low_spin, hs=high_spin, dH_and_dS=True )

# Modify the initial temperature guess, in Kelvin
Thalf = thermo.transition_temperature( ls=low_spin, hs=high_spin, guess=120.0 )

# Multiple choices
Thalf, dH, dS = thermo.transition_temperature( ls=low_spin, hs=high_spin, centers=2.0, guess=120.0, dH_and_dS=True)
```

---

## 2.3 High-Spin Population as Function of Temperature

Returns a NumPy array where the first column is the temperature, $T$, in Kelvin and the second column is the high-spin population, $`n_\mathrm{HS}`$.

- `metals: tuple(str)`. By default the contribution $`{\Delta E}_\mathrm{sco}`$ is normalized to the total number of spin conversion centers that is determined by comparing which transition metals undergo a spin transition. This list of metals defaults to `metals=('Cr', 'Mn', 'Fe', 'Co)` and can be modified by the user.
- `centers: float`. The automatic normalization can also be overruled by means of the optional keyword `centers` that fixes the total number of spin-crossover centers.
- `points: list|int`. The number of temperature points at which the high-spin population is evaluated. If not specified, a default `points=256` is used.
- `interaction: float`. The phenomenological interaction parameter, $\Gamma$, as defined in the [Slichter-Drickamer model](https://doi.org/10.1063/1.1677511). It defaults to `interaction=0.0`

```python
from pysco import thermo

# Default
T_and_nHS = thermo.high_spin_population( ls=low_spin, hs=high_spin )

# Modify the number of points
T_and_nHS = thermo.high_spin_population( ls=low_spin, hs=high_spin, points=128 )

# Specify a list of points, in Kelvin
T_and_nHS = thermo.high_spin_population( ls=low_spin, hs=high_spin, points=[100, 150, 200, 250] )

# Fix the phenomenological interaction parameter, in kJ/mol
T_and_nHS = thermo.high_spin_population( ls=low_spin, hs=high_spin, interaction=2.5 )
```
---

## 2.4 Gibbs Free Energy as Function of Temperature

Returns a NumPy array where the first column is the high-spin population, $`n_\mathrm{HS}`$, and the second column is the Gibss free energy, $G$, in kJ/mol. Here we must **specify our choice of temperature for the isotherm**.

- `metals: tuple(str)`. By default the contribution $`{\Delta E}_\mathrm{sco}`$ is normalized to the total number of spin conversion centers that is determined by comparing which transition metals undergo a spin transition. This list of metals defaults to `metals=('Cr', 'Mn', 'Fe', 'Co)` and can be modified by the user.
- `centers: float`. The automatic normalization can also be overruled by means of the optional keyword `centers` that fixes the total number of spin-crossover centers.
- `points: list|int`. The number of points to evaluate the high-spin population is evaluated. If not specified, a default `points=256` is used.
- `temperature: float`. The temperature at which the Gibbs free energy is evaluated. It defaults to `temperature=0.0` Kelvin.
- `interaction: float`. The phenomenological interaction parameter, $\Gamma$, as defined in the [Slichter-Drickamer model](https://doi.org/10.1063/1.1677511). It defaults to `interaction=0.0`

```python
from pysco import thermo

# Default
nHS_and_G = thermo.gibbs_free_energy( ls=low_spin, hs=high_spin, temperature=120.0 )

# Modify the number of points
nHS_and_G = thermo.gibbs_free_energy( ls=low_spin, hs=high_spin, temperature=120.0, points=128 )

# Specify a list of points
nHS_and_G = thermo.gibbs_free_energy( ls=low_spin, hs=high_spin, temperature=120.0, points=[0.0, 0.25, 0.5, 0.75, 1.0] )

# Fix the phenomenological interaction parameter, in kJ/mol
nHS_and_G = thermo.gibbs_free_energy( ls=low_spin, hs=high_spin, temperature=120.0, interaction=2.5 )

# Multiple choices
nHS_and_G = thermo.gibbs_free_energy( ls=low_spin, hs=high_spin, temperature=120.0, points=128, interaction=2.5 )
```
---

## 2.5 Fitting the Interaction Parameter

The phenomenological interaction parameter, $\Gamma$, was proposed by Slichter and Drickamer to account for cooperative effects in spin-crossover materials. We can consider it as an addition second-order contribution to the Gibbs free energy in the form

```math
G = (1 - n_\mathrm{HS})\,G_\mathrm{LS} + n_\mathrm{HS}\,G_\mathrm{HS} + \Gamma\,n_\mathrm{HS}\,(1 - n_\mathrm{HS}) - T\,S_\mathrm{mix}
```

where LS and HS label the low- and high-spin state, respectively, and $`S_\mathrm{mix}`$ is the mixing entropy given by $`S_\mathrm{mix} = -k_B\,N_A\,(\; n_\mathrm{HS} \ln[n_\mathrm{HS}] + (1 - n_\mathrm{HS}) \ln[1 - n_\mathrm{HS}] \;)`$. Keep in mind that we need to sample a series of spin configurations within the interval $`0 < n_\mathrm{HS} < 1`$ to compute this interaction parameter. This set of spin configurations are passed to pySCO with the keyword `ms:list` that contains the objects generated by the [read](#1-reading-outputs-from-electronic-structure-codes) module.

- `metals: tuple(str)`. By default the contribution $`{\Delta E}_\mathrm{sco}`$ is normalized to the total number of spin conversion centers that is determined by comparing which transition metals undergo a spin transition. This list of metals defaults to `metals=('Cr', 'Mn', 'Fe', 'Co)` and can be modified by the user.
- `centers: float`. The automatic normalization can also be overruled by means of the optional keyword `centers` that fixes the total number of spin-crossover centers.
- `least_squares: bool`. If set to `True`, the fitting procedure will use a least-squares approach and returns both the interaction parameter and coefficient of determination, $`R^2`$. This is the default.
- `boltzmann: bool`. If set to `True`, the fitting procedure will use a Boltzmann ensemble approach as in H. Paulsen [Magnetochemistry](https://doi.org/10.3390/magnetochemistry2010014) 2, 14 (2016). This option still is at an experimental stage, and returns both the interaction parameter and weight for each spin configuration.

```python
from pysco import thermo

# List of spin configurations
mixed_spin = [ spin_config_1, spin_config_2, spin_config_3, ... ]

# Default, in kJ/mol
interaction, R2 = thermo.interaction_parameter( ls=low_spin, hs=high_spin, ms=mixed_spin )

# Modify the metals of interest for automatic normalization
interaction, R2 = thermo.interaction_parameter( ls=low_spin, hs=high_spin, ms=mixed_spin, metals=('Fe', 'Mn') )

# Specify the number of spin conversion centers
interaction, R2 = thermo.interaction_parameter( ls=low_spin, hs=high_spin, ms=mixed_spin, centers=2.0 )

# Specify the fitting procedure
interaction, weights = thermo.interaction_parameter( ls=low_spin, hs=high_spin, ms=mixed_spin, boltzmann=True )
```
---
