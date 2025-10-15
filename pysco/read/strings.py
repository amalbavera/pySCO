#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################


#############################################################################################################################

 ######  ##        #######  ########     ###    ##        ######  
##    ## ##       ##     ## ##     ##   ## ##   ##       ##    ## 
##       ##       ##     ## ##     ##  ##   ##  ##       ##       
##   ### ##       ##     ## ########  ##     ## ##        ######  
##    ## ##       ##     ## ##     ## ######### ##             ## 
##    ## ##       ##     ## ##     ## ##     ## ##       ##    ## 
 ######  ########  #######  ########  ##     ## ########  ######

#############################################################################################################################
### Output serach strings for VASP files

vasp_total_energy_str:str     = "free  energy   TOTEN  ="
vasp_fermi_energy_str:str     = "Fermi energy:"
vasp_cell_volume_str:str      = "volume of cell :"

vasp_mag_start_str:str        = "magnetization (x)"
vasp_mag_stop_str:str         = "-----"

vasp_frequency_str:str        = "2PiTHz"
vasp_imag_freq_str:str        = "f/i"

#############################################################################################################################
### Output serach strings for GAUSSIAN files

#gaussian_single_point_str:str = "A.U. after    1 cycles"
gaussian_total_energy_str:str = "SCF Done:"
gaussian_multip_str:str       = "Multiplicity ="

gaussian_up_occ_str:str       = "Alpha  occ. eigenvalues --"
gaussian_up_vir_str:str       = "Alpha virt. eigenvalues --"
gaussian_dw_occ_str:str       = "Beta  occ. eigenvalues --"
gaussian_dw_vir_str:str       = "Beta virt. eigenvalues --"

gaussian_pop_start_str:str    = "APT charges:"
gaussian_pop_stop_str:str     = "Sum of APT charges ="

gaussian_frequency_str:str    = "Frequencies --"
gaussian_vib_temps_str:str    = "Vibrational temperatures"

gaussian_rot_sym_str:str      = "Rotational symmetry number"
gaussian_rot_temp_str:str     = "Rotational temperatures"

#############################################################################################################################
### Output serach strings for NWCHEM files

nwchem_single_point_str:str   = "NWChem Nuclear Hessian and Frequency Analysis"
nwchem_total_energy_str:str   = "Total DFT energy ="
nwchem_multip_str:str         = "Spin multiplicity:"

nwchem_up_orb_str:str         = "DFT Final Alpha Molecular Orbital Analysis"
nwchem_dw_orb_str:str         = "DFT Final Beta Molecular Orbital Analysis"
nwchem_orb_str:str            = "DFT Final Molecular Orbital Analysis"
nwchem_occ_str:str            = "Occ="

nwchem_geom_str:str           = "XYZ format geometry"
nwchem_overlap_str:str        = "alpha - beta orbital overlaps"

nwchem_frequency_str:str      = "P.Frequency"
nwchem_rot_sym_str:str        = "symmetry #  ="
nwchem_rot_temp_str:str       = "Rotational Constants"

#############################################################################################################################
### Output serach strings for ORCA files

orca_total_energy_str:str     = "FINAL SINGLE POINT ENERGY"
orca_multip_str:str           = "Multiplicity       :"
orca_general_stop_str:str     = "-----"

orca_orbital_str:str          = "ORBITAL ENERGIES"
orca_up_orb_str:str           = "SPIN UP ORBITALS"
orca_dw_orb_str:str           = "SPIN DOWN ORBITALS"

orca_pop_start_str:str        = "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS"
orca_pop_stop_str:str         = "Sum of atomic charges"

orca_frequency_str:str        = "VIBRATIONAL FREQUENCIES"
orca_rot_sym_str:str          = "Symmetry Number:"
orca_rot_const_str:str        = "Rotational constants in cm-1:"
