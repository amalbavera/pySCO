#############################################################################################################################

##       #### ########  ########     ###    ########  #### ########  ######
##        ##  ##     ## ##     ##   ## ##   ##     ##  ##  ##       ##    ##
##        ##  ##     ## ##     ##  ##   ##  ##     ##  ##  ##       ##
##        ##  ########  ########  ##     ## ########   ##  ######    ######
##        ##  ##     ## ##   ##   ######### ##   ##    ##  ##             ##
##        ##  ##     ## ##    ##  ##     ## ##    ##   ##  ##       ##    ##
######## #### ########  ##     ## ##     ## ##     ## #### ########  ######

#############################################################################################################################

import numpy             as np
import multiprocessing   as mp

from scipy.optimize import fsolve
#
### local libraries
#
from .constants import h_over_k, k_times_n
from .constants import ev_to_joule

#############################################################################################################################

##     ##  #######  ########  ##     ## ##       ########  ######
###   ### ##     ## ##     ## ##     ## ##       ##       ##    ##
#### #### ##     ## ##     ## ##     ## ##       ##       ##
## ### ## ##     ## ##     ## ##     ## ##       ######    ######
##     ## ##     ## ##     ## ##     ## ##       ##             ##
##     ## ##     ## ##     ## ##     ## ##       ##       ##    ##
##     ##  #######  ########   #######  ######## ########  ######

#############################################################################################################################
### Extract thermodynamic data from spin-states

def extract_data(ls:any=None, hs:any=None, thermo:any=None, return_nhs=False)-> list:
    
    sco    = thermo.sco(ls=ls, hs=hs, zpe=True)
    
    s_spin = thermo.s_spin(ls=ls, hs=hs)
    
    pv     = thermo.pdv(ls=ls, hs=hs)
    
    frq    = np.vstack([ls.frequencies, hs.frequencies])
    
    eig    = np.vstack([ls.orbitalenergies, hs.orbitalenergies])
    
    rot    = np.stack(( (ls.rotationalsymm,ls.rotationaltemp), (hs.rotationalsymm,hs.rotationaltemp) ))

    if return_nhs:
        nhs = thermo.nhs(ls=ls, hs=hs)
        return [sco, s_spin, pv, frq, eig, rot], nhs
    
    return [sco, s_spin, pv, frq, eig, rot]

#############################################################################################################################
### Ranges for numerical or direct methods

def set_range(sigmoid:tuple=(-10.0,10.0), temperatures=(5.0e1, 3.5e2), points:int|list=256, rtype:str="none") -> np.ndarray:
    
    if isinstance(points,list) or isinstance(points,np.ndarray):
        points      = np.asarray(points)
        points_list = np.linspace(sigmoid[0], sigmoid[1], num=points.size)

    else:
        points_list = np.linspace(sigmoid[0], sigmoid[1], num=points)

    sigmoid_points = 1.0/(1.0 + np.exp(-points_list))
    
    if rtype.lower() == "nhs":
        return sigmoid_points

    start = points[0]  if isinstance(points,np.ndarray) else temperatures[0]
    stop  = points[-1] if isinstance(points,np.ndarray) else temperatures[1]

    iteration_points = start*np.flip(sigmoid_points) + stop*sigmoid_points

    return iteration_points, 0.5*(iteration_points[-1] - iteration_points[0])

#############################################################################################################################
### Solve numeric (nhs vs T) or (G vs nhs)

def solve_numeric(nhs:float=0.0, data:any=None, idx:int=0, guess:float=298.15, interaction:float=0.0, thermo:any=None, regsolmod:bool=False) -> tuple:
    
    args = data.zipit(idx=idx, interaction=interaction, nhs=nhs)
    
    T  = fsolve(T_from_nhs, guess, args=args, xtol=1e-3)[0]
    
    if not regsolmod:
        return T, nhs
    
    dG = delta_gibbs(T=T, **data.zipit(idx=idx, interaction=interaction))
    G  = nhs*dG + nhs*(1.0 - nhs)*interaction - T*thermo.s_mix(nhs=nhs)
    
    return nhs, G

#############################################################################################################################
### Solve direct (nhs vs T) or (G vs nhs)

def solve_direct(nhs:float=0.0, T:float=298.15, data:any=None, idx:int=0, guess:float=298.15, interaction:float=0.0, thermo:any=None, regsolmod:bool=False) -> tuple:
    
    kwargs = data.zipit(idx=idx, interaction=interaction)
    
    if regsolmod:
        dG = delta_gibbs(T=guess, **kwargs)
        G  = nhs*dG + nhs*(1.0 - nhs)*interaction - guess*thermo.s_mix(nhs=nhs)

        return nhs, G
    
    dG = delta_gibbs(T=T, **kwargs)
            
    if interaction < 0.0:
        guess = 1e-8 if T <= guess else 1-1e-4
        nhs   = fsolve(nhs_from_T, guess, args=(("dG",dG), ("T",T), ("interaction",interaction)), xtol=1e-4)[0]

    else:
        nhs    = nhs_from_T(0.0, ("dG",dG), ("T",T), ("interaction",interaction))
    
    return T, nhs
    
#############################################################################################################################
### High-spin poulation vs temperature

def solve_nhs_vs_temperature(data:any=None, idx:int=0, point:float=0.0, interaction:float=0.0, guess:float=298.15, rescale:list=[0.0, 1.0], option:str="none") -> list:
    
    if   option in ("n", "num", "numeric"):
        T, nhs = solve_numeric(nhs=point, data=data, idx=idx, guess=guess, interaction=interaction)
        
    elif option in ("d", "dir", "direct"):
        T, nhs = solve_direct(T=point, data=data, idx=idx, guess=guess, interaction=interaction)
        
    return [T, rescale[0] + (rescale[1] - rescale[0])*nhs]

#############################################################################################################################
### Gibbs free energy vs High-spin poulation

def solve_gibbs_vs_nhs(data:any=None, idx:int=0, point:float=0.0, interaction:float=0.0, guess:float=298.15, thermo:any=None, option:str="numeric") -> float:
    
    if   option in ("n", "num", "numeric"):
        nhs, G = solve_numeric(nhs=point, data=data, idx=idx, guess=guess, interaction=interaction, thermo=thermo, regsolmod=True)
        
    elif option in ("d", "dir", "direct"):
        nhs, G = solve_direct(nhs=point, data=data, idx=idx, guess=guess, interaction=interaction, thermo=thermo, regsolmod=True)
        
    return [nhs, G]

#############################################################################################################################
### Calculate Boltzmann factor and microscopic interaction parameter

def solve_factor_and_interaction(ls=None, ms=None, states:any=None, guess:float=298.15, thermo:any=None, interaction:float=0.0, option:str="scf") -> np.ndarray:
        
    step = extract_data(ls=ls, hs=ms, thermo=thermo)

    sco  = step[0]
        
    nhs    = thermo.nhs(ls=ls, hs=ms)
    
    if (nhs == 1.0) or (nhs == 0.0): raise RuntimeError(f"high-spin population for {ms} is {nhs:.2f}")
    
    ms.step  = step
        
    hschoice = 0.5 if option == "fix" else nhs

    hlargs   = states.zipit(idx=0, nhs=hschoice, interaction=interaction)

    argsdict = dict(i for i in hlargs)
    
    T   = fsolve(T_from_nhs, guess, args=hlargs, xtol=1e-3)[0]
        
    KbT = k_times_n*T
        
    e_ml_hl          = sco - nhs*argsdict["E_sco"]
    
    boltzmann_factor = np.exp(-e_ml_hl/KbT)
        
    interaction_c    = e_ml_hl/(nhs*(1.0 - nhs))
        
    return np.asarray([sco, nhs, boltzmann_factor, interaction_c])

#############################################################################################################################
### Fermi-Dirac Entropy

def electronic_entropy(energy:np.ndarray=None, T:float=298.15, sigma:float=1e-2, thresh:float=700.0) -> float:

    size, KbT     = len(energy), k_times_n*T/ev_to_joule

    fermi, dummy  = np.zeros(size), np.zeros(size)

    dos = np.ones(size) - sigma*energy**2
    dos = dos*(dos>0)
    
    Ei_over_KbT = energy/KbT
    
    fermi[Ei_over_KbT<-thresh] = 1.0
    
    fermi[Ei_over_KbT>-thresh] = 0.0
    
    idx = np.where( (Ei_over_KbT>-thresh)&(Ei_over_KbT<thresh) )
    
    fermi[idx] = 1.0/(1.0 + np.exp(Ei_over_KbT[idx]))
    
    idx = np.where( (fermi>0.0)&(fermi<1.0) )
    
    dummy[idx] = dos[idx]*(fermi[idx]*np.log(fermi[idx]) + (1.0 - fermi[idx])*np.log(1.0 - fermi[idx]))
    
    idx = np.logical_not(np.isnan(dummy))

    S_elec = -k_times_n*np.trapz(energy[idx], dummy[idx])

    return S_elec

#############################################################################################################################
### Determine temperature from high-spin population

def T_from_nhs(T:float, *args:tuple) -> float:

    kwargs = dict(i for i in args)

    nhs = kwargs["nhs"]

    dG = delta_gibbs(T=T, E_sco=kwargs["E_sco"], S_spin=kwargs["S_spin"], PdV=kwargs["PdV"],
                     frq=kwargs["frq"], eig=kwargs["eig"], rot=kwargs["rot"])

    difference = k_times_n*T*np.log((1.0 - nhs)/nhs) - dG - (1.0 - 2.0*nhs)*kwargs["interaction"]

    return difference

#############################################################################################################################
### Determine high-spin population at given temperature

def nhs_from_T(nhs:float, *args:tuple) -> float:
    
    kwargs = dict(i for i in args)
    
    exp_thresh, log_thresh = 7e2, 1e-16
    
    dG, T, interaction  = kwargs["dG"], kwargs["T"], kwargs["interaction"]
    
    nhs = min(max(nhs, log_thresh), 1.0-log_thresh)
    
    if interaction: return k_times_n*T*np.log((1.0 - nhs)/nhs) - dG - (1.0 - 2.0*nhs)*interaction

    exp_arg = dG/(k_times_n*T)

    if np.isnan(exp_arg) or np.isinf(exp_arg): raise RuntimeError("high-spin population undefined")

    if -exp_thresh < exp_arg < exp_thresh: return 1.0/(1.0 + np.exp(exp_arg))
    
    nhs = 1.0 if exp_arg <= -exp_thresh else 0.0
    
    return nhs

#############################################################################################################################
### Calculate delta Gibbs at given temperature

def delta_gibbs(T:float=298.15, E_sco:float=0.0, S_spin:float=0.0, PdV:float=0.0, frq:np.ndarray=None,
                eig:np.ndarray=None, rot:np.ndarray=None, thresh:float=7e2, dH_and_dS=False) -> float:

    if T <= 0.0: return -np.inf
    if T is np.nan: return -np.inf
    
    hkv_ls   = h_over_k*frq[0,:]
    hkv_hs   = h_over_k*frq[1,:]
    
    hkvT_ls  = hkv_ls/T
    hkvT_hs  = hkv_hs/T
    
    exp_ls   = np.exp(-hkvT_ls)
    exp_hs   = np.exp(-hkvT_hs)
    
    oexp_ls  = 1.0 - exp_ls
    oexp_hs  = 1.0 - exp_hs
    
    E_vib_ls = k_times_n*hkv_ls*exp_ls/oexp_ls
    E_vib_hs = k_times_n*hkv_hs*exp_ls/oexp_hs
    
    S_vib_ls = k_times_n*(hkvT_ls)*exp_ls/oexp_ls - k_times_n*np.log(oexp_ls)
    S_vib_hs = k_times_n*(hkvT_hs)*exp_hs/oexp_hs - k_times_n*np.log(oexp_hs)
    
    S_ele_ls = electronic_entropy(energy=np.sort(eig[0,:]), T=T, thresh=thresh)
    S_ele_hs = electronic_entropy(energy=np.sort(eig[1,:]), T=T, thresh=thresh)

    S_fermi  = S_ele_hs - S_ele_ls

    dU       = E_sco + np.sum(E_vib_hs - E_vib_ls)

    dS       = S_spin + S_fermi + np.sum(S_vib_hs - S_vib_ls) + k_times_n
    
    if not any( i is None for i in rot.ravel() ):
        rot_symm_ls, rot_temp_ls = rot[0,:]
        rot_symm_hs, rot_temp_hs = rot[1,:]
                
        rot_theta_ls = rot_symm_ls*np.sqrt(rot_temp_ls*T*T*T)
        rot_theta_hs = rot_symm_hs*np.sqrt(rot_temp_hs*T*T*T)
        
        S_rot_ls = k_times_n*np.log(rot_theta_ls)
        S_rot_hs = k_times_n*np.log(rot_theta_hs)
                
        dS  += S_rot_hs - S_rot_ls

    dH = dU + PdV
    dG = dH - dS*T
        
    if not dH_and_dS:
        return dG
    
    return dG, dH, dS

#############################################################################################################################
### Calculate the transition temperature

def calculate_thalf(data:any=None, guess=298.15, interaction:float=0.0, centers:any=None) -> float:

    args   = data.zipit(idx=0, interaction=interaction, nhs=0.5)
    kwargs = dict(i for i in args)

    thalf = fsolve(T_from_nhs, guess, args=args, xtol=1e-3)[0]

    _, dH, dS = delta_gibbs(T=thalf, E_sco=kwargs["E_sco"], S_spin=kwargs["S_spin"], PdV=kwargs["PdV"],
                            frq=kwargs["frq"], eig=kwargs["eig"], rot=kwargs["rot"], dH_and_dS=True)
    return thalf, dH, dS

#############################################################################################################################
### Calculate the relative high-spin population

def calculate_nhs(option:str="direct", guess:float=298.15, points:int|list=256, data:any=None, interaction:float=0.0) -> np.ndarray:

    cpus = mp.cpu_count()

    if option == "numeric":
        points        = set_range(points=points, rtype="nhs")

    else:
        points, guess = set_range(points=points, rtype="tmp")

    for idx in range(data.crossovers):
        
        resize   = (0.0, data.nhs[idx]) if idx == 0 else (data.nhs[idx-1], data.nhs[idx])
        
        pool     = mp.Pool(cpus)
        
        mparg    = [(data, idx, point, interaction, guess, resize, option) for point in points]

        T_and_hs = pool.starmap_async(solve_nhs_vs_temperature, mparg).get()
        
        pool.close()
        
    return np.asarray(T_and_hs)

#############################################################################################################################
### Numerical or direct methods

def calculate_rsm(option:str="numeric", guess:float=298.15, points:list|int=256, data:any=None, interaction:float=0.0, thermo:any=None) -> np.ndarray:

    cpu_avail = mp.cpu_count()
    
    points    = set_range(points=points, rtype="nhs")
                    
    pool  = mp.Pool(cpu_avail)
    
    mparg = [(data, 0, point, interaction, guess, thermo, option) for point in points]

    nhs_and_gibbs = pool.starmap_async(solve_gibbs_vs_nhs, mparg).get()
                        
    pool.close()

    return np.asarray(nhs_and_gibbs)

#############################################################################################################################
### Calcualte interaction parameter for a set of configurations

def calculate_interaction(ls:any=None, hs:any=None, ms:any=None, interaction:float=0.0, states:any=None, thermo:any=None,
                          option:str="scf") -> np.ndarray:

    cpu_avail = mp.cpu_count()

    ms_states = len(ms)

    cpus  = ms_states if cpu_avail > ms_states else cpu_avail
    pool  = mp.Pool(cpus)

    mparg = [(ls, i, states, 298.15, thermo, interaction, option) for i in ms]

    interaction_data = pool.starmap_async(solve_factor_and_interaction, mparg).get()

    return np.asarray(interaction_data)

#############################################################################################################################
