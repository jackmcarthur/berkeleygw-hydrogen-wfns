import numpy as np
from scipy.special import sph_harm
import h5py as h5
import sys

# Analytical wavefunctions in momentum space, written over a provided BerkeleyGW WFN.h5 file.
# written by Jack McArthur, last edited 7/26/2024

# radial function in momentum space:
# 1s: -1/(1+p^2)^2
# 2s: (2-8p^2)/(1+4p^2)^3
# 2p: -2i*p/(1+4p^2)^3
# 3s: -3(1-30p^2+81p^4)/(1+9p^2)^4
# 3p: -12i*p*(1-9p^2)/(1+9p^2)^4
# 3d: 9p^2/(1+9p^2)^4 (not implemented)
def radial_1s(r):
    return -np.power((1.+np.power(r,2)),-2)

def radial_2s(r):
    r2 = np.power(r,2)
    return (2.-8*r2)*np.power((1.+4*r2),-3)

def radial_2p(r):
    return -2j*r*np.power((1.+4*np.power(r,2)),-3)

def radial_3s(r):
    r2 = np.power(r,2)
    return -3*(1-30*r2+81*r2*r2)*np.power(1+9*r2,-4)

def radial_3p(r):
    r2 = np.power(r,2)
    return -12j*r*(-1+9*r2)*np.power(1+9*r2,-4)

def radial_3d(r):
    r2 = np.power(r,2)
    return 9*r2*np.power(1+9*r2,-4)

def radial_fn(r, n, l):
    if (n==0 or n==1):
        return radial_1s(r)
    elif (n==2):
        if (l==0):
            return radial_2s(r)
        elif (l==1):
            return radial_2p(r)
    elif (n==3):
        if (l==0):
            return radial_3s(r)
        elif (l==1):
            return radial_3p(r)
        elif (l==2):
            return radial_3d(r)
    else:
        print('arguments not supported')
        return

# analytical *single particle* energies incl. fine structure in Ry
# it would be more correct to calculate the quasiparticle energies and use those,
# but that requires the self-interaction error to be corrected and this is simpler. 
def get_energies(nwfc, wfn_nlms):
    Enmls = np.zeros(nwfc)
    for nb in range(nwfc):
        n = wfn_nlms[nb,0]
        l = wfn_nlms[nb,1]
        m = wfn_nlms[nb,2]
        ms = 0.5 - wfn_nlms[nb,3]

        j = np.abs(l + ms)
        alpha = 1/137.035999

        Enmls[nb] = -1./(n*n) * (1 + alpha**2/(n*n) * (n/(j+0.5) - 3./4.))
    return Enmls

def modify_file(nwfc, filename, wfn_nlms):
    with h5.File(filename, 'r+') as wfnfile:
        ngk = wfnfile['mf_header/kpoints/ngk'][0]
        wfnfile['mf_header/kpoints/mnband'][()] = nwfc

        del wfnfile['wfns/coeffs'] #.resize((nwfc,2,ngk,2))
        wfnfile.create_dataset('wfns/coeffs', data=np.zeros((nwfc,2,ngk,2), dtype=np.float64))

        del wfnfile['mf_header/kpoints/occ'] #.resize((1,1,nwfc))
        wfnfile.create_dataset('mf_header/kpoints/occ', data=np.zeros((1,1,nwfc)))
        wfnfile['mf_header/kpoints/occ'][0,0,0] = 1.0
    
        del wfnfile['mf_header/kpoints/el'] #.resize((1,1,nwfc))
        wfnfile.create_dataset('mf_header/kpoints/el', data=get_energies(nwfc, wfn_nlms))
        wfnfile['mf_header/kpoints/ifmin'][()] = np.ones((1,1),dtype=int)
        wfnfile['mf_header/kpoints/ifmax'][()] = np.ones((1,1),dtype=int)

def write_wfns(nwfc, filename, wfn_nlms):
    # loop over bands, for each band write wfn in chunks
    with h5.File(filename, 'r+') as wfnfile:
        gvecs = wfnfile['mf_header/gspace/components'][()]
        blat = wfnfile['mf_header/crystal/blat'][()]
        bvec = wfnfile['mf_header/crystal/bvec'][()] * blat
        celvol = wfnfile['mf_header/crystal/celvol'][()]
        occs = wfnfile['mf_header/kpoints/occ'][()]
        ngk = wfnfile['mf_header/kpoints/ngk'][:]
        nbnd = occs.shape[2]

        chunk_size = 50000

        for nb in range(nwfc):
            n = wfn_nlms[nb,0]
            l = wfn_nlms[nb,1]
            m = wfn_nlms[nb,2]
            ms = wfn_nlms[nb,3]
            print(f'n = {n},\tl = {l},\tm = {m},\tms = {(-ms+0.5):+.1f}')

            for chunk_start in range(0,ngk[0],chunk_size):
                chunk_end = min(chunk_start+chunk_size, ngk[0])
                chsize = chunk_end - chunk_start
                tmp_vecs = np.zeros((chsize,3))
                tmp_pthetaphi = np.zeros((chsize,3))

                # cart->sph conversion
                tmp_vecs = gvecs[chunk_start:chunk_end]
                tmp_vecs = np.einsum('ij,kj->kj', bvec, tmp_vecs)
                tmp_pthetaphi[:,0] = np.linalg.norm(tmp_vecs, axis=1)
                tmp_pthetaphi[:,1] = np.arctan2(np.linalg.norm(tmp_vecs[:,0:2], axis=1), tmp_vecs[:,2])
                tmp_pthetaphi[:,2] = np.arctan2(tmp_vecs[:,1],tmp_vecs[:,0])

                # spherical harmonics
                wfn = sph_harm(m,l,tmp_pthetaphi[:,2],tmp_pthetaphi[:,1])
                # radial function
                R = radial_fn(tmp_pthetaphi[:,0], n, l)
                # wavefunction
                wfn = wfn * R
                wfnfile['wfns/coeffs'][nb,ms,chunk_start:chunk_end,0] = wfn.real
                wfnfile['wfns/coeffs'][nb,ms,chunk_start:chunk_end,1] = wfn.imag
                wfnfile['wfns/coeffs'][nb,(ms+1)%2,chunk_start:chunk_end,0] = 0.0
                wfnfile['wfns/coeffs'][nb,(ms+1)%2,chunk_start:chunk_end,1] = 0.0

        # normalize wavefunctions
        for nb in range(nwfc):
            Anml2 = 0. # sum for normalization 
            for chunk_start in range(0,ngk[0],chunk_size):
                chunk_end = min(chunk_start+chunk_size, ngk[0])
                chsize = chunk_end - chunk_start
                tmp_wfnvals = np.zeros((2,chsize,2))
                tmp_wfnvals = wfnfile['wfns/coeffs'][nb,:,chunk_start:chunk_end,:]
                Anml2 += np.sum(np.power(tmp_wfnvals,2))        
            for chunk_start in range(0,ngk[0],chunk_size):
                chunk_end = min(chunk_start+chunk_size, ngk[0])
                chsize = chunk_end - chunk_start
                wfnfile['wfns/coeffs'][nb,:,chunk_start:chunk_end,:] *= 1/np.sqrt(Anml2)

    print(f'Wrote {filename} with {nwfc} wavefunctions.')



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python write_hydrogen_wfns.py [filename] [nbnd < 20]")
        sys.exit(1)
    filename = sys.argv[1]
    nwfc = int(sys.argv[2])

    # array of quantum numbers (n,l,m,ms=0up,1down) for each wfn:
    # goes up to 3p orbitals; can be extended
    wfn_nlms = np.asarray([
        [1,0,0,0],
        [1,0,0,1],
        [2,0,0,0],
        [2,0,0,1],
        [2,1,-1,0],
        [2,1,-1,1],
        [2,1,0,0],
        [2,1,0,1],
        [2,1,1,0],
        [2,1,1,1],
        [3,0,0,0],
        [3,0,0,1],
        [3,1,-1,0],
        [3,1,-1,1],
        [3,1,0,0],
        [3,1,0,1],
        [3,1,1,0],
        [3,1,1,1],
        ],np.int32)
    
    modify_file(nwfc, filename, wfn_nlms)
    write_wfns(nwfc, filename, wfn_nlms)