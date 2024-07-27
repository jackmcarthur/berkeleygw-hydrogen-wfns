# berkeleygw-hydrogen-wfns
If you want to modify the BerkeleyGW package, it is wise to compare the outputs of new code (exchange integrals, transition matrix elements, orbital moments...) with published results. I recently wrote a routine for which the most obvious comparison was an analytical result for the hydrogen single-particle wavefunctions $\psi_{nlm\uparrow/\downarrow}$.

The DFT wavefunctions typically used as inputs in BerkeleyGW would be insufficiently accurate for this purposes due to self-interaction error, so I instead implemented the exact form of the hydrogen wavefunctions in the BerkeleyGW WFN.h5 format. As BerkeleyGW stores the wavefunctions in plane wave components, this required the momentum-space representation of the hydrogen wavefunctions:[[1]](#1)

$$\Psi_{nlm}(p,\theta_p,\phi_p)\propto (-i)^l \frac{n^l p^l}{(n^2p^2+1)^{l+2}} C_{n-l-1}^{l+1}\left(\frac{n^2p^2-1}{n^2p^2+1}\right) Y^l_m(\theta_p,\phi_p)$$

The python script in this repository can be run as `write_hydrogen_wfns.py [WFNfilename] [nbnd < 20]` and it will populate an existing BerkeleyGW WFN.h5 file with the first `nbnd` spinor wavefunctions for a hydrogen atom at the origin, in the order:
$$(n,l,m,m_s)=(1,0,0,+1/2),(1,0,0,-1/2),(2,0,0,+1/2),(2,0,0,-1/2),(2,1,-1,+1/2),(2,1,-1,-1/2),(2,1,0,+1/2)...$$.

The file `WFN_small.h5` was generated from a calculation in Quantum ESPRESSO with a 100 Ry plane-wave cutoff in a $20 \r{A} \times 20 \r{A} \times 20 \r{A}$ unit cell. These parameters provide fairly converged results up to the 3s/3p states, but it may be useful to generate another calculation with a higher cutoff or larger box for wavefunctions higher than that.


## References
<a id="1">[1]</a> 
Bransden, B. H., Joachain, C. J. (1983). Physics of atoms and molecules. pp 621. United Kingdom: Longman.
