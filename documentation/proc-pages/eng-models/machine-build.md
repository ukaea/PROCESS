# Radial and Vertical Build

Simplified scale diagrams of the vertical and horizontal cross-sections of the machine can be
output in the `2-page summary` using the utility `plot_proc.py` (currently stored in `process/process/io`).  

The coordinate system is $(R,Z)$ system, where $R$ is the radial distance from the vertical
centreline (axis) of the torus, and $Z$ is the vertical distance from the equatorial midplane.

Components are often referred to as being 'inboard' or 'outboard', which simply
means that they lie at a radius $R$ less than or greater than $R_0$,
respectively, where $R_0$ is the plasma major radius (`rmajor`).

The radial build is described in detail in the OUT.DAT file as in the example below, which lists
the major radius of each component in the midplane.  The machine is axisymmetric, except for the
TF coils which are discrete.  The variables marked "IP" below are input variables.  Those marked
"ITV" are available as iteration variables (although it is not always advisable to use them as
iteration variables).  Those marked with an asterisk (*) may or may not be input variables
depending on the switches used.  

```text
 ************************************************ Radial Build ***************************************
 
 Device centreline                            0.000           0.000                       
 Machine bore                                 2.124           2.124   (bore)              IP, ITV
 Central solenoid                             0.500           2.624   (ohcth)             IP, ITV
 CS precompression                            0.065           2.689   (precomp)           
 Gap                                          0.050           2.739   (gapoh)             IP, ITV
 TF coil inboard leg                          1.400           4.139   (tfcth)             IP
 Gap                                          0.050           4.189   (tftsgap)           IP
 Thermal shield, inboard                      0.050           4.239   (thshield_ib)       IP
 Gap                                          0.020           4.259   (gapds)             IP
 Vacuum vessel (and shielding)                0.600           4.859   (d_vv_in + shldith) IP
 Gap                                          0.020           4.879   (vvblgap)           IP
 Inboard blanket                              0.755           5.634   (blnkith)           IP*
 Inboard first wall                           0.018           5.652   (fwith)             
 Inboard scrape-off                           0.225           5.877   (scrapli)           IP, ITV
 Plasma geometric centre                      3.265           9.142   (rminor)            
 Plasma outboard edge                         3.265          12.408   (rminor)            
 Outboard scrape-off                          0.225          12.633   (scraplo)           IP, ITV
 Outboard first wall                          0.018          12.651   (fwoth)             
 Outboard blanket                             0.982          13.633   (blnkoth)           IP*
 Gap                                          0.020          13.653   (vvblgap)           IP
 Vacuum vessel (and shielding)                1.100          14.753   (d_vv_out+shldoth)  IP
 Gap                                          1.900          16.652   (gapsto)            
 Thermal shield, outboard                     0.050          16.702   (thshield_ob)       IP
 Gap                                          0.050          16.752   (tftsgap)           IP
 TF coil outboard leg                         1.400          18.152   (dr_tf_coil_outboard)          

```

The radial build is shown schematically below (click to zoom).

<img title="Radial build" src="../../images/radial-build.png" alt="Radial build">

The vertical build is described in detail in the OUT.DAT as in the following example, which lists
the vertical coordinate of each component at the point furthest from the midplane (excluding the
CS, the other PF coils and the cryostat).  The midplane is defined to be half way between the top
and bottom of the plasma.  A single-null scenario is assumed to have a lower divertor, in which
case the machine is not symmetric about the midplane.  

```text
 *********************************************** Vertical Build *************************
 
 Single null case
                                          Thickness (m)    Height (m)
 TF coil                                      1.576           9.862   (tfcth)             
 Gap                                          0.050           8.286   (tftsgap)           
 Thermal shield                               0.050           8.236   (thshield)          
 Gap                                          0.050           8.186   (vgap_vv_thermalshield)             
 Vacuum vessel (and shielding)                0.600           8.136   (d_vv_top+shldtth)  
 Gap                                          0.020           7.536   (vvblgap)           
 Top blanket                                  0.869           7.516   (blnktth)           
 Top first wall                               0.018           6.647   (fwtth)             
 Top scrape-off                               0.600           6.629   (vgaptop)           
 Plasma top                                   6.029           6.029   (rminor*kappa)      
 Midplane                                     0.000          -0.000                       
 Plasma bottom                                6.029          -6.029   (rminor*kappa)      
 Lower scrape-off                             2.002          -8.031   (vgap)              
 Divertor structure                           0.621          -8.652   (divfix)            
 Vacuum vessel (and shielding)                1.000          -9.652   (d_vv_bot+shldlth)  
 Gap                                          0.050          -9.702   (vgap_vv_thermalshield)             
 Thermal shield                               0.050          -9.752   (thshield)          
 Gap                                          0.050          -9.802   (tftsgap)           
 TF coil                                      1.576         -11.379   (tfcth)    

```

The vertical build is shown schematically below (click to zoom).  

<img title="Vertical build" src="../../images/vertical-build.png" alt="Vertical build">

Since PROCESS is essentially a 0-D code, the shape of each component is used to estimate its mass
and cost, but is not used otherwise.  The first wall, blanket, shield and vacuum vessel may be
either D-shaped in cross-section, or each may be defined by two half-ellipses. The choice between
these two possibilities is set using input parameter `fwbsshape`, which should be

- 1 for D-shaped,
- 2 for ellipses.

!!! Info "TF coil placement"
    The radial build can vary from the figures above dependant on the placement of the inboard TF
    coil leg when using the `tf_in_cs` switch. See [TF coil page](tf-coil.md)**
