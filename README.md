# ECCE Jet Properties
This code intends to study properties of jet detection in the [ECCE](https://www.ecce-eic.org/) detector at the EIC.  

## Key Properties
* Jet Energy
    * Scale
        * Using the definition (Reconstructed Jet Energy - Truth Jet Energy) / Truth Jet Energy, what is the scaling of jet energy in reconstruction?  
    * Resolution
        * Similar to the scale, what is the resolution achieved in reconstruction?  This is using the definition RMS ((Reconstructed Jet Energy - Truth Jet Energy) / Truth Jet Energy)
    * Both these metrics make reference to figure 8.43 in the [Yellow Report](https://arxiv.org/abs/2103.05419)

* Jet Spatially
    * Scale
    * Resolution
    * These metrics are functionally equivalent to those described for jet energy, only now we are interested in the jets eta and phi component.  
    * Instead of dividing the truth value, these measures are simply the mean and RMS of the difference of the values

* Efficiency
    * How good is the jet reconstruction/matching? 
    * Efficiency = (Num Matched Jets) / (Num Truth Jets), binned over energy
    
    