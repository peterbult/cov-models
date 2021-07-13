// complex covariance of a blackbody spectrum
// 
// Parameters:
//      ReIm    Switch function output between Real (1) or Imag mode (2)
//      kT      mean temperature of the time-average continuum blackbody
//      phi_a   phase of the constant temperature term
//      phi_b   phase of the linear temperature term
//      rho     correlation between the constant and linear temperature terms
//

#include <cmath>
#include <xsTypes.h>
//#include <XSstreams.h>
//#include <XSUtil/Numerics/BinarySearch.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

Real covbb_density(Real e, Real kt, const RealArray& scale)
{
    Real expkt = exp(e/kt);
    Real term1 = e*e / (expkt - 1);
    Real term2 = e*e*e * expkt / (kt*kt * (expkt-1)*(expkt-1));

    return scale[0] * term1 + scale[1] * term2;
}

extern "C"
void lmodcovbb(const RealArray& energyArray, const RealArray& parameters,
	   int spectrumNumber, RealArray& fluxArray, 
           RealArray& fluxErrArray, const string& initString)
{
    // Resize the output flux array to theV proper size
    fluxArray.resize(energyArray.size()-1);

    // Unpack the model parameters
    int real_switch = parameters[0];
    Real kt = parameters[1];
    Real phi_a = parameters[2] * M_PI / 180.0;
    Real phi_b = parameters[3] * M_PI / 180.0;
    Real rho = parameters[4];

    // Preprocess Real/Imag switch
    if (real_switch == 0) {
        if (FunctionUtility::inXFLT(spectrumNumber, "reim")) {
            real_switch = FunctionUtility::getXFLT(spectrumNumber, "reim");
        }
        else {
            std::string message = "WARNING: covbb: could not read REIM from file: assuming 1 [REAL mode]";
            FunctionUtility::xsWrite(message,1);
            real_switch = 1;
        }
    }

    // Construct a scaling vector for the Real/Imag output modes
    RealArray ReIm_scale(2);
    if (real_switch == 1) {
        ReIm_scale[0] = cos(phi_a);
        ReIm_scale[1] = cos(phi_b) * rho;
    }
    else {
        ReIm_scale[0] = sin(phi_a);
        ReIm_scale[1] = sin(phi_b) * rho;
    }

    // Set the internal normalization
    Real norm = 1.0344e-3;

    // Precompute the bin weight for Simpson's rule
    const Real h(1.0/6.0);

    // Precompute the flux density at the start of the first energy bin
    Real e0 = energyArray[0];
    Real f0 = covbb_density(e0, kt, ReIm_scale);

    // Populate the output flux array
    for (size_t i=0; i<fluxArray.size(); i++) {
        // Compute the flux density at the midpoint of the energy bin
        Real e1 = (energyArray[i] + energyArray[i+1]) / 2.0;
        Real f1 = covbb_density(e1, kt, ReIm_scale);

        // Compute the flux density at the end of the energy bin
        Real e2 = energyArray[i+1];
        Real f2 = covbb_density(e2, kt, ReIm_scale);

        // Compute the integrated flux via Simpson's rule
        fluxArray[i] = norm * h * (e2-e0) * (f0 + 4*f1 + f2);

        // Memorize f2 for use in the next iteration
        e0 = e2;
        f0 = f2;
    }
    
    return;
}
