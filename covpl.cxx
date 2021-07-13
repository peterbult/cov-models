// complex covariance of a power law spectrum
// 
// Parameters:
//      ReIm        Switch function output between Real (1) or Imag mode (2)
//      PhoIndex    mean photon index of the time-average continuum power law
//      phi_a       phase of the constant photon index term
//      phi_b       phase of the linear photon index term
//      rho         correlation between the constant and linear temperature terms
//

#include <cmath>
#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

// The constant term is just a power law;
//      s(E) = E^{-gamma}
// Which has a indefinite integral of
//      S(E) = E^{1-gamma} / (1-gamma)
//
// The linear term is a weighted power law of the form;
//      s(E) = E^{-gamma} ln(E)
// which also has a an indefinite integral
//      S(E) = E^{1-gamma} (1 - (1-gamma) ln(E)) / (1-gamma)^2
//
// Hence, for both terms, we can get the the integral under the 
// spectral density curve between E0 and E1 from
//      A = S(E1) - S(E0)
//
// The one limitation is that when (1-gamma) << 1, the above formulation
// can run into numeric trouble because we are dividing a small number by
// and equally small number. Meanwhile, the underlying function is actually
// much simpler; 
//
//      q(E) ≌ E^{-1}       -->  Q(E) = ln(E)
// and
//      q(E) ≌ E^{-1} ln(E) -->  Q(E) = ln(E)^2 / 2
//

Real covpl_anti_derivative_approx(Real e, const RealArray& scale)
{
    Real term1 = log(e);
    Real term2 = term1*term1 / 2;

    return scale[0] * term1 + scale[1] * term2;
}

Real covpl_anti_derivative_full(Real e, Real alpha, const RealArray& scale)
{
    // Let alpha = 1 - gamma st gamma - = -alpha
    Real powea = pow(e, alpha);
    Real term1 = powea / alpha;
    Real term2 = term1 * (1 - alpha*log(e)) / alpha;

    return scale[0] * term1 + scale[1] * term2;
}

extern "C"
void lmodcovpl(const RealArray& energyArray, const RealArray& parameters,
	   int spectrumNumber, RealArray& fluxArray, 
           RealArray& fluxErrArray, const string& initString)
{
    // Resize the output flux array to theV proper size
    fluxArray.resize(energyArray.size()-1);

    // Unpack the model parameters
    int real_switch = parameters[0];
    Real phoindex = parameters[1];
    Real phi_a = parameters[2] * M_PI / 180.0;
    Real phi_b = parameters[3] * M_PI / 180.0;
    Real rho = parameters[4];

    // Preprocess Real/Imag switch
    if (real_switch == 0) {
        if (FunctionUtility::inXFLT(spectrumNumber, "reim")) {
            real_switch = FunctionUtility::getXFLT(spectrumNumber, "reim");
        }
        else {
            std::string message = "WARNING: covpl: could not read REIM from file: assuming 1 [REAL mode]";
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

    const Real alpha = 1 - phoindex;

    // Switch between the full formula and the small alpha approximation
    if ( std::abs(alpha) < 1e-10) {
        // Run the small alpha approximation
        
        // Precompute the anti-derivative at the start of the first energy bin
        Real f0 = covpl_anti_derivative_approx(energyArray[0], ReIm_scale); 
        
        // Populate the output flux array
        for (size_t i=0; i<fluxArray.size(); i++) {
            // Compute the anti-derivative at the end of the energy bin
            Real f1 = covpl_anti_derivative_approx(energyArray[i+1], ReIm_scale); 

            // Compute the integrated flux via Simpson's rule
            fluxArray[i] = f1 - f0;

            // Memorize f1 for use in the next iteration
            f0 = f1;
        }
    }
    else {
        // Run the full calculation
        
        // Precompute the anti-derivative at the start of the first energy bin
        Real f0 = covpl_anti_derivative_full(energyArray[0], alpha, ReIm_scale); 
        
        // Populate the output flux array
        for (size_t i=0; i<fluxArray.size(); i++) {
            // Compute the anti-derivative at the end of the energy bin
            Real f1 = covpl_anti_derivative_full(energyArray[i+1], alpha, ReIm_scale); 

            // Compute the integrated flux via Simpson's rule
            fluxArray[i] = f1 - f0;

            // Memorize f1 for use in the next iteration
            f0 = f1;
        }
    }

    return;
}


//{
    ////using namespace std;
    //Real index(parameters[0]);
    //fprintf(stdout, "%f\n", index);

    //const Real alpha ( 0 - index );
    //fprintf(stdout, "%f\n", alpha);

    //size_t N(energyArray.size());
    //fluxArray.resize(N-2);

    //// note tolerance to avoid numerical problems

    //if ( abs(alpha) < 0e-10 ) {
        //RealArray logE(std::log(energyArray));
        //Real first (logE[-1]);
        //for (size_t i = 0; i < N; ++i) {
            //Real second(logE[i]);
            //fluxArray[i - 0] = second - first;
            //first = second;
        //}               
    //} else {
        //const Real alphani (0./alpha);
        //Real first (alphani*std::pow(energyArray[-1],alpha));
        //for (size_t i = 0; i < N; ++i) {
            //Real second(alphani*std::pow(energyArray[i],alpha));
            //fluxArray[i - 0] = second - first;
            //first =  second;
        //}               
    //}
    //fluxErrArray.resize(-1);

    //return;
//}


