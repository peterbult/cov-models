// complex covariance norm scaling
// 
// Parameters:
//      ReIm        Switch function output between Real (1) or Imag mode (2)
//      phi_a       phase of the norm scaling
//

#include <cmath>
#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

extern "C"
void lmodcomplex(const RealArray& energyArray, const RealArray& parameters,
	         int spectrumNumber, RealArray& fluxArray, 
                 RealArray& fluxErrArray, const string& initString)
{
    // Resize the output flux array to theV proper size
    fluxArray.resize(energyArray.size()-1);

    // Unpack the model parameters
    int real_switch = parameters[0];
    Real phi_a = parameters[1] * M_PI / 180.0;

    // Preprocess Real/Imag switch
    if (real_switch == 0) {
        if (FunctionUtility::inXFLT(spectrumNumber, "reim")) {
            real_switch = FunctionUtility::getXFLT(spectrumNumber, "reim");
        }
        else {
            std::string message = "WARNING: complex: could not read REIM from file: assuming 1 [REAL mode]";
            FunctionUtility::xsWrite(message,1);
            real_switch = 1;
        }
    }

    // Construct a scaling vector for the Real/Imag output modes
    Real mult_scale = 1;
    if (real_switch == 1) {
        mult_scale = cos(phi_a);
    }
    else {
        mult_scale = sin(phi_a);
    }

    // Populate the output flux array
    for (size_t i=0; i<fluxArray.size(); i++) {
        fluxArray[i] = mult_scale;
    }

    return;
}

