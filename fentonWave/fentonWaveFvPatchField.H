/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    fentonWaveFvPatchField

Description
	This is an implementation of the Fourier approximation method for finding
	stream function waves described in section 3.1.2 of the paper by J.D.
	Fenton, "Numerical Methods for Nonlinear Waves" (1999), Numerical Methods 
	for Nonlinear Waves, Advances in Coastal and Ocean Engineering, Vol. 5,
	pp241-324. Equation numbers in the code refer to this paper. At the time
	of writing the paper could be downloaded here:

	http://www.johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf	

	Usage: 
	
	The wave generator patch velocity should be specified in 0/U like

		patchName            
		{
			type            fentonWave;
			value			uniform (0 0 0);
		}

	Similarly the dynamic pressure in 0/pd and volumen fraction field in 
	0/alpha1 must look like:

		patchName            
		{
			type            fentonWave;
			value			uniform 0;
		}

	In the <case>/constant folder there should be a file named waveProperties
	containing the wave parameters:

	//	Name			Value ex.	Description
		height          0.1;		//Wave heigth
		period          20;			//Wave period (only used if uType = 1 or 2)
		length			2;			//Wave length (only used if uType = 0)
		uType           2; 			//Input type. 0: wave lenght specified. 1 or 2: wave period and current specified (default = 2).
		u1or2           0; 			//Current, Eulerian if uType = 1, Stokes if uType = 2 (default = 0).

		depth           0.3;		//Water depth
		direction       ( 1 0 0 ); 	//Direction of wave vector
		rampUpTime      0;			//Initial wave height ramp up time
		seabedHeight	-0.3; 		//Vertical coordinate of seabed on patch (default = lowest point on patch)
		waveMakerPos	-3.0; 		//Horizontal position of wave maker (default = backmost point on patch)
		phase			.1;			//Phase of wave, theta = k*x+omega*t + 2*pi*phase (default = 0.0)

		nModes          40; 		//Number of Fourier modes in expansion (default = 30)
		nHeightSteps    20; 		//Number of height steps in search for wave solution (default = 20)
		rootTolerance	1e-7; 		//Tolerance in multidimensional root finding algorithm (default = 1e-7)
		rootMaxIter		1000; 		//Maximum number of iterations in multidimensional root finding (default = 1000)
	
	This folder must also contain a file named g containing the gravitation 
	vector, e.g.:
	
		dimensions      [0 1 -2 0 0 0 0];
		value           ( 0 -9.81 0 );

	The file <case>/system/controlDict must end with the line:
	
		libs ( "libFentonWaveBC.so" "libOpenFOAM.so" );	
		
SourceFiles
    fentonWaveFvPatchField.C

Author
	Johan Roenby, DHI Water & Environment.

\*---------------------------------------------------------------------------*/

#ifndef fentonWaveFvPatchField_H
#define fentonWaveFvPatchField_H

#include "mixedFvPatchFields.H"
#include "uniformDimensionedFields.H" //gravity type of field
//#include "dimensionedScalar.H"
//#include "dimensionedVector.H"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fentonWaveFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fentonWaveFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

		//- wave property dictionary
		IOdictionary waveProp_;

		//- Wave Height
		scalar H_;

		//- Mean water level
        scalar d_;
				
		//- Current
		scalar u1or2_;

		//- Current type
		int uType_; //0 for wavelength input, 1 for wave period and "Eulerian" current, 2 for wave period and "Stokes" current
		
		//- Number of Fourier modes
        scalar N_;

		//- Number of steps to reach final wave height
		int nHsteps_;

		//- gravity vector
        vector g_;

		//- Wave vector
        vector K_;
		
		//- Vertical coordinate of seabed at patch
		scalar seabedHeight_;

		//- Horizontal position of wave maker patch
		scalar xWaveMaker_;
		
		//- Wave phase
        scalar phi_;

		//- ramp up time
		scalar rampUpTime_;		
		
		//- Water density
		scalar rho1_;

		//- Air density
		scalar rho2_;

		//- Wave number k = 2*pi/lambda
        scalar k_;
								
		//- Angular wave frequency
        scalar omega_;
		
		//- Mean fluid speed in co-moving frame
        scalar uBar_;
		
		//- Bernoulli constant
		scalar R_;
						
		//- Fourier coefficients for velocity field and pressure
        scalarField B_;
		
		//- Fourier coefficients for surface elevation
        scalarField E_;
		
		//- Residual tolerance in multidimensional root finding 
		scalar rootTolerance_;
		
		//- Max number of iteration in multidimensional root finding
		unsigned int rootMaxIter_;

		//- Take a wild guess
		scalar pi_;
	
		struct fentonParms{ double h; double lort; double u1or2; int uType; int N; };

    // Protected member functions
		
		//Initialize member data from dictionary (avoid constructor code duplication)
		void init();
		
		//Initialize member data from other patch (avoid constructor code duplication)
		void init(const fentonWaveFvPatchField<Type>&);

		//Calculate Fourier coefficients for stream function wave
		void fenton();

		//Multidimensional root finding using the GSL library
		void findRoots(scalarField &, scalar, scalar);
		
		//Vector function that is the zero vector for a stream function wave solution
		static int waveConditionVectorFunction (const gsl_vector *, void *, gsl_vector *);
						
		//Populate B_ and E_ either by calculation or by reading from file
		void getFourierCoeffs();
		
		void writeOutMembers();

public:

    //- Runtime type information
    TypeName("fentonWave");


    // Constructors

        //- Construct from patch and internal field
        fentonWaveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fentonWaveFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fentonWaveFvPatchField
        //  onto a new patch
        fentonWaveFvPatchField
        (
            const fentonWaveFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new fentonWaveFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fentonWaveFvPatchField
        (
            const fentonWaveFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new fentonWaveFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions
				
		//- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

		//- Update the coefficients associated with the patch field
        virtual void evaluate();

        //- Set field values on patch
		virtual void setField
		(
			const Field<Type>&, 
			const scalarField&, 
			const scalarField&, 
			const scalarField&
		);
		
		//- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fentonWaveFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
