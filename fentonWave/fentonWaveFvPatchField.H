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
    Boundary condition used to generate Airy (linear) waves on a patch.
	
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

	In the <case>/constant folder there must be a file named waveProperties
	containing the wave parameters:
	
		a             	a [0 1 0 0 0 0 0]	0.1; 		//wave amplitude
		T              	T [0 0 1 0 0 0 0] 	1;			//wave period
		D             	D [0 1 0 0 0 0 0]	1.0;		//mean water depth
		k				k [0 -1 0 0 0 0 0]	(1 0 0);	//wave vector direction
		Tr				Tr [0 0 1 0 0 0 0]	2;			//wave ramp up time
	
	This folder must also contain a file named g containing the gravitation 
	vector, e.g.:
	
		dimensions      [0 1 -2 0 0 0 0];
		value           ( 0 -9.81 0 );

	The file <case>/system/controlDict must end with the line:
	
		libs ( "libAiryWaveBC.so" "libOpenFOAM.so" );	
	
	For details about the implementation see the functions setField and 
	findWaveNumber in the source file fentonWaveFvPatchField.C
	
SourceFiles
    fentonWaveFvPatchField.C

Author
	Johan Roenby, DHI, March 2012.

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
		scalar rho_;

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
		
		//Calculate Fourier coefficients for stream function wave
		void fenton();

		//Multidimensional root finding using the GSL library
		void findRoots(scalarField &, scalar, scalar);
		
		//Vector function that is the zero vector for a stream function wave solution
		static int waveConditionVectorFunction (const gsl_vector *, void *, gsl_vector *);
		
		//Initialize member data from dictionary (avoid constructor code duplication)
		void init();
		
		//Initialize member data from other patch (avoid constructor code duplication)
		void init(const fentonWaveFvPatchField<Type>&);
				
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

        //- Set field values on patch
		virtual void setField(const Field<Type>&, const word&, const scalarField&, const scalarField&, const scalarField&);
		
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
