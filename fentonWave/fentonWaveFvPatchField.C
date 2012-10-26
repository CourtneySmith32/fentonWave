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

\*---------------------------------------------------------------------------*/

#include "fentonWaveFvPatchField.H"
#include "surfaceFields.H" //using lookupPatchField to get flux

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fentonWaveFvPatchField<Type>::fentonWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),	
	waveProp_
	(
        IOobject
        (
            "waveProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    N_(waveProp_.lookupOrDefault<int>("nModes",30)),
    B_(N_,0.0),
	E_(N_+1,0.0)
{
	Info<< "Construct from patch and internal field" << endl;
	init();
	writeOutMembers();

	//Initialising to fixed value boundary with zero as the value
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1.0;
}

template<class Type>
fentonWaveFvPatchField<Type>::fentonWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
	waveProp_
	(
        IOobject
        (
            "waveProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    N_(waveProp_.lookupOrDefault<int>("nModes",30)),
    B_(N_,0.0),
	E_(N_+1,0.0)
{
	Info<< "Construct from patch, internal field and dictionary" << endl;
	init();
	writeOutMembers();

    if (dict.found("value")) //Not sure why this is needed but it exists in many BC's derived from mixed, so better keep it (JRO).
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
	
}

template<class Type>
fentonWaveFvPatchField<Type>::fentonWaveFvPatchField
(
    const fentonWaveFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
	waveProp_(ptf.waveProp_),
	B_(ptf.B_),
	E_(ptf.E_)
{
	Info<< "Construct by mapping given fentonWaveFvPatchField onto a new patch" << endl;
	init(ptf);
}


template<class Type>
fentonWaveFvPatchField<Type>::fentonWaveFvPatchField
(
    const fentonWaveFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
	waveProp_(ptf.waveProp_),
	B_(ptf.B_),
	E_(ptf.E_)
{
	Info<< "Construct as copy setting internal field reference" << endl;
	init(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
template<class Type>
void fentonWaveFvPatchField<Type>::updateCoeffs()
{

	if (this->updated())
    {
        return;
    }
	
	scalarField x = (K_/mag(K_)) & (this->patch().Cf());
	x -= xWaveMaker_;
	scalarField y = (-g_/mag(g_)) & (this->patch().Cf());
	y -= seabedHeight_; //Vertical coordinate as measured from seabed
	
	scalar t = this->db().time().value();
	scalarField theta = k_*x - omega_*t + 2*pi_*phi_;

	scalarField eta(x.size(),0.0);
	forAll(E_,jj)
	{
		eta += 2*E_[jj]*cos(jj*theta);
	}
	
	scalarField u(y.size(),0.0);
	scalarField v(y.size(),0.0);
	int	jj;
	forAll(B_,ii)
	{
		jj = ii + 1;
		u += jj*B_[ii]*cosh(jj*k_*y)/cosh(jj*k_*d_)*cos(jj*theta);
		v += jj*B_[ii]*sinh(jj*k_*y)/cosh(jj*k_*d_)*sin(jj*theta);
	}
	u *= sqrt(mag(g_)/k_);
	v *= sqrt(mag(g_)/k_);
	u += omega_/k_ - uBar_;
	
	scalar fac(1.0);
	if (t < rampUpTime_)
	{
		fac = sin(pi_/2.0*t/rampUpTime_);
		u *= fac;
		v *= fac;
		eta = fac*(eta-d_) + d_;
	}

	scalarField alpha(neg(y - eta));

	const word& fieldName = this->dimensionedInternalField().name(); //For some reason this line cannot be executed in setField. Therefore fieldName is written here and given as an input parameter to setField.
	setField(this->refValue(),fieldName, alpha, u, v);
    mixedFvPatchField<Type>::updateCoeffs(); //This line simply runs fvPatchField::updateCoeffs() which simplys sets its private boolean updated_ = true;

}

//In the following we use template specialization to set field depending on its Type: 
//Inspired by http://www.parashift.com/c++-faq-lite/templates.html#faq-35.7

template<class Type>
void fentonWaveFvPatchField<Type>::setField(Field<Type> const& refVal, const word& fieldName, const scalarField& alpha, const scalarField& u, const scalarField& v)
{
//   dummy code executed if setField is called with other Type than scalar or vector
}   

// Setting volume fraction (alpha1) and pressure (pd) field
template<> 
void fentonWaveFvPatchField<scalar>::setField(Field<scalar> const& iF, const word& fieldName, const scalarField& alpha, const scalarField& u, const scalarField& v)
{
	if (fieldName == "alpha1")
	{			
		this->refValue() = alpha;
		this->refGrad() = 0.0;
		this->valueFraction() = 1.0;
				
		//Setting field to fixed value at inflow faces and zero gradient at outflow faces
		if (this->db().time().timeIndex() != 0)
		{
			const Field<scalar>& phip = this->patch().lookupPatchField
			(
				"phi",
				reinterpret_cast<const surfaceScalarField*>(NULL),
				reinterpret_cast<const scalar*>(NULL)
			);

			this->valueFraction() = 1.0 - pos(phip);
		}
		else
		{
			this->valueFraction() = 1.0;
		}

	}
	else if (fieldName == "pd" || fieldName == "p_rgh")
	{
//		scalar c = omega_/k_;
//		this->refValue() =  alpha*rho_*( R_ - 0.5*( pow(u-c,2) + pow(v,2) ) );
//		this->refGrad() = 0.0;
//		this->valueFraction() = 1.0;
		
		//So far only zero gradient implemented
		this->refValue() =  0.0;
		this->refGrad() = 0.0;
		this->valueFraction() = 0.0;
	}	
}

// Setting velocity field
template<> 
void fentonWaveFvPatchField<vector>::setField(Field<vector> const& refVal, const word& fieldName, const scalarField& alpha, const scalarField& u, const scalarField& v)
{			
	this->refValue() = alpha*( u*K_/mag(K_) + v*( -g_/mag(g_)) );
	this->refGrad() = vector::zero;
	this->valueFraction() = 1.0;
} 

// Initialising 
template<class Type>
void fentonWaveFvPatchField<Type>::init()
{

	//Populating member data
	pi_ = acos(-1.0);
	H_ = readScalar(waveProp_.lookup("height"));
	d_ = readScalar(waveProp_.lookup("depth"));
	nHsteps_ = waveProp_.lookupOrDefault<int>("nHeightSteps",20);
	uType_ = waveProp_.lookupOrDefault<int>("uType",2);
	rootTolerance_ = waveProp_.lookupOrDefault<scalar>("rootTolerance",1e-7);
	rootMaxIter_ = waveProp_.lookupOrDefault<int>("rootMaxIter",1000);

	if ( uType_ == 0 )
	{
		scalar lambda = readScalar(waveProp_.lookup("length"));
		k_ = 2*pi_/lambda;
		omega_ = 0.0; //co-moving frame
		//u1or2 is not used for anything in this case
	}
	else if ( (uType_ == 1) || (uType_ == 2) )
	{
		//k_ is calculated by fenton()
		scalar tau = readScalar(waveProp_.lookup("period"));
		omega_ = 2*pi_/tau;
		u1or2_ = waveProp_.lookupOrDefault("u1or2",0.0);
	}

	g_ = 
	uniformDimensionedVectorField
	(
		IOobject
		(
			"g",
			this->db().time().constant(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	).value();

	K_ = vector(waveProp_.lookup("direction"));
	K_ = K_/mag(K_);

	//Vertical seabed coordinate at boundary
	scalar yMin = min( (-g_/mag(g_)) & (this->patch().boundaryMesh().mesh().points())); //lowest point on patch
	seabedHeight_ = waveProp_.lookupOrDefault<scalar>("seabedHeight",yMin);

	//Horizontal wave maker position
	scalar xMin = min( K_ & (this->patch().boundaryMesh().mesh().points())); //lowest point on patch
	xWaveMaker_ = waveProp_.lookupOrDefault<scalar>("waveMakerPos",xMin);
	
    phi_ = waveProp_.lookupOrDefault<scalar>("phase",0.0);
	rampUpTime_ = readScalar(waveProp_.lookup("rampUpTime"));

	IOdictionary tp
	(
		IOobject
		(
			"transportProperties",
			this->db().time().constant(),
			this->db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
	dictionary sD(tp.subDict("phase1"));
	rho_ = (dimensionedScalar(sD.lookup("rho"))).value();	

	getFourierCoeffs();
}

template<class Type>
void fentonWaveFvPatchField<Type>::init(const fentonWaveFvPatchField<Type>& ptf)
{
	xWaveMaker_ = ptf.xWaveMaker_;
	seabedHeight_ = ptf.seabedHeight_;
	pi_ = ptf.pi_;
	H_ = ptf.H_;
	d_ = ptf.d_;
	u1or2_ = ptf.u1or2_;
	uType_ = ptf.uType_;
	g_ = ptf.g_;
	K_ = ptf.K_;
	phi_ = ptf.phi_;
	N_ = ptf.N_;
	nHsteps_ = ptf.nHsteps_;
	rampUpTime_ = ptf.rampUpTime_;
	rho_ = ptf.rho_;
	k_ = ptf.k_;
	omega_ = ptf.omega_;
	uBar_ = ptf.uBar_;
	R_ = ptf.R_;
	rootTolerance_ = ptf.rootTolerance_;
	rootMaxIter_ = ptf.rootMaxIter_;
}

template<class Type>
void fentonWaveFvPatchField<Type>::getFourierCoeffs()
{
	IOdictionary fentonFile
	(
		IOobject
		(
			this->patch().name(),
			this->db().time().constant(),
			this->db(),
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		)
	);
	
	bool sameWave = !(fentonFile.empty());

	if ( sameWave )
	{
		scalar H(readScalar(fentonFile.lookup("H")));
		scalar d(readScalar(fentonFile.lookup("d")));
		int N(readScalar(fentonFile.lookup("N")));
		int nHsteps(readScalar(fentonFile.lookup("nHsteps")));
		scalar rootTolerance(readScalar(fentonFile.lookup("tolerance")));
		unsigned int rootMaxIter(readInt(fentonFile.lookup("maxIter")));
		int uType(readInt(fentonFile.lookup("uType")));

		scalar prec(2*pow(10.0,-scalar(IOstream::precision_)));

		sameWave = sameWave && (N == N_);
		sameWave = sameWave && (nHsteps == nHsteps_);
		sameWave = sameWave && (abs(rootTolerance - rootTolerance_) < prec);
		sameWave = sameWave && (rootMaxIter == rootMaxIter_);
		sameWave = sameWave && (uType == uType_);
		sameWave = sameWave && (mag(H/d - H_/d_) < prec);

		scalar k(readScalar(fentonFile.lookup("k")));
		scalar omega(readScalar(fentonFile.lookup("omega")));
		scalar u1or2(readScalar(fentonFile.lookup("u1or2")));		
			
		if ( uType_ == 0 )
		{
			sameWave = sameWave && (mag(k/d - k_/d_) < prec);
		}
		else if ( uType_ == 1 || uType_ == 2 )
		{
			sameWave = sameWave && (mag(omega - omega_) < prec);
			sameWave = sameWave && (mag(u1or2 - u1or2_) < prec);
		}
		else
		{
			Info << "uType has to be 0, 1 or 2!" << endl;		
		}

		if ( sameWave )
		{
			k_ = k;
			omega_ = omega;
			uBar_ = readScalar(fentonFile.lookup("uBar"));  
			R_ = readScalar(fentonFile.lookup("R"));
			B_ = fentonFile.lookup("B");  
			E_ = fentonFile.lookup("E");		
		}
	}
	
	if ( ! sameWave )
	{
		fenton();
		//write parameters defining a unique wave calculation to file
		fentonFile.add("H",H_,true);
		fentonFile.add("d",d_,true);
		fentonFile.add("N",N_,true);
		fentonFile.add("nHsteps",nHsteps_,true);
		fentonFile.add("maxIter",rootMaxIter_,true);
		fentonFile.add("tolerance",rootTolerance_,true);
		fentonFile.add("uType",uType_,true);
		fentonFile.add("k",k_,true);
		fentonFile.add("omega",omega_,true);
		fentonFile.add("u1or2",u1or2_,true);
		fentonFile.add("R",R_,true);
		fentonFile.add("uBar",uBar_,true);
		fentonFile.add("B",B_,true);
		fentonFile.add("E",E_,true);
		fentonFile.writeObject(IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED);
	}
}

template<class Type>
void fentonWaveFvPatchField<Type>::fenton()
{
	Info << "Calculating stream function wave..." << endl;

	scalar lort, g = mag(g_);

	if (uType_ == 0)
	{
		//k_ set to 2*pi_/lambda in init();
		scalar lambda = 2*pi_/k_;
		lort = lambda/d_; //Nondimensionalised wave length
		omega_ = 0.0;
		//u1or2_ not used in this case
	} else
	{
		scalar tau = 2*pi_/omega_;
		lort = tau * Foam::sqrt(g/d_); //Nondimensionalised wave period T*sqrt(g/d)
		k_ = Foam::pow(omega_,2.0)/ ( g * Foam::pow( Foam::tanh( Foam::pow( omega_*Foam::sqrt(d_/g) , 3.0/2.0 ) ) , 2.0/3.0) );
		u1or2_ /= sqrt(g*d_);
	}
	
//	Info << "Specifying initial guess for variables from 1st order theory (Fenton 1999, p. 14)" << endl; 
	scalar kd = k_*d_; //Mean depth normalized by wavenumber, k*d.
	
	scalarField keta(N_ + 1); //Nondimensionalised surface elevation above sea level, keta = k*eta
    forAll(keta,ii)
    {
		keta[ii] = kd + 0.5 * k_ * (H_/nHsteps_) * Foam::cos(ii*pi_/N_);
    }
	scalar u = Foam::sqrt(Foam::tanh(kd)); //Nondimensionalized mean horizontal velocity in co-moving frame, Ubar*sqrt(k/g).
	B_ = 0.0; //Fourier coefficients
	B_[0] = 0.5 * k_ * (H_/nHsteps_) / u;
	scalar q = kd * u; //Nondimensionalized volume flux in co-moving frame, q = Q*sqrt(k^3/g).
	scalar r = kd + 0.5 * Foam::pow(u,2.0); //Nondimensionalized constant in Bernoulli equation, R*k/g.
//	Info << "g = " << g << ", k_ guess = " << k_ << ", kd guess = " << kd << ", u = " << u << ", B[0] = " << B[0] << ", q = " << q << ", r = " << r << endl;

//	Info << "Putting guessed values into x vector" << endl;
	scalarField x(2*N_+5);
    forAll(keta,ii)
    {
		x[ii] = keta[ii];
    }
    forAll(B_,ii)
    {
		x[ii + N_ + 1] = B_[ii];
    }
	x[2*N_+1] = u; 
	x[2*N_+2] = kd;
	x[2*N_+3] = q;
	x[2*N_+4] = r;

	for (int ii=1; ii<=nHsteps_; ii++) 
	{
		findRoots(x, ii*H_/nHsteps_/d_, lort);
		Info << "Iteration " << ii << " of " << nHsteps_ << " with H/d = " << ii*H_/nHsteps_/d_ << " gave kd = " << x[2*N_+2] << endl;
	}

	//Converting x vector into named variables
    forAll(keta,ii)
    {
		keta[ii] = x[ii];
    }
    forAll(B_,ii)
    {
		B_[ii] = x[ii + N_ + 1];
    }
	u = x[2*N_+1];
	kd = x[2*N_+2];
	q = x[2*N_+3];
	r = x[2*N_+4];

	k_ = kd/d_;
	K_ = K_/mag(K_)*k_;
	R_ = r*g/k_;

	//Calculating E_ using eq. 3.26 (with j running from 0 to N) and eq. 3.27 (with the RHS divided by N to correct the equation) in Fenton 1999.
	//With these corrections A_[j] = 2/N*E_{j+1} for j = [0,N-2], and A_[N-1] = E_{N}/N.
	scalarField etapp = keta/k_;
	etapp[0] *= 0.5;
	etapp[N_] *= 0.5;
	forAll(E_,jj)
	{
		forAll(etapp,mm)
		{
			E_[jj] += etapp[mm] * Foam::cos(jj*mm*pi_/N_);
		}
	}
	
	E_[0] *= 0.5;
	E_[N_] *= 0.5;
	E_ /= N_;

	uBar_ = u / Foam::sqrt(k_/g);
		
	Info << "Done calculating stream function wave" << endl;
}

template<class Type>
int fentonWaveFvPatchField<Type>::waveConditionVectorFunction(const gsl_vector * x, void *params, gsl_vector * f)
{
//    Info << "Enter waveConditionVectorFunction" << endl;
	//Reading parameters
	struct fentonParms *p = static_cast<struct fentonParms*>(params);
	double h = p->h;
	double lort = p->lort;
	double u1or2 = p->u1or2;
	int uType = p->uType;
	int N = p->N;
//	Info << "h: " << h << ", lort: " << lort << ", u1or2: " << u1or2 << ", uType: " << uType << ", N: " << N << endl;

//	Info << "Converting from x vector into named variables" << endl;
	scalarField keta(N+1);
    forAll(keta,ii)
    {
		keta[ii] = gsl_vector_get(x,ii);
    }
	scalarField B(N);
    forAll(B,ii)
    {
		int jj = ii + N + 1;
		B[ii] = gsl_vector_get(x,jj);
    }
	scalar u = gsl_vector_get(x,2*N+1); 
	scalar kd = gsl_vector_get(x,2*N+2);
	scalar q = gsl_vector_get(x,2*N+3);
	scalar r = gsl_vector_get(x,2*N+4); 

//	Info << "Computing wave conditions, N = " << N << endl;
	scalarField kinBC(N+1); //Kinematic boundary condition, (3.7) in Fenton 1999
	scalarField dynBC(N+1); //Dynamic boundary condition (3.8) in Fenton 1999
	scalar U, V, kx, ke, PI = acos(-1.0);
	int j;
	scalar meanDepthEq = - (0.5/N)*(keta[0] + keta[N]) - kd; //Mean depth equation (3.9) in Fenton 1999
	forAll(kinBC,ii)
    {
		kx = ii*PI/N;
		ke = keta[ii];
		kinBC[ii] = - u * ke + q;
		dynBC[ii] = ke - r;
		U = - u;
		V = 0.0;
		forAll(B,jj)
		{
			j = jj + 1;
			kinBC[ii] += B[jj] * (Foam::sinh(j*ke) / Foam::cosh(j*kd)) * Foam::cos(j*kx);
			U += j * B[jj] * (Foam::cosh(j*ke) / Foam::cosh(j*kd)) * Foam::cos(j*kx);
			V += j * B[jj] * (Foam::sinh(j*ke) / Foam::cosh(j*kd)) * Foam::sin(j*kx);
		}
		dynBC[ii] += 0.5*Foam::pow(U,2.0) + 0.5*Foam::pow(V,2.0);
		meanDepthEq += keta[ii]/N;
    }
	scalar waveHeightEq = keta[0] - keta[N] - kd*h; //Wave height equation (3.10) in Fenton 1999
	
	scalar lastEq = 0;
	if (uType == 0) {
		scalar l = lort;
		lastEq = kd - 2.0*PI/l; // Relative wavelength equation (3.11) in Fenton 1999
	}
	else {
		scalar t = lort;
		if (uType == 1) {
			lastEq = Foam::sqrt(kd)*u + kd*u1or2 - 2.0*PI/t; //Equation for mean current in laboratory frame (3.15) in Fenton 1999
		}
		else if (uType == 2) {
			lastEq = q/Foam::sqrt(kd) + kd*u1or2 - 2.0*PI/t; //Equation for mass transport velocity (3.16) in Fenton 1999
		}
	}	

//	Info << "Writing wave conditions into solution vector f" << endl;
	forAll(kinBC,ii)
    {
		gsl_vector_set (f, ii, kinBC[ii]);
		gsl_vector_set (f, ii + N + 1, dynBC[ii]);
    }
	gsl_vector_set (f, 2*N+2, meanDepthEq);
	gsl_vector_set (f, 2*N+3, waveHeightEq);
	gsl_vector_set (f, 2*N+4, lastEq);
	return GSL_SUCCESS;
}

template<class Type>
void fentonWaveFvPatchField<Type>::findRoots(scalarField & X, scalar h, scalar lort)
{

//	Code adapted example from:
//	http://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Multidimensional-Root-finding.html

//	Info << "Enter findRoots" << endl;
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s; 
	
	int status;
	size_t iter = 0;

	const size_t n = X.size();
	struct fentonParms p = { h, lort, u1or2_, uType_, N_};
	gsl_multiroot_function f = {&waveConditionVectorFunction, n, &p};

	gsl_vector *x = gsl_vector_alloc (n);

	forAll(X,ii)
	{
		gsl_vector_set (x, ii, X[ii]);
	}
	T = gsl_multiroot_fsolver_hybrid;
	s = gsl_multiroot_fsolver_alloc (T, n);
	gsl_multiroot_fsolver_set (s, &f, x);

	do
	 {
	   iter++;
	   status = gsl_multiroot_fsolver_iterate (s);

	   if (status)   /* check if solver is stuck */
		 break;

		status = gsl_multiroot_test_residual (s->f, rootTolerance_);
	 }
	while (status == GSL_CONTINUE && iter < rootMaxIter_);

	printf ("status = %s\n", gsl_strerror (status));

	forAll(X,ii)
	{
		X[ii] = gsl_vector_get (s->x, ii);
	}

	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);
}

// Write
template<class Type>
void fentonWaveFvPatchField<Type>::write(Ostream& os) const
{		
    fvPatchField<Type>::write(os);
    os.writeKeyword("height")
        << H_ << token::END_STATEMENT << nl;
    os.writeKeyword("depth")
        << d_ << token::END_STATEMENT << nl;
    os.writeKeyword("period")
        << 2*pi_/omega_ << token::END_STATEMENT << nl;
    os.writeKeyword("length")
        << 2*pi_/k_ << token::END_STATEMENT << nl;
    os.writeKeyword("uType")
        << uType_ << token::END_STATEMENT << nl;
    os.writeKeyword("u1or2")
        << u1or2_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}

template<class Type>
void fentonWaveFvPatchField<Type>::writeOutMembers()
{
		Info << "Wave height " << H_ << ", length " << 2*pi_/k_ << ", and period " << 2*pi_/omega_ << endl;
		Info << "Water depth: " << d_ << ", current type " << uType_ << " and current " << u1or2_ << endl; 
		Info << "Number of Fourier modes " << N_ << " and height steps " << nHsteps_ << endl;
		Info << "Root finding tolerance " << rootTolerance_ << " and max number of iterations " << rootMaxIter_ << endl;
		Info << "Gravity " << g_ << ", wave direction " << K_ << ", and ramp up time " << rampUpTime_ << endl; 
		Info << "Seabed height " << seabedHeight_ << ", wave maker position " << xWaveMaker_ << ", and wave phase " << phi_ << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
