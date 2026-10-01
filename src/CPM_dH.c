// file: cpm_dh.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
double calcdH(VOX* pv, NOD* pn, int* csize, int xt, int xs, int pick, int ttag, int stag)
{
	double dH, dHcontact, dHvol, dHstr;

	// go through all components of the hamiltonian

	// contact (aka adhesion) contribution
	dHcontact = 0;
	dHcontact = calcdHcontact(pv,xt,ttag,stag);

	// volume constraint contribution
	dHvol = 0;
	dHvol = calcdHvol(csize,ttag,stag);

	// durotaxis contribution
	dHstr = 0;
	dHstr = calcdHstrain(pn,xt,xs,pick,ttag,stag);

	dH = dHcontact+dHvol+dHstr;
	return dH;


}




////////////////////////////////////////////////////////////////////////////////
double calcdHcontact(VOX* pv, int xt, int ttag, int stag)
{
	// only access Moore neighborhood
	double dHcontact, Hcontact, Hcontactn;
	int nbs[8],n,nbtag;

	Hcontact=0; Hcontactn=0;
	nbs[0]=xt-1+NVX; nbs[1]=xt+NVX; nbs[2]=xt+1+NVX;
	nbs[7]=xt-1;                    nbs[3]=xt+1;
	nbs[6]=xt-1-NVX; nbs[5]=xt-NVX; nbs[4]=xt+1-NVX;
	for(n=0;n<8;n++)
	{
		nbtag = pv[nbs[n]].ctag;
		Hcontact += contactenergy(ttag,nbtag);
		Hcontactn += contactenergy(stag,nbtag);
	}
	dHcontact = Hcontactn-Hcontact;

	return dHcontact;
}

////////////////////////////////////////////////////////////////////////////////
double contactenergy(int tag1, int tag2)
{
	// set energy based on cell-cell or cell-medium adhesion
	double J;

	J = 0;
	if(tag1!=tag2)
	{
    	if((tag1==0)||(tag2==0))
        	J = JCM; // cell-medium
    	else
        	J = JCC; // cell-cell
	}
	return J;
}





////////////////////////////////////////////////////////////////////////////////
double calcdHvol(int* csize, int ttag, int stag)
{
	// for MPI, will need to comm csize for cells shared between ranks 
	// so that all ranks have an accurate view of the cell volume
	// or keep copy of entire csize array on each rank since it isn't a lot of memory
	
	double dHvol,dHvolA,dHvolB,V0,eV,eVn;
	int V;

	// assume 2 cells, A (ttag) and B (stag) are involved
	dHvolA=0; dHvolB=0; V0=TARGETVOLUME;
	if(ttag) // cell ttag retracts
	{
		V=csize[ttag-1]; eV=(V-V0)/V0; eVn=(V-1-V0)/V0;
		dHvolA = INELASTICITY*(eVn*eVn-eV*eV);
	}
	if(stag) // cell stag expands
	{
		V=csize[stag-1]; eV=(V-V0)/V0; eVn=(V+1-V0)/V0;
		dHvolB = INELASTICITY*(eVn*eVn-eV*eV);
	}
	dHvol = dHvolA+dHvolB;

	return dHvol;
}



////////////////////////////////////////////////////////////////////////////////
double calcdHstrain(NOD* pn, int xt, int xs, int pick, int ttag, int stag)
{
	double dHstrain;
	double q,vmx[8],vmy[8],vm[2];
	double estrains[3],L1,L2,v1[2],v2[2];
	double vmv1,vmv2;
	double E1,E2;
	int xtx,xty,xsx,xsy;

	dHstrain = 0;

	// unitvectors for move: vm
	q = SQ05;
	vmx[0]=-q; vmx[1]= 0; vmx[2]= q;
	vmx[7]=-1;            vmx[3]= 1;
	vmx[6]=-q; vmx[5]= 0; vmx[4]= q;

	vmy[0]= q; vmy[1]= 1; vmy[2]= q;
	vmy[7]= 0;            vmy[3]= 0;
	vmy[6]=-q; vmy[5]=-1; vmy[4]=-q;

	vm[0] = vmx[pick];
	vm[1] = vmy[pick];


	if(stag) // expansion
	{
		xty = xt / NVX;
		xtx = xt % NVX;
		get_estrains(pn,xtx,xty,estrains);
		L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);

		// inproducts of move vector with pr. strain vectors
		vmv1 = vm[0]*v1[0]+vm[1]*v1[1];
		vmv2 = vm[0]*v2[0]+vm[1]*v2[1];

		//dHstrain -= sige(L1)*vmv1*vmv1 + sige(L2)*vmv2*vmv2;
		//dHstrain -= sige(YOUNGS);

		E1 = YOUNGS; if(L1>0) {E1*=(1+L1/STIFFENINGSTIFF);}
		E2 = YOUNGS; if(L2>0) {E2*=(1+L2/STIFFENINGSTIFF);}
		dHstrain -= sige(E1)*vmv1*vmv1 + sige(E2)*vmv2*vmv2;
	}

	if(ttag) // retraction
	{
		xsy = xs / NVX;
		xsx = xs % NVX;
		get_estrains(pn,xsx,xsy,estrains);
		L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);

		// inproducts of move vector with pr. strain vectors
		vmv1 = vm[0]*v1[0]+vm[1]*v1[1];
		vmv2 = vm[0]*v2[0]+vm[1]*v2[1];

		//dHstrain += sige(L1)*vmv1*vmv1 + sige(L2)*vmv2*vmv2;
		//dHstrain += sige(YOUNGS);

		E1 = YOUNGS; if(L1>0) {E1*=(1+L1/STIFFENINGSTIFF);}
		E2 = YOUNGS; if(L2>0) {E2*=(1+L2/STIFFENINGSTIFF);}
		dHstrain += sige(E1)*vmv1*vmv1 + sige(E2)*vmv2*vmv2;
	}


	return dHstrain;
}



////////////////////////////////////////////////////////////////////////////////
double sige(double L)
{
	double x,sigL;

	x = STIFFSENSITIVITY*(L - THRESHOLDSTIFF);
	sigL = MAXDHSTR/(1+exp(-x)); //sigmoid function

	return sigL;
}



/*
////////////////////////////////////////////////////////////////////////////////
double calcdHstrain(NOD* pn, int xt, int xs, int pick, int ttag, int stag)
{
	double dHstrain;
	double q,vmx8[8],vmy8[8],vmx,vmy;
	double estrains[3],exx,eyy,exy,em;

	dHstrain = 0;
	//if ((old==0)||(new==0))
	{
		// unitvectors for move: vm
		q = SQ05;
		vmx8[0]=-q; vmx8[1]= 0; vmx8[2]= q;
		vmx8[7]=-1;             vmx8[3]= 1;
		vmx8[6]=-q; vmx8[5]= 0; vmx8[4]= q;

		vmy8[0]= q; vmy8[1]= 1; vmy8[2]= q;
		vmy8[7]= 0;             vmy8[3]= 0;
		vmy8[6]=-q; vmy8[5]=-1; vmy8[4]=-q;

		vmx = vmx8[pick];
		vmy = vmy8[pick];

		if(stag) // expansion
		{
			get_estrains(pn,xt,estrains);
			exx=estrains[0]; eyy=estrains[1]; exy=.5*estrains[2];
			em = vmx*(exx*vmx+exy*vmy)+vmy*(exy*vmx+eyy*vmy);
			dHstrain -= sige(em);
		}
		if(ttag) // retraction
		{
			get_estrains(pn,xs,estrains);
			exx=estrains[0]; eyy=estrains[1]; exy=.5*estrains[2];
			em = vmx*(exx*vmx+exy*vmy)+vmy*(exy*vmx+eyy*vmy);
			dHstrain += sige(em);
		}
	}
	return dHstrain;
}
*/


