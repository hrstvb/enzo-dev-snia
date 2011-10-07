extern "C" void FORTRAN_NAME(pgas2d_dual)(
	        float *dslice, float *eslice, float *geslice, float *pslice,
		float *uslice, float *vslice, float *wslice, float *eta1, 
		float *eta2, int *idim, int *jdim, int *i1, int *i2, 
		int *j1, int *j2, float *gamma, float *pmin);

extern "C" void FORTRAN_NAME(pgas2d)(
	        float *dslice, float *eslice, float *pslice,
		float *uslice, float *vslice, float *wslice, int *idim, 
		int *jdim, int *i1, int *i2, int *j1, int *j2, float *gamma, 
		float *pmin);

extern "C" void FORTRAN_NAME(calcdiss)(
	        float *dslice, float *eslice, float *uslice, float *v, 
		float *w, float *pslice, float *dx, float *dy, float *dz, 
		int *idim, int *jdim, int *kdim, int *i1, int *i2, int *j1, 
		int *j2, int *k, int *nzz, int *idir, int *dimx, int *dimy, 
		int *dimz, float *dt, float *gamma, int *idiff, int *iflatten, 
		float *diffcoef, float *flatten);

extern "C" void FORTRAN_NAME(inteuler)(
		float *dslice, float *pslice, int *gravity, float *grslice, 
		float *geslice, float *uslice, float *vslice, float *wslice, 
		float *dxi, float *flatten, int *idim, int *jdim, int *i1, 
		int *i2, int *j1, int *j2, int *idual, float *eta1, 
		float *eta2, int *isteep, int *iflatten,
		int *iconsrec, int *iposrec,
		float *dt, float *gamma, int *ipresfree, float *dls, float *drs, 
		float *pls, float *prs, float *gels, float *gers, float *uls, 
		float *urs, float *vls, float *vrs, float *wls, float *wrs, 
		int *ncolor, float *colslice, float *colls, float *colrs, 
		float *h1, float *h2, float *h3,  float *h4, float *h5,      // ok this is very verbose ... 
		float *h6, float *h7, float *h8,  float *h9, 
		float *h10, float *h11, float *h12, float *h13,  float *h14, 
		float *h15, float *h16, float *h17, float *h18,  float *h19, 
		float *h20, float *h21, float *h22, float *h23,  float *h24, 
		float *h25, float *h26, float *h27, float *h28,  float *h29, 
		float *h30, float *h31, float *h32, float *h33,  float *h34, 
		float *h35, float *h36, float *h37, float *h38,  float *h39, 
		float *h40, float *h41, float *h42, float *h43,  float *h44, 
		float *h45, float *h46, float *h47, float *h48,  float *h49, 
		float *h50, float *h51, float *h52, float *h53,  float *h54, 
		float *h55, float *h56, float *h57, float *h58,  float *h59, 
		float *h60, float *h61, float *h62, float *h63,  float *h64, 
		float *h65, float *h66, float *h67, float *h68,  float *h69, 
		float *h70, float *h71, float *h72, float *h73,  float *h74, 
		float *h75, float *h76, float *h77, float *h78,  float *h79, 
		float *h80, float *h81, float *h82, float *h83,  float *h84, 
		float *h85, float *h86, float *h87,
		float *h100, float *h101, float *h102, float *h103,  float *h104, 
		float *h105, float *h106, float *h107, float *h108,
		float *h200, float *h201
);

extern "C" void FORTRAN_NAME(twoshock)(
		float *dls, float *drs, float *pls, float *prs, 
		float *uls, float *urs, int *idim, int *jdim, int *i1,
		int *i2, int *j1, int *j2, float *dt, float *gamma, 
		float *pmin, int *ipresfree, float *pbar, float *ubar,
		int *gravity, float *grslice, int *idual, float *eta1,
		float *h1, float *h2, float *h3,  float *h4,  float *h5, 
		float *h6, float *h7, float *h8,  float *h9);

extern "C" void FORTRAN_NAME(flux_twoshock)(
           float *dslice, float *eslice, float *geslice, float *uslice, 
	   float *vslice, float *wslice, float *dx,
           float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
	   int *j1, int *j2, float *dt, float *gamma, int *idiff,
           int *idual, float *eta1, int *ifallback,
           float *dls, float *drs, float *pls, float *prs, 
	   float *gels, float *gers, 
           float *uls, float *urs, float *vls, float *vrs, 
	   float *wls, float *wrs, float *pbar, float *ubar,
           float *df, float *ef, float *uf, float *vf, float *wf, float *gef, 
	   float *ges,
           int *ncolor, float *colslice, float *colls, float *colrs, float *colf);

extern "C" void FORTRAN_NAME(flux_hll)(
           float *dslice, float *eslice, float *geslice, float *uslice, 
	   float *vslice, float *wslice, float *dx,
           float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
	   int *j1, int *j2, float *dt, float *gamma,
	   int *idiff, int *idual, float *eta1, int *ifallback,
           float *dls, float *drs, float *pls, float *prs, 
           float *uls, float *urs, float *vls, float *vrs, 
	   float *wls, float *wrs, float *gels, float *gers, 
           float *df, float *ef, float *uf, float *vf, float *wf, float *gef, 
	   float *ges,
           int *ncolor, float *colslice, float *colls, float *colrs, float *colf);

extern "C" void FORTRAN_NAME(flux_hllc)(
           float *dslice, float *eslice, float *geslice, float *uslice, 
	   float *vslice, float *wslice, float *dx,
           float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
	   int *j1, int *j2, float *dt, float *gamma,
	   int *idiff, int *idual, float *eta1, int *ifallback,
           float *dls, float *drs, float *pls, float *prs, 
           float *uls, float *urs, float *vls, float *vrs, 
	   float *wls, float *wrs, float *gels, float *gers, 
           float *df, float *ef, float *uf, float *vf, float *wf, float *gef, 
	   float *ges,
           int *ncolor, float *colslice, float *colls, float *colrs, float *colf);

extern "C" void FORTRAN_NAME(euler)(
		float *dslice, float *pslice, float *grslice, 
		float *geslice, float *uslice, float *vslice, float *wslice, 
		float *dx, float *diffcoef, int *idim, int *jdim, int *i1, 
		int *i2, int *j1, int *j2, float *dt, float *gamma, int *idiff, 
		int *gravity, int *idual, float *eta1, float *eta2, float *df, 
		float *ef, float *uf, float *vf, float *wf, float *gef,
		float *ges,
		int *ncolor, float *colslice, float *colf,
		float *h1, float *h2, float *h3,  float *h4,  float *h5, 
		float *h6, float *h7);

static  float h1[MAX_ANY_SINGLE_DIRECTION];
static  float h2[MAX_ANY_SINGLE_DIRECTION];
static  float h3[MAX_ANY_SINGLE_DIRECTION];
static  float h4[MAX_ANY_SINGLE_DIRECTION];
static  float h5[MAX_ANY_SINGLE_DIRECTION];
static  float h6[MAX_ANY_SINGLE_DIRECTION];
static  float h7[MAX_ANY_SINGLE_DIRECTION];
static  float h8[MAX_ANY_SINGLE_DIRECTION];
static  float h9[MAX_ANY_SINGLE_DIRECTION];
static  float h10[MAX_ANY_SINGLE_DIRECTION];
static  float h11[MAX_ANY_SINGLE_DIRECTION];
static  float h12[MAX_ANY_SINGLE_DIRECTION];
static  float h13[MAX_ANY_SINGLE_DIRECTION];
static  float h14[MAX_ANY_SINGLE_DIRECTION];
static  float h15[MAX_ANY_SINGLE_DIRECTION];
static  float h16[MAX_ANY_SINGLE_DIRECTION];
static  float h17[MAX_ANY_SINGLE_DIRECTION];
static  float h18[MAX_ANY_SINGLE_DIRECTION];
static  float h19[MAX_ANY_SINGLE_DIRECTION];
static  float h20[MAX_ANY_SINGLE_DIRECTION];
static  float h21[MAX_ANY_SINGLE_DIRECTION];
static  float h22[MAX_ANY_SINGLE_DIRECTION];
static  float h23[MAX_ANY_SINGLE_DIRECTION];
static  float h24[MAX_ANY_SINGLE_DIRECTION];
static  float h25[MAX_ANY_SINGLE_DIRECTION];
static  float h26[MAX_ANY_SINGLE_DIRECTION];
static  float h27[MAX_ANY_SINGLE_DIRECTION];
static  float h28[MAX_ANY_SINGLE_DIRECTION];
static  float h29[MAX_ANY_SINGLE_DIRECTION];
static  float h30[MAX_ANY_SINGLE_DIRECTION];
static  float h31[MAX_ANY_SINGLE_DIRECTION];
static  float h32[MAX_ANY_SINGLE_DIRECTION];
static  float h33[MAX_ANY_SINGLE_DIRECTION];
static  float h34[MAX_ANY_SINGLE_DIRECTION];
static  float h35[MAX_ANY_SINGLE_DIRECTION];
static  float h36[MAX_ANY_SINGLE_DIRECTION];
static  float h37[MAX_ANY_SINGLE_DIRECTION];
static  float h38[MAX_ANY_SINGLE_DIRECTION];
static  float h39[MAX_ANY_SINGLE_DIRECTION];
static  float h40[MAX_ANY_SINGLE_DIRECTION];
static  float h41[MAX_ANY_SINGLE_DIRECTION];
static  float h42[MAX_ANY_SINGLE_DIRECTION];
static  float h43[MAX_ANY_SINGLE_DIRECTION];
static  float h44[MAX_ANY_SINGLE_DIRECTION];
static  float h45[MAX_ANY_SINGLE_DIRECTION];
static  float h46[MAX_ANY_SINGLE_DIRECTION];
static  float h47[MAX_ANY_SINGLE_DIRECTION];
static  float h48[MAX_ANY_SINGLE_DIRECTION];
static  float h49[MAX_ANY_SINGLE_DIRECTION];
static  float h50[MAX_ANY_SINGLE_DIRECTION];
static  float h51[MAX_ANY_SINGLE_DIRECTION];
static  float h52[MAX_ANY_SINGLE_DIRECTION];
static  float h53[MAX_ANY_SINGLE_DIRECTION];
static  float h54[MAX_ANY_SINGLE_DIRECTION];
static  float h55[MAX_ANY_SINGLE_DIRECTION];
static  float h56[MAX_ANY_SINGLE_DIRECTION];
static  float h57[MAX_ANY_SINGLE_DIRECTION];
static  float h58[MAX_ANY_SINGLE_DIRECTION];
static  float h59[MAX_ANY_SINGLE_DIRECTION];
static  float h60[MAX_ANY_SINGLE_DIRECTION];
static  float h61[MAX_ANY_SINGLE_DIRECTION];
static  float h62[MAX_ANY_SINGLE_DIRECTION];
static  float h63[MAX_ANY_SINGLE_DIRECTION];
static  float h64[MAX_ANY_SINGLE_DIRECTION];
static  float h65[MAX_ANY_SINGLE_DIRECTION];
static  float h66[MAX_ANY_SINGLE_DIRECTION];
static  float h67[MAX_ANY_SINGLE_DIRECTION];
static  float h68[MAX_ANY_SINGLE_DIRECTION];
static  float h69[MAX_ANY_SINGLE_DIRECTION];
static  float h70[MAX_ANY_SINGLE_DIRECTION];
static  float h71[MAX_ANY_SINGLE_DIRECTION];
static  float h72[MAX_ANY_SINGLE_DIRECTION];
static  float h73[MAX_ANY_SINGLE_DIRECTION];
static  float h74[MAX_ANY_SINGLE_DIRECTION];
static  float h75[MAX_ANY_SINGLE_DIRECTION];
static  float h76[MAX_ANY_SINGLE_DIRECTION];
static  float h77[MAX_ANY_SINGLE_DIRECTION];
static  float h78[MAX_ANY_SINGLE_DIRECTION];
static  float h79[MAX_ANY_SINGLE_DIRECTION];
static  float h80[MAX_ANY_SINGLE_DIRECTION];
static  float h81[MAX_ANY_SINGLE_DIRECTION];
static  float h82[MAX_ANY_SINGLE_DIRECTION];
static  float h83[MAX_ANY_SINGLE_DIRECTION];
static  float h84[MAX_ANY_SINGLE_DIRECTION];
static  float h85[MAX_ANY_SINGLE_DIRECTION];
static  float h86[MAX_ANY_SINGLE_DIRECTION];
static  float h87[MAX_ANY_SINGLE_DIRECTION];
static  float h88[MAX_ANY_SINGLE_DIRECTION];
static  float h89[MAX_ANY_SINGLE_DIRECTION];

static  float h100[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h101[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h102[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h103[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h104[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h105[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h106[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h107[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h108[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];
static  float h109[MAX_ANY_SINGLE_DIRECTION*MAX_COLOR];

static  float h200[MAX_ANY_SINGLE_DIRECTION*25];
static  float h201[MAX_ANY_SINGLE_DIRECTION*25];
