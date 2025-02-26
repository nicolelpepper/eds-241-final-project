/*
9/15/20 
Ashley Larsen

Code to replicate analysis figures in the main text
*/

# delimit ;

capture log close;
set more off;
*set CD to wherever the dta file is stored;
cd /Users/[...];
log using NatSusFinal.log, text replace;

*****************;
*Figure 2*;
*****************;
use NatSusRepeat, clear;

eststo Pool_Hac: ols_spatial_HAC IHSInsectperHaA IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA ConleyConst, lat(centroid_y) lon(centroid_x) t(year) p(rMTR) dist(2.5) lag(0) bartlett ;
predict XB, xb;
gen resids = IHSInsectperHaA-XB;
gen ResidSq = resids^2;
gen SemiVar= 0;
replace SemiVar = resids^2 if resids>0;
ihstrans ResidSq SemiVar, prefix(IHS);	

quietly eststo Pool_VarHac: ols_spatial_HAC IHSResidSq IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(rMTR) dist(2.5) lag(0) bartlett;

quietly eststo Pool_SemiHac: ols_spatial_HAC IHSSemiVar IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(rMTR) dist(2.5) lag(0) bartlett;

drop  XB resids ResidSq SemiVar IHSResidSq IHSSemiVar;


local Fixed "TY SY CY PY";
		foreach FE in `Fixed' {;
	*default for reg2hdfespatial is bartlett;
	eststo `FE'_Hac: reg2hdfespatial IHSInsectperHaA IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA, lat(centroid_y) lon(centroid_x) t(year) p(`FE') dist(2.5) lag(0);
	predict XB, xb;
	gen resids = IHSInsectperHaA-XB;
	gen ResidSq = resids^2;
	gen SemiVar= 0;
	replace SemiVar = resids^2 if resids>0;
	ihstrans ResidSq SemiVar, prefix(IHS);

	quietly eststo `FE'_VarHac: reg2hdfespatial IHSResidSq IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(`FE') dist(2.5) lag(0) ;
	
	quietly eststo `FE'_SemiHac: reg2hdfespatial IHSSemiVar IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(`FE') dist(2.5) lag(0) ;

	
	drop XB resids ResidSq SemiVar IHSResidSq IHSSemiVar;

};

	
	eststo SyP_Hac: reg2hdfespatial  IHSInsectperHaA IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(PY) dist(2.5) lag(0) altfetime(SpYr) ;
	predict XB, xb;
	gen resids = IHSInsectperHaA-XB;
	gen ResidSq = resids^2;
	gen SemiVar= 0;
	replace SemiVar = resids^2 if resids>0;
	ihstrans ResidSq SemiVar, prefix(IHS);

	quietly eststo SyP_VarHac: reg2hdfespatial  IHSResidSq IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(PY) dist(2.5) lag(0) altfetime(SpYr);
	
	quietly eststo SyP_SemiHac: reg2hdfespatial  IHSSemiVar IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(PY) dist(2.5) lag(0) altfetime(SpYr) ;

	
	drop IHSInsectperHaA IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA XB resids ResidSq SemiVar IHSResidSq IHSSemiVar;


*LEVEL;
coefplot Pool_Hac, bylabel(Pool)  || TY_Hac, bylabel(TY) || SY_Hac, bylabel(SY) || CY_Hac, bylabel(CY) || PY_Hac, bylabel(FY) || SyP_Hac, bylabel(FS*Y) ||, bycoefs 
keep(IHSSimpsDivNearTempA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large))  subtitle(SDI, size(large) ring(0) pos(1)) ylabel(-0.5 (0.5) 0.5, nogrid) yscale(range(-0.5 0.5) noextend) legend(off);
graph save Hac_Level_SDI.gph, replace;
coefplot Pool_Hac, bylabel(Pool)  || TY_Hac, bylabel(TY) || SY_Hac, bylabel(SY) || CY_Hac, bylabel(CY) || PY_Hac, bylabel(FY) || SyP_Hac, bylabel(FS*Y) ||, bycoefs 
keep(IHSAgHaNearTempA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) 
subtitle(Extent, size(large) ring(0) pos(1)) ylabel(-0.5 (0.5) 0.5, nogrid) yscale(range(-0.5 0.5) noextend) legend(off) title(Level, size(large) ring(1) pos(12));
graph save Hac_Level_Extent.gph, replace;
coefplot Pool_Hac, bylabel(Pool)  || TY_Hac, bylabel(TY) || SY_Hac, bylabel(SY) || CY_Hac, bylabel(CY) || PY_Hac, bylabel(FY) || SyP_Hac, bylabel(FS*Y) ||, bycoefs 
 keep(IHSHaPrmtA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel(,labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) subtitle(Field Sz, size(large) ring(0) pos(1)) 
ylabel(-0.5 (0.5) 0.5, nogrid) yscale(range(-0.5 0.5) noextend) legend(off);
graph save Hac_Level_Field.gph, replace;
graph combine  Hac_Level_Extent.gph  Hac_Level_SDI.gph Hac_Level_Field.gph, xcommon row(3) scheme(lean1) imargin(vsmall) saving(Hac_Level_AI.gph, replace);	
	
*VARIANCE;
coefplot Pool_VarHac, bylabel(Pool)  || TY_VarHac, bylabel(TY) || SY_VarHac, bylabel(SY) || CY_VarHac, bylabel(CY) || PY_VarHac, bylabel(FY) || SyP_VarHac, bylabel(FS*Y) ||, bycoefs 
keep(IHSSimpsDivNearTempA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large))  subtitle("", size(large) ring(0) pos(1)) ylabel(-0.1 (0.1) 0.2, nogrid) yscale(range(-0.1 0.2) noextend) legend(off);
graph save Hac_Var_SDI.gph, replace;
coefplot Pool_VarHac, bylabel(Pool)  || TY_VarHac, bylabel(TY) || SY_VarHac, bylabel(SY) || CY_VarHac, bylabel(CY) || PY_VarHac, bylabel(FY) || SyP_VarHac, bylabel(FS*Y) ||, bycoefs 
 keep(IHSAgHaNearTempA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) 
subtitle("", size(large) ring(0) pos(1)) ylabel(-0.1 (0.1) 0.2, nogrid) yscale(range(-0.1 0.2) noextend)  legend(off) title(Variance, size(large) ring(1) pos(12));
graph save Hac_Var_Extent.gph, replace;
coefplot Pool_VarHac, bylabel(Pool)  || TY_VarHac, bylabel(TY) || SY_VarHac, bylabel(SY) || CY_VarHac, bylabel(CY) || PY_VarHac, bylabel(FY) || SyP_VarHac, bylabel(FS*Y) ||, bycoefs 
 keep(IHSHaPrmtA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel(,labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) subtitle("", 
 size(large) ring(0) pos(1)) ylabel(-0.1 (0.1) 0.2, nogrid) yscale(range(-0.1 0.2) noextend)  legend(off);
graph save Hac_Var_Field.gph, replace;
graph combine  Hac_Var_Extent.gph  Hac_Var_SDI.gph Hac_Var_Field.gph, xcommon row(3) scheme(lean1) imargin(vsmall) saving(Hac_Var_AI.gph, replace);	
	
*SEMI-VARIANCE;
coefplot Pool_SemiHac, bylabel(Pool)  || TY_SemiHac, bylabel(TY) || SY_SemiHac, bylabel(SY) || CY_SemiHac, bylabel(CY) || PY_SemiHac, bylabel(FY) || SyP_SemiHac, bylabel(FS*Y) ||, bycoefs 
keep(IHSSimpsDivNearTempA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large))  subtitle("", size(large) ring(0) pos(1)) ylabel(-0.1 (0.1) 0.1, nogrid) yscale(range(-0.1 0.1) noextend) legend(off);
graph save Hac_Semi_SDI.gph, replace;
coefplot Pool_SemiHac, bylabel(Pool)  || TY_SemiHac, bylabel(TY) || SY_SemiHac, bylabel(SY) || CY_SemiHac, bylabel(CY) || PY_SemiHac, bylabel(FY) || SyP_SemiHac, bylabel(FS*Y) ||, bycoefs 
 keep(IHSAgHaNearTempA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) 
subtitle("", size(large) ring(0) pos(1)) ylabel(-0.1 (0.1) 0.1, nogrid) yscale(range(-0.1 0.1) noextend) legend(off) title(Semi-variance, size(large) ring(1) pos(12));
graph save Hac_Semi_Extent.gph, replace;
coefplot Pool_SemiHac, bylabel(Pool)  || TY_SemiHac, bylabel(TY) || SY_SemiHac, bylabel(SY) || CY_SemiHac, bylabel(CY) || PY_SemiHac, bylabel(FY) || SyP_SemiHac, bylabel(FS*Y) ||, bycoefs 
 keep(IHSHaPrmtA) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel(,labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) subtitle("", size(large) ring(0) pos(1)) 
ylabel(-0.1 (0.1) 0.1, nogrid) yscale(range(-0.1 0.1) noextend) legend(off);
graph save Hac_Semi_Field.gph, replace;
graph combine  Hac_Semi_Extent.gph  Hac_Semi_SDI.gph Hac_Semi_Field.gph, xcommon row(3) scheme(lean1) imargin(vsmall) saving(Hac_Semi_AI.gph, replace);	
	
graph combine  Hac_Level_AI.gph  Hac_Var_AI.gph Hac_Semi_AI.gph, xcommon row(1) scheme(lean1) imargin(vsmall) saving(Hac_AI.gph, replace);	
graph export Fig2_Hac_AI.eps, replace;


********************;
*2) Figure 3;
*******************;

use NatSusRepeat, clear;

	eststo Pool_Hac_D: ols_spatial_HAC IHSInsectperHaA IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(rMTR) dist(2.5) lag(0) bartlett;

	predict XB, xb;
	gen resids = IHSInsectperHaA-XB;
	gen ResidSq = resids^2;
	
	gen SemiVar= 0;
	replace SemiVar = resids^2 if resids>0;
	ihstrans ResidSq  SemiVar, prefix(IHS);
	
	eststo Pool_Hac_Var_D: ols_spatial_HAC IHSResidSq IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(rMTR) dist(2.5) lag(0) bartlett;

	eststo Pool_Hac_Semi_D: ols_spatial_HAC IHSSemiVar IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(rMTR) dist(2.5) lag(0) bartlett;
	

drop XB resids ResidSq SemiVar IHSResidSq  IHSSemiVar ;


**FEs***;
*all include year;
local Fixed "TY SY CY PY";
	foreach FE in `Fixed' {;
		
	eststo `FE'_Hac_Lev_D: reg2hdfespatial IHSInsectperHaA IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA, lat(centroid_y) lon(centroid_x) t(year) p(`FE') dist(2.5) lag(0);

	predict XB, xb;
	gen resids = IHSInsectperHaA-XB;
	gen ResidSq = resids^2;
	gen SemiVar= 0;
	replace SemiVar = resids^2 if resids>0;
	ihstrans ResidSq  SemiVar, prefix(IHS);
	
	eststo `FE'_Hac_D: reg2hdfespatial IHSResidSq IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(`FE') dist(2.5) lag(0);

	eststo `FE'_Hac_Semi_D: reg2hdfespatial IHSSemiVar IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(`FE') dist(2.5) lag(0);
	

drop XB resids ResidSq SemiVar IHSResidSq IHSSemiVar;
};

*Sp*yr;

	eststo SyP_Hac_Lev_D: reg2hdfespatial IHSInsectperHaA IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(PY) dist(2.5) lag(0) altfetime(SpYr);

	predict XB, xb;
	gen resids = IHSInsectperHaA-XB;
	gen ResidSq = resids^2;
	gen SemiVar= 0;
	replace SemiVar = resids^2 if resids>0;
	ihstrans ResidSq SemiVar, prefix(IHS);
	
	eststo SyP_Hac_Var_D: reg2hdfespatial IHSResidSq IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(PY) dist(2.5) lag(0) altfetime(SpYr);

	eststo SyP_Hac_Semi_D: reg2hdfespatial IHSSemiVar IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500 IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 IHSHaPrmtA,  lat(centroid_y) lon(centroid_x) t(year) p(PY) dist(2.5) lag(0) altfetime(SpYr);
	 
drop XB resids ResidSq  SemiVar IHSResidSq IHSSemiVar IHSInsectperHaA;

	
*level;
coefplot SyP_Hac_Lev_D, keep(IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500) scheme(lean1) yline(0) norecycle nooffsets vertical  msize(large) subtitle(SDI,size(large) ring(0) pos(1))
 xlabel(,labsize(large) angle(45)) ylabel(,labsize(large)) ylabel(-.1 (.1) .1, nogrid) yscale(range(-.1 .1) noextend); 
graph save Hac_Lev_D_SDI.gph, replace;
coefplot   SyP_Hac_Lev_D, keep(IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 ) scheme(lean1) yline(0) vertical norecycle nooffsets msize(large) subtitle(Extent,size(large) ring(0) pos(1)) xlabel("")  ylabel(,labsize(large))  ylabel(-.1 (.1) .1, nogrid) yscale(range(-.1 .1) noextend) title(Level, size(large) ring(1) pos(12));
graph save Hac_Lev_D_Extent.gph, replace;
graph combine  Hac_Lev_D_Extent.gph Hac_Lev_D_SDI.gph , xcommon row(2) scheme(lean1) imargin(vsmall) saving(Hac_Lev_D_AI.gph, replace);
*variance;
coefplot SyP_Hac_Var_D, keep(IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500) scheme(lean1) yline(0) norecycle nooffsets vertical  msize(large) ylabel(,labsize(large)) ylabel(-.05 (.05) .05, nogrid) yscale(range(-.05 .05) noextend) xlabel(,labsize(large) angle(45)) ;
graph save Hac_Var_D_SDI.gph, replace;
coefplot  SyP_Hac_Var_D, keep(IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 ) scheme(lean1) yline(0) vertical norecycle nooffsets msize(large) xlabel("") ylabel(,labsize(large))  ylabel(-.05 (.05) .05, nogrid) yscale(range(-.05 .05) noextend) title(Variance, size(large) ring(1) pos(12));
graph save Hac_Var_D_Extent.gph, replace;
graph combine Hac_Var_D_Extent.gph Hac_Var_D_SDI.gph , xcommon  row(2) scheme(lean1) imargin(vsmall) saving(Hac_Var_D_AI.gph, replace);
*SemiVar;
coefplot  SyP_Hac_Semi_D, keep(IHSSDIA500 IHSSDIA1000 IHSSDIA1500 IHSSDIA2000 IHSSDIA2500) scheme(lean1) yline(0) norecycle nooffsets vertical  msize(large) xlabel(,labsize(large) angle(45))  ylabel(,labsize(large)) ylabel(-.05 (.05) .05, nogrid) yscale(range(-.05 .05) noextend);
graph save Hac_Semi_D_SDI.gph, replace;
coefplot  SyP_Hac_Semi_D, keep(IHSAgHaNearA500 IHSAgHaNearA1000 IHSAgHaNearA1500 IHSAgHaNearA2000 IHSAgHaNearA2500 ) xlabel("") scheme(lean1) yline(0) vertical  norecycle nooffsets msize(large) ylabel(,labsize(large)) ylabel(-.05 (.05) .05, nogrid) yscale(range(-.05 .05) noextend) title(Semi-variance, size(large) ring(1) pos(12))  ;
graph save Hac_Semi_D_Extent.gph, replace;
graph combine Hac_Semi_D_Extent.gph Hac_Semi_D_SDI.gph, xcommon row(2) scheme(lean1) imargin(vsmall) saving(Hac_Semi_D_AI.gph, replace) ;
graph combine Hac_Lev_D_AI.gph Hac_Var_D_AI.gph Hac_Semi_D_AI.gph, xcommon row(1) scheme(lean1) imargin(vsmall) saving(Hac_AI_D.gph, replace);
graph export Fig3_AI_D.eps, replace;


***************************;
******Figure 4_ nnual, individual crop*******;
****************************/;
* almond=3001, grape=29141, orange = 2006, pistachio =3011,tangerine 2008, alfalfa = 23001, cotton 29121; 


	local crops "3001 29141 2006 3011 2008  23001 29121 ";
		foreach cr in `crops' {;
	
	use NatSusRepeat, clear;
	drop if rcommodity!=`cr';
	
	hdfe IHSInsectperHaA IHSAgHaNearTempA IHSSimpsDivNearTempA IHSHaPrmtA, absorb(year) gen (Yr);
	rename YrIHSInsectperHaA IHSInsectperHaB;
	rename YrIHSHaPrmtA IHSHaPrmtB;
	rename YrIHSSimpsDivNearTempA IHSSimpsDivNearTempB;
	rename YrIHSAgHaNearTempA IHSAgHaNearTempB;
	
	label var IHSAgHaNearTempB "Extent";
	label var IHSSimpsDivNearTempB "SDI";
	label var IHSHaPrmtB "Field Sz";
	
	
	eststo Crop_`cr'_Lev: ols_spatial_HAC IHSInsectperHaB IHSAgHaNearTempB IHSSimpsDivNearTempB IHSHaPrmtB ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(CY) dist(2.5) lag(0) bartlett;
	predict XB, xb;
	gen resids = IHSInsectperHaB-XB;
	gen ResidSq = resids^2;
	gen SemiVar= 0;
	replace SemiVar = resids^2 if resids>0;
	ihstrans ResidSq SemiVar, prefix(IHS);

	quietly eststo Crop_`cr'_Var: ols_spatial_HAC IHSResidSq IHSAgHaNearTempB IHSSimpsDivNearTempB IHSHaPrmtB ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(CY) dist(2.5) lag(0) bartlett;
	
	quietly eststo Crop_`cr'_Semi: ols_spatial_HAC IHSSemiVar IHSAgHaNearTempB IHSSimpsDivNearTempB IHSHaPrmtB ConleyConst,  lat(centroid_y) lon(centroid_x) t(year) p(CY) dist(2.5) lag(0) bartlett;

	drop IHSInsectperHaB IHSAgHaNearTempB IHSSimpsDivNearTempB IHSHaPrmtB XB resids ResidSq SemiVar IHSResidSq IHSSemiVar;
	};
		
			
*level;
coefplot Crop_3001_Lev, bylabel(Alm)  || Crop_29141_Lev, bylabel(Grp) || Crop_2006_Lev, bylabel(Org) || Crop_3011_Lev, bylabel(Pis) || Crop_2008_Lev, bylabel(Tan) || Crop_23001_Lev, bylabel(Alf) || Crop_29121_Lev, bylabel(Cot)||, bycoefs keep(IHSAgHaNearTempB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) title(Level, size(large) ring(1) pos(12))
subtitle(Extent, size(large) ring(0) pos(1)) ylabel(-0.5 (0.5) 1, nogrid) yscale(range(-0.5 1) noextend) legend(off);
graph save Crop_Level_Extent.gph, replace;
coefplot Crop_3001_Lev, bylabel(Alm)  || Crop_29141_Lev, bylabel(Grp) || Crop_2006_Lev, bylabel(Org) || Crop_3011_Lev, bylabel(Pis) || Crop_2008_Lev, bylabel(Tan) || Crop_23001_Lev, bylabel(Alf) || Crop_29121_Lev, bylabel(Cot)||, 
bycoefs keep(IHSSimpsDivNearTempB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) 
subtitle(SDI, size(large) ring(0) pos(1)) ylabel(-0.5 (0.5) 1, nogrid) yscale(range(-0.5 1) noextend)legend(off);
graph save Crop_Level_SDI.gph, replace;
coefplot Crop_3001_Lev, bylabel(Alm)  || Crop_29141_Lev, bylabel(Grp) || Crop_2006_Lev, bylabel(Org) || Crop_3011_Lev, bylabel(Pis) || Crop_2008_Lev, bylabel(Tan) || Crop_23001_Lev, bylabel(Alf) || Crop_29121_Lev, bylabel(Cot)||, 
bycoefs keep(IHSHaPrmtB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel(,labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) subtitle(Field Sz, size(large) ring(0) pos(1)) 
ylabel(-0.5 (0.5) 1, nogrid) yscale(range(-0.5 1) noextend) legend(off);
graph save Crop_Level_Field.gph, replace;
graph combine  Crop_Level_Extent.gph Crop_Level_SDI.gph  Crop_Level_Field.gph, xcommon row(3) scheme(lean1) imargin(vsmall) saving(Crop_Level_AI.gph, replace);	
*Variance;
coefplot Crop_3001_Var, bylabel(Alm)  || Crop_29141_Var, bylabel(Grp) || Crop_2006_Var, bylabel(Org) || Crop_3011_Var, bylabel(Pis) || Crop_2008_Var, bylabel(Tan) || Crop_23001_Var, bylabel(Alf) || Crop_29121_Var, bylabel(Cot)||, 
bycoefs keep(IHSAgHaNearTempB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) title(Variance, size(large) ring(1) pos(12))
subtitle("", size(large) ring(0) pos(1)) ylabel(-0.4 (0.4) 0.4, nogrid) yscale(range(-0.4 0.4) noextend) legend(off);
graph save Crop_Var_Extent.gph, replace;
coefplot Crop_3001_Var, bylabel(Alm)  || Crop_29141_Var, bylabel(Grp) || Crop_2006_Var, bylabel(Org) || Crop_3011_Var, bylabel(Pis) || Crop_2008_Var, bylabel(Tan) || Crop_23001_Var, bylabel(Alf) || Crop_29121_Var, bylabel(Cot)||, bycoefs keep(IHSSimpsDivNearTempB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) 
subtitle("", size(large) ring(0) pos(1)) ylabel(-0.4 (0.4) 0.4, nogrid) yscale(range(-0.4 0.4) noextend) legend(off);
graph save Crop_Var_SDI.gph, replace;
coefplot Crop_3001_Var, bylabel(Alm)  || Crop_29141_Var, bylabel(Grp) || Crop_2006_Var, bylabel(Org) || Crop_3011_Var, bylabel(Pis) || Crop_2008_Var, bylabel(Tan) || Crop_23001_Var, bylabel(Alf) || Crop_29121_Var, bylabel(Cot)||,
bycoefs keep(IHSHaPrmtB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel(,labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) subtitle("", size(large) ring(0) pos(1)) 
ylabel(-0.4 (0.4) 0.4, nogrid) yscale(range(-0.4 0.4) noextend) legend(off);
graph save Crop_Var_Field.gph, replace;
graph combine Crop_Var_Extent.gph Crop_Var_SDI.gph  Crop_Var_Field.gph, xcommon row(3) scheme(lean1) imargin(vsmall) saving(Crop_Var_AI.gph, replace);	
*Semi-Variance;
coefplot Crop_3001_Semi, bylabel(Alm)  || Crop_29141_Semi, bylabel(Grp) || Crop_2006_Semi, bylabel(Org) || Crop_3011_Semi, bylabel(Pis) || Crop_2008_Semi, bylabel(Tan) || Crop_23001_Semi, bylabel(Alf) || Crop_29121_Semi, bylabel(Cot)||, bycoefs keep(IHSAgHaNearTempB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) title(Semi-variance, size(large) ring(1) pos(12)) subtitle("", size(large) ring(0) pos(1)) ylabel(-0.3 (0.3) 0.3, nogrid) yscale(range(-0.3 0.3) noextend) legend(off);
graph save Crop_Semi_Extent.gph, replace;
coefplot Crop_3001_Semi, bylabel(Alm)  || Crop_29141_Semi, bylabel(Grp) || Crop_2006_Semi, bylabel(Org) || Crop_3011_Semi, bylabel(Pis) || Crop_2008_Semi, bylabel(Tan) || Crop_23001_Semi, bylabel(Alf) || Crop_29121_Semi, bylabel(Cot)||, bycoefs keep(IHSSimpsDivNearTempB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel("",labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) subtitle("", size(large) ring(0) pos(1)) ylabel(-0.3 (0.3) 0.3, nogrid) yscale(range(-0.3 0.3) noextend) legend(off);
graph save Crop_Semi_SDI.gph, replace;
coefplot Crop_3001_Semi, bylabel(Alm)  || Crop_29141_Semi, bylabel(Grp) || Crop_2006_Semi, bylabel(Org) || Crop_3011_Semi, bylabel(Pis) || Crop_2008_Semi, bylabel(Tan) || Crop_23001_Semi, bylabel(Alf) || Crop_29121_Semi, bylabel(Cot)||, bycoefs keep(IHSHaPrmtB) vertical scheme(lean1) yline(0) norecycle nooffsets  msize(large) xlabel(,labsize(large) angle(45)) ylabel(,labsize(large)) byopts(row(3) legend(off)) subtitle("", size(large) ring(0) pos(1)) ylabel(-0.3 (0.3) 0.3, nogrid) yscale(range(-0.3 0.3) noextend) legend(off);
graph save Crop_Semi_Field.gph, replace;
graph combine   Crop_Semi_Extent.gph Crop_Semi_SDI.gph Crop_Semi_Field.gph, xcommon row(3) scheme(lean1) imargin(vsmall) saving(Crop_Semi_AI.gph, replace);	
graph combine Crop_Level_AI.gph Crop_Var_AI.gph Crop_Semi_AI.gph,  xcommon ycommon row(1) scheme(lean1) imargin(vsmall) saving(Crop_AI.gph, replace);
graph export Fig4_Crop_AI.eps, replace;



log close;

