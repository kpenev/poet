#include "IOColumns.h"

const std::vector<std::string> &output_column_names()
{
	static std::vector<std::string> names(OutCol::NUM_OUTPUT_QUANTITIES, "");
	if(names[OutCol::AGE]=="") {
		names[OutCol::AGE]="t";
		names[OutCol::SEMIMAJOR]="a";
		names[OutCol::CONV_INCLINATION]="convincl";
		names[OutCol::CONV_PERIAPSIS]="convperi";
		names[OutCol::RAD_INCLINATION]="radincl";
		names[OutCol::RAD_PERIAPSIS]="radperi";
		names[OutCol::WORB]="Worb";
		names[OutCol::PORB]="Porb";
		names[OutCol::LORB]="Lorb";
		names[OutCol::LCONV]="Lconv";
		names[OutCol::LRAD]="Lrad";
		names[OutCol::LTOT]="L";
		names[OutCol::ICONV]="Iconv";
		names[OutCol::IRAD]="Irad";
		names[OutCol::ITOT]="I";
		names[OutCol::WSURF]="Wsurf";
		names[OutCol::WRAD]="Wrad";
		names[OutCol::PSURF]="Psurf";
		names[OutCol::PRAD]="Prad";
		names[OutCol::EVOL_MODE]="mode";
		names[OutCol::WIND_STATE]="wind";
		names[OutCol::RSTAR]="R";
		names[OutCol::LSTAR]="Lum";
		names[OutCol::RRAD]="Rrad";
		names[OutCol::MRAD]="Mrad";
		names[OutCol::ICONV_DERIV]="DIconv";
		names[OutCol::IRAD_DERIV]="DIrad";
		names[OutCol::ITOT_DERIV]="DI";
		names[OutCol::RSTAR_DERIV]="DR";
		names[OutCol::RRAD_DERIV]="DRrad";
		names[OutCol::MRAD_DERIV]="DMrad";
		names[OutCol::ICONV_SECOND_DERIV]="DDIconv";
		names[OutCol::IRAD_SECOND_DERIV]="DDIrad";
		names[OutCol::ITOT_SECOND_DERIV]="DDI";
		names[OutCol::RRAD_SECOND_DERIV]="DDRrad";
#ifdef COLUMN_NAME_EMPHASIS
		for(int i=0; i<OutCol::NUM_OUTPUT_QUANTITIES; ++i)
			names[i]=(std::string(COLUMN_NAME_EMPHASIS)
					  +
					  names[i]
					  +
					  std::string(COLUMN_NAME_EMPHASIS));
#endif
	}
	return names;
}

